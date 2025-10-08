// Scan Algorithm adapted from https://github.com/simveit/effective_scan
// Extended to handle arbitrary size arrays by Alexey V. Medvedev https://github.com/a-v-medvedev/low_latency_cuda

#define WARP_SIZE 32
#define LOG_WARP_SIZE 5
#define WARP_MASK (WARP_SIZE - 1)

__device__ inline int lane_id(void) { return threadIdx.x & WARP_MASK; }
__device__ inline int warp_id(void) { return threadIdx.x >> LOG_WARP_SIZE; }
// Warp scan
template <typename TYPE>
__device__ __forceinline__ TYPE warp_scan(TYPE val) {
  TYPE x = val;
#pragma unroll
  for (int offset = 1; offset < WARP_SIZE; offset <<= 1) {
    TYPE y = __shfl_up_sync(0xffffffff, x, offset);
    if (lane_id() >= offset) x += y;
  }
  return x - val;
}

template <typename TYPE, int threadsPerBlock>
__device__ TYPE block_scan(int in) {
  __shared__ TYPE sdata[threadsPerBlock >> LOG_WARP_SIZE];
  // A. Exclusive scan within each warp
  TYPE warpPrefix = warp_scan<TYPE>(in);
  // B. Store in shared memory
  if (lane_id() == WARP_SIZE - 1) sdata[warp_id()] = warpPrefix + in;
  __syncthreads();
  // C. One warp scans in shared memory
  if (threadIdx.x < WARP_SIZE)
    sdata[threadIdx.x] = warp_scan<TYPE>(sdata[threadIdx.x]);
  __syncthreads();
  // D. Each thread calculates its final value
  TYPE thread_out_element = warpPrefix + sdata[warp_id()];
  return thread_out_element;
}

/*
template <typename TYPE>
__device__ void merge_blocks_small(int tid, int block_size, int block_id, int nblocks, TYPE *output, int numElements) {
  // postprocessing: merge of scan results for multiple blocks
  namespace cg = cooperative_groups;
  cg::grid_group grid = cg::this_grid();
  // round up the nblocks to power of two
  nblocks--; nblocks |= nblocks >> 1; nblocks |= nblocks >> 2; nblocks |= nblocks >> 4; nblocks |= nblocks >> 8; nblocks++;
  // we divide all blocks into groups
  // and we gradually increase the size of the group from just 2 blocks to 
  // nblocks (means: all blocks in a single group)
  for (int group_size = 2; group_size <= nblocks; group_size *= 2) {
    grid.sync();
    if (tid < numElements) {
      int half_group_size = group_size / 2; 
      // what is the position of this block in a group?
      int position_in_group = block_id % group_size;
      // all elements in the second half of the group get an addition equal to the value of 
      // the last element of the first half of the group
      if (position_in_group >= half_group_size) {
        // idx of an an element just before 2nd half of the group start 
        int idx_of_addition = (block_id - position_in_group + half_group_size) * block_size - 1;
        output[tid] += output[idx_of_addition];
      }
    }
  }
}
*/

// NOTE: If block_id == -1 we have a special case when blocks are so big that we handle them sequentially as "chunks"
// Chunk size is assumed to be equal to block size
template <typename TYPE, int maxBlocksInGrid>
__device__ void merge_blocks(int tid, int block_size, int block_id, TYPE *output, int numElements, int nchunks = 0) {
  __shared__ TYPE sdata[maxBlocksInGrid];
  // postprocessing: merge of scan results for multiple blocks
  namespace cg = cooperative_groups;
  cg::grid_group grid = cg::this_grid();
  grid.sync();
  TYPE addition = 0;
  // pre-load into shared memory the values of all elements of the preceding blocks
  if (threadIdx.x < maxBlocksInGrid) {
    int idx_of_addition = (threadIdx.x + 1) * block_size - 1;
    if (idx_of_addition < numElements) 
      sdata[threadIdx.x] = output[idx_of_addition];
  }
  __syncthreads();
  grid.sync();
  if (block_id == -1) {
    // Assume block size is as big as a total number of threads, so we 
    // have to handle blocks (let's call them "chunks") sequentially
    for (int chunk = 1; chunk < nchunks; chunk++) {
      int gtid = tid + chunk * block_size;
      if (gtid < numElements) {
        addition = 0;
        // go through the last-elmenet values of all preceding chunks, 
        // make a sum and add it to ours value
        for (int b = 0; b < chunk; b++)
          addition += sdata[b];
        output[gtid] += addition;
      }
    }
  } else {
    if (tid < numElements && block_id) {
      // go through the last-elmenet values of all preceding blocks, 
      // make a sum and add it to ours value
      for (int b = 0; b < block_id; b++)
        addition += sdata[b];
      output[tid] += addition;
    }
  }
}

// note: static assert: threadsPerBlock >= maxBlocksInGrid
template <typename TYPE, int threadsPerBlock, int maxBlocksInGrid>
__global__ void inclusive_scan(TYPE *input, TYPE *output, int numElements) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int chunk_stride = threadsPerBlock * gridDim.x;
  int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
  for (int chunk = 0; chunk < nchunks; chunk++) {
    int gtid = tid + (chunk * chunk_stride);
    TYPE val = input[gtid];
    TYPE result = block_scan<TYPE,threadsPerBlock>(val);
    if (gtid < numElements) {
      output[gtid] = result + val;
    }
    merge_blocks<TYPE, maxBlocksInGrid>(tid, threadsPerBlock, blockIdx.x, output + (chunk * chunk_stride), numElements);
  }
  if (nchunks > 1) {
    merge_blocks<TYPE, maxBlocksInGrid>(tid, chunk_stride, -1, output, numElements, nchunks);
  }
}

/*
template <typename TYPE, int threadsPerBlock>
__global__ void inclusive_scan_small(TYPE *input, TYPE *output, int numElements) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    TYPE val = input[tid];
    TYPE result = block_scan<TYPE,threadsPerBlock>(val);
    if (tid < numElements) {
      output[tid] = result + val;
    }
    merge_blocks_small<TYPE>(tid, threadsPerBlock, blockIdx.x, gridDim.x, output, numElements);
}
*/

// Assumed gridDim.x == 1
template <typename TYPE, int threadsPerBlock>
__global__ void inclusive_scan_one_block(TYPE *input, TYPE *output, int numElements) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int chunk_stride = threadsPerBlock;
    int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
    TYPE addition = 0; 
    for (int chunk = 0; chunk < nchunks; chunk++) {
      int gtid = tid + (chunk * chunk_stride);

      TYPE val = input[gtid] + (threadIdx.x ? 0 : addition);
      TYPE result = block_scan<TYPE,threadsPerBlock>(val);
      if (tid < numElements) {
        output[gtid] = result + val;
      }
      __syncthreads();
      if (chunk != nchunks - 1)
        addition = output[(chunk + 1) * chunk_stride - 1];
   }
}

