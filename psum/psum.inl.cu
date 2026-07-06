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

// Merge up to maxBlocksInGrid number of blocks, each block is threadsPerBlock size (block_size == threadsPerBlock)
// Or: if block_id == -1, merge large chunks sequentially
template <typename TYPE, int threadsPerBlock, int maxBlocksInGrid>
__device__ void merge_blocks(int tid, int block_id, int block_size, TYPE *output, int offset, int numElements, int nchunks = 0) {
  __shared__ TYPE sdata[maxBlocksInGrid];
  // postprocessing: merge of scan results for multiple blocks
  namespace cg = cooperative_groups;
  cg::grid_group grid = cg::this_grid();
  grid.sync();
  TYPE addition = 0, value = 0;
  int idx_of_addition = (threadIdx.x + 1) * block_size - 1;
  if (threadIdx.x < maxBlocksInGrid && offset + idx_of_addition < numElements) {
    value = output[offset + idx_of_addition];
  }
  addition = block_scan<TYPE,threadsPerBlock>(value);
  if (block_id == -1) {
    if (threadIdx.x < maxBlocksInGrid) sdata[threadIdx.x] = addition + value;
  } else {
    if (threadIdx.x == block_id) sdata[0] = addition;
  }
  __syncthreads();
  if (block_id == -1) {
    for (int chunk = 1; chunk < nchunks; chunk++) {
      int gtid = tid + chunk * block_size;
      addition = sdata[chunk-1];
      if (gtid < numElements) output[gtid] += addition;
    }
  } else {
    addition = sdata[0];
    int gtid = offset + tid;
    if (block_id && gtid < numElements) output[gtid] += addition;
  }
  grid.sync();
}

// note: static assert: threadsPerBlock >= maxBlocksInGrid
template <typename TYPE, int threadsPerBlock, int maxBlocksInGrid>
__global__ void inclusive_scan(TYPE *input, TYPE *output, int *idx, int numElements, int *npti) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int chunk_stride = threadsPerBlock * gridDim.x;
  int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
  for (int chunk = 0; chunk < nchunks; chunk++) {
    int gtid = tid + (chunk * chunk_stride);
    TYPE val = 0;
    if (gtid < numElements) {
      val = input[gtid];
    }
    TYPE result = block_scan<TYPE,threadsPerBlock>(val);
    if (gtid < numElements) {
      output[gtid] = result + val;
    }
    merge_blocks<TYPE, threadsPerBlock, maxBlocksInGrid>(tid, blockIdx.x, threadsPerBlock, output, chunk * chunk_stride, numElements);
  }
  if (nchunks > 1) {
    merge_blocks<TYPE, threadsPerBlock, maxBlocksInGrid>(tid, -1, chunk_stride, output, 0, numElements, nchunks);
  }
  if constexpr (std::is_same_v<TYPE, int>) {
    if (idx) {
      for (int chunk = 0; chunk < nchunks; chunk++) {
        int gtid = tid + (chunk * chunk_stride);
        if (gtid < numElements) {
          if (input[gtid] != 0) idx[output[gtid]] = gtid;
        }
      }
      if (threadIdx.x + blockIdx.x == 0 && npti) *npti = output[numElements - 1]; 
    }
  }
}

// Assumed gridDim.x == 1
template <typename TYPE, int threadsPerBlock>
__global__ void inclusive_scan_one_block(TYPE *input, TYPE *output, int *idx, int numElements, int *npti) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int chunk_stride = threadsPerBlock;
    int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
    TYPE addition = 0; 
    for (int chunk = 0; chunk < nchunks; chunk++) {
      int gtid = tid + (chunk * chunk_stride);

      TYPE val = 0;
      if (gtid < numElements) {
        val = input[gtid] + (threadIdx.x ? 0 : addition);
      }
      TYPE result = block_scan<TYPE,threadsPerBlock>(val);
      if (gtid < numElements) {
        output[gtid] = result + val;
        if constexpr (std::is_same_v<TYPE, int>) if (idx && input[gtid]) idx[output[gtid]] = gtid;  
      }
      __syncthreads();
      if (chunk != nchunks - 1)
        addition = output[(chunk + 1) * chunk_stride - 1];
   }
   if constexpr (std::is_same_v<TYPE, int>) if (threadIdx.x == 0 && npti) *npti = output[numElements - 1];
}

