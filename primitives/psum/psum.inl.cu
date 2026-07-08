// Scan Algorithm based on block_scan adapted from https://github.com/simveit/effective_scan
// Extended to handle arbitrary size arrays by Alexey V. Medvedev https://github.com/a-v-medvedev/low_latency_cuda
//
// inclusive_scan / inclusive_scan_one_block are meant to be normal inclusive scan with
// the "+" basic operation (i.e. prefix sum operation) focused on achieving low latency
// for input array size less than 1M elements

#include "block_scan.inl.cu"

// Merge as much separate sequential blocks as we have CUDA blocks in the 
// execution grid. We can't merge more in one shot as it is a limitation
// of cooperative groups API we use for barrier syncing.
// 
// We are also aware that we may go through the whole array in sequential
// chunks, so we mind the offset and the max value of a previous chunk
// if the offset is not 0.
template <typename TYPE, int threadsPerBlock>
__device__ void merge_blocks(int tid, TYPE *output, int offset, int numElements) {
  __shared__ TYPE sdata[1];
  namespace cg = cooperative_groups;
  cg::grid_group grid = cg::this_grid();
  grid.sync();
  int gtid = offset + tid;
  TYPE addition = 0;
  if (blockIdx.x) {
    TYPE value = 0;
    int idx_of_addition = (threadIdx.x + 1) * threadsPerBlock - 1;
    if (offset + idx_of_addition < numElements) {
      value = output[offset + idx_of_addition];
    }
    addition = block_scan<TYPE,threadsPerBlock>(value);
    if (threadIdx.x == blockIdx.x) {
      sdata[0] = addition;
    }
    __syncthreads();
    addition = sdata[0];
  }
  grid.sync();
  if (gtid < numElements) {
    TYPE addition_from_prev_chunk = 0;
    if (offset) addition_from_prev_chunk = output[offset - 1];
    output[gtid] += addition + addition_from_prev_chunk;
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
    TYPE val = 0;
    if (gtid < numElements) {
      val = input[gtid];
    }
    TYPE result = block_scan<TYPE,threadsPerBlock>(val) + val;
    if (gtid < numElements) {
      output[gtid] = result;
    }
    merge_blocks<TYPE, threadsPerBlock>(tid, output, chunk * chunk_stride, numElements);
  }
}

// Assumed gridDim.x == 1
template <typename TYPE, int threadsPerBlock>
__global__ void inclusive_scan_one_block(TYPE *input, TYPE *output, int numElements) {
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
      }
      __syncthreads();
      if (chunk != nchunks - 1)
        addition = output[(chunk + 1) * chunk_stride - 1];
   }
}

