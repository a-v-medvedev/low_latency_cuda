// Packloc algorithm based on block_scan (adapted from https://github.com/simveit/effective_scan)
// by Alexey V. Medvedev https://github.com/a-v-medvedev/low_latency_cuda
//
// PACKLOC is a Fortran language extension made by NVIDIA for their nvfrotran compiler.
// You use it like:
//   idx = packloc(mask, count=output_count)
// to have indeces of all true elements in mask array collected in idx, and count is a scalar integer
// for the total number of true elements in mask array.
//
// Here we implement the primitive similar to packloc focusing on low latency for
// input arrays of less than 1M elements

#include "block_scan.inl.cu"

template <typename TYPE, int threadsPerBlock, int maxBlocksInGrid>
__device__ TYPE merge_blocks_nooutput(int tid, int chunk, bool lastinchunk, TYPE input_value) {
  static __device__ TYPE latest_inputs_per_block[maxBlocksInGrid];
  static __device__ TYPE latest_output_for_previous_chunk[1];
  __shared__ TYPE sdata[1];
  namespace cg = cooperative_groups;
  cg::grid_group grid = cg::this_grid();
  if (threadIdx.x == threadsPerBlock - 1)
    latest_inputs_per_block[blockIdx.x] = input_value;
  grid.sync();

  TYPE addition = 0;
  if (blockIdx.x) {
    TYPE input_value_for_addition = 0; 
    if (threadIdx.x < maxBlocksInGrid) {
      input_value_for_addition = latest_inputs_per_block[threadIdx.x];
    }
    addition = block_scan<TYPE,threadsPerBlock>(input_value_for_addition);
    if (threadIdx.x == blockIdx.x) {
      sdata[0] = addition;
    }
    __syncthreads();
    addition = sdata[0];
  }
  TYPE addition_from_prev_chunk = 0;
  if (chunk) addition_from_prev_chunk = latest_output_for_previous_chunk[0];
  input_value += addition + addition_from_prev_chunk;
  grid.sync();
  if (lastinchunk) latest_output_for_previous_chunk[0] = input_value;    
  return input_value;
}

// NOTE: we assume Fortran style for indexing, that means the index of the first element is "1", not "0"
// therefore the lowest index value for idx[] is 1.
// note: static assert: threadsPerBlock >= maxBlocksInGrid
template <typename TYPEIN, typename TYPEOUT, int threadsPerBlock, int maxBlocksInGrid>
__global__ void packloc(TYPEIN *input, TYPEOUT *idx, int *n, int numElements) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int chunk_stride = threadsPerBlock * gridDim.x;
  int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
  TYPEOUT result = 0;
  for (int chunk = 0; chunk < nchunks; chunk++) {
    int gtid = tid + (chunk * chunk_stride);
    TYPEOUT val = 0;
    if (gtid < numElements) {
      val = (TYPEOUT)(input[gtid]?1:0);
    }
    result = block_scan<TYPEOUT,threadsPerBlock>(val) + val;
    result = merge_blocks_nooutput<TYPEOUT, threadsPerBlock,maxBlocksInGrid>(tid, chunk, tid == chunk_stride - 1, result);
    if constexpr (std::is_same_v<TYPEOUT, int>) {
      if (gtid < numElements && input[gtid] != 0) idx[result - 1] = gtid + 1;
    }
  }
  if constexpr (std::is_same_v<TYPEOUT, int>) {
    if (tid + ((nchunks - 1) * chunk_stride) == numElements - 1) *n = result;
  }
}

// Assumed gridDim.x == 1
template <typename TYPEIN, typename TYPEOUT, int threadsPerBlock>
__global__ void packloc_one_block(TYPEIN *input, TYPEOUT *idx, int *n, int numElements) {
    __shared__ TYPEOUT sdata[2];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int chunk_stride = threadsPerBlock;
    int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
    sdata[0] = 0;
    for (int chunk = 0; chunk < nchunks; chunk++) {
      int gtid = tid + (chunk * chunk_stride);
      TYPEOUT val = 0;
      if (gtid < numElements) {
        val = (TYPEOUT)(input[gtid]?1:0) + (threadIdx.x ? 0 : sdata[0]);
      }
      TYPEOUT result = block_scan<TYPEOUT,threadsPerBlock>(val) + val;
      if (gtid < numElements) {
        if constexpr (std::is_same_v<TYPEOUT, int>) if (input[gtid]) idx[result - 1] = gtid + 1;
      }
      if (tid == chunk_stride - 1) sdata[0] = result;
      if (chunk == nchunks - 1 && gtid == numElements - 1) sdata[1] = result;
      __syncthreads();
   }
   if constexpr (std::is_same_v<TYPEOUT, int>) if (threadIdx.x == 0) *n = sdata[1];
}

