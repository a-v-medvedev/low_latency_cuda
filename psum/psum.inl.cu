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
template <typename TYPEIN, typename TYPEOUT, int threadsPerBlock, int maxBlocksInGrid>
__global__ void inclusive_scan(TYPEIN *input, TYPEOUT *output, int numElements) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int chunk_stride = threadsPerBlock * gridDim.x;
  int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
  for (int chunk = 0; chunk < nchunks; chunk++) {
    int gtid = tid + (chunk * chunk_stride);
    TYPEIN val = 0;
    if (gtid < numElements) {
      val = input[gtid];
    }
    TYPEOUT result = block_scan<TYPEOUT,threadsPerBlock>((TYPEOUT)val);
    if (gtid < numElements) {
      output[gtid] = result + val;
    }
    merge_blocks<TYPEOUT, threadsPerBlock>(tid, output, chunk * chunk_stride, numElements);
  }
}

// Assumed gridDim.x == 1
template <typename TYPEIN, typename TYPEOUT, int threadsPerBlock>
__global__ void inclusive_scan_one_block(TYPEIN *input, TYPEOUT *output, int numElements) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int chunk_stride = threadsPerBlock;
    int nchunks = (numElements + chunk_stride - 1) / chunk_stride;
    TYPEOUT addition = 0; 
    for (int chunk = 0; chunk < nchunks; chunk++) {
      int gtid = tid + (chunk * chunk_stride);

      TYPEOUT val = 0;
      if (gtid < numElements) {
        val = (TYPEOUT)input[gtid] + (threadIdx.x ? 0 : addition);
      }
      TYPEOUT result = block_scan<TYPEOUT,threadsPerBlock>(val);
      if (gtid < numElements) {
        output[gtid] = result + val;
      }
      __syncthreads();
      if (chunk != nchunks - 1)
        addition = output[(chunk + 1) * chunk_stride - 1];
   }
}

//--- packloc:

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
    TYPE input_value_for_addition = latest_inputs_per_block[threadIdx.x];
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
      if (gtid < numElements && input[gtid] != 0) idx[result] = gtid;
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
    TYPEOUT addition = 0; 
    for (int chunk = 0; chunk < nchunks; chunk++) {
      int gtid = tid + (chunk * chunk_stride);
      TYPEOUT val = 0;
      if (gtid < numElements) {
        val = (TYPEOUT)(input[gtid]?1:0) + (threadIdx.x ? 0 : addition);
      }
      TYPEOUT result = block_scan<TYPEOUT,threadsPerBlock>(val);
      if (gtid < numElements) {
        if constexpr (std::is_same_v<TYPEOUT, int>) if (input[gtid]) idx[result + val] = gtid;
      }
      if (tid == chunk_stride - 1) sdata[0] = result + val;
      if (chunk == nchunks - 1 && gtid == numElements - 1) sdata[1] = result + val;
      __syncthreads();
   }
   if constexpr (std::is_same_v<TYPEOUT, int>) if (threadIdx.x == 0) *n = sdata[1];
}

