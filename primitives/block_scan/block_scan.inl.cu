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

