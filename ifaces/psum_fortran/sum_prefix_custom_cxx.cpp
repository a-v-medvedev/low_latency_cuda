#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <cooperative_groups.h>
#include <cuda_runtime.h>
#include <semaphore>
#include <atomic>

#include "psum.inl.cu"

#define PTRCAST(ptr, type) thrust::device_ptr<type>(ptr)

static int maxNumBlocksPerDevice = 0;

// NOTE: async GPU kernels execution, expected to sync the stream outside
template <typename TYPE>
void scan_wrapper_cxx(TYPE *input, TYPE *output, int numElements, void *stream_void) { 
  cudaStream_t stream = reinterpret_cast<cudaStream_t>(stream_void); 
  if (numElements > 5*1024*1024) { 
    thrust::plus<TYPE> binary_op; 
    thrust::inclusive_scan(thrust::cuda::par.on(stream), PTRCAST(input, TYPE), PTRCAST(input, TYPE) + 
        numElements, PTRCAST(output, TYPE), binary_op); 
  } else { 
    const int threadsPerBlock = 1024; 
    const int maxBlocksInGrid = 256; 
    int numBlocks = (numElements + threadsPerBlock - 1) / threadsPerBlock; 
    if (numBlocks < 7) { 
      inclusive_scan_one_block<TYPE,threadsPerBlock><<<1,threadsPerBlock,0,stream>>>(input, output, numElements); 
    } else { 
        if (!maxNumBlocksPerDevice) { 
          cudaDeviceProp deviceProp; 
          cudaGetDeviceProperties(&deviceProp, 0); 
          int numBlocksPerSm = 0; 
          cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, inclusive_scan<TYPE,threadsPerBlock,maxBlocksInGrid>, threadsPerBlock, 0); 
          maxNumBlocksPerDevice = deviceProp.multiProcessorCount * numBlocksPerSm; 
        } 
        int maxNumBlocks = std::min(maxNumBlocksPerDevice, maxBlocksInGrid); 
        numBlocks = std::min(numBlocks, maxNumBlocks); 
        void *args[] = {&input, &output, &numElements}; 
        cudaLaunchCooperativeKernel((void *)inclusive_scan<TYPE,threadsPerBlock,maxBlocksInGrid>, numBlocks, threadsPerBlock, args, 0, stream); 
    } 
  } 
}

#define DECLARE_SCAN_WRAPPER(TYPE) \
void scan_##TYPE##_wrapper(TYPE *input, TYPE *output, int numElements, void *stream_void) { \
  scan_wrapper_cxx<TYPE>(input, output, numElements, stream_void); \
}

extern "C" {
DECLARE_SCAN_WRAPPER(int)
DECLARE_SCAN_WRAPPER(float)
DECLARE_SCAN_WRAPPER(double)
}

