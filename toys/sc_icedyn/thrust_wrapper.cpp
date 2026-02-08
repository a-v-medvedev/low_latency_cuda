#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <cooperative_groups.h>
#include <cuda_runtime.h>

#include "../../psum/psum.inl.cu"

#define PTRCAST(ptr, type) thrust::device_ptr<type>(ptr)

#define DEFINE_SCAN_WRAPPER(TYPE) \
int scan_##TYPE##_wrapper(TYPE *input, TYPE *output, int *idx, int numElements, void *stream_void) { \
  cudaStream_t stream = reinterpret_cast<cudaStream_t>(stream_void); \
  int retval = 0; \
  if (numElements > 5*1024*1024 || (numElements > 256 * 1024 && numElements < 512 * 1024)) { \
    thrust::plus<TYPE> binary_op; \
    thrust::inclusive_scan(thrust::cuda::par.on(stream), PTRCAST(input, TYPE), PTRCAST(input, TYPE) + \
        numElements, PTRCAST(output, TYPE), binary_op); \
  } else { \
    retval = 1; \
    const int threadsPerBlock = 1024; \
    const int maxBlocksInGrid = 256; \
    int numBlocks = (numElements + threadsPerBlock - 1)/threadsPerBlock; \
    void *args[] = {&input, &output, &idx, &numElements}; \
    if (numBlocks < 7) { \
      if (idx == 0) { \
        inclusive_scan_one_block<TYPE,threadsPerBlock,false><<<1,threadsPerBlock,0,stream>>>(input, output, 0, numElements); \
      } else { \
        inclusive_scan_one_block<TYPE,threadsPerBlock,true><<<1,threadsPerBlock,0,stream>>>(input, output, idx, numElements); \
      } \
    } else { \
        if (!maxNumBlocksPerDevice) { \
          cudaDeviceProp deviceProp; \
          cudaGetDeviceProperties(&deviceProp, 0); \
          int numBlocksPerSm = 0; \
          cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, inclusive_scan<int,threadsPerBlock,maxBlocksInGrid,false>, threadsPerBlock, 0); \
          maxNumBlocksPerDevice = deviceProp.multiProcessorCount * numBlocksPerSm; \
        } \
        int maxNumBlocks = std::min(maxNumBlocksPerDevice, maxBlocksInGrid); \
        numBlocks = std::min(numBlocks, maxNumBlocks); \
        if (idx == 0) { \
          cudaLaunchCooperativeKernel((void *)inclusive_scan<TYPE,threadsPerBlock,maxBlocksInGrid,false>, numBlocks, threadsPerBlock, args, 0, stream); \
        } else { \
          cudaLaunchCooperativeKernel((void *)inclusive_scan<TYPE,threadsPerBlock,maxBlocksInGrid,true>, numBlocks, threadsPerBlock, args, 0, stream); \
        } \
    } \
  } \
  return retval; \
}

static int maxNumBlocksPerDevice = 0;

extern "C" {
DEFINE_SCAN_WRAPPER(int)
//DEFINE_SCAN_WRAPPER(float)
//DEFINE_SCAN_WRAPPER(double)
}

