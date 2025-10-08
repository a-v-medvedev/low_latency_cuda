#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <cooperative_groups.h>
#include <cuda_runtime.h>

#include "../../psum/psum.inl.cu"

#define PTRCAST(ptr, type) thrust::device_ptr<type>(ptr)

#define DEFINE_SCAN_WRAPPER(TYPE) \
void scan_##TYPE##_wrapper(TYPE *input, TYPE *output, int numElements) { \
  if (numElements > 5*1024*1024 || (numElements > 256 * 1024 && numElements < 512 * 1024)) { \
    thrust::plus<TYPE> binary_op; \
    thrust::inclusive_scan(PTRCAST(input, TYPE), PTRCAST(input, TYPE) + \
        numElements, PTRCAST(output, TYPE), binary_op); \
  } else { \
    const int threadsPerBlock = 1024; \
    const int maxBlocksInGrid = 256; \
    int numBlocks = (numElements + threadsPerBlock - 1)/threadsPerBlock; \
    void *args[] = {&input, &output, &numElements}; \
    if (numBlocks < 7) { \
      inclusive_scan_one_block<TYPE,threadsPerBlock><<<1,threadsPerBlock>>>(input, output, numElements); \
    } else { \
        if (!maxNumBlocksPerDevice) { \
          cudaDeviceProp deviceProp; \
          cudaGetDeviceProperties(&deviceProp, 0); \
          int numBlocksPerSm = 0; \
          cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, inclusive_scan<int,threadsPerBlock,maxBlocksInGrid>, threadsPerBlock, 0); \
          maxNumBlocksPerDevice = deviceProp.multiProcessorCount * numBlocksPerSm; \
        } \
        int maxNumBlocks = std::min(maxNumBlocksPerDevice, maxBlocksInGrid); \
        numBlocks = std::min(numBlocks, maxNumBlocks); \
        cudaLaunchCooperativeKernel((void *)inclusive_scan<TYPE,threadsPerBlock,maxBlocksInGrid>, numBlocks, threadsPerBlock, args, 0, 0); \
    } \
  } \
}

static int maxNumBlocksPerDevice = 0;

extern "C" {
DEFINE_SCAN_WRAPPER(int)
DEFINE_SCAN_WRAPPER(float)
DEFINE_SCAN_WRAPPER(double)
}

