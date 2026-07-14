#include <cooperative_groups.h>
#include <cuda_runtime.h>
#include <semaphore>
#include <unordered_map>

#include "packloc.inl.cu"

static int maxNumBlocksPerDevice = 0;
static int *nptis = nullptr;
static int *nptis_dev = nullptr;
static unsigned int cuda_streams_counter = 0;
static std::unordered_map<cudaStream_t, int> stream_to_id;
static std::binary_semaphore sem{1};

// NOTE: we do sync the CUDA stream inside in order to have the correct npti value in the CPU code flow
template <typename TYPEIN, typename TYPEOUT>
int packloc_wrapper_cxx(TYPEIN *input, TYPEOUT *idx, int numElements, void *stream_void) { 
  cudaStream_t stream = reinterpret_cast<cudaStream_t>(stream_void); 
  const int threadsPerBlock = 1024; 
  const int maxBlocksInGrid = 256; 
  while (nptis == nullptr) {
    sem.acquire();
    if (nptis == nullptr) {
      cudaHostAlloc(&nptis, sizeof(int) * 32, cudaHostAllocMapped);
      cudaHostGetDevicePointer(&nptis_dev, nptis, 0);
    }
    sem.release();
  }
  int my_npti_idx = 0;
  auto it = stream_to_id.find(stream);
  if (it != stream_to_id.end()) {
    my_npti_idx = it->second;
  } else {
    sem.acquire();
    my_npti_idx = cuda_streams_counter++;
    stream_to_id.emplace(stream, my_npti_idx);
    sem.release();
  }
  int *my_npti_dev = &nptis_dev[my_npti_idx];
  int numBlocks = (numElements + threadsPerBlock - 1) / threadsPerBlock;
  if (numBlocks < 7) { 
    packloc_one_block<TYPEIN,TYPEOUT,threadsPerBlock><<<1,threadsPerBlock,0,stream>>>(input, idx, my_npti_dev, numElements); 
  } else { 
    if (!maxNumBlocksPerDevice) { 
      cudaDeviceProp deviceProp; 
      cudaGetDeviceProperties(&deviceProp, 0); 
      int numBlocksPerSm = 0; 
      cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, packloc<TYPEIN,TYPEOUT,threadsPerBlock,maxBlocksInGrid>, threadsPerBlock, 0); 
      maxNumBlocksPerDevice = deviceProp.multiProcessorCount * numBlocksPerSm; 
    } 
    int maxNumBlocks = std::min(maxNumBlocksPerDevice, maxBlocksInGrid); 
    numBlocks = std::min(numBlocks, maxNumBlocks); 
    void *args[] = {&input, &idx, &my_npti_dev, &numElements}; 
    cudaLaunchCooperativeKernel((void *)packloc<TYPEIN,TYPEOUT,threadsPerBlock,maxBlocksInGrid>, numBlocks, threadsPerBlock, args, 0, stream); 
  }
  cudaStreamSynchronize(stream);
  int npti = nptis[my_npti_idx]; 
  return npti;
}

#define DECLARE_PACKLOC_WRAPPER(TYPEIN,TYPEOUT) \
int packloc_##TYPEIN##_##TYPEOUT##_wrapper(TYPEIN *input, TYPEOUT *idx, int numElements, void *stream_void) { \
  return packloc_wrapper_cxx<TYPEIN,TYPEOUT>(input, idx, numElements, stream_void); \
}

extern "C" {
DECLARE_PACKLOC_WRAPPER(char,int)
DECLARE_PACKLOC_WRAPPER(int,int)
}

