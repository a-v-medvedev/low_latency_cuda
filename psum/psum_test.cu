#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <chrono>

// Scan Algorithm adapted from https://github.com/simveit/effective_scan
// Extended to handle arbitrary size arrays by Alexey V. Medvedev https://github.com/a-v-medvedev/low_latency_cuda

/*
$ nvcc --extended-lambda -arch=native -O3 -o psum_test_custom psum_test.cu
$ nvcc --extended-lambda -arch=native -O3 -DWITH_THRUST_SCAN -o psum_test_thrust psum_test.cu
$ nvfortran -cuda -acc=gpu -O3 -o psum_test_fortran psum_test_fortran.f90 -cudalib=cutensor
$ ./psum_test_custom > custom.log
$ ./psum_test_thrust > thrust.log
$ ./psum_test_fortran > fortran.log
$ paste custom.log thrust.log fortran.log > compare.txt
$ cat compare.txt | sed 's/[ \t]i=[^ ]* / /g;s/[iusec]*=//g;/^[ \t]*1[ \t]/d' > table.txt
$ cat table.txt | awk '{if (NF==4) printf "%10d %10.6f %10.6f %10.6f -- %5.1f %5.1f %5.1f\n", $1, $1 / $2 / 1024, $1 / $3 / 1024, $1 / $4 / 1024, $2, $3, $4 }' > table_pretty.txt

General observations:
- we can reach 2.4..6.6 usec latency for arrays of less than 100K elements
- we can have BW in the diapason: 15..30 Gigatransfers per second for 100K..5M elements (type: int 4 bytes)
- thrust is better in two diapasons:
  - [256K; 512K) elements
  - [5M; Inf) elements
- BW smoothly grows after 5M with the thrust version up to ~160 Gigatranfers per second
- we also have the Fortran code with the SUM_PREFIX intrinsic doing the same: it is generally worse 
  besides one diapason: [1.6M; 10M] where it significantly outpeforms thrust code

     bytes   bw GINTps custom   bw GINTps thrust        usec custom  usec thrust 
         2   0.000824           0.000148         --       2.4          13.2
         3   0.001237           0.000223         --       2.4          13.1
       101   0.041545           0.007421         --       2.4          13.3
       941   0.388248           0.069123         --       2.4          13.3
      1031   0.342893           0.074800         --       2.9          13.5
      2151   0.616333           0.151515         --       3.4          13.9
      3121   0.782022           0.218851         --       3.9          13.9
      4551   1.002444           0.322898         --       4.4          13.8
      6051   1.197645           0.428003         --       4.9          13.8
      6651   1.242918           0.469216         --       5.2          13.8
     10691   1.814717           0.756021         --       5.8          13.8
     20791   3.502451           1.464914         --       5.8          13.9
     30411   5.073154           2.118126         --       5.9          14.0
     48951   8.008663           3.393945         --       6.0          14.1
    104881  15.549241           6.982265         --       6.6          14.7
    153521  16.992276          10.085627         --       8.8          14.9
    204321  21.181765          13.040470         --       9.4          15.3
    247221  24.101703          15.453291         --      10.0          15.6
>   271941  13.132597          17.138908         --      20.2          15.5
>   481731  23.155015          29.107810         --      20.3          16.2
    705281  24.938481           5.348110         --      27.6         128.8
    938711  26.383180           7.105563         --      34.7         129.0
   1374351  26.966840           9.423814         --      49.8         142.4
   2012161  30.366265          14.550174         --      64.7         135.1
   3240571  30.983162          23.248752         --     102.1         136.1
   4313181  29.552311          25.610086         --     142.5         164.5
>  5218931  29.550718          32.524648         --     172.5         156.7
>  6946391  27.901061          38.785506         --     243.1         174.9
>  8405121  27.583849          43.926608         --     297.6         186.9
> 10170191  28.047294          52.491027         --     354.1         189.2
> 11187201  27.934749          56.213023         --     391.1         194.3
> 12305921  28.267826          54.899502         --     425.1         218.9
> 13536511  28.242034          60.558198         --     468.1         218.3
> 14890161  28.129868          66.997663         --     516.9         217.0
> 16379171  28.277706          70.142449         --     565.6         228.0
*/

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

//#define WITH_THRUST_SCAN
//#define WITH_CHECK
#define INITIAL_NCYCLES 10000

auto make_zip_begin(const thrust::counting_iterator<int> &counter, const thrust::device_vector<int> &v) {
  return thrust::make_zip_iterator(thrust::make_tuple(counter, v.begin()));
}

auto make_zip_end(const thrust::counting_iterator<int> &counter, const thrust::device_vector<int> &v) {
  return thrust::make_zip_iterator(thrust::make_tuple(counter + v.size(), v.end()));
}

#define PTRCAST(ptr, type) thrust::device_ptr<type>(ptr)

int ncycles = INITIAL_NCYCLES;
int test(const unsigned int N) {
  if (N < 2)
      return 0;
  thrust::device_vector<int> x(N), y(N);
  thrust::sequence(x.begin(), x.end(), 1);
  struct zipit {
    const thrust::counting_iterator<int> &_counter;
    const thrust::device_vector<int> &_v;
    zipit(const thrust::counting_iterator<int> &counter, const thrust::device_vector<int> &v) : _v(v), _counter(counter) {}
    auto begin() { return thrust::make_zip_iterator(thrust::make_tuple(_counter, _v.begin())); }
    auto end() { return thrust::make_zip_iterator(thrust::make_tuple(_counter + _v.size(), _v.end())); }
  };
  thrust::counting_iterator<int> counter(0);
  thrust::transform(zipit(counter, x).begin(), zipit(counter, x).end(),
    x.begin(),
    [] __host__ __device__ (thrust::tuple<int,int> t) {
        auto [idx, val] = t;
        return (idx % 2 == 1) ? -val : val;
  });
  thrust::fill(y.begin(), y.end(), 0);
  int *input = thrust::raw_pointer_cast(thrust::device_ptr<int>(x.data()));
  int *output = thrust::raw_pointer_cast(thrust::device_ptr<int>(y.data()));
  size_t numElements = x.size();
 
  using namespace std::chrono;
  auto t1 = high_resolution_clock::now();  
#if !defined WITH_THRUST_SCAN
  {
      const int threadsPerBlock = 1024; 
      static int maxNumBlocksPerDevice = 0;
      const int maxBlocksPerGrid = 256; 
      int numBlocks = (N + threadsPerBlock - 1) / threadsPerBlock;
      void *args[] = {&input, &output, &numElements};
      if (numBlocks < 7) {
        for (int j = 0; j < ncycles; j++) {
          inclusive_scan_one_block<int,threadsPerBlock><<<1,threadsPerBlock>>>(input, output, numElements);
        }
      //} else if (numBlocks < 4) { 
      //  for (int j = 0; j < ncycles; j++) {
      //    cudaLaunchCooperativeKernel((void *)inclusive_scan_small<int,threadsPerBlock>, numBlocks, threadsPerBlock, args, 0, 0);
      //  }
      } else {
        if (!maxNumBlocksPerDevice) {
          cudaDeviceProp deviceProp;
          cudaGetDeviceProperties(&deviceProp, 0);
          int numBlocksPerSm = 0;
          cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, inclusive_scan<int,threadsPerBlock,maxBlocksPerGrid>, threadsPerBlock, 0);
          maxNumBlocksPerDevice = deviceProp.multiProcessorCount * numBlocksPerSm;
        }
        int maxNumBlocks = std::min(maxNumBlocksPerDevice, maxBlocksPerGrid);
        // we limit the number of blocks to the maximum available due to cooperative kernels limitations.
        // the kernel will automatically detect that we have more elements to process than threads in grid, 
        // and will do sequential execution in chunks
        numBlocks = std::min(numBlocks, maxNumBlocks); 
        for (int j = 0; j < ncycles; j++) {
          cudaLaunchCooperativeKernel((void *)inclusive_scan<int,threadsPerBlock,maxBlocksPerGrid>, numBlocks, threadsPerBlock, args, 0, 0);
        } 
      } 
  }
#else
  {
      thrust::plus<int> binary_op; 
      for (int j = 0; j < ncycles; j++) {
        thrust::inclusive_scan(PTRCAST(input, int), PTRCAST(input, int) + numElements, PTRCAST(output, int), binary_op); 
      }
  }
#endif
  cudaDeviceSynchronize();
  auto t2 = high_resolution_clock::now();  
  auto usec = duration_cast<microseconds>(t2 - t1).count();
#if defined WITH_CHECK
  thrust::host_vector<int> xh = x;
  thrust::host_vector<int> yh = y;
  int cpu_psum = xh[0];
  for (size_t i = 1; i < yh.size(); i++) {
    cpu_psum += xh[i];
    if (yh[i] != cpu_psum) {
      std::cout << ">> 1: i=" << i << " " << yh[i] << " " << cpu_psum << std::endl;
      return 1;
    }
  }
  //std::cout << ">> " << yh.back() << " " << cpu_psum << std::endl;
#endif
  return -(int)usec;
}

int main(int argc, char **argv)
{
  //for (unsigned i = 20*1024*1024 + 1; i <= 1024*1024*1024;) {
  for (unsigned i = 1; i <= 16*1024*1024;) {
  //for (unsigned i = 1; i <= 32*1024;) {
    if (i > 16*1024) ncycles = INITIAL_NCYCLES / 10;
    if (i > 1024*1024) ncycles = INITIAL_NCYCLES / 100;
    int result = test(i);
    int usec = 0;
    if (result < 0) {
      usec = -result;
      result = 0;
    }
    switch (result) {
      case 0: break;
      case 1:
      case 2:
      case 3: std::cout << "ERROR (" << result << "): " << i << std::endl; break;
    }
    if (!result && ((i < 10) || (i % 10 == 1))) {
      std::cout << "i=" << i << " usec=" << (float)usec/(float)ncycles << std::endl;
    }
    if (result)
      break;

    int incr = (i/100)*10;
    if (!incr && i < 11)
      incr = 1; 
    i += (incr ? incr : 10);
  }
  return 0;
}


