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
- we can reach 2.5..6.2 usec latency for arrays of less than 100K elements
- we can have BW in the diapason: 15..40 Gigatransfers per second for 100K..5M elements (type: int 4 bytes)
- thrust is better than our code starting from 5M elements
- BW smoothly grows after 5M with the thrust version up to ~160 Gigatranfers per second
- we also have the Fortran code with the SUM_PREFIX intrinsic doing the same: it is generally worse 
  besides one diapason: [3M; 10M] where it significantly outpeforms both ours code and thrust code

     bytes   bw GINTps     bw GINTps     bw GINTps         latency usec   latency usec  latency usec
             custom        thrust        SUM_PREFIX        custom         thrust        SUM_PREFIX

         2   0.000768      0.000152      0.000148   --     2.5            12.8          13.2
         3   0.001155      0.000228      0.000222   --     2.5            12.8          13.2
       101   0.038515      0.007607      0.007081   --     2.6            13.0          13.9
       941   0.352519      0.070776      0.044074   --     2.6            13.0          20.9
      1031   0.321611      0.077514      0.046313   --     3.1            13.0          21.7
      2151   0.576403      0.155469      0.068334   --     3.6            13.5          30.7
      3121   0.733768      0.223881      0.078634   --     4.2            13.6          38.8
      4551   0.951454      0.330803      0.138324   --     4.7            13.4          32.1
      6051   1.138263      0.437762      0.153445   --     5.2            13.5          38.5
     10691   1.714553      0.775560      0.334951   --     6.1            13.5          31.2
     20791   3.353768      1.495780      0.649927   --     6.1            13.6          31.2
     30411   4.885383      2.170448      0.922019   --     6.1            13.7          32.2
     48951   7.885799      3.458023      1.434255   --     6.1            13.8          33.3
    104881  16.506503      7.065594      2.206913   --     6.2            14.5          46.4
    153521  21.371754     10.116252      3.802253   --     7.0            14.8          39.4
    204321  27.682051     13.278248      5.504337   --     7.2            15.0          36.2
    247221  32.104622     15.663839      5.606752   --     7.5            15.4          43.1
    271941  23.050723     17.173266      6.863980   --    11.5            15.5          38.7
    481731  37.917339     29.404365     11.087448   --    12.4            16.0          42.4
    705281  40.045990      5.260894     11.634307   --    17.2           130.9          59.2
    938711  40.623503      5.821453     19.537723   --    22.6           157.5          46.9
   1374351  40.438073     11.396278     21.626485   --    33.2           117.8          62.1
   2012161  43.666688     16.530672     37.321956   --    45.0           118.9          52.6
>> 3240571  43.250241     25.653535     57.960075   --    73.2           123.4          54.6
>> 4313181  38.643035     30.122941     64.346025   --   109.0           139.8          65.5
>> 5218931  38.613625     29.709194     60.862339   --   132.0           171.6          83.7
>> 6946391  38.110028     40.099219     54.565516   --   178.0           169.2         124.3
>> 8405121  38.101128     44.179590     49.749233   --   215.4           185.8         165.0
> 10170191  38.571700     55.051423     50.059613   --   257.5           180.4         198.4
> 11187201  38.528005     56.883271     50.588076   --   283.6           192.1         216.0
> 12305921  38.914257     60.971593     51.464610   --   308.8           197.1         233.5
> 13536511  38.786600     61.053247     51.891066   --   340.8           216.5         254.8
> 14890161  38.655855     66.516504     51.643189   --   376.2           218.6         281.6
> 16379171  38.621026     73.484101     52.515872   --   414.2           217.7         304.6
*/

#include "psum.inl.cu"

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
      int numBlocks = (N + threadsPerBlock - 1) / threadsPerBlock;
      if (numBlocks < 7) {
        for (int j = 0; j < ncycles; j++) {
          inclusive_scan_one_block<int,threadsPerBlock><<<1,threadsPerBlock>>>(input, output, numElements);
        }
      //} else if (numBlocks < 4) { 
      //  for (int j = 0; j < ncycles; j++) {
      //    cudaLaunchCooperativeKernel((void *)inclusive_scan_small<int,threadsPerBlock>, numBlocks, threadsPerBlock, args, 0, 0);
      //  }
      } else {
        static int maxNumBlocksPerDevice = 0;
        const int maxBlocksPerGrid = 256; 
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
        void *args[] = {&input, &output, &numElements};
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


