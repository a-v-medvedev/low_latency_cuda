# Collection of small CUDA primitives that may provide lower latency than standard library functions

## Custom inclusive scan

Inspired by the work done here: `https://github.com/simveit/effective_scan` by Simon Veitner. Described in the blog post: https://veitner.bearblog.dev/making-prefix-sum-really-fast

Trying to beat `thrust::inclusive_scan` on small inputs. Currently, the custom code is significantly faster on test arrays consisting of 4-byte integers with up to 5M elements. 

Another competing solution is the `SUM_PREFIX` subroutine from the `cuTENSOR` Fortran library (is it Fortran-only??). In the case of `SUM_PREFIX`, it is faster than our solution for array sizes in the range [3M; 10M] elements.

But the arrays with fewer than 1M elements are our real target, and for such cases we have 2x to 5x lower latency compared to both `thrust::inclusive_scan` and `SUM_PREFIX` of `cuTENSOR`.

## Custom packloc

Another primitive we are addressing is the `PACKLOC` intrinsic specific to NVFORTRAN. We implemented a custom version of the same idea, also focusing on input arrays with fewer than 1M elements. Indeed, we observe up to 2.5x lower latency compared to NVIDIA's PACKLOC from cuTENSOR. The custom packloc code also uses the `block_scan` implementation from the `effective_scan` project mentioned above.

The additional features of the Custom packloc are:
- it can run in a specific CUDA stream, and this helps to eliminate the stream synchronization overhead in some cases
- it can take the array size as an argument, in that case, only the initial part of the array is processed (which is important for our practical application in the NEMO Ocean Modeling code)
- this code can potentially be ported to the AMD ecosystem using HIP.

## Toy code: `sc_icedyn`

This Fortran code uses the Custom inclusive scan and Custom packloc primitives described above to mock the array compaction pattern found in the NEMO v.4 and NEMO v.5 Sea-Ice (SI3) module (see https://nemo-ocean.eu). We try to compare different GPU-specific solutions with the original CPU-centric serial code using the toy mock test cases. Currently, the Custom packloc-based solution delivers the best performance.

