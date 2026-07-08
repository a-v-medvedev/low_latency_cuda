# Collection of small CUDA primitives that may show better latency than standard library functions

## Custom inclusive scan

Inspired by the work done here: `https://github.com/simveit/effective_scan` by Simon Veitner. Described in the blog post: https://veitner.bearblog.dev/making-prefix-sum-really-fast

Trying to beat `thrust::inclusive_scan` on small inputs. Currecntly, it is really faster on test arrays consisting of 4-byte integers with up to 5M elements. Another competing solution is the `SUM_PREFIX` subroutine from `cutensor` Fortran library (is it ortran-only??). In case of `SUM_PREFIX`, it is faster than our solution for the diapason of [3M; 10M] elements.

But the real target for us is arrays with less than 1M elements, and here we have 2x to 5x better latency value compared both to `thrust::inclusive_scan` and `SUM_PREFIX` of `cutensor`.

## Custom packloc

Another primitive we are addressing is the `PACKLOC` intrinsic specific to nvfortran. We tried to have a custom implementation of the same idea, also focusing on input arrays with less than 1M elements. Indeed we can have up to 2.5x less latency compared to NVIDIA's PACKLOC from cutensor. The custom packloc code also uses the `block_scan` of the `https://github.com/simveit/effective_scan` mentioned above as a building block.

## Toy code: `sc_icedyn`

The Fortran code that uses both Custom inclusive scan and Custom packloc primitives mentioned above that mocks the specific array compaction pattern as found in the NEMO v.4 and NEMO v.5 Sea-Ice (SI3) module (see: https://nemo-ocean.eu). We try to compare different GPU-specific solutions with the original CPU-centric serial code using the toy mocking-code tests. As for now, the Custom packloc based solution is the best.

