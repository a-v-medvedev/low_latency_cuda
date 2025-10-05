# Collection of small CUDA primitives that may show better latency than standard library functions

## Custom inclusive scan

Inspired by the work done here: `https://github.com/simveit/effective_scan` by Simon Veitner. Described in the blog post: https://veitner.bearblog.dev/making-prefix-sum-really-fast

Trying to beat `thrust::inclusive_scan` on small inputs.


