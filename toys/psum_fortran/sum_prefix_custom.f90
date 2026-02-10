module sum_prefix_custom
   implicit none

#define DEVICE_ATTR , device

   interface sum_prefix_custom
      module procedure scan_int, scan_float, scan_double 
   end interface

   interface
   subroutine scan_int_wrapper(input, output, N, stream) bind(C, name="scan_int_wrapper")
      use iso_c_binding
      use openacc
      integer(c_int) DEVICE_ATTR :: input(*)
      integer(c_int) DEVICE_ATTR :: output(*)
      integer(c_int), value :: N
      integer(acc_handle_kind), value :: stream
   end subroutine

   subroutine scan_float_wrapper(input, output, N, stream) bind(C, name="scan_float_wrapper")
      use iso_c_binding
      use openacc
      real(c_float) DEVICE_ATTR :: input(*)
      real(c_float) DEVICE_ATTR :: output(*)
      integer(c_int) :: N
      integer(acc_handle_kind), value :: stream
   end subroutine

   subroutine scan_double_wrapper(input, output, N, stream) bind(C, name="scan_double_wrapper")
      use iso_c_binding
      use openacc
      real(c_double) DEVICE_ATTR :: input(*)
      real(c_double) DEVICE_ATTR :: output(*)
      integer(c_int) :: N
      integer(acc_handle_kind), value :: stream
   end subroutine
   end interface

contains

   subroutine scan_int(input, output)
      USE openacc
#ifdef WITH_CUTENSOREX
      use cutensorex, only: SUM_PREFIX
#endif
      integer, intent(in) DEVICE_ATTR :: input(:)
      integer, intent(out) DEVICE_ATTR :: output(:)
      INTEGER(acc_handle_kind) :: stream
      integer(c_int) :: N
#ifdef WITH_CUTENSOREX
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
#endif
        N = size(input)
        stream = acc_get_cuda_stream(1)   
        call scan_int_wrapper(input, output, N, stream)
#ifdef WITH_CUTENSOREX
      endif
#endif
   end subroutine

   subroutine scan_float(input, output)
      use iso_fortran_env, only : real32, real64
      USE openacc
#ifdef WITH_CUTENSOREX
      use cutensorex, only: SUM_PREFIX
#endif
      real(real32), intent(in) DEVICE_ATTR :: input(:)
      real(real32), intent(out) DEVICE_ATTR :: output(:)
      INTEGER(acc_handle_kind) :: stream
#ifdef WITH_CUTENSOREX
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
#endif
        stream = acc_get_cuda_stream(1)
        call scan_float_wrapper(input, output, size(input), stream)
#ifdef WITH_CUTENSOREX
      endif
#endif
   end subroutine

   subroutine scan_double(input, output)
      use iso_fortran_env, only : real32, real64
      USE openacc
#ifdef WITH_CUTENSOREX
      use cutensorex, only: SUM_PREFIX
#endif
      real(real64), intent(in) DEVICE_ATTR :: input(:)
      real(real64), intent(out) DEVICE_ATTR :: output(:)
      INTEGER(acc_handle_kind) :: stream
#ifdef WITH_CUTENSOREX
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
#endif
        stream = acc_get_cuda_stream(1)
        call scan_double_wrapper(input, output, size(input), stream)
#ifdef WITH_CUTENSOREX
      endif
#endif
   end subroutine

end module sum_prefix_custom
