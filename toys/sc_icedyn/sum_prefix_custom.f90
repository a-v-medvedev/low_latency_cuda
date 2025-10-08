module sum_prefix_custom
   implicit none

#define DEVICE_ATTR device

   interface sum_prefix_custom
   subroutine scan_int_wrapper(input, output, N) bind(C, name="scan_int_wrapper")
      use iso_c_binding
      use openacc
      integer(c_int), DEVICE_ATTR :: input(*)
      integer(c_int), DEVICE_ATTR :: output(*)
      integer(c_int) :: N
   end subroutine

   subroutine scan_float_wrapper(input, output, N) bind(C, name="scan_float_wrapper")
      use iso_c_binding
      use openacc
      real(c_float), DEVICE_ATTR :: input(*)
      real(c_float), DEVICE_ATTR :: output(*)
      integer(c_int) :: N
   end subroutine

   subroutine scan_double_wrapper(input, output, N) bind(C, name="scan_double_wrapper")
      use iso_c_binding
      use openacc
      real(c_double), DEVICE_ATTR :: input(*)
      real(c_double), DEVICE_ATTR :: output(*)
      integer(c_int) :: N
   end subroutine
   end interface

contains

   subroutine scan_int(input, output)
      use iso_c_binding
      use cutensorex, only: SUM_PREFIX
      integer(c_int), intent(in), DEVICE_ATTR:: input(:)
      integer(c_int), intent(out), DEVICE_ATTR:: output(:)
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
        call scan_int_wrapper(input, output, size(input))
      endif
   end subroutine

   subroutine scan_float(input, output)
      use iso_c_binding
      use cutensorex, only: SUM_PREFIX
      real(c_float), intent(in), DEVICE_ATTR:: input(:)
      real(c_float), intent(out), DEVICE_ATTR:: output(:)
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
        call scan_float_wrapper(input, output, size(input))
      endif
   end subroutine

   subroutine scan_double(input, output)
      use iso_c_binding
      use cutensorex, only: SUM_PREFIX
      real(c_double), intent(in), DEVICE_ATTR:: input(:)
      real(c_double), intent(out), DEVICE_ATTR:: output(:)
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
        call scan_double_wrapper(input, output, size(input))
      endif
   end subroutine


end module sum_prefix_custom
