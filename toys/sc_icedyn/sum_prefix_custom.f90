module sum_prefix_custom
   implicit none

#define DEVICE_ATTR device
   interface sum_prefix_custom
   function scan_int_wrapper(input, output, nptidx, N, stream) bind(C, name="scan_int_wrapper")
      use iso_c_binding
      use openacc
      integer(c_int) :: scan_int_wrapper
      integer(c_int), DEVICE_ATTR :: input(*)
      integer(c_int), DEVICE_ATTR :: output(*)
      integer(c_int), DEVICE_ATTR :: nptidx(*)
      integer(c_int) :: N
      integer(acc_handle_kind), value :: stream
   end function

!!   function scan_float_wrapper(input, output, nptidx, N, stream) bind(C, name="scan_float_wrapper")
!!      use iso_c_binding
!!      use openacc
!!      integer(c_int) :: scan_float_wrapper
!!      real(c_float), DEVICE_ATTR :: input(*)
!!      real(c_float), DEVICE_ATTR :: output(*)
!!      integer(c_int), DEVICE_ATTR :: nptidx(*)
!!      integer(c_int) :: N
!!      integer(acc_handle_kind), value :: stream
!!   end function
!!
!!   function scan_double_wrapper(input, output, nptidx, N, stream) bind(C, name="scan_double_wrapper")
!!      use iso_c_binding
!!      use openacc
!!      integer(c_int) :: scan_double_wrapper
!!      real(c_double), DEVICE_ATTR :: input(*)
!!      real(c_double), DEVICE_ATTR :: output(*)
!!      integer(c_int), DEVICE_ATTR :: nptidx(*)
!!      integer(c_int) :: N
!!      integer(acc_handle_kind), value :: stream
!!   end function
   end interface
contains

   subroutine scan_int(input, output, nptidx, nptidxdone)
      use iso_c_binding
      USE openacc
      use cutensorex, only: SUM_PREFIX
      integer(c_int), intent(in), DEVICE_ATTR:: input(:)
      integer(c_int), intent(out), DEVICE_ATTR:: output(:)
      integer(c_int), intent(out), DEVICE_ATTR:: nptidx(:)
      logical, intent(inout) :: nptidxdone
      INTEGER(acc_handle_kind) :: stream
      integer(c_int) :: retval 
      if (size(input) > 1600000 .and. size(input) < 10000000) then
        output = SUM_PREFIX(input)
      else
        stream = acc_get_cuda_stream(acc_async_sync)   
        retval = scan_int_wrapper(input, output, nptidx, size(input), stream)
        if (retval == 1) nptidxdone = .true.
      endif
   end subroutine

!!   subroutine scan_float(input, output, nptidx, nptidxdone)
!!      use iso_c_binding
!!      USE openacc
!!      use cutensorex, only: SUM_PREFIX
!!      real(c_float), intent(in), DEVICE_ATTR:: input(:)
!!      real(c_float), intent(out), DEVICE_ATTR:: output(:)
!!      integer(c_int), intent(out), DEVICE_ATTR:: nptidx(:)
!!      logical, intent(inout) :: nptidxdone
!!      INTEGER(acc_handle_kind) :: stream
!!      integer(c_int) :: retval 
!!      if (size(input) > 1600000 .and. size(input) < 10000000) then
!!        output = SUM_PREFIX(input)
!!      else
!!        stream = acc_get_cuda_stream(acc_async_sync)
!!        retval = scan_float_wrapper(input, output, nptidx, size(input), stream)
!!        if (retval == 1) nptidxdone = .true.
!!      endif
!!   end subroutine
!!
!!   subroutine scan_double(input, output, nptidx, nptidxdone)
!!      use iso_c_binding
!!      USE openacc
!!      use cutensorex, only: SUM_PREFIX
!!      real(c_double), intent(in), DEVICE_ATTR:: input(:)
!!      real(c_double), intent(out), DEVICE_ATTR:: output(:)
!!      integer(c_int), intent(out), DEVICE_ATTR:: nptidx(:)
!!      logical, intent(inout) :: nptidxdone
!!      INTEGER(acc_handle_kind) :: stream
!!      integer(c_int) :: retval 
!!      if (size(input) > 1600000 .and. size(input) < 10000000) then
!!        output = SUM_PREFIX(input)
!!      else
!!        stream = acc_get_cuda_stream(acc_async_sync)
!!        retval = scan_double_wrapper(input, output, nptidx, size(input), stream)
!!        if (retval == 1) nptidxdone = .true.
!!      endif
!!   end subroutine

end module sum_prefix_custom
