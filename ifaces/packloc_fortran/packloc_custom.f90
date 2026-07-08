module packloc_custom
   implicit none

#define DEVICE_ATTR , device

   interface packloc_custom
      module procedure packloc_int, packloc_char
   end interface

   interface
   function packloc_char_int_wrapper(input, idx, N, stream) bind(C, name="packloc_char_int_wrapper")
      use iso_c_binding
      use openacc
      integer(c_int)              :: packloc_char_int_wrapper
      logical*1      DEVICE_ATTR  :: input(*)
      integer(c_int) DEVICE_ATTR  :: idx(*)
      integer(c_int), value       :: N
      integer(acc_handle_kind), value :: stream
   end function

   function packloc_int_int_wrapper(input, idx, N, stream) bind(C, name="packloc_int_int_wrapper")
      use iso_c_binding
      use openacc
      integer(c_int)              :: packloc_int_int_wrapper
      logical        DEVICE_ATTR  :: input(*)
      integer(c_int) DEVICE_ATTR  :: idx(*)
      integer(c_int), value       :: N
      integer(acc_handle_kind), value :: stream
   end function
   end interface

contains

   subroutine packloc_int(input, idx, npti)
      USE openacc
      logical, intent(in)  DEVICE_ATTR  :: input(:)
      integer, intent(out) DEVICE_ATTR  :: idx(:)
      integer, intent(out)              :: npti
      INTEGER(acc_handle_kind)          :: stream
      integer(c_int)                    :: N
      N = size(input)
      stream = acc_get_cuda_stream(1)  !! NOTE: "1" means acc stream No.1 
      npti = packloc_int_int_wrapper(input, idx, N, stream)
   end subroutine

   subroutine packloc_char(input, idx, npti)
      USE openacc
      logical*1, intent(in) DEVICE_ATTR  :: input(:)
      integer, intent(out)  DEVICE_ATTR  :: idx(:)
      integer, intent(out)               :: npti
      INTEGER(acc_handle_kind)           :: stream
      integer(c_int)                     :: N
      N = size(input)
      stream = acc_get_cuda_stream(1)  !! NOTE: "1" means acc stream No.1 
      npti = packloc_char_int_wrapper(input, idx, N, stream)
   end subroutine

end module packloc_custom
