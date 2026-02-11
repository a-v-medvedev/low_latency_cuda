#define ongpu
 
module print_norms
ongpu use openacc
      use iso_fortran_env, only : real32, real64
      implicit none
      integer, parameter :: sp = real32
      integer, parameter :: dp = real64

      private

#define PRINT_CALL_GPU  print '(A, I3, A, I3, 3A, E20.15)', "norm: rank: ", mpprank, " i: ", i, " ", arrname, "@norm: ", l2_norm(arr)
#define DUMP_CALL_GPU  print '(A, I3, A, I3, 3A)', "dump: rank: ", mpprank, " i: ", i, " ", arrname, " not implemented"

   public :: print_norms_gpu, dump_array_gpu, l2_norm

   interface print_norms_gpu
      module procedure print_norms_1d_sp_gpu, print_norms_2d_sp_gpu, print_norms_3d_sp_gpu, &
                     & print_norms_4d_sp_gpu, print_norms_5d_sp_gpu, &
                     & print_norms_1d_dp_gpu, print_norms_2d_dp_gpu, print_norms_3d_dp_gpu, &
                     & print_norms_4d_dp_gpu, print_norms_5d_dp_gpu
   end interface print_norms_gpu

   interface dump_array_gpu
      module procedure dump_array_1d_sp_gpu, dump_array_2d_sp_gpu, dump_array_3d_sp_gpu, &
                     & dump_array_4d_sp_gpu, dump_array_5d_sp_gpu, &
                     & dump_array_1d_dp_gpu, dump_array_2d_dp_gpu, dump_array_3d_dp_gpu, &
                     & dump_array_4d_dp_gpu, dump_array_5d_dp_gpu
   end interface dump_array_gpu

   interface l2_norm
      module procedure l2_norm_1d_sp, l2_norm_2d_sp, l2_norm_3d_sp, l2_norm_4d_sp, l2_norm_5d_sp, &
                     & l2_norm_1d_dp, l2_norm_2d_dp, l2_norm_3d_dp, l2_norm_4d_dp, l2_norm_5d_dp
   end interface l2_norm

   contains

   function l2_norm_1d_sp(arr) result(norm) 
      real(sp), dimension(:), intent(in) :: arr       
      real(dp)                           :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
      norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_2d_sp(arr) result(norm) 
      real(sp), dimension(:,:), intent(in) :: arr       
      real(dp)                             :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
      norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_3d_sp(arr) result(norm) 
      real(sp), dimension(:,:,:), intent(in) :: arr       
      real(dp)                               :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_4d_sp(arr) result(norm) 
      real(sp), dimension(:,:,:,:), intent(in) :: arr       
      real(dp)                                 :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_5d_sp(arr) result(norm) 
      real(sp), dimension(:,:,:,:,:), intent(in) :: arr       
      real(dp)                                   :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_1d_dp(arr) result(norm) 
      real(dp), dimension(:), intent(in) :: arr       
      real(dp)                           :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_2d_dp(arr) result(norm) 
      real(dp), dimension(:,:), intent(in) :: arr       
      real(dp)                             :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_3d_dp(arr) result(norm) 
      real(dp), dimension(:,:,:), intent(in) :: arr       
      real(dp)                               :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_4d_dp(arr) result(norm) 
      real(dp), dimension(:,:,:,:), intent(in) :: arr       
      real(dp)                                 :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   function l2_norm_5d_dp(arr) result(norm) 
      real(dp), dimension(:,:,:,:,:), intent(in) :: arr       
      real(dp)                                   :: norm 
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
     norm = sqrt(sum(reshape(arr, [size(arr)])**2)) 
   end function

   subroutine print_norms_1d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU 
   end subroutine

   subroutine print_norms_2d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU 
   end subroutine

   subroutine print_norms_3d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU 
   end subroutine

   subroutine print_norms_4d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU 
   end subroutine

   subroutine print_norms_5d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU 
   end subroutine

   subroutine print_norms_1d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU
   end subroutine

   subroutine print_norms_2d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU
   end subroutine

   subroutine print_norms_3d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU
   end subroutine

   subroutine print_norms_4d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU
   end subroutine

   subroutine print_norms_5d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      PRINT_CALL_GPU
   end subroutine


   subroutine dump_array_1d_sp_gpu(mpprank, counter, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: counter
      real(sp), dimension(:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      INTEGER :: i
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
      do i=1,size(arr,1);
         print '(A, I3, A, I3, 3A, I, A, E20.15)', "dump: rank: ", mpprank, " i: ", counter, " ", arrname, " @idx: ", i, " @value: ", arr(i)
      end do
   end subroutine

   subroutine dump_array_2d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      DUMP_CALL_GPU 
   end subroutine

   subroutine dump_array_3d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      DUMP_CALL_GPU 
   end subroutine

   subroutine dump_array_4d_sp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(sp), dimension(:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      DUMP_CALL_GPU 
   end subroutine

   subroutine dump_array_5d_sp_gpu(mpprank, counter, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: counter
      real(sp), dimension(:,:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      INTEGER :: i, j, k, l, f
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
      do f=1,size(arr,5); do l=1,size(arr,4); do k=1,size(arr,3); do j=1,size(arr,2); do i=1,size(arr,1);
         print '(A, I3, A, I3, 3A, 5I, A, E20.15)', "dump: rank: ", mpprank, " i: ", counter, " ", arrname, " @idx: ", i,j,k,l,f, " @value: ", arr(i,j,k,l,f)
      end do; end do; end do; end do; end do
   end subroutine

   subroutine dump_array_1d_dp_gpu(mpprank, counter, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: counter
      real(dp), dimension(:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      INTEGER :: i
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
      do i=1,size(arr,1);
         print '(A, I3, A, I3, 3A, I, A, E20.15)', "dump: rank: ", mpprank, " i: ", counter, " ", arrname, " @idx: ", i, " @value: ", arr(i)
      end do
   end subroutine

   subroutine dump_array_2d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      DUMP_CALL_GPU
   end subroutine

   subroutine dump_array_3d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      DUMP_CALL_GPU
   end subroutine

   subroutine dump_array_4d_dp_gpu(mpprank, i, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: i
      real(dp), dimension(:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname 
      DUMP_CALL_GPU
   end subroutine

   subroutine dump_array_5d_dp_gpu(mpprank, counter, arr, arrname) 
      implicit none 
      integer, intent(in)       :: mpprank 
      integer, intent(in)       :: counter
      real(dp), dimension(:,:,:,:,:), intent(in)  :: arr   
      character(len=*), intent(in)  :: arrname
      INTEGER :: i, j, k, l, f
ongpu if (acc_is_present(arr)) then
ongpu !$acc update host(arr)
ongpu endif
      do f=1,size(arr,5); do l=1,size(arr,4); do k=1,size(arr,3); do j=1,size(arr,2); do i=1,size(arr,1);
         print '(A, I3, A, I3, 3A, 5I, A, E20.15)', "dump: rank: ", mpprank, " i: ", counter, " ", arrname, " @idx: ", i,j,k,l,f, " @value: ", arr(i,j,k,l,f)
      end do; end do; end do; end do; end do
   end subroutine

end module print_norms
