!! nvfortran -cuda -acc=gpu -O3 -o sc_test sc_test.f90 -cudalib=cutensor

PROGRAM SCTEST
   use iso_fortran_env, only : real32, real64
   use print_norms
   implicit none
   integer, parameter :: sp = real32
   integer, parameter :: dp = real64
   INTEGER               :: jpi, jpj
   INTEGER               :: ncycles = 1, ncycles_tab = 20
   REAL(dp), ALLOCATABLE :: x_(:,:), y_(:,:)
   INTEGER               :: threshold
   INTEGER               :: input

   call get_cmdline_arg(1, input)

   select case (input)
   case(1) 
      write (*,*) "-- generate"
      jpi = 524; jpj = 704
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = generate(x_, y_, 6)
   case(2)
      write (*,*) "-- read from file: orca025_idx_135_352__14395.txt" 
      jpi = 135; jpj = 352
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_135_352__14395.txt', 14395)
   case(3)
      write (*,*) "-- read from file: orca025_idx_364_526__1326.txt" 
      jpi = 364; jpj = 526
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_364_526__1326.txt', 1326)
   case(4)
      write (*,*) "-- read from file: orca025_idx_364_526__17176.txt" 
      jpi = 364; jpj = 526
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_364_526__17176.txt', 17176)
   case(5)
      write (*,*) "-- read from file: orca025_idx_364_526__20196.txt" 
      jpi = 364; jpj = 526
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_364_526__20196.txt', 20196)
   case(6)
      write (*,*) "-- read from file: orca025_idx_364_526__22811.txt" 
      jpi = 364; jpj = 526
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_364_526__22811.txt', 22811)
   case(7)
      write (*,*) "-- read from file: orca025_idx_364_526__25868.txt" 
      jpi = 364; jpj = 526
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_364_526__25868.txt', 25868)
   case(8)
      write (*,*) "-- read from file: orca025_idx_364_526__34783.txt" 
      jpi = 364; jpj = 526
      allocate(x_(jpi,jpj),y_(jpi,jpj))
      threshold = read_file(x_, y_, 'orca025_idx_364_526__34783.txt', 34783)
   case default
      stop "Unknown input index. Stop."
   end select

   write (*,*) "jpi: ", jpi
   write (*,*) "jpj: ", jpj
   write (*,*) "L2_NORM x: ", l2_norm(x_)
   write (*,*) "L2_NORM y: ", l2_norm(y_)

   !! warm-ups
   call test1(x_,y_)
   !! NOTE: test2 is slow and not interesting
   !!call test2(x_,y_)
   call test3(x_,y_)
#if defined WITH_CUTENSOREX
   call test4(x_,y_)
   call test5(x_,y_)
   call test6(x_, y_)
#endif
   call test7(x_, y_)
   call test8(x_, y_)
   call test9(x_, y_)
  
   !! actual tests 
   ncycles = 10
   write (*,*) "-- test1:"
   call test1(x_, y_)
   !! NOTE: test2 is slow and not interesting
   !!write (*,*) "-- test2:"
   !!call test2(x_, y_)
   write (*,*) "-- test3:"
   call test3(x_, y_)
#if defined WITH_CUTENSOREX
   write (*,*) "-- test4:"
   call test4(x_, y_)
   write (*,*) "-- test5:"
   call test5(x_, y_)
   write (*,*) "-- test6:"
   call test6(x_, y_)
#endif
   write (*,*) "-- test7:"
   call test7(x_, y_)
!!   write (*,*) "-- test8:"
!!   call test8(x_, y_)
!!   write (*,*) "-- test9:"
!!   call test9(x_, y_)

   deallocate(x_, y_)

CONTAINS 
   subroutine get_cmdline_arg(idx, value)
      implicit none
      integer, intent(in)  :: idx
      integer, intent(out) :: value

      character(len=64) :: arg
      integer :: ios

      if (command_argument_count() < idx) then
         value = 1
      else
         call get_command_argument(idx, arg)

         read(arg, *, iostat=ios) value
         if (ios /= 0) stop "Invalid integer argument. Stop."
      end if
   end subroutine

   FUNCTION read_file(x, y, file_name, n)
      REAL(dp)     :: x(:,:), y(:,:)
      CHARACTER(*) :: file_name
      INTEGER      :: n
      REAL(dp) :: read_file
      INTEGER, ALLOCATABLE :: idx(:)
      INTEGER      :: ios, i, j, p, npti
      REAL(dp)     :: thr = 1000.0_dp

      allocate(idx(n))

      open(10, file=file_name, status='old', action='read')
      read(10, *, iostat=ios) idx
      if (ios /= 0) then
         write (*,*) 'File reading error: ', file_name
         stop "File reading error. Stop."
      endif
 
      npti = 1
      do j=1,size(x,2); do i=1,size(x,1)
         p = (j - 1) * size(x,1) + i
         if (idx(npti) == p) then
            x(i,j) = thr + 1.0_dp
            npti = npti + 1
         endif
         y(i,j) = j
      enddo; enddo

      write (*,*) 'Values in file count=', count(x>thr)
      write (*,*) 'Success:', (count(x>thr) == n)

      deallocate(idx)
      read_file = thr 
   END FUNCTION


   FUNCTION generate(x, y, n)
      REAL(dp) :: x(:,:), y(:,:)
      INTEGER  :: n
      REAL(dp) :: generate
      INTEGER  :: FIB(10) = (/1, 1, 2, 3, 5, 8, 13, 21, 34, 55/)
      INTEGER  :: THR(10) = (/10, 10, 30, 50, 1000, 100000, 1000000, 10000000, 1000000000, 1000000000/)
      INTEGER  :: ji, jj, a1 = 13, b1 = 8, c1 = 1, a2 = 1, b2 = 1, c2 = 1
      a1 = FIB(n)
      b1 = FIB(n-1)
      DO jj = 1, jpj 
         DO ji = 1, jpi
            if (MOD((jj - 1) * jpi + ji, a1) == 0) then
               c1 = b1; b1 = a1; a1 = b1 + c1
               a2 = 1; b2 = 1; c2 = 1
            endif
            if (a2 > a1) then
              a2 = 1; b2 = 1; c2 = 1
            else
              c2 = b2; b2 = a2; a2 = b2 + c2
            endif
            x(ji,jj) = a2
            y(ji,jj) = jj  
         END DO 
      END DO 
      generate = real(THR(n))
   END FUNCTION

   ! original CPU code that creates a map (the nptidx array) and then does the 2d->1d compaction 
   ! according to the map using the tab_2d_1d() subroutine
   SUBROUTINE test1(x, y)
      REAL(dp) :: x(:,:), y(:,:)
      INTEGER :: npti 
      INTEGER, ALLOCATABLE :: nptidx(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      REAL(dp) :: tmp1, tmp2, tmp3
      INTEGER :: j, k, m, count_rate, count_start, count_end
      ALLOCATE(nptidx(jpi*jpj))
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size
      nptidx = 0
      z1 = 0; z2 = 0; z3 = 0; z4 = 0; z5 = 0; z6 = 0; z7 = 0; z8 = 0; z9 = 0; z10 = 0;
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         npti = 0
         DO jj = 1, jpj 
            DO ji = 1, jpi 
               IF ( x(ji,jj) > threshold ) THEN 
                  npti = npti + 1 
                  nptidx( npti ) = (jj - 1) * jpi + ji 
               ENDIF 
            END DO 
         END DO 
         DO m = 1, ncycles_tab
            call tab_2d_1d(npti, nptidx, z1, y)
            call tab_2d_1d(npti, nptidx, z2, y)
            call tab_2d_1d(npti, nptidx, z3, y)
            call tab_2d_1d(npti, nptidx, z4, y)
            call tab_2d_1d(npti, nptidx, z5, y)
            call tab_2d_1d(npti, nptidx, z6, y)
            call tab_2d_1d(npti, nptidx, z7, y)
            call tab_2d_1d(npti, nptidx, z8, y)
            call tab_2d_1d(npti, nptidx, z9, y)
            call tab_2d_1d(npti, nptidx, z10, y)
         END DO
         DO m = 1, ncycles_tab
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
         END DO
         DO m = 1, ncycles_tab
            call tab_1d_2d(npti, nptidx, z1, y)
            call tab_1d_2d(npti, nptidx, z2, y)
            call tab_1d_2d(npti, nptidx, z3, y)
            call tab_1d_2d(npti, nptidx, z4, y)
            call tab_1d_2d(npti, nptidx, z5, y)
            call tab_1d_2d(npti, nptidx, z6, y)
            call tab_1d_2d(npti, nptidx, z7, y)
            call tab_1d_2d(npti, nptidx, z8, y)
            call tab_1d_2d(npti, nptidx, z9, y)
            call tab_1d_2d(npti, nptidx, z10, y)
         END DO
      end do
      call system_clock(count_end)
      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx)
         write (*,*) "L2_NORM z: ", l2_norm(z10)
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)
   END SUBROUTINE

   ! the CPU code doing the same that creates a 1D mask array using reshape intrinsic, then uses this mask for the
   ! pack intrinsic to do actual compaction
   ! NOTE: this code is slow, so it is left here only for historic reasons
   SUBROUTINE test2(x, y)
      REAL(dp) :: x(:,:), y(:,:)
      INTEGER :: i, npti 
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      INTEGER :: k, m, count_rate, count_start, count_end

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size

      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         mask    = reshape(x > threshold, [jpi*jpj])
         npti    = count(mask)
         nptidx  = pack([(i, i=1,jpi*jpj)], mask)
         DO m = 1, ncycles_tab
            z1 = pack(reshape(y, [jpi*jpj]), mask)
            z2 = pack(reshape(y, [jpi*jpj]), mask)
            z3 = pack(reshape(y, [jpi*jpj]), mask)
            z4 = pack(reshape(y, [jpi*jpj]), mask)
            z5 = pack(reshape(y, [jpi*jpj]), mask)
            z6 = pack(reshape(y, [jpi*jpj]), mask)
            z7 = pack(reshape(y, [jpi*jpj]), mask)
            z8 = pack(reshape(y, [jpi*jpj]), mask)
            z9 = pack(reshape(y, [jpi*jpj]), mask)
            z10= pack(reshape(y, [jpi*jpj]), mask)
         END DO
      end do
      call system_clock(count_end)
      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx)
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)
   END SUBROUTINE

   ! the GPU version of test1 -- calls the packloc_custom() custom subroutine defined below, also calls gpu versions of tab_*() 
   ! subroutines. Each tab_*() subroutine here is a separate GPU kernels annotated with async(1)
   SUBROUTINE test3(x, y)
      REAL(dp) :: x(:,:), y(:,:)
      INTEGER :: npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      INTEGER :: j, k, m, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size
      nptidx = 0
      z1 = 0; z2 = 0; z3 = 0; z4 = 0; z5 = 0; z6 = 0; z7 = 0; z8 = 0; z9 = 0; z10 = 0;

      !$acc data copy(nptidx,z10,x,y) create(z1,z2,z3,z4,z5,z6,z7,z8,z9)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         call packloc_custom(x, threshold, nptidx, npti)
         DO m = 1, ncycles_tab
            call tab_2d_1d_gpu(npti, nptidx, z1, y)
            call tab_2d_1d_gpu(npti, nptidx, z2, y)
            call tab_2d_1d_gpu(npti, nptidx, z3, y)
            call tab_2d_1d_gpu(npti, nptidx, z4, y)
            call tab_2d_1d_gpu(npti, nptidx, z5, y)
            call tab_2d_1d_gpu(npti, nptidx, z6, y)
            call tab_2d_1d_gpu(npti, nptidx, z7, y)
            call tab_2d_1d_gpu(npti, nptidx, z8, y)
            call tab_2d_1d_gpu(npti, nptidx, z9, y)
            call tab_2d_1d_gpu(npti, nptidx, z10, y)
         END DO
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,tmp1,tmp2,tmp3) async(1)
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
            !$acc end parallel loop
         END DO
         DO m = 1, ncycles_tab
            call tab_1d_2d_gpu(npti, nptidx, z1, y)
            call tab_1d_2d_gpu(npti, nptidx, z2, y)
            call tab_1d_2d_gpu(npti, nptidx, z3, y)
            call tab_1d_2d_gpu(npti, nptidx, z4, y)
            call tab_1d_2d_gpu(npti, nptidx, z5, y)
            call tab_1d_2d_gpu(npti, nptidx, z6, y)
            call tab_1d_2d_gpu(npti, nptidx, z7, y)
            call tab_1d_2d_gpu(npti, nptidx, z8, y)
            call tab_1d_2d_gpu(npti, nptidx, z9, y)
            call tab_1d_2d_gpu(npti, nptidx, z10, y)
         END DO
         !$acc wait(1)
      end do
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)
   END SUBROUTINE

#if defined WITH_CUTENSOREX 
   ! the GPU version of test1 -- expected to call packloc(), pack() and unpack() from the cutensorex, uses pointer ll_condition_1d instead of reshape
   SUBROUTINE test4(x, y)
      USE cutensorex, only: pack, unpack, packloc
      REAL(dp), INTENT(in) :: x(:,:)
      REAL(dp), INTENT(out) :: y(:,:)
      LOGICAL, ALLOCATABLE :: ll_condition(:,:)
      TARGET :: ll_condition
      INTEGER :: npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      LOGICAL, POINTER :: ll_condition_1d(:)
      INTEGER :: j, k, m, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size
      allocate(ll_condition(size(x,1),size(x,2)), source=.false.)

      !$acc data copy(x,y,nptidx,z10) create(z2,z3,z4,z5,z6,z7,z8,z9,z1,ll_condition)
      ll_condition_1d(1:jpi*jpj) => ll_condition(:,:)
      !$acc parallel loop gang vector collapse(2)
      do j=1,size(x,2) ; do k=1,size(x,1) ; ll_condition(k,j) = (x(k,j) > threshold) ; end do ; end do
      !$acc end parallel loop 
      !!!$acc kernels
      !!where (x > threshold)  ll_condition = .true.
      !!!$acc end kernels 
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc host_data use_device(ll_condition_1d, nptidx)
         nptidx  = packloc(ll_condition_1d, count=npti)
         !$acc end host_data
         DO m = 1, ncycles_tab
            !$acc host_data use_device(x, y, z1)
            z1 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z2)
            z2 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z3)
            z3 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z4)
            z4 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z5)
            z5 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z6)
            z6 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z7)
            z7 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z8)
            z8 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z9)
            z9 = pack(y, ll_condition)
            !$acc end host_data
            !$acc host_data use_device(x, y, z10)
            z10= pack(y, ll_condition)
            !$acc end host_data
         END DO
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,tmp1,tmp2,tmp3) async(1)
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
            !$acc end parallel loop
         END DO
         !$acc wait(1)
         DO m = 1, ncycles_tab
            !$acc host_data use_device(x, y, z1)
            y = unpack(z1, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z2)
            y = unpack(z2, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z3)
            y = unpack(z3, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z4)
            y = unpack(z4, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z5)
            y = unpack(z5, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z6)
            y = unpack(z6, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z7)
            y = unpack(z7, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z8)
            y = unpack(z8, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z9)
            y = unpack(z9, ll_condition, y)
            !$acc end host_data
            !$acc host_data use_device(x, y, z10)
            y = unpack(z10, ll_condition, y)
            !$acc end host_data
         END DO
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)
   END SUBROUTINE
 
   ! variation of test4 -- still use packloc from tensorex, but replace pack() calls with GPU versions of tab_*() calls, as in test3
   SUBROUTINE test5(x, y)
      USE cutensorex, only: pack, packloc
      REAL(dp) :: x(:,:), y(:,:)
      LOGICAL, ALLOCATABLE :: ll_condition(:,:)
      TARGET :: ll_condition
      INTEGER :: npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      LOGICAL, POINTER :: ll_condition_1d(:)
      INTEGER :: j, k, m, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(ll_condition(jpi,jpj))
      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size

      !$acc data copy(x,y,nptidx,z10) create(z2,z3,z4,z5,z6,z7,z8,z9,z1,ll_condition)
      ll_condition_1d(1:jpi*jpj) => ll_condition(:,:)
      !$acc parallel loop gang vector collapse(2)
      do j=1,size(x,2) ; do k=1,size(x,1) ; ll_condition(k,j) = (x(k,j) > threshold) ; end do ; end do
      !$acc end parallel loop 
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc host_data use_device(ll_condition_1d, nptidx)
         nptidx  = packloc(ll_condition_1d, count=npti)
         !$acc end host_data
         DO m = 1, ncycles_tab
            call tab_2d_1d_gpu(npti, nptidx, z1, y)
            call tab_2d_1d_gpu(npti, nptidx, z2, y)
            call tab_2d_1d_gpu(npti, nptidx, z3, y)
            call tab_2d_1d_gpu(npti, nptidx, z4, y)
            call tab_2d_1d_gpu(npti, nptidx, z5, y)
            call tab_2d_1d_gpu(npti, nptidx, z6, y)
            call tab_2d_1d_gpu(npti, nptidx, z7, y)
            call tab_2d_1d_gpu(npti, nptidx, z8, y)
            call tab_2d_1d_gpu(npti, nptidx, z9, y)
            call tab_2d_1d_gpu(npti, nptidx, z10, y)
         END DO
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,tmp1,tmp2,tmp3) async(1)
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
            !$acc end parallel loop
         END DO
         DO m = 1, ncycles_tab
            call tab_1d_2d_gpu(npti, nptidx, z1, y)
            call tab_1d_2d_gpu(npti, nptidx, z2, y)
            call tab_1d_2d_gpu(npti, nptidx, z3, y)
            call tab_1d_2d_gpu(npti, nptidx, z4, y)
            call tab_1d_2d_gpu(npti, nptidx, z5, y)
            call tab_1d_2d_gpu(npti, nptidx, z6, y)
            call tab_1d_2d_gpu(npti, nptidx, z7, y)
            call tab_1d_2d_gpu(npti, nptidx, z8, y)
            call tab_1d_2d_gpu(npti, nptidx, z9, y)
            call tab_1d_2d_gpu(npti, nptidx, z10, y)
         END DO
         !$acc wait(1)
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, ll_condition)
   END SUBROUTINE

   ! variation of test5 -- replace simple GPU versions of tab_*() calls with a single-kernel approach: here we have a complex loop
   ! that calls the devece routines for tab_*(), all in parallel
   SUBROUTINE test6(x, y)
      USE cutensorex, only: pack, packloc
      REAL(dp) :: x(:,:), y(:,:)
      TARGET :: x
      REAL(dp), POINTER :: x1d(:)
      INTEGER :: npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      LOGICAL, ALLOCATABLE :: ll_condition_1d(:)
      INTEGER :: j, k, m, p, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(ll_condition_1d(jpi*jpj),source=.false.)
      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size
     
      x1d(1:jpi*jpj) => x(:,:)
      do k=1,size(ll_condition_1d); ll_condition_1d(k)=(x1d(k) > threshold); enddo

      !$acc data copy(x,y,ll_condition_1d,nptidx,z10) create(z2,z3,z4,z5,z6,z7,z8,z9,z1)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc host_data use_device(ll_condition_1d, nptidx)
         nptidx  = packloc(ll_condition_1d, count=npti)
         !$acc end host_data
         DO m = 1, ncycles_tab
            !$acc parallel loop gang(dim:2) num_gangs(256,10) async(1)
            do p=1,10 
            select case(p)
            case(1)  ; call tab_2d_1d_device(npti, nptidx, z1, y)
            case(2)  ; call tab_2d_1d_device(npti, nptidx, z2, y)
            case(3)  ; call tab_2d_1d_device(npti, nptidx, z3, y)
            case(4)  ; call tab_2d_1d_device(npti, nptidx, z4, y)
            case(5)  ; call tab_2d_1d_device(npti, nptidx, z5, y)
            case(6)  ; call tab_2d_1d_device(npti, nptidx, z6, y)
            case(7)  ; call tab_2d_1d_device(npti, nptidx, z7, y)
            case(8)  ; call tab_2d_1d_device(npti, nptidx, z8, y)
            case(9)  ; call tab_2d_1d_device(npti, nptidx, z9, y)
            case(10) ; call tab_2d_1d_device(npti, nptidx, z10, y)
            end select
            end do
         END DO
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,tmp1,tmp2,tmp3) async(1)
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
            !$acc end parallel loop
         END DO
         DO m = 1, ncycles_tab
         !$acc parallel loop gang(dim:2) num_gangs(256,10) async(1)
            do p=1,10
            select case(p)
            case(1)  ; call tab_1d_2d_device(npti, nptidx, z1, y)
            case(2)  ; call tab_1d_2d_device(npti, nptidx, z2, y)
            case(3)  ; call tab_1d_2d_device(npti, nptidx, z3, y)
            case(4)  ; call tab_1d_2d_device(npti, nptidx, z4, y)
            case(5)  ; call tab_1d_2d_device(npti, nptidx, z5, y)
            case(6)  ; call tab_1d_2d_device(npti, nptidx, z6, y)
            case(7)  ; call tab_1d_2d_device(npti, nptidx, z7, y)
            case(8)  ; call tab_1d_2d_device(npti, nptidx, z8, y)
            case(9)  ; call tab_1d_2d_device(npti, nptidx, z9, y)
            case(10) ; call tab_1d_2d_device(npti, nptidx, z10, y)
            end select
            end do
         END DO
         !$acc wait(1)
      end do
      call system_clock(count_end)
      !! NOTE: the extra copy because the previous copies may overwrite y in random order 
      !! This is for having a correct y norm
      call tab_1d_2d_gpu(npti, nptidx, z10, y)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, ll_condition_1d)
   END SUBROUTINE
 #endif

   ! variation of test6 -- replace packloc with packloc_custom
   SUBROUTINE test7(x, y)
      REAL(dp) :: x(:,:), y(:,:)
      TARGET :: x
      INTEGER :: npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      INTEGER :: j, k, m, p, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size

      !$acc data copy(x,y,nptidx,z10) create(z2,z3,z4,z5,z6,z7,z8,z9,z1)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         call packloc_custom(x, threshold, nptidx, npti)
         DO m = 1, ncycles_tab
            !$acc parallel loop gang(dim:2) num_gangs(256,10) async(1)
            do p=1,10
            select case(p)
            case(1)  ; call tab_2d_1d_device(npti, nptidx, z1, y)
            case(2)  ; call tab_2d_1d_device(npti, nptidx, z2, y)
            case(3)  ; call tab_2d_1d_device(npti, nptidx, z3, y)
            case(4)  ; call tab_2d_1d_device(npti, nptidx, z4, y)
            case(5)  ; call tab_2d_1d_device(npti, nptidx, z5, y)
            case(6)  ; call tab_2d_1d_device(npti, nptidx, z6, y)
            case(7)  ; call tab_2d_1d_device(npti, nptidx, z7, y)
            case(8)  ; call tab_2d_1d_device(npti, nptidx, z8, y)
            case(9)  ; call tab_2d_1d_device(npti, nptidx, z9, y)
            case(10) ; call tab_2d_1d_device(npti, nptidx, z10, y)
            end select
            end do
         END DO
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,tmp1,tmp2,tmp3) async(1)
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
            !$acc end parallel loop
         END DO
         DO m = 1, ncycles_tab
         !$acc parallel loop gang(dim:2) num_gangs(256,10) async(1)
            do p=1,10
            select case(p)
            case(1)  ; call tab_1d_2d_device(npti, nptidx, z1, y)
            case(2)  ; call tab_1d_2d_device(npti, nptidx, z2, y)
            case(3)  ; call tab_1d_2d_device(npti, nptidx, z3, y)
            case(4)  ; call tab_1d_2d_device(npti, nptidx, z4, y)
            case(5)  ; call tab_1d_2d_device(npti, nptidx, z5, y)
            case(6)  ; call tab_1d_2d_device(npti, nptidx, z6, y)
            case(7)  ; call tab_1d_2d_device(npti, nptidx, z7, y)
            case(8)  ; call tab_1d_2d_device(npti, nptidx, z8, y)
            case(9)  ; call tab_1d_2d_device(npti, nptidx, z9, y)
            case(10) ; call tab_1d_2d_device(npti, nptidx, z10, y)
            end select
            end do
         END DO
         !$acc wait(1)
      end do
      call system_clock(count_end)
      !! NOTE: the extra copy because the previous copies may overwrite y in random order 
      !! This is for having a correct y norm
      call tab_1d_2d_gpu(npti, nptidx, z10, y)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)
   END SUBROUTINE

   ! variation of test7 -- use logical array for calculation kernels instead of compacted data. Measure only calculation kernels
   ! time 
   SUBROUTINE test8(x, y)
      REAL(dp) :: x(:,:), y(:,:)
      LOGICAL, ALLOCATABLE :: ll_condition(:,:)
      TARGET :: x
      REAL(dp), ALLOCATABLE :: y1(:,:), y2(:,:), y3(:,:), y4(:,:), y5(:,:), y6(:,:), y7(:,:), y8(:,:), y9(:,:), y10(:,:)

      INTEGER :: j, k, m, p, npti, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(ll_condition(jpi,jpj))
      allocate(y1(jpi,jpj),y2(jpi,jpj),y3(jpi,jpj),y4(jpi,jpj),y5(jpi,jpj),y6(jpi,jpj),y7(jpi,jpj),y8(jpi,jpj),y9(jpi,jpj),y10(jpi,jpj))
      y1 = y; y2 = y; y3 = y; y4 = y; y5 = y; y6 = y; y7 = y; y8 = y; y9 = y; y10 = y

      !$acc data copy(y10) copyin(x,y2,y3,y4,y5,y6,y7,y8,y9,y1) create(ll_condition)
      !$acc parallel loop gang vector collapse(2) async(1)
      do j=1,size(x,2) ; do k=1,size(x,1) ; ll_condition(k,j) = (x(k,j) > threshold) ; end do ; end do
      !$acc end parallel loop 
      !$acc kernels
      npti = count(ll_condition)
      !$acc end kernels  
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do p = 1, ncycles
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,k,tmp1,tmp2,tmp3) collapse(2) async(1)
            do j=1,size(x,2) ; do k=1,size(x,1) 
               if (ll_condition(k,j)) then
                  tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
                  if (j>1)  tmp1 = real(y4(k,j)) / real(y5(k,j-1))
                  if (j>1)  tmp2 = real(y6(k,j)) / real(y7(k,j-1))
                  if (j>1)  tmp3 = 3.0_dp * real(y8(k,j)) / real(y9(k,j-1))
                  y10(k,j) = int(real(y10(k,j)) + real(y10(k,j)) / real(y2(k,j)) + real(y3(k,j)) * EXP(-0.1_dp * tmp1) + real(y1(k,j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
               endif
            end do; end do
            !$acc end parallel loop
         END DO
         !$acc wait(1)
      end do
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "L2_NORM y: ", l2_norm(y1)
         write (*,*) "time calc only: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
      DEALLOCATE(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, ll_condition)
   END SUBROUTINE

   ! variation of test7 -- original 1d scheme, but measure time for calculation kernels only
   SUBROUTINE test9(x, y)
      REAL(dp) :: x(:,:), y(:,:)
      TARGET :: x
      INTEGER :: npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      REAL(dp), ALLOCATABLE :: z1(:), z2(:), z3(:), z4(:), z5(:), z6(:), z7(:), z8(:), z9(:), z10(:)
      INTEGER :: j, k, m, p, count_rate, count_start, count_end
      REAL(dp) :: tmp1, tmp2, tmp3

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z1(jpi*jpj),z2(jpi*jpj),z3(jpi*jpj),z4(jpi*jpj),z5(jpi*jpj),z6(jpi*jpj),z7(jpi*jpj),z8(jpi*jpj),z9(jpi*jpj),z10(jpi*jpj))       ! largest possible size

      !$acc data copy(x,y,nptidx,z10) create(z1,z2,z3,z4,z5,z6,z7,z8,z9)
      call packloc_custom(x, threshold, nptidx, npti)


      !$acc parallel loop gang(dim:2) num_gangs(2048) async(1)
      do p=1,10
      select case(p)
      case(1)  ; call tab_2d_1d_device(npti, nptidx, z1, y)
      case(2)  ; call tab_2d_1d_device(npti, nptidx, z2, y)
      case(3)  ; call tab_2d_1d_device(npti, nptidx, z3, y)
      case(4)  ; call tab_2d_1d_device(npti, nptidx, z4, y)
      case(5)  ; call tab_2d_1d_device(npti, nptidx, z5, y)
      case(6)  ; call tab_2d_1d_device(npti, nptidx, z6, y)
      case(7)  ; call tab_2d_1d_device(npti, nptidx, z7, y)
      case(8)  ; call tab_2d_1d_device(npti, nptidx, z8, y)
      case(9)  ; call tab_2d_1d_device(npti, nptidx, z9, y)
      case(10) ; call tab_2d_1d_device(npti, nptidx, z10, y)
      end select
      end do
      !$acc wait(1)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         DO m = 1, ncycles_tab
            !$acc parallel loop gang vector default(present) private(j,tmp1,tmp2,tmp3) async(1)
            do j=1,npti
               tmp1 = 1.0_dp; tmp2 = 2.0_dp; tmp3 = 3.0_dp
               if (j>1)  tmp1 = real(z4(j)) / real(z5(j-1))
               if (j>1)  tmp2 = real(z6(j)) / real(z7(j-1))
               if (j>1)  tmp3 = 3.0_dp * real(z8(j)) / real(z9(j-1))
               z10(j) = int(real(z10(j)) + real(z10(j)) / real(z2(j)) + real(z3(j)) * EXP(-0.1_dp * tmp1) + real(z1(j)) * EXP(-0.1_dp * tmp2) * SQRT(tmp3))
            end do
            !$acc end parallel loop
         END DO
         !$acc wait(1)
      end do
      call system_clock(count_end)
      !$acc parallel loop gang(dim:2) num_gangs(2048) async(1)
      do p=1,10
      select case(p)
      case(1)  ; call tab_1d_2d_device(npti, nptidx, z1, y)
      case(2)  ; call tab_1d_2d_device(npti, nptidx, z2, y)
      case(3)  ; call tab_1d_2d_device(npti, nptidx, z3, y)
      case(4)  ; call tab_1d_2d_device(npti, nptidx, z4, y)
      case(5)  ; call tab_1d_2d_device(npti, nptidx, z5, y)
      case(6)  ; call tab_1d_2d_device(npti, nptidx, z6, y)
      case(7)  ; call tab_1d_2d_device(npti, nptidx, z7, y)
      case(8)  ; call tab_1d_2d_device(npti, nptidx, z8, y)
      case(9)  ; call tab_1d_2d_device(npti, nptidx, z9, y)
      case(10) ; call tab_1d_2d_device(npti, nptidx, z10, y)
      end select
      end do
      !$acc wait(1)
      !! NOTE: the extra copy because the previous copies may overwrite y in random order 
      !! This is for having a correct y norm
      call tab_1d_2d_gpu(npti, nptidx, z10, y)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "L2_NORM z: ", l2_norm(z10(1:npti))
         write (*,*) "L2_NORM y: ", l2_norm(y)
         write (*,*) "time calc only: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
      DEALLOCATE(nptidx, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10)
   END SUBROUTINE


   SUBROUTINE tab_2d_1d( ndim1d, tab_ind, tab1d, tab2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d   ! 1d size
      INTEGER , DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind  ! input index
      REAL(dp), DIMENSION(:,:), INTENT(in   ) ::   tab2d    ! input 2D field
      REAL(dp), DIMENSION(ndim1d) , INTENT(out) ::   tab1d    ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , size(tab2d,1) ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / size(tab2d,1) + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO
   END SUBROUTINE tab_2d_1d

   SUBROUTINE tab_2d_1d_gpu( ndim1d, tab_ind, tab1d, tab2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   ndim1d   ! 1d size
      INTEGER, DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind  ! input index
      REAL(dp), DIMENSION(:,:), INTENT(in   ) ::   tab2d    ! input 2D field
      REAL(dp), DIMENSION(ndim1d) , INTENT(out)   ::   tab1d    ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      !$acc parallel loop gang vector default(present) private(jid,jjd) async(1)
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , size(tab2d,1) ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / size(tab2d,1) + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO
      !$acc end parallel loop
   END SUBROUTINE tab_2d_1d_gpu

   SUBROUTINE tab_2d_1d_device( ndim1d, tab_ind, tab1d, tab2d )
      !$acc routine gang(dim:1) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   ndim1d   ! 1d size
      INTEGER, DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind  ! input index
      REAL(dp), DIMENSION(:,:), INTENT(in   ) ::   tab2d    ! input 2D field
      REAL(dp), DIMENSION(ndim1d) , INTENT(out)   ::   tab1d    ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      !$acc loop gang(dim:1) vector private(jid,jjd)
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , size(tab2d,1) ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / size(tab2d,1) + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO
      !$acc end loop
   END SUBROUTINE tab_2d_1d_device

   SUBROUTINE tab_1d_2d( ndim1d, tab_ind, tab1d, tab2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_1d_2d  ***
      !! restore ji,jj dimensions from a packed form using tab_ind as a map
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d    ! 1D size
      INTEGER, DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind   ! input index
      REAL(dp), DIMENSION(:) ,      INTENT(in   ) ::   tab1d     ! input 1D field
      REAL(dp), DIMENSION(:,:), INTENT(out)   ::   tab2d     ! output 2D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1 ,  size(tab2d,1) ) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / size(tab2d,1)  + 1
         tab2d(jid, jjd) = tab1d( jn)
      END DO
   END SUBROUTINE tab_1d_2d

   SUBROUTINE tab_1d_2d_gpu( ndim1d, tab_ind, tab1d, tab2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_1d_2d  ***
      !! restore ji,jj dimensions from a packed form using tab_ind as a map
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d    ! 1D size
      INTEGER, DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind   ! input index
      REAL(dp), DIMENSION(:) ,      INTENT(in   ) ::   tab1d     ! input 1D field
      REAL(dp), DIMENSION(:,:), INTENT(out)   ::   tab2d     ! output 2D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      !$acc parallel loop gang vector default(present) private(jid,jjd) async(1)
      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1 ,  size(tab2d,1) ) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / size(tab2d,1)  + 1
         tab2d(jid, jjd) = tab1d( jn)
      END DO
      !$acc end parallel loop
   END SUBROUTINE tab_1d_2d_gpu

   SUBROUTINE tab_1d_2d_device( ndim1d, tab_ind, tab1d, tab2d )
      !$acc routine gang(dim:1) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_1d_2d  ***
      !! restore ji,jj dimensions from a packed form using tab_ind as a map
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d    ! 1D size
      INTEGER, DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind   ! input index
      REAL(dp), DIMENSION(:) ,      INTENT(in   ) ::   tab1d     ! input 1D field
      REAL(dp), DIMENSION(:,:), INTENT(out)   ::   tab2d     ! output 2D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      !$acc loop gang(dim:1) vector private(jid,jjd)
      DO jn = 1, ndim1d
         jid             = MOD( tab_ind(jn) - 1 ,  size(tab2d,1) ) + 1
         jjd             =    ( tab_ind(jn) - 1 ) / size(tab2d,1)  + 1
         tab2d(jid, jjd) = tab1d( jn)
      END DO
      !$acc end loop
   END SUBROUTINE tab_1d_2d_device

   SUBROUTINE packloc_custom(x, threshold, nptidx, npti)
      USE sum_prefix_custom, only: sum_prefix_custom
      REAL(dp), INTENT(in) :: x(:,:)
      INTEGER, INTENT(in) :: threshold
      INTEGER, INTENT(inout) :: nptidx(:)
      INTEGER, INTENT(inout) :: npti 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: scan_idxflags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: scan_idxoffsets
      INTEGER :: scan_idx
      INTEGER :: scan_offset
      INTEGER :: sz, jpi, jpj, ji, jj

      npti = 0 
      sz = size(x)
      jpi = size(x, 1)
      jpj = size(x, 2)

      ALLOCATE(scan_idxflags(sz))
      ALLOCATE(scan_idxoffsets(sz))
      !$acc data create(scan_idxflags,scan_idxoffsets)
      
      !$acc parallel loop collapse(2) private(scan_idx) present(x) async(1)
      DO jj = 1, jpj
         DO ji = 1, jpi
           scan_idx = (jj - 1) * jpi + ji
           IF ( x(ji,jj) > threshold ) THEN
              scan_idxflags(scan_idx) = 1
           else  
              scan_idxflags(scan_idx) = 0
           ENDIF
         END DO
      END DO
      !$acc end parallel loop

      !$acc host_data use_device(scan_idxflags, scan_idxoffsets)
      CALL sum_prefix_custom(scan_idxflags, scan_idxoffsets)
      !$acc end host_data
     
      !$acc parallel loop collapse(2) private(scan_idx,scan_offset) default(present) async(1)
      DO jj = 1, jpj
         DO ji = 1, jpi
            scan_idx = (jj - 1) * jpi + ji
            IF (scan_idxflags(scan_idx) == 1) THEN
              scan_offset = scan_idxoffsets(scan_idx)
              nptidx( scan_offset ) = scan_idx
            ENDIF
         END DO
      END DO
      !$acc end parallel loop
     
      !$acc kernels async(1)
      npti = scan_idxoffsets(sz)
      !$acc end kernels

      !$acc wait(1)  

      !$acc end data
      DEALLOCATE(scan_idxflags)
      DEALLOCATE(scan_idxoffsets)
            
   END SUBROUTINE

!!---------------------------------------------------------------------------------------------
#if 0
   !! Backup versions of sum_prefix_custom 
   SUBROUTINE sum_prefix_custom(input, output)
      INTEGER, INTENT(in) :: input(:)
      INTEGER, INTENT(out) :: output(:)
      INTEGER :: i
  
      !$acc update host(input)
      output(1) = input(1) 
      do i=2,size(input) ; output(i) = output(i - 1) + input(i); enddo 
      !$acc update device(output)
   END SUBROUTINE

   SUBROUTINE sum_prefix_custom(input, output)
      use cutensorex, only: sum_prefix
      INTEGER, INTENT(in) :: input(:)
      INTEGER, INTENT(out) :: output(:)
      !$acc data present(input,output)
      !$acc host_data use_device(input, output) 
      output = SUM_PREFIX(input)
      !$acc end host_data 
      !$acc end data
   END SUBROUTINE
#endif

END PROGRAM


!!#define  TAB_SET_BEGIN(N)   !$acc parallel loop gang(dim:2) num_gangs(256,N) async(1)  \
!!     do p=1,N \
!!      select case(p) 
!!
!!#define TAB_2D_1D(N, a, b)   case(N)  ; call tab_2d_1d_device(npti, nptidx, a, b)
!!
!!#define TAB_SET_END   end select;       end do
!!
!!
!!    TAB_SET_BEGIN(10)
!!         TAB_2D_1D(1, z1, y)
!!         TAB_2D_1D(2, z2, y)
!!
!!
!!      TAB_SET_END

