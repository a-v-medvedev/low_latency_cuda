!! nvfortran -cuda -acc=gpu -O3 -o sc_test sc_test.f90 -cudalib=cutensor

PROGRAM SCTEST
   IMPLICIT NONE
   INTEGER, PARAMETER :: jpi = 999, jpj = 1899
   !!INTEGER, PARAMETER :: jpi = 5, jpj = 3
   INTEGER :: ncycles = 1
   INTEGER :: x_(jpi, jpj), y_(jpi, jpj)
   INTEGER :: threshold

   write (*,*) "-- generate"
   threshold = generate(x_, y_, 6)
   write (*,*) "SUM x: ", SUM(x_)
   write (*,*) "SUM y: ", SUM(y_)

   call test1(x_,y_)
   call test2(x_,y_)
   call test3(x_,y_)
   call test4(x_,y_)
   
   ncycles = 2
   write (*,*) "-- test1:"
   call test1(x_, y_)
   write (*,*) "-- test2:"
   !!call test2(x_, y_)
   write (*,*) "-- test3:"
   call test3(x_, y_)
   write (*,*) "-- test4:"
   call test4(x_, y_)
   write (*,*) "-- test5:"
   call test5(x_, y_)


CONTAINS 

   FUNCTION generate(x, y, n)
      INTEGER :: x(:,:), y(:,:)
      INTEGER :: n
      INTEGER :: generate
      INTEGER :: FIB(10) = (/1, 1, 2, 3, 5, 8, 13, 21, 34, 55/)
      INTEGER :: THR(10) = (/10, 10, 30, 50, 1000, 100000, 1000000, 10000000, 1000000000, 1000000000/)
      INTEGER :: ji, jj, a1 = 13, b1 = 8, c1 = 1, a2 = 1, b2 = 1, c2 = 1
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
      generate = THR(n)
   END FUNCTION

   ! original CPU code that creates a map (the nptidx array) and then does the 2d->1d compaction 
   ! according to the map using the tab_2d_1d() subroutine
   SUBROUTINE test1(x, y)
      INTEGER :: x(:,:), y(:,:)
      INTEGER :: npti 
      INTEGER, ALLOCATABLE :: nptidx(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER :: k, m, count_rate, count_start, count_end
      ALLOCATE(nptidx(jpi*jpj), z(jpi*jpj))
      nptidx = 0
      z = 0
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
         DO m = 1, 10
            CALL tab_2d_1d( npti, nptidx, z, y)
         END DO
      end do
      call system_clock(count_end)
      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx)
         write (*,*) "SUM z: ", SUM(z)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
      DEALLOCATE(nptidx, z)
   END SUBROUTINE

   ! the CPU code doing the same that creates a 1D mask array using reshape intrinsic, then uses this mask for the
   ! pack intrinsic to do actual compaction
   SUBROUTINE test2(x, y)
      INTEGER :: x(:,:), y(:,:)
      INTEGER :: i, npti 
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER :: k, m, count_rate, count_start, count_end
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         mask    = reshape(x > threshold, [jpi*jpj])
         npti    = count(mask)
         nptidx  = pack([(i, i=1,jpi*jpj)], mask)
         DO m = 1, 10
            z = pack(reshape(y, [jpi*jpj]), mask)
         END DO
      end do
      call system_clock(count_end)
      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx)
         write (*,*) "SUM z: ", SUM(z)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
   END SUBROUTINE

   ! the GPU version of test2 (?) -- expected to call pack() and reshape() from the cutensorex
   SUBROUTINE test3(x, y)
      INTEGER :: x(:,:), y(:,:)
      INTEGER :: i, npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER :: k, m, count_rate, count_start, count_end

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z(jpi*jpj))       ! largest possible size
      nptidx = 0
      z = 0

      !$acc data copyin(x,y) copyout(nptidx,z)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         call packloc_custom(x, threshold, nptidx, npti)
         DO m = 1, 10
            call tab_2d_1d_gpu(npti, nptidx, z, y)
         END DO
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "SUM z: ", SUM(z(1:npti))
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      deallocate(z)
      deallocate(nptidx)
   END SUBROUTINE
 
   ! the GPU version of test2 (?) -- expected to call pack() and reshape() from the cutensorex
   SUBROUTINE test4(x, y)
      USE cutensorex, only: pack, packloc
      INTEGER :: x(:,:), y(:,:)
      TARGET :: x
      INTEGER :: i, npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER, POINTER :: x1d(:)
      INTEGER :: k, m, count_rate, count_start, count_end

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z(jpi*jpj))       ! largest possible size

      !$acc data copyin(x,y) copyout(nptidx,z)
      x1d(1:jpi*jpj) => x(:,:)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc host_data use_device(x1d, nptidx)
         nptidx  = packloc(x1d > threshold, count=npti)
         !$acc end host_data
         DO m = 1, 10
            !$acc host_data use_device(x, y, z)
            z = pack(y, x > threshold)
            !$acc end host_data
         END DO
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "SUM z: ", SUM(z(1:npti))
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      deallocate(z)
      deallocate(nptidx)
   END SUBROUTINE
 
   ! the GPU version of test2 (?) -- expected to call pack() and reshape() from the cutensorex
   SUBROUTINE test5(x, y)
      USE cutensorex, only: pack, packloc
      INTEGER :: x(:,:), y(:,:)
      TARGET :: x
      INTEGER :: i, npti
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER, POINTER :: x1d(:)
      INTEGER :: k, m, count_rate, count_start, count_end

      allocate(nptidx(jpi*jpj))  ! largest possible size
      allocate(z(jpi*jpj))       ! largest possible size

      !$acc data copyin(x,y) copyout(nptidx,z)
      x1d(1:jpi*jpj) => x(:,:)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc host_data use_device(x1d, nptidx)
         nptidx  = packloc(x1d > threshold, count=npti)
         !$acc end host_data
         DO m = 1, 10
            call tab_2d_1d_gpu(npti, nptidx, z, y)
         END DO
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "npti: ", npti
         write (*,*) "SUM nptidx: ", SUM(nptidx(1:npti))
         write (*,*) "SUM z: ", SUM(z(1:npti))
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif 
      deallocate(z)
      deallocate(nptidx)
   END SUBROUTINE
 
   SUBROUTINE tab_2d_1d( ndim1d, tab_ind, tab1d, tab2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d   ! 1d size
      INTEGER , DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind  ! input index
      INTEGER, DIMENSION(jpi,jpj), INTENT(in   ) ::   tab2d    ! input 2D field
      INTEGER, DIMENSION(ndim1d) , INTENT(out) ::   tab1d    ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , jpi ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / jpi + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO
   END SUBROUTINE tab_2d_1d

   SUBROUTINE tab_2d_1d_gpu(ndim1d, tab_ind, tab1d, tab2d )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d   ! 1d size
      INTEGER , DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind  ! input index
      INTEGER, DIMENSION(jpi,jpj), INTENT(in   ) ::   tab2d    ! input 2D field
      INTEGER, DIMENSION(ndim1d) , INTENT(out) ::   tab1d    ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      !$acc data present(tab1d,tab2d,tab_ind)
      !$acc parallel loop private(jid,jjd)
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , jpi ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / jpi + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO
      !$acc end parallel loop
      !$acc end data
   END SUBROUTINE tab_2d_1d_gpu

   SUBROUTINE packloc_custom(x, threshold, nptidx, npti)
      INTEGER, INTENT(in) :: x(:,:)
      INTEGER, INTENT(in) :: threshold
      INTEGER, INTENT(inout) :: nptidx(:)
      INTEGER, INTENT(inout) :: npti 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: scan_idxflags
      INTEGER, ALLOCATABLE, DIMENSION(:) :: scan_idxoffsets
      INTEGER :: scan_idx
      INTEGER :: scan_offset
      INTEGER :: sz, jpi, jpj, ji, jj
      LOGICAL :: nptidxdone = .false.

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

      !!$acc host_data use_device(scan_idxflags, scan_idxoffsets)
      !scan_idxoffsets = SUM_PREFIX(scan_idxflags) 
      !!$acc end host_data
      CALL sum_prefix_custom(scan_idxflags, scan_idxoffsets, nptidx, nptidxdone)
     
      if (.not. nptidxdone) then
      !$acc parallel loop collapse(2) private(scan_idx,scan_offset) present(nptidx) async(1)
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
      endif
     
      !$acc kernels async(1)
      npti = scan_idxoffsets(sz)
      !$acc end kernels

      !$acc end data
      DEALLOCATE(scan_idxflags)
      DEALLOCATE(scan_idxoffsets)
            
   END SUBROUTINE

!!   SUBROUTINE sum_prefix_custom(input, output, nptidx, nptidxdone)
!!      use cutensorex, only: sum_prefix
!!      INTEGER, INTENT(in) :: input(:)
!!      INTEGER, INTENT(out) :: output(:)
!!      INTEGER, ALLOCATABLE :: y(:)
!!      !$acc data present(input,output)
!!      !!!$acc host_data use_device(input, y) 
!!      output = SUM_PREFIX(input)
!!      !!!$acc end host_data 
!!      !$acc end data
!!   END SUBROUTINE

END PROGRAM

