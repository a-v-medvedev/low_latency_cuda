!! nvfortran -cuda -acc=gpu -O3 -o sc_test sc_test.f90 -cudalib=cutensor

PROGRAM SCTEST
   IMPLICIT NONE
   INTEGER, PARAMETER :: jpi = 999, jpj = 1899
   INTEGER, PARAMETER :: ncycles = 1
   INTEGER :: x_(jpi, jpj), y_(jpi, jpj)
   INTEGER :: threshold

   write (*,*) "-- generate"
   threshold = generate(x_, y_, 6)
   write (*,*) "SUM x: ", SUM(x_)
   write (*,*) "SUM y: ", SUM(y_)

   write (*,*) "-- test1:"
   call test1(x_, y_)
   write (*,*) "-- test2:"
   call test2(x_, y_)
   write (*,*) "-- test3:"
   call test3(x_, y_)

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
      LOGICAL, ALLOCATABLE :: mask(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER :: k, m, count_rate, count_start, count_end
      ALLOCATE(nptidx(jpi*jpj), z(jpi*jpj))
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         npti = 0; nptidx(:) = 0 
         z = 0
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
      write (*,*) "npti: ", npti
      write (*,*) "SUM nptidx: ", SUM(nptidx)
      write (*,*) "SUM z: ", SUM(z)
      write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
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
         if (.not. allocated(z)) allocate(z(npti))
         DO m = 1, 10
            z = pack(reshape(y, [jpi*jpj]), mask)
         END DO
      end do
      call system_clock(count_end)
      write (*,*) "npti: ", npti
      write (*,*) "SUM nptidx: ", SUM(nptidx)
      write (*,*) "SUM z: ", SUM(z)
      write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      deallocate(z)
   END SUBROUTINE

   ! the GPU version of test2 (?) -- expected to call pack() and reshape() from the cutensorex
   SUBROUTINE test3(x, y)
      USE cutensorex, only: pack, reshape, count
      INTEGER :: x(:,:), y(:,:)
      INTEGER :: i, npti 
      INTEGER, ALLOCATABLE :: nptidx(:)
      LOGICAL, ALLOCATABLE :: mask(:)
      INTEGER, ALLOCATABLE :: z(:)
      INTEGER, ALLOCATABLE :: y1d(:)
      INTEGER :: k, m, count_rate, count_start, count_end
      !$acc data copyin(x,y)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc host_data use_device(x)
         mask    = reshape(x > threshold, [jpi*jpj])
         !$acc end host_data  
         !$acc host_data use_device(mask)
         npti    = count(mask)
         !$acc end host_data
         !$acc host_data use_device(mask)
         nptidx  = pack([(i, i=1,jpi*jpj)], mask)
         !$acc end host_data
         if (.not. allocated(z)) allocate(z(npti))
         DO m = 1, 10
            !$acc host_data use_device(y)  
            y1d = reshape(y, [jpi*jpj])
            !$acc end host_data
            !$acc host_data use_device(y1d,mask)
            z = pack(y1d, mask)
            !$acc end host_data
         END DO
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data
      write (*,*) "npti: ", npti
      write (*,*) "SUM nptidx: ", SUM(nptidx)
      write (*,*) "SUM z: ", SUM(z)
      write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      deallocate(z)
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

END PROGRAM



