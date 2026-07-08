!! nvfortran -cuda -acc=gpu -O3 -o packloc_test_fortran packloc_test_fortran.f90 -cudalib=cutensor

!!#define WITH_CHECK

PROGRAM packloc_test
   IMPLICIT NONE
   INTEGER, PARAMETER :: INITIAL_NCYCLES = 10000
   INTEGER :: i = 1, imax = 16*1024*1024
   INTEGER :: ncycles, usec, incr
   ncycles = INITIAL_NCYCLES
   do 
      if (i > imax) exit 
      if (i > 16*1024) ncycles = INITIAL_NCYCLES / 10;
      if (i > 1024*1024) ncycles = INITIAL_NCYCLES / 100;
      usec = test(i)
      if (i < 10 .or. mod(i,10) == 1) then
         write (*,'(A,I0,A,F0.2)') "i=", i, " usec=", real(usec)/real(ncycles)
      endif
      incr = (i/100)*10
      if (incr == 0 .and. i < 11) incr = 1
      if (incr == 0) incr = 10
      i = i + incr 
   end do

CONTAINS

   FUNCTION TEST(sz)
#if defined WITH_PACKLOC_CUTENSOREX
      use cutensorex, only: packloc
#endif
      IMPLICIT NONE
      INTEGER :: sz
      INTEGER :: test
      INTEGER :: i, k, threshold = 6, n, n_cpu
      INTEGER :: count_rate, count_start, count_end
      INTEGER, ALLOCATABLE :: x(:), y(:), y_cpu(:)
      LOGICAL, ALLOCATABLE :: condition(:) 
      ALLOCATE(X(sz), source=0)
      ALLOCATE(Y(sz), source=0)
      ALLOCATE(Y_CPU(sz), source=0)
      ALLOCATE(condition(sz), source=.false.)

      x = [( merge(i, -i, mod(i,2)==1), i=1,sz )]
      do i=1,sz; if (x(i) > 0) x(i) = mod(x(i),10); end do 
      do i=1,sz; condition(i) = (x(i) > threshold); end do

#if defined WITH_CHECK
      call packloc_cpu(condition,y_cpu,n_cpu)
#endif

      !$acc data copy(y,condition)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
#if defined WITH_PACKLOC_CUTENSOREX
         !$acc host_data use_device(condition,y)
         y = packloc(condition,count=n)
         !$acc end host_data
#elif defined WITH_PACKLOC_CUSTOM_OLD
         call packloc_custom_old(condition,y,n)
#elif defined WITH_PACKLOC_CUSTOM_NEW
         call packloc_custom_new(condition,y,n)
#endif
      end do
      call cudaDeviceSynchronize()
      call system_clock(count_end)
      !$acc end data 

      test = INT(real(count_end - count_start) / real(count_rate) * 1e6)

#if defined WITH_CHECK
      if (n /= n_cpu) then
         write(*,'(A,I0,A,I0,A,I0)') ">> i=", sz, " comparison of n failed, n on gpu: ", n, "; on cpu: ", n_cpu
         stop 0 
      end if
      do i=1,n_cpu
         if (y(i) /= y_cpu(i)) then
            write (*,'(A,I0,A,I0,A,I0,A,I0)') ">> i=", sz, " comparison failed for elem=", i, " gpu: ", y(i), " cpu: ", y_cpu(i); 
            stop 0
         end if
      end do
#endif

      DEALLOCATE(X, Y, Y_CPU, condition) 

   END FUNCTION

   SUBROUTINE packloc_cpu(condition, y, n)
       IMPLICIT NONE
       LOGICAL, INTENT(in)  :: condition(:)
       INTEGER, INTENT(out) :: y(:)
       INTEGER, INTENT(out) :: n
       INTEGER              :: sz
       sz = size(condition)
       n = 0
       do i=1,sz
          if (condition(i)) then
             n = n + 1
             y(n) = i
          endif
       end do
   END SUBROUTINE

   SUBROUTINE packloc_custom_old(condition, y, n)
       USE sum_prefix_custom, only: sum_prefix_custom
       IMPLICIT NONE
       LOGICAL, INTENT(in)  :: condition(:)
       INTEGER, INTENT(out) :: y(:)
       INTEGER, INTENT(out) :: n
       INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: scan_idxflags
       INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: scan_idxoffsets
       INTEGER :: sz, i
 
       n = 0
       sz = size(condition)
       ALLOCATE(scan_idxflags(sz),source=0)
       ALLOCATE(scan_idxoffsets(sz),source=0)
       !$acc data create(scan_idxflags,scan_idxoffsets)
 
       !$acc parallel loop default(present) async(1)
       do i=1,sz; if (condition(i)) scan_idxflags=1; enddo
       !$acc end parallel loop  

      !$acc host_data use_device(scan_idxflags, scan_idxoffsets)
      CALL sum_prefix_custom(scan_idxflags, scan_idxoffsets)
      !$acc end host_data

      !$acc parallel loop default(present) async(1)
      do i=1,sz; if (scan_idxflags(i) == 1) y(scan_idxoffsets(i)) = i; enddo
      !$acc end parallel loop

      !$acc kernels async(1)
      n = scan_idxoffsets(sz)
      !$acc end kernels

      !$acc wait(1)  
 
      !$acc end data
      DEALLOCATE(scan_idxflags)
      DEALLOCATE(scan_idxoffsets)
   END SUBROUTINE

   SUBROUTINE packloc_custom_new(condition, y, n)
       USE sum_prefix_custom, only: packloc_custom
       IMPLICIT NONE
       LOGICAL, INTENT(in)  :: condition(:)
       INTEGER, INTENT(out) :: y(:)
       INTEGER, INTENT(out) :: n
       INTEGER :: sz, i
 
       n = 0
       sz = size(condition)
       !$acc host_data use_device(condition, y)
       CALL packloc_custom(condition, y, n)
       !$acc end host_data
    END SUBROUTINE


END PROGRAM
