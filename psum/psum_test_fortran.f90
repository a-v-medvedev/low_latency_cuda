!! nvfortran -cuda -acc=gpu -O3 -o psum_fortran_test psum_fortran_test.f90 -cudalib=cutensor

PROGRAM psum_test
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

FUNCTION TEST(N)
   use cutensorex, only: SUM_PREFIX 
   INTEGER :: N
   INTEGER :: test
   INTEGER :: i, k 
   INTEGER :: count_rate, count_start, count_end
   INTEGER, ALLOCATABLE :: x(:), y(:)
   ALLOCATE(X(N), Y(N))
   x = [( merge(i, -i, mod(i,2)==1), i=1,N )]
   !$acc data copy(x,y)
   call system_clock(count_rate = count_rate)
   call system_clock(count_start)
   do k = 1, ncycles
      !!!$acc host_data use_device(x,y)
      y = SUM_PREFIX(x)
   end do
   call cudaDeviceSynchronize()
   call system_clock(count_end)
   !$acc end data 
   test = INT(real(count_end - count_start) / real(count_rate) * 1e6)
END FUNCTION

END PROGRAM
