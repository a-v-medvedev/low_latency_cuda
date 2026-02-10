PROGRAM SCTEST
   IMPLICIT NONE
   INTEGER :: ncycles = 1, N = 20240, chksum = 0
   INTEGER, ALLOCATABLE :: input(:), output(:)
   INTEGER :: threshold

   allocate(input(N),output(N))

   !$acc data create(input, output)

   call generate(input,chksum)

   ! warmup
   call test1(input, output, chksum)
   call test2(input, output, chksum)
   call test3(input, output, chksum)
   call cleanup(output)

   ncycles = 30
   write (*,*) "-- test1:"
   call test1(input, output, chksum)
   call cleanup(output)
   write (*,*) "-- test2:"
   call test2(input, output, chksum)
   call cleanup(output)
   write (*,*) "-- test3:"
   call test3(input, output, chksum)
   call cleanup(output)

   !$acc end data

   deallocate(input,output)

CONTAINS 

   SUBROUTINE generate(input, chksum)
      INTEGER, INTENT(inout) :: input(:) 
      INTEGER, INTENT(out)   :: chksum
      INTEGER :: i, total = 0

      !$acc parallel loop gang vector default(present)
      do i=1,N ; input(i)=i; if (mod(i, 2) == 0) input(i) = -input(i); if (mod(i, 97) == 0) input(i) = input(i) * 2; enddo
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present) reduction(+:total)
      do i=1,N ; total = total + input(i); end do
      !$acc end parallel loop
      chksum = total
   END SUBROUTINE

   SUBROUTINE cleanup(output)
      INTEGER, INTENT(inout) :: output(:)
      INTEGER :: i

      !$acc parallel loop gang vector default(present)
      do i=1,N ; output(i)=0; enddo
      !$acc end parallel loop
   END SUBROUTINE

   ! Naive CPU code 
   SUBROUTINE test1(x, y, chksum)
      INTEGER, INTENT(inout) :: x(:), y(:), chksum
      INTEGER :: i, k, count_rate, count_start, count_end
    
      !$acc data present(x,y)
      !$acc update host(x)
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         y(1) = x(1)
         do i=2,size(x) ; y(i) = y(i - 1) + x(i) ; enddo
      enddo
      call system_clock(count_end)
      !$acc update device(y)
      !$acc end data

      if (ncycles /= 1) then
         write (*,*) "output(last): ", y(size(y))
         write (*,*) "test: ", (y(size(y)) == chksum)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
   END SUBROUTINE

   ! SUM_PREFIX GPU code 
   SUBROUTINE test2(x, y, chksum)
#ifdef WITH_CUTENSOREX
      use cutensorex, only: sum_prefix
      INTEGER, INTENT(inout) :: x(:), y(:), chksum
      INTEGER :: i, k, count_rate, count_start, count_end
    
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc data present(x,y)
         !$acc host_data use_device(x,y) 
         y = SUM_PREFIX(x)
         !$acc end host_data 
         !$acc end data
      enddo
      call system_clock(count_end)

      if (ncycles /= 1) then
         write (*,*) "output(last): ", y(size(y))
         write (*,*) "test: ", (y(size(y)) == chksum)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
#else
      write (*,*) ">> NOTE: test2 is not compiled in since it uses CUTENSOREX CUDA library"
#endif
   END SUBROUTINE

   ! "CALL sum_prefix_custom()" GPU code -- fancy combination of SUM_PREFIX, thrust::inclusive_scan and custom C++ kernel
   ! this option is expected to be fastest (or equal to test2) for array sizes more than 10k elements
   SUBROUTINE test3(x, y, chksum)
      USE sum_prefix_custom, only: sum_prefix_custom
      INTEGER, INTENT(inout) :: x(:), y(:), chksum
      INTEGER :: i, k, count_rate, count_start, count_end
    
      call system_clock(count_rate = count_rate)
      call system_clock(count_start)
      do k = 1, ncycles
         !$acc data present(x,y)
         !$acc host_data use_device(x,y) 
         CALL sum_prefix_custom(x, y)
         !$acc end host_data 
         !$acc end data
      enddo
      !$acc wait(1)
      call system_clock(count_end)

      if (ncycles /= 1) then
         write (*,*) "output(last): ", y(size(y))
         write (*,*) "test: ", (y(size(y)) == chksum)
         write (*,*) "time: ", INT(real(count_end - count_start) / real(count_rate) / ncycles * 1e6)
      endif
   END SUBROUTINE

END PROGRAM

