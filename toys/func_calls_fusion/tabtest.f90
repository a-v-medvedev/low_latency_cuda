program tabtest
use iceglobals
use icetab
implicit none

integer, parameter         :: npti = 10
integer, dimension(npti)   :: map
real, dimension(npti)      :: t1d_1, t1d_2, t1d_3, t1d_4, t1d_5, t1d_6, t1d_7, t1d_8, t1d_9
real, dimension(jpi,jpj)   :: t2d_1, t2d_2, t2d_3, t2d_4, t2d_5, t2d_6, t2d_7, t2d_8, t2d_9
integer                    :: i 

t1d_1 = 0; t2d_1 = 1
t1d_2 = 0; t2d_2 = 1
t1d_3 = 0; t2d_3 = 1
t1d_4 = 0; t2d_4 = 1
t1d_5 = 0; t2d_5 = 1
t1d_6 = 0; t2d_6 = 1
t1d_7 = 0; t2d_7 = 1
t1d_8 = 0; t2d_8 = 1
t1d_9 = 0; t2d_9 = 1

t2d_1 = reshape([( merge(i, -i, mod(i,2)==1), i=1,jpi*jpj )], (/jpi, jpj/))

map = (/0, 2, 3, 5, 8, 13, 21, 34, 55, 89/)

!$acc data copy(map, t1d_1, t2d_1,  t1d_2, t2d_2, t1d_3, t2d_3) async(1)
!$acc data copy(     t1d_4, t2d_4,  t1d_5, t2d_5, t1d_6, t2d_6) async(1)
!$acc data copy(     t1d_7, t2d_7,  t1d_8, t2d_8, t1d_9, t2d_9) async(1)

!$acc kernels async(1)
call tab_2d_1d_gpu(10, map, t1d_1, t2d_1)
call tab_2d_1d_gpu(10, map, t1d_2, t2d_2)
call tab_2d_1d_gpu(10, map, t1d_3, t2d_3)
!$acc end kernels

!$acc end data
!$acc end data
!$acc end data
!$acc wait(1)

print *,t1d_1
print *,t1d_2
print *,t1d_3

end program

