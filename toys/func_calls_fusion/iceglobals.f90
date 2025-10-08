module iceglobals
implicit none
integer,public :: jpi=100,jpj=100,jpl=100
!$acc declare create(jpi,jpj,jpl)
end module
