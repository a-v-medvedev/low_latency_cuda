MODULE icetab
 use iceglobals
 implicit none
 PUBLIC   tab_1d_2d_gpu, tab_2d_1d_gpu
CONTAINS
 SUBROUTINE tab_1d_2d_gpu( ndim1d, tab_ind, tab1d, tab2d )
 !$acc routine worker
 !!----------------------------------------------------------------------
 !!                  ***  ROUTINE tab_2d_1d  ***
 !!----------------------------------------------------------------------
 INTEGER                     , INTENT(in) ::   ndim1d    ! 1D size
 INTEGER , DIMENSION(ndim1d) , INTENT(in) ::   tab_ind   ! input index
 REAL, DIMENSION(ndim1d) , INTENT(in)     ::   tab1d     ! input 1D field
 REAL, DIMENSION(jpi,jpj), INTENT(inout)  ::   tab2d     ! output 2D field
 !
 INTEGER ::   jn , jid, jjd
 !!----------------------------------------------------------------------
 if(ndim1d<=0) then; return; ENDIF
 !$acc loop private(jn)
 DO jn = 1, ndim1d
    jid             = MOD( tab_ind(jn) - 1 ,  jpi ) + 1
    jjd             =    ( tab_ind(jn) - 1 ) / jpi  + 1
    tab2d(jid, jjd) = tab1d( jn)
 END DO
 !$acc end loop
 END SUBROUTINE tab_1d_2d_gpu

   SUBROUTINE tab_2d_1d_gpu( ndim1d, tab_ind, tab1d, tab2d )
   !$acc routine worker
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tab_2d_1d  ***
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ndim1d   ! 1d size
      INTEGER , DIMENSION(ndim1d) , INTENT(in   ) ::   tab_ind  ! input index
      REAL, DIMENSION(jpi,jpj), INTENT(in   ) ::   tab2d    ! input 2D field
      REAL, DIMENSION(ndim1d) , INTENT(inout) ::   tab1d    ! output 1D field
      !
      INTEGER ::   jn , jid, jjd
      !!----------------------------------------------------------------------
      !$acc loop private(jn)
      DO jn = 1, ndim1d
         jid        = MOD( tab_ind(jn) - 1 , jpi ) + 1
         jjd        =    ( tab_ind(jn) - 1 ) / jpi + 1
         tab1d( jn) = tab2d( jid, jjd)
      END DO
      !$acc end loop
   END SUBROUTINE tab_2d_1d_gpu

end module

