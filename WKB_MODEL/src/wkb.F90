!>@file   wkb.F90
!>@brief Wave model main program
!>@author R. Benshila (CNRS)
!>@version 0.0001

!> Main program calling :
!> - initialisation
!> - forcing
!> - steady state equilibrium if needed
!> - time integrtion
!> - outputs 

PROGRAM wkb

   USE par_wkb
   USE wkbutil
   USE wkbini
   USE wkbstp
   USE wkbmod
   USE wkbcdf
   USE rw_wkb
   USE wkbbry
   
   IMPLICIT NONE
     
   CALL wkb_nam
   
   wstp=1
   CALL wkb_ini

   winfo=1   
   iwave=1
   thwave=1.D+10  
   wnew=1 
!   IF ( .NOT. ln_rst ) THEN       ! WKB ray steady mode
     DO WHILE ( iwave .LE. nitermax .AND. thwave .GE. eps ) 
         wstp=wnew
         wnew=wstp+1
         IF (wnew.ge.3) wnew=1
         CALL wkb_stp
         CALL wkb_diag 
         iwave=iwave+1
         thwave=MAX(av_wac,av_wkn)
      END DO 
!  ENDIF
   iwave=2
      !   CALL wkb_wri
   DO nbstp= nit000, nitend         ! WKB ray equation time stepping
      IF( .NOT. ln_anabry) CALL wkb_bry(nbstp)
      wstp=wnew
      wnew=wstp+1
      IF (wnew.ge.3) wnew=1
      CALL wkb_stp 
      IF( MOD(nbstp,nstock) == 0) THEN
         WRITE(*,*) 'write output at kstp = ', nbstp
         CALL wkb_wri
      ENDIF
   END DO 

   !

!call get_bry_wkb
!          call set_bry_wkb (tile)

END PROGRAM
