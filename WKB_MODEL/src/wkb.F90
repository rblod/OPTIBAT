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
   
   IMPLICIT NONE
   
   INTEGER :: iif
  
   CALL wkb_nam
   
   wstp=1
   CALL wkb_ini

   winfo=1   
   iwave=1
   thwave=1.D+10  
   wnew=1 
   IF ( .NOT. ln_rst ) THEN       ! WKB ray steady mode
     DO WHILE ( iwave .LE. nitermax .AND. thwave .GE. eps ) 
         wstp=wnew
         wnew=wstp+1
         IF (wnew.ge.3) wnew=1
         CALL wkb_stp
         CALL wkb_diag 
         iwave=iwave+1
         thwave=MAX(av_wac,av_wkn)
      END DO 
   ENDIF
   iwave=2
   DO iif= nit000, nitend         ! WKB ray equation time stepping
      wstp=wnew
      wnew=wstp+1
      IF (wnew.ge.3) wnew=1
      CALL wkb_stp 
   END DO 

   !
   CALL wkb_wri

!call get_bry_wkb
!          call set_bry_wkb (tile)

END PROGRAM
