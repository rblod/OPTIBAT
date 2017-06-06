!>@file   wkbcst.F90
!>@brief Wave model physical constants
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE wkbcst
	
   IMPLICIT NONE
   PUBLIC
	
   INTEGER, PUBLIC, PARAMETER ::   sp = SELECTED_REAL_KIND( 6, 37)   !: single precision (real 4)
   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                           !: working precision
   INTEGER, PUBLIC, PARAMETER ::   lc = 256                          !: Lenght of Character strings   

   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any early
                                            !  output can be collected; do not change
   INTEGER ::   numnam          =   -1      !: logical unit for reference namelist
   LOGICAL ::   lwp      = .TRUE.           !: boolean : true on the 1st processor only .OR. ln_ctl

   REAL(wp), PARAMETER :: g = 9.81
   REAL(wp), PARAMETER :: pi = 3.14159265358979323846
   REAL(wp), PARAMETER :: deg2rad = pi/180.
   REAL(wp), PARAMETER :: sec2day = 1./86400.

   REAL, PARAMETER ::  khmax=20.D0     ! deep-water limit for k*h
        		
END MODULE wkbcst
