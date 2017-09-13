!>@file   wkbbry.F90
!>@brief Wave boundary conditions
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE wkbbry

   USE netcdf
   USE par_wkb
   USE wkbmod
   USE rw_wkb

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC wkb_bry
   
   CONTAINS
   
   !> Update bundary conditions :\
   !> forcing offshore and "internal" conditions. 
   !> We can deal with periodicity or clamped conditions.
   
   SUBROUTINE wkb_bry(kstp)
   	INTEGER, INTENT(in) :: kstp
   	!
   	INTEGER :: jj,kread
   	REAL(wp) :: cff,cdir,wamp,cfrq,khd,kh,dd,hh,kw,cosw,sinw,hrm 
   	      
      IF( MOD (real(kstp,wp),bry_frq/rdt) == 0 ) THEN
         wkxbry_west_dt(:,2) = wkxbry_west_dt(:,1)
         wkebry_west_dt(:,2) = wkebry_west_dt(:,1)
         wacbry_west_dt(:,2) = wacbry_west_dt(:,1)
         hbry_west_dt(:,2) = hbry_west_dt(:,1)
         !
         kread = ((kstp-1)*rdt+bry_frq)/bry_frq +1
         write(*,*) 'read forcing at kstp = ',  kstp 
         write(*,*) 'read record for kread = ', kread 
         write(*,*)
   	   ! bry_tide
   	   ! bry_period
   	   ! bry_hs
   	   ! bry_dir
   	   CALL wkb_read_bry(kread)
   	   !
   	   DO jj=1, jpj 
         dd = h(1,jj) + hbry_west(jj,1)
         cdir=deg2rad*dirbry_west(jj,1)
         wamp=0.5*0.7*hsbry_west(jj,1)
         cfrq=2*pi/perbry_west(jj,1)
         khd = dd*cfrq*cfrq/g
         kh = sqrt(    khd*khd + khd/(1.0 + khd*(0.6666666666 &
         +khd*(0.3555555555 + khd*(0.1608465608 &
         +khd*(0.0632098765 + khd*(0.0217540484 &
         +khd*0.0065407983)))))) )
         kw=kh/dd
         cosw=cos(cdir)
         sinw=sin(cdir)
         hrm=wamp*2.0
         wacbry_west_dt(:,1)=0.125*g*(hrm**2)/cfrq
         wkxbry_west_dt(:,1)=kw*cosw
         wkebry_west_dt(:,1)=kw*sinw
       END DO
   	ENDIF
      cff= MOD (real(kstp),bry_frq/rdt)
      cff=1 
      !
   	wkxbry_west(:) = cff*wkxbry_west_dt(:,1) + (1-cff)*wkxbry_west_dt(:,2)
   	wkebry_west(:) = cff*wkebry_west_dt(:,1) + (1-cff)*wkebry_west_dt(:,2)
   	wacbry_west(:) = cff*wacbry_west_dt(:,1) + (1-cff)*wacbry_west_dt(:,2)
   	hbry_west(:,1)   = cff*hbry_west_dt(:,1)   + (1-cff)*hbry_west_dt(:,2)
   	
   
   END SUBROUTINE wkb_bry

END MODULE wkbbry
