!>@file   wkbini.F90
!>@brief Wave model initialisation
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE wkbini

   USE par_wkb
   USE wkbmod
   USE wkbutil
   USE wkbcdf
   USE rw_wkb
   USE wkbbry
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC wkb_nam, wkb_ini
   
   CONTAINS
   
   SUBROUTINE wkb_nam
   
      INTEGER ::   ios
 
      NAMELIST/namwkb/ ln_rst, nit000, nitend, nstock,                     &
                       jpi, jpj, rdx, rdt, nitermax, eps,                  &
                       cn_dirin, cn_dirout, cn_filein, cn_fileout,         &
                       ln_perio, ln_anabry,                                &
                       ln_brywest, ln_bryeast, ln_brysouth, ln_brynorth,   &
                       wkb_amp, wkb_prd, wkb_dir,                          &  
                       wkb_roller, wkb_gam, wkb_btg, wkb_tide, wkb_rsb,    &
                       cn_bryin, bry_frq

      CALL ctl_opn( numnam, 'namelist', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. )
      REWIND( numnam )              ! Namelist namctl in reference namelist : Control prints & Benchmark
      READ  ( numnam, namwkb, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namwkb in namelist', .TRUE. )
   
      IF(lwp) THEN                            ! open listing units
         !
         CALL ctl_opn( numout, 'wkb.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
         !
         WRITE(numout,*)
         WRITE(numout,*) '                        CNRS - INP - IRD    '
         WRITE(numout,*) '                        The OptiBat team   '
         WRITE(numout,*) '                     version 0.0001  (2016) '
         WRITE(numout,*)
         WRITE(numout,*)
      ENDIF
         !
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'wkb_nam: Control prints'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namwkb'
         WRITE(numout,*) '      X-dimension                 jpi          = ', jpi
         WRITE(numout,*) '      Y-dimension                 jpj          = ', jpj
         WRITE(numout,*) '      Restart                     ln_rst       = ', ln_rst
         WRITE(numout,*) '      Frist time steps            nit000       = ', nit000
         WRITE(numout,*) '      Last of time steps          nitend       = ', nitend
         WRITE(numout,*) '      Output frequency            nstock       = ', nstock
         WRITE(numout,*) '      Resolution                  rdx          = ', rdt
         WRITE(numout,*) '      Time step                   rdt          = ', rdt
         WRITE(numout,*) '      Number of iteration         nitermax     = ', nitermax
         WRITE(numout,*) '      Input directory             cn_dirin     = ', TRIM(cn_dirin)
         WRITE(numout,*) '      Output directory            cn_dirout    = ', TRIM(cn_dirout)
         WRITE(numout,*) '      Input file                  cn_filein    = ', TRIM(cn_filein)
         WRITE(numout,*) '      Output file                 cn_fileout   = ', TRIM(cn_fileout)
         WRITE(numout,*) '      Periodicity (0/1)           ln_perio     = ', ln_perio
         WRITE(numout,*) '      Analytical forcing (0/1)    ln_anabry    = ', ln_anabry
         IF( ln_anabry) THEN
            WRITE(numout,*) '      Wave amplitude              wkb_amp      = ', wkb_amp
            WRITE(numout,*) '      Wave period                 wkb_prd      = ', wkb_prd
            WRITE(numout,*) '      Wave direction              wkb_dir      = ', wkb_dir
         ELSE   
            WRITE(numout,*) '      Forcing file                cn_bryin     = ', cn_bryin
            WRITE(numout,*) '      Forcing freq                bry_frq      = ', bry_frq
         ENDIF
      ENDIF   
   
   END SUBROUTINE wkb_nam

   SUBROUTINE wkb_ini
      INTEGER :: ierr
      !!!!!!!!!!
      INTEGER :: istr,iend,jstr,jend, i, j, k, &
     &          istrR,iendR,jstrR,jendR  
      REAL(wp) ::  sbc, kh,  nw, &
     &     cosw, sinw, cw, peg,  &
     &     inv_k, khn, cff1, cff2, &
     &     c_roller, wkbbry(3)
      REAL(wp) :: wamp, wh, cfrq, cdir, Btg, gamw, khd, kw, kr, ks, &
     &     ho, dd, co, cgo, dsup, cfrc,wkb_rsb, zeps
      REAL(wp), DIMENSION(jpi,jpj) :: Dstp !, zeta
            
      !!!!!!!!!!
      istr=2 ; iend=jpi-1 ; jstr=2 ; jend=jpj-1
      istrR=1 ; iendR=jpi ; jstrR=1 ; jendR=jpj
      wkb_rsb=0.1
      !!!!!!!!!!
      zeps=1.e-10

      ierr = wkb_alloc() 
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'wkb_alloc : unable to allocate standard wave arrays' )

      CALL wkb_read_ini
      
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'wkb_ini: Read Bathymetry'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) 'Max depth = ', MAXVAL(h)
      ENDIF 
      
      IF(.not. ln_anabry) THEN
         CALL wkb_bry(0)
         zeta(:,:) = MAXVAL(hbry_west)
      ELSE
         zeta(:,:) = wkb_tide
      ENDIF
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'wkb_ini: Water level'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) 'Max level = ', MAXVAL(zeta)
      ENDIF 
      
      DO j=jstrR,jendR
         DO i=istrR,iendR
            IF (zeta(i,j) .lt. Dcrit-h(i,j)) THEN
              zeta(i,j)=Dcrit-h(i,j)
            ENDIF
         END DO
      END DO      

      sbc=3.0/16.0*sqrt(pi)*g*(wkb_btg**3)/(2.0*pi)
      c_roller=g*wkb_rsb  ! for roller dissipation term, epsilon_r
      cfrc = 0.5/sqrt(pi)
      peg=8.0/g

!      IF( .NOT. ln_rst ) THEN
         IF (ln_brywest) THEN
            ho=h(1,1) + zeta (1,1) !+ wkb_tide         !+wkb_tide       ! offshore depth
         ELSE IF (ln_bryeast) THEN
            ho=h(jpi,1) + zeta(jpi,1)! + wkb_tide      !+wkb_tide       ! offshore depth
         ELSE IF (ln_brynorth) THEN
            ho=h(jpi,jpj)+ zeta(jpi,jpj)! +wkb_tide  !+wkb_tide       ! offshore depth
         ELSE IF (ln_brysouth) THEN
            ho=h(1,1)+ zeta (1,1)!+wkb_tide          !+wkb_tide       ! offshore depth
         ENDIF
   
   
         IF( ln_anabry ) THEN
            wamp = wkb_amp*0.7/2.        ! wave amplitude (m)
            cfrq = 2.0*pi/wkb_prd        ! peak wave frequency (rad/s)
            cdir = wkb_dir*deg2rad       ! wave direction rad
         ELSE 
            CALL wkb_bry(0)
            IF ( ln_brywest ) THEN
               wkbbry(1)=wacbry_west(1) 
               wkbbry(2)=wkxbry_west(1)
               wkbbry(3)=wkebry_west(1)
            ELSE IF ( ln_bryeast ) THEN
               wkbbry(1)=wacbry_east(1) 
               wkbbry(2)=wkxbry_east(1)
               wkbbry(3)=wkebry_east(1)
            ENDIF
            !MPI
            cff1=sqrt(wkbbry(2)**2+wkbbry(3)**2)
            inv_k =1.0/cff1
            khn =cff1*ho
            cfrq =sqrt(g*cff1*tanh(khn))
            cosw = inv_k*wkbbry(2)          
            sinw = inv_k*wkbbry(3)           
            cdir=ATAN2(sinw,cosw)
            wamp=sqrt(2*wkbbry(1)/g*cfrq)
         ENDIF   ! ln_anabry     
         Btg  = wkb_btg               ! B parameter
         gamw = wkb_gam               ! gamma paramemer (Hrms/h ratio)

         khd = ho*cfrq*cfrq/g
         kh  = sqrt(    khd*khd + khd/(1.0 + khd*(0.6666666666   &
        &                 +khd*(0.3555555555 + khd*(0.1608465608   &
        &                 +khd*(0.0632098765 + khd*(0.0217540484   &
        &                               +khd*0.0065407983)))))) )
         co=sqrt(g/kh*ho*tanh(kh))
         cgo=co*0.5*(1.+2.*kh/sinh(2.*kh))
         DO j=jstrR,jendR
            DO i=istrR,iendR
               Dstp(i,j)=h(i,j)+zeta(i,j)
               dsup=0.0
               DO k=1,4
                  dd = Dstp(i,j) + dsup         !+ wkb_tide
                  khd = dd*cfrq*cfrq/g
                  kh = sqrt(    khd*khd + khd/(1.0 + khd*(0.6666666666  &
        &                      +khd*(0.3555555555 + khd*(0.1608465608   &
        &                      +khd*(0.0632098765 + khd*(0.0217540484   &
        &                                      +khd*0.0065407983)))))) )
                  kw=kh/max(dd,zeps)
                  cw=sqrt(g/kw*tanh(kh))
                  nw=0.5*(1.+2.*kh/sinh(2.*kh))  ! n=Cg/C
                  ks=sqrt(cgo/(2.*nw*cw))        ! shoaling coefficient
                  kr=sqrt(cos(cdir)/cos(cdir))   ! refraction coefficient
                  cosw=cos(cdir)
                  sinw=sin(cdir)
                  cff1=gamw*dd
                  cff2=2.0*wamp*ks*kr
                  wh = min(cff1,cff2)
                  dsup=-0.125*(wh**2)*kw/sinh(2.*kw*dd) ! wave set-up
               END DO
               DO k=1,2
                 hrm(i,j,k)=wamp*2.0
                 frq(i,j,k)=cfrq
                 wcg(i,j,k)=nw*cw/kw
                 wac(i,j,k)=0.125*g*(hrm(i,j,k)**2)/frq(i,j,k)
                 wkx(i,j,k)=kw*cosw
                 wke(i,j,k)=kw*sinw
                 wsb(i,j,k)=sbc/(Dstp(i,j)**5)*(hrm(i,j,k)**7)
                 wvn(i,j,k)=kw
                 war(i,j,k)=0.0
                 wsr(i,j,k)=0.0
                 wcr(i,j,k)=cw/kw
                 wfc(i,j,k)=0.0
               END DO
               wdrx(i,j)=cosw      ! cosine wave direction (xi)
               wdre(i,j)=sinw      ! sine wave direction (eta)
               sup(i,j)=dsup       ! wave setup
               wdsp(i,j)=0.0       ! sbc/(h(i,j)**5)*(hrm(i,j,k)**7)/cfrq 
                                   ! S (\ep_b/rho, wave dissipation)
            END DO
         END DO
 !     ELSE      !ln_rst
!         DO j=jstrR,jendR
!            DO i=istrR,iendR
!               kk=max(sqrt(wkx(i,j,wstp)**2+wke(i,j,wstp)**2),zeps)
!               kh=min(kk*h(i,j),khmax)
!               nw=0.5*(1.0+2.0*kh/max(sinh(2.0*kh),zeps))
!               wvn(i,j,wstp)=kk
!               frq(i,j,wstp)=sqrt(g*kk*tanh(kh))
!!               hrm(i,j,wstp)=sqrt(peg*max(wac(i,j,wstp),0.0)*frq(i,j,wstp))
!               wcr(i,j,wstp)=frq(i,j,wstp)/max(kk**2,zeps)  ! c/k
!               wsr(i,j,wstp)=c_roller*war(i,j,wstp)*kk   &  ! ep_r/rho/sigma
!     &                      /max(frq(i,j,wstp)**2,zeps)
!            END DO
!         END DO

!      ENDIF   ! ln_rst

     IF ( ln_anabry ) THEN
         IF( ln_bryeast ) THEN
            DO j=jstr-1,jend+1
               wacbry_east(j)=wac(iendR,j,1)
               wkxbry_east(j)=wkx(iendR,j,1)
               wkebry_east(j)=wke(iendR,j,1)
            END DO
            !
         ELSE IF( ln_brywest ) THEN
            DO j=jstr-1,jend+1
               wacbry_west(j)=wac(istr-1,j,1)
               wkxbry_west(j)=wkx(istr-1,j,1)
               wkebry_west(j)=wke(istr-1,j,1)
            END DO
         ELSE IF( ln_brysouth ) THEN
            DO i=istr-1,iend+1
               wacbry_south(j)=wac(i,jstr-1,1)
               wkxbry_south(j)=wkx(i,jstr-1,1)
               wkebry_south(j)=wke(i,jstr-1,1)
            END DO
         ELSE IF( ln_brynorth ) THEN
            DO i=istr-1,iend+1
                wacbry_north(i)=wac(i,jendR,1)
                wkxbry_north(i)=wkx(i,jendR,1)
                wkebry_north(i)=wke(i,jendR,1)
            END DO
         ENDIF   ! bry east
      ENDIF   ! anabry  
                     
   END SUBROUTINE wkb_ini

END MODULE wkbini
