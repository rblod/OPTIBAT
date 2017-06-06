!>@file   wkbstp.F90
!>@brief Wave model time integration
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE wkbstp

   USE par_wkb
   USE wkbmod
   USE wkbbry
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC wkb_stp, wkb_diag
   
CONTAINS

   !> Time integration of the follwoing equations : \n
   !> \f$
!>\frac{\partial \mathcal{A}}{\partial t} + \nabla \cdot (\bf{c_g}\mathcal{A})=-\frac{\varepsilon_\textit{b}+\varepsilon_\textit{d} }{\sigma}
!>\f$
!> 
!>\f$
!>\frac{\partial \bf{k}}{\partial t} + {\bf{c_g}} \cdot \nabla {\bf{k}} =-\frac{k \sigma}{\sinh2kD}\nabla D
!>\f$
      
   SUBROUTINE wkb_stp
   
         !!!!!!!!!!
      INTEGER ::  &
     &          istrR,iendR,jstrR,jendR  
      integer istr,iend,jstr,jend, i, j, wprv, wbak
      real(wp), dimension(jpi,jpj) :: cgx,cge,acx,ace, &
     &                                             Dstp,rkx,rke,rwac,rwar,zeta
      real(wp) ::sbc, kh, nw, a0, ax1, ax2, ay1, ay2, &
     &     cxu, cev, peg,  &
     &     inv_k, khn, cff, cff1, w1, w2,  &
     &     dtwave, cfrc, inv_f, urms, abot, fw, c_roller, wkb_rsb, zob !!!!
      real(wp), parameter :: zeps=1.e-10 
            
      !!!!!!!!!!
      istr=2 ; iend=jpi-1 ; jstr=2 ; jend=jpj-1
      wkb_rsb=0.1 ; zob=0.01
      istrR=1 ; iendR=jpi ; jstrR=1 ; jendR=jpj
      !!!!!!!!!!
             cgx=0.;cge=0.;acx=0.;ace=0. ;rwar=0.        
           Dstp=0.;rkx=0.;rke=0.;rwac=0.

      dtwave=rdt
      wprv=wstp
      if (iwave.eq.1) then              ! Option 1: Adams-Bashforth2
        wbak=wstp                       ! The first wave step is solved
        w1=1.                           ! with forward Euler scheme.
        w2=0.
      else
        wbak=wstp-1
        if (wbak.lt.1) wbak=wstp
        w1=1.5
        w2=-0.5
      endif
 !
! Initializing parameters for breaking and roller models
!
!      if (iwave.eq.1) then

#undef WAVE_BREAK_TG86
# ifdef WAVE_BREAK_TG86
      sbc=3.0/16.0*sqrt(pi)*g*(wkb_btg**3)/(wkb_gam**4)/(2.0*pi)
# elif defined WAVE_BREAK_TG86A
      sbc=3.0/16.0*sqrt(pi)*g*(wkb_btg**3)/(wkb_gam**2)/(2.0*pi)
# elif defined WAVE_BREAK_CT93
      sbc=3.0/16.0*sqrt(pi)*g*(wkb_btg**3)/(2.0*pi)
# elif defined WAVE_BREAK_R93
      sbc=0.25*g*(wkb_btg**3)/(2.0*pi)
# else
      sbc=3.0/16.0*sqrt(pi)*g*(wkb_btg**3)/(2.0*pi)
# endif
      c_roller=g*wkb_rsb  ! for roller dissipation term, epsilon_r
      cfrc = 0.5/sqrt(pi)
      peg=8.0/g
!      endif 

!
!
! Prepair for multiple-time step algorithms
!
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          rkx(i,j) =w1*wkx(i,j,wstp)+w2*wkx(i,j,wbak)
          rke(i,j) =w1*wke(i,j,wstp)+w2*wke(i,j,wbak)
          rwac(i,j)=w1*wac(i,j,wstp)+w2*wac(i,j,wbak)
          rwar(i,j)=w1*war(i,j,wstp)+w2*war(i,j,wbak)
        enddo
      enddo

      zeta(:,:)=0.
      DO j=jstrR,jendR
         DO i=istrR,iendR
            IF (zeta(i,j) .lt. Dcrit-h(i,j)) THEN
              zeta(i,j)=Dcrit-h(i,j)
            ENDIF
         END DO
      END DO      
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          Dstp(i,j)=h(i,j)+zeta(i,j)
          cgx(i,j)=wcg(i,j,wstp)*rkx(i,j)
          cge(i,j)=wcg(i,j,wstp)*rke(i,j)
        enddo
      enddo

!
!---------------------------------------------------------------------
! Solve primary wavenumber equation
!---------------------------------------------------------------------
!
      cff=0.5
      do j=jstr,jend
        do i=istr,iend
          cxu=cgx(i,j)
          cev=cge(i,j)
          ax1=cxu+abs(cxu)                ! 1st order upwind scheme
          ax2=cxu-abs(cxu)
          ay1=cev+abs(cev)
          ay2=cev-abs(cev)
          kh =min(wvn(i,j,wstp)*Dstp(i,j),khmax)
          a0 =frq(i,j,wstp)*wvn(i,j,wstp)/max(sinh(2.0*kh),zeps)   
          wkx(i,j,wnew) =wkx(i,j,wprv) -dtwave   *1./rdx                        & 
     &                *( a0*0.5*(Dstp(i+1,j)-Dstp(i-1,j))    &
     &               +0.5*( ax1*(rkx(i  ,j)-rkx(i-1,j))    &
     &                     +ax2*(rkx(i+1,j)-rkx(i  ,j))    &
     &                     +ay1*(rkx(i,j  )-rkx(i,j-1))    &
     &                     +ay2*(rkx(i,j+1)-rkx(i,j  )) )  &
     &                                                            )
          wke(i,j,wnew) =wke(i,j,wprv) -dtwave  *1./rdx                      &
     &                *( a0*0.5*(Dstp(i,j+1)-Dstp(i,j-1))   &
     &               +0.5*( ax1*(rke(i  ,j)-rke(i-1,j))   &
     &                     +ax2*(rke(i+1,j)-rke(i  ,j))   &
     &                     +ay1*(rke(i,j  )-rke(i,j-1))   &
     &                     +ay2*(rke(i,j+1)-rke(i,j  )) ) &
     &                                                            )
        enddo
      enddo               
      wkx(:,1,wnew)=wkx(:,jpj-1,wnew)
      wkx(:,jpj,wnew)=wkx(:,2,wnew)
      wkx(jpi,:,wnew)=wkx(jpi-1,:,wnew)
      wke(:,1,wnew)=wke(:,jpj-1,wnew)
      wke(:,jpj,wnew)=wke(:,2,wnew)
      wke(jpi,:,wnew)=wke(jpi-1,:,wnew)
 !
!---------------------------------------------------------------------
! Solve primary action balance equation
!---------------------------------------------------------------------
    !
      do j=jstr-1,jend                    ! 1st order upwind scheme
        do i=istr-1,iend
          cxu=0.5*(cgx(i+1,j)+cgx(i,j))
          cev=0.5*(cge(i,j+1)+cge(i,j))
          acx(i,j)=0.5*( (cxu-abs(cxu))*rwac(i+1,j)   &
     &                              +(cxu+abs(cxu))*rwac(i,j) )
          ace(i,j)=0.5*( (cev-abs(cev))*rwac(i,j+1)   &
     &                              +(cev+abs(cev))*rwac(i,j) )
        enddo
      enddo           ! <--- discard cgx,cge,rwac
!
      do j=jstr,jend
        do i=istr,iend
          wac(i,j,wnew) =wac(i,j,wprv)-dtwave*(1./rdx  &
     &                 *(acx(i,j)-acx(i-1,j)+ace(i,j)-ace(i,j-1) &
     &                          ) +wsb(i,j,wstp)+wfc(i,j,wstp) )
        enddo
      enddo           ! <--- discard acx,ace

      wac(:,1,wnew)=wac(:,jpj-1,wnew)
      wac(:,jpj,wnew)=wac(:,2,wnew)
      wac(jpi,:,wnew)=wac(jpi-1,:,wnew)
      wac(1,:,wnew)=wacbry_west(:)
      
!----------------------------------------------------------------------
! Solve roller action balance equation
!---------------------------------------------------------------------
!
      do j=jstr-1,jend+1                  ! phase velocity (not Cg)
        do i=istr-1,iend+1
          cgx(i,j)=wcr(i,j,wstp)*rkx(i,j)
          cge(i,j)=wcr(i,j,wstp)*rke(i,j)
        enddo
      enddo           ! <--- discard rkx,rke

      do j=jstr-1,jend                    ! 1st order upwind scheme
        do i=istr-1,iend
          cxu=0.5*(cgx(i+1,j)+cgx(i,j))
          cev=0.5*(cge(i,j+1)+cge(i,j))
          acx(i,j)=0.5*( (cxu-abs(cxu))*rwar(i+1,j)   &
     &                              +(cxu+abs(cxu))*rwar(i,j) )
          ace(i,j)=0.5*( (cev-abs(cev))*rwar(i,j+1)   &
     &                              +(cev+abs(cev))*rwar(i,j) )
        enddo
      enddo           ! <--- discard cgx,cge,rwar

      do j=jstr,jend
        do i=istr,iend
          war(i,j,wnew) =war(i,j,wprv)-dtwave*(1./rdx   &
     &               *(                              &
     &                   acx(i,j)-acx(i-1,j)+ace(i,j)-ace(i,j-1)  &
                      ) &
     &                -wkb_roller*wsb(i,j,wstp) )!+wsr(i,j,wstp) )
        enddo
      enddo           ! <--- discard acx,ace
      war(:,1,wnew)=war(:,jpj-1,wnew)
      war(:,jpj,wnew)=war(:,2,wnew)
      war(jpi,:,wnew)=war(jpi-1,:,wnew)
      war(1,:,wnew)=0.

!---------------------------------------------------------------------
! Set boundary conditions for the WKB wave quantities
! For now, only gradient and clamped BCs are available.
!---------------------------------------------------------------------
!
      call wkb_bry
!
!---------------------------------------------------------------------
! Estimate wave-associated variables for MRL_WCI/BBL routines
!---------------------------------------------------------------------
!
      do j=jstrR,jendR
        do i=istrR,iendR
          wvn(i,j,wnew)=max(sqrt(wkx(i,j,wnew)**2+wke(i,j,wnew)**2),zeps)
          inv_k =1.0/wvn(i,j,wnew)
          khn =min(wvn(i,j,wnew)*Dstp(i,j),khmax)
          nw =0.5*(1.0+2.0*khn/max(sinh(2.0*khn),zeps))
          frq(i,j,wnew) =sqrt(g*wvn(i,j,wnew)*tanh(khn))
          wcg(i,j,wnew) =frq(i,j,wnew)*(inv_k**2)*nw  ! cg/k
          hrm(i,j,wnew) =sqrt(peg*max(wac(i,j,wnew),0.0)*frq(i,j,wnew))
#undef WAVE_BREAK_TG86
# ifdef WAVE_BREAK_TG86
          wsb(i,j,wnew) =sbc*((1.0/Dstp(i,j))**5)*(hrm(i,j,wnew)**7) ! ep_b/rho/sigma
# elif defined WAVE_BREAK_TG86A
          cff1=hrm(i,j,wnew)/(wkb_gam*Dstp(i,j))
          wsb(i,j,wnew)=sbc/(Dstp(i,j)**3)*(hrm(i,j,wnew)**5) &
     &                            *(1.0-(1.0+peg**2)**(-2.5))
# elif defined WAVE_BREAK_CT93
          cff1=hrm(i,j,wnew)/(wkb_gam*Dstp(i,j))
          wsb(i,j,wnew) =sbc/Dstp(i,j)*(hrm(i,j,wnew)**3)
     &      *( 1.0+tanh(8.0*(cff1-1.0)) )*( 1.0-(1.0+cff1**2)**(-2.5) )

#  elif defined WAVE_BREAK_R93
            cff1=hrm(i,j,wnew)/(wkb_gam*Dstp(i,j))
            wsb(i,j,wstp)=sbc/Dstp(i,j)*hrm(i,j,wnew)**3*cff1**4 &
     &       *(1 -exp (-cff1**5 ) )
# else
          cff1=hrm(i,j,wnew)/(wkb_gam*Dstp(i,j))
          wsb(i,j,wnew) =sbc/Dstp(i,j)*(hrm(i,j,wnew)**3) &
     &      *( 1.0+tanh(8.0*(cff1-1.0)) )*( 1.0-(1.0+cff1**2)**(-2.5) )
# endif
          wcr(i,j,wnew) = frq(i,j,wnew)*(inv_k**2)   ! c/k
          wsr(i,j,wnew) = c_roller*war(i,j,wnew)    & ! ep_r/rho/sigma
     &                    *wvn(i,j,wnew)/max(frq(i,j,wnew)**2,zeps)
     
 !         write(*,*) c_roller, maxval(wvn), maxval(frq),maxval(war)
          inv_f= 1.0/max(frq(i,j,wnew),eps)
          urms = 0.5*frq(i,j,wnew)*hrm(i,j,wnew)/max(sinh(khn),eps)
          abot = inv_f*urms
          fw   = 1.39*(Zob/max(abot,eps))**0.52 ! Soulsby (1995)
          wfc(i,j,wnew) = cfrc*fw*inv_f*urms**3 ! ep_d/rho/sigma
        enddo
      enddo
      
   END SUBROUTINE wkb_stp 
   
      SUBROUTINE wkb_diag
      integer  i,nsubs,ie
      real(wp), dimension(2) :: my_ww
      integer max_check_line
      parameter (max_check_line=128)
      character check_line*(max_check_line), tstring*18
!
! Diagnose difference of wave action density and modulus wavenumber
! between two neighboring time steps to ensure the steady state.
! First, compute thread-integrated deviations, my_wac and my_wkn.
!
      if (mod(iwave-1,winfo) .eq. 0) then
        my_ww(1) = maxval(abs(wac(:,:,wnew)-wac(:,:,wstp)))
        my_ww(2) = maxval(abs(wvn(:,:,wnew)-wvn(:,:,wstp)))
          nsubs=1
!
! Compute global summation for OMP and MPI situations
!
!        if (tile_count.eq.0) then       ! Initialize global sums
          av_wac=my_ww(1)               ! for multithreaded shared
          av_wkn=my_ww(2)               ! memory summation.
!        else                            ! Perform global summation
!          av_wac=av_wac+my_ww(1)        ! among the threads within
!          av_wkn=av_wkn+my_ww(2)        ! each MPI process.
!        endif

            if (first_time.eq.0) then
              first_time=1
              write(*,'(5x,A,1x,A,3x,A,3x,A,2x,A)') &
     &         'STEP','time[DAYS]','DIFF_ACTION','DIFF_WAVENUMBER'
            endif
            write(check_line,'(2(1PE16.10,X),I3)') av_wac, av_wkn
            ie=max_check_line
            do while (check_line(ie:ie).eq.' ' .and. ie.gt.0)
              ie=ie-1
            enddo                                 ! Suppress FORTRAN
            i=0                                   ! floating point Es
            do while (i.lt.ie)                    ! to shorten the
              i=i+1                               ! diagnostic line.
              if (check_line(i:i).eq.'E' .or.&
     &            check_line(i:i).eq.'e') then
                check_line(i:ie-1)=check_line(i+1:ie)
                check_line(ie:ie)=' '
                ie=ie-1
              elseif (ichar(check_line(i:i)).lt.48 .or. &
     &                ichar(check_line(i:i)).gt.57) then
                if (check_line(i:i).ne.' ' .and. &
     &              check_line(i:i).ne.'+'  .and. &! Set may_day_flag
     &              check_line(i:i).ne.'-'  .and.& ! to terminate the
     &              check_line(i:i).ne.'.') then  ! run in the case
     !             may_day_flag=1                  ! of floating point
                endif                             ! exception
              endif
            enddo

            write(tstring,'(F18.8)') real(iwave-1)*rdt*sec2day
            i=1
            do while (i.lt.18 .and. tstring(i:i).eq.' ')
              i=i+1
            enddo
            write(*,'(I9,1x,A,1x,A)') iwave-1, tstring(i:i+9), &
     &                                       check_line(1:ie)


          if (winfo.eq. 1   .and.  iwave.gt.   8) winfo=2
          if (winfo.eq. 2   .and.  iwave.gt.  16) winfo=4
          if (winfo.eq. 4   .and.  iwave.gt.  64) winfo=8
          if (winfo.eq. 8   .and.  iwave.gt. 256) winfo=16
          if (winfo.eq.16   .and.  iwave.gt. 512) winfo=32

      endif              ! <-- mod(iwave-1,winfo).eq.0
   
   END SUBROUTINE wkb_diag
  

END MODULE wkbstp
