module m_fixhycom_eco_metno
!Ehouarn: March 2011
!fixanalysis: remapping of tracers after physical analysis!
!use of remapping subroutines embedded in hycom to interpolate!
!biogeochemical tracer on the analysis grid (dp)!
!Remapping is realized after correction of negative anlaysis dp!

contains
      subroutine hybgen_weno_coefs(s,dp,lc,ci,kk,ks,thin)
      implicit none

      integer kk,ks
      logical lc(kk)
      real    s(kk,ks),dp(kk),ci(kk,ks,2),thin
!
!-----------------------------------------------------------------------
!  1) coefficents for remaping from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the 
!             interfacial values 
!
!     REFERENCE?
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       dp    - initial layer thicknesses (>=thin)
!       lc    - use PCM for selected layers
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       ci    - coefficents for hybgen_weno_remap
!                ci.1 is value at interface above
!                ci.2 is value at interface below
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
!-----------------------------------------------------------------------
!
      real, parameter :: dsmll=1.0e-8
!
      integer j,i
      real    q,q01,q02,q001,q002
      real    qdpjm(kk),qdpjmjp(kk),dpjm2jp(kk)
      real    zw(kk+1,3)

      !compute grid metrics
      do j=2,kk-1
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
        qdpjmjp(j) = 1.0/(dp(j-1) +     dp(j) + dp(j+1))
        dpjm2jp(j) =      dp(j-1) + 2.0*dp(j) + dp(j+1)
      enddo !j
      j=kk
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
!
      do i= 1,ks
        do j=2,kk
          zw(j,3) = qdpjm(j)*(s(j,i)-s(j-1,i))
        enddo !j
          j = 1  !PCM first layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
          else
            q001 = dp(j)*zw(j+1,3)
            q002 = dp(j)*zw(j,  3)
            if (q001*q002 < 0.0) then
              q001 = 0.0
              q002 = 0.0
            endif
            q01 = dpjm2jp(j)*zw(j+1,3)
            q02 = dpjm2jp(j)*zw(j,  3)
            if     (abs(q001) > abs(q02)) then
              q001 = q02
            endif
            if     (abs(q002) > abs(q01)) then
              q002 = q01
            endif
            q    = (q001-q002)*qdpjmjp(j)
            q001 = q001-q*dp(j+1)
            q002 = q002+q*dp(j-1)

            ci(j,i,2) = s(j,i)+q001
            ci(j,i,1) = s(j,i)-q002
            zw(  j,1) = (2.0*q001-q002)**2
            zw(  j,2) = (2.0*q002-q001)**2
          endif  !PCM:WEND
        enddo !j
          j = kk  !PCM last layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0

        do j=2,kk
          q002 = max(zw(j-1,2),dsmll)
          q001 = max(zw(j,  1),dsmll)
          zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
        enddo !j
          zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
          zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

        do j=2,kk-1
          if     (.not.(lc(j) .or. dp(j).le.thin)) then  !don't use PCM
            q01  = zw(j+1,3)-s(j,i)
            q02  = s(j,i)-zw(j,3)
            q001 = 2.0*q01
            q002 = 2.0*q02
            if     (q01*q02 < 0.0) then
              q01 = 0.0
              q02 = 0.0
            elseif (abs(q01) > abs(q002)) then
              q01 = q002
            elseif (abs(q02) > abs(q001)) then
              q02 = q001
            endif
            ci(j,i,1) = s(j,i)-q02
            ci(j,i,2) = s(j,i)+q01
          endif  !PCM:WEND
        enddo !j
      enddo !i
      return
      end subroutine hybgen_weno_coefs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hybgen_weno_remap(si,pi,dpi,ci,&
                                  so,po,ki,ko,ks,thin)
      implicit none
!
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,2),&
             so(ko,ks),po(ko+1),thin
!
!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the 
!             interfacial values 
!             the output is the average of the interpolation
!             profile across each output cell.
!
!     REFERENCE?
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       pi    - initial layer interface depths (non-negative)
!                  pi(   1) is the surface
!                  pi(ki+1) is the bathymetry
!                  pi(k+1) >= pi(k)
!       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
!       ci    - coefficents from hybgen_weno_coefs
!                ci.1 is value at interface above
!                ci.2 is value at interface below
!       ki    - number of  input layers
!       ko    - number of output layers
!       ks    - number of fields
!       po    - target interface depths (non-negative)
!                  po(   1) is the surface
!                  po(ko+1) is the bathymetry (== pi(ki+1))
!                  po(k+1) >= po(k)
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       so    - scalar fields in po-layer space
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
      integer i,k,l,lb,lt
      real    dpb,dpt,qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
!
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
!       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
!
! ---       thin or bottomed layer, values taken from layer above
!
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else

!         form layer averages.
!
!         if     (pi(lb).gt.zt) then
!           write(lp,*) 'bad lb = ',lb
!           stop
!         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            dpt = pi(lt+1)-zt
            dpb = zb-pi(lb)
            qt1 = xt*(xt-1.0)
            qt2 = qt1+xt
            qt0 = 1.0-qt1-qt2
            qb1 = (xb-1.0)**2
            qb2 = qb1-1.0+xb
            qb0 = 1.0-qb1-qb2
            do i= 1,ks
              o = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(qt0*(si(lt,i)  -o) + &
                       qt1*(ci(lt,i,1)-o) + &
                       qt2*(ci(lt,i,2)-o)  )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              sz  = sz + dpb*(qb0*(si(lb,i)  -o) + &
                             qb1*(ci(lb,i,1)-o) + &
                             qb2*(ci(lb,i,2)-o)  )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          else !single layer
            qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
            qt2 = qt1 - 1.0 + (xb+xt)
            qt0 = 1.0 - qt1 - qt2
            do i= 1,ks
              sz=qt0*(si(lt,i)  -o) + &
                qt1*(ci(lt,i,1)-o) + &
                qt2*(ci(lt,i,2)-o) 
              so(k,i) = o + sz
            enddo !i
          endif !layers
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_weno_remap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      subroutine trcr_int2char(cfield,itrac) 
      !function returns the name of the tracer
      !char => integer
      implicit none   
      character(len=8)::cfield     
      integer :: itrac
      
      
      select case (itrac)
         case (1)
           cfield='sco212'
	 case (2)
           cfield='alkali'
	 case (3)
           cfield='phosph'
	 case (4)
           cfield='oxygen'
	 case (5)
           cfield='gasnit'
	 case (6)
           cfield='ano3'
	 case (7)
           cfield='silica'
	 case (8)
           cfield='doc'
	 case (9)
           cfield='poc'
	 case (10)
           cfield='hi'
	 case (11)
           cfield='co3'
	 case (12)
           cfield='phyto'
	 case (13)
           cfield='grazer'
	 case (14)
	   cfield='calciu'
	 case (15)
	   cfield='opal'
	 case (16)
	   cfield='n2o'
	 case (17)
	   cfield='dms'
	 case (18)
	   cfield='fdust'
	 case (19)
	   cfield='iron'
	 case (20)
	   cfield='akw3'
	 case (21)
	   cfield='akb3'
	 case (22)
	   cfield='ak13'
	 case (23)
	   cfield='ak23'
	 case (24)
	   cfield='aksp'
	 case (25)
	   cfield='satoxy'
	 case default
           print *,'tracer unknown',itrac
	   stop
      end select
      
      end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    subroutine conservation_mxlayer(nz,ntrc,trcf,trca,dpf,dpa,epsil)
    implicit none
    integer::nz,ntrc
    real,dimension(nz,ntrc)::trcf,trca
    real,dimension(nz)::dpf,dpa,mask
    real::epsil
    integer::k,m
    
    mask(:)=0.
    do k=1,nz
      if(dpf(k).gt.epsil)then
        mask(k)=1.     
      endif
    enddo
    
    trca(:,:)=0.
    do m=1,ntrc
      trca(1,m)=trcf(1,m)*dpf(1)/dpa(1)
      do k=2,nz
        trca(2,m)=trca(2,m)+mask(k)*trcf(k,m)*dpf(k)
      enddo
      trca(2,m)=trca(2,m)/dpa(2)
      do k=3,nz
        trca(k,m)=trca(2,m)
      enddo
    enddo
    
    
    end subroutine
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1   

    subroutine conservation_wc(nz,trcf,trca,dpf,dpa,epsil)
    implicit none
    integer::nz,ntrc
    real,dimension(nz)::trcf,trca
    real,dimension(nz)::dpf,dpa
    real::epsil
    integer::k,m,l,cpt,k0,p
    real::dptot
    integer,dimension(nz)::mask
    
    k=1
    mask(:)=1
    do while (k.le.nz)
      if(dpf(k).gt.0.)then
        if(dpa(k).gt.0.)then
!      if(dpf(k).gt.epsil)then
!        if(dpa(k).gt.epsil)then
	  trca(k)=trcf(k)*dpf(k)/dpa(k)
	  mask(k)=1	  
	endif
	k=k+1
      else
        m=k
        do while (dpf(m)==0.)
!        do while (dpf(m).le.epsil)
	  m=m+1
	  if (m.ge.nz)then
	    m=nz
	    exit
	  endif  
	enddo
	!m first non null layer 
	cpt=0
	dptot=0.
	do l=k,m
	  if(dpa(l).gt.0.)then
!	  if(dpa(l).gt.epsil)then
	    cpt=cpt+1
	    dptot=dptot+dpa(l)
	  endif
	enddo
	
	if(m.ge.nz)then
	  !all forecasted layers are null below k-1
 	  do l=k-1,m
             if(dpa(l).gt.0.)then
!             if(dpa(l).gt.epsil)then
!	        trca(l)=trcf(k-1)*dpf(k-1)/(dpa(l)*real(cpt+1))
                trca(l)=trcf(k-1)*dpf(k-1)/(dptot+dpa(k-1))
		mask(k-1)=1
	     else
	        trca(l)=trca(l-1)
	     endif 
	  enddo		
	else
   	  do l=k,m
             if(dpa(l).gt.0.)then
!             if(dpa(l).gt.epsil)then
!	        trca(l)=trcf(m)*dpf(m)/(dpa(l)*real(cpt))
	        trca(l)=trcf(m)*dpf(m)/dptot
		mask(m)=1
	     endif 
	  enddo
	endif  

	k=m+1
      endif
    enddo
  
 !!!!!!!!!!!!!! upward to fill out null analysis layers  
   
    k=nz
    do while (dpa(k)==0.)
      k=k-1
    enddo 
    
    do while (k.ge.3)
      if(dpa(k)==0.)then
        m=k
        do while (dpa(m)==0.)
!        do while (dpa(m).le.epsil)
	  m=m-1
	  if (m==2)then
	    exit
	  endif  	  
	enddo
	if(m==3) then
	  if(dpf(k+1).gt.0)then
	    cpt=0
	    dptot=0.
	    do l=k+1,m,-1!m+1,-1
	      if(dpa(l).gt.0.)then
!	     if(dpa(l).gt.epsil)then
	        cpt=cpt+1
		dptot=dptot+dpa(l)
	     endif
	    enddo
	    !trca(m+1:k+1)=0.
	    trca(m:k+1)=0.
	    do l=k+1,m,-1!m+1,-1
	      if(dpa(l).gt.0.)then
	        do p=k+1,m,-1
!	    if(dpf(l).gt.epsil)then
!	          trca(l)=trca(l)+trcf(p)*dpf(p)/(dpa(l)*real(cpt))
	           if(mask(p)==1)then	    
	             trca(l)=trca(l)+trcf(p)*dpf(p)/dptot
	           endif	  
	        enddo
	      endif
	    enddo
	  endif 
	  do l=k,m,-1!m+1,-1
	    if(dpa(l)==0.)then
!	    if(dpa(l).gt.epsil)then
	      trca(l)=trca(l+1)    
	    endif
	  enddo	
	elseif((m.gt.3).or.(m==2)) then
	  if(dpf(k+1).gt.0)then
	    cpt=0
	    dptot=0.
	    do l=k+1,m+1,-1!m+1,-1
	      if(dpa(l).gt.0.)then
!	     if(dpa(l).gt.epsil)then
	        cpt=cpt+1
		dptot=dptot+dpa(l)
	     endif
	    enddo
	    !trca(m+1:k+1)=0.
	    trca(m+1:k+1)=0.
	    do l=k+1,m+1,-1!m+1,-1
	      if(dpa(l).gt.0.)then
	        do p=k+1,m+1,-1
!	    if(dpf(l).gt.epsil)then
!	          trca(l)=trca(l)+trcf(p)*dpf(p)/(dpa(l)*real(cpt))
                  if(mask(p)==1)then
	   	    trca(l)=trca(l)+trcf(p)*dpf(p)/dptot
		  endif   
	        enddo
	      endif 
	    enddo
	  endif 
	  do l=k,m+1,-1!m+1,-1
	    if(dpa(l)==0.)then
!	    if(dpa(l).gt.epsil)then
	      trca(l)=trca(l+1)    
	    endif
	  enddo		
!	else
!	  do l=k,m+1,-1
!	    if(dpa(l)==0.)then
!	  if(dpa(l).gt.epsil)then
!	      trca(l)=trca(l+1)    
!	    endif
!	  enddo
	endif   
!	do l=k,m+1,-1
!	  if(dpa(l)==0.)then
!	  if(dpa(l).gt.epsil)then
!	    trca(l)=trca(l+1)    
!	  endif
!	enddo	
	k=m-1
      else
	k=k-1
      endif
    enddo  
   
    end subroutine
    
!!!!!!!!!!!!!!!!!!    
    subroutine conservation_wc2(nz,trcf,trca,dpf,dpa,epsil)
    implicit none
    integer::nz,ntrc
    real,dimension(nz)::trcf,trca
    real,dimension(nz)::dpf,dpa
    real::epsil
    integer::k,m,l,cpt,k0
    
    k=1
    do while (k.le.nz)
      
      if(dpf(k).gt.0.)then
        if(dpa(k).gt.0.)then
!      if(dpf(k).gt.epsil)then
!        if(dpa(k).gt.epsil)then
	  trca(k)=trcf(k)*dpf(k)/dpa(k)
	else
	  trca(k)=trca(k-1)	  
	endif
	k=k+1
      else
        m=k
        do while (dpf(m)==0.)
!        do while (dpf(m).le.epsil)
	  m=m+1
	  if (m.ge.nz)then
	    m=nz
	    exit
	  endif  
	enddo
	!m first non null layer 
	cpt=0
	do l=k,m
	  if(dpa(l).gt.0.)then
!	  if(dpa(l).gt.epsil)then
	    cpt=cpt+1
	  endif
	enddo
	
	if(m.ge.nz)then
	  !all forecasted layers are null below k-1
 	  do l=k-1,m
             if(dpa(l).gt.0.)then
!             if(dpa(l).gt.epsil)then
	        trca(l)=trcf(k-1)*dpf(k-1)/(dpa(l)*real(cpt+1))
	     else
	        trca(l)=trca(l-1)
	     endif 
	  enddo
	
	
	else
   	  do l=k,m
             if(dpa(l).gt.0.)then
!             if(dpa(l).gt.epsil)then
	        trca(l)=trcf(m)*dpf(m)/(dpa(l)*real(cpt))
	     else
	       trca(l)=trca(l-1)
	     endif 
	  enddo

	endif  

	k=m+1
      endif
    
    enddo
   
    end subroutine
    
    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine conservation_wc3(nz,trcf,trca,dpf,dpa,epsil,verbose)
    implicit none
    integer::nz,ntrc
    real,dimension(nz)::trcf,trca
    real,dimension(nz)::dpf,dpa
    real::epsil
    integer::k,m,l,cpt,k0,p
    real::dptot,trctot
    logical::verbose
    
    k=1
    do while (k.le.nz)   
!      if((dpf(k).gt.0.).and.(dpa(k).gt.0.))then
      if((dpf(k).gt.epsil).and.(dpa(k).gt.epsil))then      
	trca(k)=trcf(k)*dpf(k)/dpa(k)
        k=k+1      
      else
	m=k
	dptot=0.
	trctot=0.
!	do while((dpf(m)==0.).or.(dpa(m)==0.))  	  
	do while((dpf(m).le.epsil).or.(dpa(m).le.epsil))  
          if(dpa(m).gt.epsil) dptot=dptot+dpa(m)
	  if(dpf(m).gt.epsil) trctot=trctot+trcf(m)*dpf(m)
	  m=m+1
!	  if(verbose)then
!	    print*,'m ',m-1+nz,dptot,trctot
!	  endif
	  if(m.ge.nz)then
	    m=nz
	    exit
	  endif  
        enddo
	if(m.lt.nz)then
	  dptot=dptot+dpa(m)
	  trctot=trctot+trcf(m)*dpf(m)
	endif  
!	if(verbose)then
!	    print*,'m ',m+nz,dptot,trctot
!	endif
	
	if(dptot==0.)then
	  trca(k-1)=(trctot+trca(k-1)*dpa(k-1))/dpa(k-1)
	  do l=k,m
	    trca(l)=trca(l-1)
	  enddo
	elseif(trctot==0.)then  
	  do l=k-1,m
!	    if(dpa(l).gt.0)then
	    if(dpa(l).gt.epsil)then
	      trca(l)=trcf(k-1)*dpf(k-1)/(dptot+dpa(k-1))
	      k0=l
	    endif
	  enddo 
	  do l=k0+1,nz
	    trca(k)=trca(k-1)
	  enddo
	  do l=k0,k,-1
!	    if(dpa(l)==0.)then
	    if(dpa(l).le.epsil)then
	      trca(l)=trca(l+1)
	    endif
	  enddo 
	else	
	  do l=k,m
!	    if(dpa(l).gt.0.)then
	    if(dpa(l).gt.epsil)then
	      trca(l)=trctot/dptot
	      k0=l
	    endif
	  enddo
	  if(m==nz)then
	    do l=k0+1,nz
	      trca(k)=trca(k-1)
	    enddo
	    do l=k0,k,-1
!	      if(dpa(l)==0.)then
	      if(dpa(l).le.epsil)then
	        trca(l)=trca(l+1)
	      endif
	    enddo  
	  else
	    do l=m,k,-1
!	      if(dpa(l)==0.)then
	      if(dpa(l).le.epsil)then
	        trca(l)=trca(l+1)
	      endif
	    enddo  
	  endif
	  
	endif
	k=m+1
	
      endif
    enddo
  
 
   
    end subroutine
    
end module
