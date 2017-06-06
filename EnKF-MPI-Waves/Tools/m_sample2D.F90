module m_sample2D
! This routine samples pseudo random fields with improved independence
! or orthogonality.  This is done by first drawing a large sample
! and construct the final sample using the dominant singular vectors
! of the large sample.

contains
subroutine sample2D(A2,nx,ny,nrens,nre,rh,samp_fix)
#if defined (QMPI)
   use qmpi, only : master, stop_mpi
#else
   use qmpi_fake
#endif
   use m_pseudo2D
   use m_fixsample
   implicit none
   integer, intent(in)     ::  nx,ny
   integer, intent(in)     ::  nrens
   integer, intent(in)     ::  nre
   real,    intent(in)     ::  rh
   logical, intent(in)     ::  samp_fix
   real,    intent(out)    ::  A2(nx,ny,nrens)

   integer n1,n2,n
   integer ns,msx,i,j,nsx,reclA
   integer lwork,ierr,isize,iens
   integer ishape3(3)
   character(len=3) tag3
   real summ,summ2, pi2
   real, allocatable, dimension(:,:,:) :: A,A3,A0,UU
   real, allocatable, dimension(:,:)   :: U,VT,VT1,mean,var, A1
   real, allocatable, dimension(:)     :: sig,work

   logical :: debug=.false.
   
   !ehouarn: x-periodicity in NorESM! 
   logical :: xper=.true.
    
#ifdef SGI
   n1=nx+2*nint(rh); if (mod(n1,2) == 1 ) n1=n1+1
   n2=ny+2*nint(rh); if (mod(n2,2) == 1 ) n2=n2+1
#else
   !ehouarn: x-periodicity
   if(xper)then  
      n1 = nx
   else
      n1 = ceiling(log(real(nx)+2.*rh)/log(2.))
      n1 = 2**n1
   endif
   n2 = ceiling(log(real(ny)+2.*rh)/log(2.))
   n2 = 2**n2
   if (master) print*, 'FFT size',n1,n2
#endif

   n=nx*ny
   ns=nre*nrens
   msx=min(ns,n)
   nsx=min(nrens,n)
   pi2 = 2. * 3.14159253589

!   print *,'sample2D: samp_fix= ',samp_fix

   ! Standard Monte Carlo sampling
   if (.true.) then
   !if (nre == 1) then
      if(master) print *,'NB using standard sampling in sample2D!!'

      if (master) print *,'calling pseudo2d'
      !print *,nx,ny,nrens,n1,n2
      call pseudo2D(A2,nx,ny,nrens,rh,n1,n2)
      if (master) print *,'pseudo2d done'




   ! Start with oversized ensemble of ns members
   elseif (nre > 1) then
      lwork=2*max(3*ns+max(n,ns),5*ns)
      allocate(work(lwork))

      if (master) print*, 'Allocate oversized ensemble'
      allocate(A(nx,ny,ns))
      if (master) print*, 'oversized ensemble allocated'
      call pseudo2D(A ,nx,ny,ns   ,rh,n1,n2)
      if (master) print *,'pseudo2D done'


      ! make an orthogonal VT1 used as linear combination for final ensemble
      if (master) print*, ' make an orthogonal VT1 used as linear combination for final ensemble'
      allocate (A0(nsx,nsx,1), A1(nsx,nsx), U(nsx,nsx), sig(nsx), VT1(nsx,nsx) )

#ifdef SGI
      n1=nsx+2*nint(rh); if (mod(n1,2) == 1 ) n1=n1+1
      n2=nsx+2*nint(rh); if (mod(n2,2) == 1 ) n2=n2+1
#else
      n1 = ceiling(log(real(nsx)+2.*rh)/log(2.))
      n1 = 2**n1
      n2 = n1
      if (master) print*, 'FFT size',n1,n2
#endif

      !print*, ' Re-pseudo2D to generate random orth. matrix'
      !call pseudo2D(A0(:,:,1),nsx,nsx,1,rh,n1,n2)
      call random_number(A0(:,:,1))
      call random_number(A1(:,:))
      A0(:,:,1) = sqrt(-2.*log(A0(:,:,1))) * cos(pi2*A1)
!$OMP CRITICAL
      if (master) print*, ' DGESVD (singular value) '
      call dgesvd('N', 'S', nsx, nsx, A0, nsx, sig, U, nsx, VT1, nsx, work, lwork, ierr)
!$OMP END CRITICAL
      if (ierr /= 0) print *, 'ierr',ierr
      deallocate(A0, sig, U,A1)



      ! Compute SVD of oversized ensemble
      allocate( U(n,msx), sig(msx), VT(msx,msx) )
!$OMP CRITICAL
      if (master) print*, ' Compute SVD of oversized ensemble'
      call dgesvd('S', 'N', n, ns, A, n, sig, U, n, VT, msx, work, lwork, ierr)
!$OMP END CRITICAL
      if (ierr /= 0) print *, 'ierr',ierr
      !print *,'max/min A:',maxval(A),minval(A)
      !print *,'max/min U:',maxval(U),minval(U)

      write(tag3,'(i3.3)')ns

!      inquire(unit=11,opened=lopen)
!      open(11,file='sigma_'//tag3//'.dat')
!         summ=0.0
!         do i=1,ns
!            summ=summ+sig(i)**2
!            write(10,'(i4,3e12.4)')i,sig(i)/sig(1),sig(i)**2/sig(1)**2,summ/real(n*ns)
!         enddo
!      close(11)


      ! Generate first min(nrens,n) members
      if (master) print*, ' Generate first min(nrens,n) members'
      ! KAL -- U is size msx !!
      ishape3=(/nx,ny,nsx/)
      allocate(UU(nx,ny,nsx))
      !KALishape3=(/nx,ny,msx/)
      !KALallocate(UU(nx,ny,msx))
      !print *,size(U)
      !print *,size(U(:,1:nsx))
      !print *,size(UU)
      !print *,minval(U),maxval(U)

      ! KAL 
      UU=reshape(U,ishape3)
      !KALUU=reshape(U(:,1:nsx),ishape3)


      !print *,'aft reshape'
      A2=0.0
      do j=1,nsx
         !print *,'m_sample2D test ',j
         do i=1,nsx
         !KALdo i=1,msx
            A2(:,:,j)=A2(:,:,j)+UU(:,:,i)*sig(i)/sqrt(real(nre))*VT1(i,j)
         enddo
      enddo
      deallocate(U, UU, VT, sig, VT1,A)
      !call stop_mpi()


      if (debug) then
         ! SVD of new ensemble
         allocate (U(n,nsx))
         allocate (sig(nsx))
         allocate (VT(nsx,nsx))
         sig=0.0
!$OMP CRITICAL
         call dgesvd('S', 'S', n, nsx, A2, n, sig, U, n, VT, nsx, work, lwork, ierr)
!$OMP END CRITICAL
         if (ierr /= 0) print *, 'ierr',ierr

!$OMP CRITICAL
         open(10,file='sigma2_'//tag3//'.dat')
            summ=0.0
            do i=1,nsx
               summ=summ+sig(i)**2
               write(10,'(i4,3e12.4)')i,sig(i)/sig(1),sig(i)**2/sig(1)**2,summ/real(n*ns)
            enddo
         close(10)
!$OMP END CRITICAL
         deallocate(U, VT, sig)
         call stop_mpi()
      endif
      
      deallocate(work)
   else
      if (master) print *,'invalid value for nre=',nre
      call stop_mpi()
   endif





   ! subtract mean and correct variance
   
   !print *,'max/min A2 before fixsample:',maxval(A2),minval(A2)
   if (samp_fix) call fixsample(A2,nx,ny,nrens)
   !print *,'max/min A2 after  fixsample:',maxval(A2),minval(A2)

   if (debug) then
      allocate(mean(nx,ny))
      allocate(var(nx,ny))
      mean=0.0
      do j=1,nrens
         mean(:,:)=mean(:,:)+A2(:,:,j)
      enddo
      mean=(1.0/real(nrens))*mean


      var=0.0
      do iens=1,nrens
         do j=1,ny
         do i=1,nx
            var(i,j)=var(i,j)+A2(i,j,iens)**2
         enddo
         enddo
      enddo
      var=(1.0/real(nrens-1))*var

!$OMP CRITICAL
      open(10,file='check.dat')
         do j=1,ny
         do i=1,nx
            write(10,'(2i5,2g13.5)')i,j,mean(i,j),var(i,j)
         enddo
         enddo
      close(10)
!$OMP END CRITICAL

      deallocate(mean)
      deallocate(var)

      call stop_mpi()
   endif

end subroutine sample2D
end module m_sample2D
