module m_param_ensemble
contains

  
   subroutine perturb_param(work,nrens,idm,jdm,rh,var,distrib)
#if defined (QMPI)
   use qmpi, only : master, stop_mpi
#else
   use qmpi_fake, only : master, stop_mpi
#endif
   use m_sample2D
   implicit none
   integer:: nrens,idm,jdm
   real,dimension(idm,jdm,nrens)::work
   character(len=4)::distrib
   real:: rh,var
   integer::i,j,k,npoints
   real:: variance
   real,dimension(idm,jdm)::ave
   real::z
   
   call sample2D(work,idm,jdm,nrens,10,rh,.true.) 
   
   ave(:,:)=0.0
   variance=0.
   npoints=0
   
   do k=1,nrens
      do j=1,jdm
      do i=1,idm 
        ave(i,j)=ave(i,j)+work(i,j,k)
      enddo
      enddo
   enddo
   ave=ave/real(nrens)
   
   do k=1,nrens
      do j=1,jdm
      do i=1,idm
          work(i,j,k)=work(i,j,k)-ave(i,j)
          variance = variance + work(i,j,k)**2.
          work(i,j,k)=var*work(i,j,k)
          npoints=npoints+1    
      enddo
      enddo
   enddo
   
   variance=sqrt(variance/real(npoints))
   
   
   if(trim(distrib)=='logn')then
    do j=1,jdm
    do i=1,idm
        !work(i,j,:)=mp+work(i,j,:)/variance
      do k=1,nrens
	work(i,j,k)=exp(work(i,j,k)/variance-(var**2)/2)
      enddo
    enddo
    enddo
   elseif(trim(distrib)=='norm')then
    do j=1,jdm
    do i=1,idm
      do k=1,nrens
        work(i,j,k)=work(i,j,k)/variance
      enddo
    enddo
    enddo   
   else
      print*,'to be done'
   endif
   
   
   end subroutine perturb_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
end module m_param_ensemble
