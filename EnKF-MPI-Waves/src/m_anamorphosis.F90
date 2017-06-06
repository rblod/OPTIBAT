!-----------------------------------------------------------------------
!BOP
!
! !MODULE: m_anamorphosis
!
! !INTERFACE:
   MODULE m_anamorphosis
!
! !DESCRIPTION: 
!   !module including the experimental anamorphosis routines for the EnKF! 

! !USE:
!   List of mudules used (can be empty)
   IMPLICIT NONE
!


! !PUBLIC MEMBER FUNCTIONS:
!
! !PUBLIC DATA MEMBERS:
!
    type anamorphosis
        character(len=3) smooth !nature of smoothness of the empirical anamorphosis!
        character(len=3) id !nature of the physical variable!
        integer::n !size of the sample!
        real, dimension(:), allocatable::y  !gaussian values invG(i/N)!
        real, dimension(:), allocatable::z  !sorted physical values!
        real, dimension(:), allocatable::phi !anamorphosis coefficients!
        real::Gmin,Gmax !bounds in gaussian space!
        real::Zmin,Zmax !bounds in real space!
        real:: alpha !linear interpolation: translation coefficient!
        integer:: m !0:linear tail, 1:nonlinear tail! 
        integer::nsamp,samp !size of the sampling after clustering!
        real::eps !threshold (lowest relevant value)!
	integer::ifirst
	real::a1,a2 !param of the left tail: nonlinear function!
	logical::obs
    end type anamorphosis

    real,parameter,private::onem=9806 !NORWECOM
    character(len=*), parameter, public :: infile_ana='analysisfields_ana.in'
    integer,save :: numfields_ana
    character(len=3), dimension(:), save, allocatable:: fieldnames_ana
    type(anamorphosis),dimension(:),save,allocatable::ana_enkf
    integer,dimension(:),save, allocatable::obs_bias
      
    integer::nana_stat,nana_dyn,nana_par
    integer::nana_statP
    parameter(nana_stat=241,nana_statP=250,nana_dyn=223,nana_par=49)
    real::gp1,gm1
    parameter(gm1=0.1587,gp1=0.8413)
   
    character(len=8),private::biomodel_ana='NORWECOM'
    integer, parameter,private:: nzchl=40!15 
    
    integer,parameter,private::nens_stat=9
    
    
! !REVISION HISTORY:
!  Author(s): NAME Ehouarn
!  01Jan2000: Ver. ?.?.? (???): 
!
! !LOCAL VARIABLES:
!
! !BUGS
!
!EOP
!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
  
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!         intrinsec functions     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      integer function ana_get_ind(ana,nana,char3) 
      implicit none
      integer :: nana
      type(anamorphosis),dimension(nana)::ana
      character(len=3) :: char3
      integer::k

      do k=1,nana
        if(trim(ana(k)%id)==trim(char3))then
          ana_get_ind=k
          return
       endif  
      enddo
      ana_get_ind=-1
   
      end function
   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!         reading of analysisfields_ana.in     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function get_nrfields_ana()
      implicit none
      integer :: ios,first,last
      logical :: ex
      character(len=3) :: char3

      inquire(exist=ex,file=infile_ana)
      if (.not. ex) then
        print *,'Could not find '//infile_ana
      end if

      open(103,status='old',form='formatted',file=infile_ana)
      ios=0
      get_nrfields_ana=0
      do while (ios==0)
        read(103,100,iostat=ios) char3
        if (ios==0) get_nrfields_ana=get_nrfields_ana+1
      end do
      close(103)
  100 format (a3,2i3)
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   
      subroutine get_analysisfields_ana()
      implicit none
      integer :: first,last,k,nfld,ios
      logical :: ex
      character(len=3) :: char3

      numfields_ana=get_nrfields_ana()
      print *,'numfields_ana is ',numfields_ana
      if (numfields_ana<=0 .or.numfields_ana > 20) then !
         print *,'numfields_ana is higher than max allowed setting or = 0'
      end if
      allocate(fieldnames_ana(numfields_ana))

      inquire(exist=ex,file=infile_ana)
      if (.not. ex) then
        print *,'Could not find '//infile_ana
      end if

      open(103,status='old',form='formatted',file=infile_ana)
      ios=0
      nfld=0
      do while (ios==0)
        read(103,100,iostat=ios) char3
        if (ios==0) then
          fieldnames_ana (nfld+1)=char3
          nfld=nfld+1
        end if
      end do
      close(103)
  100 format (a3,2i3)

      if (nfld/=numfields_ana) then
        print *,'An error occured when reading '//infile_ana
      end if

      ! List fields used in analysis
      do k=1,numfields_ana
        print *,fieldnames_ana(k)
      end do

      end subroutine      
      
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!     physical data    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   
  subroutine ana_data_z_stat(ana,work,nrens,nz,obs,obstot,nwork,m0)  
  use mod_observations
  use m_bio_estimation
  use m_sort2
  use netcdf
  implicit none
  type(anamorphosis)::ana
  integer::nrens,m0,obstot,nz,nspec,nwork
  type(observations)::obs(obstot)
  real,dimension(nwork)::work 
  integer, dimension(nwork)::itemp
  integer::cptz,i,j,k,it,nzz

  character(len=8) :: cfield
  logical::ex
  integer::ncid,varid,ierr,nindex,ntime
  real,dimension(:,:),allocatable::fld
  real,dimension(:,:,:),allocatable::fld2
  integer::tincr,syear,mp
  integer,dimension(nens_stat)::dstat
    
  cptz=1
  
  !getting obs data!
  ierr=nf90_open('observation.nc',nf90_nowrite,ncid)
  ierr=nf90_inq_dimid(ncid,'time',varid)
  if(ierr/=nf90_noerr)stop
  ierr=nf90_inquire_dimension(ncid,varid,len=ntime) 
  if(ierr/=nf90_noerr)stop
  
  syear=365
  tincr=30!45
  print*,'stat tincr',tincr
  if(tincr.ge.m0)then
    mp=1
  else
    mp=m0-tincr+1!eho 190811
  endif
  
  allocate(fld(1:nz,1:ntime))
  do k=1,obstot         
    ierr=nf90_inq_varid(ncid,trim(ana%id),varid)
    if(ierr/=nf90_noerr)stop
    ierr=nf90_get_var(ncid,varid,fld)
    do j=mp,m0
      work(cptz)=fld(obs(k)%k,j)
      cptz=cptz+1
    enddo
  enddo
  deallocate(fld)
  ierr=nf90_close(ncid)
  
  !second set of osbervation!
  ierr=nf90_open('observation2.nc',nf90_nowrite,ncid)

  allocate(fld(1:nz,1:ntime))
  do k=1,obstot         
    ierr=nf90_inq_varid(ncid,trim(ana%id),varid)
    if(ierr/=nf90_noerr)stop
    ierr=nf90_get_var(ncid,varid,fld)
    do j=mp,min(m0+tincr-1,ntime)
      work(cptz)=fld(obs(k)%k,j)
      cptz=cptz+1
    enddo
  enddo
  deallocate(fld)
  ierr=nf90_close(ncid)
  
  !getting model data!
 
  cfield=trim(ana%id)
  ierr=NF90_OPEN('Ana_data_set.nc' ,NF90_NOWRITE,ncid)
  ierr=nf90_inq_dimid(ncid,'time',varid)
  if(ierr/=nf90_noerr)stop
  ierr=nf90_inquire_dimension(ncid,varid,len=ntime) 
  if(ierr/=nf90_noerr)stop

  allocate(fld2(1:nrens,1:nz,1:ntime))
  ierr=nf90_inq_varid(ncid,trim(cfield),varid)
  if(ierr/=nf90_noerr)stop
  ierr=nf90_get_var(ncid,varid,fld2)
  if(ierr/=nf90_noerr)stop
  ierr=nf90_close(ncid)
  
  
  if(m0.gt.4*syear)then
    mp=m0-2*syear
  elseif(m0.gt.3*syear)then
    mp=m0-syear
  elseif(m0.gt.2*syear)then
    mp=m0!-syear
  elseif(m0.gt.syear)then
    mp=m0+syear
  else
    mp=m0+2*syear
  endif
  print*,'date m',mp,m0
  
  call compute_dstat(dstat,nens_stat)
  if(trim(ana%id)=='chl')then
    nzz=nz-nzchl+1
  else
    nzz=1
  endif
  
  do k=1,nens_stat
    do j=nz,nzz,-1
      do i=mp-tincr,mp+tincr
       work(cptz)=fld2(dstat(k),j,i)
       cptz=cptz+1
       if(fld2(dstat(k),j,i).lt.0.) print*,'negative param',m0,trim(cfield),k
      enddo   
    enddo
  enddo
  deallocate(fld2)
  
  !print*,'cptz',cptz,nwork
  call quick_sort(work,itemp)
   
  if(ana%obs)then
  open(unit=10,FILE='work_chl.txt',FORM='formatted',STATUS='unknown')
  do k=1,nwork
    write(10,*)k,work(k)
  enddo
  close(10)
  endif
  
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_dstat(dstat,nens)
  implicit none
  integer::nens
  integer,dimension(nens)::dstat
  integer::k
  
  do k=1,nens
    dstat(k)=10*k+1
  enddo
  
  end subroutine
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
  subroutine ana_data_z_dyn(ana,work,nrens,nz,nspec,ccda,obs,obstot,nwork,m0,dt,nstep)  
  use mod_observations
  use m_bio_estimation
  use m_sort2
  use netcdf
  implicit none
  type(anamorphosis)::ana
  integer::nrens,m0,obstot,nz,nspec,nstep
  type(observations)::obs(obstot)
  real,dimension(nwork)::work 
  real,dimension(1:nspec,0:nz,1:nrens)::ccda
  real::dt
  
  integer::cptz,i,j,k,it

  character(len=8) :: cfield
  logical::ex
  integer::ncid,varid,ierr,nindex,ntime,var_id 
  integer::nwork
  integer, dimension(nwork)::itemp
  real,dimension(:,:),allocatable::fld
  real,dimension(:),allocatable::tmpo
  character(len=3) :: cipx
  integer::tincr,syear,mp,idvar
  
  
  ierr=nf90_open('observation.nc',nf90_nowrite,ncid)
  ierr=nf90_inq_dimid(ncid,'time',var_id)
  if(ierr/=nf90_noerr)stop
  ierr=nf90_inquire_dimension(ncid,var_id,len=ntime) 
  if(ierr/=nf90_noerr)stop
  
  
  syear=365
!  if(m0.gt.syear)then
!    mp=m0-syear
!  else
!    mp=m0
!  endif
  
  !tincr=floor(45.*3600.*24./(dt*real(nstep))) !fixed 3 months window
  tincr=30!45
  print*,'dyn tincr',tincr
  if(tincr.ge.m0)then
    mp=1
  else
    mp=m0-tincr+1
  endif
  
  cptz=1
  idvar=get_ind_varmodel(ana%id,biomodel_ana)
  
  if(idvar.gt.0)then     
    do k=1,nrens
      do j=1,nz
        work(cptz)=ccda(idvar,j,k)
        cptz=cptz+1
      enddo
    enddo
    do k=1,obstot
      if(ana%obs)then
         allocate(fld(1:nz,1:ntime))
	 ierr=nf90_inq_varid(ncid,trim(ana%id),var_id)
         if(ierr/=nf90_noerr)stop
         ierr=nf90_get_var(ncid,var_id,fld)
	 do j=mp,m0
	   work(cptz)=fld(obs(k)%k,j)
	   cptz=cptz+1
	 enddo
	 deallocate(fld)
      endif
    enddo
  elseif(trim(ana%id)=='chl')then
   do k=1,nrens
      do j=nz,nz-nzchl+1,-1
#if defined (CHLA)
        work(cptz)=(ccda(get_ind_varmodel('chf',biomodel_ana),j,k)&
	            +ccda(get_ind_varmodel('chd',biomodel_ana),j,k))
#else
        work(cptz)=(ccda(get_ind_varmodel('fla',biomodel_ana),j,k)&
	            +ccda(get_ind_varmodel('dia',biomodel_ana),j,k))&
		    /n2chla
#endif
        cptz=cptz+1
      enddo
    enddo
    do k=1,obstot      
       allocate(fld(1:nz,1:ntime))
       ierr=nf90_inq_varid(ncid,trim(ana%id),var_id)
       if(ierr/=nf90_noerr)stop
       ierr=nf90_get_var(ncid,var_id,fld)
       do j=mp,m0
          work(cptz)=fld(obs(k)%k,j)
	  cptz=cptz+1
       enddo
       deallocate(fld)
    enddo
  else
    !cfield=trim(ana%id)
    call ana_get_name_parmodel(trim(ana%id),cfield)
    call get_parmodel(work,nrens,cfield,biomodel_ana) 
  endif  
 
  ierr=nf90_close(ncid)
  
  call quick_sort(work,itemp)
  
  
  if(trim(ana%id)=='chl')then
    open(unit=104,FILE='work.txt',FORM='formatted',STATUS='unknown')
    do k=1,nwork
      write(104,*)k,work(k)
    enddo
    close(104)
  endif
  
  end subroutine


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!    smoothness param   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   
   subroutine ana_smooth_param(ana)
   implicit none
   type(anamorphosis)::ana
   character(len=80) ::memfile
   logical::ex
   
   memfile='infile_ana.'//trim(ana%id)//'.in'
   
   inquire(exist=ex,file=memfile)
   if (.not. ex) then
      print *,'Could not find '//memfile
   end if
   
   open(105,status='old',form='formatted',file=memfile)
   read(105,*) ana%smooth
   read(105,*) ana%Gmin
   read(105,*) ana%Gmax
   read(105,*) ana%Zmin   
   read(105,*) ana%Zmax
   read(105,*) ana%alpha
   read(105,*) ana%m
   read(105,*) ana%nsamp
   read(105,*) ana%eps 
   read(105,*) ana%obs
   close(105)

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!    sampling   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
   subroutine ana_samp_z(ana)
   implicit none
   type(anamorphosis)::ana
   !to improve the clustering process: now only a value each samp value!
   
   if(mod(ana%n-2,ana%nsamp)==0)then
      ana%samp=(ana%n-2)/ana%nsamp+1
   else
      ana%samp=(ana%n-2)/ana%nsamp+2
   endif
!   print *,'ana_samp',ana%samp 
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   
   subroutine ana_gauss_y(ana)
   implicit none
   type(anamorphosis)::ana
   integer::k
   real::tmp
   
   allocate(ana%y(ana%samp))
   do k=1,ana%samp-1
      tmp=real(1+(k-1)*ana%nsamp)/real(ana%n)
      call gaussinv(tmp,ana%y(k))
   enddo
   tmp=real(ana%n-1)/real(ana%n)
   call gaussinv(tmp,ana%y(ana%samp))
     
   end subroutine
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!    phi   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      
   
   subroutine ana_cpt_phi_eco3(ana,work,nwork)
   implicit none
   integer::nwork
   type(anamorphosis)::ana
   real,dimension(nwork)::work
   
   logical::alert
   integer::k,kz,kzp1,ksed,k0,kl
   real::y,tmp
   character(len=80) :: file
   
   alert=.false.
   
   if((ana%smooth=='lin').or.(ana%smooth=='li2'))then
      !linear interpolation!
      allocate(ana%phi(ana%samp+1)) 
      
      tmp=ana%alpha
      !left tile!               
      !y=tmp*ana%y(2)+(1.-tmp)*ana%y(1) !ehouarn 21/12
      y=tmp*ana%y(3)+(1.-tmp)*ana%y(2)
      ana%phi(2)=(work(1+ana%nsamp)-work(1))/(y-ana%y(1))     
      ana%z(2)=work(1+ana%nsamp)
      
      !left tile redefining Gmin
      ana%z(1)=work(1)
      ana%phi(1)=ana%phi(2)
      if(ana%phi(1).ne.0.)then
        ana%Gmin=ana%y(1)-(work(1)-ana%Zmin)/ana%phi(1)
      endif
      
      !approx!
      do k=3,ana%samp-1
       y=tmp*ana%y(k+1)+(1.-2.*tmp)*ana%y(k)+(tmp-1.)*ana%y(k-1)
       kz=1+(k-2)*ana%nsamp
       kzp1=1+(k-1)*ana%nsamp
       ana%phi(k)=(work(kzp1)-work(kz))/y
       
       ana%z(k)=work(kzp1)
      enddo
      
      y=ana%y(ana%samp)-(tmp*ana%y(ana%samp)+(1.-tmp)*ana%y(ana%samp-1))
      kz=1+(ana%samp-2)*ana%nsamp
      ana%phi(ana%samp)=(work(ana%n)-work(kz))/y
      ana%z(ana%samp)=work(ana%n)
      
      !right tile!
      y=ana%y(ana%samp)
      ana%phi(ana%samp+1)=ana%phi(ana%samp)
      if(ana%phi(ana%samp+1).ne.0.)then
        ana%Gmax=(ana%Zmax-work(ana%n))/ana%phi(ana%samp+1)+y
      endif
      
      k=1
      do while (k.lt.ana%samp+1)
        if(ana%phi(k).ne.0.)then
	  k=k+1
	else
	  do k0=k,ana%samp+1
	    if(ana%phi(k0).ne.0.)exit  	    
	  enddo
	  !if((k0==ana%samp+2).and.(ana%phi(k0-1)==0.))then !test redondant
	  if(k0==ana%samp+2)then
	    !no larger values..    
	    if(k==1)then
	      print*,'alert: anamorposis null ',trim(ana%id)
	      stop
	    else
	      !k>=3
	      !y=ana%y(ana%samp)-(tmp*ana%y(k-2)+(1.-tmp)*ana%y(k-3))
              !do ksed=k-1,k0-1
	      !  ana%phi(ksed)=(ana%z(ana%samp)-ana%z(k-2))/y
	      !enddo
	      if(ana%z(ana%samp)==ana%Zmax)then
	        !special case: compute Gmax first
		y=real(ana%samp)/(real(ana%samp)+1.)
	        call gaussinv(y,ana%Gmax)
		
		y=ana%Gmax-(tmp*ana%y(k-1)+(1.-tmp)*ana%y(k-2))
	        do ksed=k-1,k0-1
	          ana%phi(ksed)=(ana%Zmax-ana%z(k-2))/y		
	        enddo
	      
	        do ksed=k-1,ana%samp-1
                  y=tmp*(ana%y(ksed+1)-ana%y(k-2))&
		                   +(1.-tmp)*(ana%y(ksed)-ana%y(k-3))
		  ana%z(ksed)=ana%phi(ksed)*y+ana%z(k-2)
	        enddo
		ana%z(ana%samp)=ana%phi(ana%samp)*ana%y(ana%samp)+ana%z(k-2)
		
	      else	      
	        y=ana%y(ana%samp)-(tmp*ana%y(k-1)+(1.-tmp)*ana%y(k-2))
                do ksed=k-1,k0-1
	          ana%phi(ksed)=(ana%z(ana%samp)-ana%z(k-2))/y		
	        enddo
	        ana%Gmax=(ana%Zmax-ana%z(k-2))/ana%phi(ana%samp+1)&
	                +tmp*ana%y(k-1)+(1.-tmp)*ana%y(k-2)
	        do ksed=k-1,ana%samp-1
                  y=tmp*(ana%y(ksed+1)-ana%y(k-2))&
		                   +(1.-tmp)*(ana%y(ksed)-ana%y(k-3))
		  ana%z(ksed)=ana%phi(ksed)*y+ana%z(k-2)
	        enddo
	      endif
	    endif 
	  else
	    if(k==1)then
	      !y=tmp*ana%y(k0)+(1.-tmp)*ana%y(k0-1)-ana%y(1)
	      !do ksed=1,k0
	      !  ana%phi(ksed)=(ana%z(k0)-ana%z(1))/y
	      !enddo
	      
	      if(ana%z(1)==ana%Zmin)then
	        !special case: need to define Gmin first
	        !Zmin=z1=...=z(k0-1)
		y=1./(real(ana%samp)+1.)
	        call gaussinv(y,ana%Gmin)
		
		y=tmp*ana%y(k0+1)+(1.-tmp)*ana%y(k0)-ana%Gmin
		
		do ksed=1,k0
	          ana%phi(ksed)=(ana%z(k0)-ana%Gmin)/y
	        enddo
		do ksed=1,k0-1
	       	  ana%z(ksed)=ana%Zmin+ana%phi(ksed)*&
		         (tmp*ana%y(ksed+1)+(1.-tmp)*ana%y(ksed)-ana%Gmin)
	        enddo
	      else	      	      
	        y=tmp*ana%y(k0+1)+(1.-tmp)*ana%y(k0)-ana%y(1)
	      
	        do ksed=1,k0
	          ana%phi(ksed)=(ana%z(k0)-ana%z(1))/y
	        enddo
	        do ksed=2,k0-1
	       	  ana%z(ksed)=ana%z(1)+ana%phi(ksed)*&
		           (tmp*ana%y(ksed+1)+(1.-tmp)*ana%y(ksed)-ana%y(1))
	        enddo
	        !ana%z(1)=ana%phi(1)*ana%y(1)+ana%z(1)	      
	        ana%Gmin=ana%y(1)-(ana%z(1)-ana%Zmin)/ana%phi(1)
	      endif
	    
	    else
	      !y=tmp*ana%y(k0)+(1.-tmp)*ana%y(k0-1)-tmp*ana%y(k-2)+(1.-tmp)*ana%y(k-3)
	      !do ksed=k-1,k0
	      !  ana%phi(ksed)=(ana%z(k0)-ana%z(k-2))/y
	      !enddo
	      y=tmp*ana%y(k0+1)+(1.-tmp)*ana%y(k0)-tmp*ana%y(k)+(1.-tmp)*ana%y(k-1)
	      do ksed=k,k0
	        ana%phi(ksed)=(ana%z(k0)-ana%z(k-1))/y
	      enddo
	      do ksed=k,k0-1
	        ana%z(ksed)=ana%z(k-1)+ana%phi(ksed)*&
		  (tmp*(ana%y(ksed+1)-ana%y(k))+(1.-tmp)*(ana%y(ksed)-ana%y(k-1)))
	      enddo
	      	      
	    endif
	  endif
	  k=k0
	endif
      enddo
      
   endif
   
      
   file='ana_'//trim(ana%id)//'_phi.txt'
   open(unit=111,FILE=file,FORM='formatted',STATUS='unknown')
   file='ana_'//trim(ana%id)//'_y.txt'
   open(unit=131,FILE=file,FORM='formatted',STATUS='unknown')
   file='ana_'//trim(ana%id)//'_zs.txt'
   open(unit=121,FILE=file,FORM='formatted',STATUS='unknown')
   
   write(131,*)ana%Gmin
   write(121,*)ana%Zmin
   do k=1,ana%samp
      write(111,*)ana%phi(k)
      write(131,*)ana%y(k)
      write(121,*)ana%z(k)
      if(ana%phi(k)==0.)then
        print*,'alert phi=0',trim(ana%id),k
	alert=.true.
      endif
   enddo   
   write(131,*)ana%Gmax
   write(111,*)ana%phi(ana%samp+1)   
   write(121,*)ana%Zmax
   
   close(111)
   close(131)
   close(121)   
   
   if((ana%phi(ana%samp)==0.).and.(ana%phi(ana%samp+1)==0.))then
      print*,'alert phi=0',trim(ana%id),ana%samp+1
      alert=.true.
   endif

   
   end subroutine   
 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!    gaussian distribution   !!!!!
	   !!!!!   issued from A. Hollard  !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
   subroutine pdgauss(y,out)
   !  fonction de repartition gaussienne
   implicit none


   real,intent(in)::y
   real,intent(out)::out

   real,dimension(5)::b
   real             ::p,t,tn,u,x
   integer::n

   p=0.2316419
   b(1)=0.319381530
   b(2)=-0.356563782
   b(3)=1.781477937
   b(4)=-1.821255978
   b(5)=1.330274429


   if (y < 0.) then 
     x=-y
   else
     x=y
   endif
! then we proceed with the calculus
   t=1./(1.+p*x)
   u=0.
   tn=1.
   do n=1,5
     tn=tn*t
     u=u+tn*b(n)
   enddo
   u=1.-0.3989423*u*dexp(-0.5*x*x)
! then we give the general value of pdgauss
   out=u
   if (y < 0.) then
     out=1.-u
   endif
 
   end subroutine
   
   subroutine gaussinv(u,gau)
   !  fonction de repartition de gauss inverse
   !      u= probability
   !  we called in this section the previous function pdgauss
   !       
   implicit none

   real,intent(in)::u
   real,intent(out)::gau
   integer::n
   real::v,x,y,x0,y0,x1,y1,eps

!  check out the validity of the arguments

   if ((u <= 0.).or.(u>=1.)) then
     print *,'there is a problem with your probability'
     stop
   endif

!  restriction to the case >0.5
   v=u
   if (u < 0.5) then
     v=1.-u
   endif

!  precision of the calculus 
   eps=0.000001

!  initialisation of the dichotomy  
   x0=0.
   y0=0.5
   x1=0.
   y1=0.
        
   do while (y1 <= v)
     x1=x1+1.
     call pdgauss(x1,y1)
   enddo
   	
!  dichotomie
   do while ((dabs(x1-x0) > eps).or.(dabs(y1-y0) > eps))
     x=(x0+x1)/2.
     call pdgauss(x,y)
     if (y > v)then
       x1=x
       y1=y
     else
       x0=x
       y0=y
     endif
   enddo

!  end of the dichotomy
   x=(x0+x1)/2.

!  return to the genral case
   gau=x
   if (u < 0.5) then
     gau=-x
   end if

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!   anamorphosis functions in EnKF   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ana_exp_eco22(y,z,ana)
  use m_sort2
  implicit none
  real :: y,z
  type(anamorphosis)::ana
  integer::ind,i
  real::b,tmp
  real,dimension(:),allocatable::yt
  real::d

#if defined (ANA_LOG)
    if(y.le.-40.)then
      z=ana%Zmin
      print*,'ltail ana_exp_eco2',trim(ana%id),z 
      return
    elseif(y.ge.log(ana%Zmax))then 
      z=ana%Zmax
      print*,'rtail ana_exp_eco2',trim(ana%id),z
      return
    endif  
    
    z=exp(y)
#else

!    if(y.le.ana%Gmin)then     
!      if(ana%m==0)then
!        !linear tail
!	z=ana%Zmin
!      else
!        call ana_ltail(ana,y,z)
!      endif
!      return
!    elseif(y.ge.ana%Gmax)then       
!      if(ana%m==0)then
!        !linear tail
!	z=ana%Zmax 
!      else
!        call ana_rtail(ana,y,z)
!      endif      
!      return
!    endif  
    
    if(trim(ana%id)=='chl')then
      if(y.le.-40.)then
        z=ana%Zmin
        print*,'ltail ana_exp_eco2',trim(ana%id),z 
        return
      elseif(y.ge.log(ana%Zmax))then 
        z=ana%Zmax
        print*,'rtail ana_exp_eco2',trim(ana%id),z
        return
       endif  
       print*,'connard log'
      z=exp(y)
    else
    
    if(ana%m==0)then
      if(y.le.ana%Gmin)then           
        !linear tail
        z=ana%Zmin
        return 
      elseif(y.ge.ana%Gmax)then            
        z=ana%Zmax 
        return            
      endif 
      
    else
      if(y.le.ana%y(1))then           
        !linear tail
        call ana_ltail(ana,y,z)
        return
      elseif(y.ge.ana%y(ana%samp))then            
        call ana_rtail(ana,y,z)         
        return
      endif 
    endif 
    
    
    tmp=ana%alpha      
    allocate(yt(1:ana%samp+2))
    yt(1)=ana%Gmin
    yt(2)=ana%y(1)
    do i=2,ana%samp-1
      yt(i+1)=tmp*ana%y(i+1)+(1.-tmp)*ana%y(i)
    enddo
    yt(ana%samp+1)=ana%y(ana%samp)
    yt(ana%samp+2)=ana%Gmax
    
    call search_value(yt,y,ana%samp+2,ind)
    
    if(ind==ana%samp+1)then
      b=ana%Zmax-ana%phi(ana%samp+1)*ana%Gmax
    elseif(ind==ana%samp)then  
      b=ana%z(ind)-ana%phi(ind)*ana%y(ana%samp)    
    elseif(ind==1)then
      b=ana%z(1)-ana%phi(1)*ana%y(1)
    else
      !b=ana%z(ind)-ana%phi(ind)*yt(ind) !21/12
      b=ana%z(ind)-ana%phi(ind)*yt(ind+1)
    endif
    
    z=ana%phi(ind)*y+b

    if(ana%m==1)then
    !nonlinear tail
      if(ind==ana%samp+1)then
        call ana_rtail(ana,y,z)
        !print*,'rtail ana_exp_eco2',trim(ana%id),y     
      elseif(ind==1)then
        call ana_ltail(ana,y,z)
        !print*,'ltail ana_exp_eco2',trim(ana%id),y	   
      endif 
    endif
    
    deallocate(yt)
    
    endif !ana%chl=log
    
!LOG!    
#endif 
   
  end subroutine  
   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine ana_exp_inv_eco22(z,y,ana)
  use m_sort2
  implicit none
  real :: y,z
  type(anamorphosis)::ana
  integer::ind,i
  real::b,tmp
  real,dimension(:),allocatable::zt
  real::d

#if defined (ANA_LOG)
    if(z.le.ana%Zmin)then
      y=-40.     
      return
    elseif(z.ge.ana%Zmax)then !infinite rtail
      y=log(ana%Zmax)     
      return
    endif

    y=log(z)
 
#else
    if(trim(ana%id)=='chl')then
      if(z.le.ana%Zmin)then
        y=-40.     
        return
      elseif(z.ge.ana%Zmax)then !infinite rtail
        y=log(ana%Zmax)     
        return
      endif
      print*,'connard2 log'
      y=log(z)  
    else
    
     
    if(ana%m==0)then
      if(z.le.ana%Zmin)then
        y=ana%Gmin
        return
      elseif(z.ge.ana%Zmax)then 
        y=ana%Gmax  
        return
      endif  
    else
      if((z.le.ana%Zmin).or.(z.ge.ana%Zmax))then
        print*,'alert z out of bounds ',trim(ana%id),z
        print*,'ana_bounds: ',ana%Zmin,ana%Zmax
        stop
      endif   
    endif
    
    tmp=ana%alpha
    allocate(zt(1:ana%samp+2))
    zt(1)=ana%Zmin
    do i=1,ana%samp
      zt(i+1)=ana%z(i)
    enddo
    zt(ana%samp+2)=ana%Zmax  
    
    call search_value(zt,z,ana%samp+2,ind)
    
    if(ind==ana%samp+1)then
      b=ana%Zmax-ana%phi(ind)*ana%Gmax
    elseif(ind==ana%samp)then  
      b=ana%z(ind)-ana%phi(ind)*ana%y(ind)    
    elseif(ind==1)then
      b=zt(2)-ana%phi(1)*ana%y(1)
    else
      b=ana%z(ind)-ana%phi(ind)*(tmp*ana%y(ind+1)+(1.-tmp)*ana%y(ind))
    endif
    
    y=(z-b)/ana%phi(ind)

    if(ana%m==1)then !eho2605
      if(ind==ana%samp+1)then
        call ana_rtail_inv(ana,y,z)
        !print*,'rtail ana_exp_inv',trim(ana%id),z
        !print*,'rtail ana_exp_inv',y,ana%y(ana%samp)
      elseif(ind==1)then
        call ana_ltail_inv(ana,y,z)
        !print*,'ltail ana_exp_inv',trim(ana%id),z
        !print*,'ltail ana_exp_inv',y,ana%y(1)  
      endif 
    endif
      
    deallocate(zt)
    
    
    endif !ana%chl=log
!LOG!
#endif    
  
  end subroutine

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   subroutine ana_chi_exp_obs(obs,nrobs,ana)
   !real exp:coeff_var=1. until 1999 011!
   !real exp: coeff_var=0.49 from 1999 019 to 1999 078!
   
   use mod_observations
   implicit none
   integer :: nrobs       ! Number of measurements
   type(observations) :: obs(nrobs) 
   type(anamorphosis)::ana
   integer::k
   real::d

   do k=1,nrobs
      call ana_exp_inv_eco22(obs(k)%d,d,ana)
      obs(k)%d=d
       !obs(k)%var=0.09!(0.1*obs(k)%d)**2!(0.3*obs(k)%d)**2!0.004!abs(d*0.015)
   enddo      
       
   end subroutine 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
      subroutine ana_varobs(obs,ana,nrens)
      !subroutine ana_varobs(obs,nrobs,ana,nrens,nx,ny,meandx)
      use mod_observations
      use m_random
      implicit none
      integer::nrobs,nrens!,nx,ny
      type(observations):: obs
      type(anamorphosis)::ana
      real::varo!,meandx
      integer::k,iens
      real::md,d
      real,dimension(nrens)::work,work2
      
      !obs(k)%var=1.!0.09
      obs%var=0.0818!0.0625!0.25
       !test: obs< seawif limits!
      !if((obs%d.le.ana%z(2)).or.(obs%d.ge.ana%z(ana%samp-1)))then
      !  print*,'dsig: obs out of xf-range'
      !  obs%var=1.
      !endif
      !if(obs%d.le.0.01)then
      !  print*,'dsig: obs lower than 0.01'
      !  obs%var=4.
      !endif
      
      call randn(nrens,work2)
      do iens=1,nrens
        if(obs%d==0.)then
          print*,'obs null',iens
          obs%d=1.e-8
        endif
        work(iens)=obs%d*exp(sqrt(obs%var)*&
	                 (work2(iens)-sqrt(obs%var)/2.))
	
        call ana_exp_inv_eco22(work(iens),d,ana)
        work(iens)=d	
      enddo	
      call ana_scatter(work,nrens,'obs')
      
      md=sum(work)/real(nrens)
      varo=0.
      do iens=1,nrens
        varo=varo+(work(iens)-md)**2
      enddo	
      obs%var=varo/real(nrens-1)
      !log!eho070212
      obs%var=0.09
      print*,'dsig',md,sqrt(obs%var)
      
      print*,'avant',obs%d
      call ana_exp_inv_eco22(obs%d,d,ana)
      obs%d=d
      print*,'apres',obs%d

      end  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   subroutine ana_chi_exp_S(S,nrobs,nrens,iens,ana)   
   implicit none
   integer :: nrobs,nrens,iens        ! Number of measurements
   real, dimension(1:nrobs,1:nrens) :: S
   type(anamorphosis)::ana
   integer::k
   real:: d
     
   if(ana%obs)then
      do k=1,nrobs
       call ana_exp_inv_eco22(S(k,iens),d,ana) 
       S(k,iens)=d
      enddo 
   else
      write(*,*)'ana_chi_exp_S: you are not in the observation space'
   endif  
       
   end subroutine  
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   
   subroutine ana_phi_exp(fld,nz,ana)   
   implicit none
   integer :: nz      
   type(anamorphosis)::ana
   real, dimension(1:nz) :: fld 
   integer :: k
   real:: d
  
   do k=1,nz
      call ana_exp_inv_eco22(fld(k),d,ana)       
      fld(k)=d
   enddo
        
   end subroutine
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
   subroutine ana_invphi_exp(fld,nz,ana)   
   implicit none
   integer :: nz     ! Number of measurements
   type(anamorphosis)::ana
   real, dimension(1:nz) :: fld
   integer :: k
   real ::d
   
   do k=1,nz
      call ana_exp_eco22(fld(k),d,ana)       
      fld(k)=d
   enddo
       
   end subroutine

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine test_ana(fld,nz,ana,mode)   
   implicit none
   integer :: nz,ka,mode     ! Number of measurements
   type(anamorphosis)::ana
   real, dimension(1:nz) :: fld
   integer :: k
   real ::d,d2
   character(len=80) :: file
   
   file='test_ana_'//trim(ana%id)//'.txt'
   open(unit=107,FILE=file,FORM='formatted',STATUS='unknown',position='append')
   
   if(mode==0)then   
      do k=1,nz
        write(107,*)'test_ana: avant ',fld(k)
        call ana_exp_inv_eco22(fld(k),d,ana) 
        write(107,*)'test_ana: gauss ',d
        call ana_exp_eco22(d,d2,ana) 
        write(107,*)'test_ana: apres ',d2
      enddo
    else
      do k=1,nz
        write(107,*)'test_ana: avant ',fld(k)
        call ana_exp_eco22(fld(k),d,ana) 
        write(107,*)'test_ana: reel ',d
        call ana_exp_inv_eco22(d,d2,ana) 
        write(107,*)'test_ana: apres ',d2
      enddo    
    endif
    
    
   close(107) 
       
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

   /* subroutine anamorphosis_stat(m,ana,ixp)
    use mod_dimensions
    implicit none
    integer:: m,ixp
    type(anamorphosis), dimension(:),allocatable::ana  
    !real,dimension(99*kdim)::zt
    real,dimension(:),allocatable::zt
    integer::k
    
    !reading of the data which will be anamorphosed: analysisfields_ana.in!
    call get_analysisfields_ana()
   
    !allocation of the anamorphosis structure!
    allocate(ana(numfields_ana))
    
    do k=1,numfields_ana
      ana(k)%id=trim(fieldnames_ana(k))
      print *,'----------   ',ana(k)%id,'   ----------'
      
      !reading of the smoothing param of the anamorphosis!
      call ana_smooth_param(ana(k))

      if(ana(k)%obs)then
        !ana(k)%n=nrens*kdim+nrmes_z
        ana(k)%n=90*9*nz+obstot*(floor(45.*3600.*24./(dt*real(nstep)))+1)
        ana(k)%nsamp=9
      else
        if(isstate(ana(k)%id,biomodel_ana))then
          ana(k)%n=90*9*nz
          ana(k)%nsamp=9
        else !parameters!
          ana(k)%n=nrens
	  ana(k)%nsamp=2
        endif
      endif
      allocate(zt(ana(k)%n))
      
      if(mod(ana(k)%n-2,ana(k)%nsamp)==0)then
        ana(k)%samp=(ana(k)%n-2)/ana(k)%nsamp+1
      else
        ana(k)%samp=(ana(k)%n-2)/ana(k)%nsamp+2
      endif
      allocate(ana(k)%z(ana(k)%samp+1))
      
      !ana(k)%nsamp=ana(k)%n/(ana(k)%samp-1)
      print*,'ana nsamp',ana(k)%nsamp,ana(k)%samp
      
      !reading and sorting of the data!
      !call ana_data_z_stat(ana(k),zt,m,9,40,ana(k)%n,ixp)
      !call ana_data_z_stat(ana(k),zt,m,4,91,ana(k)%n,ixp)!4 years,4 days
       call ana_data_z_stat(ana(k),zt,m,9,91,ana(k)%n,ixp)!9 years , 4days
      !computation of yi=invG(i/n)!
      print *,'gauss Y'
      call ana_gauss_y(ana(k))
      
      !smooth anamorphosis!
      print *,'phi'
      call ana_cpt_phi_eco1d2(ana(k),zt,ana(k)%n)   
      
      deallocate(zt)
    enddo
    
    end*/
    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               
    subroutine anamorphosis_dyn(ccda,nrens,nz,nspec,nrobs,ana,obs,m0,dt,nstep)
    use m_bio_estimation
    use mod_observations
    implicit none
    integer:: nrens,m0,nspec,nz,nstep,nrobs
    real,dimension(1:nspec,0:nz,1:nrens)::ccda
    type(anamorphosis), dimension(:),allocatable::ana  
    !real,dimension(nrens*kdim)::zt
    real::dt
    real,dimension(:),allocatable::zt
    integer::k,m
    type(observations), intent(in) :: obs(nrobs)
    type(observations), allocatable :: d(:)
    integer::obstot
    integer::obschl
    
    !reading of the data which will be anamorphosed: analysisfields_ana.in!
    call get_analysisfields_ana()
   
    !allocation of the anamorphosis structure!
    allocate(ana(numfields_ana))
    
    !observation!
    obstot=0
    obschl=0
    do m=1,nrobs
      if (obs(m)%action) obstot=obstot+1
      if (trim(obs(m)%ch)=='chl') obschl=obschl+1
    enddo
    allocate(d(obstot))
    obstot=0
    do m=1,nrobs
      if (obs(m)%action) then
         obstot=obstot+1
         d(obstot)=obs(m)
      endif 
    enddo
    
    do k=1,numfields_ana
      ana(k)%id=trim(fieldnames_ana(k))
      print *,'----------   ',ana(k)%id,'   ----------'
      
      !reading of the smoothing param of the anamorphosis!
      call ana_smooth_param(ana(k))
      if(trim(ana(k)%id)=='chl')then
        if(m0.gt.30)then
          ana(k)%n=nrens*nzchl+obschl*30!45!(floor(45.*3600.*24./(dt*real(nstep)))+1)
        else
	  ana(k)%n=nrens*nzchl+obschl*m0
	endif
	ana(k)%nsamp=3    
      elseif(ana(k)%obs)then
        !ana(k)%n=nrens*kdim+nrmes_z
        ana(k)%n=nrens*nz+(obstot-obschl)*30!45!(floor(45.*3600.*24./(dt*real(nstep)))+1)
        ana(k)%nsamp=9         
      else
        if(isstate(ana(k)%id,biomodel_ana))then
          ana(k)%n=nrens*nz
          ana(k)%nsamp=9
        else !parameters!
          ana(k)%n=nrens
	  ana(k)%nsamp=2
        endif
      endif
      allocate(zt(ana(k)%n))
      
!      if((trim(ana(k)%id)=='P').or.(trim(ana(k)%id)=='N')&
!                                 .or.(trim(ana(k)%id)=='H'))then     
        !allocate(ana(k)%z(nana_dyn))
!        allocate(ana(k)%z(nana_dyn+1))!eco1d2
!	ana(k)%samp=nana_dyn
!      else
        !allocate(ana(k)%z(nana_par))
!        allocate(ana(k)%z(nana_par+1))!eco1d2
!	ana(k)%samp=nana_par
!      endif
      
      !ana(k)%nsamp=ana(k)%n/(ana(k)%samp-1)      
      if(mod(ana(k)%n-2,ana(k)%nsamp)==0)then
        ana(k)%samp=(ana(k)%n-2)/ana(k)%nsamp+1
      else
        ana(k)%samp=(ana(k)%n-2)/ana(k)%nsamp+2
      endif
      allocate(ana(k)%z(ana(k)%samp))
      
      !reading and sorting of the data!
      call ana_data_z_dyn(ana(k),zt,nrens,nz,nspec,ccda,d,obstot,ana(k)%n,m0,dt,nstep)
   
      !computation of yi=invG(i/n)!
      print *,'gauss Y'
      call ana_gauss_y(ana(k))
      
      !smooth anamorphosis!
      print *,'phi'
      call ana_cpt_phi_eco2(ana(k),zt,ana(k)%n)   
      
      deallocate(zt)
    enddo
    
    end   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine anamorphosis_partial_dyn(ccda,nrens,nz,nspec,nrobs,ana,obs,m0,dt,nstep)
    use m_bio_estimation
    use mod_observations
    use netcdf
    implicit none
    integer:: nrens,m0,nspec,nz,nstep,nrobs
    real,dimension(1:nspec,0:nz,1:nrens)::ccda
    type(anamorphosis), dimension(:),allocatable::ana  
    !real,dimension(nrens*kdim)::zt
    real::dt
    real,dimension(:),allocatable::zt
    integer::k,m
    type(observations), intent(in) :: obs(nrobs)
    type(observations), allocatable :: d(:)
    integer::obstot
    integer::obschl
    integer::ierr,ntime,varid,ncid
    
    !reading of the data which will be anamorphosed: analysisfields_ana.in!
    call get_analysisfields_ana()
   
    !allocation of the anamorphosis structure!
    allocate(ana(numfields_ana))
    
    !observation!
    obstot=0
    obschl=0
    do m=1,nrobs
      if (obs(m)%action) obstot=obstot+1
      if (trim(obs(m)%ch)=='chl') obschl=obschl+1
    enddo
    allocate(d(obstot))
    obstot=0
    do m=1,nrobs
      if (obs(m)%action) then
         obstot=obstot+1
         d(obstot)=obs(m)
      endif 
    enddo
    
    do k=1,numfields_ana
      ana(k)%id=trim(fieldnames_ana(k))
      print *,'----------   ',ana(k)%id,'   ----------'
      
      !reading of the smoothing param of the anamorphosis!
      call ana_smooth_param(ana(k))
      
      if(trim(ana(k)%id)=='chl')then      
        !ana(k)%n=9*27*kdim+41*nrmes_z !obs every 4 days
        !ana(k)%nsamp=9
	ierr=nf90_open('observation.nc',nf90_nowrite,ncid)
        ierr=nf90_inq_dimid(ncid,'time',varid)
        if(ierr/=nf90_noerr)stop
        ierr=nf90_inquire_dimension(ncid,varid,len=ntime) 
        if(ierr/=nf90_noerr)stop
	ierr=nf90_close(ncid)
	
        if(m0.gt.30)then
          !ana(k)%n=nens_stat*61*nzchl+obschl*30
          ana(k)%n=nens_stat*61*nzchl+obschl*(60+min(29,ntime-m0+1))!obs2
	else
	  !ana(k)%n=nens_stat*61*nzchl+obschl*m0
          ana(k)%n=nens_stat*61*nzchl+obschl*(2*m0+29) !obs2
	endif
        ana(k)%nsamp=10!3    
      elseif(ana(k)%obs)then
        !ana(k)%n=nrens*kdim+nrmes_z
        ana(k)%n=nens_stat*61*nz+(obstot-obschl)*30!45
        ana(k)%nsamp=50!9       
      else
        if(isstate(ana(k)%id,biomodel_ana))then
          ana(k)%n=nrens*nz
          ana(k)%nsamp=50!9
        else !parameters!
          ana(k)%n=nrens
	  ana(k)%nsamp=5
        endif
      endif
      allocate(zt(ana(k)%n))
      
      if(mod(ana(k)%n-2,ana(k)%nsamp)==0)then
        ana(k)%samp=(ana(k)%n-2)/ana(k)%nsamp+1 !eho21/12
      else
        ana(k)%samp=(ana(k)%n-2)/ana(k)%nsamp+2
      endif
      allocate(ana(k)%z(ana(k)%samp))
      
      !reading and sorting of the data!
      if(ana(k)%obs)then
	call ana_data_z_stat(ana(k),zt,nrens,nz,d,obstot,ana(k)%n,m0)  
      else
        call ana_data_z_dyn(ana(k),zt,nrens,nz,nspec,ccda,d,obstot,ana(k)%n,m0,dt,nstep)
      endif
      
      !computation of yi=invG(i/n)!
      print *,'gauss Y'
      call ana_gauss_y(ana(k))
      
      !smooth anamorphosis!
      print *,'phi'
      call ana_cpt_phi_eco3(ana(k),zt,ana(k)%n)   
      
      deallocate(zt)
    enddo
    
    end   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ana_ltail(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      real, parameter :: pi=3.1415927
     
      if(y.ge.0.) then
        print*,'alert ltail'
	stop
      endif     
      
      d1=ana%phi(1)/ana%z(1)
      z=ana%z(1)*exp(d1*(y-ana%y(1)))
           
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      
      subroutine ana_ltail_inv(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      real, parameter :: pi=3.1415927
      
      if(z==0.) then
        print*,'alert ltail_inv'
	stop
      endif           
      
      d1=ana%phi(1)/ana%z(1)
      y=log(z/ana%z(1))/d1+ana%y(1)
      
      if(y.ge.0)then
        print*,'alert ltail_inv: y pos',trim(ana%id)
	stop
      endif
      
           
      end
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ana_rtail(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      real, parameter :: pi=3.1415927
      
      if(y.le.0.) then
        print*,'alert rtail',trim(ana%id)
	stop
      endif

      
      if(ana%z(ana%samp).ge.ana%Zmax)then
        print*,'alert rtail: z >Zmax ',trim(ana%id)
	stop	
      endif
      
      d1=ana%phi(ana%samp+1)/(ana%Zmax-ana%z(ana%samp))
      z=(ana%Zmax-(ana%Zmax-ana%z(ana%samp))*exp(d1*(ana%y(ana%samp)-y)))
      
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      
      subroutine ana_rtail_inv(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      real, parameter :: pi=3.1415927
      
      if(ana%z(ana%samp)==0.)then
        print*,'alert rtail_inv'
	stop
      endif      
      
      if(ana%z(ana%samp).ge.ana%Zmax)then
        print*,'alert rtail_inv: y neg ',trim(ana%id)      
	stop	
      endif
      
      d1=ana%phi(ana%samp+1)/(ana%Zmax-ana%z(ana%samp))
      y=ana%y(ana%samp)-log((ana%Zmax-z)/(ana%Zmax-ana%z(ana%samp)))/d1
                              
      if(y.le.0.)then
        print*,'alert rtail_inv: y neg 2 ',trim(ana%id)
	stop
      endif
           
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      subroutine ana_scatter(work,kdim,char3)
      implicit none
      integer::kdim
      real,dimension(1:kdim)::work
      character(len=3)::char3
      character(len=80)::memfile
      integer::k
  
      memfile='scatter_'//trim(char3)//'.txt'
      open(unit=420,FILE=trim(memfile),FORM='formatted',STATUS='unknown')          
      do k=1,kdim
        write(420,*)work(k)
      enddo

      close(420)    

      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! interaction with biomodel     !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    subroutine ana_get_id_parmodel(parname,cfld)
    implicit none 
    character(len=8)::parname
    character(len=3)::cfld
    
    if(trim(biomodel_ana)=='NORWECOM')then
      if(parname(1:3)=='phi')then
        select case(trim(parname))
	  case('phimes1')
	    cfld='ms1'
	  case('phimes2')
	    cfld='ms2' 
	  case('phimes3')
	    cfld='ms3' 
	  case('phimic1')
	    cfld='mc1' 
	  case('phimic2')
	    cfld='mc2' 
	  case('phimic3')
	    cfld='mc3' 
	  case default
	    print*,'ana: unknown preferences ', trim(parname) 
	end select
      else
        cfld=parname(1:3)
      endif
      !print*,'ana_get_id_parmodel: ',trim(parname),cfld
    else
      print*,'not implemented'
      stop
    endif
        
    end subroutine
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
         
    subroutine ana_get_name_parmodel(cfld,parname)
    implicit none 
    character(len=8)::parname
    character(len=3)::cfld
    
    if(trim(biomodel_ana)=='NORWECOM')then
      select case(trim(cfld))
	case('ms1')
	  parname='phimes1'
	case('ms2')
	  parname='phimes2'
	case('ms3')
	  parname='phimes3'
	case('mc1')
	  parname='phimic1'
	case('mc2')
	  parname='phimic2'
	case('mc3')
	  parname='phimic3'
	case default
	  parname=trim(cfld)
      end select
    else
      print*,'not implemented'
      stop
    endif
        
    end subroutine
     
      
   END MODULE m_anamorphosis

!-----------------------------------------------------------------------

!Copyright (C) 2000 - NAME
