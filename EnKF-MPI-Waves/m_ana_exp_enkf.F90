module m_ana_exp_enkf
!Ehouarn!
!module including the experimental anamorphosis routines for the EnKF! 
      type anamorphosis
        character(len=3) smooth !nature of smoothness of the empirical anamorphosis!
        character(len=8) id !nature of the physical variable!
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
   
      character(len=*), parameter :: infile_ana='analysisfields_ana.in'
      integer,save :: numfields_ana
      character(len=8), dimension(:), save, allocatable:: fieldnames_ana   
      type(anamorphosis),dimension(:),save,allocatable::ana_enkf
      integer,dimension(:),save, allocatable::obs_bias
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!         intrinsec functions     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      integer function ana_get_ind(ana,nana,char3) 
      implicit none
      integer :: nana
      type(anamorphosis),dimension(nana)::ana
      character(len=8) :: char3     
      integer::k

      do k=1,nana
        if(trim(ana(k)%id)==trim(char3))then
          ana_get_ind=k
          return
       endif  
      enddo
      ana_get_ind=-1
   
      end function
      
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      subroutine init_ana_exp_eco()     
      implicit none
      integer::k,samp,j
      character(len=80) :: file
      real::y,tmp
   
      !reading analysfields_ana.in!
      call get_analysisfields_ana()
   
      allocate(ana_enkf(numfields_ana))
   
      open(unit=50,FILE='ana_size.txt',FORM='formatted',STATUS='unknown')
   
      do k=1,numfields_ana
        ana_enkf(k)%id=fieldnames_ana(k)
        call ana_smooth_param(ana_enkf(k))
        read(50,*)samp
        ana_enkf(k)%samp=samp
        allocate(ana_enkf(k)%z(samp))
        allocate(ana_enkf(k)%y(samp))
        allocate(ana_enkf(k)%phi(samp+1))
        file='ana_'//trim(ana_enkf(k)%id)//'_phi.dat'
        open(unit=11,FILE=file,FORM='unformatted',STATUS='unknown')
        read(11)ana_enkf(k)%phi
        close(11)
        file='ana_'//trim(ana_enkf(k)%id)//'_y.dat'
        open(unit=21,FILE=file,FORM='unformatted',STATUS='unknown')
        read(21)ana_enkf(k)%y
        close(21)
        file='ana_'//trim(ana_enkf(k)%id)//'_zs.dat'
        open(unit=31,FILE=file,FORM='unformatted',STATUS='unknown')
        read(31)ana_enkf(k)%z
        close(31)

        ana_enkf(k)%Gmax=(ana_enkf(k)%Zmax-ana_enkf(k)%z(samp))&
                              /ana_enkf(k)%phi(samp+1)+ana_enkf(k)%y(samp)!y!passe par Y(n-1)
        
	ana_enkf(k)%Gmin=ana_enkf(k)%y(1)-(ana_enkf(k)%z(1)-ana_enkf(k)%Zmin)/ana_enkf(k)%phi(1)
	print *,'GMAX ',ana_enkf(k)%id,ana_enkf(k)%Gmax
	print *,'Gmin ',ana_enkf(k)%id,ana_enkf(k)%Gmin
	
      enddo
      close(50)
     
      end subroutine
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!         reading of analysisfields_ana.in     !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function get_nrfields_ana()
#if defined (QMPI)
      use qmpi
#else
      use qmpi_fake
#endif   
      implicit none
      integer :: ios,first,last
      logical :: ex
      character(len=8) :: char3

      inquire(exist=ex,file=infile_ana)
      if (.not. ex) then
        if (master) print *,'Could not find '//infile_ana
        call stop_mpi()
      end if

      open(10,status='old',form='formatted',file=infile_ana)
      ios=0
      get_nrfields_ana=0
      do while (ios==0)
        read(10,100,iostat=ios) char3
        if (ios==0) get_nrfields_ana=get_nrfields_ana+1
      end do
      close(10)
  100 format (a8)
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
   
      subroutine get_analysisfields_ana()
#if defined (QMPI)
      use qmpi
#else
      use qmpi_fake
#endif   
      implicit none
      integer :: first,last,k,nfld,ios
      logical :: ex
      character(len=8) :: char3

      numfields_ana=get_nrfields_ana()
      if (master) print *,'numfields_ana is ',numfields_ana
      if (numfields_ana<=0 .or.numfields_ana > 20) then !
         if (master) print *,'numfields_ana is higher than max allowed setting or = 0'
         call stop_mpi()
      end if
      allocate(fieldnames_ana(numfields_ana))

      inquire(exist=ex,file=infile_ana)
      if (.not. ex) then
        if (master) print *,'Could not find '//infile_ana
        call stop_mpi()
      end if

      open(10,status='old',form='formatted',file=infile_ana)
      ios=0
      nfld=0
      do while (ios==0)
        read(10,100,iostat=ios) char3
        if (ios==0) then
          fieldnames_ana (nfld+1)=char3
          nfld=nfld+1
        end if
      end do
      close(10)
  100 format (a8)

      if (nfld/=numfields_ana) then
        if (master) print *,'An error occured when reading '//infile_ana
        call stop_mpi()
      end if

      ! List fields used in analysis
      do k=1,numfields_ana
        if (master) print *,fieldnames_ana(k)
      end do

      end subroutine
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !!!!!        allocation of the structure    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      subroutine ana_size_sample(ana,nx,ny,nz,nens,depths)  
      use m_get_micom_fld
      implicit none
      type(anamorphosis)::ana  
      integer::nx,ny,nz, nens
      integer::i,j,k,m,ks,l,ks2
      real, dimension(nx,ny):: fld_dp,depths,pb
      character(len=3) :: cmem
      character(len=80) :: memfile
      character(len=8) :: cfield
      real,dimension(:),allocatable::phi
      integer,dimension(1:nz)::assign
      integer::dpnull
      real, dimension(:,:,:),allocatable::fld3d
      real::d
      
      ana%n=0
      if (ana%obs)then
        !values from the observations (could be merged with the forecast, let see..)     
        print*,'Twin experiments: forecast at the surface'
	cfield='dp'
	do m=1,nens
	  write(cmem,'(i3.3)') m
          memfile='forecast'//cmem 
	  call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 1, 1, nx, ny)
            do j=1,ny
              do i=1,nx
	        if((fld_dp(i,j).gt.0.).and.(fld_dp(i,j).lt.1.e10))then
	        !to be tuned..
	          ana%n=ana%n+1
	        endif
	      enddo
            enddo
	enddo
      else
	!values for the forecast ensemble    
        if(trim(ana%id)=='dpphi')then
!	  do j=1,ny
!	  do i=1,nx
!	    if(depths(i,j).gt.0.) then
!	      ana%n=ana%n+(nz-1)*nens
!	    endif	
!          enddo	
!	  enddo 
         allocate(fld3d(nx,ny,nz))
         allocate(phi(nz-1))
      
         do m=1,nens
           write(cmem,'(i3.3)') m
           memfile='forecast'//cmem 
       
           do l=1,2
	   
	   cfield='pb'
	   call get_micom_fld_new(trim(memfile), pb, m, cfield,&
                 l, 1, nx, ny)
	   do k=1+(l-1)*nz,nz+(l-1)*nz
	     if(k.gt.nz)then
	       ks=k-nz
	     else
	       ks=k
	     endif	     
	     
!	     if(k.lt.nz)then
!	       ks2=k+nz
!	     else
!	       ks2=k
!	     endif
	     
	     cfield='dp'        
             call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 k, 1, nx, ny)		     
	     
	     do j=1,ny
	     do i=1,nx
	       if(depths(i,j).gt.0.) then
	         fld3d(i,j,ks)=abs(fld_dp(i,j))/pb(i,j)
	       endif	
	     enddo
	     enddo 
	   enddo
	
	   !numerical corrections: sum(fld3d,3)=1
	   do j=1,ny
	   do i=1,nx
	     if(depths(i,j).gt.0.)then
	       d=sum(fld3d(i,j,:)) 
	       fld3d(i,j,:)=fld3d(i,j,:)/d
	     endif
	   enddo
	   enddo

	
	   do j=1,ny
	   do i=1,nx
	     if(depths(i,j).gt.0.) then
	       call ana_pi2phi(nz,fld3d(i,j,:),phi,assign,dpnull)
	       do k=1,nz-1	       
	         if((phi(k).gt.1.e-6).and.(phi(k).lt.(1.-1.e-6)))then
		 !if((phi(k).gt.1.e-8).and.(phi(k).lt.0.99999999))then
	          ana%n=ana%n+1
		 endif  
	       enddo
	     endif	
	   enddo
           enddo	 	
           
	   enddo !l=1,2
	   
	  enddo
          deallocate(fld3d,phi)

        else
	cfield='dp'
        do m=1,nens
          write(cmem,'(i3.3)') m
          memfile='forecast'//cmem 
          do k=1,nz           
	    call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 k, 1, nx, ny)            
	    do j=1,ny
              do i=1,nx
	        if((k.ge.3).and.(fld_dp(i,j).ge.9806.).and.(fld_dp(i,j).lt.1.e10))then
	        !to be tuned..
	          ana%n=ana%n+1
	        endif
	      enddo
            enddo	        
          enddo       
        enddo   
        endif
      
      endif
           
      write(*,*) trim(ana%id), ' - size of the sample: ',ana%n
   
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!     physical data    !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   
  subroutine ana_data_z(ana,zt, nx,ny,nz,nens,nsize,depths)  
  use m_get_micom_fld
  use nfw_mod
  implicit none
  type(anamorphosis)::ana
  integer::nx,ny,nz,nens,nsize
  real,dimension(1:nsize)::zt
  integer::m,cptz,i,j,k,ks,l,ks2
  real, dimension(nx,ny):: fld_dp,fld,pb, depths
  real, dimension(:,:,:),allocatable::fld3d!,fldtmp
  character(len=3) :: cind
  character(len=3) :: cmem
  character(len=80) :: memfile,memfilephy
  character(len=8) :: cfield
  logical::ex
  integer,dimension(nx,ny,2*nz)::assign
  real, dimension(nx,ny,2*(nz-1))::phi
  integer,dimension(nx,ny,2)::dpnull
  real::d
  integer::vid_phi,vid_a,ncid,dimnx,dimny,dimnz,dimnzm1,vid_depths,vid_dp1,dimndt
  integer,dimension(nz)::atmp
  real, dimension(nz-1)::ptmp
  integer::dptmp
  
  cptz=1
  
  if (ana%obs)then
    !values from the observations (could be merged with the forecast, let see..)     
    print*,'Twin experiments: forecast at the surface'
    do m=1,nens
      write(cmem,'(i3.3)') m
      memfile='forecast'//cmem 
      cfield='dp'
      call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 1, 1, nx, ny)
      cfield=trim(ana%id)	
      call get_micom_fld_new(trim(memfile), fld, m, cfield,&
                 1, 1, nx, ny)	  
      do j=1,ny
        do i=1,nx
	  if((fld_dp(i,j).gt.0.).and.(fld_dp(i,j).lt.1.e10))then
	     zt(cptz)=fld(i,j)
	      cptz=cptz+1
	   endif
	 enddo
       enddo
    enddo
  else
    !values from the forecast ensemble    

    if(trim(ana%id)=='dpphi')then
      !special case: needs to get the wohle water column before transformation!
      allocate(fld3d(nx,ny,nz))
      !allocate(fldtmp(nx,ny,2*nz))
      
      do m=1,nens
        write(cmem,'(i3.3)') m
        memfile='forecast'//cmem 
       
        phi(:,:,:)=0.
        assign(:,:,:)=0
        !fldtmp(:,:,:)=0.
        dpnull(:,:,:)=0
       
        do l=1,2
       
        cfield='pb'
	call get_micom_fld_new(trim(memfile), pb, m, cfield,&
                 l, 1, nx, ny)
	do k=1+(l-1)*nz,nz+(l-1)*nz
	  if(k.gt.nz)then
	     ks=k-nz
	   else
	     ks=k
	  endif
	  
!	  if(k.lt.nz)then
!	     ks2=k+nz
!	   else
!	     ks2=k
!	  endif
	  
	  cfield='dp'        
          call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 k, 1, nx, ny)	 

	  do j=1,ny
	    do i=1,nx
	      if(depths(i,j).gt.0.) then
	        fld3d(i,j,ks)=abs(fld_dp(i,j))/pb(i,j)
		!fldtmp(i,j,k)=fld_dp(i,j)
	      endif	
	    enddo
	  enddo 
	enddo
	
	!numerical corrections: sum(fld3d,3)=1
	do j=1,ny
	do i=1,nx
	  if(depths(i,j).gt.0.)then
	    d=sum(fld3d(i,j,:)) 
	    fld3d(i,j,:)=fld3d(i,j,:)/d
	  endif
	enddo
	enddo

	do j=1,ny
	  do i=1,nx
	    if(depths(i,j).gt.0.) then
	      call ana_pi2phi(nz,fld3d(i,j,:),ptmp,atmp,dptmp)
	      phi(i,j,1+(l-1)*(nz-1):nz-1+(l-1)*(nz-1))=ptmp(1:nz-1)
	      assign(i,j,1+(l-1)*nz:nz+(l-1)*nz)=atmp(1:nz)
	      dpnull(i,j,l)=dptmp
!	      call ana_pi2phi(nz,fld3d(i,j,:),phi(i,j,1+(l-1)*(nz-1):nz-1+(l-1)*(nz-1))&	                         
!			           	  ,assign(i,j,1+(l-1)*nz:nz+(l-1)*nz))
	      do k=1+(l-1)*(nz-1),nz-1+(l-1)*(nz-1)       
	        !if(k.gt.nz)then
	        !ks=k-nz
	        !else
	        ! ks=k
	        !endif
		!fldtmp(i,j,k)=fld3d(i,j,ks)
		if((phi(i,j,k).gt.1.e-6).and.(phi(i,j,k).lt.(1.-1.e-6)))then
		!if((phi(i,j,k).gt.1.e-8).and.(phi(i,j,k).lt.0.99999999))then
	   	  zt(cptz)=phi(i,j,k)
	          cptz=cptz+1
		endif  
	      enddo
	    endif	
	  enddo
        enddo	 	
        
	enddo !l=1,2
	
        write(cmem,'(i3.3)') m
        memfilephy='forecast_phi'//cmem//'.nc' 
        call nfw_create(memfilephy, nf_write, ncid)

        call nfw_def_dim(memfilephy, ncid, 'nx', nx, dimnx)
        call nfw_def_dim(memfilephy, ncid, 'ny', ny, dimny)
        call nfw_def_dim(memfilephy, ncid, 'nz', 2*nz, dimnz)
        call nfw_def_dim(memfilephy, ncid, 'nphi', 2*(nz-1), dimnzm1)
	call nfw_def_dim(memfilephy, ncid, 'ndt', 2, dimndt)

        call nfw_def_var(memfilephy, ncid, 'dpphi', nf_double, 3,(/dimnx,dimny,dimnzm1/), vid_phi)
        call nfw_def_var(memfilephy, ncid, 'assign', nf_int, 3,(/dimnx,dimny,dimnz/), vid_a)
        call nfw_def_var(memfilephy, ncid, 'depths', nf_double, 2,(/dimnx,dimny/), vid_depths)
	!call nfw_def_var(memfilephy, ncid, 'dp1', nf_double, 3,(/dimnx,dimny,dimnz/), vid_dp1)
	call nfw_def_var(memfilephy, ncid, 'dpnull', nf_int, 3,(/dimnx,dimny,dimndt/), vid_dp1)
	
        call nfw_enddef(memfilephy, ncid)

        call nfw_put_var_double(memfilephy, ncid, vid_phi, phi)
        call nfw_put_var_int(memfilephy, ncid, vid_a, assign)
        call nfw_put_var_double(memfilephy, ncid, vid_depths, depths)
	!call nfw_put_var_double(memfilephy, ncid, vid_dp1, fldtmp)
	call nfw_put_var_int(memfilephy, ncid, vid_dp1, dpnull)
	
        call nfw_close(memfilephy, ncid)
     
      enddo
      deallocate(fld3d)

    else
      do m=1,nens
        write(cmem,'(i3.3)') m
        memfile='forecast'//cmem 
 
        do k=1,nz  
	  cfield='dp'        
          call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 k, 1, nx, ny)		
	  cfield=trim(ana%id)	
          call get_micom_fld_new(trim(memfile), fld, m, cfield,&
                 k, 1, nx, ny) 	
	  do j=1,ny
            do i=1,nx
	      if((k.ge.3).and.(fld_dp(i,j).ge.9806.).and.(fld_dp(i,j).lt.1.e10))then
	      !to be tuned..
	        zt(cptz)=fld(i,j)
	        cptz=cptz+1
	      endif
	    enddo
          enddo 	         
        enddo       
      enddo     
    endif
  
  endif
  end subroutine


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
           !!!!!    smoothness param   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   
   subroutine ana_smooth_param(ana)
#if defined (QMPI)
   use qmpi
#else
   use qmpi_fake
#endif      
   !use qmpi
   !use qmpi_fake
   implicit none
   type(anamorphosis)::ana
   character(len=80) ::memfile
   logical::ex
   
   memfile='infile_ana.'//trim(ana%id)//'.in'
   
   inquire(exist=ex,file=memfile)
   if (.not. ex) then
      print *,'Could not find '//memfile
      call stop_mpi()
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
	  if(k0==ana%samp+2)then
	    !no larger values..    
	    if(k==1)then
	      print*,'alert: anamorposis null ',trim(ana%id)
	      stop
	    else
	      !k>=3	     
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
	        ana%Gmin=ana%y(1)-(ana%z(1)-ana%Zmin)/ana%phi(1)
	      endif
	    
	    else
	      y=tmp*ana%y(k0+1)+(1.-tmp)*ana%y(k0)-tmp*ana%y(k)-(1.-tmp)*ana%y(k-1)
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
   
   file='ana_'//trim(ana%id)//'_phi.dat'
   open(unit=111,FILE=file,FORM='unformatted',STATUS='unknown')
   write(111)ana%phi
   close(111)
   file='ana_'//trim(ana%id)//'_y.dat'
   open(unit=121,FILE=file,FORM='unformatted',STATUS='unknown')
   write(121)ana%y
   close(121)
   file='ana_'//trim(ana%id)//'_zs.dat'
   open(unit=131,FILE=file,FORM='unformatted',STATUS='unknown')
   write(131)ana%z
   close(131)
   
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
   u=1.-0.3989423*u*exp(-0.5*x*x)
! then we give the general value of pdgauss
   out=u
   if (y < 0.) then
     out=1.-u
   endif
 
   end subroutine
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
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
   do while ((abs(x1-x0) > eps).or.(abs(y1-y0) > eps))
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
      elseif(ind==1)then
        call ana_ltail_inv(ana,y,z)
      endif 
    endif
      
    deallocate(zt)

!LOG!
#endif    
  
  end subroutine

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   subroutine ana_chi_exp_obs(obs,nrobs,ana)
   
   use mod_measurement
   implicit none
   integer :: nrobs       ! Number of measurements
   type(measurement) :: obs(nrobs) 
   type(anamorphosis)::ana
   integer::k
   real::d

   do k=1,nrobs
      call ana_exp_inv_eco22(obs(k)%d,d,ana)
      obs(k)%d=d
   enddo      
       
   end subroutine 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
      subroutine ana_varobs(obs,ana,nrens)
      !subroutine ana_varobs(obs,nrobs,ana,nrens,nx,ny,meandx)
      use mod_measurement
      use m_random
      implicit none
      integer::nrobs,nrens!,nx,ny
      type(measurement):: obs
      type(anamorphosis)::ana
      real::varo!,meandx
      integer::k,iens
      real::md,d
      real,dimension(nrens)::work,work2
      
      obs%var=0.0818
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
      !call ana_scatter(work,nrens,'obs')
      
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

      end subroutine
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
   
   subroutine ana_phi_exp_2D(fld,nx,ny,depth,ana)   
   implicit none
   integer :: nx,ny   
   type(anamorphosis)::ana
   real, dimension(nx,ny) :: fld,depth 
   integer :: i,j
   real:: d
  
    do j=1,ny
      do i=1,nx
        if(depth(i,j).gt.0.)then
          call ana_exp_inv_eco22(fld(i,j),d,ana)    
          fld(i,j)=d
	endif
      enddo
    enddo
        
   end subroutine
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
   subroutine ana_invphi_exp_2D(fld,nx,ny,depth,ana)   
   implicit none
   integer :: nx,ny   
   type(anamorphosis)::ana
   real, dimension(nx,ny) :: fld,depth
   integer :: i,j
   real ::d
   
    do j=1,ny
      do i=1,nx
        if(depth(i,j).gt.0.)then
          call ana_exp_eco22(fld(i,j),d,ana)       
          fld(i,j)=d
	endif  
      enddo
    enddo
       
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
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
   integer :: nz    
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
   
  subroutine nl_left_tail(y,z,ana)
  !gaussian left tail!
  implicit none
  type(anamorphosis)::ana
  real::y,z
   
  z=exp(-(y-ana%a1)**2/ana%a2)
        
  end subroutine 
   
  subroutine nl_left_tail_inv(z,y,ana)
  !gaussian left tail!
  implicit none
  type(anamorphosis)::ana
  real::y,z
  
  y=ana%a1-sqrt(-ana%a2*log(z))
        
  end subroutine 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ana_ltail(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      !real, parameter :: pi=3.1415927
     
      if(y.ge.0.) then
        print*,'alert ltail'
	stop
      endif

      d1=ana%phi(1)/ana%z(1)
      z=ana%z(1)*exp(d1*(y-ana%y(1)))
           
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      
      subroutine ana_ltail_inv(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      !real, parameter :: pi=3.1415927
      
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
      
           
      end subroutine
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ana_rtail(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      !real, parameter :: pi=3.1415927
      
      if(y.le.0.) then
        print*,'alert rtail'
	stop
      endif
      
      if(ana%z(ana%samp).ge.ana%Zmax)then
        print*,'alert rtail: y neg ',trim(ana%id)
	!z=ana%Zmax
	!return
	stop
      endif
      
      d1=ana%phi(ana%samp+1)/(ana%Zmax-ana%z(ana%samp))
      z=(ana%Zmax-(ana%Zmax-ana%z(ana%samp))*exp(d1*(ana%y(ana%samp)-y)))
      
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      
      subroutine ana_rtail_inv(ana,y,z)
      implicit none
      type(anamorphosis)::ana
      real::y,z
      real::d1,d2
      !real, parameter :: pi=3.1415927
      
      if(ana%z(ana%samp)==0.)then
        print*,'alert rtail_inv',trim(ana%id)
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
           
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine test_ana(fld,nx,ny,depth,ana,mode)   
   implicit none
   integer :: nx,ny,ka,mode     ! Number of measurements
   type(anamorphosis)::ana   
   real, dimension(nx,ny) :: fld, depth
   integer :: i,j
   real ::d,d2
   character(len=80) :: file
   
   file='test_ana_'//trim(ana%id)//'.txt'
   open(unit=10,FILE=file,FORM='formatted',STATUS='unknown',position='append')
   
   if(mode==0)then       
      do j=1,ny
        do i=1,nx
	  if (depth(i,j).gt.0.)then
          !write(10,*)'test_ana: avant ',fld(i,j)
          call ana_exp_inv_eco22(fld(i,j),d,ana) 
          !write(10,*)'test_ana: gauss ',d
          call ana_exp_eco22(d,d2,ana) 
          !write(10,*)'test_ana: apres ',d2
	  d=abs(d2-fld(i,j))
	  if(d.ge.1.e-8)then
	    write(10,*)fld(i,j),fld(i,j)-d2
	    write(10,*)'********************'
	  endif
	  endif
        enddo
      enddo	
    else
      do j=1,ny
        do i=1,nx
	  if (depth(i,j).gt.0.)then
          write(10,*)'test_ana: avant ',fld(i,j)
          call ana_exp_eco22(fld(i,j),d,ana) 
          write(10,*)'test_ana: reel ',d
          call ana_exp_inv_eco22(d,d2,ana) 
          write(10,*)'test_ana: apres ',d2
	  write(10,*)'********************'
	  endif
        enddo
      enddo 
    endif  
    
   close(10) 
       
   end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! subroutines for the sum-to-one constraint  !!

    subroutine ana_pi2phi2(nz,z,phi,assign)
      implicit none
      integer::nz
      real, dimension(1:nz)::z,pi_tmp
      real, dimension(1:nz-1)::phi
      integer,dimension(1:nz)::assign
      real::d,d2
      integer::k,i,j,cpt1,cpt2
      !real, parameter :: pi=3.141592653589
      real::pi2
      
      if(nz.lt.3)then
        print*,'ana_pi2phi: nz<3'
	stop
      endif
      
      pi2=acos(0.)
      
      do k=1,nz
        z(k)=abs(z(k))
      enddo
      
      !assign the variable (dp) in the cov: z(k)=0 at the end
      !assumption: dp(1) and dp(2) cannot be null    
      !sigma 2 coordinate: most likely to find null layers in the subsurface 
      if ((z(1)==0.).or.(z(2)==0.))then
        print*,'alert: fist dp null'
	stop
      endif
      cpt1=1
      cpt2=nz
      do k=1,2
        pi_tmp(k)=z(k)
        assign(k)=k
	cpt1=cpt1+1
      enddo
      /*do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).gt.1.e-3)then
	  pi_tmp(cpt1)=z(k)
	  assign(cpt1)=k
	  cpt1=cpt1+1
	else
	  pi_tmp(cpt2)=0.
	  assign(cpt2)=k
	  cpt2=cpt2-1
	endif	
      enddo*/	
      do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).gt.1.e-3)then
	  pi_tmp(cpt1)=z(k)
	  assign(cpt1)=k
	  cpt1=cpt1+1	
	endif	
      enddo
      do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).le.1.e-3)then
	  pi_tmp(cpt1)=z(k)
	  assign(cpt1)=k
	  cpt1=cpt1+1	
	endif	
      enddo
      
      if (cpt2+1-cpt1/=0)then
        print*,'alert: cpt1,cpt2 ',cpt1,cpt2
	stop
      endif
      
      phi(:)=0.
      phi(1)=(1./pi2)*acos(sqrt(pi_tmp(1)))
      do k=2,nz-2
        d=1.
	do j=1,k-1
	  d=d*sqrt((sin(pi2*phi(j)))**2)
	enddo
	if(d==0.)then
          phi(k)=0.!1.!0. !a verifier
	else
	  d2=min(1.,sqrt(pi_tmp(k))/d)
	  phi(k)=(1./pi2)*acos(d2)
	 !pb with numerical accuracy... have to truncate..!
	  !phi(k)=(1./pi2)*acos(sqrt(z(k)/d))
	endif 
	
!	d=1.
!        if(phi(k-1)==1.)then
!	  !cos(pi2*phi(k-1))=0=>pi(k-1)=0=>pi(k)=0
!	  phi(k)=1.  
!	elseif(phi(k-1)==0.)then
!	  !sin(pi2*phi(k-1))=0=>pi(k)=0
!	  phi(k)=1. 
!	else
!	  d=sqrt(pi_tmp(k-1))*abs(tan(pi2*phi(k-1)))
!	  d2=min(1.,sqrt(pi_tmp(k))/d)
!	  phi(k)=(1./pi2)*acos(d2)	  
!	endif 

      enddo
!      if((pi_tmp(nz-1).gt.0.).and.(phi(nz-2).gt.0.))then
!        !d2=min(1.,sqrt(z(nz)/z(nz-1)))
!        phi(nz-1)=(1./pi2)*atan(sqrt(pi_tmp(nz)/pi_tmp(nz-1)))  
!      else
!        phi(nz-1)=1.!pi2
!      endif
!      if((pi_tmp(nz-1)==0.).and.(pi_tmp(nz)==0.))then
!        phi(nz-1)=0.
!      endif
      
      if(phi(nz-2).gt.0.)then
        if(pi_tmp(nz-1)==0.)then
	   phi(nz-1)=0.
	else
	  phi(nz-1)=(1./pi2)*atan(sqrt(pi_tmp(nz)/pi_tmp(nz-1)))
	endif
      else
        phi(nz-1)=0.
      endif
      
      

   
   
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ana_phi2pi2(nz,z,phi,assign)
      implicit none
      integer::nz
      real, dimension(1:nz)::z,pi_tmp
      real, dimension(1:nz-1)::phi
      integer,dimension(1:nz)::assign
      real::d
      integer::k,i,j
      !real, parameter :: pi=3.141592653589
      real::pi2
      
      if(nz.lt.3)then
        print*,'ana_phi2pi: nz<3'
	stop
      endif
     
      pi2=acos(0.)
      
      z(:)=0.
      pi_tmp(1)=(cos(pi2*phi(1)))**2
      do k=2,nz-1
        d=1.
	do j=1,k-1
	  d=d*(sin(pi2*phi(j)))**2
	enddo
        pi_tmp(k)=d*(cos(pi2*phi(k)))**2
      enddo
      
      if((phi(nz-1)==1.).or.(phi(nz-1)==0.))then
        pi_tmp(nz-1)=0.
	pi_tmp(nz)=0.
      else
        pi_tmp(nz)=pi_tmp(nz-1)*(tan(pi2*phi(nz-1)))**2 
      endif
      
!      if(phi(nz-1).lt.1.)then
!        pi_tmp(nz)=pi_tmp(nz-1)*(tan(pi2*phi(nz-1)))**2 
!        !pi_tmp(nz)=d*(sin(pi2*phi(nz-1)))**2
!      else
!       pi_tmp(nz)=0.
!      endif 
          
      do k=1,nz
        !z(k)=pi_tmp(assign(k))	
        z(assign(k))=pi_tmp(k)	
      enddo
         		   
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    subroutine ana_pi2phi3(nz,z,phi,assign)
      implicit none
      integer::nz
      real, dimension(1:nz)::z,pi_tmp
      real, dimension(1:nz-1)::phi
      integer,dimension(1:nz)::assign
      real::d,d2
      integer::k,i,j,cpt1,cpt2
      !real, parameter :: pi=3.141592653589
      real::pi2
      
      if(nz.lt.3)then
        print*,'ana_pi2phi: nz<3'
	stop
      endif
      
      pi2=acos(0.)
      
      do k=1,nz
        z(k)=abs(z(k))
      enddo
      
      !assign the variable (dp) in the cov: z(k)=0 at the end
      !assumption: dp(1) and dp(2) cannot be null    
      !sigma 2 coordinate: most likely to find null layers in the subsurface 
      if ((z(1)==0.).or.(z(2)==0.))then
        print*,'alert: fist dp null'
	stop
      endif
    
      pi_tmp(nz)=z(1)
      assign(nz)=1
      pi_tmp(nz-1)=z(2)
      assign(nz-1)=2
      cpt1=1
      cpt2=nz-2
   
      do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).gt.1.e-3)then
	  pi_tmp(cpt2)=z(k)
	  assign(cpt2)=k
	  cpt2=cpt2-1	
	endif	
      enddo
      do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).le.1.e-3)then
	  pi_tmp(cpt2)=z(k)
	  assign(cpt2)=k
	  cpt2=cpt2-1	
	endif	
      enddo
      
      if (cpt2+1-cpt1/=0)then
        print*,'alert: cpt1,cpt2 ',cpt1,cpt2
	stop
      endif
      
      phi(:)=0.
      phi(1)=(1./pi2)*acos(sqrt(pi_tmp(1)))
      do k=2,nz-2
        d=1.
	do j=1,k-1
	  d=d*sqrt((sin(pi2*phi(j)))**2)
	enddo
	if(d==0.)then
         print*,'alert d=0'
	 stop
	else
	  d2=min(1.,sqrt(pi_tmp(k))/d)
	  phi(k)=(1./pi2)*acos(d2)
	 !pb with numerical accuracy... have to truncate..!
	  !phi(k)=(1./pi2)*acos(sqrt(z(k)/d))
	endif 

      enddo
      phi(nz-1)=(1./pi2)*atan(sqrt(pi_tmp(nz)/pi_tmp(nz-1)))

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ana_phi2pi3(nz,z,phi,assign)
      implicit none
      integer::nz
      real, dimension(1:nz)::z,pi_tmp
      real, dimension(1:nz-1)::phi
      integer,dimension(1:nz)::assign
      real::d
      integer::k,i,j
      !real, parameter :: pi=3.141592653589
      real::pi2
      
      if(nz.lt.3)then
        print*,'ana_phi2pi: nz<3'
	stop
      endif
     
      pi2=acos(0.)
      
      z(:)=0.
      pi_tmp(1)=(cos(pi2*phi(1)))**2
      do k=2,nz-1
        d=1.
	do j=1,k-1
	  d=d*(sin(pi2*phi(j)))**2
	enddo
        pi_tmp(k)=d*(cos(pi2*phi(k)))**2
      enddo
      if(phi(nz-1).lt.1.)then
        pi_tmp(nz)=pi_tmp(nz-1)*(tan(pi2*phi(nz-1)))**2 
        !pi_tmp(nz)=d*(sin(pi2*phi(nz-1)))**2
      else
        pi_tmp(nz)=0.
      endif 
      !pi_tmp(nz)=d*(sin(pi2*phi(nz-1)))**2
          
      do k=1,nz
        !z(k)=pi_tmp(assign(k))	
        z(assign(k))=pi_tmp(k)	
      enddo
         		   
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ana_pi2phi(nz,z,phi,assign,dpnull)
      implicit none
      integer::nz,dpnull
      real, dimension(1:nz)::z,pi_tmp
      real, dimension(1:nz-1)::phi
      integer,dimension(1:nz)::assign
      real::d,d2
      integer::k,i,j,cpt1,cpt2
      !real, parameter :: pi=3.141592653589
      real::pi2
      
      if(nz.lt.3)then
        print*,'ana_pi2phi: nz<3'
	stop
      endif
      
      pi2=acos(0.)
      
      do k=1,nz
        z(k)=abs(z(k))
      enddo
      
      !assign the variable (dp) in the cov: z(k)=0 at the end
      !assumption: dp(1) and dp(2) cannot be null    
      !sigma 2 coordinate: most likely to find null layers in the subsurface 
      if ((z(1)==0.).or.(z(2)==0.))then
        print*,'alert: fist dp null'
	stop
      endif
      cpt1=1
      do k=1,2
        pi_tmp(k)=z(k)
        assign(k)=k
	cpt1=cpt1+1
      enddo
      do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).gt.1.e-3)then
	  pi_tmp(cpt1)=z(k)
	  assign(cpt1)=k
	  cpt1=cpt1+1	
	endif	
      enddo
      cpt2=cpt1-1
      dpnull=cpt2
      do k=3,nz
      !do k=nz,3,-1
        !if(z(k).gt.0.)then
	if(z(k).le.1.e-3)then
	  pi_tmp(cpt1)=z(k)
	  assign(cpt1)=k
	  cpt1=cpt1+1	
	endif	
      enddo

      
      phi(:)=0.
      phi(1)=(1./pi2)*acos(sqrt(pi_tmp(1)))      
      select case (cpt2)
        case(2)
	   !all layers below the surface (2 first layers) are null
         ! print*,'only 2 layers not null..'	  
	
	case(3)
	  !only one layer not null below the surface
	  phi(2)=(1./pi2)*atan(sqrt(pi_tmp(3)/pi_tmp(2)))
	
	case default
          !regular case
          do k=2,cpt2-2
            d=1.
	    do j=1,k-1
	      d=d*sqrt((sin(pi2*phi(j)))**2)
	    enddo
	    if(d==0.)then
              phi(k)=0.!1.!0. !a verifier
	    else
	      d2=min(1.,sqrt(pi_tmp(k))/d)
	      phi(k)=(1./pi2)*acos(d2)
	      !pb with numerical accuracy... have to truncate..!
	    endif 	
          enddo
	  if(pi_tmp(cpt2-1)==0.)then
	    print*,'alert ana_pi2phi'
	    stop
	  endif
          phi(cpt2-1)=(1./pi2)*atan(sqrt(pi_tmp(cpt2)/pi_tmp(cpt2-1))) 
      end select

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ana_phi2pi(nz,z,phi,assign,dpnull)
      implicit none
      integer::nz,dpnull
      real, dimension(1:nz)::z,pi_tmp
      real, dimension(1:nz-1)::phi
      integer,dimension(1:nz)::assign
      real::d
      integer::k,cpt2,i,j
      !real, parameter :: pi=3.141592653589
      real::pi2
      
      if(nz.lt.3)then
        print*,'ana_phi2pi: nz<3'
	stop
      endif
     
      pi2=acos(0.)
      
!     cpt2=-1
!     do k=2,nz
!       if(assign(k).lt.assign(k-1))then
!	  cpt2=k-1
!	  exit	  
!	endif	
!      enddo
      
!      if(cpt2==-1)then
!        do k=2,nz-2
!          if((phi(k)==0.).and.(phi(k+1)==0.))then
!	    cpt2=k-1
!	    exit	  
!	  endif	
!        enddo
!      endif      
!      print*,'cpt2= ',cpt2
      
      
      z(:)=0.
      pi_tmp(:)=0.
      pi_tmp(1)=(cos(pi2*phi(1)))**2
      select case (dpnull)
        case(2)
	   !all layers below the surface (2 first layers) are null
          !print*,'only 2 layers not null..'	  
	  pi_tmp(2)=(sin(pi2*phi(1)))**2
	case(3)
	  !only one layer not null below the surface
	  pi_tmp(2)=(sin(pi2*phi(1))*cos(pi2*phi(2)))**2
	  pi_tmp(3)=(sin(pi2*phi(1))*sin(pi2*phi(2)))**2
	case default
          !regular case
        do k=2,dpnull-1
          d=1.
	  do j=1,k-1
	    d=d*(sin(pi2*phi(j)))**2
	  enddo
          pi_tmp(k)=d*(cos(pi2*phi(k)))**2
        enddo
        if((phi(dpnull-1)==1.).or.(phi(dpnull-1)==0.))then
          pi_tmp(dpnull-1)=0.
	  pi_tmp(dpnull)=0.
        else
          pi_tmp(dpnull)=pi_tmp(dpnull-1)*(tan(pi2*phi(dpnull-1)))**2 
        endif
     
      end select
     
      do k=1,dpnull
        !z(k)=pi_tmp(assign(k))	
        z(assign(k))=pi_tmp(k)	
      enddo
         		   
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine test_phi2pi(nz)
      implicit none
      integer::nz,dpnull
      real, dimension(1:nz)::z,z2  
      real, dimension(1:nz-1)::work,phi
      integer,dimension(1:nz)::assign
      integer::i,j,k

      !open(unit=10,FILE='test_phi2pi_phi.txt',FORM='formatted',STATUS='unknown')!,position='append')
      open(unit=11,FILE='test_phi2pi_z.txt',FORM='formatted',STATUS='unknown')!,position='append')
      
      do j=1,100
        
	do k=1,nz-1
	  call random_number(work(k))
	enddo
	assign(1)=1
	assign(2)=2
	do k=nz,3,-1
          assign(nz-k+3)=k
	enddo
	
        call ana_phi2pi(nz,z2,work,assign,dpnull)		
	
	call ana_pi2phi(nz,z2,phi,assign,dpnull)     
        call ana_phi2pi(nz,z,phi,assign,dpnull)
	
	
        do i=1,nz
	  if(abs(z(i)-z2(i)).gt.1e-8)then
            !write(10,*) 'randn ',work(i),phi(i)
	    write(11,*) 'z ',z2(i),z(i)-z2(i)
	  endif  
        enddo
        !write(10,*)'************************'
        !write(11,*) 'z ',z2(nz),z(nz)-z2(nz)
        write(11,*)'************************'
      
      enddo
      
      !close(10)
      close(11)
      
    end subroutine
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
  subroutine test_pi2phi_micom(nx,ny,nz,nens,depths)  
  use m_get_micom_fld
  use nfw_mod
  implicit none
  integer::nx,ny,nz,nens,nsize
  integer::m,cptz,i,j,k,ks,l,dptmp
  real, dimension(nx,ny):: fld_dp,fld,pb, depths
  real, dimension(:,:,:),allocatable::fld3d,fldtmp
  character(len=3) :: cind
  character(len=3) :: cmem
  character(len=80) :: memfile,memfilephy
  character(len=8) :: cfield
  logical::ex
  integer,dimension(nx,ny,2*nz)::assign
  real, dimension(nx,ny,2*(nz-1))::phi
  real::d
  integer::vid_phi,vid_a,ncid,dimnx,dimny,dimnz,dimnzm1,vid_depths!,vid_dp1
  integer,dimension(nz)::atmp
  real, dimension(nz-1)::ptmp
  
  
  cptz=1

  allocate(fld3d(nx,ny,nz))
  allocate(fldtmp(nx,ny,nz))

  write(cmem,'(i3.3)') nens
  memfile='forecast'//cmem 
       
   phi(:,:,:)=0.
   assign(:,:,:)=0
   !fldtmp(:,:,:)=0.
         
   do l=1,1!2
       
      cfield='pb'
      call get_micom_fld_new(trim(memfile), pb, m, cfield,&
                 l, 1, nx, ny)
      do k=1+(l-1)*nz,nz+(l-1)*nz
	cfield='dp'        
        call get_micom_fld_new(trim(memfile), fld_dp, m, cfield,&
                 k, 1, nx, ny)	 
	 
	if(k.gt.nz)then
          ks=k-nz
        else
	  ks=k
        endif
	  
	do j=1,ny
	  do i=1,nx
	    if(depths(i,j).gt.0.) then
	      fld3d(i,j,ks)=abs(fld_dp(i,j))/pb(i,j)
	    endif	
	  enddo
	enddo 
      enddo
	
	!numerical corrections: sum(fld3d,3)=1
      do j=1,ny
        do i=1,nx
	  if(depths(i,j).gt.0.)then
	    d=sum(fld3d(i,j,:)) 
	    fld3d(i,j,:)=fld3d(i,j,:)/d
	  endif
        enddo
      enddo

      do j=1,ny
	do i=1,nx
	  if(depths(i,j).gt.0.) then
	    call ana_pi2phi(nz,fld3d(i,j,:),ptmp,atmp,dptmp)
	    phi(i,j,1+(l-1)*(nz-1):nz-1+(l-1)*(nz-1))=ptmp(1:nz-1)
	    assign(i,j,1+(l-1)*nz:nz+(l-1)*nz)=atmp(1:nz)

	    call ana_phi2pi(nz,fldtmp(i,j,:),ptmp,atmp,dptmp) 
	    d=sum(fldtmp(i,j,:))  
            fldtmp(i,j,:)=(fldtmp(i,j,:)/d-fld3d(i,j,:))*pb(i,j)
	    fld3d(i,j,:)=fld3d(i,j,:)*pb(i,j)
	  endif	
	enddo
      enddo	 	
    enddo !l=1,2
	
        
    memfilephy='test_phi2pi'//cmem//'.nc' 
    call nfw_create(memfilephy, nf_write, ncid)

    call nfw_def_dim(memfilephy, ncid, 'nx', nx, dimnx)
    call nfw_def_dim(memfilephy, ncid, 'ny', ny, dimny)
    call nfw_def_dim(memfilephy, ncid, 'nz', nz, dimnz)


    call nfw_def_var(memfilephy, ncid, 'diff', nf_double, 3,(/dimnx,dimny,dimnz/), vid_phi)
    call nfw_def_var(memfilephy, ncid, 'dp', nf_double, 3,(/dimnx,dimny,dimnz/), vid_a)
    call nfw_def_var(memfilephy, ncid, 'depths', nf_double, 2,(/dimnx,dimny/), vid_depths)
	
    call nfw_enddef(memfilephy, ncid)

    call nfw_put_var_double(memfilephy, ncid, vid_phi, fldtmp)
    call nfw_put_var_double(memfilephy, ncid, vid_a, fld3d)
    call nfw_put_var_double(memfilephy, ncid, vid_depths, depths)
	
    call nfw_close(memfilephy, ncid)
     

    deallocate(fld3d,fldtmp)


  end subroutine
   
end module m_ana_exp_enkf 
