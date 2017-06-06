program anamorphosis_state
!Ehouarn!
!prog that generates samples file!
   use qmpi_fake
   use mod_raw_io
   use m_parse_blkdat
   use m_get_micom_nrens
   use m_get_micom_grid
   use m_get_micom_dim
   use m_get_micom_fld
   use m_ana_exp_enkf
   use m_sort2
   use netcdf
   use nfw_mod
   implicit none

   integer*4, external :: iargc

   integer:: imem                ! ensemble member
   character(len=80) :: ecotemplate,fname
   character(len=80) :: memfile,memfilephy
   logical          :: ex
   character(len=8) :: cfld,ctmp,method
   character(len=3) :: cstart,cind,ck
   integer          :: idm,jdm,kdm,klay
   integer :: i,k,fnd,j,m,cptz,nchl,neco,nrens
   integer :: dimx,var_id,var1d(1),ncid,ierr,dimnx,dimny,dimnz
   real::rdummy
   !type(anamorphosis), dimension(:),allocatable::ana   
   type(anamorphosis)::ana
   real,dimension(:,:),allocatable::modlat, modlon,depths,fld
   integer, dimension(:), allocatable::itemp
   real,dimension(:),allocatable::zt
   real:: mindx,meandx,latmax,radius
   logical::toosmall,zero_phi
   
   if (iargc()==1) then
      call getarg(1,ctmp)
      read(ctmp,*) kdm  
   else
      print *,'usage: anamorphosis kdim'
      call exit(1)
   endif

   ! Get dimensions from blkdat
   call get_micom_dim(idm,jdm)
   !call parse_blkdat('kdm   ','integer',rdummy,kdm)
   
   !get model grid!
   allocate(modlon(idm,jdm))
   allocate(modlat(idm,jdm))
   allocate(depths(idm,jdm)) 
   allocate(fld(idm,jdm))
   call get_micom_grid(modlon, modlat, depths, mindx, meandx, idm, jdm)
   deallocate(modlon,modlat)
   
   !reading of the data which will be anamorphosed: analysisfields_ana.in!
   call get_analysisfields_ana()
   
   !allocation of the anamorphosis structure!
   !allocate(ana(numfields_ana))
   
   !size of the ensmeble!
   nrens = get_micom_nrens(idm, jdm)
   print*,'nrens= ', nrens
		       
   open(unit=50,FILE='ana_size.txt',FORM='formatted',STATUS='unknown')
   
   !reading of the data and storage in ana!   
   do k=1,numfields_ana
      !id of the data!
      ana%id=trim(fieldnames_ana(k))
      print *,'----------   ',ana%id,'   ----------'
     
      !reading of the smoothing param of the anamorphosis!
      call ana_smooth_param(ana)
     
      !size of the sample!
      print *,'----------   nb sample   ----------'
      call ana_size_sample(ana,idm,jdm,kdm,nrens,depths)
      allocate(zt(ana%n))
      
      if(mod(ana%n-2,ana%nsamp)==0)then
        ana%samp=(ana%n-2)/ana%nsamp+1 !eho21/12
      else
        ana%samp=(ana%n-2)/ana%nsamp+2
      endif
      allocate(ana%z(ana%samp))
      
      !reading the data!   
      call ana_data_z(ana,zt,idm,jdm,kdm,nrens,ana%n,depths)
      !do i=1,ana%n
      !print*,'connard',zt(i)
      !enddo
      
      !sorting the array Z!
      print *,'sort Z'
      allocate(itemp(ana%n))
      call quick_sort(zt,itemp)
      deallocate(itemp)

      !sampling for the computation of invG!
      !call ana_samp_z(ana(k))
      
      !computation of yi=invG(i/n)!
      print *,'gauss Y'
      call ana_gauss_y(ana)
     
      !smooth anamorphosis!
      print *,'phi'
      call ana_cpt_phi_eco3(ana,zt,ana%n)
      
      !print size of the y,phy array!
      write(50,*)ana%samp
      
      !cleaning
      deallocate(zt,ana%y,ana%z,ana%phi)
      
   enddo!numfields_ana!
   deallocate(depths)
   
   close(50)
   
#if defined(AIX)
   call exit_(0)
#else
   call exit(0)
#endif

end program
