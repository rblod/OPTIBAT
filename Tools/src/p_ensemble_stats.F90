program ensemble_stats
   use m_ensemble_stats
   use nfw_mod    
   implicit none

   integer imem, nrens, nrfields, af     
   character(len=*), parameter :: infile='stats_fields.in'
   character(len=3), dimension(:), allocatable:: fieldsnames
   character(len=8) :: cfld, ctmp,cmem  
   integer::ncid, field_id,dim1_id,dim2_id,ncid2
   character(len=80) :: fname1,fname2
   integer ::idm,jdm,k,kens
   real, allocatable, dimension(:,:,:)     :: fld1
   real, allocatable, dimension(:,:)     :: fld2,fld3
   integer,dimension(2)::nc,ns
   
   if (command_argument_count() == 2) then
      call getarg(1,ctmp)
      read(ctmp,*) nrens
      call getarg(2,ctmp)
      read(ctmp,*) af
   else
      print *,'ensemble_stats'
      print *
      print *,'usage: '
      print *,'  ensemble_stats nens af'
      print *,'   "nens" is the ensemble size'
      print *,'   "af" = 0 if forecast, 1 if analysis'
      call exit(1)
   endif

   
   ! fields that will be analyzed
   nrfields=get_nrfields_stats(infile)
   allocate(fieldsnames(nrfields))
   call get_fields_stats(infile,nrfields,fieldsnames)
  
   ! netcdf file
   ! get the dimensions
   if(af==0)then
      fname2='forecast001.nc'
   else
      fname2='analysis001.nc'
   end if
   call nfw_open(trim(fname2), nf_nowrite, ncid2)
   call nfw_inq_dimid(trim(fname2), ncid2, 'xi_rho', dim1_id)
   call nfw_inq_dimlen(trim(fname2), ncid2, dim1_id, idm)
   call nfw_inq_dimid(trim(fname2), ncid2, 'eta_rho', dim2_id)
   call nfw_inq_dimlen(trim(fname2), ncid2, dim2_id, jdm)
   call nfw_close(trim(fname2), ncid2)
   
   allocate(fld1(idm,jdm,nrens)) 
   allocate(fld2(idm,jdm))
   allocate(fld3(idm,jdm))
   
   ns(1)=1
   ns(2)=1
   nc(1)=idm
   nc(2)=jdm
   
   ! create the file
   fname1='Ensemble_stats.nc'
   call nfw_create(trim(fname1),nf_write, ncid)
   call nfw_def_dim(trim(fname1),ncid,'xi_rho',idm,dim1_id)
   call nfw_def_dim(trim(fname1),ncid,'eta_rho',jdm,dim2_id)
   call nfw_enddef(trim(fname1),ncid)
   
   do k=1,nrfields
   	cfld=trim(fieldsnames(k))
	print*,'field ',trim(cfld)
	
	fld1(:,:,:)=0
	fld2(:,:)=0  ! ensemble mean
        
	do kens=1,nrens
	   write(cmem,'(i3.3)')kens
	   if(af==0)then
	      fname2='forecast'//trim(cmem)//'.nc'
	   else
	      fname2='analysis'//trim(cmem)//'.nc'
	   end if
	   
	   call nfw_open(trim(fname2), nf_nowrite, ncid2)
	   call nfw_inq_varid(trim(fname2), ncid2,cfld,field_id)
           call nfw_get_vara_double(trim(fname2), ncid2, field_id, ns, nc, fld3)
	   call nfw_close(trim(fname2), ncid2)
	   
	   fld2(:,:)=fld2(:,:)+fld3(:,:)
	   fld1(:,:,kens)=fld3(:,:)
	enddo
	fld2(:,:)=fld2(:,:)/real(nrens)
	
	fld3(:,:)=0 ! ensemble spread
	! second loop, not prone to rounding errors
	do kens=1,nrens
	   fld3(:,:)=fld3(:,:)+(fld1(:,:,kens)-fld2(:,:))**2
	enddo
	fld3(:,:)=sqrt(fld3(:,:)/(real(nrens)-1.))
	
	! storing the variables in the netcdf files
	ctmp=trim(cfld)//'_mean'
	call nfw_redef(trim(fname1),ncid)
	call nfw_def_var(trim(fname1),ncid,trim(ctmp),nf_float,2,(/dim1_id,dim2_id/),field_id)
	call nfw_enddef(trim(fname1),ncid)
	call nfw_put_var_double(trim(fname1),ncid,field_id,fld2(1:idm,1:jdm))
	
	ctmp=trim(cfld)//'_std'
	call nfw_redef(trim(fname1),ncid)
	call nfw_def_var(trim(fname1),ncid,trim(ctmp),nf_float,2,(/dim1_id,dim2_id/),field_id)
	call nfw_enddef(trim(fname1),ncid)
	call nfw_put_var_double(trim(fname1),ncid,field_id,fld3(1:idm,1:jdm))
		
   enddo
   
   call nfw_close(trim(fname1),ncid)
   
end program
