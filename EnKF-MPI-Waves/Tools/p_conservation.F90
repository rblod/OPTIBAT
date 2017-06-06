program conservation
!Ehouarn!
!prog converting eco restart files .ab to netcd format!
#if defined (QMPI)
      use qmpi
#else
      use qmpi_fake
#endif   
 !  use mod_raw_io
   use m_get_micom_grid
   use m_get_micom_dim
   use mod_eosfun
   use netcdf
   use nfw_mod
   implicit none

   integer*4, external :: iargc

   integer imem                  ! ensemble member
   character(len=80) :: memfile,memfilenc,memfile2
   logical          :: ex
   character(len=8) :: cfld,ctmp,cfld2
   character(len=3) :: cmem
   integer          :: idm,jdm,kdm,nrens,ana
   real,dimension(:,:),allocatable::depths,modlon,modlat
   real,dimension(:,:,:),allocatable::ma,mf,dpa,dpf,diff,mtf,mta,difft,fldt
   real,dimension(:,:,:,:),allocatable:: dp,saln,temp
   !real, allocatable, dimension(:,:,:) ::  nit,pho,sil,det,oxy,sis,fld_dp,fla,dia,chl,detp,mes,mic
   integer :: i,k,fnd,j,fnd2,kens
   integer :: dimx,var_id,var1d(1),ncid,ierr,dimnx,dimny,dimnz,cpt
   real::rdummy   
   real::onem
   parameter(onem=98060.)
   real::sig2,mindx,meandx
   real::mmf,mma,mdpf,mdpa,mmtf,mmta
   
   real, parameter ::c1= 9.77093E+00
   real, parameter ::   c2= -2.26493E-02
   real, parameter ::   c3= 7.89879E-01
   real, parameter ::   c4= -6.43205E-03
   real, parameter ::  c5= -2.62983E-03
   real, parameter ::  c6= 2.75835E-05
   real, parameter ::  c7= 3.48658E-05
   integer::z_id,vDP_ID,vPBOT_ID
   integer,dimension(:),allocatable::nc,ns
   
   if (iargc()==2) then
      call getarg(1,ctmp)
      read(ctmp,*) nrens
      call getarg(2,cfld2)
   else
      print *,'usage: ensemble_onservation nrens tracer'
      call exit(1)
   endif
   
   ! Get dimensions from blkdat
   call get_micom_dim(idm,jdm)

   allocate(mf  (idm,jdm,nrens+1))
   allocate(ma  (idm,jdm,nrens+1))  
   allocate(diff  (idm,jdm,nrens+1))  
   allocate(mtf  (idm,jdm,nrens+1))
   allocate(mta  (idm,jdm,nrens+1))  
   allocate(difft  (idm,jdm,nrens+1))
   allocate(dpf  (idm,jdm,nrens+1))
   allocate(dpa  (idm,jdm,nrens+1)) 
  
  

   allocate(modlon(idm,jdm))
   allocate(modlat(idm,jdm))
   allocate(depths(idm,jdm)) 
   
   call get_micom_grid(modlon, modlat, depths, mindx, meandx, idm,jdm)
   deallocate(modlat,modlon)

   call eosini
   
   mf(:,:,:)=0.
   ma(:,:,:)=0. 
   dpf(:,:,:)=0.
   dpa(:,:,:)=0.
   mtf(:,:,:)=0.
   mta(:,:,:)=0. 
   print*,'connard'
   do kens=1,nrens  
      write(cmem,'(i3.3)') kens
      
      memfile='forecast'//cmem

      call nfw_open(trim(memfile)//'.nc', nf_write, ncid)
      call nfw_inq_dimid(trim(memfile)//'.nc', ncid, 'kk', z_ID)
      call nfw_inq_dimlen(trim(memfile)//'.nc', ncid, z_ID, kdm)
  
      allocate(saln(idm,jdm,2*kdm,1 ))
      allocate(dp (idm,jdm,2*kdm,1))
      allocate(temp(idm,jdm,2*kdm,1 ))
   
      !Reading dp 
      print *,'Reading dp and pb'
      allocate(ns(4))
      allocate(nc(4))
      ns(1)=1
      ns(2)=1
      ns(3)=1
      ns(4)=1
      nc(1)=idm
      nc(2)=jdm
      nc(3)=2*kdm
      nc(4)=1
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'dp',vDP_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vDP_ID, ns, nc, dp)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'temp',vPBOT_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vPBOT_ID, ns, nc, temp)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'saln',vPBOT_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vPBOT_ID, ns, nc, saln)
      call nfw_close(trim(memfile)//'.nc',ncid) 
      deallocate(ns,nc)
     
      memfile='forecastECO'//cmem
      allocate(fldt(idm,jdm,kdm))
      call nfw_open(trim(memfile)//'.nc', nf_write, ncid)
      allocate(ns(3))
      allocate(nc(3))
      ns(1)=1
      ns(2)=1
      ns(3)=1
      nc(1)=idm
      nc(2)=jdm
      nc(3)=kdm
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld2),vDP_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vDP_ID, ns, nc, fldt)
      call nfw_close(trim(memfile)//'.nc',ncid) 
      deallocate(ns,nc)
      
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
	if(depths(i,j).gt.0.)then
	  !sig2=c1+c3*saln(i,j,k,1)+temp(i,j,k,1)&
	  !       *(c2+c5*saln(i,j,k,1)+temp(i,j,k,1)*(c4+c7*saln(i,j,k,1)+c6*temp(i,j,k,1)))	 
	  sig2=sig(temp(i,j,k+kdm,1),saln(i,j,k+kdm,1))   
	  mf(i,j,kens)=mf(i,j,kens)+dp(i,j,k+kdm,1)/(sig2*onem)
	  dpf(i,j,kens)=dpf(i,j,kens)+dp(i,j,k+kdm,1)/onem
	    
	  mtf(i,j,kens)=mtf(i,j,kens)+fldt(i,j,k)*dp(i,j,k+kdm,1)/onem
	  if((i==9).and.(j==99))then
	    print*,'forecast',k+kdm,fldt(i,j,k)*dp(i,j,k+kdm,1)/onem
	    print*,'forecast',fldt(i,j,k),dp(i,j,k+kdm,1)/onem
	    print*,'forecast',mtf(i,j,kens),dpf(i,j,kens)
	    print*,'forecast ****************'
	  endif
	endif
      enddo
      enddo
      enddo
      
      !!!!!!!!!!!!!!!!!1
      
      memfile='analysis'//cmem
      call nfw_open(trim(memfile)//'.nc', nf_write, ncid)
      call nfw_inq_dimid(trim(memfile)//'.nc', ncid, 'kk', z_ID)
      call nfw_inq_dimlen(trim(memfile)//'.nc', ncid, z_ID, kdm)
   
      !Reading dp 
      print *,'Reading dp and pb'
      allocate(ns(4))
      allocate(nc(4))
      ns(1)=1
      ns(2)=1
      ns(3)=1
      ns(4)=1
      nc(1)=idm
      nc(2)=jdm
      nc(3)=2*kdm
      nc(4)=1
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'dp',vDP_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vDP_ID, ns, nc, dp)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'temp',vPBOT_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vPBOT_ID, ns, nc, temp)
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'saln',vPBOT_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vPBOT_ID, ns, nc, saln)
      call nfw_close(trim(memfile)//'.nc',ncid) 
      deallocate(ns,nc)
     
      memfile='analysisECO'//cmem
      allocate(fldt(idm,jdm,kdm))
      call nfw_open(trim(memfile)//'.nc', nf_write, ncid)
      allocate(ns(3))
      allocate(nc(3))
      ns(1)=1
      ns(2)=1
      ns(3)=1
      nc(1)=idm
      nc(2)=jdm
      nc(3)=kdm
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld2),vDP_ID)
      call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vDP_ID, ns, nc, fldt)
      call nfw_close(trim(memfile)//'.nc',ncid) 
      
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
	if(depths(i,j).gt.0.)then
	  !sig2=c1+c3*saln(i,j,k,1)+temp(i,j,k,1)&
	  !       *(c2+c5*saln(i,j,k,1)+temp(i,j,k,1)*(c4+c7*saln(i,j,k,1)+c6*temp(i,j,k,1)))	    
	  sig2=sig(temp(i,j,k+kdm,1),saln(i,j,k+kdm,1)) 
	  ma(i,j,kens)=ma(i,j,kens)+dp(i,j,k+kdm,1)/(sig2*onem)
	  dpa(i,j,kens)=dpa(i,j,kens)+dp(i,j,k+kdm,1)/onem
	    
	  mta(i,j,kens)=mta(i,j,kens)+fldt(i,j,k)*dp(i,j,k+kdm,1)/onem
	  if((i==9).and.(j==99))then
	    print*,'analysis',k+kdm,fldt(i,j,k)*dp(i,j,k+kdm,1)/onem
	    print*,'analysis',fldt(i,j,k),dp(i,j,k+kdm,1)/onem
	    print*,'analysis',mta(i,j,kens),dpa(i,j,kens)
	    print*,'analysis ****************'
	  endif
	endif
      enddo
      enddo
      enddo

    
    enddo             
   
   mmf=0.
   mma=0.
   mdpf=0.
   mdpa=0.
   cpt=0
   diff(:,:,:)=0.
   mmtf=0.
   mmta=0.
   difft(:,:,:)=0.
   do kens=1,nrens
      do j=1,jdm
      do i=1,idm
        if(depths(i,j).gt.0.)then
	  mf(i,j,nrens+1)=mf(i,j,nrens+1)+mf(i,j,kens)/real(nrens)
	  dpf(i,j,nrens+1)=dpf(i,j,nrens+1)+dpf(i,j,kens)/real(nrens)
	  mmf=mmf+mf(i,j,kens)
	  mdpf=mdpf+dpf(i,j,kens)
	  
	  ma(i,j,nrens+1)=ma(i,j,nrens+1)+ma(i,j,kens)/real(nrens)
	  dpa(i,j,nrens+1)=dpa(i,j,nrens+1)+dpa(i,j,kens)/real(nrens)
	  mma=mma+ma(i,j,kens)
	  mdpa=mdpa+dpa(i,j,kens)
	  
	  !diff(i,j,kens)=(ma(i,j,kens)-mf(i,j,kens))*100./mf(i,j,kens)
	  !diff(i,j,nrens+1)=diff(i,j,nrens+1)+&
	  !                (ma(i,j,kens)-mf(i,j,kens))/real(nrens)
	  
	  diff(i,j,kens)=(dpa(i,j,kens)-dpf(i,j,kens))*100./dpf(i,j,kens)
	  diff(i,j,nrens+1)=diff(i,j,nrens+1)+&
	                  (dpa(i,j,kens)-dpf(i,j,kens))/real(nrens)
	   
	  !mtf(i,j,kens)=mtf(i,j,kens)/dpf(i,j,kens)
	  mtf(i,j,nrens+1)=mtf(i,j,nrens+1)+mtf(i,j,kens)/real(nrens)
	  mmtf=mmtf+mtf(i,j,kens)
	  
	  !mta(i,j,kens)=mta(i,j,kens)/dpa(i,j,kens)
	  mta(i,j,nrens+1)=mta(i,j,nrens+1)+mta(i,j,kens)/real(nrens)
	  mmta=mmta+mta(i,j,kens)
	  
	  difft(i,j,kens)=(mta(i,j,kens)-mtf(i,j,kens))*100./mtf(i,j,kens)
	  difft(i,j,nrens+1)=difft(i,j,nrens+1)+&
	                  (mta(i,j,kens)-mtf(i,j,kens))/real(nrens)
	  		  
	  cpt=cpt+1
	endif
      enddo
      enddo
   enddo 
   
   print*,'mmf,mdpf,mma,mdpa'
   print*,mmf/real(cpt),mdpf/real(cpt),mma/real(cpt),mdpa/real(cpt)
   
   
   print*,'mmtf,mmta'
   print*,mmtf/real(cpt),mmta/real(cpt)
   print*,'difft ', (mmtf-mmta)/mmtf
   
   
   memfilenc='conservation.nc'
   call nfw_create(trim(memfilenc), nf_write, ncid)

   call nfw_def_dim(trim(memfilenc),ncid,'nx',idm,dimnx)
   call nfw_def_dim(trim(memfilenc),ncid,'ny',jdm,dimny)
   call nfw_def_dim(trim(memfilenc),ncid,'nens',nrens+1,dimnz)
   call nfw_enddef(trim(memfilenc), ncid)

   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'mass_forecast', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,mf(1:idm,1:jdm,1:nrens+1))   
     
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'vol_forecast', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,dpf(1:idm,1:jdm,1:nrens+1))  
     
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'mass_analysis', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,ma(1:idm,1:jdm,1:nrens+1))   
   
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'vol_analysis', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,dpa(1:idm,1:jdm,1:nrens+1)) 
   
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'mass_diff', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,diff(1:idm,1:jdm,1:nrens+1)) 
  
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'tracer_forecast', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,mtf(1:idm,1:jdm,1:nrens+1))
   
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'tracer_analysis', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,mta(1:idm,1:jdm,1:nrens+1))
   
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'tracer_diff', nf_float,3,(/dimnx,dimny,dimnz/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,difft(1:idm,1:jdm,1:nrens+1)) 
   
   call nfw_redef(trim(memfilenc),ncid)
   call nfw_def_var(trim(memfilenc), ncid, 'depths', nf_float,2,(/dimnx,dimny/), var_id)
   call nfw_enddef(trim(memfilenc), ncid)
   call nfw_put_var_double(trim(memfilenc), ncid,var_id,depths(1:idm,1:jdm)) 
   
   call nfw_close(trim(memfilenc), ncid)
   
   deallocate(depths)
   deallocate(mf,ma,dpf,dpa)
   
#if defined(AIX)
   call exit_(0)
#else
   call exit(0)
#endif

end program
