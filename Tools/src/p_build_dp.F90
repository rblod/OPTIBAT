! File:          p_build_dp.F90
!
program build_dp
use netcdf
use nfw_mod
use m_ana_exp_enkf
   implicit none

!   integer*4, external :: iargc
   real, parameter :: onem=9806.

   integer imem                  ! ensemble member
   character(len=80) :: oldfile,newfile, forfile, testfile
   logical          :: ex
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc,cmem
   integer          :: tlevel, vlevel, nproc
   integer          :: nx,ny,nz
   real, allocatable:: fld(:,:)
   real, allocatable, dimension(:,:)     :: depths,dp1,dp2
   real, allocatable, dimension(:,:,:,:)   :: pb,dp,kfpla,pbu,pbv,tmp
   real, allocatable, dimension(:,:,:)   :: phi,dp_ana
   integer, allocatable, dimension(:,:,:)   :: assign,dpnull
   real, allocatable,dimension(:)::dp_tmp
   integer, dimension(:,:),allocatable::ins,jns,inw,jnw
   
   real, parameter  :: epsil=1.e-6!11

   integer,parameter :: numfields=2
   integer :: ios,ios2
   integer :: i,j,k
   real :: dpsum,d

   integer :: ncid, x_ID, y_ID, z_ID, fld_id
   integer ::vid1,vid2,vid_depths,dimnz,dimnx,dimny
   integer :: ncid2,ins_ID,jns_ID,inw_ID,jnw_ID,vPBU_ID,vPBV_ID,vTMP_ID
   !real, allocatable :: press(:)
   integer, allocatable :: ns(:), nc(:)
   
   
    print*,'build_dp: recontructing dp from analysis_phi???.nc'

   if (iargc()==1 ) then
      call getarg(1,ctmp)
      read(ctmp,*) imem
      write(cmem,'(i3.3)') imem
   else
      print *,'"build_dp" -- Compute dp from phi (hyperspherical formulation)'
      print *
      print *,'usage: '
      print *,'   build_dp ensemble_member'
      print *,'   "ensemble_member" is the ensemble member'
      call exit(1)
   endif
  
   oldfile='analysis_phi'//cmem//'.nc'
   newfile='analysis'//cmem//'.nc'
   forfile='forecast_phi'//cmem//'.nc'
   
   
   !!!!!!!!!! forecast : assign !!!!!!
   call nfw_open(forfile, nf_nowrite, ncid)
   call nfw_inq_dimid(forfile, ncid, 'nx', x_ID)
   call nfw_inq_dimid(forfile, ncid, 'ny', y_ID)
   call nfw_inq_dimid(forfile, ncid, 'nz', z_ID)

   call nfw_inq_dimlen(forfile, ncid, x_ID, nx)
   call nfw_inq_dimlen(forfile, ncid, y_ID, ny)
   call nfw_inq_dimlen(forfile, ncid, z_ID, nz)
   
   nz=nz/2
   allocate(assign(nx,ny,2*nz))
   allocate(depths(nx,ny))
   allocate(dp_tmp(nz))
   allocate(dpnull(nx,ny,2))
      
   cfld='assign'
   call nfw_inq_varid(forfile, ncid,trim(cfld),fld_id)
   call nfw_get_vara_int(forfile, ncid,fld_id,(/1,1,1/),(/nx,ny,2*nz/),assign)
   
   cfld='dpnull'
   call nfw_inq_varid(forfile, ncid,trim(cfld),fld_id)
   call nfw_get_vara_int(forfile, ncid,fld_id,(/1,1,1/),(/nx,ny,2/),dpnull)
   
   cfld='depths'
   call nfw_inq_varid(forfile, ncid,trim(cfld),fld_id)
   call nfw_get_vara_double(forfile, ncid,fld_id,(/1,1/),(/nx,ny/),depths)
  
   call nfw_close(forfile,ncid)
   
   !!!!!!!!! analysis : phi  !!!!!!!!!!!!!
   
   allocate(phi(nx,ny,2*(nz-1)))
   call nfw_open(oldfile, nf_nowrite, ncid)
   
   cfld='dpphi'
   call nfw_inq_varid(oldfile, ncid,trim(cfld),fld_id)
   call nfw_get_vara_double(oldfile, ncid,fld_id,(/1,1,1/),(/nx,ny,2*(nz-1)/),phi)
   call nfw_close(oldfile,ncid)
   
   !!!!!!!!!!! computation of dp !!!!!!!!!!
   allocate(dp(nx,ny,2*nz,1))
   dp(:,:,:,:)=-1.
   
   allocate(pb(nx,ny,2,1)) 
   allocate(kfpla(nx,ny,2,1)) 
   pb(:,:,:,:)=0.
   kfpla(:,:,:,:)=0.
   
   call nfw_open(newfile, nf_write, ncid)
   cfld='pb'
   call nfw_inq_varid(newfile, ncid,trim(cfld),fld_id)
   call nfw_get_vara_double(newfile, ncid,fld_id,(/1,1,1,1/),(/nx,ny,2,1/),pb)
   
   do j=1,ny
      do i=1,nx
        if(depths(i,j).gt.0.)then
	   call ana_phi2pi(nz,dp_tmp,phi(i,j,1:nz-1),assign(i,j,1:nz),dpnull(i,j,1))
	   d=sum(dp_tmp)
	   dp(i,j,1:nz,1)=pb(i,j,1,1)*dp_tmp(1:nz)/d
	   !call ana_phi2pi(nz,dp_tmp,phi(i,j,nz:2*nz-2),assign(i,j,nz+1:2*nz))
	   !d=sum(dp_tmp)
	   dp(i,j,nz+1:2*nz,1)=dp(i,j,1:nz,1)
	   pb(i,j,2,1)=pb(i,j,1,1)
	endif
      enddo    
    enddo    

    cfld='dp'
    call nfw_inq_varid(newfile, ncid,trim(cfld),fld_id)
    call nfw_put_vara_double(newfile, ncid, fld_id,(/1,1,1,1/),(/nx,ny,2*nz,1/),dp)
    
    deallocate(assign,dp_tmp,phi,dpnull)
    
    !!!!!!!!!!!!!!!!!!! computation of kfpla!!!!!!!!!!!!!!!!
    
    do j=1,ny
      do i=1,nx
          if(depths(i,j).gt.0.)then
          k=3
          do while (dp(i,j,k,1).lt.epsil)
            k=k+1
            if (k.gt.nz) exit
          enddo
          kfpla(i,j,1,1)=real(k)
	  !k=35
          !do while (dp(i,j,k,1).lt.epsil)
          !  k=k+1
          !  if (k.gt.2*nz) exit
          !enddo
	  kfpla(i,j,2,1)=kfpla(i,j,1,1)
	  endif
      enddo
    enddo
    
    cfld='kfpla'
    call nfw_inq_varid(newfile, ncid,trim(cfld),fld_id)
    call nfw_put_vara_double(newfile, ncid, fld_id,(/1,1,1,1/),(/nx,ny,2,1/),kfpla)
    deallocate(kfpla)
    
    !!!!!!!!!!!!!!!! computation of dpu, dpv !!!!!!!!!!!!!!!!!!!!!
    
   call nfw_inq_varid(newfile, ncid,'pbu',vPBU_ID)
   call nfw_inq_varid(newfile, ncid,'pbv',vPBV_ID)

   allocate(pbu (nx,ny,2,1    ))
   allocate(pbv (nx,ny,2,1    ))
   call nfw_get_vara_double(trim(newfile)//'.nc', ncid, vPBU_ID,(/1,1,1,1/),(/nx,ny,2,1/),pbu)
   call nfw_get_vara_double(trim(newfile)//'.nc', ncid, vPBV_ID,(/1,1,1,1/),(/nx,ny,2,1/),pbv)
   allocate(ins(nx,ny))
   allocate(jns(nx,ny))
   allocate(inw(nx,ny))
   allocate(jnw(nx,ny))

    call nfw_open('grid.nc', nf_write, ncid2)
    call nfw_inq_varid('grid.nc', ncid2,'jns',jns_ID)
    call nfw_inq_varid('grid.nc', ncid2,'ins',ins_ID)
    call nfw_inq_varid('grid.nc', ncid2,'jnw',jnw_ID)
    call nfw_inq_varid('grid.nc', ncid2,'inw',inw_ID)

    call nfw_get_vara_int('grid.nc', ncid2, jns_ID, (/1,1/),(/nx,ny/),jns)
    call nfw_get_vara_int('grid.nc', ncid2, ins_ID,(/1,1/),(/nx,ny/),ins)
    call nfw_get_vara_int('grid.nc', ncid2, jnw_ID,(/1,1/),(/nx,ny/),jnw)
    call nfw_get_vara_int('grid.nc', ncid2, inw_ID,(/1,1/),(/nx,ny/),inw)
    
    call nfw_close('grid.nc',ncid2)
   
   !Here we do not have acess to ifu and ilu, so we need to make some tricks to
   !make the ponts close to the land equal to 0
    do j=1,ny
      do i=1,nx
         if ( pbu(i,j,1,1) .ne. 0) then
            pbu(i,j,1,1)=min(pb(i,j,1,1),pb(inw(i,j),jnw(i,j),1,1))
            pbu(i,j,2,1)=pbu(i,j,1,1)
	 endif
         if ( pbv(i,j,1,1) .ne. 0) then
            pbv(i,j,1,1)=min(pb(i,j,1,1),pb(ins(i,j),jns(i,j),1,1))          
            pbv(i,j,2,1)=pbv(i,j,1,1)
	 endif
	 
	 !if ( pbu(i,j,2,1) .ne. 0) then
         !   pbu(i,j,2,1)=min(pb(i,j,2,1),pb(inw(i,j),jnw(i,j),2,1))
         !endif
         !if ( pbv(i,j,2,1) .ne. 0) then
          !  pbv(i,j,2,1)=min(pb(i,j,2,1),pb(ins(i,j),jns(i,j),2,1))
         !endif
      enddo
   enddo
    
   call nfw_put_vara_double(newfile, ncid,vPBU_ID,(/1,1,1,1/),(/nx,ny,2,1/),pbu)
   call nfw_put_vara_double(newfile, ncid,vPBV_ID,(/1,1,1,1/),(/nx,ny,2,1/),pbv)    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 
  
   allocate(ns(4))
   allocate(nc(4))
   allocate(tmp   (nx,ny,2*nz,1))
  
   ns(1)=1
   ns(2)=1
   ns(3)=1
   ns(4)=1
   nc(1)=nx
   nc(2)=ny
   nc(3)=2*nz
   nc(4)=1
   !temp
   call nfw_inq_varid(newfile, ncid,'temp',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !saln
   call nfw_inq_varid(newfile, ncid,'saln',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !u
   call nfw_inq_varid(newfile, ncid,'u',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !v
   call nfw_inq_varid(newfile, ncid,'v',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !sigma
   call nfw_inq_varid(newfile, ncid,'sigma',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !pgfx
   call nfw_inq_varid(newfile, ncid,'pgfx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !pgfy
   call nfw_inq_varid(newfile, ncid,'pgfy',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !uflx
   call nfw_inq_varid(newfile, ncid,'uflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !vflx
   call nfw_inq_varid(newfile, ncid,'vflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !utflx
   call nfw_inq_varid(newfile, ncid,'utflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !vtflx
   call nfw_inq_varid(newfile, ncid,'vtflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !usflx
   call nfw_inq_varid(newfile, ncid,'usflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !vsflx
   call nfw_inq_varid(newfile, ncid,'vsflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,nz+1:2*nz,1)=tmp(:,:,1:nz,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !!!!! 2D variables
   nc(3)=2
   !ub
   call nfw_inq_varid(newfile, ncid,'ub',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !vb
   call nfw_inq_varid(newfile, ncid,'vb',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !ubflx
   call nfw_inq_varid(newfile, ncid,'ubflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !vb
   call nfw_inq_varid(newfile, ncid,'vbflx',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !pvtrop
   call nfw_inq_varid(newfile, ncid,'pvtrop',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !pgfxm
   call nfw_inq_varid(newfile, ncid,'pgfxm',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !pgfym
   call nfw_inq_varid(newfile, ncid,'pgfym',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !xixm
   call nfw_inq_varid(newfile, ncid,'xixm',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !xiym
   call nfw_inq_varid(newfile, ncid,'xiym',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !xixp
   call nfw_inq_varid(newfile, ncid,'xixp',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   !xiyp
   call nfw_inq_varid(newfile, ncid,'xiyp',vTMP_ID)
   call nfw_get_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(newfile, ncid, vTMP_ID, ns, nc, tmp)
  
   call nfw_close(newfile,ncid)
  
   deallocate(pbu,pbv,ins,jns,inw,jnw)
   deallocate(ns,nc,tmp)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !check the sum of dp
    /*allocate(dp1(nx,ny))
    allocate(dp2(nx,ny))
    dp1(:,:)=-1.
    dp2(:,:)=-1.
    
    do j=1,ny
      do i=1,nx
        if(depths(i,j).gt.0.)then	   
	   dp1(i,j)=sum(dp(i,j,1:nz,1))/pb(i,j,1,1)
	   dp2(i,j)=sum(dp(i,j,nz+1:2*nz,1))/pb(i,j,2,1)   
	endif
      enddo    
    enddo 
    
       
    testfile='test_analysis_'//cmem//'.nc' 
    call nfw_create(testfile, nf_write, ncid)

    call nfw_def_dim(testfile, ncid, 'nx', nx, dimnx)
    call nfw_def_dim(testfile, ncid, 'ny', ny, dimny)
    call nfw_def_dim(testfile, ncid, 'nz', nz, dimnz)

    call nfw_def_var(testfile, ncid, 'depths', nf_double, 2,(/dimnx,dimny/), vid_depths)    
    call nfw_def_var(testfile, ncid, 'dp1', nf_double, 2,(/dimnx,dimny/), vid1)    
    call nfw_def_var(testfile, ncid, 'dp2', nf_double, 2,(/dimnx,dimny/), vid2)
    
    call nfw_enddef(testfile, ncid)

    call nfw_put_var_double(testfile, ncid, vid1, dp1)
    call nfw_put_var_double(testfile, ncid, vid2, dp2)
    call nfw_put_var_double(testfile, ncid, vid_depths, depths)
	
    call nfw_close(testfile, ncid)
    deallocate(dp1,dp2)*/
  
   !!!!!!!!!!
    deallocate(pb,dp,depths)
    
    
     
end program
