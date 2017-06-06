! File:          p_fixenkf.F90
!
! Created:       Francois counillon
!
! Last modified: 24/08/2010
!
! Purpose:       Fixes EnKF output.
!
! Description:  
!            This program recompute the kfpla from the model output
!                         ensure that the sum of dp = pb
!            Input is the mem
!
! Modifications:
!
program fixenkf
use netcdf
use mod_eosfun
use nfw_mod
use m_fixhycom_eco_metno
   implicit none

   integer*4, external :: iargc
   real, parameter :: onem=98060.

   integer imem                  ! ensemble member
   character(len=80) :: oldfile,newfile, char80
   logical          :: ex
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc,cmem
   integer          :: tlevel, vlevel, nproc
   integer          :: idm,jdm,kdm
   real, allocatable:: fld(:,:)
   real, allocatable, dimension(:,:)     :: depths,modlon,modlat
   real, allocatable, dimension(:,:,:)     :: ficem,hicem
   real, allocatable, dimension(:,:,:,:)   :: dp, pb, kfpla,pbu,pbv,tmp,saln,temp
   integer, allocatable, dimension(:,:)   :: jns,ins,jnn,inn,inw,jnw,ine,jne
   real, parameter  :: epsil=1.e-12

   integer,parameter :: numfields=2
   integer :: ios,ios2
   integer :: i,j,k
   real :: dpsum
   integer, allocatable :: ns(:), nc(:)
   integer, allocatable :: ns2(:), nc2(:),ns3(:), nc3(:)
   integer :: ncid, x_ID, y_ID, z_ID, vDP_ID, vPBOT_ID, vKFP_ID
   integer :: vPBMN_ID, vPBP_ID, vPBU_ID, vPBV_ID ,vPBUP_ID,vPBVP_ID,vTMP_ID
   integer :: vFICEM_ID, vHICEM_ID
   integer :: ncid2, jns_ID, ins_ID, inw_ID, jnw_ID,jnn_ID, inn_ID, ine_ID, jne_ID
   real, allocatable :: press(:)  
   
#if defined (ECO)   
    real,dimension(:,:,:), allocatable::cfi,pressa,prsf,dpf 
    character(len=80) :: restfor
    real,dimension(:,:,:,:), allocatable::fldeco,fldeco2,dpfor,pbfor,kffor
    integer::ntracr,ktrcr,vfld_ID
    real::dpthin
    character(2)::ctrcr
    logical, dimension(:,:,:), allocatable::lcm
    !real, dimension(:,:,:), allocatable::temp,sal
    integer::kisop,nlptrc
    logical::verbose
#endif
    integer::k1,k2,k0,m 
    real::tmpfld
     
   if (iargc()==1 ) then
      call getarg(1,ctmp)
      read(ctmp,*) imem
      write(cmem,'(i3.3)') imem
   else
      print *,'"fixmycom" -- A crude routine to correct restart files obvious errors and complete diagnostic variable'
      print *
      print *,'usage: '
      print *,'   fixmicom ensemble_member'
      print *,'   "ensemble_member" is the ensemble member'
      call exit(1)
   endif

   oldfile='analysis'//cmem//'.nc'
   print *, 'fixenkf file:',oldfile
   ! Get dimensions from blkdat
   inquire(exist=ex,file=trim(oldfile))
   if (.not.ex) then
      write(*,*) 'Can not find '//'analysis'//cmem//'.nc'
      stop '(EnKF_postprocess)'
   end if
   ! Reading the restart file
   call nfw_open(trim(oldfile), nf_write, ncid)
   ! Get dimension id in netcdf file ...
   !nb total of data
   call nfw_inq_dimid(trim(oldfile), ncid, 'x', x_ID)
   call nfw_inq_dimid(trim(oldfile), ncid, 'y', y_ID)
   call nfw_inq_dimid(trim(oldfile), ncid, 'kk', z_ID)
   !nb total of track
   call nfw_inq_dimlen(trim(oldfile), ncid, x_ID, idm)
   call nfw_inq_dimlen(trim(oldfile), ncid, y_ID, jdm)
   call nfw_inq_dimlen(trim(oldfile), ncid, z_ID, kdm)

#if defined (ECO)   
   !ehouarn: bio
   k0=kdm
   k1=2
   k2=1
#else
   k0=0
   k1=1
   k2=2
#endif 
  
   print*,'k0,k1,k2 ', k0,k1,k2
  
  ! print *, 'The model dimension is :',idm,jdm,kdm
   allocate(pb (idm,jdm,2,1 ))
   allocate(dp   (idm,jdm,2*kdm,1))
   allocate(kfpla(idm,jdm,2,1 ))
   allocate(press(kdm+1))
#if defined (ECO)   
   allocate(pressa(idm,jdm,kdm+1))
#endif   
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
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'dp',vDP_ID)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'kfpla',vKFP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vDP_ID, ns, nc, dp)
   !Reading pb 
   nc(3)=2
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pb',vPBOT_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vPBOT_ID, ns, nc, pb)
!   !Reading kfpla
!   nc(3)=2
!   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'kfpla',vKFP_ID)
!   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vKFP_ID, ns, nc, kfpla)

   ! DP correction
   do j=1,jdm
   do i=1,idm
   !only if not land mask
     if (dp(i,j,1+k0,1)<100000000000.) then
      !!! Move negative layers to neighbouring layers.
      dpsum=dp(i,j,1+k0,1)
      do k = 1, kdm-1
        if (k<3) then
           ! The first two layers cannot be zero, set them to 1 m
           if (dp(i,j,k+k0,1)<=0.) then
            dp(i,j,k+1+k0,1) = dp(i,j,k+1+k0,1) + (dp(i,j,k+k0,1)-98060.0)
            dp(i,j,k+k0,1  ) = 98060.0
           end if
        else
         dp(i,j,k+1+k0,1) = dp(i,j,k+1+k0,1) + min(0.0,dp(i,j,k+k0,1))
         dp(i,j,k+k0,1  ) = max(dp(i,j,k+k0,1),0.0)
        endif
         dpsum=dpsum+ dp(i,j,k+1+k0,1) 
      end do
      !!! Go backwards to fix lowermost layer.
      do k = kdm, 4, -1
         dp(i,j,k-1+k0,1) = dp(i,j,k-1+k0,1) + min(0.0,dp(i,j,k+k0,1))
         dp(i,j,k+k0,1)   =   max(dp(i,j,k+k0,1),0.0)
      end do

#if defined (ECO)
       !ehouarn 
      pressa(i,j,1) = 0.0         
      do k = 1, kdm-1
         pressa(i,j,k+1) = pressa(i,j,k) + dp(i,j,k+k0,1)
         pressa(i,j,k+1) = min(pb(i,j,k1,1),pressa(i,j,k+1))
      end do
      pressa(i,j,kdm+1) = pb(i,j,k1,1)
#endif

     endif
   end do
   end do


!
   ! Compute the Kpfla
#if defined(ECO)
   do j=1,jdm
      do i=1,idm
          k=kdm+3
          dpsum=0.
          do while (dp(i,j,k,1).lt.epsil)
            dpsum=dpsum+dp(i,j,k,1)
            dp(i,j,k-kdm,1)=0.
            dp(i,j,k,1)=0.
            k=k+1
            if (k.gt.2*kdm) exit
          enddo
          if (k.gt.2*kdm) then
            dp(i,j,2+kdm,1)=dp(i,j,2+kdm,1)+dpsum
            dp(i,j,2,1)=dp(i,j,2+kdm,1)
          else
            dp(i,j,k,1)=dp(i,j,k,1)+dpsum
            dp(i,j,k-kdm,1)=dp(i,j,k,1)
          endif
          kfpla(i,j,1,1)=real(k-kdm)
          kfpla(i,j,2,1)=real(k-kdm)
      enddo
   enddo
   dp(:,:,1:kdm,1)=dp(:,:,kdm+1:2*kdm,1)
#else
   do j=1,jdm
      do i=1,idm
          k=3
          dpsum=0.
          do while (dp(i,j,k,1).lt.epsil)
            dpsum=dpsum+dp(i,j,k,1)
            dp(i,j,k,1)=0.
            dp(i,j,k+kdm,1)=0.
            k=k+1
            if (k.gt.kdm) exit
          enddo
          if (k.gt.kdm) then
            dp(i,j,2,1)=dp(i,j,2,1)+dpsum
            dp(i,j,2+kdm,1)=dp(i,j,2,1)
          else
            dp(i,j,k,1)=dp(i,j,k,1)+dpsum
            dp(i,j,k+kdm,1)=dp(i,j,k,1)
          endif
          kfpla(i,j,1,1)=real(k)
          kfpla(i,j,2,1)=real(k)
      enddo
   enddo
   dp(:,:,kdm+1:2*kdm,1)=dp(:,:,1:kdm,1)
#endif 
 
   pb(:,:,k2,1)=pb(:,:,k1,1)
   kfpla(:,:,k2,1)=kfpla(:,:,k1,1)


!
   ! Now we should be finished with dp pbot and kfpla
   nc(3)=2*kdm
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vDP_ID, ns, nc, dp)
   nc(3)=2
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vPBOT_ID, ns, nc, pb)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vKFP_ID, ns, nc, kfpla)
  
#ifndef ECO 
    !ehouarn: need these variables the in bio conservation  
   deallocate(dp,kfpla)
#endif
   
   !Now need to dump the 15 variables that have different name, but that almost
   !do the same stuff. Set them all equal!
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pb_mn',vPBMN_ID)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pb_p',vPBP_ID)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pbu',vPBU_ID)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pbv',vPBV_ID)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pbu_p',vPBUP_ID)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pbv_p',vPBVP_ID)
   allocate(pbu (idm,jdm,2,1    ))
   allocate(pbv (idm,jdm,2,1    ))
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vPBU_ID, ns, nc, pbu)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vPBV_ID, ns, nc, pbv)
   allocate(ins(idm,jdm))
   allocate(jns(idm,jdm))
   allocate(inn(idm,jdm))
   allocate(jnn(idm,jdm))
   allocate(inw(idm,jdm))
   allocate(jnw(idm,jdm))
   allocate(ine(idm,jdm))
   allocate(jne(idm,jdm))
  allocate(ns3(3))
  allocate(nc3(3))
  ns3(1)=1
  ns3(2)=1
  ns3(3)=1
  nc3(1)=idm
  nc3(2)=jdm
  nc3(3)=1
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vPBP_ID, ns3, nc3, pb(:,:,1,1))
   call nfw_open('grid.nc', nf_write, ncid2)
   call nfw_inq_varid('grid.nc', ncid2,'jns',jns_ID)
   call nfw_inq_varid('grid.nc', ncid2,'ins',ins_ID)
   call nfw_inq_varid('grid.nc', ncid2,'jnw',jnw_ID)
   call nfw_inq_varid('grid.nc', ncid2,'inw',inw_ID)
   allocate(ns2(2))
   allocate(nc2(2))
   ns2(1)=1
   ns2(2)=1
   nc2(1)=idm
   nc2(2)=jdm
   call nfw_get_vara_int('grid.nc', ncid2, jns_ID, ns2, nc2, jns)
   call nfw_get_vara_int('grid.nc', ncid2, ins_ID, ns2, nc2, ins)
   call nfw_get_vara_int('grid.nc', ncid2, jnw_ID, ns2, nc2, jnw)
   call nfw_get_vara_int('grid.nc', ncid2, inw_ID, ns2, nc2, inw)
   !Here we do not have acess to ifu and ilu, so we need to make some tricks to
   !make the ponts close to the land equal to 0
   do j=1,jdm
      do i=1,idm
         if ( pbu(i,j,k1,1) .ne. 0) then
            pbu(i,j,k1,1)=min(pb(i,j,k1,1),pb(inw(i,j),jnw(i,j),k1,1))
            pbu(i,j,k2,1)=pbu(i,j,k1,1)
         endif
         if ( pbv(i,j,k1,1) .ne. 0) then
            pbv(i,j,k1,1)=min(pb(i,j,k1,1),pb(ins(i,j),jns(i,j),k1,1))
            pbv(i,j,k2,1)=pbv(i,j,k1,1)
         endif
      enddo
   enddo
   nc(3)=2
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vPBU_ID, ns, nc, pbu)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vPBV_ID, ns, nc, pbv)
   deallocate(pbu,pbv,ins,jns,inw,jnw)
   !Copy to the second time level , T,S,pgfx,pgfy, [u,v]flx,,
   ns(1)=1
   ns(2)=1
   ns(3)=1
   ns(4)=1
   nc(1)=idm
   nc(2)=jdm
   nc(3)=2*kdm
   nc(4)=1
   allocate(temp   (idm,jdm,2*kdm,1))
   allocate(saln   (idm,jdm,2*kdm,1))
   allocate(tmp   (idm,jdm,2*kdm,1))
   !temp
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'temp',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, temp)
#if defined (ECO)
   temp(:,:,1:kdm,1)=temp(:,:,kdm+1:2*kdm,1)
#else
   temp(:,:,kdm+1:2*kdm,1)=temp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, temp)
   !saln
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'saln',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, saln)
#if defined (ECO)
   saln(:,:,1:kdm,1)=saln(:,:,kdm+1:2*kdm,1)
#else
   saln(:,:,kdm+1:2*kdm,1)=saln(:,:,1:kdm,1)
#endif
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, saln)
   !sigma
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'sigma',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !recompute sigma from T and S
   call eosini
   do j=1,jdm
   do i=1,idm
   do k=1,kdm
     if (temp(i,j,k+k0,1)<100000.) then
     tmp(i,j,k+k0,1)= sig(temp(i,j,k+k0,1),saln(i,j,k+k0,1))
     endif
   enddo
   enddo
   enddo
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else   
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   deallocate(temp,saln)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'u',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !v
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'v',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif 
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfy
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfy',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !uflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'uflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif 
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif 
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !utflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'utflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else   
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vtflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vtflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else 
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif 
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !usflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'usflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else 
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif    
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vsflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vsflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
#if defined (ECO)
   tmp(:,:,1:kdm,1)=tmp(:,:,kdm+1:2*kdm,1)
#else 
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
#endif   
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !!!!! 2D variables
   nc(3)=2
   !ub
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'ub',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)  
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vb
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vb',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !ubflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'ubflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vb
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vbflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pvtrop
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pvtrop',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfxm
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfxm',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfym
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfym',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xixm
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xixm',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xiym
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xiym',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xixp
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xixp',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xiyp
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xiyp',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,k2,1)=tmp(:,:,k1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   call nfw_close(trim(oldfile)//'.nc', ncid)
   
#if defined (ECO)   
   !ehouarn: remapping of the biogeoch. tracers after correction of the dp
   ! => mass conservation
   
   !get the old vertical grid
   restfor='forecast'//cmem//'.nc'
   allocate(dpfor(idm,jdm,2*kdm,1))
   allocate(pbfor(idm,jdm,2,1))
   allocate(kffor(idm,jdm,2,1 ))
   allocate(dpf(idm,jdm,kdm))
   allocate(prsf(idm,jdm,kdm+1))
   allocate(lcm(idm,jdm,kdm))
   
   inquire(exist=ex,file=trim(restfor))
   if (.not.ex) then
      write(*,*) 'Can not find '//'forecast'//cmem//'.nc'
      stop '(EnKF_postprocess)'
   end if
   ! Reading the restart file
   call nfw_open(trim(restfor), nf_write, ncid)
   ns(1)=1
   ns(2)=1
   ns(3)=1
   ns(4)=1
   nc(1)=idm
   nc(2)=jdm
   nc(3)=2*kdm
   nc(4)=1
   call nfw_inq_varid(trim(restfor),ncid,'dp',vDP_ID)
   call nfw_get_vara_double(trim(restfor),ncid, vDP_ID, ns, nc, dpfor)  
   nc(3)=2
   call nfw_inq_varid(trim(restfor),ncid,'kfpla',vDP_ID)
   call nfw_get_vara_double(trim(restfor),ncid, vDP_ID, ns, nc, kffor)
   call nfw_inq_varid(trim(restfor),ncid,'pb',vDP_ID)
   call nfw_get_vara_double(trim(restfor),ncid, vDP_ID, ns, nc, pbfor)
    
   call nfw_close(trim(restfor), ncid)
   
   dpthin=0.01
   do j=1,jdm
   do i=1,idm
      prsf(i,j,1) = 0.0
      do k = 1, kdm-1

!         prsf(i,j,k+1) = prsf(i,j,k) + max(dpfor(i,j,k+kdm,1),dpthin)
!	 prsf(i,j,k+1) =min(pbfor(i,j,2,1),prsf(i,j,k+1))
	 
!	 dpf(i,j,k)=max(dpfor(i,j,k+kdm,1),dpthin) 

         prsf(i,j,k+1) = prsf(i,j,k) + dpfor(i,j,k+kdm,1)
	 prsf(i,j,k+1) =min(pbfor(i,j,2,1),prsf(i,j,k+1))
	 
	 dpf(i,j,k)=max(dpfor(i,j,k+kdm,1),dpthin) 

	 !if(k.le.kffor(i,j,1,1))then
	 if(k.le.2)then
	 ! a regader, peut etre que k=2 suffit en sigma2
	   lcm(i,j,k)=.false.
	 else
	   lcm(i,j,k)= dpfor(i,j,k+kdm,1).le.dpthin
	 endif
      end do
      dpf(i,j,kdm)=max(dpfor(i,j,2*kdm,1),dpthin) 
      lcm(i,j,kdm)= dpfor(i,j,2*kdm,1).le.dpthin
      prsf(i,j,kdm+1)=pbfor(i,j,2,1) 

   enddo
   enddo
   
   !get the tracers and remap them
   restfor='analysisECO'//cmem//'.nc'

   inquire(exist=ex,file=trim(restfor))
   if (.not.ex) then
      write(*,*) 'Can not find '//'forecast'//cmem//'.nc'
      stop '(EnKF_postprocess)'
   end if
   ! Reading the restart file
   call nfw_open(trim(restfor), nf_write, ncid)
   ntracr=25
   ns3(1)=1
   ns3(2)=1
   ns3(3)=1
   nc3(1)=idm
   nc3(2)=jdm
   nc3(3)=kdm
   nlptrc=5
   allocate(fldeco(idm,jdm,kdm,nlptrc))
   allocate(fldeco2(idm,jdm,kdm,nlptrc))
   allocate(cfi(kdm,nlptrc,2))
   do ktrcr=1,ntracr,nlptrc
      do k=0,nlptrc-1
        call trcr_int2char(cfld,ktrcr+k) 
        call nfw_inq_varid(trim(restfor),ncid,trim(cfld),vfld_ID)
        call nfw_get_vara_double(trim(restfor),ncid, vfld_ID, ns3, nc3,fldeco(:,:,:,k+1))	
      enddo
      fldeco2(:,:,:,:)=fldeco(:,:,:,:)
      
      do j=1,jdm
      do i=1,idm
        if(dpfor(i,j,1+kdm,1)<100000000000.)then
          if(kfpla(i,j,2,1).gt.kdm)then
	    !special case: all the layers are null below the mixed layer
	    !issues with the interpolation
	    call conservation_mxlayer(kdm,nlptrc,fldeco(i,j,:,:),fldeco2(i,j,:,:),dpfor(i,j,kdm+1:2*kdm,1),&
	                            dp(i,j,kdm+1:2*kdm,1),epsil)
	    
	  else
!	    do k=1,nlptrc
!	      verbose=(i==9).and.(j==99).and.(k==1).and.(ktrcr==1)
!	      call conservation_wc3(kdm,fldeco(i,j,:,k),fldeco2(i,j,:,k),dpfor(i,j,kdm+1:2*kdm,1),&
!	          dp(i,j,kdm+1:2*kdm,1),epsil,verbose)
!            enddo
	    !computation of the weno coefficients	
	    call hybgen_weno_coefs(fldeco(i,j,:,:),dpf(i,j,1:kdm),lcm(i,j,:),cfi,kdm,nlptrc,dpthin)
	    !interpolation
	    call hybgen_weno_remap(fldeco(i,j,:,:),prsf(i,j,:),dpfor(i,j,kdm+1:2*kdm,1),cfi,fldeco2(i,j,:,:),&
                        pressa(i,j,:),kdm,kdm,nlptrc,dpthin)
	    !mixed layer	    
!	    call conservation_mxlayer(2,nlptrc,fldeco(i,j,1:2,:),fldeco2(i,j,1:2,:),dpfor(i,j,kdm+1:kdm+2,1),&
!	                            dp(i,j,kdm+1:kdm+2,1),epsil)
				    
	  endif	  
        endif 
      enddo
      enddo
      
      do k=0,nlptrc-1
        call trcr_int2char(cfld,ktrcr+k) 
        call nfw_inq_varid(trim(restfor),ncid,trim(cfld),vfld_ID)
	call nfw_put_vara_double(trim(restfor),ncid, vfld_ID, ns3, nc3,fldeco2(:,:,:,k+1))
      enddo
   enddo
   deallocate(dpf,lcm,cfi,prsf,dpfor,fldeco,fldeco2, pressa,kffor,pbfor)
   deallocate(dp,kfpla)
#endif  
 
end program
