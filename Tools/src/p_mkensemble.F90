program mkensemble
! use netcdf
use mod_eosfun
use nfw_mod
use m_pseudo2D
   implicit none

   integer*4, external :: iargc
   real, parameter :: onem=9806.

   integer imem                  ! ensemble member
   character(len=80) :: oldfile,newfile, char80
   logical          :: ex
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc,cmem
   integer          :: tlevel, vlevel, nproc
   integer          :: idm,jdm,kdm
   real, allocatable:: fld(:,:)
   real, allocatable, dimension(:,:)     :: depths,modlon,modlat,fld2d,ranfld,accranfld,fld0
   real, allocatable, dimension(:,:,:)     :: ficem,hicem
   real, allocatable, dimension(:,:,:,:)   :: dp, pb, kfpla,pbu,pbv,tmp,saln,temp
   integer, allocatable, dimension(:,:)   :: jns,ins,jnn,inn,inw,jnw,ine,jne
   real, parameter  :: epsil=1.e-11

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
   real::rv,sd_d,rh,alp2,bet
   integer::n1,n2

!   if (iargc()==1 ) then
   if (command_argument_count() /= 1) then
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
   
   open(11,file='mkensemble.in',action='read', status='old')
   read(11,*) sd_d   ! Std dev of log(d), no unit, 0.1 = 10%
   read(11,*) rv     ! Vertical correlation range (in nb layers)
   read(11,*) rh     ! horizontal correlation range (in nb of grid cells)
   close(11)
   
   oldfile='analysis'//cmem//'.nc'
   print *, 'mkensemble file:',oldfile
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
  ! print *, 'The model dimension is :',idm,jdm,kdm
   allocate(pb (idm,jdm,2,1 ))
   allocate(dp   (idm,jdm,2*kdm,1))
   allocate(kfpla(idm,jdm,2,1 ))
   allocate(press(kdm+1))
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
   

   !  Generates the vertical correlation of the ensemble.
   alp2=exp(-1.0/rv)
   print *,'Alp=',alp2
   bet=sqrt(1.0-alp2**2) ! keeps var(enstmp)=var(ensmem)
   print *,'Bet=',bet

   ! original Hycom state kept as first member in record 1 - 
   ! FFT dimensions
   n1=2**(ceiling(log(real(idm))/log(2.)))
   n2=2**(ceiling(log(real(jdm))/log(2.)))

   ! Process dp layers -- dp must be the same for both time steps,
   ! otherwise temporal gradients may cause problems
   allocate(fld2d(idm,jdm))
   allocate(ranfld(idm,jdm))
   allocate(accranfld(idm,jdm))
   allocate(fld0(idm,jdm))
   print '(a)','-perturbing layers '
   fld0(:,:)=dp(:,:,1,1)
   do k=1,kdm 
    ! 3D vertically correlated
    fld2d(:,:)=dp(:,:,k,1)
    where(fld0>100000000000.) fld2d=0.
    call pseudo2D(ranfld,idm,jdm,1,rh,n1,n2)
    if (k==1) then
       call pseudo2D(accranfld,idm,jdm,1,rh,n1,n2)
    else
       accranfld=alp2*accranfld+bet*ranfld
    end if
    !print *,'dp test unperturbed',k,minval(fld2d)/onem,maxval(fld2d)/onem
    fld2d=fld2d*exp(accranfld*sd_d-0.5*sd_d**2)
    ! We can afford one 3D variable (about 60 MB for a 800x880x22 grid)
    dp(:,:,k,1)=fld2d
   end do
   
   ! DP correction
   do j=1,jdm
   do i=1,idm
   !only if not land mask
     if (dp(i,j,1,1)<100000000000.) then
      !!! Move negative layers to neighbouring layers.
      dpsum=dp(i,j,1,1)
      do k = 1, kdm-1
        if (k<3) then
           ! The first two layers cannot be zero, set them to 1 m
           if (dp(i,j,k,1)<=0.) then
            dp(i,j,k+1,1) = dp(i,j,k+1,1) + (dp(i,j,k,1)-98060.0)
            dp(i,j,k,1  ) = 98060.0
           end if
        else
         dp(i,j,k+1,1) = dp(i,j,k+1,1) + min(0.0,dp(i,j,k,1))
         dp(i,j,k,1  ) = max(dp(i,j,k,1),0.0)
        endif
         dpsum=dpsum+ dp(i,j,k+1,1) 
      end do
      !!! Go backwards to fix lowermost layer.
      do k = kdm, 4, -1
         dp(i,j,k-1,1) = dp(i,j,k-1,1) + min(0.0,dp(i,j,k,1))
         dp(i,j,k,1)   =   max(dp(i,j,k,1),0.0)
      end do

!      !!! No layers below the sea bed.
!      press(  1) = 0.0         
!      do k = 1, kdm-1
!         press(k+1) = press(k) + dp(i,j,k,1)
!         press(k+1) = min(pb(i,j,1,1),press(k+1))
!      end do
!      press(kdm+1) = pb(i,j,1,1)
!
!      dpsum=0
!      do k = 1, kdm
!         dp(i,j,k,1) = press(k+1) - press(k)
!         dpsum=dpsum+dp(i,j,k,1)
!      end do
!!
!      if (dpsum-pb(i,j,1,1)>1 .and. pb(i,j,1,1)<100000000000 ) then
!        print *,'Inconsistency at point:',i,j,dpsum,pb(i,j,1,1)
!      endif
     endif
   end do
   end do


!
   ! Compute the Kpfla
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
   pb(:,:,2,1)=pb(:,:,1,1)
   kfpla(:,:,2,1)=kfpla(:,:,1,1)


!
   ! Now we should be finished with dp pbot and kfpla
   nc(3)=2*kdm
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vDP_ID, ns, nc, dp)
   nc(3)=2
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vPBOT_ID, ns, nc, pb)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vKFP_ID, ns, nc, kfpla)
   deallocate(dp,kfpla)
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
         if ( pbu(i,j,1,1) .ne. 0) then
            pbu(i,j,1,1)=min(pb(i,j,1,1),pb(inw(i,j),jnw(i,j),1,1))
            pbu(i,j,2,1)=pbu(i,j,1,1)
         endif
         if ( pbv(i,j,1,1) .ne. 0) then
            pbv(i,j,1,1)=min(pb(i,j,1,1),pb(ins(i,j),jns(i,j),1,1))
            pbv(i,j,2,1)=pbv(i,j,1,1)
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
   temp(:,:,kdm+1:2*kdm,1)=temp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, temp)
   !saln
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'saln',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, saln)
   saln(:,:,kdm+1:2*kdm,1)=saln(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, saln)
   !sigma
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'sigma',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !recompute sigma from T and S
   call eosini
   do j=1,jdm
   do i=1,idm
   do k=1,kdm
     if (temp(i,j,k,1)<100000.) then
     tmp(i,j,k,1)= sig(temp(i,j,k,1),saln(i,j,k,1))
     endif
   enddo
   enddo
   enddo
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   deallocate(temp,saln)
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'u',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !v
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'v',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfy
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfy',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !uflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'uflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !utflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'utflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vtflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vtflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !usflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'usflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vsflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vsflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,kdm+1:2*kdm,1)=tmp(:,:,1:kdm,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !!!!! 2D variables
   nc(3)=2
   !ub
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'ub',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vb
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vb',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !ubflx
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'ubflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !vb
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'vbflx',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pvtrop
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pvtrop',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfxm
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfxm',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !pgfym
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'pgfym',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xixm
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xixm',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xiym
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xiym',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xixp
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xixp',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   !xiyp
   call nfw_inq_varid(trim(oldfile)//'.nc', ncid,'xiyp',vTMP_ID)
   call nfw_get_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   tmp(:,:,2,1)=tmp(:,:,1,1)
   call nfw_put_vara_double(trim(oldfile)//'.nc', ncid, vTMP_ID, ns, nc, tmp)
   call nfw_close(trim(oldfile)//'.nc', ncid)

end program mkensemble


