program pert_param
! use netcdf
use mod_eosfun
use nfw_mod
use m_sample2D
use m_param_ensemble
   implicit none

   integer imem                  ! ensemble member
   character(len=80) :: oldfile,newfile, char80
   logical          :: ex
   character(len=8) :: cfld, ctmp
   character(len=4) :: distrib
    character(len=3) :: cmem
   integer          :: idm,jdm,nrens
   real, allocatable, dimension(:,:)     :: fld,fld0
   real, allocatable, dimension(:,:,:) :: work

   integer :: n1,n2
   integer :: i,j,k,x_ID,y_ID,h_ID,ncid
   integer, allocatable :: ns(:), nc(:)
   integer, allocatable :: ns2(:), nc2(:),ns3(:), nc3(:)
 
   real::rv,sd_d,rh

!   if (iargc()==1 ) then
   if (command_argument_count() == 1) then
      call getarg(1,ctmp)
      read(ctmp,*) nrens
   else
      print *,'pert_param'
      print *
      print *,'usage: '
      print *,'   pert_param ensemble_member'
      print *,'   "ensemble_member" is the ensemble member'
      call exit(1)
   endif
   
   open(11,file='pert_param.in',action='read', status='old')
   read(11,*) sd_d   ! Std dev of log(d), no unit, 0.1 = 10%
!   read(11,*) rv     ! Vertical correlation range (in nb layers)
   read(11,*) rh     ! horizontal correlation range (in nb of grid cells)
   read(11,*) distrib  
   close(11) 
   
   oldfile='forecast001.nc'
   print *, 'perturb_param file:',oldfile
   ! Get dimensions from blkdat
   inquire(exist=ex,file=trim(oldfile))
   if (.not.ex) then
      write(*,*) 'Can not find '//'forecast'//cmem//'.nc'
      stop '(pertub_param)'
   end if
   ! Reading the restart file
   call nfw_open(trim(oldfile), nf_write, ncid)
   ! Get dimension id in netcdf file ...
   !nb total of data
   call nfw_inq_dimid(trim(oldfile), ncid, 'xi_rho', x_ID)
   call nfw_inq_dimid(trim(oldfile), ncid, 'eta_rho', y_ID)
   !nb total of track
   call nfw_inq_dimlen(trim(oldfile), ncid, x_ID, idm)
   call nfw_inq_dimlen(trim(oldfile), ncid, y_ID, jdm)
   call nfw_close(trim(oldfile), ncid)
   
   allocate(fld (idm,jdm))
   allocate(work(idm,jdm,nrens))   
   allocate(ns(2))
   allocate(nc(2))
   ns(1)=1
   ns(2)=1
   nc(1)=idm
   nc(2)=jdm
   
   
   print*,'dimenssions:',idm,jdm,nrens
   call perturb_param(work,nrens,idm,jdm,rh,sd_d,trim(distrib))
   
   do k=1,nrens
      write(cmem,'(i3.3)') k
     
      oldfile='forecast'//cmem//'.nc'
      print *, 'perturb_param file:',trim(oldfile)
      ! Get dimensions from blkdat
      inquire(exist=ex,file=trim(oldfile))
      if (.not.ex) then
        write(*,*) 'Can not find '//'forecast'//cmem//'.nc'
        stop '(pertub_param)'
      end if
      ! Reading the restart file
      call nfw_open(trim(oldfile), nf_write, ncid)
      print *,'Reading h'
      call nfw_inq_varid(trim(oldfile), ncid,'h',h_ID)
      call nfw_get_vara_double(trim(oldfile), ncid, h_ID, ns, nc, fld)
      
      if(trim(distrib)=='logn')then
        !fld(:,:)=max(0.2,fld(:,:)*work(:,:,k))
	fld(:,:)=fld(:,:)*work(:,:,k)
	!print*, work(:,:,k)
      elseif(trim(distrib)=='norm')then
        fld(:,:)=fld(:,:)+work(:,:,k)	
      end if	

!
    ! Now we should be finished with dp pbot and kfpla
      call nfw_put_vara_double(trim(oldfile), ncid, h_ID, ns, nc, fld)
      call nfw_close(trim(oldfile), ncid)
   
  enddo
   



end program pert_param


