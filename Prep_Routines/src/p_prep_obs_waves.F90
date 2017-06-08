program p_prep_obs_waves
  use mod_measurement
  use m_uobs
  use nfw_mod
  use m_set_random_seed2
  use m_write_wet_file
  use m_random
  implicit none

  integer, parameter :: STRLEN = 512

  type (measurement), allocatable :: obs(:)

  integer :: nx, ny,nmonth
  integer, parameter :: maxobs = 5000000
  character(STRLEN) :: fname, fnamehdr, dataformat, producer,fnamehdr2
  character(len=3) :: form
  character(len=5) :: obstype

  integer :: nrobs
  integer :: grpoints, k
  real :: factor, var
  real :: rdummy

  logical :: data_eq_obs

  ! superobs
  logical :: dosuperob
  logical :: is3d
  integer :: nrsobs
  type(measurement), allocatable :: sobs(:)

  integer :: i,j
  integer :: nthisobs
  integer, allocatable, dimension(:) :: thisobs
  real,dimension(:),allocatable::noise1d,xi,eta
  real,dimension(:,:),allocatable :: fld,depths
  integer::nechx,nechy

  !ehouarn: length of the assimilation window
  integer, parameter::lasswin=1!120!3!12
  logical ex
  integer:: ncid, x_ID, y_ID,l_ID,nl,var_id,cpt,ntype_obs
  character(len=3) :: type_obs(2)
  
  data_eq_obs = .false.

  open(10, file = 'infile.data')
  read(10, '(a)') producer
  read(10, '(a)') obstype
  read(10, '(a)') fnamehdr
  read(10, '(a)') fname
  close(10)

  read(fnamehdr,'(i2)')nechx
  read(fname,'(i2)')nechy

  print *, 'Data producer: ', trim(producer)
  print *, 'Data to be processed for TOPAZ: ', trim(obstype)
  print *, 'Filenames to be processed are: "', trim(fnamehdr), '" "', trim(fname), '"'
  print *, 'Result of processing is stored in temporary file "observation.uf"'

  ! Get grid dimensions 

  inquire(file='shoreface_obs.nc',exist=ex)
  if (ex) then
    ! Reading the grid file
    call nfw_open('shoreface_obs.nc', nf_nowrite, ncid)
    ! Get dimension id in netcdf file ...
    call nfw_inq_dimid('shoreface_obs.nc', ncid, 'xi_rho', x_ID)
    call nfw_inq_dimid('shoreface_obs.nc', ncid, 'eta_rho', y_ID)
!    call nfw_inq_dimid('shoreface_obs.nc', ncid, 'time', l_ID)
    !Get the dimension
    call nfw_inq_dimlen('shoreface_obs.nc', ncid, x_ID, nx)
    call nfw_inq_dimlen('shoreface_obs.nc', ncid, y_ID, ny)
!    call nfw_inq_dimlen('shoreface_obs.nc', ncid, l_ID, nl)
    
    !get the field
    !allocate(fld(ny,nx),xi(nx),eta(ny))
    !allocate(depths(ny,nx))
    allocate(fld(nx,ny),xi(nx),eta(ny))
    allocate(depths(nx,ny))
    ! depths
    call nfw_inq_varid('shoreface_obs.nc', ncid,'h',var_id)
    call nfw_get_var_double('shoreface_obs.nc', ncid, var_id, depths)
    
    if(trim(obstype)=='EPB') then
      call nfw_inq_varid('shoreface_obs.nc', ncid,'epb',var_id)
    elseif(trim(obstype)=='OCG') then
      call nfw_inq_varid('shoreface_obs.nc', ncid,'Cg',var_id)
    end if  
    call nfw_get_var_double('shoreface_obs.nc', ncid, var_id, fld)
    
    call nfw_inq_varid('shoreface_obs.nc', ncid,'xi_rho',var_id)
    call nfw_get_var_double('shoreface_obs.nc', ncid, var_id, xi)
    
    call nfw_inq_varid('shoreface_obs.nc', ncid,'eta_rho',var_id)
    call nfw_get_var_double('shoreface_obs.nc', ncid, var_id, eta)
    
    call nfw_close('shoreface_obs.nc',ncid)
    
    print *, 'The model dimension is :',nx,ny
  else
      stop 'ERROR: file .nc is missing'
  endif
   

  dosuperob = .false.
  is3d = .false.

  ! Fill the "data" array by calling subroutines specific for the producer
  ! and observation type
  

  dosuperob = .false.
!  grpoints=(floor(real(nx)/real(nechx))-1)*(floor(real(ny)/real(nechy))-1)
  grpoints=nx*ny
  print*,'nb obs in shoreface: ',grpoints
  allocate(obs(grpoints))
  allocate(noise1d(grpoints))
  call set_random_seed3 
  call randn(grpoints,noise1d)
  
  cpt=0   
  
!    do i=1+nechx,nx-nechx,nechx
!      do j=1+nechy ,ny-nechy ,nechy 
    do i=1,nx
       do j=1 ,ny
       print*,trim(obstype),fld(i,j)    
       cpt=cpt+1
!       obs(cpt)%lat=10+20*(eta(j)-1)
!       obs(cpt)%lon=10+20*(xi(i)-1)
!!!       obs(cpt)%lat=eta(j)
!!!       obs(cpt)%lon=xi(i)
       obs(cpt)%lat=xi(i)
       obs(cpt)%lon=eta(j)
       obs(cpt)%depth=0
       obs(cpt) % var=0.01
       obs(cpt)%d=fld(i,j)+noise1d(cpt)*obs(cpt) % var
       obs(cpt) %date=0
       obs(cpt) %ipiv=i
       obs(cpt) %jpiv=j
       obs(cpt)%orig_id =0
 !      obs(cpt)%i_orig_grid = -1
 !      obs(cpt)%j_orig_grid = -1
       obs(cpt)%h = 1
       obs(cpt)%status = .true.
 !      obs(cpt)%status = .false.	  
       obs(cpt)%a1 = 1
       obs(cpt)%a2 = 0
       obs(cpt)%a3 = 0
       obs(cpt)%a4 = 0
       obs(cpt)%ns = 0
       obs(cpt)%date = 0
       obs(cpt)%id = trim(obstype)
      end do
    end do  
    
   ! do j=1,nx*ny
   !    print*, obs(j)%d,fld(obs(j) %jpiv,obs(j) %ipiv)
   !    print*,obs(j) %jpiv,obs(j) %ipiv
   ! end do
    
  deallocate(fld) 
  deallocate(noise1d)
  data_eq_obs = .true.
  print *,'Data read'
  
  nrobs=grpoints
  
  if (nrobs .ge. maxobs) then
     print *, 'max No. of data reached, increase it!'
     stop 'ERROR: p_prep_obs'
  elseif (nrobs .le. 1) then
     print *, 'less than one observation in the whole dataset'
     !PS 4/9/2011 stop 'ERROR: p_prep_obs: Not worth the money'
  end if

  ! Write data to the binary file "observations.uf"
  !
  call write_wet_file(obs, nrobs)

  call uobs_get(obs(1 : nrobs) % id, nrobs, .true.)
  allocate(thisobs(nrobs))
  do i = 1, nuobs
     nthisobs = 0
     do k = 1, nrobs
        if (trim(unique_obs(i)) == trim(obs(k) % id)) then
           nthisobs = nthisobs + 1
           thisobs(nthisobs) = k
        end if
     end do

     if (nthisobs > 0) then
        call obs2nc(nthisobs, obs(thisobs(1 : nthisobs)))
     end if
  end do
  deallocate(thisobs)

  print *, 'Last observation:'
  print '(a)','   obs       var    id      lon   lat  depth   ipiv  jpiv   nsup'//&
         '  4-bilin-coeffs    active  orig (i,j)   dp    age orig_id'
  print '(2g10.2,a6,3f6.1,3i6,4f5.1,l5,2i7,f7.1,2i5)', obs(nrobs)
  
  call nc_obs2d(obs,nrobs,nx,ny,depths)
  
  deallocate(obs,depths)

  print *, 'prep_obs: end of processing'
end program p_prep_obs_waves


subroutine obs2nc(nobs, obs)
  use mod_measurement
  use nfw_mod
  implicit none

  integer, parameter :: STRLEN = 512

  integer, intent(in) :: nobs
  type(measurement), intent(in) :: obs(nobs)

  character(STRLEN) :: fname
  integer :: ncid, obsdimid(1), lon_id, lat_id, depth_id, d_id, var_id, age_id
  integer :: n_id, ipiv_id, jpiv_id
  integer :: n(nobs)

  ! Create netcdf file of observations
  !
  write(fname, '(a, a, a)') 'observations-', trim(obs(1) % id), '.nc'
  print *, 'dumping observations to "', trim(fname), '"'

  call nfw_create(fname, nf_clobber, ncid)

  call nfw_def_dim(fname, ncid, 'nobs', nobs, obsdimid(1))
  call nfw_def_var(fname, ncid, 'lon', nf_float, 1, obsdimid(1), lon_id)
  call nfw_def_var(fname, ncid,  'lat', nf_float, 1, obsdimid(1), lat_id)
  call nfw_def_var(fname, ncid, 'depth', nf_float, 1, obsdimid(1), depth_id)
  call nfw_def_var(fname, ncid, 'd', nf_float, 1, obsdimid(1), d_id)
  call nfw_def_var(fname, ncid, 'var', nf_float, 1, obsdimid(1), var_id)
  call nfw_def_var(fname, ncid, 'age', nf_int, 1, obsdimid(1), age_id)
  call nfw_def_var(fname, ncid, 'n', nf_int, 1, obsdimid(1), n_id)
  call nfw_def_var(fname, ncid, 'ipiv', nf_int, 1, obsdimid(1), ipiv_id)
  call nfw_def_var(fname, ncid, 'jpiv', nf_int, 1, obsdimid(1), jpiv_id)
  call nfw_enddef(fname, ncid)

  call nfw_put_var_double(fname, ncid, lon_id, obs(1:nobs) % lon)
  call nfw_put_var_double(fname, ncid, lat_id, obs(1:nobs) % lat)
  call nfw_put_var_double(fname, ncid, depth_id, obs(1:nobs) % depth)
  call nfw_put_var_double(fname, ncid, d_id, obs(1:nobs) % d)
  call nfw_put_var_double(fname, ncid, var_id, obs(1:nobs) % var)
  call nfw_put_var_int(fname, ncid, age_id, obs(1:nobs) % date)
  call nfw_put_var_int(fname, ncid, ipiv_id, obs(1:nobs) % ipiv)
  call nfw_put_var_int(fname, ncid, jpiv_id, obs(1:nobs) % jpiv)
  n = int(obs(1:nobs) % h)
  call nfw_put_var_int(fname, ncid, n_id, n)
  
  call nfw_close(fname, ncid)
end subroutine obs2nc

 
  subroutine nc_obs2d(obs,nobs,nx,ny,depths)
    use mod_measurement
    !use netcdf
    use nfw_mod
    implicit none
    integer::nobs,nx,ny
    type (measurement)::obs(nobs)
    real,dimension(nx,ny)::depths
    real,dimension(nx,ny)::obs2d
    integer::k,ierr,ncid,dimnx,dimny,var_id,var_id2
    character(len=80) :: fname   
    
    fname='Observations_2D.nc'
    
    obs2d(:,:)=-1.
    do k=1,nobs
      obs2d(obs(k)%ipiv, obs(k)%jpiv)=obs(k)%d 
    enddo
    
    
    call nfw_create(trim(fname), nf_write, ncid)
        
    call nfw_def_dim(trim(fname), ncid, 'nx', nx, dimnx)
    call nfw_def_dim(trim(fname), ncid, 'ny', ny, dimny)    
    !ierr=NF90_DEF_DIM(ncid,'nx',nx,dimnx)
    !ierr=NF90_DEF_DIM(ncid,'ny',ny,dimny)

    !ierr=NF90_REDEF(ncid)
    !ierr=NF90_DEF_VAR(ncid,'CHLA',NF90_Float,(/dimnx,dimny/),var_id)
    !ierr=NF90_ENDDEF(ncid)
    !ierr=NF90_PUT_VAR(ncid,var_id,obs2d(1:nx,1:ny))
    
    call nfw_def_var(trim(fname), ncid, 'obs', nf_float,2,(/dimnx,dimny/), var_id)
    call nfw_def_var(trim(fname), ncid, 'depths', nf_float,2,(/dimnx,dimny/), var_id2)
    call nfw_enddef(trim(fname), ncid)

    call nfw_put_var_double(trim(fname), ncid, var_id,obs2d(1:nx,1:ny))
    call nfw_put_var_double(trim(fname), ncid, var_id2,depths(1:nx,1:ny))
    !ierr=NF90_REDEF(ncid)
    !ierr=NF90_DEF_VAR(ncid,'depths',NF90_Float,(/dimnx,dimny/),var_id)
    !ierr=NF90_ENDDEF(ncid)
    !ierr=NF90_PUT_VAR(ncid,var_id,depths(1:nx,1:ny))

    !ierr=NF90_CLOSE(ncid)
    call nfw_close(trim(fname), ncid)
    
 end subroutine nc_obs2d
 

