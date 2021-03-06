program p_prep_obs_OOPS
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
  real,dimension(:),allocatable::noise1d
  real,dimension(:,:,:),allocatable :: fld
  integer::nechx,nechy

  !ehouarn: length of the assimilation window
  integer, parameter::lasswin=1!120!3!12
  logical ex
  integer:: ncid, x_ID, y_ID,l_ID,nl,var_id,cpt
  
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

  inquire(file='observation-OOPS-QG.nc',exist=ex)
  if (ex) then
    ! Reading the grid file
    call nfw_open('observation-OOPS-QG.nc', nf_nowrite, ncid)
    ! Get dimension id in netcdf file ...
    call nfw_inq_dimid('observation-OOPS-QG.nc', ncid, 'nx', x_ID)
    call nfw_inq_dimid('observation-OOPS-QG.nc', ncid, 'ny', y_ID)
    call nfw_inq_dimid('observation-OOPS-QG.nc', ncid, 'nl', l_ID)
    !Get the dimension
    call nfw_inq_dimlen('observation-OOPS-QG.nc', ncid, x_ID, nx)
    call nfw_inq_dimlen('observation-OOPS-QG.nc', ncid, y_ID, ny)
    call nfw_inq_dimlen('observation-OOPS-QG.nc', ncid, l_ID, nl)
    
    !get the field
    allocate(fld(nl,ny,nx))
    call nfw_inq_varid('observation-OOPS-QG.nc', ncid,'stream',var_id)
    call nfw_get_var_double('observation-OOPS-QG.nc', ncid, var_id, fld)
    
    call nfw_close('observation-OOPS-QG.nc',ncid)
    
    print *, 'The model dimension is :',nx,ny,nl
  else
      stop 'ERROR: file .nc is missing'
  endif
   

  dosuperob = .false.
  is3d = .false.

  ! Fill the "data" array by calling subroutines specific for the producer
  ! and observation type

  dosuperob = .false.
  grpoints=(floor(real(nx)/real(nechx))-1)*(floor(real(ny)/real(nechy))-1)
  print*,'nb obs in OOPS: ',grpoints
  allocate(obs(grpoints))
  allocate(noise1d(grpoints))
  call set_random_seed3 
  call randn(grpoints,noise1d)
  
  cpt=0      
  do i=1+nechx,nx-nechx,nechx
    do j=1+nechy ,ny-nechy ,nechy 
       cpt=cpt+1
       obs(cpt)%lat=0
       obs(cpt)%lon=0
       obs(cpt)%depth=0
       obs(cpt) % var=0.01
       obs(cpt)%d=fld(1,j,i)+noise1d(cpt)*obs(cpt) % var
       obs(cpt) %date=0
       obs(cpt) %ipiv=i
       obs(cpt) %jpiv=j
       obs(cpt)%orig_id =0
       obs(cpt)%i_orig_grid = -1
       obs(cpt)%j_orig_grid = -1
       obs(cpt)%h = 1
       obs(cpt)%status = .true.
       obs(cpt)%a1 = 1
       obs(cpt)%a2 = 0
       obs(cpt)%a3 = 0
       obs(cpt)%a4 = 0
       obs(cpt)%ns = 0
       obs(cpt)%date = 0
       obs(cpt)%id ='STR'
    end do
  end do   
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
  
  !call nc_obs2d(obs,nrobs,nx,ny,depths)
  
  deallocate(obs)

  print *, 'prep_obs: end of processing'
end program p_prep_obs_OOPS


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


