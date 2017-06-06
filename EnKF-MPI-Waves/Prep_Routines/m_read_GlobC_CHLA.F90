module m_read_GlobC_CHLA
! Ehouarn !
! Reads CHLA1 data and grid from the GlobColour files !
contains

  subroutine read_CHLA(fname,gr,data)
  use mod_measurement
  use mod_grid
  !use m_datatest
  use m_spherdist
  use netcdf
  !use m_nf90_err
  use nfw_mod
  implicit none

  type (measurement),  intent(inout) :: data(:)
  type (grid),         intent(inout) :: gr ! CLS measurement grid
  character(len=80),   intent(in) :: fname

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, GridDepth_ID, LatLon_ID

! Variable ids
  integer :: vLatLon_ID, vNbLatitudes_ID, vNbLongitudes_ID, vLatLonMin_ID, vLatLonStep_ID, &
           vGrid0001_ID,vflagsid

! Array dimensions
  integer :: NbLatitudes, NbLongitudes, GridDepth

! Data arrays
  real,allocatable :: chla(:,:), chla_var(:,:), lon(:),lat(:)

! utilitary
  integer ncid, ijmax(2)
  real undef(1),undef_lat, undef_lon
  integer i, j,k
  logical valid
  real, parameter :: epsilon = 0.01  ! test for undefined values

! Open file
  !call nf90_err(NF90_OPEN(trim(fname),NF90_NOCLOBBER,ncid))
  call nfw_open(trim(fname), nf_nowrite, ncid)
  
! Get dimension id in netcdf file ...
  !call nf90_err(nf90_Inq_Dimid(ncid,'lat',NbLatitudes_ID))
  !call nf90_err(nf90_Inq_Dimid(ncid,'lon',NbLongitudes_ID))
  call nfw_inq_dimid(trim(fname), ncid, 'lat', NbLatitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'lon', NbLongitudes_ID)
  print*,'How far do you go'

! Get dimension length from id
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
  call nfw_inq_dimlen(trim(fname),ncid,NbLatitudes_ID,NbLatitudes)
  call nfw_inq_dimlen(trim(fname),ncid,NbLongitudes_ID,NbLongitudes)
  print*, 'Dimensions:', NbLatitudes, NbLongitudes

! State which variable you want here.. Available vars are shown when you do
  allocate(lon(NbLongitudes))
  allocate(lat(NbLatitudes))
  allocate(chla(NbLongitudes,NbLatitudes))
  allocate(chla_var(NbLongitudes,NbLatitudes))
  !allocate(chla(NbLatitudes,NbLongitudes))
  !allocate(chla_var(NbLatitudes,NbLongitudes))
  
! Variable ids in netcdf file
  !call nf90_err(nf90_inq_varid(ncid,'lat' ,vNbLatitudes_ID),'lat')
  !call nf90_err(nf90_inq_varid(ncid,'lon' ,vNbLongitudes_ID),'lon')
  !call nf90_err(nf90_inq_varid(ncid,'CHL1_mean' ,vGrid0001_ID),'CHL1_mean')
  !call nf90_err(nf90_inq_varid(ncid,'CHL1_flags' ,vflagsid),'CHL1_flags')
  call nfw_inq_varid(trim(fname),ncid,'lat' ,vNbLatitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'lon' ,vNbLongitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'CHL1_mean' ,vGrid0001_ID)
  call nfw_inq_varid(trim(fname),ncid,'CHL1_flags' ,vflagsid)

! Variable _FillValue attributes
  !call nf90_err(nf90_get_att(ncid,vGrid0001_ID ,'_FillValue',undef))
  call nfw_get_att_double(trim(fname),ncid,vGrid0001_ID ,'_FillValue',undef(1))
  print*, 'Undefined values are ',undef(1)
  gr%undef = undef(1)

! actual variable values (for dimensions of var -- see ncdump, or improve this program)
! NB: note that index dimensions are different between fortran and C internals. 
! "ncdump" gives C internals.
  print *,'test'
  !call nf90_err(nf90_get_var(ncid,vNbLongitudes_ID  ,lon))
  call nfw_get_var_double(trim(fname), ncid,vNbLongitudes_ID  ,lon)
  print *,'Range Lon', minval(lon), maxval(lon)
  !call nf90_err(nf90_get_var(ncid,vNbLatitudes_ID   ,lat))
  call nfw_get_var_double(trim(fname),ncid,vNbLatitudes_ID   ,lat)
  print *,'Range Lat', minval(lat), maxval(lat) 
  !call nf90_err(nf90_get_var(ncid,vGrid0001_ID,chla))
  call nfw_get_var_double(trim(fname),ncid,vGrid0001_ID,chla)
  print *,'Range CHLA ', minval(chla), maxval(chla)
  !call nf90_err(nf90_get_var(ncid,vflagsid,chla_var))
  call nfw_get_var_double(trim(fname),ncid,vflagsid,chla_var)
  print *,'Range CHLA_flags ', minval(chla_var), maxval(chla_var)

  !call nf90_err (nf90_close(ncid))
  call nfw_close(trim(fname), ncid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fill the data(:) vector
  
  print*,'connard',chla(629,182)
  
  do j=1,NbLongitudes ! gr%ny
  do i=1,NbLatitudes ! gr%nx
     k=(j-1)*gr%nx+i
    
     data(k)%id = 'CHLA'
    data(k)%d = chla(j,i)
    !data(k)%d = chla(i,j)

     data(k)%ipiv = i
     data(k)%jpiv = j

     data(k)%lat=lat(i)
     data(k)%lon=ang180(lon(j))

!LB: Data support is assumed = a square grid cell
!support diameter in meters stored in %a1 (tricky, isn't it ?)
!     data(k)%a1 = spherdist(lon(j)-.5*gr%dx,lat(i)-.5*gr%dy, &
!                            lon(j)+.5*gr%dx,lat(i)+.5*gr%dy)
     
     data(k)%a1 = 1
     data(k)%a2 = 0
     data(k)%a3 = 0
     data(k)%a4 = 0
     data(k)%ns = 0
     
!     data(k)%ns = 1
 
    data(k)%var = (0.35*chla(j,i))**2 ! 10 percent of the value
      !data(k)%var = (0.35*chla(i,j))**2 

     data(k)%depth=0.0

      valid =  abs( (chla(j,i)-undef(1))   / undef(1) )  > epsilon  
      !valid =  abs( (chla(i,j)-undef(1))   / undef(1) )  > epsilon  

     data(k)%status = valid

  enddo
  enddo
  print*, 'Number of data read:', k, gridpoints(gr)

!CONSISTENCY...
  !call datatest(gr,data,k,'consistency_d_'//trim(data(1)%id))

  deallocate(lat,lon,chla)
   
end subroutine read_CHLA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_CHLA_grid(filename,gr)
  !use mod_dimensions
  use mod_grid
  use netcdf
  use nfw_mod
  !use m_nf90_err
  implicit none
  
  character(len=80), intent(in) :: filename
  type(grid),        intent(out) :: gr

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, LatLon_ID


! Array dimensions
  integer :: LatLon, NbLatitudes, NbLongitudes
  real, allocatable :: dlatlon(:),lat(:),lon(:)

!variables ids
  integer vLatLonMin_ID, vLatLonStep_ID,vNbLatitudes_ID,vNbLongitudes_ID

  integer :: ncid,ierr
  real undef,scale,error_variance 

  gr = default_grid

! Open file
!filename='sst_topaz_19510.nc'
  !call nf90_err(NF90_OPEN(trim(filename),NF90_NOCLOBBER,ncid))
  call nfw_open(trim(filename), nf_nowrite, ncid)
  
! Get dimension id in netcdf file ...
  !call nf90_err(nf90_Inq_Dimid(ncid,'lat',NbLatitudes_ID))
  !call nf90_err(nf90_Inq_Dimid(ncid,'lon',NbLongitudes_ID))
  call nfw_inq_dimid(trim(filename), ncid, 'lat', NbLatitudes_ID)
  call nfw_inq_dimid(trim(filename), ncid, 'lon', NbLongitudes_ID)
  
  print*,'How far do you go'
! Get dimension length from id
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
  call nfw_inq_dimlen(trim(filename),ncid,NbLatitudes_ID,NbLatitudes)
  call nfw_inq_dimlen(trim(filename),ncid,NbLongitudes_ID,NbLongitudes)
  print*, 'Dimensions:', NbLatitudes, NbLongitudes

  allocate(dlatlon(2))  ! dx and dy

! Variables in NetCDF file
  !call nf90_err(nf90_inq_varid(ncid,'lat' ,vNbLatitudes_ID),'lat')
  !call nf90_err(nf90_inq_varid(ncid,'lon' ,vNbLongitudes_ID),'lon')
  call nfw_inq_varid(trim(filename),ncid,'lat' ,vNbLatitudes_ID)
  call nfw_inq_varid(trim(filename),ncid,'lon' ,vNbLongitudes_ID)
  
  allocate(lon(NbLongitudes))
  allocate(lat(NbLatitudes))
  !call nf90_err(nf90_get_var(ncid,vNbLongitudes_ID,lon))
  !call nf90_err(nf90_get_var(ncid,vNbLatitudes_ID,lat))
  call nfw_get_var_double(trim(filename), ncid,vNbLongitudes_ID  ,lon)
  call nfw_get_var_double(trim(filename),ncid,vNbLatitudes_ID   ,lat)
  
  !call nf90_err(nf90_get_att(ncid,NF90_GLOBAL ,'lat_step',dlatlon(1)))
  !call nf90_err(nf90_get_att(ncid,NF90_GLOBAL ,'lon_step',dlatlon(2)))
  call nfw_get_att_double(trim(filename),ncid,nf_global ,'lat_step',dlatlon(1))
  call nfw_get_att_double(trim(filename),ncid,nf_global ,'lon_step',dlatlon(2))
 
  !call nf90_err (nf90_close(ncid))
  call nfw_close(trim(filename), ncid)
  
  gr%nx=NbLatitudes
  gr%ny=NbLongitudes
  gr%x0=   lat(NbLatitudes)
  gr%y0=   lon(1)
!  gr%dx=   0.179
  gr%dx=   dlatlon(1)
  gr%dy=   dlatlon(2)
  gr%reg = .true.
  gr%order = 2
  gr%ux = 'deg'
  gr%uy = 'deg'
  gr%set = .true.

  deallocate(lat,lon,dlatlon)

 end subroutine read_CHLA_grid
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
end module m_read_GlobC_CHLA
