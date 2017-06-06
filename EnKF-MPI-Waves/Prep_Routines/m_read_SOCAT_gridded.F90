module m_read_SOCAT_gridded
! Ehouarn !
! Reads CHLA1 data and grid from the GlobColour files !
contains

  subroutine read_SOCAT_FCO2(fname,nmonth,gr,data)
  use mod_measurement
  use mod_grid
  !use m_datatest
  use m_spherdist
  use netcdf
  !use m_nf90_err
  use nfw_mod
  implicit none

  type (measurement),  intent(inout) :: data(:)
  character(len=80),   intent(in) :: fname
  type(grid),        intent(out) :: gr
  integer::nmonth

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, GridDepth_ID, LatLon_ID, NbTime_id

! Variable ids
  integer :: vLatLon_ID, vNbLatitudes_ID, vNbLongitudes_ID, vLatLonMin_ID, vLatLonStep_ID, &
           vGrid0001_ID,vflagsid

! Array dimensions
  integer :: NbLatitudes, NbLongitudes, GridDepth, NbTime

! Data arrays
  real,allocatable :: FCO2(:,:,:), FCO2_var(:,:,:), lon(:),lat(:)

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
  call nfw_inq_dimid(trim(fname), ncid, 'YLAT', NbLatitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'XLON', NbLongitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'TMNTH', Nbtime_ID)
  print*,'How far do you go'

! Get dimension length from id
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
  call nfw_inq_dimlen(trim(fname),ncid,NbLatitudes_ID,NbLatitudes)
  call nfw_inq_dimlen(trim(fname),ncid,NbLongitudes_ID,NbLongitudes)
  call nfw_inq_dimlen(trim(fname),ncid,Nbtime_ID,NbTime)
  print*, 'Dimensions:', NbLatitudes, NbLongitudes,NbTime

! State which variable you want here.. Available vars are shown when you do
  allocate(lon(NbLongitudes))
  allocate(lat(NbLatitudes))
  !allocate(FCO2(NbTime,NbLongitudes,NbLatitudes))
  !allocate(FCO2_var(NbTime,NbLongitudes,NbLatitudes))
  !allocate(FCO2(NbTime,NbLatitudes,NbLongitudes))
  !allocate(FCO2_var(NbTime,NbLatitudes,NbLongitudes))
   allocate(FCO2(NbLongitudes,NbLatitudes, NbTime))
  allocate(FCO2_var(NbLongitudes,NbLatitudes, NbTime))
  
! Variable ids in netcdf file
  !call nf90_err(nf90_inq_varid(ncid,'lat' ,vNbLatitudes_ID),'lat')
  !call nf90_err(nf90_inq_varid(ncid,'lon' ,vNbLongitudes_ID),'lon')
  !call nf90_err(nf90_inq_varid(ncid,'CHL1_mean' ,vGrid0001_ID),'CHL1_mean')
  !call nf90_err(nf90_inq_varid(ncid,'CHL1_flags' ,vflagsid),'CHL1_flags')
  call nfw_inq_varid(trim(fname),ncid,'YLAT' ,vNbLatitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'XLON' ,vNbLongitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'FCO2_AVE_UNWTD' ,vGrid0001_ID)
  call nfw_inq_varid(trim(fname),ncid,'FCO2_STD_UNWTD' ,vflagsid)

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
  call nfw_get_var_double(trim(fname),ncid,vGrid0001_ID,FCO2)
  print *,'Range FCO2 ', minval(FCO2), maxval(FCO2)
  !call nf90_err(nf90_get_var(ncid,vflagsid,chla_var))
  call nfw_get_var_double(trim(fname),ncid,vflagsid,FCO2_var)
  print *,'Range FCO2_var ', minval(FCO2_var), maxval(FCO2_var)

  !call nf90_err (nf90_close(ncid))
  call nfw_close(trim(fname), ncid)  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fill the data(:) vector
  k=0
  do j=1,NbLongitudes ! gr%ny
  do i=1,NbLatitudes ! gr%nx
    
    k=k+1
    data(k)%id = 'FCO2'
    data(k)%d = FCO2(j,i,nmonth)
    !data(k)%d = chla(i,j)

     data(k)%ipiv = i
     data(k)%jpiv = j

     data(k)%lat=lat(i)
     data(k)%lon=ang180(lon(j))

!LB: Data support is assumed = a square grid cell
!support diameter in meters stored in %a1 (tricky, isn't it ?)
      !data(k)%a1 = -1.0e10 
      !data(k)%a2 = -1.0e10 
      !data(k)%a3 = -1.0e10 
      !data(k)%a4 = -1.0e10 
      !data(k)%ns = 0
     data(k)%a1 = spherdist(lon(j)-.5*gr%dx,lat(i)-.5*gr%dy, &
                            lon(j)+.5*gr%dx,lat(i)+.5*gr%dy)
     data(k)%ns = 1
    
     data(k)%var = FCO2_var(j,i,nmonth)

     data(k)%depth=0.0

      valid =  abs( (FCO2(j,i,nmonth)-undef(1))   / undef(1) )  > epsilon  
      !valid =  abs( (chla(i,j)-undef(1))   / undef(1) )  > epsilon  

     data(k)%status = valid
  enddo
  enddo
  print*, 'Number of data read:', k, gridpoints(gr)

!CONSISTENCY...
  !call datatest(gr,data,k,'consistency_d_'//trim(data(1)%id))

  deallocate(lat,lon,FCO2,FCO2_var)
   
end subroutine read_SOCAT_FCO2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine read_SOCAT_grid(fname,gr)
  use mod_grid
  !use m_datatest
  use m_spherdist
  use netcdf
  !use m_nf90_err
  use nfw_mod
  implicit none


  character(len=80),   intent(in) :: fname
  type(grid),        intent(out) :: gr
  integer::nmonth

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, GridDepth_ID, LatLon_ID, NbTime_id

! Variable ids
  integer :: vLatLon_ID, vNbLatitudes_ID, vNbLongitudes_ID, vLatLonMin_ID, vLatLonStep_ID, &
           vGrid0001_ID,vflagsid

! Array dimensions
  integer :: NbLatitudes, NbLongitudes, GridDepth, NbTime

  real,dimension(:),allocatable::lon,lat

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
  call nfw_inq_dimid(trim(fname), ncid, 'YLAT', NbLatitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'XLON', NbLongitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'TMNTH', Nbtime_ID)
  print*,'How far do you go'

! Get dimension length from id
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
  call nfw_inq_dimlen(trim(fname),ncid,NbLatitudes_ID,NbLatitudes)
  call nfw_inq_dimlen(trim(fname),ncid,NbLongitudes_ID,NbLongitudes)
  call nfw_inq_dimlen(trim(fname),ncid,Nbtime_ID,NbTime)
  print*, 'Dimensions:', NbLatitudes, NbLongitudes,NbTime
 
  call nfw_inq_varid(trim(fname),ncid,'YLAT' ,vNbLatitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'XLON' ,vNbLongitudes_ID)
  
  allocate(lon(NbLongitudes))
  allocate(lat(NbLatitudes))
  call nfw_get_var_double(trim(fname), ncid,vNbLongitudes_ID  ,lon)
  call nfw_get_var_double(trim(fname),ncid,vNbLatitudes_ID   ,lat)
 
  !call nf90_err (nf90_close(ncid))
  call nfw_close(trim(fname), ncid)  
  
  !grid
  gr%nx=NbLatitudes
  gr%ny=NbLongitudes
  gr%x0=lat(1)!0
  gr%y0=lon(1)!0
  gr%dx=1.
  gr%dy=1.
  gr%reg = .true.
  gr%order = 2
  gr%ux = 'deg'
  gr%uy = 'deg'
  gr%set = .true.

end subroutine read_SOCAT_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_SOCAT_FCO2_asynchro(fname,nmonth,gr,data,nwind)
  use mod_measurement
  use mod_grid
  !use m_datatest
  use m_spherdist
  use netcdf
  !use m_nf90_err
  use nfw_mod
  implicit none

  type (measurement),  intent(inout) :: data(:)
  character(len=80),   intent(in) :: fname
  type(grid),        intent(out) :: gr
  integer::nmonth,nwind

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, GridDepth_ID, LatLon_ID, NbTime_id

! Variable ids
  integer :: vLatLon_ID, vNbLatitudes_ID, vNbLongitudes_ID, vLatLonMin_ID, vLatLonStep_ID, &
           vGrid0001_ID,vflagsid,vlatoff_id,vlonoff_id

! Array dimensions
  integer :: NbLatitudes, NbLongitudes, GridDepth, NbTime

! Data arrays
  real,allocatable :: FCO2(:,:,:), FCO2_var(:,:,:), lon(:),lat(:)
  real, dimension(:,:,:),allocatable::lonoff,latoff 
! utilitary
  integer ncid, ijmax(2)
  real undef(1),undef_lat, undef_lon
  integer i, j,k,l
  logical valid
  real, parameter :: epsilon = 0.01  ! test for undefined values

  real, dimension(:,:),allocatable::obsday
  integer::vday_id
  
! Open file
  !call nf90_err(NF90_OPEN(trim(fname),NF90_NOCLOBBER,ncid))
  call nfw_open(trim(fname), nf_nowrite, ncid)
  
! Get dimension id in netcdf file ...
  !call nf90_err(nf90_Inq_Dimid(ncid,'lat',NbLatitudes_ID))
  !call nf90_err(nf90_Inq_Dimid(ncid,'lon',NbLongitudes_ID))
  call nfw_inq_dimid(trim(fname), ncid, 'YLAT', NbLatitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'XLON', NbLongitudes_ID)
  call nfw_inq_dimid(trim(fname), ncid, 'TMNTH', Nbtime_ID)
  print*,'How far do you go'

! Get dimension length from id
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  !call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
  call nfw_inq_dimlen(trim(fname),ncid,NbLatitudes_ID,NbLatitudes)
  call nfw_inq_dimlen(trim(fname),ncid,NbLongitudes_ID,NbLongitudes)
  call nfw_inq_dimlen(trim(fname),ncid,Nbtime_ID,NbTime)
  print*, 'Dimensions:', NbLatitudes, NbLongitudes,NbTime

! State which variable you want here.. Available vars are shown when you do
  allocate(lon(NbLongitudes))
  allocate(lat(NbLatitudes))
  !allocate(FCO2(NbTime,NbLongitudes,NbLatitudes))
  !allocate(FCO2_var(NbTime,NbLongitudes,NbLatitudes))
  !allocate(FCO2(NbTime,NbLatitudes,NbLongitudes))
  !allocate(FCO2_var(NbTime,NbLatitudes,NbLongitudes))
  allocate(FCO2(NbLongitudes,NbLatitudes, NbTime))
  allocate(FCO2_var(NbLongitudes,NbLatitudes, NbTime))
  allocate(lonoff(NbLongitudes,NbLatitudes, NbTime))
  allocate(latoff(NbLongitudes,NbLatitudes, NbTime))
  
  
! Variable ids in netcdf file
  !call nf90_err(nf90_inq_varid(ncid,'lat' ,vNbLatitudes_ID),'lat')
  !call nf90_err(nf90_inq_varid(ncid,'lon' ,vNbLongitudes_ID),'lon')
  !call nf90_err(nf90_inq_varid(ncid,'CHL1_mean' ,vGrid0001_ID),'CHL1_mean')
  !call nf90_err(nf90_inq_varid(ncid,'CHL1_flags' ,vflagsid),'CHL1_flags')
  call nfw_inq_varid(trim(fname),ncid,'YLAT' ,vNbLatitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'XLON' ,vNbLongitudes_ID)
  call nfw_inq_varid(trim(fname),ncid,'FCO2_AVE_UNWTD' ,vGrid0001_ID)
  call nfw_inq_varid(trim(fname),ncid,'FCO2_STD_UNWTD' ,vflagsid)
  call nfw_inq_varid(trim(fname),ncid,'LON_OFFSET_UNWTD' ,vlonoff_id)
  call nfw_inq_varid(trim(fname),ncid,'LAT_OFFSET_UNWTD' ,vlatoff_id)

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
  call nfw_get_var_double(trim(fname),ncid,vGrid0001_ID,FCO2)
  print *,'Range FCO2 ', minval(FCO2), maxval(FCO2)
  !call nf90_err(nf90_get_var(ncid,vflagsid,chla_var))
  call nfw_get_var_double(trim(fname),ncid,vflagsid,FCO2_var)
  print *,'Range FCO2_var ', minval(FCO2_var), maxval(FCO2_var)

  call nfw_get_var_double(trim(fname),ncid,vlatoff_id,latoff)
  print *,'Range Lat_offset ', minval(latoff), maxval(latoff)
  call nfw_get_var_double(trim(fname),ncid,vlonoff_id,lonoff)
  print *,'Range Lon_offset ', minval(lonoff), maxval(lonoff)
  
  
  
  allocate(obsday(2,NbTime))
  call nfw_inq_varid(trim(fname),ncid,'TMNTH_bnds' ,vday_id)
  call nfw_get_var_double(trim(fname), ncid,vday_id  ,obsday)
  
  !call nf90_err (nf90_close(ncid))
  call nfw_close(trim(fname), ncid)  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fill the data(:) vector
  k=0
  do l=0,nwind-1
  do j=1,NbLongitudes ! gr%ny
  do i=1,NbLatitudes ! gr%nx
    
    k=k+1
    data(k)%id = 'AFCO2'
    data(k)%d = FCO2(j,i,nmonth+l)
    !data(k)%d = chla(i,j)

     data(k)%ipiv = i
     data(k)%jpiv = j

     data(k)%lat=lat(i)+latoff(j,i,nmonth+l)
     data(k)%lon=ang180(lon(j)+lonoff(j,i,nmonth+l))

!LB: Data support is assumed = a square grid cell
!support diameter in meters stored in %a1 (tricky, isn't it ?)
!     data(k)%a1 = spherdist(lon(j)-.5*gr%dx,lat(i)-.5*gr%dy, &
!                            lon(j)+.5*gr%dx,lat(i)+.5*gr%dy)

     data(k)%a1 = 1
     data(k)%a2 = 0
     data(k)%a3 = 0
     data(k)%a4 = 0
     
     data(k)%ns = 0

     data(k)%date = l+1
     data(k)%var = FCO2_var(j,i,nmonth+l)

     data(k)%depth=0.0

      valid =  abs( (FCO2(j,i,nmonth+l)-undef(1))   / undef(1) )  > epsilon  
      !valid =  abs( (chla(i,j)-undef(1))   / undef(1) )  > epsilon  

     data(k)%status = valid
  enddo
  enddo
  print*,'day',l,obsday(1,nmonth+l)
  enddo
  print*, 'Number of data read:', k, gridpoints(gr)*nwind

!CONSISTENCY...
  !call datatest(gr,data,k,'consistency_d_'//trim(data(1)%id))

  deallocate(lat,lon,FCO2,FCO2_var,lonoff,latoff)
  deallocate(obsday)
   
  end subroutine read_SOCAT_FCO2_asynchro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111  
  
  subroutine nc_obs2d(obs,nobs,nx,ny,depths)
    use mod_measurement
    use netcdf
    use nfw_mod
    implicit none
    integer::nobs,nx,ny
    type (measurement)::obs(nobs)
    real,dimension(nx,ny)::depths
    !real,dimension(nx,ny)::obs2d
    real,dimension(:,:,:),allocatable::obs2d
    integer::k,ierr,ncid,dimnx,dimny,var_id,var_id2,i,j,dimnz,cpt
    character(len=80) :: fname   
    integer::a0
    character(len=8)::cfld
    
    fname='Observations_2D.nc'
    
    cpt=1
    a0=obs(1)%date
 
    do k=1,nobs
      if(obs(k)%date.gt.a0)then
        cpt=cpt+1
	a0=obs(k)%date
      endif
    enddo

    allocate(obs2d(nx,ny,cpt))
    obs2d(:,:,:)=-2.
    do k=1,nobs
      obs2d(obs(k)%ipiv, obs(k)%jpiv,obs(k)%date)=obs(k)%d  
    enddo
    do j=1,ny
      do i=1,nx
        if(depths(i,j).le.0.)then
	  obs2d(i,j,:)=-2.
	endif
      enddo
    enddo
    
    call nfw_create(trim(fname), nf_write, ncid)
        
    call nfw_def_dim(trim(fname), ncid, 'nx', nx, dimnx)
    call nfw_def_dim(trim(fname), ncid, 'ny', ny, dimny)
    call nfw_def_dim(trim(fname), ncid, 'nmonth', cpt, dimnz)
    !ierr=NF90_DEF_DIM(ncid,'nx',nx,dimnx)
    !ierr=NF90_DEF_DIM(ncid,'ny',ny,dimny)

    !ierr=NF90_REDEF(ncid)
    !ierr=NF90_DEF_VAR(ncid,'CHLA',NF90_Float,(/dimnx,dimny/),var_id)
    !ierr=NF90_ENDDEF(ncid)
    !ierr=NF90_PUT_VAR(ncid,var_id,obs2d(1:nx,1:ny))
    
    cfld=trim(obs(1)%id)
    call nfw_def_var(trim(fname), ncid, trim(cfld), nf_float,3,(/dimnx,dimny,dimnz/), var_id)
    call nfw_def_var(trim(fname), ncid, 'depths', nf_float,2,(/dimnx,dimny/), var_id2)
    call nfw_enddef(trim(fname), ncid)

    call nfw_put_var_double(trim(fname), ncid, var_id,obs2d(1:nx,1:ny,1:cpt))
    call nfw_put_var_double(trim(fname), ncid, var_id2,depths(1:nx,1:ny))
    !ierr=NF90_REDEF(ncid)
    !ierr=NF90_DEF_VAR(ncid,'depths',NF90_Float,(/dimnx,dimny/),var_id)
    !ierr=NF90_ENDDEF(ncid)
    !ierr=NF90_PUT_VAR(ncid,var_id,depths(1:nx,1:ny))

    !ierr=NF90_CLOSE(ncid)
    call nfw_close(trim(fname), ncid)
    deallocate(obs2d)
    
    
 end subroutine nc_obs2d
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer function get_nobs_SOCAT_TW(nx,ny,nwnd)
   use mod_measurement
   use nfw_mod
   implicit none  
   integer :: nobs,nwnd,nx,ny
   character(80) :: fname,memfile
   integer :: ncid, var_id,k,j
   type (measurement), allocatable :: data(:)
   integer,dimension(:),allocatable::itemp 
   real::att(2)
   real,dimension(nx,ny)::fld
   character(3)::cmonth
 
   nobs=-1
   memfile='observations-AFCO2.nc'
   call nfw_open(memfile, nf_nowrite, ncid)
    
   call nfw_inq_dimid(memfile, ncid, 'nobs', var_id)
   call nfw_inq_dimlen(memfile,ncid,var_id,nobs)
   
   allocate(data(nobs))
   allocate(itemp(nobs))
   
   call nfw_inq_varid(memfile,ncid,'ipiv' ,var_id)
   call nfw_get_var_int(memfile, ncid, var_id, itemp)
   data(1:nobs) % ipiv=itemp(1:nobs)
  
   call nfw_inq_varid(memfile,ncid,'jpiv' ,var_id)
   call nfw_get_var_int(memfile, ncid, var_id, itemp)
   data(1:nobs) % jpiv=itemp(1:nobs)
  
   call nfw_inq_varid(memfile,ncid,'age' ,var_id)
   call nfw_get_var_int(memfile, ncid, var_id, itemp)
   data(1:nobs) % date=itemp(1:nobs)
   call nfw_close(memfile, ncid)
   
   get_nobs_SOCAT_TW=0
   do k=1,nwnd
      write(cmonth,'(i3.3)')k
      memfile='monthly_ref_T'//cmonth
      
      call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
    
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'fice',var_id)
      call nfw_get_var_double(trim(memfile)//'.nc', ncid, var_id, fld)
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, var_id, 'add_offset', att(2))
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, var_id, 'scale_factor', att(1))
     
      call nfw_close(trim(memfile)//'.nc', ncid)
      
      
      do j=1,nobs
	if((data(j)%date==k).and.(fld(data(j)%ipiv,data(j)%jpiv)*att(1)+att(2)==0.))then
	  get_nobs_SOCAT_TW=get_nobs_SOCAT_TW+1
	endif  
      enddo
       
   enddo

   deallocate(data,itemp)
   print*,'FCO2 in ice: ', 1.-real(get_nobs_SOCAT_TW)/real(nobs)

   return
   
 end function
 
 subroutine get_obs_SOCAT_TW(obs,nobs,nx,ny,nwnd)
   use mod_measurement
   use nfw_mod
   implicit none  
   integer :: nobs,nwnd,nx,ny
   type(measurement):: obs(nobs)
   character(80) :: fname,memfile
   integer :: ncid, var_id,ntobs
   real,dimension(:),allocatable::rtemp 
   integer,dimension(:),allocatable::itemp 
   type (measurement), allocatable :: data(:)
   real,dimension(nx,ny)::fld
   character(3)::cmonth
   real::att(2)
   integer::j,k,cpt
   
   fname='observations-AFCO2.nc'
   call nfw_open(fname, nf_nowrite, ncid)
   call nfw_inq_dimid(memfile, ncid, 'nobs', var_id)
   call nfw_inq_dimlen(memfile,ncid,var_id,ntobs)
   
   allocate(data(ntobs))
   allocate(rtemp(ntobs))
   allocate(itemp(ntobs))
   
   call nfw_inq_varid(fname,ncid,'lat' ,var_id)
   call nfw_get_var_double(fname, ncid, var_id, rtemp)
   data(1:ntobs) % lat=rtemp(1:ntobs)
   
   call nfw_inq_varid(fname,ncid,'lon' ,var_id)
   call nfw_get_var_double(fname, ncid, var_id, rtemp)
   data(1:ntobs) % lon=rtemp(1:ntobs)
   
   call nfw_inq_varid(fname,ncid,'depth' ,var_id)
   call nfw_get_var_double(fname, ncid, var_id, rtemp)
   data(1:ntobs) % depth=rtemp(1:ntobs)
   
   call nfw_inq_varid(fname,ncid,'d' ,var_id)
   call nfw_get_var_double(fname, ncid, var_id, rtemp)
   data(1:ntobs) % d=rtemp(1:ntobs)
   
   call nfw_inq_varid(fname,ncid,'var' ,var_id)
   call nfw_get_var_double(fname, ncid, var_id, rtemp)
   data(1:ntobs) % var=rtemp(1:ntobs)
   
   call nfw_inq_varid(fname,ncid,'age' ,var_id)
   call nfw_get_var_int(fname, ncid, var_id, itemp)
   data(1:ntobs) % date=itemp(1:ntobs)
  
   call nfw_inq_varid(fname,ncid,'ipiv' ,var_id)
   call nfw_get_var_int(fname, ncid, var_id, itemp)
   data(1:ntobs) % ipiv=itemp(1:ntobs)
  
   call nfw_inq_varid(fname,ncid,'jpiv' ,var_id)
   call nfw_get_var_int(fname, ncid, var_id, itemp)
   data(1:ntobs) % jpiv=itemp(1:ntobs)
   
   call nfw_close(fname, ncid)
   
   deallocate(rtemp,itemp)
   
   cpt=1
   do k=1,nwnd
      write(cmonth,'(i3.3)')k
      memfile='monthly_ref_T'//cmonth
      
      call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
    
      call nfw_inq_varid(trim(memfile)//'.nc', ncid,'fice',var_id)
      call nfw_get_var_double(trim(memfile)//'.nc', ncid, var_id, fld)
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, var_id, 'add_offset', att(2))
      call nfw_get_att_double(trim(memfile)//'.nc', ncid, var_id, 'scale_factor', att(1))
     
      call nfw_close(trim(memfile)//'.nc', ncid)
      
      
      do j=1,ntobs
	if((data(j)%date==k).and.(fld(data(j)%ipiv,data(j)%jpiv)*att(1)+att(2)==0.))then
	  obs(cpt)%lat=data(j) % lat
	  obs(cpt)%lon=data(j) % lon
	  obs(cpt)%depth=data(j) % depth
	  obs(cpt)%d=data(j) % d
	  obs(cpt) % var=data(j) % var
	  obs(cpt) %date=data(j) % date
	  obs(cpt) %ipiv=data(j) % ipiv
	  obs(cpt) %jpiv=data(j) % jpiv
	  cpt=cpt+1
	endif  
      enddo
       
   enddo
   
   deallocate(data)
   
 end subroutine
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 
   subroutine compute_noise_SOCAT(noise,obs,nobs)
   use mod_measurement
   use m_set_random_seed2
   use m_random
   implicit none
   integer :: nobs
   type(measurement):: obs(nobs)
   real, dimension(1:nobs)::noise
   integer::i 
   real::var 
    
   call set_random_seed3
   call randn(nobs,noise)
   
   print*,'noise',noise(5:10)
   var=0.
   do i=1,nobs
      var=var+obs(i)%var
   enddo
   var=var/real(nobs)
   noise(:)=noise(:)*var
   obs(:)%var=var
   
   end subroutine
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11  
   
   subroutine read_hamocc_AFCO2(fname,obs,nx,ny,nobs,noise,nrobs)
   use mod_measurement
   use nfw_mod
   implicit none
   integer::nx,ny,nwnd,nobs,nrobs
   type(measurement):: obs(nobs)
   real,dimension(nobs)::noise
   character(80) :: fname,memfile
   integer::d0,k,vFIELD_ID,ncid
   character(3)::cmonth
   real,dimension(nx,ny)::fld
   real::att(2)
   
   d0=obs(1)%date
   
   do k=1,nobs
      write(cmonth,'(i3.3)')obs(k)%date
      memfile=trim(fname)//cmonth

      if((obs(k)%date.gt.d0).or.(k==1))then
        call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
        call nfw_inq_varid(trim(memfile)//'.nc', ncid,'pco2',vFIELD_ID)
	call nfw_get_var_double(trim(memfile)//'.nc', ncid, vFIELD_ID, fld)
        call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'add_offset', att(2))
        call nfw_get_att_double(trim(memfile)//'.nc', ncid, vFIELD_ID, 'scale_factor', att(1))
	call nfw_close(trim(memfile)//'.nc', ncid)
	
	!call get_micom_fld_eco(trim(memfile),fld,iens,'pco2',1,1,nx,ny)
	d0=obs(k)%date
      endif    
      
      obs(k)%a1 = 1.
      obs(k)%a2 = 0.
      obs(k)%a3 = 0.
      obs(k)%a4 = 0.
     
      obs(k)%ns = 1!0

      obs(k)%id ='AFCO2'
      obs(k)%h = 1
      obs(k)%d =fld(obs(k)%ipiv,obs(k)%jpiv)*att(1)+att(2)+noise(k)
      obs(k)%status=.true.
      !if(obs(k)%date.gt.1)obs(k)%status=.false.
      obs(k)%i_orig_grid = -1
      obs(k)%j_orig_grid = -1
   enddo
   nrobs=nobs
   
   end subroutine	
!!!!!!!!!!!!!!!!!!!!!!!!! 
 
end module m_read_SOCAT_gridded
