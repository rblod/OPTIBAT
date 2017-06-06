module m_read_OCCCI
! Ehouarn !
! Reads CHLA1 data and grid from the GlobColour files !
contains

  subroutine read_OCCCI_grid(fname,gr)
  use mod_measurement
  use mod_grid
  !use m_datatest
  use m_spherdist
  use netcdf
  !use m_nf90_err
  use nfw_mod
  implicit none


  type (grid),         intent(inout) :: gr ! CLS measurement grid
  character(len=80),   intent(in) :: fname

! Array dimensions
  integer :: blist_id, bin_list
  
! utilitary
  integer ncid
  logical :: ex
  
  print *, 'read_OC-CCI_grid'
  
  gr = default_grid
  gr%nx=0
  
  inquire(file=trim(fname),exist=ex)
  if(ex) then
      call nfw_open(fname, nf_nowrite, ncid)
      print *, '  found "', trim(fname), '"...'

      ! Get dimension id in netcdf file ...
      call nfw_inq_dimid(fname, ncid, 'bin_list',blist_id)
     
      ! Get dimension length from id
      call nfw_inq_dimlen(fname, ncid, blist_id, bin_list)
      call nfw_close(fname, ncid)

      gr%nx=bin_list
      print*,'gr%nx',gr%nx
      gr%ny=1
      gr%x0=0
      gr%y0=0
      gr%dx=0.1
      gr%dy=0.1
      gr%reg = .false.
      gr%order = 1
      gr%ux = 'm'
      gr%uy = 'm'
      gr%set = .true.
  endif
   
  end subroutine read_OCCCI_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_OCCCI_CHLA_v0(fname,gr,data)
    use mod_measurement
    use mod_grid
    use nfw_mod
    implicit none

    character(len=80), intent(in) :: fname
    type(grid), intent(inout) :: gr ! CLS measurement grid
    type(measurement), intent(inout) :: data(:)

    integer :: bindex_ID, blist_ID
    integer :: vNbpoint_ID, vLongitude_ID, vLatitude_ID, vBegindate_ID, vchla_ID, &
               vchla_rms_id,vbextent_id,vbstartnum_ID

    ! array dimensions
    integer :: bin_index,bin_list

    ! data arrays
    real(8), allocatable :: vchla(:), vlat(:),vchla_rms(:),vchla_bias(:)
    integer, allocatable ::bi_extent(:),vstart(:)
    logical, allocatable :: isgood(:)

    integer :: ncid
    real(8), dimension(1) :: undef_chla, undef_chla_rms, undef_chla_bias!, undef_begindate
    real(8) :: varsat
    integer, dimension(1) :: undef_nbpoint
    integer :: i, j, k, nobs, obsid, sid, age
    real(8), parameter :: EPS = 0.01  ! test for undefined values
    logical :: ex
    real::lat,lon 
     
     
    print *, 'read_OC-CCI_CHLA'

    nobs = 0
  

    inquire(file = trim(fname), exist = ex)
    if (.not. ex) then
        stop
    end if

    ! Reading the observation file of satellite
    call nfw_open(fname, nf_nowrite, ncid)  
    call nfw_inq_dimid(fname, ncid, 'bin_index', bindex_ID)
    call nfw_inq_dimid(fname, ncid, 'bin_list', blist_ID)
          
    call nfw_inq_dimlen(fname, ncid, bindex_ID, bin_index)
    call nfw_inq_dimlen(fname, ncid, blist_ID, bin_list)
    print '(1x, a, 3i8)', '    dimensions (# bin_index, # bin_list):',bin_index,bin_list
       
    allocate(vlat(bin_index))
    allocate(vstart(bin_index),bi_extent(bin_index))
    allocate(vchla(bin_list),vchla_rms(bin_list))!,vchla_bias(bin_list))   
    
    !allocate(vnbpoint(ntracks), vbegindate(ncycles, ntracks))
   ! allocate(isgood(bin_list))
           
    ! Variable ids in netcdf file
    !call nfw_inq_varid(fname, ncid, 'bi_vsize', vLatitude_ID)
    call nfw_inq_varid(fname, ncid,'bi_hsize', vLongitude_ID)
    !call nfw_inq_varid(fname, ncid,'BeginDates', vBegindate_ID)
    !call nfw_inq_varid(fname, ncid,'NbPoints', vNbpoint_ID)
    call nfw_inq_varid(fname, ncid,'chlor_a', vchla_id)
    call nfw_inq_varid(fname, ncid,'chlor_a_rms_uncertainty', vchla_rms_id)
    !call nfw_inq_varid(fname, ncid,'chlor_a_bias_uncertainty', vchla_bias_id)
    call nfw_inq_varid(fname, ncid,'bi_extent', vbextent_ID)
    call nfw_inq_varid(fname, ncid,'bi_start_num', vbstartnum_ID)
       
    ! Variable _FillValue attributes
    call nfw_get_att_double(fname, ncid, vchla_ID, '_FillValue', undef_chla(1))
    call nfw_get_att_double(fname, ncid, vchla_rms_ID, '_FillValue', undef_chla_rms(1))
    !call nfw_get_att_double(fname, ncid,vchla_bias_ID, '_FillValue', undef_chla_bias(1))
    gr % undef = undef_chla(1)
          
    !call nfw_get_var_double(fname, ncid, vLongitude_ID, vlon)
    call nfw_get_var_double(fname, ncid, vLatitude_ID, vlat)
    call nfw_get_var_double(fname, ncid, vchla_ID, vchla)
    call nfw_get_var_double(fname, ncid, vchla_rms_ID, vchla_rms)
    !call nfw_get_var_double(fname, ncid, vchla_ID, vchla_bias)
    call nfw_get_var_int(fname, ncid, vbstartnum_ID, vstart)
    
    !lon = ang180(lon)
    !vlon = vlon * 1.e-06
    !vlat = vlat * 1.e-06
    !print '(1x, a, 2f10.2)', '    range Lon = ', minval(vlon), maxval(vlon)
    !print '(1x, a, 2f10.2)', '    range Lat = ', minval(vlat), maxval(vlat)
    print '(1x, a, 2f10.2)', '    range CHLA = ', minval(vchla), maxval(vchla)
    print '(1x, a, 2f10.2)', '    range CHLA RMS= ', minval(vchla_rms), maxval(vchla_rms)
    !print '(1x, a, 2f10.2)', '    range CHLA BIAS= ', minval(vchla_bias), maxval(vchla_bias)
          
    call nfw_get_var_int(fname, ncid,vbextent_ID,bi_extent)
    
    call nfw_close(fname, ncid)
    
    nobs=0
    do i=1,bin_list!index
      j=get_row_index_v0(bin_index,bi_extent,vstart,i)
         
      if ((vchla(i) == undef_chla(1)).or.(vchla_rms(i) == undef_chla_rms(1))) then
        cycle
      end if
      nobs = nobs + 1	
      data(nobs) % id = 'CHLA'
      data(nobs) % d = vchla(i)              
      data(nobs) % ipiv = -1 ! to be filled
      data(nobs) % jpiv = -1
      lat=(real(j)+0.5)*180./real(bin_index)-90.
      data(nobs) % lat = lat
      lon=(real(i-(vstart(j)-1))+0.5)*vlat(i)-180.
      data(nobs) % lon = lon
      data(nobs) % a1 = -1.0e10 ! to be filled
      data(nobs) % a2 = -1.0e10
      data(nobs) % a3 = -1.0e10
      data(nobs) % a4 = -1.0e10
      data(nobs) % ns = 0
      data(nobs) % var = vchla_rms(i)
      data(nobs) % date = 0
      data(nobs) % depth = 0.0
      data(nobs) % status = .true.
   enddo


   print*, '    # of obs read so far = ', nobs,bin_list
   deallocate(vlat, vchla,vchla_rms,bi_extent,vstart)

    gr % nx = nobs
  end subroutine read_OCCCI_CHLA_v0
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function get_row_index_v0(bin_index,number_of_bins,start_index,bin_target)
   implicit none
   integer::bin_index,bin_target
   integer,dimension(bin_index)::number_of_bins,start_index
   integer::m
   
   get_row_index_v0=0
   do m=1,bin_index
      if((number_of_bins(m)+start_index(m)).gt.bin_target)then
        exit 
      else
        get_row_index_v0=get_row_index_v0+1 
      endif
   enddo
   get_row_index_v0=get_row_index_v0-1
   return
   
   end function
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 
   integer function get_row_index(bin_index,start_index,bin_target)
   implicit none
   integer::bin_index,bin_target
   integer,dimension(bin_index)::start_index
   integer::m
   
   get_row_index=0
   do while (get_row_index .lt. bin_index-1)
      if(start_index(get_row_index+1).gt.bin_target)then
	return
      else
        get_row_index=get_row_index+1 
      endif
   enddo
   return
   
   end function
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!   
 
  subroutine read_OCCCI_CHLA_binned(fname,gr,data)
    use mod_measurement
    use mod_grid
    use nfw_mod
    implicit none

    character(len=80), intent(in) :: fname
    type(grid), intent(inout) :: gr ! CLS measurement grid
    type(measurement), intent(inout) :: data(:)

    integer :: bindex_ID, blist_ID
    integer :: vNbpoint_ID, vLongitude_ID, vLatitude_ID, vBegindate_ID, vchla_ID, &
               vchla_rms_id,vbextent_id,vbstartnum_ID,vbimax_id,vbegin_id,vbl_id

    ! array dimensions
    integer :: bin_index,bin_list

    ! data arrays
    real(8), allocatable :: vchla(:), vlat(:),vchla_rms(:),vchla_bias(:)
    integer, allocatable ::bi_extent(:),vstart(:),vbimax(:),vbegin(:),vbl_bin_num(:)
    logical, allocatable :: isgood(:)

    integer :: ncid
    real(8), dimension(1) :: undef_chla, undef_chla_rms, undef_chla_bias!, undef_begindate
    real(8) :: varsat
    integer, dimension(1) :: undef_nbpoint
    integer :: i, j, k, nobs, obsid, sid, age,j2
    real(8), parameter :: EPS = 0.01  ! test for undefined values
    logical :: ex
    real::lat,lon 
     
     
    print *, 'read_OC-CCI_CHLA'

    nobs = 0
  

    inquire(file = trim(fname), exist = ex)
    if (.not. ex) then
        stop
    end if

    ! Reading the observation file of satellite
    call nfw_open(fname, nf_nowrite, ncid)  
    call nfw_inq_dimid(fname, ncid, 'bin_index', bindex_ID)
    call nfw_inq_dimid(fname, ncid, 'bin_list', blist_ID)
          
    call nfw_inq_dimlen(fname, ncid, bindex_ID, bin_index)
    call nfw_inq_dimlen(fname, ncid, blist_ID, bin_list)
    print '(1x, a, 3i8)', '    dimensions (# bin_index, # bin_list):',bin_index,bin_list
       
    !allocate(vlat(bin_index))
    allocate(vstart(bin_index),bi_extent(bin_index),vbimax(bin_index),vbegin(bin_index))
    allocate(vchla(bin_list),vchla_rms(bin_list))!,vchla_bias(bin_list))   
    allocate(vbl_bin_num(bin_list))
 
           
    ! Variable ids in netcdf file
   
    call nfw_inq_varid(fname, ncid,'bi_max', vbimax_ID)
    call nfw_inq_varid(fname, ncid,'chlor_a', vchla_id)
    call nfw_inq_varid(fname, ncid,'chlor_a_rms_uncertainty', vchla_rms_id)
    call nfw_inq_varid(fname, ncid,'bi_extent', vbextent_ID)
    !call nfw_inq_varid(fname, ncid,'bi_start_num', vbstartnum_ID)
    call nfw_inq_varid(fname, ncid,'bi_start_num', vbstartnum_ID)
    call nfw_inq_varid(fname, ncid,'bi_begin', vbegin_ID)
    call nfw_inq_varid(fname, ncid,'bl_bin_num', vbl_ID)
       
    ! Variable _FillValue attributes
    call nfw_get_att_double(fname, ncid, vchla_ID, '_FillValue', undef_chla(1))
    call nfw_get_att_double(fname, ncid, vchla_rms_ID, '_FillValue', undef_chla_rms(1))
    gr % undef = undef_chla(1)
          
    !call nfw_get_var_double(fname, ncid, vLongitude_ID, vlon)
    !call nfw_get_var_double(fname, ncid, vLatitude_ID, vlat)
    call nfw_get_var_double(fname, ncid, vchla_ID, vchla)
    call nfw_get_var_double(fname, ncid, vchla_rms_ID, vchla_rms)
    !call nfw_get_var_double(fname, ncid, vchla_ID, vchla_bias)
    call nfw_get_var_int(fname, ncid, vbstartnum_ID, vstart)
    call nfw_get_var_int(fname, ncid, vbimax_ID, vbimax)
    call nfw_get_var_int(fname, ncid, vbegin_ID, vbegin)
    call nfw_get_var_int(fname, ncid, vbl_ID, vbl_bin_num)
    
    !lon = ang180(lon)
    !vlon = vlon * 1.e-06
    !vlat = vlat * 1.e-06
    !print '(1x, a, 2f10.2)', '    range Lon = ', minval(vlon), maxval(vlon)
    !print '(1x, a, 2f10.2)', '    range Lat = ', minval(vlat), maxval(vlat)
    print '(1x, a, 2f10.2)', '    range CHLA = ', minval(vchla), maxval(vchla)
    print '(1x, a, 2f10.2)', '    range CHLA RMS= ', minval(vchla_rms), maxval(vchla_rms)
    !print '(1x, a, 2f10.2)', '    range CHLA BIAS= ', minval(vchla_bias), maxval(vchla_bias)
          
    call nfw_get_var_int(fname, ncid,vbextent_ID,bi_extent)
    
    call nfw_close(fname, ncid)
    
    nobs=0
    do i=1,bin_list!index
      
      j=get_row_index(bin_index,vbegin,vbl_bin_num(i))
    
      if ((vchla(i) == undef_chla(1)).or.(vchla_rms(i) == undef_chla_rms(1))) then
        cycle
      end if
      nobs = nobs + 1	
      data(nobs) % id = 'CHLA'
      data(nobs) % d = vchla(i)              
      data(nobs) % ipiv = -1 ! to be filled
      data(nobs) % jpiv = -1
      lat=(real(j)+0.5)*180./real(bin_index)-90.
      data(nobs) % lat = lat
      !print*,'i,j, ',i,j
      !lon=(real(j)+0.5)*360./real(vbimax(i))-180.
      !lon=(real(vbl_bin_num(i))-real(vstart(j))+0.5)*360./real(vbimax(j))-180.
      lon=(real(vbl_bin_num(i))-real(vstart(j))+0.5)*360./real(vbimax(j))-180.
      !lon=(real(j2)-real(vstart(j))+0.5)*360./real(vbimax(j))-180.
      data(nobs) % lon = ang180(lon)
      data(nobs) % a1 = -1.0e10 ! to be filled
      data(nobs) % a2 = -1.0e10
      data(nobs) % a3 = -1.0e10
      data(nobs) % a4 = -1.0e10
      data(nobs) % ns = 0
      data(nobs) % var = vchla_rms(i)
      !print*,' obs,var ', data(nobs) % d, data(nobs) % var
      !print*,' lat,lon ', data(nobs) % lat, data(nobs) % lon
      data(nobs) % date = 0
      data(nobs) % depth = 0.0
      data(nobs) % status = .true.
      print*,'lat,lon',lat,ang180(lon)
   
   enddo


   print*, '    # of obs read so far = ', nobs,bin_list
   deallocate(vbimax, vchla,vchla_rms,bi_extent,vstart,vbl_bin_num)

    gr % nx = nobs
  end subroutine read_OCCCI_CHLA_binned  
  
!!!!!!!!!!!!!!!!1 

  subroutine read_OCCCI_singrid(fname,gr)
  use mod_measurement
  use mod_grid
  !use m_datatest
  use m_spherdist
  use netcdf
  !use m_nf90_err
  use nfw_mod
  implicit none


  type (grid),         intent(inout) :: gr ! CLS measurement grid
  character(len=80),   intent(in) :: fname

! Array dimensions
  integer :: sing_id, singrid
  
! utilitary
  integer ncid
  logical :: ex
  
  print *, 'read_OC-CCI_grid'
  
  gr = default_grid
  gr%nx=0
  
  inquire(file=trim(fname),exist=ex)
  if(ex) then
      call nfw_open(fname, nf_nowrite, ncid)
      print *, '  found "', trim(fname), '"...'

      ! Get dimension id in netcdf file ...
      call nfw_inq_dimid(fname, ncid, 'sin_grid',sing_id)
     
      ! Get dimension length from id
      call nfw_inq_dimlen(fname, ncid, sing_id, singrid)
      call nfw_close(fname, ncid)

      gr%nx=singrid
      print*,'gr%nx',gr%nx
      gr%ny=1
      gr%x0=0
      gr%y0=0
      gr%dx=0.1
      gr%dy=0.1
      gr%reg = .false.
      gr%order = 1
      gr%ux = 'm'
      gr%uy = 'm'
      gr%set = .true.
  endif
   
  end subroutine read_OCCCI_singrid

 !!!!!!!!!!!!!!!!!!!!!!!!!   
 
  subroutine read_OCCCI_CHLA_sin(fname,gr,data)
    use mod_measurement
    use mod_grid
    use nfw_mod
    implicit none

    character(len=80), intent(in) :: fname
    type(grid), intent(inout) :: gr ! CLS measurement grid
    type(measurement), intent(inout) :: data(:)

    integer :: bindex_ID, blist_ID
    integer :: vNbpoint_ID, vLongitude_ID, vLatitude_ID, vBegindate_ID, vchla_ID, &
               vchla_rms_id,vbextent_id,vbstartnum_ID,vbimax_id,vbegin_id,vbl_id

    ! array dimensions
    integer :: bin_index,bin_list

    ! data arrays
    real(8), allocatable :: vchla(:), vlat(:),vchla_rms(:),vlon(:)

    integer :: ncid
    real(8), dimension(1) :: undef_chla, undef_chla_rms, undef_chla_bias!, undef_begindate
    real(8) :: varsat
    integer, dimension(1) :: undef_nbpoint
    integer :: i, j, k, nobs, obsid, sid, age,sin_grid,vlon_id,vlat_id,sgrid_id
    real(8), parameter :: EPS = 0.01  ! test for undefined values
    logical :: ex
    real::lat,lon 
     
     
    print *, 'read_OC-CCI_CHLA'

    nobs = 0
  

    inquire(file = trim(fname), exist = ex)
    if (.not. ex) then
        stop
    end if

    ! Reading the observation file of satellite
    call nfw_open(fname, nf_nowrite, ncid)  
    call nfw_inq_dimid(fname, ncid, 'sin_grid', sgrid_ID)

          
    call nfw_inq_dimlen(fname, ncid, sgrid_ID, sin_grid)
    print '(1x, a, 3i8)', '    dimensions (# bin_index, # bin_list):',sin_grid
       
    allocate(vchla(sin_grid),vchla_rms(sin_grid),vlat(sin_grid),vlon(sin_grid))  

 
           
    ! Variable ids in netcdf file
   
    call nfw_inq_varid(fname, ncid,'chlor_a', vchla_id)
   ! call nfw_inq_varid(fname, ncid,'chlor_a_rms_uncertainty', vchla_rms_id)
    call nfw_inq_varid(fname, ncid,'latitude', vlat_ID)
    call nfw_inq_varid(fname, ncid,'longitude', vlon_ID)
       
    ! Variable _FillValue attributes
    call nfw_get_att_double(fname, ncid, vchla_ID, '_FillValue', undef_chla(1))
   ! call nfw_get_att_double(fname, ncid, vchla_rms_ID, '_FillValue', undef_chla_rms(1))
    gr % undef = undef_chla(1)
          
    call nfw_get_var_double(fname, ncid, vlon_ID, vlon)
    call nfw_get_var_double(fname, ncid, vlat_ID, vlat)
    call nfw_get_var_double(fname, ncid, vchla_ID, vchla)
   ! call nfw_get_var_double(fname, ncid, vchla_rms_ID, vchla_rms)
   
    call nfw_close(fname, ncid)

    print '(1x, a, 2f10.2)', '    range CHLA = ', minval(vchla), maxval(vchla)
    print '(1x, a, 2f10.2)', '    range CHLA RMS= ', minval(vchla_rms), maxval(vchla_rms)
          
    nobs=0
    do i=1,sin_grid

    !  if ((vchla(i) == undef_chla(1)).or.(vchla_rms(i) == undef_chla_rms(1))) then
      if ((vchla(i) == undef_chla(1))) then
        cycle
      end if
      nobs = nobs + 1	
      data(nobs) % id = 'CHLA'
      data(nobs) % d = vchla(i)              
      data(nobs) % ipiv = -1 ! to be filled
      data(nobs) % jpiv = -1
      data(nobs) % lat = vlat(i)
      data(nobs) % lon = ang180(vlon(i))
      data(nobs) % a1 = -1.0e10 ! to be filled
      data(nobs) % a2 = -1.0e10
      data(nobs) % a3 = -1.0e10
      data(nobs) % a4 = -1.0e10
      data(nobs) % ns = 0
    !  data(nobs) % var = vchla_rms(i)
      !print*,' obs,var ', data(nobs) % d, data(nobs) % var
      !print*,' lat,lon ', data(nobs) % lat, data(nobs) % lon
      data(nobs) % date = 1!0
      data(nobs) % depth = 0.0
      data(nobs) % status = .true.

   
   enddo


   print*, '    # of obs read so far = ', nobs,sin_grid
   deallocate(vchla,vchla_rms,vlon,vlat)

    gr % nx = nobs
  end subroutine read_OCCCI_CHLA_sin  
  
!!!!!!!!!!!!!!!!1 

end module m_read_OCCCI
