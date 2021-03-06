module m_read_metno_icec

contains

  subroutine read_metno_icec(fname, data, gr)
    use nfw_mod
    use mod_measurement
    use mod_grid
    implicit none

    character(*), intent(in) :: fname
    type (measurement), allocatable, intent(out) :: data(:)
    type(grid), intent(out) :: gr

    logical :: ex
    integer :: ncid
    integer :: xc_id, yc_id
    integer :: nx, ny
    integer :: lon_id, lat_id, icec_id, std_id, flag_id
    real, allocatable :: lon(:,:), lat(:,:), icec(:,:), std(:, :)
    integer, allocatable :: flag(:,:)

    integer :: i, j, nobs

    print *, 'reading "', trim(fname), '"...'

    inquire(file = trim(fname), exist = ex)
    if (.not. ex) then
       print *, 'ERROR: file "', trim(fname), '" not found'
       stop
    end if

    call nfw_open(fname, nf_nowrite, ncid)
    call nfw_inq_dimid(fname, ncid, 'xc', xc_id)
    call nfw_inq_dimid(fname, ncid, 'yc', yc_id)
    call nfw_inq_dimlen(fname, ncid, xc_id, nx)
    call nfw_inq_dimlen(fname, ncid, yc_id, ny)
    print *, '  nx = ', nx
    print *, '  ny = ', ny
    allocate(lon(nx, ny))
    allocate(lat(nx, ny))
    allocate(icec(nx, ny))
    allocate(std(nx, ny))
    allocate(flag(nx, ny))
    call nfw_inq_varid(fname, ncid, 'lon', lon_id)
    call nfw_inq_varid(fname, ncid, 'lat', lat_id)
    call nfw_inq_varid(fname, ncid, 'ice_conc', icec_id)
    call nfw_inq_varid(fname, ncid, 'standard_error', std_id)
    call nfw_inq_varid(fname, ncid, 'status_flag', flag_id)
    call nfw_get_var_double(fname, ncid, lon_id, lon)
    call nfw_get_var_double(fname, ncid, lat_id, lat)
    call nfw_get_var_double(fname, ncid, icec_id, icec)
    call nfw_get_var_double(fname, ncid, std_id, std)
    call nfw_get_var_int(fname, ncid, flag_id, flag)
    call nfw_close(fname, ncid)

    print *, 'filling the measurements array...'

    allocate(data(nx * ny))

    ! 0.995 is the max allowed by the model
    where (99.5d0 <= icec .and. icec <= 100.0d0)
       icec = 99.5d0
    end where

    nobs = 0
    do j = 1, ny
       do i = 1, nx
          nobs = nobs + 1
          if (flag(i, j) /= 0) then
             data(nobs) % status = .false.
             cycle
          end if
          data(nobs) % id = 'ICEC'
          data(nobs) % d = icec(i, j) * 1d-4
          data(nobs) % var = max(1d-8 * std(i, j) ** 2, 0.01d0 + (0.5d0 - abs(0.5d0 - data(nobs) % d)) ** 2)
          data(nobs) % ipiv = i
          data(nobs) % jpiv = j
          data(nobs) % lon = lon(i, j)
          data(nobs) % lat = lat(i, j)
          data(nobs) % a1 = 1e10
          data(nobs) % a2 = 1e10
          data(nobs) % a3 = 1e10
          data(nobs) % a4 = 1e10
          data(nobs) % ns = 1
          data(nobs) % date = 0
          data(nobs) % depth = 0.0
          data(nobs) % status = .true.
       end do
    end do
    print *, '  ', nobs, 'primary ICEC observations'
    print *, '  ', minval(data % d), ' <= icec <= ', maxval(data % d)

    gr = default_grid
    gr % nx = nx
    gr % ny = ny
    gr%reg = .true.
    gr % order = 2
    gr%ux = '10 km'
    gr%uy = '10 km'
    gr%set = .true.

    deallocate(lat, lon, icec, std, flag)
  end subroutine read_metno_icec

end module m_read_metno_icec
