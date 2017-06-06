! Bilinear coeffisients are calculated witin this program for all the 
! observation points.
! Only wet points are stored in the observation.uf file

module m_get_def_wet_point

  implicit none

  integer, parameter, private :: STRLEN = 512
  character(STRLEN), parameter, private :: MEANSSHFNAME = "meanssh.uf"
  
  private read_mean_ssh
  private land_nearby

contains 

  subroutine get_def_wet_point(obs, data, gr, depths, modlat, modlon, nrobs, nx, ny)
    ! Program converts to a general format readable by ENKF (observations.uf)
    use mod_measurement
    use mod_grid
    ! Functions to be used
    use m_oldtonew
    use m_bilincoeff
    use m_pivotp_micom
    use m_confmap
    use m_spherdist

    integer, intent(in) :: nx, ny
    type (measurement), intent(in) :: data(:)
    type (measurement), intent(inout)   :: obs(:)
    type (grid), intent(in) :: gr      ! observations grid
    real, dimension(nx, ny), intent(in)  ::  depths, modlat, modlon
    integer, intent(out) :: nrobs
    integer, parameter :: maxobs = 1441 * 722 !2*400*600 ! maximum number of observations

    real, dimension(nx, ny) :: mean_ssh
    integer k, imin, imax, jmin, jmax
    integer ipiv, jpiv, nsupport, nsmin, nsmax
    real :: x0, y0
    real wetsill, griddiag, mingridsize, minobssize
    real, dimension(nx,ny) :: min_r, max_r, itw, jtw, its, jts, itn, jtn, ite, jte

    logical wet

    ! gr = default_grid
    nrobs = 0; 
    nsmin = maxobs; 
    nsmax = 0
    mingridsize = 1.E+10; 
    minobssize = 1.E+10    ! in meters
    call ini_pivotp(modlon,modlat, nx, ny, min_r, max_r, itw, jtw, itn, jtn, its, jts, ite, jte)
    ipiv=1
    jpiv=1


    !Calculate pivot points
    !Find wet points (all neigbours in water)
    !Find the points with defined data value
    !Put the data into the obs data structture
    !Compute bilinear coefficients
    


    do k = 1, gridpoints(gr)
       if (data(k) % id .eq. 'SLA' .or. data(k) % id .eq. 'sla' .or. &
            data(k) % id.eq. 'SSH' .or. data(k)%id .eq. 'ssh' .or.&
            data(k)%id.eq.'TSLA') then
          call read_mean_ssh(mean_ssh, nx, ny)
          wetsill = 200.   ! Discarding data in shallow waters
       else
          wetsill=10.
       endif
       call pivotp_micom(data(k)%lon, data(k)%lat, modlon, modlat, ipiv, jpiv, &
          nx, ny, min_r, max_r,itw, jtw, itn, jtn, its, jts, ite, jte)
       if (depths(ipiv, jpiv) < wetsill ) cycle
       wet = data(k) % status ! Discards inconsistent/Fill values
       if (data(k) % id .eq. 'SLA' .or. data(k) % id .eq. 'sla' .or.&
            data(k) % id .eq. 'TSLA') then
          wet = wet .and. (mean_ssh(ipiv, jpiv) < 990)
          wet = wet .and. .not. land_nearby(nx, ny, mean_ssh, modlon, modlat,&
               ipiv, jpiv, data(k) % lon, data(k) % lat)
       endif

       if(.not. undefined(data(k) % d, gr) .and. wet) then
          nrobs = nrobs + 1
          obs(nrobs) = data(k)
          obs(nrobs) % ipiv = ipiv
          obs(nrobs) % jpiv=  jpiv
          obs(nrobs) % status = .true. ! Wet obs
          obs(nrobs) % ns = 0    ! point measurements have zero support
	  !ehouarn: call bilicoeff	 
       endif
    end do
    print*, 'Number of defined and wet observations: nrobs ', nrobs
    print*, 'Support (in nb of cells) between: ', nsmin, ' and ', nsmax
    print '(2(a,f8.3),a)', ' Minimum obs support: ', 0.001 * minobssize, &
         'km, min grid diagonal: ', 0.001 * mingridsize, ' km' 
  end subroutine get_def_wet_point


  subroutine read_mean_ssh(mean_ssh, nx, ny)
    integer, intent(in) :: nx, ny
    real, intent(out):: mean_ssh(nx, ny)
    logical :: exists

    inquire(file = trim(MEANSSHFNAME), exist = exists)
    if (.not. exists) then
       print *,'ERROR: read_mean_ssh(): file "', trim(MEANSSHFNAME), '" not found'
       stop
    end if
       
    open (10, file = trim(MEANSSHFNAME), status = 'unknown',form = 'unformatted')
    read (10) mean_ssh
    close (10)
  end subroutine read_mean_ssh


  logical function land_nearby(nx, ny, mean_ssh, modlon, modlat, ipiv, jpiv, obslon, obslat)
    use m_spherdist
    implicit none
    real, parameter :: Dis0 = 50.0d0
    integer, intent (in) :: nx, ny, ipiv, jpiv
    real, dimension(nx,ny), intent(in) :: mean_ssh, modlon, modlat
    real, intent (in) :: obslon,obslat 
    integer :: ii, jj, ncells
    real :: griddist

    land_nearby = .false.
    ncells = ceiling(Dis0 / spherdist(modlon(ipiv, jpiv), modlat(ipiv, jpiv),&
         modlon(ipiv, jpiv + 1), modlat(ipiv, jpiv + 1)))
    do jj = max(jpiv - ncells, 1), min(jpiv + ncells, ny)
       do ii = max(ipiv - ncells, 1), min(ipiv + ncells, nx)
          griddist = spherdist(modlon(ii, jj), modlat(ii, jj), obslon, obslat)
          if (mean_ssh(ipiv,jpiv) < 990 .and. griddist < Dis0) then
             land_nearby = .true.
             return
          end if
       enddo
    enddo
  end function land_nearby

end module m_get_def_wet_point
