module m_get_waves_dim
contains 
subroutine get_waves_dim(nx,ny)
  ! use netcdf
   use nfw_mod

   implicit none
   integer, intent(out) :: nx,ny
   integer:: ncid, x_ID, y_ID

   character(len=7) tag7    
   logical ex
   
   print*,'connard1'
   inquire(file='forecast001.nc',exist=ex)
   if (ex) then
     ! Reading the grid file
      call nfw_open('forecast001.nc', nf_nowrite, ncid)
     ! Get dimension id in netcdf file ...
     call nfw_inq_dimid('forecast001.nc', ncid, 'xi_rho', x_ID)
     call nfw_inq_dimid('forecast001.nc', ncid, 'eta_rho', y_ID)
     !Get the dimension
     call nfw_inq_dimlen('forecast001.nc', ncid, x_ID, nx)
     call nfw_inq_dimlen('forecast001.nc', ncid, y_ID, ny)
     print *, 'The model dimension is :',nx,ny
     call nfw_close('forecast001.nc',ncid)
   else
      stop 'ERROR: file grid.nc is missing'
   endif
end subroutine  get_waves_dim
end module  m_get_waves_dim
