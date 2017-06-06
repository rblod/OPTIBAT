module m_get_oops_dim
contains 
subroutine get_oops_dim(nx,ny)
  ! use netcdf
   use nfw_mod

   implicit none
   integer, intent(out) :: nx,ny
   integer:: ncid, x_ID, y_ID

   character(len=7) tag7    
   logical ex

   inquire(file='forecast001.nc',exist=ex)
   if (ex) then
     ! Reading the grid file
      call nfw_open('forecast001.nc', nf_nowrite, ncid)
     ! Get dimension id in netcdf file ...
     call nfw_inq_dimid('forecast001.nc', ncid, 'nx', x_ID)
     call nfw_inq_dimid('forecast001.nc', ncid, 'ny', y_ID)
     !Get the dimension
     call nfw_inq_dimlen('forecast001.nc', ncid, x_ID, nx)
     call nfw_inq_dimlen('forecast001.nc', ncid, y_ID, ny)
     print *, 'The model dimension is :',nx,ny
   else
      stop 'ERROR: file grid.nc is missing'
   endif
end subroutine  get_oops_dim
end module  m_get_oops_dim
