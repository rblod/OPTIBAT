module m_get_waves_grid
contains 
subroutine get_waves_grid(modlon, modlat, depths, mindx, meandx, nx, ny)  
   use nfw_mod
   implicit none
   integer, intent(in) :: nx,ny
   real, dimension(nx,ny), intent(out) :: modlon,modlat,depths
   real,intent(out)   :: mindx,meandx
   real, dimension(:),allocatable :: eta,xi
   real::dx,dy
   integer :: i
   logical ::  ex
    integer:: ncid, x_ID, y_ID,var_id
   
    inquire(file='forecast001.nc',exist=ex)
    if (ex) then
      ! Reading the grid file
        allocate(eta(ny),xi(nx))
        call nfw_open('forecast001.nc', nf_nowrite, ncid)
     
        call nfw_inq_varid('forecast001.nc', ncid,'xi_rho',var_id)
        call nfw_get_var_double('forecast001.nc', ncid, var_id, xi)
    
        call nfw_inq_varid('forecast001.nc', ncid,'eta_rho',var_id)
        call nfw_get_var_double('forecast001.nc', ncid, var_id, eta)
        call nfw_close('forecast001.nc',ncid)

      do i=1,nx
!        modlon(i,:)=10+20*(xi(i)-1)
!!!!        modlon(i,:)=xi(i)
         modlat(i,:)=xi(i)
      enddo	
      do i=1,ny
!        modlat(:,i)=10+20*(eta(i)-1)
!!!        modlat(:,i)=eta(i)
        modlon(:,i)=eta(i)
      enddo
     
      mindx = 20
      meandx=mindx
      depths(:,:)=100.
      deallocate(xi,eta)
   else
      print*, 'ERROR: file grid.nc is missing' 
      stop 
   endif
end subroutine  get_waves_grid
end module  m_get_waves_grid
