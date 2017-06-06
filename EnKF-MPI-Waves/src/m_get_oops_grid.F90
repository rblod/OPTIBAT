module m_get_oops_grid
contains 
subroutine get_oops_grid(modlon, modlat, depths, mindx, meandx, nx, ny)
   implicit none
   integer, intent(in) :: nx,ny
   real, dimension(nx,ny), intent(out) :: modlon,modlat,depths
   real,intent(out)   :: mindx,meandx
   real::dx,dy
   integer :: i
   logical ::  ex
   
   inquire(file='grid.txt',exist=ex)
   if (ex) then
      ! Reading the grid file
      open(unit=50,FILE='grid.txt',FORM='formatted',STATUS='unknown')
      read(50,*)dx
      read(50,*)dy
      close(50)
      
      do i=1,nx
        modlon(i,:)=dx/2+(i-1)*dx
      enddo	
      do i=1,ny
        modlat(:,i)=dy/2+(i-1)*dy
      enddo
     
      mindx = min(dx,dy)
      meandx=mindx
      depths(:,:)=100.
   else
      print*, 'ERROR: file grid.nc is missing' 
      stop 
   endif
end subroutine  get_oops_grid
end module  m_get_oops_grid
