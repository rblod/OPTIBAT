module m_put_waves_fld
#ifdef QMPI
use qmpi, only : master
#else
use qmpi_fake, only : master
#endif
!use netcdf
use nfw_mod
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! KAL - This is for the new file type
subroutine put_waves_fld(memfile,fld,iens,cfld,vlevel,tlevel,nx,ny)
   implicit none
   integer, intent(in) :: nx,ny
   integer,                intent(in)  :: iens   ! Ensemble member to read
   real, dimension(nx,ny), intent(in)  :: fld    ! output fld
   character(len=*),       intent(in)  :: memfile! base name of input files
   character(len=9),       intent(in)  :: cfld   ! name of fld
   integer,                intent(in)  :: tlevel ! time level
   integer,                intent(in)  :: vlevel ! vertical level

   integer :: ncid, vFIELD_ID
   integer, allocatable :: ns(:), nc(:)
   logical :: ex


   if (master .and. iens.eq.1) then
      print *,'Dumping  ',trim(cfld),vlevel
   endif
   inquire(file=trim(memfile)//'.nc',exist=ex)
     ! Reading the observation file of satellite
     call nfw_open(trim(memfile)//'.nc', nf_write, ncid)
     call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld),vFIELD_ID)

     
      allocate(ns(2))
      allocate(nc(2))
      ns(1)=1
      ns(2)=1
	 
      nc(1)=nx
      nc(2)=ny
      
      call nfw_put_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, fld)
     
     call nfw_close(trim(memfile)//'.nc', ncid)
end subroutine



end module m_put_waves_fld


