module m_get_waves_fld
! KAL -- This routine reads one of the fields from the model, specified
! KAL -- by name, vertical level and time level 
! KAL -- This routine is really only effective for the new restart files.
!use netcdf
use nfw_mod
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine get_waves_fld_new(memfile,fld,iens,cfld,vlevel,tlevel,nx,ny)
   use mod_raw_io
#if defined (QMPI)
   use qmpi, only : qmpi_proc_num
#else
   use qmpi_fake
#endif
   implicit none
   integer,      intent(in)            :: nx,ny  ! Grid dimension
   integer,      intent(in)            :: iens   ! Ensemble member to read
   real, dimension(nx,ny), intent(out) :: fld    ! output fld
   character(len=*), intent(in)        :: memfile! base name of input files
   character(len=*), intent(in)        :: cfld   ! name of fld
   integer, intent(in)                 :: tlevel ! time level
   integer, intent(in)                 :: vlevel ! vertical level

   real, dimension(nx,ny) :: readfld
   integer ::  ncid, vFIELD_ID!,ex
   integer, allocatable :: ns(:), nc(:)
   real::att(2)
   logical ::ex
   
   inquire(file=trim(memfile)//'.nc',exist=ex)
   if (ex) then
     ! Reading the observation file of satellite
     call nfw_open(trim(memfile)//'.nc', nf_nowrite, ncid)
     call nfw_inq_varid(trim(memfile)//'.nc', ncid,trim(cfld),vFIELD_ID)
 
     allocate(ns(2))
     allocate(nc(2))
     ns(1)=1
     ns(2)=1
     	
     nc(1)=nx
     nc(2)=ny
     
     call nfw_get_vara_double(trim(memfile)//'.nc', ncid, vFIELD_ID, ns, nc, readfld)
     fld=readfld(:,:)
    
     call nfw_close(trim(memfile)//'.nc', ncid)
  else
     print*, 'ERROR: forecast file is missing '//trim(memfile)//'.nc'
     stop 
  endif

end subroutine
end module m_get_waves_fld


