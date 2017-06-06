module m_read_OOPS_gridded
! Ehouarn !
! Reads CHLA1 data and grid from the GlobColour files !
#ifdef OOPS
contains

 integer function get_nobs_OOPS_QG(nx,ny,nechx,nechy)
   use mod_measurement
   use nfw_mod
   implicit none  
   integer :: nx,ny,nechx,nechy
   
   character(80) :: memfile
   integer :: ncid, var_id,k,j
   
   real,dimension(nx,ny)::fld
   character(3)::cmonth
 
   nobs=-1
   memfile='observations-OOPS-QG.nc'
   call nfw_open(memfile, nf_nowrite, ncid)
   
   call nfw_inq_varid(memfile,ncid,'stream' ,var_id)
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

#endif
end
