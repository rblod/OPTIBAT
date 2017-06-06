module m_pivotp_micom
  use netcdf
  use nfw_mod
  use m_spherdist
  implicit none

contains
  ! F. Counillon (adapted from an algorithm of Mats Bentsen) 
  ! This subroutine search the pivot point for a given observations
  ! The search is linear moving toward the neigboring grid cell that minimize
  ! the distance to the obs. This search is in the worst case in O(n). The input
  ! ipiv, jpiv corresponds to the pivot point from the previous search. If there
  ! is a kind of order in the way the observation are given the search will be
  ! very fast.
  !
  subroutine pivotp_micom(lon, lat,modlon,modlat, ipiv, jpiv, nx, ny, min_r, max_r, &
  itw, jtw, its, jts, itn, jtn, ite, jte)
   real, intent(in) ::  lon, lat
   integer, intent(in) :: nx, ny
   real, intent(in), dimension(nx,ny) :: modlon,modlat, min_r, max_r, itw,jtw, &
   its, jts, itn, jtn, ite, jte
   integer, intent(inout) :: ipiv, jpiv
   real*8 :: min_d, d
   integer :: i, j, ito, jto

   min_d = spherdist(modlon(ipiv,jpiv), modlat(ipiv,jpiv), lon, lat)
   do while (min_d > min_r(ipiv,jpiv))

      ito = ipiv
      jto = jpiv

      i = itw(ito,jto)
      j = jtw(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = ite(ito,jto)
      j = jte(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = its(ito,jto)
      j = jts(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif
      i = itn(ito,jto)
      j = jtn(ito,jto)
      d = spherdist(modlon(i,j), modlat(i,j), lon, lat)
      if (d.lt.min_d) then
        ipiv = i
        jpiv = j
        min_d = d
      endif

      if (ipiv == ito .and. jpiv == jto) exit

   enddo
end subroutine pivotp_micom

subroutine ini_pivotp(modlon,modlat, nx, ny, min_r, max_r, itw, jtw, itn, jtn, its, jts, ite, jte)
   integer, intent(in) :: nx, ny
   real, intent(in), dimension(nx,ny) :: modlon,modlat
   real, intent(out), dimension(nx,ny):: min_r, max_r, itw, jtw, &
      itn, jtn, its, jts, ite, jte
   integer :: ncid,vITW_ID, vJTW_ID, vITE_ID ,vJTE_ID ,vITS_ID,  &
      vJTS_ID, vITN_ID, vJTN_ID, vVCLON_ID, vVCLAT_ID, vUCLON_ID &
      , vUCLAT_ID, vTCLON_ID, vTCLAT_ID ,ii,jj

   real, dimension(nx,ny,4) :: vclon, vclat, uclon, uclat, tclon, tclat 
   call nfw_open('grid.nc', nf_nowrite, ncid)

   call nfw_inq_varid('grid.nc', ncid,'inw' ,vITW_ID)
   call nfw_inq_varid('grid.nc', ncid,'jnw' ,vJTW_ID)
   call nfw_inq_varid('grid.nc', ncid,'ine' ,vITE_ID)
   call nfw_inq_varid('grid.nc', ncid,'jne' ,vJTE_ID)
   call nfw_inq_varid('grid.nc', ncid,'ins' ,vITS_ID)
   call nfw_inq_varid('grid.nc', ncid,'jns' ,vJTS_ID)
   call nfw_inq_varid('grid.nc', ncid,'inn' ,vITN_ID)
   call nfw_inq_varid('grid.nc', ncid,'jnn' ,vJTN_ID)
   call nfw_inq_varid('grid.nc', ncid,'vclon' ,vVCLON_ID)
   call nfw_inq_varid('grid.nc', ncid,'vclat' ,vVCLAT_ID)
   call nfw_inq_varid('grid.nc', ncid,'uclon' ,vUCLON_ID)
   call nfw_inq_varid('grid.nc', ncid,'uclat' ,vUCLAT_ID)
   call nfw_inq_varid('grid.nc', ncid,'pclon' ,vTCLON_ID)
   call nfw_inq_varid('grid.nc', ncid,'pclat' ,vTCLAT_ID)
   call nfw_get_var_double('grid.nc', ncid, vITW_ID, itw)
   call nfw_get_var_double('grid.nc', ncid, vJTW_ID, jtw)
   call nfw_get_var_double('grid.nc', ncid, vITS_ID, its)
   call nfw_get_var_double('grid.nc', ncid, vJTS_ID, jts)
   call nfw_get_var_double('grid.nc', ncid, vITN_ID, itn)
   call nfw_get_var_double('grid.nc', ncid, vJTN_ID, jtn)
   call nfw_get_var_double('grid.nc', ncid, vITE_ID, ite)
   call nfw_get_var_double('grid.nc', ncid, vJTE_ID, jte)
   call nfw_get_var_double('grid.nc', ncid, vVCLON_ID, vclon)
   call nfw_get_var_double('grid.nc', ncid, vVCLAT_ID, vclat)
   call nfw_get_var_double('grid.nc', ncid, vUCLON_ID, uclon)
   call nfw_get_var_double('grid.nc', ncid, vUCLAT_ID, uclat)
   call nfw_get_var_double('grid.nc', ncid, vTCLON_ID, tclon)
   call nfw_get_var_double('grid.nc', ncid, vTCLAT_ID, tclat)
   do jj = 1, ny
         do ii = 1, nx
            min_r(ii,jj) =                                   &
               min(spherdist(vclon(ii,jj,4), vclat(ii,jj,4), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(vclon(ii,jj,3), vclat(ii,jj,3), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(uclon(ii,jj,2), uclat(ii,jj,2), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(uclon(ii,jj,3), uclat(ii,jj,3), &
                              modlon(ii,jj), modlat(ii,jj)))
            max_r(ii,jj) =                                   &
               max(spherdist(tclon(ii,jj,1), tclat(ii,jj,1), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(tclon(ii,jj,2), tclat(ii,jj,2), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(tclon(ii,jj,3), tclat(ii,jj,3), &
                              modlon(ii,jj), modlat(ii,jj)), &
                   spherdist(tclon(ii,jj,4), tclat(ii,jj,4), &
                              modlon(ii,jj), modlat(ii,jj)))
         enddo
   enddo
end subroutine ini_pivotp



end module m_pivotp_micom
