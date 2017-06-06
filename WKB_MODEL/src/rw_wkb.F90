!>@file   rw_wkb.F90
!>@brief Wave outputs
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE rw_wkb

   USE netcdf
   USE par_wkb
   USE wkbmod
   USE wkbcdf
    
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC :: wkb_wri, wkb_read
#define toto   
CONTAINS

   !> @brief Write wkb outputs
   
   SUBROUTINE wkb_wri
      INTEGER :: status, ncid, ji
      CHARACTER(len=lc),DIMENSION(2) :: dimnames
      CHARACTER(lc) :: clname

   
      clname=TRIM(cn_dirout)//TRIM(cn_fileout)
      status = nf90_create(clname,NF90_WRITE,ncid)
      status = nf90_close(ncid)

      dimnames = (/ 'xi_rho ','eta_rho' /)
      CALL Write_Ncdf_dim(dimnames(1), clname, jpi)
      CALL Write_Ncdf_dim(dimnames(2), clname, jpj)

      CALL Write_Ncdf_var('xi_rho' , dimnames(1),clname, (/ (ji*rdx , ji=1,jpi)/), 'float')
      CALL Write_Ncdf_var('eta_rho', dimnames(2),clname, (/ (ji*rdx , ji=1,jpj)/), 'float')
      CALL Write_Ncdf_var('h'  , dimnames, clname, h  (:,:)     ,'double')
      CALL Write_Ncdf_var('epb', dimnames, clname, wsb(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wac', dimnames, clname, wac(:,:,wnew), 'double')
      CALL Write_Ncdf_var('war', dimnames, clname, war(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wkx', dimnames, clname, wkx(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wky', dimnames, clname, wke(:,:,wnew), 'double')
      CALL Write_Ncdf_var('Cg' , dimnames, clname, wcg(:,:,wnew), 'double')
#ifdef toto
      CALL Write_Ncdf_var('hrm' , dimnames, clname, hrm(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wvn' , dimnames, clname, wvn(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wcr' , dimnames, clname, wcr(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wsr' , dimnames, clname, wsr(:,:,wnew), 'double')
      CALL Write_Ncdf_var('frq' , dimnames, clname, frq(:,:,wnew), 'double')
      CALL Write_Ncdf_var('wfc' , dimnames, clname, wfc(:,:,wnew), 'double')
#endif
   END SUBROUTINE wkb_wri

   !> @brief Read wkb outputs

   SUBROUTINE wkb_read
      CHARACTER(lc) :: clname
      
      clname=TRIM(cn_dirin)//TRIM(cn_filein)
      
      CALL Read_Ncdf_var('h', clname, h(:,:))
      IF ( ln_rst ) THEN
         CALL Read_Ncdf_var('epb', clname, wsb(:,:,wstp))
         CALL Read_Ncdf_var('wac', clname, wac(:,:,wstp))
         CALL Read_Ncdf_var('war', clname, war(:,:,wstp))
         CALL Read_Ncdf_var('wkx', clname, wkx(:,:,wstp))
         CALL Read_Ncdf_var('wky', clname, wke(:,:,wstp))
         CALL Read_Ncdf_var('Cg',  clname, wcg(:,:,wstp))
#ifdef toto
      CALL Read_Ncdf_var('hrm' , clname, hrm(:,:,wstp) )
      CALL Read_Ncdf_var('wvn' , clname, wvn(:,:,wstp) )
      CALL Read_Ncdf_var('wcr' , clname, wcr(:,:,wstp) )
      CALL Read_Ncdf_var('wsr' , clname, wsr(:,:,wstp) )
      CALL Read_Ncdf_var('frq' , clname, frq(:,:,wstp) )
      CALL Read_Ncdf_var('wfc' , clname, wfc(:,:,wstp))
      ENDIF   
#endif

   END SUBROUTINE wkb_read

END MODULE rw_wkb