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
   
   PUBLIC :: wkb_wri, wkb_read_ini, wkb_read_bry
   
   LOGICAL :: file_exist=.FALSE.
   INTEGER :: time_out
   
#define toto 
  
CONTAINS

   !> @brief Write wkb outputs
   
   SUBROUTINE wkb_wri
      INTEGER :: status, ncid, ji
      CHARACTER(len=lc),DIMENSION(3) :: dimnames
      CHARACTER(lc) :: clname
      REAL(kind=8), DIMENSION(jpi,jpj,1) :: ztab
  
      clname=TRIM(cn_dirout)//TRIM(cn_fileout)
      
      IF(.NOT. file_exist)THEN
         status = nf90_create(clname,NF90_WRITE,ncid)
         status = nf90_close(ncid)

         dimnames = (/ 'xi_rho      ','eta_rho     ', 'time_counter' /)
         CALL Write_Ncdf_dim(dimnames(1), clname, jpi)
         CALL Write_Ncdf_dim(dimnames(2), clname, jpj)
         CALL Write_Ncdf_dim(dimnames(3), clname, 0)
         time_out=1
         file_exist=.true.
      ENDIF
       
      CALL Write_Ncdf_var('xi_rho' , dimnames(1),clname, (/ (ji*rdx , ji=1,jpi)/), 'float')
      CALL Write_Ncdf_var('eta_rho', dimnames(2),clname, (/ (ji*rdx , ji=1,jpj)/), 'float')
      CALL Write_Ncdf_var('time_counter', (/dimnames(3)/),clname, (/real(nbstp)/),time_out,'float')
      !CALL Write_Ncdf_var('h'  , dimnames, clname, h  (:,:), 'double')
      ztab(:,:,1) = h(:,:) 
      CALL Write_Ncdf_var('h'  , dimnames, clname, ztab(:,:,1), time_out     ,'double')
      CALL Write_Ncdf_var('epb', dimnames, clname, wsb(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wac', dimnames, clname, wac(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('war', dimnames, clname, war(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wkx', dimnames, clname, wkx(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wky', dimnames, clname, wke(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('Cg' , dimnames, clname, wcg(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('zeta',dimnames,clname,zeta(:,:), time_out, 'double')
#ifdef toto
      CALL Write_Ncdf_var('hrm' , dimnames, clname, hrm(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wvn' , dimnames, clname, wvn(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wcr' , dimnames, clname, wcr(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wsr' , dimnames, clname, wsr(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('frq' , dimnames, clname, frq(:,:,wnew), time_out, 'double')
      CALL Write_Ncdf_var('wfc' , dimnames, clname, wfc(:,:,wnew), time_out, 'double')
#endif
       time_out=time_out+1
   
   END SUBROUTINE wkb_wri

   !> @brief Read wkb outputs

   SUBROUTINE wkb_read_ini
      CHARACTER(lc) :: clname
      
      clname=TRIM(cn_dirin)//TRIM(cn_filein)
      
      CALL Read_Ncdf_var('h', clname, h(:,:))

      clname=TRIM(cn_dirin)//TRIM(cn_rstin)
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

   END SUBROUTINE wkb_read_ini

   SUBROUTINE wkb_read_bry(kread)
   	INTEGER, INTENT(in) :: kread
      CHARACTER(lc) :: clname
      
      clname=TRIM(cn_dirin)//TRIM(cn_bryin)
        
      IF(ln_brywest) THEN
         CALL Read_Ncdf_var('tide_west'  , clname, hbry_west , kread)
         hbry_west_dt(:,1)=hbry_west(:,1)
         CALL Read_Ncdf_var('period_west', clname, perbry_west, kread)
         CALL Read_Ncdf_var('hs_west'    , clname, hsbry_west, kread)
         CALL Read_Ncdf_var('dir_west'   , clname, dirbry_west, kread)
      ENDIF
      
   END SUBROUTINE wkb_read_bry

END MODULE rw_wkb
