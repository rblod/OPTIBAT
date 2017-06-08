!>@file   par_wkb.F90
!>@brief Wave physics : parameters
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE par_wkb

   USE wkbcst
   
   IMPLICIT NONE
   PUBLIC
 
   LOGICAL :: ln_rst
   INTEGER :: jpi, jpj
   INTEGER :: nitermax, nit000, nitend, nstock, nbstp
   REAL(wp) :: rdx, rdt
   REAL(wp) :: eps

   CHARACTER(lc) :: cn_dirin, cn_dirout
   CHARACTER(lc) :: cn_filein, cn_fileout, cn_bryin

   LOGICAL :: ln_perio
   LOGICAL :: ln_anabry=.TRUE.
   LOGICAL :: ln_brywest, ln_bryeast, ln_brysouth, ln_brynorth

   REAL(wp) :: wkb_amp, wkb_prd, wkb_dir
   REAL(wp) :: wkb_roller, wkb_gam, wkb_btg, wkb_tide, wkb_rsb
   REAL(wp) :: bry_frq
      
   INTEGER :: wstp 
   INTEGER :: wnew   
   INTEGER :: iwave, winfo
   INTEGER ::first_time=0
   REAl(wp) :: thwave,av_wac,av_wkn  !!!QUAD
   
   INTEGER :: narea = 1  
   
   REAl(wp) :: Dcrit=0.2  
      
END MODULE par_wkb
