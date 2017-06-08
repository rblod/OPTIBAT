!>@file   wkbmod.F90
!>@brief Wave arrays delcaration and allocation
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE wkbmod

   USE par_wkb
   USE wkbutil
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC wkb_alloc
   
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: h

      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wkx
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wke
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wac
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: hrm
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: frq
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wcg
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wsb
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wvn
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wfc
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: war
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wcr
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: wsr

      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wdrx
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wdre
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wdsp
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wdrg
      REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: sup

     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: wacbry_west , warbry_west , wkxbry_west , wkebry_west
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: wacbry_east , warbry_east , wkxbry_east , wkebry_east
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: wacbry_north, warbry_north, wkxbry_north, wkebry_north
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: wacbry_south, warbry_south, wkxbry_south, wkebry_south

     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wacbry_west_dt , warbry_west_dt , wkxbry_west_dt , wkebry_west_dt
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wacbry_east_dt , warbry_east_dt , wkxbry_east_dt , wkebry_east_dt
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wacbry_north_dt, warbry_north_dt, wkxbry_north_dt, wkebry_north_dt
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: wacbry_south_dt, warbry_south_dt, wkxbry_south_dt, wkebry_south_dt
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)  ::  hbry_west  , hbry_east  , hbry_north  , hbry_south
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)  ::  hsbry_west , hsbry_east , hsbry_north , hsbry_south
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)  ::  perbry_west, perbry_east, perbry_north, perbry_south
     REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)  ::  dirbry_west, dirbry_east, dirbry_north, dirbry_south
     
CONTAINS
	
   INTEGER FUNCTION wkb_alloc()

      INTEGER :: ierr(5)

      ALLOCATE( wkx(jpi,jpj,2), wke(jpi,jpj,2), wac(jpi,jpj,2),   &
                hrm(jpi,jpj,2), frq(jpi,jpj,2), wcg(jpi,jpj,2),   &
                wsb(jpi,jpj,2), wvn(jpi,jpj,2), wfc(jpi,jpj,2),   &
                war(jpi,jpj,2), wcr(jpi,jpj,2), wsr(jpi,jpj,2),   &
                h(jpi,jpj)                    , STAT=ierr(1)      )

      ALLOCATE( wdrx(jpi,jpj), wdre(jpi,jpj),   &
                wdsp(jpi,jpj), wdrg(jpi,jpj),   &
                sup(jpi,jpj),                      STAT=ierr(2)    )
                                                
      ALLOCATE( wacbry_west(jpj) , warbry_west(jpj) , wkxbry_west(jpj) , wkebry_west(jpj) ,   &
                wacbry_east(jpj) , warbry_east(jpj) , wkxbry_east(jpj) , wkebry_east(jpj) ,   &
                wacbry_north(jpi), warbry_north(jpi), wkxbry_north(jpi), wkebry_north(jpi),   &
                wacbry_south(jpi), warbry_south(jpi), wkxbry_south(jpi), wkebry_south(jpi),   &
                                                 STAT=ierr(3)    )
 
      ALLOCATE( wacbry_west_dt (jpj,2), warbry_west_dt (jpj,2), wkxbry_west_dt (jpj,2), wkebry_west_dt (jpj,2) ,  &
                wacbry_east_dt (jpj,2), warbry_east_dt (jpj,2), wkxbry_east_dt (jpj,2), wkebry_east_dt (jpj,2) ,  &
                wacbry_north_dt(jpi,2), warbry_north_dt(jpi,2), wkxbry_north_dt(jpi,2), wkebry_north_dt(jpi,2),   &
                wacbry_south_dt(jpi,2), warbry_south_dt(jpi,2), wkxbry_south_dt(jpi,2), wkebry_south_dt(jpi,2),   &
                STAT=ierr(4) )
      ALLOCATE( hbry_west  (jpj,1), hbry_east  (jpj,1), hbry_north  (jpi,1), hbry_south  (jpi,1), &
                hsbry_west (jpj,1), hsbry_east (jpj,1), hsbry_north (jpi,1), hsbry_south (jpi,1), &
                perbry_west(jpj,1), perbry_east(jpj,1), perbry_north(jpi,1), perbry_south(jpi,1), &
                dirbry_west(jpj,1), dirbry_east(jpj,1), dirbry_north(jpi,1), dirbry_south(jpi,1), &
                STAT=ierr(5) )
                   
      wkb_alloc = MAXVAL( ierr )
      IF( wkb_alloc /= 0 )   CALL ctl_warn('wkb_alloc: failed to allocate arrays')
      
   END FUNCTION wkb_alloc

END MODULE wkbmod
