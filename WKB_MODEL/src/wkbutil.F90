!>@file wkbutil.F90
!>@brief Library for misc. utilities (control print, etc ...)
!>@author R. Benshila (CNRS)
!>@version 0.0001

MODULE wkbutil

   USE par_wkb

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam

   INTEGER       ::   nstop = 0             !: error flag (=number of reason for a premature stop run)
   INTEGER       ::   nwarn = 0             !: warning flag (=number of warning found during the run)
   CHARACTER(lc) ::   cform_err = "(/,' ===>>> : E R R O R',     /,'         ===========',/)"       !:
   CHARACTER(lc) ::   cform_war = "(/,' ===>>> : W A R N I N G', /,'         ===============',/)"   !:

   CONTAINS
   
      SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5 ,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_opa  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the error number (nstop) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nstop = nstop + 1
      IF(lwp) THEN
         WRITE(numout,cform_err)
         IF( PRESENT(cd1 ) )   WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) )   WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) )   WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) )   WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) )   WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) )   WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) )   WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) )   WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) )   WRITE(numout,*) cd9
         IF( PRESENT(cd10) )   WRITE(numout,*) cd10
      ENDIF
      CALL FLUSH(numout    )
      !
      IF( cd1 == 'STOP' ) THEN
         IF(lwp) WRITE(numout,*)  'huge E-R-R-O-R : immediate stop'
         STOP
         !CALL mppstop()
      ENDIF
      !
   END SUBROUTINE ctl_stop


   SUBROUTINE ctl_warn( cd1, cd2, cd3, cd4, cd5,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_warn  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the warning number (nwarn) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nwarn = nwarn + 1
      IF(lwp) THEN
         WRITE(numout,cform_war)
         IF( PRESENT(cd1 ) ) WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) ) WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) ) WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) ) WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) ) WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) ) WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) ) WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) ) WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) ) WRITE(numout,*) cd9
         IF( PRESENT(cd10) ) WRITE(numout,*) cd10
      ENDIF
      CALL FLUSH(numout)
      !
   END SUBROUTINE ctl_warn


   SUBROUTINE ctl_opn( knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_opn  ***
      !!
      !! ** Purpose :   Open file and check if required file is available.
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(  out) ::   knum      ! logical unit to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdfile    ! file name to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdstat    ! disposition specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdform    ! formatting specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdacce    ! access specifier
      INTEGER          , INTENT(in   ) ::   klengh    ! record length
      INTEGER          , INTENT(in   ) ::   kout      ! number of logical units for write
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      INTEGER, OPTIONAL, INTENT(in   ) ::   karea     ! proc number
      !!
      CHARACTER(len=80) ::   clfile
      INTEGER           ::   iost
      !!----------------------------------------------------------------------

      ! adapt filename
      ! ----------------
      clfile = TRIM(cdfile)
      IF( PRESENT( karea ) ) THEN
         IF( karea > 1 )   WRITE(clfile, "(a,'_',i4.4)") TRIM(clfile), karea-1
      ENDIF
      knum=get_unit()

      iost=0
      IF( cdacce(1:6) == 'DIRECT' )  THEN
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
      ELSE
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat             , ERR=100, IOSTAT=iost )
      ENDIF
      IF( iost == 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*) '     file   : ', clfile,' open ok'
            WRITE(kout,*) '     unit   = ', knum
            WRITE(kout,*) '     status = ', cdstat
            WRITE(kout,*) '     form   = ', cdform
            WRITE(kout,*) '     access = ', cdacce
            WRITE(kout,*)
         ENDIF
      ENDIF
100   CONTINUE
      IF( iost /= 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) ' ===>>>> : bad opening file: ', clfile
            WRITE(kout,*) ' =======   ===  '
            WRITE(kout,*) '           unit   = ', knum
            WRITE(kout,*) '           status = ', cdstat
            WRITE(kout,*) '           form   = ', cdform
            WRITE(kout,*) '           access = ', cdacce
            WRITE(kout,*) '           iostat = ', iost
            WRITE(kout,*) '           we stop. verify the file '
            WRITE(kout,*)
         ENDIF
         STOP 'ctl_opn bad opening'
      ENDIF

   END SUBROUTINE ctl_opn

   SUBROUTINE ctl_nam ( kios, cdnam, ldwp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_nam  ***
      !!
      !! ** Purpose :   Informations when error while reading a namelist
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(inout) ::   kios      ! IO status after reading the namelist
      CHARACTER(len=*) , INTENT(in   ) ::   cdnam     ! group name of namelist for which error occurs
      CHARACTER(len=4)                 ::   clios     ! string to convert iostat in character for print
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      !!----------------------------------------------------------------------

      ! 
      ! ----------------
      WRITE (clios, '(I4.0)') kios, ldwp
      IF( kios < 0 ) THEN         
         CALL ctl_warn( 'W A R N I N G:  end of record or file while reading namelist ' &
 &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF

      IF( kios > 0 ) THEN
         CALL ctl_stop( 'E R R O R :   misspelled variable in namelist ' &
 &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      kios = 0
      RETURN
      
   END SUBROUTINE ctl_nam

   INTEGER FUNCTION get_unit()
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  get_unit  ***
      !!
      !! ** Purpose :   return the index of an unused logical unit
      !!----------------------------------------------------------------------
      LOGICAL :: llopn
      !!----------------------------------------------------------------------
      !
      get_unit = 15   ! choose a unit that is big enough then it is not already used in NEMO
      llopn = .TRUE.
      DO WHILE( (get_unit < 998) .AND. llopn )
         get_unit = get_unit + 1
         INQUIRE( unit = get_unit, opened = llopn )
      END DO
      IF( (get_unit == 999) .AND. llopn ) THEN
         CALL ctl_stop( 'get_unit: All logical units until 999 are used...' )
         get_unit = -1
      ENDIF
      !
   END FUNCTION get_unit

   
END MODULE wkbutil
