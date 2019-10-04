MODULE IPE_Wrapper
   
  USE ESMF
  USE IPE_Precision
  USE IPE_Model_Class

  IMPLICIT NONE

  TYPE( IPE_Model ) :: ipe
  INTEGER           :: lps, lpe, mps, mpe

  PUBLIC

CONTAINS
  
  SUBROUTINE Initialize_IPE ( clock, rc )

    TYPE(ESMF_Clock)     :: clock
    INTEGER, intent(out) :: rc
    ! Local
    LOGICAL               :: init_success
    TYPE(ESMF_VM)         :: vm_IPE
    TYPE(ESMF_Time)       :: currTime
    INTEGER               :: year, month, day, hour, minute, second
    INTEGER               :: mpicomm, mpierr
    INTEGER(ESMF_KIND_I4) :: mm, dd, h, m, s
    INTEGER(ESMF_KIND_I4) :: yystop, mmstop, ddstop, hstop, mstop, sstop
    LOGICAL               :: file_exists
    CHARACTER(LEN=30)     :: init_file

    rc = ESMF_SUCCESS

    CALL ESMF_VMGetCurrent(vm=vm_IPE, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    ! Obtain communicator
    CALL ESMF_VMGet(vm=vm_IPE, mpiCommunicator=mpicomm, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    ! Build IPE
    PRINT*, ' IPE : Initializing IPE'
    CALL ipe % Build( init_success, mpi_comm = mpicomm )
    PRINT*, ' IPE : Done Initializing IPE'

    IF( .not. init_success )THEN
      CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error building IPE", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      RETURN
    ENDIF

    ! Set IPE run mode to coupled
    ipe % forcing % coupled = .true.

    ! Set local (lp,mp) bounds
    lps = 1
    lpe = ipe % grid % NLP
    mps = ipe % grid % mp_low
    mpe = ipe % grid % mp_high

    ! Set IPE internal clock
    CALL ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    CALL ESMF_TimeGet(currTime, yy=year,mm=month,dd=day,h=hour,m=minute,s=second,rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    CALL ipe % time_tracker % Set_Date( year, month, day )
    CALL ipe % time_tracker % Set_HourMinute( hour, minute )
    ipe % time_tracker % elapsed_sec = 0.0_prec

    init_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5"
    INQUIRE( FILE = TRIM(init_file), EXIST = file_exists )

    IF( file_exists )THEN
      CALL ipe % Initialize( init_file, init_success )
      IF ( .not. init_success ) THEN
        CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error initializing IPE", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        RETURN
      ENDIF
    ELSE
      PRINT*, ' IPE : File not found '//TRIM(init_file)
      CALL ESMF_LogSetError(ESMF_RC_FILE_OPEN, msg="IPE: File "//trim(init_file)//" not found", &
        line=__LINE__, file=__FILE__, rcToReturn=rc)
      RETURN
    ENDIF

  END SUBROUTINE Initialize_IPE

  
  SUBROUTINE Update_IPE( clock, rc )

    TYPE(ESMF_Clock)     :: clock
    INTEGER, intent(out) :: rc
     ! Local
    TYPE(ESMF_Time)         :: currTime, startTime
    TYPE(ESMF_TimeInterval) :: elapsedTime
    REAL(ESMF_KIND_R8)      :: elapsed_s
    INTEGER                 :: error
    INTEGER                 :: year, month, day, hour, minute, second
    CHARACTER(LEN=30)       :: hdf5_file

    rc = ESMF_SUCCESS

    CALL ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    CALL ESMF_TimeGet(currtime, yy=year, mm=month, dd=day, h=hour, m=minute, s=second, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    CALL ipe % time_tracker % Set_Date( year, month, day )
    CALL ipe % time_tracker % Set_HourMinute( hour, minute )

    elapsedTime = currTime - startTime
    CALL ESMF_TimeIntervalGet(elapsedTime, s_r8=elapsed_s, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    ipe % time_tracker % elapsed_sec = REAL(elapsed_s, KIND=prec)

    CALL ipe % Update( )

    IF( MOD( ipe % time_tracker % elapsed_sec, ipe % parameters % file_output_frequency ) == 0.0_prec )THEN

      hdf5_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5"
      CALL ipe % Write_to_HDF5( hdf5_file, error )
      IF ( error /= 0 ) THEN
        CALL ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg="Error writing to HDF5 file "//hdf5_file, &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        RETURN
      ENDIF

    ENDIF

  END SUBROUTINE Update_IPE

  
  SUBROUTINE Finalize_IPE( rc )

    INTEGER, intent(out) :: rc  
 
    rc = ESMF_SUCCESS

    CALL ESMF_LogWrite("Finalize_IPE start:", ESMF_LOGMSG_INFO, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

    CALL ipe % Trash( )

    CALL ESMF_LogWrite("Finalize_IPE end:", ESMF_LOGMSG_INFO, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__,FILE=__FILE__) ) RETURN  ! bail out

 END SUBROUTINE Finalize_IPE

END MODULE IPE_Wrapper

