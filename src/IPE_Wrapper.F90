MODULE IPE_Wrapper
   
  USE ESMF
  USE IPE_Precision
  USE IPE_Model_Class
  USE ipe_error_module

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: IPE_Model

  PUBLIC :: Initialize_IPE
  PUBLIC :: Update_IPE
  PUBLIC :: Finalize_IPE

CONTAINS
  
  SUBROUTINE Initialize_IPE ( ipe, clock, vm, rc )

    TYPE(IPE_Model)                      :: ipe
    TYPE(ESMF_Clock)                     :: clock
    TYPE(ESMF_VM), OPTIONAL              :: vm
    INTEGER,       OPTIONAL, INTENT(OUT) :: rc
    ! Local
    LOGICAL               :: init_success
    LOGICAL               :: file_exists
    INTEGER               :: localrc
    INTEGER               :: year, month, day, hour, minute, second
    INTEGER               :: mpicomm, mpierr
    INTEGER(ESMF_KIND_I4) :: mm, dd, h, m, s
    INTEGER(ESMF_KIND_I4) :: yystop, mmstop, ddstop, hstop, mstop, sstop
    CHARACTER(LEN=30)     :: init_file
    TYPE(ESMF_VM)         :: vm_IPE
    TYPE(ESMF_Time)       :: currTime

    ! begin
    IF (PRESENT(rc)) rc = ESMF_SUCCESS

    IF (PRESENT(vm)) THEN
      vm_IPE = vm
    ELSE
      CALL ESMF_VMGetCurrent(vm=vm_IPE, rc=localrc)
      IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        FILE=__FILE__, &
        rcToReturn=rc) ) RETURN  ! bail out
    END IF

    ! Obtain communicator
    CALL ESMF_VMGet(vm=vm_IPE, mpiCommunicator=mpicomm, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    ! Build IPE
    CALL ipe % Build(mpi_comm=mpicomm, rc=localrc)
    IF( localrc /= IPE_SUCCESS ) THEN
      CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error building IPE", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      RETURN
    ENDIF

    ! Set IPE run mode to coupled
    ipe % forcing % coupled = .true.

    ! Set IPE internal clock
    CALL ESMF_ClockGet(clock, currTime=currTime, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL ipe % time_tracker % Set_Date( year, month, day )
    CALL ipe % time_tracker % Set_HourMinute( hour, minute )
    ipe % time_tracker % elapsed_sec = 0.0_prec

    init_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5"
    INQUIRE( FILE = TRIM(init_file), EXIST = file_exists, IOSTAT = localrc )
    IF( localrc /= 0 ) THEN
      CALL ESMF_LogSetError(ESMF_RC_FILE_UNEXPECTED, &
        msg="Failed to inquire file "//init_file, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      RETURN
    ENDIF

    IF( file_exists )THEN
      CALL ipe % Initialize( init_file, rc=localrc )
      IF ( localrc /= IPE_SUCCESS ) THEN
        CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error initializing IPE", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        RETURN
      ENDIF
    ELSE
      PRINT*, ' IPE : File not found '//TRIM(init_file)
      CALL ESMF_LogSetError(ESMF_RC_FILE_OPEN, &
        msg="IPE: File "//trim(init_file)//" not found", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      RETURN
    ENDIF

  END SUBROUTINE Initialize_IPE

  
  SUBROUTINE Update_IPE( ipe, clock, rc )

    TYPE(IPE_Model)                :: ipe
    TYPE(ESMF_Clock)               :: clock
    INTEGER, OPTIONAL, INTENT(OUT) :: rc
     ! Local
    TYPE(ESMF_Time)         :: currTime, startTime
    TYPE(ESMF_TimeInterval) :: elapsedTime
    REAL(ESMF_KIND_R8)      :: elapsed_s
    INTEGER                 :: localrc
    INTEGER                 :: error
    INTEGER                 :: year, month, day, hour, minute, second
    CHARACTER(LEN=30)       :: hdf5_file

    ! begin
    IF (PRESENT(rc)) rc = ESMF_SUCCESS

    CALL ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL ESMF_TimeGet(currtime, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL ipe % time_tracker % Set_Date( year, month, day )
    CALL ipe % time_tracker % Set_HourMinute( hour, minute )

    elapsedTime = currTime - startTime
    CALL ESMF_TimeIntervalGet(elapsedTime, s_r8=elapsed_s, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    ipe % time_tracker % elapsed_sec = REAL(elapsed_s, KIND=prec)

    CALL ipe % Update( rc=localrc )
    IF( localrc /= IPE_SUCCESS ) THEN
      CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error updating IPE", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      RETURN
    ENDIF

    IF( MOD( ipe % time_tracker % elapsed_sec, ipe % parameters % file_output_frequency ) == 0.0_prec )THEN

      hdf5_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5"
      CALL ipe % Write( hdf5_file, rc=localrc )
      IF( localrc /= IPE_SUCCESS ) THEN
        CALL ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg="Error writing to HDF5 file "//hdf5_file, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        RETURN
      ENDIF

    ENDIF

  END SUBROUTINE Update_IPE

  
  SUBROUTINE Finalize_IPE( ipe, rc )

    TYPE(IPE_Model)                :: ipe
    INTEGER, OPTIONAL, INTENT(OUT) :: rc

    ! local variables
    integer :: localrc
 
    ! begin
    IF (PRESENT(rc)) rc = ESMF_SUCCESS

    CALL ESMF_LogWrite("Finalize_IPE start:", ESMF_LOGMSG_INFO, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL ipe % Trash( rc=localrc )
    IF( localrc /= IPE_SUCCESS ) THEN
      CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error finalizing IPE", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      RETURN
    ENDIF

    CALL ESMF_LogWrite("Finalize_IPE end:", ESMF_LOGMSG_INFO, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

 END SUBROUTINE Finalize_IPE

END MODULE IPE_Wrapper
