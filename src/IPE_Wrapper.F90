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
    LOGICAL               :: file_exists
    INTEGER               :: localrc
    INTEGER               :: mpicomm
    CHARACTER(LEN=30)     :: init_file
    TYPE(ESMF_VM)         :: vm_IPE
    TYPE(ESMF_Time)       :: currTime, startTime

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

    ! Set IPE internal clock
    CALL ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL IPE_SetClock( ipe, currTime, startTime, localrc )
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    ! Set IPE run mode to coupled
    ipe % forcing % coupled = .true.

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
    TYPE(ESMF_TimeInterval) :: timeStep
    INTEGER                 :: localrc
    CHARACTER(LEN=30)       :: hdf5_file

    ! begin
    IF (PRESENT(rc)) rc = ESMF_SUCCESS

    CALL ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
      timeStep=timeStep, rc=localrc)
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL IPE_SetClock( ipe, currTime, startTime, localrc )
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

    CALL ipe % Update( rc=localrc )
    IF( localrc /= IPE_SUCCESS ) THEN
      CALL ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="Error updating IPE", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      RETURN
    ENDIF

    currTime = currTime + timeStep
    CALL IPE_SetClock( ipe, currTime, startTime, localrc )
    IF( ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__, &
      rcToReturn=rc) ) RETURN  ! bail out

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


  SUBROUTINE IPE_SetClock( ipe, currTime, startTime, rc )

    TYPE(IPE_Model)              :: ipe
    TYPE(ESMF_Time), INTENT(IN)  :: currTime
    TYPE(ESMF_Time), INTENT(IN)  :: startTime
    INTEGER,         INTENT(OUT) :: rc

    ! local variables
    INTEGER                 :: year, month, day, hour, minute, second
    REAL(ESMF_KIND_R8)      :: elapsed_s
    TYPE(ESMF_TimeInterval) :: elapsedTime


    ! begin
    rc = ESMF_SUCCESS

    CALL ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__) ) RETURN  ! bail out

    CALL ipe % time_tracker % Set_Date( year, month, day )
    CALL ipe % time_tracker % Set_HourMinute( hour, minute )

    elapsedTime = currTime - startTime
    CALL ESMF_TimeIntervalGet(elapsedTime, s_r8=elapsed_s, rc=rc)
    IF( ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      FILE=__FILE__) ) RETURN  ! bail out

    ipe % time_tracker % elapsed_sec = REAL(elapsed_s, KIND=prec)

  END SUBROUTINE IPE_SetClock

END MODULE IPE_Wrapper
