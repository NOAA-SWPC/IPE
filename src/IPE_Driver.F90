PROGRAM IPE_Driver

USE IPE_Precision
USE IPE_Model_Class

IMPLICIT NONE

  TYPE( IPE_Model ) :: ipe
  LOGICAL           :: init_success, fileExists
  REAL(prec)        :: t0_ipe, t1_ipe
  REAL(prec)        :: t2, t1
  INTEGER           :: i, rc, stat
  CHARACTER(200)    :: init_file

    CALL ipe % Build( rc=rc )
    IF ( ipe_error_check( rc, msg="IPE model initialization unsuccessful", &
      line=__LINE__, file=__FILE__ ) ) CALL ipe % Trash()

    init_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5"
    INQUIRE( FILE = TRIM(init_file), EXIST = fileExists, iostat=stat )
    IF ( ipe_iostatus_check( rc, msg="Error inquiring about IPE initial state file "//init_file, &
      line=__LINE__, file=__FILE__ ) ) CALL ipe % Trash()

    IF ( ipe_status_check( fileExists, msg="IPE initial state file not found: "//init_file, &
      line=__LINE__, file=__FILE__ ) ) CALL ipe % Trash()

    CALL ipe % Initialize( init_file, rc=rc )
    IF ( ipe_iostatus_check( rc, msg="Error initializing IPE", &
      line=__LINE__, file=__FILE__ ) ) CALL ipe % Trash()

    CALL ipe % time_tracker % Update( ipe % parameters % start_time )

    DO i = 1, ipe % parameters % n_model_updates

      IF( ipe % mpi_layer % rank_id == 0 ) PRINT*, 'Starting time loop at time :', ipe % time_tracker % DateStamp( )

      CALL CPU_TIME(t1)
      t0_ipe = ipe % parameters % start_time + REAL(i-1,prec)*ipe % parameters % file_output_frequency
      t1_ipe = ipe % parameters % start_time + REAL(i,prec)*ipe % parameters % file_output_frequency

      CALL ipe % Update( t0_ipe, t1_ipe, rc=rc )
      IF ( ipe_iostatus_check( rc, msg="Error updating IPE model", &
        line=__LINE__, file=__FILE__ ) ) CALL ipe % Trash()

      CALL CPU_TIME(t2)

      CALL ipe % Write( rc=rc )
      IF ( ipe_iostatus_check( rc, msg="Error writing IPE output file", &
        line=__LINE__, file=__FILE__ ) ) CALL ipe % Trash()
      IF( ipe % mpi_layer % rank_id == 0 )THEN
        write(6,*) '***********************************************'
        write(6,*) '*                                             *'
        write(6,*) '*       Total Update Time : ', t2-t1,'        *'
        write(6,*) '*                                             *'
        write(6,*) '***********************************************'
      ENDIF

    ENDDO

    CALL ipe % Trash( )

END PROGRAM IPE_Driver
