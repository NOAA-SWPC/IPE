PROGRAM IPE_Driver

USE IPE_Precision
USE IPE_Model_Class

IMPLICIT NONE


  TYPE( IPE_Model ) :: ipe
  LOGICAL           :: init_success, fileExists
  REAL(prec)        :: t0_ipe, t1_ipe
  REAL(prec)        :: t2, t1
  INTEGER           :: i, rc
  CHARACTER(200)    :: init_file

    CALL ipe % Build( init_success )

    IF( init_success )THEN

      init_file = "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5"
      INQUIRE( FILE = TRIM(init_file), EXIST = fileExists )

        IF( fileExists )THEN
          CALL ipe % Initialize( init_file, init_success )
        ELSE
          PRINT*, ' IPE_Model_Class : Build : File not found : '//init_file
        ENDIF

      CALL ipe % time_tracker % Update( ipe % parameters % start_time )

      DO i = 1, ipe % parameters % n_model_updates

        IF( ipe % mpi_layer % rank_id == 0 ) PRINT*, 'Starting time loop at time :', ipe % time_tracker % DateStamp( )

        CALL CPU_TIME(t1)
        t0_ipe = ipe % parameters % start_time + REAL(i-1,prec)*ipe % parameters % file_output_frequency
        t1_ipe = ipe % parameters % start_time + REAL(i,prec)*ipe % parameters % file_output_frequency

        CALL ipe % Update( t0_ipe, t1_ipe )

        CALL CPU_TIME(t2)

#ifndef HAVE_MPI
!        CALL ipe % eldyn % Write_NetCDF( ipe % time_tracker, "eldyn.apex."//ipe % time_tracker % DateStamp( )//".nc", rc )
#endif

        CALL ipe % Write_to_HDF5( "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".h5", rc )
        IF( ipe % mpi_layer % rank_id == 0 )THEN
          write(6,*) '***********************************************'
          write(6,*) '*                                             *'
          write(6,*) '*       Total Update Time : ', t2-t1,'        *'
          write(6,*) '*                                             *'
          write(6,*) '***********************************************'
        ENDIF


      ENDDO

      CALL ipe % Trash( )

    ELSE

      PRINT*, '  IPE Model initialization unsuccessful'

    ENDIF


END PROGRAM IPE_Driver
