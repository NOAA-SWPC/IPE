MODULE IPE_Model_Parameters_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Common_Routines
  USE IPE_MPI_Layer_Class


  IMPLICIT NONE

  TYPE IPE_Model_Parameters

    ! SpaceManagement
    CHARACTER(200) :: grid_file
    INTEGER        :: NLP
    INTEGER        :: NMP
    INTEGER        :: NPTS2D
    INTEGER        :: nFluxTube

    ! TimeStepping
    REAL(prec)    :: time_step
    REAL(prec)    :: start_time
    REAL(prec)    :: end_time
    REAL(prec)    :: msis_time_step
    CHARACTER(12) :: initial_timestamp

    ! Forcing
    REAL(prec)     :: solar_forcing_time_step
    INTEGER        :: f107_kp_size
    INTEGER        :: f107_kp_interval
    INTEGER        :: f107_kp_skip_size
    INTEGER        :: f107_kp_data_size
    INTEGER        :: f107_kp_read_in_start
    LOGICAL        :: use_f107_kp_file
    CHARACTER(200) :: f107_kp_file

    ! IPECAP
    REAL(prec) :: mesh_height_min
    REAL(prec) :: mesh_height_max
    INTEGER    :: mesh_fill
    INTEGER    :: mesh_write
    CHARACTER(200) :: mesh_write_file

    ! >> Fixed parameters
    REAL(prec) :: f107
    INTEGER    :: f107_flag
    REAL(prec) :: f107_81day_avg
    REAL(prec) :: kp
    INTEGER    :: kp_flag
    REAL(prec) :: kp_1day_avg
    REAL(prec) :: ap
    INTEGER    :: ap_flag
    REAL(prec) :: ap_1day_avg
    REAL(prec) :: nhemi_power
    INTEGER    :: nhemi_power_index
    REAL(prec) :: shemi_power
    INTEGER    :: shemi_power_index
    REAL(prec) :: solarwind_By
    REAL(prec) :: solarwind_angle
    REAL(prec) :: solarwind_velocity
    REAL(prec) :: solarwind_Bz
    REAL(prec) :: solarwind_density

    !FileIO
    LOGICAL :: read_apex_neutrals
    LOGICAL :: read_geographic_neutrals
    LOGICAL :: write_apex_neutrals
    LOGICAL :: write_geographic_neutrals
    LOGICAL :: write_geographic_eldyn
    LOGICAL :: write_apex_eldyn
    REAL(prec) :: file_output_frequency

    !ElDyn
    LOGICAL :: dynamo_efield

    INTEGER :: n_model_updates

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model_Parameters

  END TYPE IPE_Model_Parameters


CONTAINS

  SUBROUTINE Build_IPE_Model_Parameters( params, mpi_layer, read_success )

    IMPLICIT NONE

    CLASS( IPE_Model_Parameters ), INTENT(out) :: params
    TYPE( IPE_MPI_Layer ),         INTENT(in)  :: mpi_layer
    LOGICAL,                       INTENT(out) :: read_success

    ! Local
    CHARACTER(200) :: grid_file
    INTEGER        :: fUnit, ierr, iostatus
    LOGICAL        :: fileExists
    INTEGER        :: NLP, NMP, NPTS2D, nFluxTube
    REAL(prec)     :: solar_forcing_time_step
    REAL(prec)     :: time_step, start_time, end_time, msis_time_step
    CHARACTER(12)  :: initial_timestamp
    INTEGER        :: f107_kp_size
    INTEGER        :: f107_kp_interval
    INTEGER        :: f107_kp_skip_size
    INTEGER        :: f107_kp_data_size
    INTEGER        :: f107_kp_read_in_start
    CHARACTER(200) :: f107_kp_file
    LOGICAL        :: read_apex_neutrals
    LOGICAL        :: read_geographic_neutrals
    LOGICAL        :: write_apex_neutrals
    LOGICAL        :: write_geographic_neutrals
    LOGICAL        :: write_geographic_eldyn
    LOGICAL        :: write_apex_eldyn
    LOGICAL        :: dynamo_efield
    REAL(prec)     :: file_output_frequency
    REAL(prec)     :: mesh_height_min
    REAL(prec)     :: mesh_height_max
    INTEGER        :: mesh_fill
    INTEGER        :: mesh_write
    CHARACTER(200) :: mesh_write_file
    ! >> Fixed parameters
    REAL(prec) :: f107
    REAL(prec) :: f107_81day_avg
    INTEGER    :: f107_flag
    REAL(prec) :: kp
    INTEGER    :: kp_flag
    REAL(prec) :: kp_1day_avg
    REAL(prec) :: ap
    INTEGER    :: ap_flag
    REAL(prec) :: ap_1day_avg
    REAL(prec) :: nhemi_power
    INTEGER    :: nhemi_power_index
    REAL(prec) :: shemi_power
    INTEGER    :: shemi_power_index
    REAL(prec) :: solarwind_By
    REAL(prec) :: solarwind_angle
    REAL(prec) :: solarwind_velocity
    REAL(prec) :: solarwind_Bz
    REAL(prec) :: solarwind_density

    ! Communication buffers
    CHARACTER(LEN=200), DIMENSION( 4) :: sbuf
    INTEGER,            DIMENSION(24) :: ibuf
    REAL(prec),         DIMENSION(21) :: rbuf


    NAMELIST / SpaceManagement / grid_file, NLP, NMP, NPTS2D, nFluxTube
    NAMELIST / TimeStepping    / time_step, start_time, end_time, msis_time_step, initial_timestamp
    NAMELIST / Forcing         / solar_forcing_time_step, f107_kp_size, f107_kp_interval, f107_kp_skip_size, &
                                 f107_kp_data_size, f107_kp_read_in_start, f107_kp_file, f107, f107_flag, f107_81day_avg, &
                                 kp, kp_flag, kp_flag, kp_1day_avg, ap, ap_flag, ap_1day_avg, nhemi_power, &
                                 nhemi_power_index, shemi_power, shemi_power_index, solarwind_By, solarwind_angle, &
                                 solarwind_velocity, solarwind_Bz, solarwind_density
    NAMELIST / FileIO          / read_apex_neutrals, read_geographic_neutrals, write_apex_neutrals, write_geographic_neutrals, &
                                 write_geographic_eldyn, write_apex_eldyn, file_output_frequency
    NAMELIST / IPECAP          / mesh_height_min, mesh_height_max, mesh_fill, mesh_write, mesh_write_file
    NAMELIST / ElDyn           / dynamo_efield

    read_success = .FALSE.

    ! Default Parameters !

    ! SpaceManagement
    grid_file = './IPE_Grid.h5'
    NLP              = 170
    NMP              = 80
    NPTS2D           = 44514
    nFluxTube        = 1115

    ! TimeStepping !
    time_step   = 180.0_prec
    start_time  = 0.0_prec
    end_time    = 360.0_prec
    msis_time_step = 900.0_prec
    initial_timestamp = "201303160000"

    ! Forcing !
    solar_forcing_time_step = 60.0_prec
    f107_kp_size      = 1
    f107_kp_interval  = 60
    f107_kp_skip_size = 0
    f107_kp_data_size = 1
    f107_kp_read_in_start = 0
    f107_kp_file      = ''

    ! Default settings
    f107              = 120.0_prec
    f107_flag         = 1
    f107_81day_avg    = 120.0_prec
    kp                = 3.0_prec
    kp_flag           = 3
    kp_1day_avg       = 3.0_prec
    ap                = 0.0_prec
    ap_1day_avg       = 0.0_prec
    nhemi_power       = 21.0_prec
    nhemi_power_index = 6
    shemi_power       = 21.0_prec
    shemi_power_index = 6

    solarwind_angle    = 0.0_prec
    solarwind_velocity = 400.0_prec
    solarwind_Bz       = -5.0_prec
    solarwind_By       = 0.0_prec
    solarwind_density  = 3.0_prec

    ! FileIO !
    read_apex_neutrals        = .TRUE.
    read_geographic_neutrals  = .TRUE.
    write_apex_neutrals       = .TRUE.
    write_geographic_neutrals = .TRUE.
    write_geographic_eldyn    = .TRUE.
    write_apex_eldyn          = .TRUE.
    file_output_frequency     = 180.0_prec

    ! IPECAP !
    mesh_height_min =   0.
    mesh_height_max = 782.
    mesh_fill  = 1
    mesh_write = 0
    mesh_write_file = 'ipemesh'

    ! ElDyn !
    dynamo_efield          = .TRUE.

    ! Initialize buffers
    sbuf = ""
    ibuf = 0
    rbuf = 0._prec

    ! Read in namelist parameters
    IF ( mpi_layer % rank_id == 0 ) THEN

      INQUIRE( FILE = 'IPE.inp', EXIST = fileExists, IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      IF( .not. fileExists ) THEN
        OPEN( UNIT = NEWUNIT(fUnit), FILE = 'IPE.inp', ACTION = 'WRITE', DELIM = 'APOSTROPHE', IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        WRITE( UNIT = fUnit, NML = SpaceManagement, IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        WRITE( UNIT = fUnit, NML = TimeStepping,    IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        WRITE( UNIT = fUnit, NML = Forcing,         IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        WRITE( UNIT = fUnit, NML = FileIO,          IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        WRITE( UNIT = fUnit, NML = IPECAP,          IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        WRITE( UNIT = fUnit, NML = ElDyn,           IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        CLOSE( UNIT = fUnit, IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN

        PRINT*, ' '
        PRINT*, '  Module IPE_Model_Parameters_Class.F90 : S/R Build_IPE_Model_Parameters : '
        PRINT*, '    IPE.inp not found. A sample IPE.inp namelist file has been'
        PRINT*, '    generated for you in your current directory.'
      END IF

      OPEN( UNIT = NewUnit(fUnit), FILE = 'IPE.inp', ACTION = 'READ', IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      REWIND( UNIT = fUnit, IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      READ( UNIT = fUnit, NML = SpaceManagement, IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      READ( UNIT = fUnit, NML = TimeStepping,    IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      READ( UNIT = fUnit, NML = Forcing,         IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      READ( UNIT = fUnit, NML = FileIO,          IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      READ( UNIT = fUnit, NML = IPECAP,          IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      READ( UNIT = fUnit, NML = ElDyn,           IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      CLOSE( fUnit, IOSTAT = iostatus )
      IF ( iostatus /= 0 ) RETURN

      IF ( LEN_TRIM( f107_kp_file ) > 0 ) THEN
        INQUIRE( FILE = TRIM( f107_kp_file ), EXIST = params % use_f107_kp_file, IOSTAT = iostatus )
        IF ( iostatus /= 0 ) RETURN
      ENDIF

      ! prepare buffers
      ! -- strings
      sbuf = (/ grid_file, initial_timestamp, f107_kp_file, mesh_write_file /)
      ! -- integers
      ibuf(1:16) = (/ NLP, NMP, NPTS2D, nFluxTube, f107_kp_size, f107_kp_interval, &
                      f107_kp_skip_size, f107_kp_data_size, f107_kp_read_in_start, mesh_fill, mesh_write, &
                      f107_flag, kp_flag, ap_flag, nhemi_power_index, shemi_power_index /)
      ! -- logicals
      IF ( read_apex_neutrals        ) ibuf(17) = 1
      IF ( read_geographic_neutrals  ) ibuf(18) = 1
      IF ( write_apex_neutrals       ) ibuf(19) = 1
      IF ( write_geographic_neutrals ) ibuf(20) = 1
      IF ( write_geographic_eldyn    ) ibuf(21) = 1
      IF ( write_apex_eldyn          ) ibuf(22) = 1
      IF ( params % use_f107_kp_file ) ibuf(23) = 1
      IF ( dynamo_efield             ) ibuf(24) = 1

      ! -- reals
      rbuf = (/ time_step, start_time, end_time, msis_time_step, solar_forcing_time_step, &
                mesh_height_min, mesh_height_max, f107, f107_81day_avg, kp, kp_1day_avg, ap,   &
                ap_1day_avg, nhemi_power, shemi_power, solarwind_By, solarwind_angle,          &
                solarwind_velocity, solarwind_Bz, solarwind_density, file_output_frequency /)

    ENDIF

#ifdef HAVE_MPI
    CALL MPI_BCAST( sbuf, 200 * size(sbuf), MPI_CHAR, 0, mpi_layer % mpi_communicator, ierr )
#endif

    params % grid_file  = sbuf(1)
    params % initial_timestamp = sbuf(2)
    params % f107_kp_file      = sbuf(3)
    params % mesh_write_file   = sbuf(4)

#ifdef HAVE_MPI
    CALL MPI_BCAST( ibuf, size(ibuf), MPI_INTEGER, 0, mpi_layer % mpi_communicator, ierr )
#endif

    params % NLP                   = ibuf(1)
    params % NMP                   = ibuf(2)
    params % NPTS2D                = ibuf(3)
    params % nFluxTube             = ibuf(4)
    params % f107_kp_size          = ibuf(5)
    params % f107_kp_interval      = ibuf(6)
    params % f107_kp_skip_size     = ibuf(7)
    params % f107_kp_data_size     = ibuf(8)
    params % f107_kp_read_in_start = ibuf(9)
    params % mesh_fill             = ibuf(10)
    params % mesh_write            = ibuf(11)
    params % f107_flag             = ibuf(12)
    params % kp_flag               = ibuf(13)
    params % ap_flag               = ibuf(14)
    params % nhemi_power_index     = ibuf(15)
    params % shemi_power_index     = ibuf(16)
    params % read_apex_neutrals        = ( ibuf(17) == 1 )
    params % read_geographic_neutrals  = ( ibuf(18) == 1 )
    params % write_apex_neutrals       = ( ibuf(19) == 1 )
    params % write_geographic_neutrals = ( ibuf(20) == 1 )
    params % write_geographic_eldyn    = ( ibuf(21) == 1 )
    params % write_apex_eldyn          = ( ibuf(22) == 1 )
    params % use_f107_kp_file          = ( ibuf(23) == 1 )
    params % dynamo_efield             = ( ibuf(24) == 1 )

#ifdef HAVE_MPI
    CALL MPI_BCAST( rbuf, size(rbuf), mpi_layer % mpi_prec, 0, mpi_layer % mpi_communicator, ierr )
#endif

    params % time_step               = rbuf(1)
    params % start_time              = rbuf(2)
    params % end_time                = rbuf(3)
    params % msis_time_step          = rbuf(4)
    params % solar_forcing_time_step = rbuf(5)
    params % mesh_height_min         = km_to_m * rbuf(6)
    params % mesh_height_max         = km_to_m * rbuf(7)
    params % f107                    = rbuf(8)
    params % f107_81day_avg          = rbuf(9)
    params % kp                      = rbuf(10)
    params % kp_1day_avg             = rbuf(11)
    params % ap                      = rbuf(12)
    params % ap_1day_avg             = rbuf(13)
    params % nhemi_power             = rbuf(14)
    params % shemi_power             = rbuf(15)
    params % solarwind_By            = rbuf(16)
    params % solarwind_angle         = rbuf(17)
    params % solarwind_velocity      = rbuf(18)
    params % solarwind_Bz            = rbuf(19)
    params % solarwind_density       = rbuf(20)
    params % file_output_frequency   = rbuf(21)

    params % n_model_updates = INT( ( params % end_time - params % start_time ) / params % file_output_frequency )

    read_success = .TRUE.

  END SUBROUTINE Build_IPE_Model_Parameters


END MODULE IPE_Model_Parameters_Class
