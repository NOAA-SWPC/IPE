MODULE IPE_Model_Parameters_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Common_Routines
  USE IPE_MPI_Layer_Class
  USE ipe_error_module


  IMPLICIT NONE

  TYPE IPE_Model_Parameters

    ! SpaceManagement
    CHARACTER(200) :: grid_file

    ! TimeStepping
    REAL(prec)    :: time_step
    REAL(prec)    :: start_time
    REAL(prec)    :: end_time
    REAL(prec)    :: msis_time_step
    CHARACTER(12) :: initial_timestamp

    ! Forcing
    REAL(prec)     :: solar_forcing_time_step
    INTEGER        :: ifp_realtime_interval
    LOGICAL        :: use_ifp_file
    CHARACTER(200) :: ifp_file

    ! IPECAP
    REAL(prec) :: mesh_height_min
    REAL(prec) :: mesh_height_max
    INTEGER    :: mesh_fill
    INTEGER    :: mesh_write
    CHARACTER(200) :: mesh_write_file

    ! >> Fixed parameters
    REAL(prec) :: f107
    REAL(prec) :: f107_81day_avg
    REAL(prec) :: kp
    REAL(prec) :: kp_1day_avg
    REAL(prec) :: ap
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

    ! >> Operations
    REAL(prec) :: colfac
    REAL(prec) :: offset1_deg
    REAL(prec) :: offset2_deg
    INTEGER    :: potential_model
    REAL(prec) :: hpeq
    INTEGER    :: transport_highlat_lp
    INTEGER    :: perp_transport_max_lp  
    REAL(prec) :: vertical_wind_limit

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model_Parameters

  END TYPE IPE_Model_Parameters


CONTAINS

  SUBROUTINE Build_IPE_Model_Parameters( parameters, mpi_layer, rc )

    IMPLICIT NONE

    CLASS( IPE_Model_Parameters ), INTENT(out) :: parameters
    TYPE( IPE_MPI_Layer ),         INTENT(in)  :: mpi_layer
    INTEGER, OPTIONAL,             INTENT(out) :: rc

    ! Local
    CHARACTER(200) :: grid_file
    INTEGER        :: fUnit, ierr, iostatus
    LOGICAL        :: fileExists
    REAL(prec)     :: solar_forcing_time_step
    REAL(prec)     :: time_step, start_time, end_time, msis_time_step
    CHARACTER(12)  :: initial_timestamp
    INTEGER        :: ifp_realtime_interval
    CHARACTER(200) :: ifp_file
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
    REAL(prec) :: kp
    REAL(prec) :: kp_1day_avg
    REAL(prec) :: ap
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
    ! >> Ops
    REAL(prec) :: colfac
    REAL(prec) :: offset1_deg
    REAL(prec) :: offset2_deg
    INTEGER    :: potential_model
    REAL(prec) :: hpeq
    INTEGER    :: transport_highlat_lp
    INTEGER    :: perp_transport_max_lp  
    REAL(prec) :: vertical_wind_limit

    ! Communication buffers
    CHARACTER(LEN=200), DIMENSION( 4) :: sbuf
    INTEGER,            DIMENSION(16) :: ibuf
    REAL(prec),         DIMENSION(26) :: rbuf


    NAMELIST / SpaceManagement / grid_file
    NAMELIST / TimeStepping    / time_step, start_time, end_time, msis_time_step, initial_timestamp
    NAMELIST / Forcing         / solar_forcing_time_step, ifp_realtime_interval, &
                                 ifp_file, f107, f107_81day_avg, &
                                 kp, kp_1day_avg, ap, ap_1day_avg, nhemi_power, &
                                 nhemi_power_index, shemi_power, shemi_power_index, solarwind_By, solarwind_angle, &
                                 solarwind_velocity, solarwind_Bz, solarwind_density
    NAMELIST / FileIO          / read_apex_neutrals, read_geographic_neutrals, write_apex_neutrals, write_geographic_neutrals, &
                                 write_geographic_eldyn, write_apex_eldyn, file_output_frequency
    NAMELIST / IPECAP          / mesh_height_min, mesh_height_max, mesh_fill, mesh_write, mesh_write_file
    NAMELIST / ElDyn           / dynamo_efield
    NAMELIST / OPERATIONAL     / colfac, offset1_deg, offset2_deg, potential_model, hpeq, &
                                 transport_highlat_lp, perp_transport_max_lp, vertical_wind_limit                       

    ! Begin
    IF (PRESENT(rc)) rc = IPE_SUCCESS

    ! Default Parameters !

    ! SpaceManagement
    grid_file = './IPE_Grid.h5'

    ! TimeStepping !
    time_step   = 180.0_prec
    start_time  = 0.0_prec
    end_time    = 360.0_prec
    msis_time_step = 900.0_prec
    initial_timestamp = "201303160000"

    ! Forcing !
    solar_forcing_time_step = 60.0_prec
    ifp_realtime_interval = -1
    ifp_file      = ''

    ! Default settings
    f107              = 120.0_prec
    f107_81day_avg    = 120.0_prec
    kp                = 3.0_prec
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
    dynamo_efield       = .TRUE.

    ! Operational
    colfac               = 1.3_prec
    offset1_deg          = 5.0_prec
    offset2_deg          = 20.0_prec
    potential_model      = 2
    hpeq                 = 0.0_prec
    transport_highlat_lp = 30
    perp_transport_max_lp   = 151
    vertical_wind_limit  = 100.0_prec

    ! Initialize buffers
    sbuf = ""
    ibuf = 0
    rbuf = 0._prec

    ! Read in namelist parameters
    IF ( mpi_layer % rank_id == 0 ) THEN

      INQUIRE( FILE = 'IPE.inp', EXIST = fileExists, IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, msg="Failed to retrieve IPE.inp status", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
      IF ( ipe_status_check( fileExists, msg="Input file IPE.inp not found", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      OPEN( UNIT = NewUnit(fUnit), FILE = 'IPE.inp', ACTION = 'READ', IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, msg="Unable to open file IPE.inp for reading", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      REWIND( UNIT = fUnit, IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = SpaceManagement, IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = TimeStepping,    IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = Forcing,         IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = FileIO,          IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = IPECAP,          IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = ElDyn,           IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ( UNIT = fUnit, NML = OPERATIONAL,     IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      CLOSE( fUnit, IOSTAT = iostatus )
      IF ( ipe_iostatus_check( iostatus, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      IF ( LEN_TRIM( ifp_file ) > 0 ) THEN
        INQUIRE( FILE = TRIM( ifp_file ), EXIST = parameters % use_ifp_file, IOSTAT = iostatus )
        IF ( ipe_iostatus_check( iostatus, msg="Failed to inquire file "//ifp_file, &
          line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
      ENDIF

      ! prepare buffers
      ! -- strings
      sbuf = (/ grid_file, initial_timestamp, ifp_file, mesh_write_file /)
      ! -- integers
      ibuf(1:5) = (/ ifp_realtime_interval, mesh_fill, mesh_write, &
                      nhemi_power_index, shemi_power_index /)
      ! -- logicals
      IF ( read_apex_neutrals        ) ibuf(6) = 1
      IF ( read_geographic_neutrals  ) ibuf(7) = 1
      IF ( write_apex_neutrals       ) ibuf(8) = 1
      IF ( write_geographic_neutrals ) ibuf(9) = 1
      IF ( write_geographic_eldyn    ) ibuf(10) = 1
      IF ( write_apex_eldyn          ) ibuf(11) = 1
      IF ( parameters % use_ifp_file     ) ibuf(12) = 1
      IF ( dynamo_efield             ) ibuf(13) = 1
      ! -- integers for operations
      ibuf(14) = potential_model
      ibuf(15) = transport_highlat_lp
      ibuf(16) = perp_transport_max_lp

      ! -- reals
      rbuf = (/ time_step, start_time, end_time, msis_time_step, solar_forcing_time_step, &
                mesh_height_min, mesh_height_max, f107, f107_81day_avg, kp, kp_1day_avg, ap,   &
                ap_1day_avg, nhemi_power, shemi_power, solarwind_By, solarwind_angle,          &
                solarwind_velocity, solarwind_Bz, solarwind_density, file_output_frequency, &
                colfac, offset1_deg, offset2_deg, hpeq, vertical_wind_limit /)

    ENDIF

#ifdef HAVE_MPI
    CALL MPI_BCAST( sbuf, 200 * size(sbuf), MPI_CHAR, 0, mpi_layer % mpi_communicator, ierr )
    IF ( ipe_status_check( ierr == 0, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
#endif

    parameters % grid_file         = sbuf(1)
    parameters % initial_timestamp = sbuf(2)
    parameters % ifp_file          = sbuf(3)
    parameters % mesh_write_file   = sbuf(4)

#ifdef HAVE_MPI
    CALL MPI_BCAST( ibuf, size(ibuf), MPI_INTEGER, 0, mpi_layer % mpi_communicator, ierr )
    IF ( ipe_status_check( ierr == 0, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
#endif

    parameters % ifp_realtime_interval     = ibuf(1)
    parameters % mesh_fill                 = ibuf(2)
    parameters % mesh_write                = ibuf(3)
    parameters % nhemi_power_index         = ibuf(4)
    parameters % shemi_power_index         = ibuf(5)
    parameters % read_apex_neutrals        = ( ibuf(6) == 1 )
    parameters % read_geographic_neutrals  = ( ibuf(7) == 1 )
    parameters % write_apex_neutrals       = ( ibuf(8) == 1 )
    parameters % write_geographic_neutrals = ( ibuf(9) == 1 )
    parameters % write_geographic_eldyn    = ( ibuf(10) == 1 )
    parameters % write_apex_eldyn          = ( ibuf(11) == 1 )
    parameters % use_ifp_file              = ( ibuf(12) == 1 )
    parameters % dynamo_efield             = ( ibuf(13) == 1 )
    parameters % potential_model           = ibuf(14)
    parameters % transport_highlat_lp      = ibuf(15)
    parameters % perp_transport_max_lp     = ibuf(16)

#ifdef HAVE_MPI
    CALL MPI_BCAST( rbuf, size(rbuf), mpi_layer % mpi_prec, 0, mpi_layer % mpi_communicator, ierr )
    IF ( ipe_status_check( ierr == 0, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
#endif

    parameters % time_step               = rbuf(1)
    parameters % start_time              = rbuf(2)
    parameters % end_time                = rbuf(3)
    parameters % msis_time_step          = rbuf(4)
    parameters % solar_forcing_time_step = rbuf(5)
    parameters % mesh_height_min         = km_to_m * rbuf(6)
    parameters % mesh_height_max         = km_to_m * rbuf(7)
    parameters % f107                    = rbuf(8)
    parameters % f107_81day_avg          = rbuf(9)
    parameters % kp                      = rbuf(10)
    parameters % kp_1day_avg             = rbuf(11)
    parameters % ap                      = rbuf(12)
    parameters % ap_1day_avg             = rbuf(13)
    parameters % nhemi_power             = rbuf(14)
    parameters % shemi_power             = rbuf(15)
    parameters % solarwind_By            = rbuf(16)
    parameters % solarwind_angle         = rbuf(17)
    parameters % solarwind_velocity      = rbuf(18)
    parameters % solarwind_Bz            = rbuf(19)
    parameters % solarwind_density       = rbuf(20)
    parameters % file_output_frequency   = rbuf(21)
    parameters % colfac                  = rbuf(22)
    parameters % offset1_deg             = rbuf(23)
    parameters % offset2_deg             = rbuf(24)
    parameters % hpeq                    = rbuf(25)
    parameters % vertical_wind_limit     = rbuf(26)

    parameters % n_model_updates = INT( ( parameters % end_time - parameters % start_time ) / parameters % file_output_frequency )

  END SUBROUTINE Build_IPE_Model_Parameters


END MODULE IPE_Model_Parameters_Class
