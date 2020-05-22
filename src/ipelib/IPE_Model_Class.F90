#include "IPE_Macros.inc"

MODULE IPE_Model_Class

  USE IPE_Precision
  USE IPE_Model_Parameters_Class
  USE IPE_Time_Class
  USE IPE_Grid_Class
  USE IPE_Neutrals_Class
  USE IPE_Forcing_Class
  USE IPE_Plasma_Class
  USE IPE_Electrodynamics_Class
  USE IPE_MPI_Layer_Class

  USE COMIO


  IMPLICIT NONE


  ! The IPE_Model serves as a wrapper for all of the underlying attributes.
  ! This class should be used to orchestrate model setup, updating, and
  ! breakdown. At this level, we define the API for interacting with the
  ! deeper attributes within IPE.

  TYPE IPE_Model

    TYPE( IPE_Time )             :: time_tracker
    TYPE( IPE_Model_Parameters ) :: parameters
    TYPE( IPE_Grid )             :: grid
    TYPE( IPE_Forcing )          :: forcing
    TYPE( IPE_Neutrals )         :: neutrals
    TYPE( IPE_Plasma )           :: plasma
    TYPE( IPE_Electrodynamics )  :: eldyn
    TYPE( IPE_MPI_Layer )        :: mpi_layer
    CLASS( COMIO_T ), POINTER    :: io => NULL()

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model
      PROCEDURE :: Trash => Trash_IPE_Model

      PROCEDURE :: Update     => Update_IPE_Model
      PROCEDURE :: Initialize => Initialize_IPE_Model
      PROCEDURE :: Write      => Write_IPE_State
      PROCEDURE :: Read       => Read_IPE_State

  END TYPE IPE_Model


  REAL(prec), PARAMETER :: fillValue = -999999.9_prec

CONTAINS

  SUBROUTINE Build_IPE_Model( ipe, init_success, mpi_comm )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    LOGICAL,            INTENT(out)   :: init_success
    INTEGER, OPTIONAL,  INTENT(in)    :: mpi_comm

    ! Local
    INTEGER        :: error
    CHARACTER(200) :: init_file
    LOGICAL        :: fileExists

!   INTEGER, DIMENSION(3) :: global_dims, local_dims, local_start

    init_success = .false.

    CALL ipe % mpi_layer % Initialize( comm = mpi_comm )

    CALL ipe % parameters % Build( ipe % mpi_layer, init_success )
    IF ( .NOT. init_success ) RETURN

    init_success = .false.

    CALL ipe % mpi_layer % Set_Domain( ipe % parameters % NLP, &
                                       ipe % parameters % NMP, error )
    IF ( error /= 0 ) RETURN

    ! Initialize I/O
    IF ( ipe % mpi_layer % enabled ) THEN
      ipe % io => COMIO_T(fmt=COMIO_FMT_HDF5, &
                          comm=ipe % mpi_layer % mpi_communicator, &
                          info=ipe % mpi_layer % mpi_info)
      IF (ipe % io % err % check(msg="Failed to initialize I/O layer", &
        file=__FILE__, line=__LINE__)) RETURN
    ELSE
      ipe % io => COMIO_T(fmt=COMIO_FMT_HDF5)
      IF (ipe % io % err % check(msg="Failed to initialize I/O layer", &
        file=__FILE__, line=__LINE__)) RETURN
    END IF

    ! Initialize the clock
    CALL ipe % time_tracker % Build( ipe % parameters % initial_timestamp, ipe % mpi_layer % rank_id )

    ! ////// Forcing ////// !

    CALL ipe % forcing % Build( ipe % parameters, error )
    IF ( error /= 0 ) RETURN

    ! ////// grid ////// !


    CALL ipe % grid % Create( ipe % mpi_layer, ipe % parameters )

    CALL ipe % grid % ReadFile( ipe % io, "IPE_Grid.h5", error )
    IF ( error /= 0 ) RETURN

    CALL ipe % grid % Initialize( ipe % parameters, error )
    IF ( error /= 0 ) RETURN

    ! ////// neutrals ////// !

    CALL ipe % neutrals % Build( nFluxtube       = ipe % grid % nFluxTube, &
                                 NLP             = ipe % grid % NLP, &
                                 NMP             = ipe % grid % NMP, &
                                 mp_low          = ipe % mpi_layer % mp_low, &
                                 mp_high         = ipe % mpi_layer % mp_high )

    ! ////// electric field /////// !

    CALL ipe % eldyn % Build( nFluxTube = ipe % grid % nFluxTube,&
                              NLP       = ipe % grid % NLP, &
                              NMP       = ipe % grid % NMP, &
                              dynamo    = ipe % parameters % dynamo_efield, &
                              mp_low    = ipe % mpi_layer % mp_low, &
                              mp_high   = ipe % mpi_layer % mp_high, &
                              halo      = ipe % mpi_layer % mp_halo_size )


    ! ////// plasma ////// !

    CALL ipe % plasma % Build( nFluxTube = ipe % grid % nFluxTube, &
                               NLP       = ipe % grid % NLP, &
                               NMP       = ipe % grid % NMP, &
                               mp_low    = ipe % mpi_layer % mp_low, &
                               mp_high   = ipe % mpi_layer % mp_high, &
                               halo      = ipe % mpi_layer % mp_halo_size )


    ! Setup data decomposition for common disk I/O
!   global_dims = (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % grid % NMP /)
!   local_dims  = (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /)
!   local_start = (/ 1, 1, ipe % mpi_layer % mp_low /)

!   CALL ipe % io % domain(global_dims, local_start, local_dims)
!   IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
!     file=__FILE__, line=__LINE__)) RETURN

    init_success = .true.

  END SUBROUTINE Build_IPE_Model


  SUBROUTINE Trash_IPE_Model( ipe )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe

    CALL ipe % forcing  % Trash( )
    CALL ipe % grid     % Trash( )
    CALL ipe % neutrals % Trash( )
    CALL ipe % plasma   % Trash( )
    CALL ipe % eldyn    % Trash( )

    CALL ipe % io       % Shutdown( )

    CALL ipe % mpi_layer % Finalize( )

  END SUBROUTINE Trash_IPE_Model


  SUBROUTINE Update_IPE_Model( ipe, t0, t1 )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    REAL(prec), OPTIONAL, INTENT(in)  :: t0, t1

    ! Local
    INTEGER    :: i, nSteps

    IF( PRESENT( t0 ) .AND. PRESENT(t1) )THEN

      nSteps = INT((t1-t0)/ipe % parameters % time_step)

    ELSE

      nSteps = 1

    ENDIF

    DO i = 1, nSteps

      call ipe % forcing % Update_Current_Index( ipe % parameters, ipe % time_tracker % elapsed_sec )

      CALL ipe % neutrals % Update( ipe % parameters, ipe % grid, ipe % time_tracker, ipe % forcing, ipe % mpi_layer )

      CALL ipe % eldyn % Update( ipe % grid, &
                                 ipe % forcing, &
                                 ipe % time_tracker, &
                                 ipe % plasma, &
                                 ipe % mpi_layer )

      CALL ipe % plasma % Update( ipe % grid, &
                                  ipe % neutrals, &
                                  ipe % forcing, &
                                  ipe % time_tracker, &
                                  ipe % mpi_layer, &
                                  ipe % eldyn % v_ExB_apex, &
                                  ipe % parameters % time_step )

      ! Update the timer
      CALL ipe % time_tracker % Increment( ipe % parameters % time_step )

    ENDDO

  END SUBROUTINE Update_IPE_Model

  SUBROUTINE Initialize_IPE_Model( ipe, filename, success )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    LOGICAL,            INTENT(out)   :: success

    ! Local
    INTEGER    :: error

    success = .false.

    CALL ipe % Read( filename, error )
    IF ( error /= 0 ) RETURN

    CALL ipe % plasma % Calculate_Field_Line_Integrals( ipe % grid, ipe % neutrals, ipe % mpi_layer )

    success = .true.

  END SUBROUTINE Initialize_IPE_Model


  SUBROUTINE Read_IPE_State( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    INTEGER, PARAMETER :: num_ion_densities = 9
    CHARACTER(LEN=*), DIMENSION(num_ion_densities), PARAMETER :: ion_densities = &
      (/ &
        "/apex/o_plus_density   ", &
        "/apex/h_plus_density   ", &
        "/apex/he_plus_density  ", &
        "/apex/n_plus_density   ", &
        "/apex/no_plus_density  ", &
        "/apex/o2_plus_density  ", &
        "/apex/n2_plus_density  ", &
        "/apex/o_plus_2D_density", &
        "/apex/o_plus_2P_density"  &
      /)

    INTEGER, PARAMETER :: num_ion_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_ion_velocities), PARAMETER :: ion_velocities = &
      (/ &
        "/apex/o_plus_velocity ", &
        "/apex/h_plus_velocity ", &
        "/apex/he_plus_velocity"  &
      /)

    INTEGER, PARAMETER :: num_plasma_datasets = 2
    CHARACTER(LEN=*), DIMENSION(num_plasma_datasets), PARAMETER :: plasma_datasets = &
      (/ &
        "/apex/ion_temperature     ", &
        "/apex/electron_temperature"  &
      /)

    INTEGER, PARAMETER :: num_apex_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_apex_velocities), PARAMETER :: apex_velocities = &
      (/ &
        "/apex/neutral_apex1_velocity", &
        "/apex/neutral_apex2_velocity", &
        "/apex/neutral_apex3_velocity"  &
      /)


    INTEGER, PARAMETER :: num_geo_datasets = 3
    CHARACTER(LEN=*), DIMENSION(num_geo_datasets), PARAMETER :: geo_datasets = &
      (/ &
        "/apex/neutral_geographic_velocity1", &
        "/apex/neutral_geographic_velocity2", &
        "/apex/neutral_geographic_velocity3"  &
      /)

    INTEGER :: item

    ! Begin
    error = -1

    IF( ipe % mpi_layer % rank_id == 0 )THEN
      PRINT *, '  Reading output file : '//TRIM(filename)
    ENDIF

    ! Open HDF5 input file
    CALL ipe % io % open(filename, "r")
    IF (ipe % io % err % check(msg="Failed to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Setup common data decomposition for datasets
    CALL ipe % io % domain( (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % grid % NMP /), &
      (/ 1, 1, ipe % mpi_layer % mp_low /), &
      (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /) )
    IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
      file=__FILE__, line=__LINE__)) RETURN

    ! Read individual datasets
    ! -- Ion densities
    DO item = 1, num_ion_densities
      CALL ipe % io % read(ion_densities(item), &
        ipe % plasma % ion_densities(item,:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset "//ion_densities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Ion velocities
    DO item = 1, num_ion_velocities
      CALL ipe % io % read(ion_velocities(item), &
        ipe % plasma % ion_velocities(item,:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset "//ion_velocities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Plasma temperatures
    DO item = 1, num_plasma_datasets
      SELECT CASE (item)
      CASE (1)
        CALL ipe % io % read(plasma_datasets(item), &
          ipe % plasma % ion_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      CASE (2)
        CALL ipe % io % read(plasma_datasets(item), &
          ipe % plasma % electron_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END SELECT
    END DO

    ! -- Compute plasma electron density from ion densities
    ipe % plasma % electron_density2(:,:,ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high) = &
      SUM(ipe % plasma % ion_densities(1:9, :, :, ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high), dim=1)

    ! -- Read neutral datasets if requested
    IF( ipe % parameters % read_apex_neutrals )THEN
      CALL ipe % io % read("/apex/o_density", &
        ipe % neutrals % oxygen(:,:,                &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/o_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/h_density", &
        ipe % neutrals % hydrogen(:,:,              &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/h_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/he_density", &
        ipe % neutrals % helium(:,:,                 &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/he_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/n_density", &
        ipe % neutrals % nitrogen(:,:,              &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/n_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/o2_density", &
        ipe % neutrals % molecular_oxygen(:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/o2_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/n2_density", &
        ipe % neutrals % molecular_nitrogen(:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/n2_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/neutral_temperature", &
        ipe % neutrals % temperature(:,:,                     &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/neutral_temperature", &
        file=__FILE__, line=__LINE__)) RETURN

      ! Read neutral velocities on apex grid
      DO item = 1, num_apex_velocities
        CALL ipe % io % read(apex_velocities(item), &
          ipe % neutrals % velocity_apex(item,:,:,         &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//apex_velocities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Read neutral velocities on geographic grid, if requested
    IF( ipe % parameters % read_geographic_neutrals )THEN
      DO item = 1, num_geo_datasets
        CALL ipe % io % read(geo_datasets(item), &
          ipe % neutrals % velocity_geographic(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//geo_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Close HDF5 file
    CALL ipe % io % close()
    IF (ipe % io % err % check(msg="Failed to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    error = 0

  END SUBROUTINE Read_IPE_State

  SUBROUTINE Write_IPE_State( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    INTEGER, PARAMETER :: num_groups = 1
    CHARACTER(LEN=*), DIMENSION(num_groups),   PARAMETER :: groups = (/ "apex" /)

    INTEGER, PARAMETER :: num_ion_densities = 9
    CHARACTER(LEN=*), DIMENSION(num_ion_densities), PARAMETER :: ion_densities = &
      (/ &
        "/apex/o_plus_density   ", &
        "/apex/h_plus_density   ", &
        "/apex/he_plus_density  ", &
        "/apex/n_plus_density   ", &
        "/apex/no_plus_density  ", &
        "/apex/o2_plus_density  ", &
        "/apex/n2_plus_density  ", &
        "/apex/o_plus_2D_density", &
        "/apex/o_plus_2P_density"  &
      /)

    INTEGER, PARAMETER :: num_ion_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_ion_velocities), PARAMETER :: ion_velocities = &
      (/ &
        "/apex/o_plus_velocity ", &
        "/apex/h_plus_velocity ", &
        "/apex/he_plus_velocity"  &
      /)

    INTEGER, PARAMETER :: num_plasma_datasets = 2
    CHARACTER(LEN=*), DIMENSION(num_plasma_datasets), PARAMETER :: plasma_datasets = &
      (/ &
        "/apex/ion_temperature     ", &
        "/apex/electron_temperature"  &
      /)

    INTEGER, PARAMETER :: num_apex_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_apex_velocities), PARAMETER :: apex_velocities = &
      (/ &
        "/apex/neutral_apex1_velocity", &
        "/apex/neutral_apex2_velocity", &
        "/apex/neutral_apex3_velocity"  &
      /)

    INTEGER, PARAMETER :: num_geo_datasets = 3
    CHARACTER(LEN=*), DIMENSION(num_geo_datasets), PARAMETER :: geo_datasets = &
      (/ &
        "/apex/neutral_geographic_velocity1", &
        "/apex/neutral_geographic_velocity2", &
        "/apex/neutral_geographic_velocity3"  &
      /)

    INTEGER :: item
    CHARACTER(LEN=28) :: dset_name

    ! Begin

    error = -1

    IF( ipe % mpi_layer % rank_id == 0 )THEN
      PRINT *, '  Writing output file : '//TRIM(filename)
    ENDIF

    ! Create HDF5 input file
    CALL ipe % io % open(filename, "c")
    IF (ipe % io % err % check(msg="Unable to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Setup common data decomposition for datasets
    CALL ipe % io % domain( (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % grid % NMP /), &
      (/ 1, 1, ipe % mpi_layer % mp_low /), &
      (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /) )
    IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
      file=__FILE__, line=__LINE__)) RETURN
    ! Create groups
!   DO item = 1, num_groups
!     CALL ipe % io % grp_build(groups(item))
!     IF (ipe % io % err % check(msg="Unable to create group "//groups(item), &
!       file=__FILE__, line=__LINE__)) RETURN
!   END DO

    ! Write individual datasets
    ! -- Ion densities
    DO item = 1, num_ion_densities
      CALL ipe % io % write(ion_densities(item), &
        ipe % plasma % ion_densities(item,:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//ion_densities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Ion velocities
    DO item = 1, num_ion_velocities
      CALL ipe % io % write(ion_velocities(item), &
        ipe % plasma % ion_velocities(item,:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//ion_velocities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Plasma temperatures
    DO item = 1, num_plasma_datasets
      SELECT CASE (item)
      CASE (1)
        CALL ipe % io % write(plasma_datasets(item), &
          ipe % plasma % ion_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      CASE (2)
        CALL ipe % io % write(plasma_datasets(item), &
          ipe % plasma % electron_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END SELECT
    END DO

    ! -- Read neutral datasets if requested
    IF( ipe % parameters % write_apex_neutrals )THEN

      dset_name = "/apex/o_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % oxygen(:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/h_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % hydrogen(:,:,      &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/he_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % helium(:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/n_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % nitrogen(:,:,      &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/o2_density"
      CALL ipe % io % write(dset_name,    &
        ipe % neutrals % molecular_oxygen(:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/n2_density"
      CALL ipe % io % write(dset_name,      &
        ipe % neutrals % molecular_nitrogen(:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/neutral_temperature"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % temperature(:,:,   &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      ! Read neutral velocities on apex grid
      DO item = 1, num_apex_velocities
        CALL ipe % io % write(apex_velocities(item), &
          ipe % neutrals % velocity_apex(item,:,:,         &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//apex_velocities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Read neutral velocities on geographic grid, if requested
    IF( ipe % parameters % write_geographic_neutrals )THEN
      DO item = 1, num_geo_datasets
        CALL ipe % io % write(geo_datasets(item), &
          ipe % neutrals % velocity_geographic(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//geo_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Close HDF5 file
    CALL ipe % io % close()
    IF (ipe % io % err % check(msg="Unable to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    error = 0

  END SUBROUTINE Write_IPE_State


END MODULE IPE_Model_Class
