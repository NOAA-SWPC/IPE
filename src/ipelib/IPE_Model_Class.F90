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
USE IPE_IO_Class

#ifdef HAVE_NETCDF
USE netcdf
#endif
USE HDF5


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
    TYPE( IPE_IO )               :: io

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model
      PROCEDURE :: Trash => Trash_IPE_Model

      PROCEDURE :: Update => Update_IPE_Model
      PROCEDURE :: Initialize => Initialize_IPE_Model

      PROCEDURE :: Write_Geographic_NetCDF_IPE
      PROCEDURE :: Write_NMF2_TEC_NetCDF_IPE
      PROCEDURE :: Write_Electron_Density_NetCDF_IPE
      PROCEDURE :: Write_Electron_Density_NetCDF_IPE_integer
      PROCEDURE :: Write_Reduced_Geographic_NetCDF_IPE
      PROCEDURE :: Write_to_HDF5
      PROCEDURE :: Read_from_HDF5
      PROCEDURE :: Read_from_Electron_Density_HDF5

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: Geographic_Interpolation

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

    INTEGER, DIMENSION(3) :: global_dims, local_dims, local_offsets

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
      CALL ipe % io % build(comm=ipe % mpi_layer % mpi_communicator, &
                            info=ipe % mpi_layer % mpi_info)
      IF (ipe % io % error_check(msg="Failed to initialize I/O layer", &
        file=__FILE__, line=__LINE__)) RETURN
    ELSE
      CALL ipe % io % build()
      IF (ipe % io % error_check(msg="Failed to initialize I/O layer", &
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
    global_dims = (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % grid % NMP /)
    local_dims  = (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /)
    local_offsets = (/ 0, 0, ipe % mpi_layer % mp_low - 1 /)

    CALL ipe % io % domain_set(global_dims, local_offsets, local_dims)
    IF (ipe % io % error_check(msg="Failed to setup I/O data decomposition", &
      file=__FILE__, line=__LINE__)) RETURN

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

    CALL ipe % io       % Trash( )

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


  SUBROUTINE Geographic_Interpolation( ipe )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe

    IF( ipe % parameters % write_geographic_eldyn )THEN
      CALL ipe % eldyn % Interpolate_to_GeographicGrid( ipe % grid )
    ENDIF

    CALL ipe % plasma % Interpolate_to_GeographicGrid( ipe % grid )
    CALL ipe % neutrals % Interpolate_to_GeographicGrid( ipe % grid )

  END SUBROUTINE Geographic_Interpolation


  SUBROUTINE Initialize_IPE_Model( ipe, filename, success )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    LOGICAL,            INTENT(out)   :: success

    ! Local
    INTEGER    :: error

    success = .false.

    CALL ipe % Read_from_HDF5( filename, error )
    IF ( error /= 0 ) RETURN

    CALL ipe % plasma % Calculate_Field_Line_Integrals( ipe % grid, ipe % neutrals, ipe % mpi_layer )

    success = .true.

  END SUBROUTINE Initialize_IPE_Model


  SUBROUTINE Read_from_HDF5( ipe, filename, error )

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
    CALL ipe % io % file_open(filename, "r")
    IF (ipe % io % error_check(msg="Failed to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Read individual datasets
    ! -- Ion densities
    DO item = 1, num_ion_densities
      CALL ipe % io % dset_read(ion_densities(item), &
        ipe % plasma % ion_densities(item,:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset "//ion_densities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Ion velocities
    DO item = 1, num_ion_velocities
      CALL ipe % io % dset_read(ion_velocities(item), &
        ipe % plasma % ion_velocities(item,:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset "//ion_velocities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Plasma temperatures
    DO item = 1, num_plasma_datasets
      SELECT CASE (item)
      CASE (1)
        CALL ipe % io % dset_read(plasma_datasets(item), &
          ipe % plasma % ion_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Failed to read dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      CASE (2)
        CALL ipe % io % dset_read(plasma_datasets(item), &
          ipe % plasma % electron_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Failed to read dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END SELECT
    END DO

    ! -- Compute plasma electron density from ion densities
    ipe % plasma % electron_density2(:,:,ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high) = &
      SUM(ipe % plasma % ion_densities(1:9, :, :, ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high), dim=1)

    ! -- Read neutral datasets if requested
    IF( ipe % parameters % read_apex_neutrals )THEN
      CALL ipe % io % dset_read("/apex/o_density", &
        ipe % neutrals % oxygen(:,:,                &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/o_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % dset_read("/apex/h_density", &
        ipe % neutrals % hydrogen(:,:,              &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/h_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % dset_read("/apex/he_density", &
        ipe % neutrals % helium(:,:,                 &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/he_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % dset_read("/apex/n_density", &
        ipe % neutrals % nitrogen(:,:,              &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/n_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % dset_read("/apex/o2_density", &
        ipe % neutrals % molecular_oxygen(:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/o2_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % dset_read("/apex/n2_density", &
        ipe % neutrals % molecular_nitrogen(:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/n2_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % dset_read("/apex/neutral_temperature", &
        ipe % neutrals % temperature(:,:,                     &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Failed to read dataset /apex/neutral_temperature", &
        file=__FILE__, line=__LINE__)) RETURN

      ! Read neutral velocities on apex grid
      DO item = 1, num_apex_velocities
        CALL ipe % io % dset_read(apex_velocities(item), &
          ipe % neutrals % velocity_apex(item,:,:,         &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Failed to read dataset "//apex_velocities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Read neutral velocities on geographic grid, if requested
    IF( ipe % parameters % read_geographic_neutrals )THEN
      DO item = 1, num_geo_datasets
        CALL ipe % io % dset_read(geo_datasets(item), &
          ipe % neutrals % velocity_geographic(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Failed to read dataset "//geo_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Close HDF5 file
    CALL ipe % io % file_close()
    IF (ipe % io % error_check(msg="Failed to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    error = 0

  END SUBROUTINE Read_from_HDF5

  SUBROUTINE Write_to_HDF5( ipe, filename, error )

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
    CALL ipe % io % file_open(filename, "c")
    IF (ipe % io % error_check(msg="Unable to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Create groups
    DO item = 1, num_groups
      CALL ipe % io % grp_build(groups(item))
      IF (ipe % io % error_check(msg="Unable to create group "//groups(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO

    ! Write individual datasets
    ! -- Ion densities
    DO item = 1, num_ion_densities
      CALL ipe % io % dset_write(ion_densities(item), &
        ipe % plasma % ion_densities(item,:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//ion_densities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Ion velocities
    DO item = 1, num_ion_velocities
      CALL ipe % io % dset_write(ion_velocities(item), &
        ipe % plasma % ion_velocities(item,:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//ion_velocities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Plasma temperatures
    DO item = 1, num_plasma_datasets
      SELECT CASE (item)
      CASE (1)
        CALL ipe % io % dset_write(plasma_datasets(item), &
          ipe % plasma % ion_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Unable to write dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      CASE (2)
        CALL ipe % io % dset_write(plasma_datasets(item), &
          ipe % plasma % electron_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Unable to write dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END SELECT
    END DO

    ! -- Read neutral datasets if requested
    IF( ipe % parameters % write_apex_neutrals )THEN

      dset_name = "/apex/o_density"
      CALL ipe % io % dset_write(dset_name, &
        ipe % neutrals % oxygen(:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/h_density"
      CALL ipe % io % dset_write(dset_name, &
        ipe % neutrals % hydrogen(:,:,      &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/he_density"
      CALL ipe % io % dset_write(dset_name, &
        ipe % neutrals % helium(:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/n_density"
      CALL ipe % io % dset_write(dset_name, &
        ipe % neutrals % nitrogen(:,:,      &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/o2_density"
      CALL ipe % io % dset_write(dset_name,    &
        ipe % neutrals % molecular_oxygen(:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/n2_density"
      CALL ipe % io % dset_write(dset_name,      &
        ipe % neutrals % molecular_nitrogen(:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/neutral_temperature"
      CALL ipe % io % dset_write(dset_name, &
        ipe % neutrals % temperature(:,:,   &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % error_check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      ! Read neutral velocities on apex grid
      DO item = 1, num_apex_velocities
        CALL ipe % io % dset_write(apex_velocities(item), &
          ipe % neutrals % velocity_apex(item,:,:,         &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Unable to write dataset "//apex_velocities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Read neutral velocities on geographic grid, if requested
    IF( ipe % parameters % write_geographic_neutrals )THEN
      DO item = 1, num_geo_datasets
        CALL ipe % io % dset_write(geo_datasets(item), &
          ipe % neutrals % velocity_geographic(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % error_check(msg="Unable to write dataset "//geo_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Close HDF5 file
    CALL ipe % io % file_close()
    IF (ipe % io % error_check(msg="Unable to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    error = 0

  END SUBROUTINE Write_to_HDF5

  SUBROUTINE Read_from_Electron_Density_HDF5( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local

    ! Begin
    error = -1

    IF( ipe % mpi_layer % rank_id == 0 )THEN
      PRINT*, '  Reading electron density file : '//TRIM(filename)
    ENDIF

    ! Open HDF5 input file
    CALL ipe % io % file_open(filename, "r")
    IF (ipe % io % error_check(msg="Unable to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Read dataset
    CALL ipe % io % dset_read("/apex/electron_density", &
      ipe % plasma % electron_density2(:,:,             &
      ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
    IF (ipe % io % error_check(msg="Unable to read dataset /apex/electron_density", &
      file=__FILE__, line=__LINE__)) RETURN

    ! -- Close file
    CALL ipe % io % file_close()
    IF (ipe % io % error_check(msg="Unable to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    error = 0

  END SUBROUTINE Read_from_Electron_Density_HDF5


  SUBROUTINE Write_NetCDF_IPE( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(in)  :: ipe
    CHARACTER(*),       INTENT(in)  :: filename
    INTEGER,            INTENT(out) :: error

    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid, time_varid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid, ve1_varid, ve2_varid
    INTEGER :: op_varid, hp_varid, hep_varid, np_varid, n2p_varid, o2p_varid, nop_varid
    INTEGER :: op2d_varid, op2p_varid, ion_temp_varid, phi_varid, mhd_phi_varid, hc_varid, pc_varid, bc_varid, electron_temp_varid
    INTEGER :: opv_varid(1:3), hpv_varid(1:3), hepv_varid(1:3), npv_varid(1:3), n2pv_varid(1:3), o2pv_varid(1:3), nopv_varid(1:3)
    INTEGER :: ion_vel_op_varid, ion_vel_hp_varid, ion_vel_hep_varid
    INTEGER :: recStart(1:4), recCount(1:4)


    error = 0

#ifdef HAVE_NETCDF
    recStart = (/ 1, 1, 1, 1 /)
    recCount = (/ ipe % grid % nFluxTube, ipe % grid % NLP, ipe % grid % NMP, 1 /)


    time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    IF( prec == sp )THEN
      NF90_PREC = NF90_FLOAT
    ELSE
      NF90_PREC = NF90_DOUBLE
    ENDIF

    _ncCheck( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
    _ncCheck( nf90_put_att( ncid, NF90_GLOBAL, "Version", "1.1") )

    _ncCheck( nf90_def_dim( ncid, "s", ipe % grid % nFluxTube, z_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "lp", ipe % grid % NLP, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "mp", ipe % grid % NMP, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )


    _ncCheck( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )

    IF( ipe % parameters % write_apex_neutrals )THEN

      ! Neutrals
      _ncCheck( nf90_def_var( ncid, "helium", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , helium_varid ) )
      _ncCheck( nf90_put_att( ncid, helium_varid, "long_name", "Neutral Helium Density" ) )
      _ncCheck( nf90_put_att( ncid, helium_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "oxygen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , oxygen_varid ) )
      _ncCheck( nf90_put_att( ncid, oxygen_varid, "long_name", "Neutral Oxygen Density" ) )
      _ncCheck( nf90_put_att( ncid, oxygen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "molecular_oxygen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , molecular_oxygen_varid ) )
      _ncCheck( nf90_put_att( ncid, molecular_oxygen_varid, "long_name", "Neutral Molecular Oxygen Density" ) )
      _ncCheck( nf90_put_att( ncid, molecular_oxygen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "molecular_nitrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , molecular_nitrogen_varid ) )
      _ncCheck( nf90_put_att( ncid, molecular_nitrogen_varid, "long_name", "Neutral Molecular Nitrogen Density" ) )
      _ncCheck( nf90_put_att( ncid, molecular_nitrogen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "nitrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , nitrogen_varid ) )
      _ncCheck( nf90_put_att( ncid, nitrogen_varid, "long_name", "Neutral Nitrogen Density" ) )
      _ncCheck( nf90_put_att( ncid, nitrogen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "hydrogen", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , hydrogen_varid ) )
      _ncCheck( nf90_put_att( ncid, hydrogen_varid, "long_name", "Neutral Hydrogen Density" ) )
      _ncCheck( nf90_put_att( ncid, hydrogen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "temperature", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , temperature_varid ) )
      _ncCheck( nf90_put_att( ncid, temperature_varid, "long_name", "Thermosphere Temperature" ) )
      _ncCheck( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

      _ncCheck( nf90_def_var( ncid, "u", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , u_varid ) )
      _ncCheck( nf90_put_att( ncid, u_varid, "long_name", "Zonal Velocity" ) )
      _ncCheck( nf90_put_att( ncid, u_varid, "units", "m s^{-1}" ) )

      _ncCheck( nf90_def_var( ncid, "v", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , v_varid ) )
      _ncCheck( nf90_put_att( ncid, v_varid, "long_name", "Meridional Velocity" ) )
      _ncCheck( nf90_put_att( ncid, v_varid, "units", "m s^{-1}" ) )

      _ncCheck( nf90_def_var( ncid, "w", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , w_varid ) )
      _ncCheck( nf90_put_att( ncid, w_varid, "long_name", "Radial Velocity" ) )
      _ncCheck( nf90_put_att( ncid, w_varid, "units", "m s^{-1}" ) )

    ENDIF

    IF( ipe % parameters % write_apex_eldyn )THEN

      _ncCheck( nf90_def_var( ncid, "phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , phi_varid ) )
      _ncCheck( nf90_put_att( ncid, phi_varid, "long_name", "Electric Potential" ) )
      _ncCheck( nf90_put_att( ncid, phi_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "mhd_phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) ,mhd_phi_varid ) )
      _ncCheck( nf90_put_att( ncid, mhd_phi_varid, "long_name", "Electric Potential - MHD Component" ) )
      _ncCheck( nf90_put_att( ncid, mhd_phi_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "v_exb1", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) ,ve1_varid ) )
      _ncCheck( nf90_put_att( ncid, ve1_varid, "long_name", "ExB Drift Velocity (Magnetic Colatitude component)" ) )
      _ncCheck( nf90_put_att( ncid, ve1_varid, "units", "m/s" ) )

      _ncCheck( nf90_def_var( ncid, "v_exb2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) ,ve2_varid ) )
      _ncCheck( nf90_put_att( ncid, ve2_varid, "long_name", "ExB Drift Velocity (Magnetic Longitude component)" ) )
      _ncCheck( nf90_put_att( ncid, ve2_varid, "units", "m/s" ) )

      _ncCheck( nf90_def_var( ncid, "hc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hc_varid ) )
      _ncCheck( nf90_put_att( ncid, hc_varid, "long_name", "Hall Conductivity" ) )
      _ncCheck( nf90_put_att( ncid, hc_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "pc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , pc_varid ) )
      _ncCheck( nf90_put_att( ncid, pc_varid, "long_name", "Pedersen Conductivity" ) )
      _ncCheck( nf90_put_att( ncid, pc_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "bc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , bc_varid ) )
      _ncCheck( nf90_put_att( ncid, bc_varid, "long_name", "Magnetic Field Aligned Conductivity" ) )
      _ncCheck( nf90_put_att( ncid, bc_varid, "units", "[Unknown]" ) )

    ENDIF

    ! Plasma
    _ncCheck( nf90_def_var( ncid, "O+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , op_varid ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "long_name", "Atomic oxygen ion number density (ground state)" ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "H+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , hp_varid ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "long_name", "Hydrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "He+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , hep_varid ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "long_name", "Helium ion number density" ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "N+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , np_varid ) )
    _ncCheck( nf90_put_att( ncid, np_varid, "long_name", "Nitrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, np_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "NO+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , nop_varid ) )
    _ncCheck( nf90_put_att( ncid, nop_varid, "long_name", "Nitrosonium ion number density" ) )
    _ncCheck( nf90_put_att( ncid, nop_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "O2+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , o2p_varid ) )
    _ncCheck( nf90_put_att( ncid, o2p_varid, "long_name", "Molecular Oxygen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, o2p_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "N2+", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , n2p_varid ) )
    _ncCheck( nf90_put_att( ncid, n2p_varid, "long_name", "Molecular Nitrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, n2p_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "O+(2D)", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , op2d_varid ) )
    _ncCheck( nf90_put_att( ncid, op2d_varid, "long_name", "Atomic oxygen ion number density (first excited state)" ) )
    _ncCheck( nf90_put_att( ncid, op2d_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "O+(2P)", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , op2p_varid ) )
    _ncCheck( nf90_put_att( ncid, op2p_varid, "long_name", "Atomic oxygen ion number density (second excited state)" ) )
    _ncCheck( nf90_put_att( ncid, op2p_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "ion_temp", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_temp_varid ) )
    _ncCheck( nf90_put_att( ncid, ion_temp_varid, "long_name", "Ion temperature" ) )
    _ncCheck( nf90_put_att( ncid, ion_temp_varid, "units", "K" ) )

    _ncCheck( nf90_def_var( ncid, "e_temp", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) ,electron_temp_varid ) )
    _ncCheck( nf90_put_att( ncid, electron_temp_varid, "long_name", "Electron temperature" ) )
    _ncCheck( nf90_put_att( ncid, electron_temp_varid, "units", "K" ) )

    _ncCheck( nf90_def_var( ncid, "ion_vel_op", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_vel_op_varid ) )
    _ncCheck( nf90_put_att( ncid, ion_vel_op_varid, "long_name", "O+ Velocity (B Parallel)" ) )
    _ncCheck( nf90_put_att( ncid, ion_vel_op_varid, "units", "?/s" ) )

    _ncCheck( nf90_def_var( ncid, "ion_vel_hp", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_vel_hp_varid ) )
    _ncCheck( nf90_put_att( ncid, ion_vel_hp_varid, "long_name", "H+ Velocity (B Parallel)" ) )
    _ncCheck( nf90_put_att( ncid, ion_vel_hp_varid, "units", "?/s" ) )

    _ncCheck( nf90_def_var( ncid, "ion_vel_hep", NF90_PREC, (/ z_dimid, x_dimid, y_dimid, time_dimid /) , ion_vel_hep_varid ) )
    _ncCheck( nf90_put_att( ncid, ion_vel_hep_varid, "long_name", "He+ Velocity (B Parallel)" ) )
    _ncCheck( nf90_put_att( ncid, ion_vel_hep_varid, "units", "?/s" ) )


    _ncCheck( nf90_enddef(ncid) )

    _ncCheck( nf90_put_var( ncid, time_varid, time ) )

    IF( ipe % parameters % write_apex_neutrals )THEN

      _ncCheck( nf90_put_var( ncid, helium_varid, ipe % neutrals % helium, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, oxygen_varid, ipe % neutrals % oxygen, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, molecular_oxygen_varid, ipe % neutrals % molecular_oxygen, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, molecular_nitrogen_varid, ipe % neutrals % molecular_nitrogen, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, nitrogen_varid, ipe % neutrals % nitrogen, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, hydrogen_varid, ipe % neutrals % hydrogen, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, temperature_varid, ipe % neutrals % temperature, recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, u_varid, ipe % neutrals % velocity_geographic(1,:,:,:), recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, v_varid, ipe % neutrals % velocity_geographic(2,:,:,:), recStart, recCount ) )
      _ncCheck( nf90_put_var( ncid, w_varid, ipe % neutrals % velocity_geographic(3,:,:,:), recStart, recCount ) )

    ENDIF

!      IF( ipe % parameters % write_apex_eldyn )THEN
!        _ncCheck( nf90_put_var( ncid, phi_varid, ipe % eldyn % electric_potential ) )
!        _ncCheck( nf90_put_var( ncid, mhd_phi_varid, ipe % eldyn % mhd_electric_potential ) )
!        _ncCheck( nf90_put_var( ncid, ve1_varid, ipe % eldyn % v_ExB_apex(1,:,:) ) )
!        _ncCheck( nf90_put_var( ncid, ve2_varid, ipe % eldyn % v_ExB_apex(2,:,:) ) )
!        _ncCheck( nf90_put_var( ncid, hc_varid, ipe % eldyn % hall_conductivity ) )
!        _ncCheck( nf90_put_var( ncid, pc_varid, ipe % eldyn % pedersen_conductivity ) )
!        _ncCheck( nf90_put_var( ncid, bc_varid, ipe % eldyn % b_parallel_conductivity ) )
!      ENDIF
!STOP

    _ncCheck( nf90_put_var( ncid, op_varid, ipe % plasma % ion_densities(1,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, hp_varid, ipe % plasma % ion_densities(2,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, hep_varid, ipe % plasma % ion_densities(3,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, np_varid, ipe % plasma % ion_densities(4,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, nop_varid, ipe % plasma % ion_densities(5,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, o2p_varid, ipe % plasma % ion_densities(6,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, n2p_varid, ipe % plasma % ion_densities(7,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, op2d_varid, ipe % plasma % ion_densities(8,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, op2p_varid, ipe % plasma % ion_densities(9,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, ion_temp_varid, ipe % plasma % ion_temperature, recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, electron_temp_varid, ipe % plasma % electron_temperature, recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, ion_vel_op_varid, ipe % plasma % ion_velocities(1,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, ion_vel_hp_varid, ipe % plasma % ion_velocities(2,:,:,:), recStart, recCount ) )
    _ncCheck( nf90_put_var( ncid, ion_vel_hep_varid, ipe % plasma % ion_velocities(3,:,:,:), recStart, recCount ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_NetCDF_IPE


  SUBROUTINE Write_Geographic_NetCDF_IPE( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid
    INTEGER :: op_varid, hp_varid, hep_varid, np_varid, n2p_varid, o2p_varid, nop_varid
    INTEGER :: op2d_varid, op2p_varid, phi_varid, mhd_phi_varid, exbu_varid, exbv_varid
    INTEGER :: ion_temp_varid, e_temp_varid, e_varid, hc_varid, pc_varid, bc_varid
    INTEGER :: ion_rate_varid, O_rate_varid, O2_rate_varid, N2_rate_varid
    INTEGER :: tec_varid, nmf2_varid, hmf2_varid


    error = 0

#ifdef HAVE_NETCDF
    CALL ipe % Geographic_Interpolation( )

    ! Time from reference time is calculated here
    time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    IF( prec == sp )THEN
      NF90_PREC = NF90_FLOAT
    ELSE
      NF90_PREC = NF90_DOUBLE
    ENDIF

    _ncCheck( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
    _ncCheck( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )

    _ncCheck( nf90_def_dim( ncid, "Z", nheights_geo, z_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

    _ncCheck( nf90_def_var( ncid, "Z", NF90_PREC, z_dimid, z_varid ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "long_name", "Altitude" ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "units", "km" ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, z_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "latitude", NF90_PREC, y_dimid, y_varid ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, y_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "longitude", NF90_PREC, x_dimid, x_varid ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, x_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )
   IF( ipe % parameters % write_geographic_neutrals )THEN

      _ncCheck( nf90_def_var( ncid, "helium", NF90_PREC, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , helium_varid ) )
      _ncCheck( nf90_put_att( ncid, helium_varid, "long_name", "NeutralHelium Density" ) )
      _ncCheck( nf90_put_att( ncid, helium_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "oxygen", NF90_PREC, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , oxygen_varid ) )
      _ncCheck( nf90_put_att( ncid, oxygen_varid, "long_name", "NeutralOxygen Density" ) )
      _ncCheck( nf90_put_att( ncid, oxygen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "molecular_oxygen", NF90_PREC, (/x_dimid, y_dimid, z_dimid, time_dimid /) , molecular_oxygen_varid ) )
      _ncCheck( nf90_put_att( ncid, molecular_oxygen_varid, "long_name","Neutral Molecular Oxygen Density" ) )
      _ncCheck( nf90_put_att( ncid, molecular_oxygen_varid, "units", "kgm^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "molecular_nitrogen", NF90_PREC, (/x_dimid, y_dimid, z_dimid, time_dimid /) , molecular_nitrogen_varid ) )
      _ncCheck( nf90_put_att( ncid, molecular_nitrogen_varid, "long_name","Neutral Molecular Nitrogen Density" ) )
      _ncCheck( nf90_put_att( ncid, molecular_nitrogen_varid, "units", "kgm^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "nitrogen", NF90_PREC, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , nitrogen_varid ) )
      _ncCheck( nf90_put_att( ncid, nitrogen_varid, "long_name", "NeutralNitrogen Density" ) )
      _ncCheck( nf90_put_att( ncid, nitrogen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "hydrogen", NF90_PREC, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , hydrogen_varid ) )
      _ncCheck( nf90_put_att( ncid, hydrogen_varid, "long_name", "NeutralHydrogen Density" ) )
      _ncCheck( nf90_put_att( ncid, hydrogen_varid, "units", "kg m^{-3}" ) )

      _ncCheck( nf90_def_var( ncid, "temperature", NF90_PREC, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , temperature_varid ) )
      _ncCheck( nf90_put_att( ncid, temperature_varid, "long_name","Thermosphere Temperature" ) )
      _ncCheck( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

      _ncCheck( nf90_def_var( ncid, "u", NF90_PREC, (/ x_dimid, y_dimid,z_dimid, time_dimid /) , u_varid ) )
      _ncCheck( nf90_put_att( ncid, u_varid, "long_name", "Apex1 Velocity" ))
      _ncCheck( nf90_put_att( ncid, u_varid, "units", "m s^{-1}" ) )

      _ncCheck( nf90_def_var( ncid, "v", NF90_PREC, (/ x_dimid, y_dimid,z_dimid, time_dimid /) , v_varid ) )
      _ncCheck( nf90_put_att( ncid, v_varid, "long_name", "Apex2 Velocity" ) )
      _ncCheck( nf90_put_att( ncid, v_varid, "units", "m s^{-1}" ) )

      _ncCheck( nf90_def_var( ncid, "w", NF90_PREC, (/ x_dimid, y_dimid,z_dimid, time_dimid /) , w_varid ) )
      _ncCheck( nf90_put_att( ncid, w_varid, "long_name", "Apex3 Velocity") )
      _ncCheck( nf90_put_att( ncid, w_varid, "units", "m s^{-1}" ) )

    ENDIF

    IF( ipe % parameters % write_geographic_eldyn )THEN

      _ncCheck( nf90_def_var( ncid, "phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , phi_varid ) )
      _ncCheck( nf90_put_att( ncid, phi_varid, "long_name", "Electric Potential" ) )
      _ncCheck( nf90_put_att( ncid, phi_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "mhd_phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , mhd_phi_varid ) )
      _ncCheck( nf90_put_att( ncid, mhd_phi_varid, "long_name", "Electric Potential - MHD Component" ) )
      _ncCheck( nf90_put_att( ncid, mhd_phi_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "exb_u", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , exbu_varid ) )
      _ncCheck( nf90_put_att( ncid, exbu_varid, "long_name", "Zonal component of ExB drift velocity" ) )
      _ncCheck( nf90_put_att( ncid, exbu_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "exb_v", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , exbv_varid ) )
      _ncCheck( nf90_put_att( ncid, exbv_varid, "long_name", "Meridional component of ExB drift velocity" ) )
      _ncCheck( nf90_put_att( ncid, exbv_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "hc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hc_varid ) )
      _ncCheck( nf90_put_att( ncid, hc_varid, "long_name", "Hall Conductivity" ) )
      _ncCheck( nf90_put_att( ncid, hc_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "pc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , pc_varid ) )
      _ncCheck( nf90_put_att( ncid, pc_varid, "long_name", "Pedersen Conductivity" ) )
      _ncCheck( nf90_put_att( ncid, pc_varid, "units", "[Unknown]" ) )

      _ncCheck( nf90_def_var( ncid, "bc", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , bc_varid ) )
      _ncCheck( nf90_put_att( ncid, bc_varid, "long_name", "Magnetic Field Aligned Conductivity" ) )
      _ncCheck( nf90_put_att( ncid, bc_varid, "units", "[Unknown]" ) )

    ENDIF

    ! Plasma
    _ncCheck( nf90_def_var( ncid, "O+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op_varid ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "long_name", "Atomic oxygen ion number density (ground state)" ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "H+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hp_varid ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "long_name", "Hydrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "He+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hep_varid ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "long_name", "Helium ion number density" ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "N+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , np_varid ) )
    _ncCheck( nf90_put_att( ncid, np_varid, "long_name", "Nitrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, np_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "NO+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , nop_varid ) )
    _ncCheck( nf90_put_att( ncid, nop_varid, "long_name", "Nitrosonium ion number density" ) )
    _ncCheck( nf90_put_att( ncid, nop_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "O2+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , o2p_varid ) )
    _ncCheck( nf90_put_att( ncid, o2p_varid, "long_name", "Molecular Oxygen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, o2p_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "N2+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , n2p_varid ) )
    _ncCheck( nf90_put_att( ncid, n2p_varid, "long_name", "Molecular Nitrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, n2p_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "O+(2D)", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op2d_varid ) )
    _ncCheck( nf90_put_att( ncid, op2d_varid, "long_name", "Atomic oxygen ion number density (first excited state)" ) )
    _ncCheck( nf90_put_att( ncid, op2d_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "O+(2P)", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op2p_varid ) )
    _ncCheck( nf90_put_att( ncid, op2p_varid, "long_name", "Atomic oxygen ion number density  (second excited state)" ) )
    _ncCheck( nf90_put_att( ncid, op2p_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "ion_temp", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , ion_temp_varid ) )
    _ncCheck( nf90_put_att( ncid, ion_temp_varid, "long_name", "Ion temperature" ) )
    _ncCheck( nf90_put_att( ncid, ion_temp_varid, "units", "K" ) )

    _ncCheck( nf90_def_var( ncid, "e", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , e_varid ) )
    _ncCheck( nf90_put_att( ncid, e_varid, "long_name", "Electron number density" ) )
    _ncCheck( nf90_put_att( ncid, e_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "e_temp", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , e_temp_varid ) )
    _ncCheck( nf90_put_att( ncid, e_temp_varid, "long_name", "Electron temperature" ) )
    _ncCheck( nf90_put_att( ncid, e_temp_varid, "units", "K" ) )

    _ncCheck( nf90_def_var( ncid, "tec", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , tec_varid ) )
    _ncCheck( nf90_put_att( ncid, tec_varid, "long_name", "Total Electron Content" ) )
    _ncCheck( nf90_put_att( ncid, tec_varid, "units", " TECU" ) )

    _ncCheck( nf90_def_var( ncid, "nmf2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , nmf2_varid ) )
    _ncCheck( nf90_put_att( ncid, nmf2_varid, "long_name", "Peak Electron Number Density" ) )
    _ncCheck( nf90_put_att( ncid, nmf2_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "hmf2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hmf2_varid ) )
    _ncCheck( nf90_put_att( ncid, hmf2_varid, "long_name", "Height of Peak Electron Number Density" ) )
    _ncCheck( nf90_put_att( ncid, hmf2_varid, "units", " km" ) )

    _ncCheck( nf90_def_var( ncid, "aur_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , ion_rate_varid ) )
    _ncCheck( nf90_put_att( ncid, ion_rate_varid, "long_name", "Total Ionization Rate from Auroral Precipitation" ) )
    _ncCheck( nf90_put_att( ncid, ion_rate_varid, "units", "[Unknown]" ) )

    _ncCheck( nf90_def_var( ncid, "O+_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , O_rate_varid ) )
    _ncCheck( nf90_put_att( ncid, O_rate_varid, "long_name", "Oxygen Ionization Rate from Auroral Precipitation" ) )
    _ncCheck( nf90_put_att( ncid, O_rate_varid, "units", "[Unknown]" ) )

    _ncCheck( nf90_def_var( ncid, "O2+_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , O2_rate_varid ) )
    _ncCheck( nf90_put_att( ncid, O2_rate_varid, "long_name", "Molecular Oxygen Ionization Rate from Auroral Precipitation" ) )
    _ncCheck( nf90_put_att( ncid, O2_rate_varid, "units", "[Unknown]" ) )

    _ncCheck( nf90_def_var( ncid, "N2+_precip", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , N2_rate_varid ) )
    _ncCheck( nf90_put_att( ncid, N2_rate_varid, "long_name", "Molecular Nitrogen Ionization Rate from Auroral Precipitation" ) )
    _ncCheck( nf90_put_att( ncid, N2_rate_varid, "units", "[Unknown]" ) )

    _ncCheck( nf90_enddef(ncid) )

    _ncCheck( nf90_put_var( ncid, time_varid, time ) )

    IF( ipe % parameters % write_geographic_neutrals )THEN

      _ncCheck( nf90_put_var( ncid, x_varid, ipe % grid % longitude_geo ) )
      _ncCheck( nf90_put_var( ncid, y_varid, ipe % grid % latitude_geo ) )
      _ncCheck( nf90_put_var( ncid, z_varid, ipe % grid % altitude_geo ) )
      _ncCheck( nf90_put_var( ncid, helium_varid, ipe % neutrals %geo_helium ) )
      _ncCheck( nf90_put_var( ncid, oxygen_varid, ipe % neutrals %geo_oxygen ) )
      _ncCheck( nf90_put_var( ncid, molecular_oxygen_varid, ipe % neutrals %geo_molecular_oxygen ) )
      _ncCheck( nf90_put_var( ncid, molecular_nitrogen_varid, ipe % neutrals% geo_molecular_nitrogen ) )
      _ncCheck( nf90_put_var( ncid, nitrogen_varid, ipe % neutrals %geo_nitrogen ) )
      _ncCheck( nf90_put_var( ncid, hydrogen_varid, ipe % neutrals %geo_hydrogen ) )
      _ncCheck( nf90_put_var( ncid, temperature_varid, ipe % neutrals %geo_temperature ) )
      _ncCheck( nf90_put_var( ncid, u_varid, ipe % neutrals %geo_velocity(1,:,:,:) ) )
      _ncCheck( nf90_put_var( ncid, v_varid, ipe % neutrals %geo_velocity(2,:,:,:) ) )
      _ncCheck( nf90_put_var( ncid, w_varid, ipe % neutrals %geo_velocity(3,:,:,:) ) )

    ENDIF

    IF( ipe % parameters % write_geographic_eldyn )THEN

      _ncCheck( nf90_put_var( ncid, phi_varid, ipe % eldyn % geo_electric_potential ) )
      _ncCheck( nf90_put_var( ncid, mhd_phi_varid, ipe % eldyn % geo_mhd_electric_potential ) )
      _ncCheck( nf90_put_var( ncid, exbu_varid, ipe % eldyn % geo_v_ExB_geographic(1,:,:) ) )
      _ncCheck( nf90_put_var( ncid, exbv_varid, ipe % eldyn % geo_v_ExB_geographic(2,:,:) ) )
      _ncCheck( nf90_put_var( ncid, hc_varid, ipe % eldyn % geo_hall_conductivity ) )
      _ncCheck( nf90_put_var( ncid, pc_varid, ipe % eldyn % geo_pedersen_conductivity ) )
      _ncCheck( nf90_put_var( ncid, bc_varid, ipe % eldyn % geo_b_parallel_conductivity ) )

    ENDIF

    _ncCheck( nf90_put_var( ncid, op_varid,  ipe % plasma % geo_ion_densities(1,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, hp_varid,  ipe % plasma % geo_ion_densities(2,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, hep_varid, ipe % plasma % geo_ion_densities(3,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, np_varid,  ipe % plasma % geo_ion_densities(4,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, nop_varid, ipe % plasma % geo_ion_densities(5,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, o2p_varid, ipe % plasma % geo_ion_densities(6,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, n2p_varid, ipe % plasma % geo_ion_densities(7,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, op2d_varid,  ipe % plasma % geo_ion_densities(8,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, op2p_varid,  ipe % plasma % geo_ion_densities(9,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, ion_temp_varid, ipe % plasma % geo_ion_temperature ) )
    _ncCheck( nf90_put_var( ncid, e_temp_varid, ipe % plasma % geo_electron_temperature ) )
    _ncCheck( nf90_put_var( ncid, e_varid, ipe % plasma % geo_electron_density ) )
    _ncCheck( nf90_put_var( ncid, tec_varid, ipe % plasma % geo_tec ) )
    _ncCheck( nf90_put_var( ncid, nmf2_varid, ipe % plasma % geo_nmf2 ) )
    _ncCheck( nf90_put_var( ncid, hmf2_varid, ipe % plasma % geo_hmf2 ) )
    _ncCheck( nf90_put_var( ncid, ion_rate_varid, ipe % plasma % geo_ionization_rates(1,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, O_rate_varid, ipe % plasma % geo_ionization_rates(2,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, O2_rate_varid, ipe % plasma % geo_ionization_rates(3,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, N2_rate_varid, ipe % plasma % geo_ionization_rates(4,:,:,:) ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_Geographic_NetCDF_IPE


  SUBROUTINE Write_NMF2_TEC_NetCDF_IPE( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: tec_varid, nmf2_varid, hmf2_varid


    error = 0

#ifdef HAVE_NETCDF
    CALL ipe % Geographic_Interpolation( )

    ! Time from reference time is calculated here
    time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    IF( prec == sp )THEN
      NF90_PREC = NF90_FLOAT
    ELSE
      NF90_PREC = NF90_DOUBLE
    ENDIF

    _ncCheck( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
    _ncCheck( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )

    _ncCheck( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )


    _ncCheck( nf90_def_var( ncid, "latitude", NF90_PREC, y_dimid, y_varid ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, y_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "longitude", NF90_PREC, x_dimid, x_varid ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, x_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )

    ! Plasma

    _ncCheck( nf90_def_var( ncid, "tec", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , tec_varid ) )
    _ncCheck( nf90_put_att( ncid, tec_varid, "long_name", "Total Electron Content" ) )
    _ncCheck( nf90_put_att( ncid, tec_varid, "units", " TECU" ) )

    _ncCheck( nf90_def_var( ncid, "nmf2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , nmf2_varid ) )
    _ncCheck( nf90_put_att( ncid, nmf2_varid, "long_name", "Peak Electron Number Density" ) )
    _ncCheck( nf90_put_att( ncid, nmf2_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "hmf2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hmf2_varid ) )
    _ncCheck( nf90_put_att( ncid, hmf2_varid, "long_name", "Height of Peak Electron Number Density" ) )
    _ncCheck( nf90_put_att( ncid, hmf2_varid, "units", " km" ) )

    _ncCheck( nf90_enddef(ncid) )

    _ncCheck( nf90_put_var( ncid, time_varid, time ) )

    _ncCheck( nf90_put_var( ncid, tec_varid, ipe % plasma % geo_tec ) )
    _ncCheck( nf90_put_var( ncid, nmf2_varid, ipe % plasma % geo_nmf2 ) )
    _ncCheck( nf90_put_var( ncid, hmf2_varid, ipe % plasma % geo_hmf2 ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_NMF2_TEC_NetCDF_IPE

  SUBROUTINE Write_Electron_density_NetCDF_IPE( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    REAL(prec) :: time
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: e_varid
    INTEGER :: op_varid, hp_varid, hep_varid, np_varid, n2p_varid, o2p_varid, nop_varid
    INTEGER :: op2d_varid, op2p_varid
    INTEGER :: tec_varid, nmf2_varid, hmf2_varid
    INTEGER :: helium_varid, oxygen_varid, molecular_oxygen_varid
    INTEGER :: molecular_nitrogen_varid, nitrogen_varid, hydrogen_varid
    INTEGER :: temperature_varid, u_varid, v_varid, w_varid

    REAL(4) :: geo_ion_densities_real4(9,nlon_geo,nlat_geo,nheights_geo)

    REAL(4) :: geo_helium_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_oxygen_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_molecular_oxygen_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_molecular_nitrogen_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_nitrogen_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_hydrogen_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_temperature_real4(nlon_geo,nlat_geo,nheights_geo)
    REAL(4) :: geo_velocity_real4(3,nlon_geo,nlat_geo,nheights_geo)

    error = 0

#ifdef HAVE_NETCDF
    CALL ipe % Geographic_Interpolation( )

    geo_ion_densities_real4(:,:,:,:)    = real(ipe % plasma % geo_ion_densities(:,:,:,:))

    geo_helium_real4(:,:,:)             = real(ipe % neutrals % geo_helium(:,:,:))
    geo_oxygen_real4(:,:,:)             = real(ipe % neutrals % geo_oxygen(:,:,:))
    geo_molecular_oxygen_real4(:,:,:)   = real(ipe % neutrals % geo_molecular_oxygen(:,:,:))
    geo_molecular_nitrogen_real4(:,:,:) = real(ipe % neutrals % geo_molecular_nitrogen(:,:,:))
    geo_nitrogen_real4(:,:,:)           = real(ipe % neutrals % geo_nitrogen(:,:,:))
    geo_hydrogen_real4(:,:,:)           = real(ipe % neutrals % geo_hydrogen(:,:,:))
    geo_temperature_real4(:,:,:)        = real(ipe % neutrals % geo_temperature(:,:,:))
    geo_velocity_real4(:,:,:,:)         = real(ipe % neutrals % geo_velocity(:,:,:,:))

    time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    print *, 'open',TRIM(filename)
    _ncCheck( nf90_create('./test', NF90_NETCDF4, ncid))
    print *, 'open',TRIM(filename)
    _ncCheck( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )


    print *, 'def_dim'
    _ncCheck( nf90_def_dim( ncid, "altitude", nheights_geo, z_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )


    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "altitude", NF90_FLOAT, z_dimid, z_varid ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "long_name", "Altitude" ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "units", "km" ) )


    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "latitude", NF90_FLOAT, y_dimid, y_varid ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )


    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "longitude", NF90_FLOAT, x_dimid, x_varid ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "time", NF90_FLOAT, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "o_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op_varid ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "long_name", "O+ number density" ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "h_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hp_varid ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "long_name", "H+ number density" ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "he_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hep_varid ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "long_name", "He+ number density" ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "n_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , np_varid ) )
    _ncCheck( nf90_put_att( ncid, np_varid, "long_name", "N+ number density" ) )
    _ncCheck( nf90_put_att( ncid, np_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "no_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , nop_varid ) )
    _ncCheck( nf90_put_att( ncid, nop_varid, "long_name", "NO+ number density" ) )
    _ncCheck( nf90_put_att( ncid, nop_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "o2_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , o2p_varid ) )
    _ncCheck( nf90_put_att( ncid, o2p_varid, "long_name", "O2+ number density" ) )
    _ncCheck( nf90_put_att( ncid, o2p_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "n2_plus", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , n2p_varid ) )
    _ncCheck( nf90_put_att( ncid, n2p_varid, "long_name", "N2+ number density" ) )
    _ncCheck( nf90_put_att( ncid, n2p_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "o_plus_2d", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op2d_varid ) )
    _ncCheck( nf90_put_att( ncid, op2d_varid, "long_name", "O+(2D) number density" ) )
    _ncCheck( nf90_put_att( ncid, op2d_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "o_plus_2p", NF90_FLOAT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op2p_varid ) )
    _ncCheck( nf90_put_att( ncid, op2p_varid, "long_name", "O+(2P) number density" ) )
    _ncCheck( nf90_put_att( ncid, op2p_varid, "units", " m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "helium", NF90_FLOAT, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , helium_varid ) )
    _ncCheck( nf90_put_att( ncid, helium_varid, "long_name", "Neutral He Density" ) )
    _ncCheck( nf90_put_att( ncid, helium_varid, "units", "kg m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "oxygen", NF90_FLOAT, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , oxygen_varid ) )
    _ncCheck( nf90_put_att( ncid, oxygen_varid, "long_name", "Neutral O Density" ) )
    _ncCheck( nf90_put_att( ncid, oxygen_varid, "units", "kg m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "molecular_oxygen", NF90_FLOAT, (/x_dimid, y_dimid, z_dimid, time_dimid /) , molecular_oxygen_varid ) )
    _ncCheck( nf90_put_att( ncid, molecular_oxygen_varid, "long_name","Neutral O2 Density" ) )
    _ncCheck( nf90_put_att( ncid, molecular_oxygen_varid, "units", "kgm^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "molecular_nitrogen", NF90_FLOAT, (/x_dimid, y_dimid, z_dimid, time_dimid /) , molecular_nitrogen_varid ) )
    _ncCheck( nf90_put_att( ncid, molecular_nitrogen_varid, "long_name","Neutral N2 Density" ) )
    _ncCheck( nf90_put_att( ncid, molecular_nitrogen_varid, "units", "kgm^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "nitrogen", NF90_FLOAT, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , nitrogen_varid ) )
    _ncCheck( nf90_put_att( ncid, nitrogen_varid, "long_name", "Neutral N Density" ) )
    _ncCheck( nf90_put_att( ncid, nitrogen_varid, "units", "kg m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "hydrogen", NF90_FLOAT, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , hydrogen_varid ) )
    _ncCheck( nf90_put_att( ncid, hydrogen_varid, "long_name", "Neutral H Density" ) )
    _ncCheck( nf90_put_att( ncid, hydrogen_varid, "units", "kg m^{-3}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "temperature", NF90_FLOAT, (/ x_dimid,y_dimid, z_dimid, time_dimid /) , temperature_varid ) )
    _ncCheck( nf90_put_att( ncid, temperature_varid, "long_name","Thermosphere Temperature" ) )
    _ncCheck( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "u", NF90_FLOAT, (/ x_dimid, y_dimid,z_dimid, time_dimid /) , u_varid ) )
    _ncCheck( nf90_put_att( ncid, u_varid, "long_name", "Apex1 Velocity" ))
    _ncCheck( nf90_put_att( ncid, u_varid, "units", "m s^{-1}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "v", NF90_FLOAT, (/ x_dimid, y_dimid,z_dimid, time_dimid /) , v_varid ) )
    _ncCheck( nf90_put_att( ncid, v_varid, "long_name", "Apex2 Velocity" ) )
    _ncCheck( nf90_put_att( ncid, v_varid, "units", "m s^{-1}" ) )

    print *, 'def_var'
    _ncCheck( nf90_def_var( ncid, "w", NF90_FLOAT, (/ x_dimid, y_dimid,z_dimid, time_dimid /) , w_varid ) )
    _ncCheck( nf90_put_att( ncid, w_varid, "long_name", "Apex3 Velocity") )
    _ncCheck( nf90_put_att( ncid, w_varid, "units", "m s^{-1}" ) )

    print *, 'enddef'
    _ncCheck( nf90_enddef(ncid) )

    print *, 'put_dim'
    _ncCheck( nf90_put_var( ncid, x_varid, real(ipe % grid % longitude_geo )) )
    _ncCheck( nf90_put_var( ncid, y_varid, real(ipe % grid % latitude_geo )) )
    _ncCheck( nf90_put_var( ncid, z_varid, real(ipe % grid % altitude_geo )) )
    _ncCheck( nf90_put_var( ncid, time_dimid, time ) )

    print *, 'put_var'
    _ncCheck( nf90_put_var( ncid, op_varid,   geo_ion_densities_real4(1,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, hp_varid,   geo_ion_densities_real4(2,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, hep_varid,  geo_ion_densities_real4(3,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, np_varid,   geo_ion_densities_real4(4,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, nop_varid,  geo_ion_densities_real4(5,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, o2p_varid,  geo_ion_densities_real4(6,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, n2p_varid,  geo_ion_densities_real4(7,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, op2d_varid, geo_ion_densities_real4(8,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, op2p_varid, geo_ion_densities_real4(9,:,:,:) ) )

    print *, 'put_var'
    _ncCheck( nf90_put_var( ncid, helium_varid, geo_helium_real4))
    _ncCheck( nf90_put_var( ncid, oxygen_varid, geo_oxygen_real4 ))
    _ncCheck( nf90_put_var( ncid, molecular_oxygen_varid, geo_molecular_oxygen_real4))
    _ncCheck( nf90_put_var( ncid, molecular_nitrogen_varid, geo_molecular_nitrogen_real4))
    _ncCheck( nf90_put_var( ncid, nitrogen_varid, geo_nitrogen_real4 ))
    _ncCheck( nf90_put_var( ncid, hydrogen_varid, geo_hydrogen_real4 ))
    _ncCheck( nf90_put_var( ncid, temperature_varid, geo_temperature_real4))
    _ncCheck( nf90_put_var( ncid, u_varid, geo_velocity_real4(1,:,:,:) ))
    _ncCheck( nf90_put_var( ncid, v_varid, geo_velocity_real4(2,:,:,:) ))
    _ncCheck( nf90_put_var( ncid, w_varid, geo_velocity_real4(3,:,:,:) ))

    print *, 'close'
    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_Electron_density_NetCDF_IPE

  SUBROUTINE Write_Electron_density_NetCDF_IPE_integer( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    REAL(prec) :: time
    REAL(prec) :: e_max,e_min
    REAL(prec) :: i_max,i_min
    REAL(prec) :: ti_max,ti_min
    REAL(prec) :: te_max,te_min
    REAL(prec) :: e_den(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    REAL(prec) :: i_den(1:9,1:nlon_geo,1:nlat_geo,1:nheights_geo)
    REAL(prec) :: i_den_single_ion(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    REAL(prec) :: ti(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    REAL(prec) :: te(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    REAL(prec), PARAMETER :: intmax = 32767.0_prec
    REAL(prec), PARAMETER :: intmin = -32768.0_prec
    REAL(prec) :: minmax_e(1:2)
    REAL(prec) :: minmax_i(1:9,1:2)
    REAL(prec) :: minmax_ti(1:2)
    REAL(prec) :: minmax_te(1:2)
    INTEGER*2 :: e_den_int(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    INTEGER*2 :: i_den_int(1:9,1:nlon_geo,1:nlat_geo,1:nheights_geo)
    INTEGER*2 :: ti_int(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    INTEGER*2 :: te_int(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    INTEGER*2 :: lon_int(1:nlon_geo)
    INTEGER*2 :: lat_int(1:nlat_geo)
    INTEGER*2 :: ht_int(1:nheights_geo)
    INTEGER :: NF90_PREC
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid, minmax_dimid, ions_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid, minmax_varid, minmax_i_varid
    INTEGER :: minmax_ti_varid, minmax_te_varid
    INTEGER :: e_varid, i_varid, ti_varid, te_varid
    INTEGER :: tec_varid, nmf2_varid, hmf2_varid
    INTEGER :: ilon, ilat, ih , ion


    error = 0

#ifdef HAVE_NETCDF
    CALL ipe % Geographic_Interpolation( )

    ! Time from reference time is calculated here
    time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    IF( prec == sp )THEN
      NF90_PREC = NF90_FLOAT
    ELSE
      NF90_PREC = NF90_DOUBLE
    ENDIF

    e_den(1:nlon_geo,1:nlat_geo,1:nheights_geo) = ipe % plasma % geo_electron_density(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    i_den(1:9,1:nlon_geo,1:nlat_geo,1:nheights_geo) = ipe % plasma % geo_ion_densities(1:9,1:nlon_geo,1:nlat_geo,1:nheights_geo)
    ti(1:nlon_geo,1:nlat_geo,1:nheights_geo) = ipe % plasma % geo_ion_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo)
    te(1:nlon_geo,1:nlat_geo,1:nheights_geo) = ipe % plasma % geo_electron_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo)

    lon_int(1:nlon_geo) = nint(ipe % grid % longitude_geo(1:nlon_geo))
    lat_int(1:nlat_geo) = nint(ipe % grid % latitude_geo(1:nlat_geo))
    ht_int(1:nheights_geo) = nint(ipe % grid % altitude_geo(1:nheights_geo))

    e_max = maxval(e_den)
    e_min = minval(e_den)
    minmax_e(1) = e_min
    minmax_e(2) = e_max
    do ilon = 1 , nlon_geo
    do ilat = 1 , nlat_geo
    do ih = 1 , nheights_geo
    e_den_int(ilon,ilat,ih) = nint((((e_den(ilon,ilat,ih) - e_min) / (e_max - e_min)) * (intmax - intmin)) + intmin)
    enddo
    enddo
    enddo

    ti_max = maxval(ti)
    ti_min = minval(ti)
    minmax_ti(1) = ti_min
    minmax_ti(2) = ti_max
    do ilon = 1 , nlon_geo
    do ilat = 1 , nlat_geo
    do ih = 1 , nheights_geo
    ti_int(ilon,ilat,ih) = nint((((ti(ilon,ilat,ih) - ti_min) / (ti_max - ti_min)) * (intmax - intmin)) + intmin)
    enddo
    enddo
    enddo

    te_max = maxval(te)
    te_min = minval(te)
    minmax_te(1) = te_min
    minmax_te(2) = te_max
    do ilon = 1 , nlon_geo
    do ilat = 1 , nlat_geo
    do ih = 1 , nheights_geo
    te_int(ilon,ilat,ih) = nint((((te(ilon,ilat,ih) - te_min) / (te_max - te_min)) * (intmax - intmin)) + intmin)
    enddo
    enddo
    enddo

    do ion = 1 , 9
    i_den_single_ion(1:nlon_geo,1:nlat_geo,1:nheights_geo) = ipe % plasma % geo_ion_densities(ion,1:nlon_geo,1:nlat_geo,1:nheights_geo)
    i_max = maxval(i_den_single_ion)
    i_min = minval(i_den_single_ion)
    minmax_i(ion,1) = i_min
    minmax_i(ion,2) = i_max
    do ilon = 1 , nlon_geo
    do ilat = 1 , nlat_geo
    do ih = 1 , nheights_geo
    i_den_int(ion,ilon,ilat,ih) = nint((((i_den(ion,ilon,ilat,ih) - i_min) / (i_max - i_min)) * (intmax - intmin)) + intmin)
    enddo
    enddo
    enddo
    enddo

    _ncCheck( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
    _ncCheck( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )

    _ncCheck( nf90_def_dim( ncid, "Z", nheights_geo, z_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "ions", 9, ions_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "minmax", 2, minmax_dimid ) )

    _ncCheck( nf90_def_var( ncid, "Z", NF90_SHORT, z_dimid, z_varid ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "long_name", "Altitude" ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "units", "km" ) )

    _ncCheck( nf90_def_var( ncid, "latitude", NF90_SHORT, y_dimid, y_varid ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )

    _ncCheck( nf90_def_var( ncid, "longitude", NF90_SHORT, x_dimid, x_varid ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )

    _ncCheck( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "minmax", NF90_PREC, minmax_dimid, minmax_varid ) )
    _ncCheck( nf90_def_var( ncid, "minmax_ti", NF90_PREC, minmax_dimid, minmax_ti_varid ) )
    _ncCheck( nf90_def_var( ncid, "minmax_te", NF90_PREC, minmax_dimid, minmax_te_varid ) )
    _ncCheck( nf90_def_var( ncid, "minmax_i", NF90_PREC, (/ ions_dimid, minmax_dimid /), minmax_i_varid ) )

    ! Plasma

    _ncCheck( nf90_def_var( ncid, "e", NF90_SHORT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , e_varid ) )
    _ncCheck( nf90_put_att( ncid, e_varid, "long_name", "Electron number density" ) )
    _ncCheck( nf90_put_att( ncid, e_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "ti", NF90_SHORT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , ti_varid ) )
    _ncCheck( nf90_put_att( ncid, ti_varid, "long_name", "Ion temperature" ) )
    _ncCheck( nf90_put_att( ncid, ti_varid, "units", " K" ) )

    _ncCheck( nf90_def_var( ncid, "te", NF90_SHORT, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , te_varid ) )
    _ncCheck( nf90_put_att( ncid, te_varid, "long_name", "Electron temperature" ) )
    _ncCheck( nf90_put_att( ncid, te_varid, "units", " K" ) )

    _ncCheck( nf90_def_var( ncid, "i", NF90_SHORT, (/ ions_dimid, x_dimid, y_dimid, z_dimid, time_dimid /) , i_varid ) )
    _ncCheck( nf90_put_att( ncid, i_varid, "long_name", "Ions number density" ) )
    _ncCheck( nf90_put_att( ncid, i_varid, "units", " m^{-3}" ) )

!     _ncCheck( nf90_def_var( ncid, "tec", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , tec_varid ) )
!     _ncCheck( nf90_put_att( ncid, tec_varid, "long_name", "Total Electron Content" ) )
!     _ncCheck( nf90_put_att( ncid, tec_varid, "units", " TECU" ) )

!     _ncCheck( nf90_def_var( ncid, "nmf2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , nmf2_varid ) )
!     _ncCheck( nf90_put_att( ncid, nmf2_varid, "long_name", "Peak Electron Number Density" ) )
!     _ncCheck( nf90_put_att( ncid, nmf2_varid, "units", " m^{-3}" ) )

!     _ncCheck( nf90_def_var( ncid, "hmf2", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) , hmf2_varid ) )
!     _ncCheck( nf90_put_att( ncid, hmf2_varid, "long_name", "Height of Peak Electron Number Density" ) )
!     _ncCheck( nf90_put_att( ncid, hmf2_varid, "units", " km" ) )

    _ncCheck( nf90_enddef(ncid) )

    _ncCheck( nf90_put_var( ncid, x_varid, lon_int ) )
    _ncCheck( nf90_put_var( ncid, y_varid, lat_int ) )
    _ncCheck( nf90_put_var( ncid, z_varid, ht_int ) )
    _ncCheck( nf90_put_var( ncid, time_varid, time ) )

!     _ncCheck( nf90_put_var( ncid, tec_varid, ipe % plasma % geo_tec ) )
!     _ncCheck( nf90_put_var( ncid, nmf2_varid, ipe % plasma % geo_nmf2 ) )
!     _ncCheck( nf90_put_var( ncid, hmf2_varid, ipe % plasma % geo_hmf2 ) )

    _ncCheck( nf90_put_var( ncid, minmax_varid, minmax_e ) )
    _ncCheck( nf90_put_var( ncid, minmax_i_varid, minmax_i ) )
    _ncCheck( nf90_put_var( ncid, minmax_ti_varid, minmax_ti ) )
    _ncCheck( nf90_put_var( ncid, minmax_te_varid, minmax_te ) )
    _ncCheck( nf90_put_var( ncid, e_varid, e_den_int ) )
    _ncCheck( nf90_put_var( ncid, i_varid, i_den_int ) )
    _ncCheck( nf90_put_var( ncid, ti_varid, ti_int ) )
    _ncCheck( nf90_put_var( ncid, te_varid, te_int ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_Electron_density_NetCDF_IPE_integer


  SUBROUTINE Write_Reduced_Geographic_NetCDF_IPE( ipe, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER,            INTENT(out)   :: error

    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: op_varid, hp_varid, hep_varid
    INTEGER :: e_varid


    error = 0

#ifdef HAVE_NETCDF
    CALL ipe % Geographic_Interpolation( )

    ! Time from reference time is calculated here
    time = ipe % time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    IF( prec == sp )THEN
      NF90_PREC = NF90_FLOAT
    ELSE
      NF90_PREC = NF90_DOUBLE
    ENDIF

    _ncCheck( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))
    _ncCheck( nf90_put_att( ncid, NF90_GLOBAL, "Version", 1.0) )

    _ncCheck( nf90_def_dim( ncid, "Z", nheights_geo, z_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "longitude", nlon_geo, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "latitude", nlat_geo, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

    _ncCheck( nf90_def_var( ncid, "Z", NF90_PREC, z_dimid, z_varid ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "long_name", "Altitude" ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "units", "km" ) )
    _ncCheck( nf90_put_att( ncid, z_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, z_varid, "missing_value", fillValue) )


    _ncCheck( nf90_def_var( ncid, "latitude", NF90_PREC, y_dimid, y_varid ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )
    _ncCheck( nf90_put_att( ncid, y_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, y_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "longitude", NF90_PREC, x_dimid, x_varid ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )
    _ncCheck( nf90_put_att( ncid, x_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, x_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
    _ncCheck( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )

    _ncCheck( nf90_def_var( ncid, "O+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , op_varid ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "long_name", "Atomic oxygen ion number density (ground state)" ) )
    _ncCheck( nf90_put_att( ncid, op_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "H+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hp_varid ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "long_name", "Hydrogen ion number density" ) )
    _ncCheck( nf90_put_att( ncid, hp_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "He+", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , hep_varid ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "long_name", "Helium ion number density" ) )
    _ncCheck( nf90_put_att( ncid, hep_varid, "units", " m^{-3}" ) )

    _ncCheck( nf90_def_var( ncid, "e", NF90_PREC, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , e_varid ) )
    _ncCheck( nf90_put_att( ncid, e_varid, "long_name", "Electron number density" ) )
    _ncCheck( nf90_put_att( ncid, e_varid, "units", " m^{-3}" ) )
    _ncCheck( nf90_enddef(ncid) )

    _ncCheck( nf90_put_var( ncid, x_varid, ipe % grid % longitude_geo ) )
    _ncCheck( nf90_put_var( ncid, y_varid, ipe % grid % latitude_geo ) )
    _ncCheck( nf90_put_var( ncid, z_varid, ipe % grid % altitude_geo ) )
    _ncCheck( nf90_put_var( ncid, time_varid, time ) )
    _ncCheck( nf90_put_var( ncid, op_varid,  ipe % plasma % geo_ion_densities(1,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, hp_varid,  ipe % plasma % geo_ion_densities(2,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, hep_varid, ipe % plasma % geo_ion_densities(3,:,:,:) ) )
    _ncCheck( nf90_put_var( ncid, e_varid, ipe % plasma % geo_electron_density ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_Reduced_Geographic_NetCDF_IPE

END MODULE IPE_Model_Class
