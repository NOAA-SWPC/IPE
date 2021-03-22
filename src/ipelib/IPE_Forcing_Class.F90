MODULE IPE_Forcing_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines
USE IPE_Model_Parameters_Class
USE IPE_MPI_Layer_Class
USE ipe_error_module
USE comio

IMPLICIT NONE

  TYPE IPE_Forcing

    LOGICAL                 :: coupled

    REAL(prec)              :: dt
    REAL(prec)              :: current_time
    REAL(prec)              :: start_time
    INTEGER                 :: current_index
    INTEGER                 :: max_read_index
    INTEGER                 :: ifp_interval
    INTEGER                 :: ifp_skip
    !
    REAL(prec), ALLOCATABLE :: f107(:)
    REAL(prec), ALLOCATABLE :: f107_81day_avg(:)
    REAL(prec), ALLOCATABLE :: kp(:)
    REAL(prec), ALLOCATABLE :: kp_1day_avg(:)
    REAL(prec), ALLOCATABLE :: ap(:)
    REAL(prec), ALLOCATABLE :: ap_1day_avg(:)
    REAL(prec), ALLOCATABLE :: nhemi_power(:)
    INTEGER, ALLOCATABLE    :: nhemi_power_index(:)
    REAL(prec), ALLOCATABLE :: shemi_power(:)
    INTEGER, ALLOCATABLE    :: shemi_power_index(:)
    ! Solar wind drivers
    REAL(prec), ALLOCATABLE :: solarwind_angle(:)
    REAL(prec), ALLOCATABLE :: solarwind_velocity(:)
    REAL(prec), ALLOCATABLE :: solarwind_Bz(:)
    REAL(prec), ALLOCATABLE :: solarwind_By(:)
    REAL(prec), ALLOCATABLE :: solarwind_Bt(:)
    REAL(prec), ALLOCATABLE :: solarwind_density(:)

    ! Tiros
    REAL(prec), ALLOCATABLE :: emaps(:,:,:)
    REAL(prec), ALLOCATABLE :: cmaps(:,:,:)
    REAL(prec), ALLOCATABLE :: djspectra(:,:)

    !
    REAL(prec) :: sun_longitude

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Forcing
      PROCEDURE :: Trash => Trash_IPE_Forcing

      PROCEDURE :: Update_Current_Index
      PROCEDURE :: GetAP
      PROCEDURE :: GetKP

      PROCEDURE, PRIVATE :: Allocate_IFP_IPE_Forcing
      PROCEDURE, PRIVATE :: Check_Write_Lock
      PROCEDURE, PRIVATE :: Manage_Write_Lock

      PROCEDURE :: Read_IFP_IPE_Forcing
      PROCEDURE :: Read_Tiros_IPE_Forcing

  END TYPE IPE_Forcing

CONTAINS

  SUBROUTINE Build_IPE_Forcing( forcing, parameters, mpi_layer, io, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ),          INTENT(out) :: forcing
    CLASS( IPE_Model_Parameters ), INTENT(in)  :: parameters
    CLASS( IPE_MPI_Layer ),        INTENT(in)  :: mpi_layer
    CLASS( COMIO_t ),              INTENT(in)  :: io
    INTEGER, OPTIONAL,             INTENT(out) :: rc

    ! Local
    INTEGER :: localrc, stat
    REAL(preC) :: dt

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    forcing % coupled = .false.

    dt            = parameters % solar_forcing_time_step

    forcing % dt            = parameters % solar_forcing_time_step
    forcing % start_time    = 0.0_prec
    forcing % current_time  = 0.0_prec
    forcing % current_index = 1

    ALLOCATE( forcing % emaps(1:maps_ipe_size(1),1:maps_ipe_size(2),1:maps_ipe_size(3)), &
              forcing % cmaps(1:maps_ipe_size(1),1:maps_ipe_size(2),1:maps_ipe_size(3)), &
              forcing % djspectra(1:n_flux_ipe,1:n_bands_ipe), &
              stat = stat )

    forcing % emaps     = 0.0_prec
    forcing % cmaps     = 0.0_prec
    forcing % djspectra = 0.0_prec
    CALL forcing % Read_IFP_IPE_Forcing( parameters, &
                                         mpi_layer, &
                                         io, &
                                         localrc )
    IF( ipe_error_check( localrc, msg="call to Read_IFP_IPE_Forcing failed", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    CALL forcing % Read_Tiros_IPE_Forcing( localrc )
    IF( ipe_error_check( localrc, msg="call to Read_Tiros_IPE_Forcing failed", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Build_IPE_Forcing

  SUBROUTINE Allocate_IFP_IPE_Forcing ( forcing, n_time_levels, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    INTEGER,              INTENT(in)    :: n_time_levels
    INTEGER, OPTIONAL,    INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    if ( allocated ( forcing % f107 ) ) call forcing % trash()

    ALLOCATE( forcing % f107(n_time_levels), &
              forcing % f107_81day_avg(n_time_levels), &
              forcing % kp(n_time_levels), &
              forcing % kp_1day_avg(n_time_levels), &
              forcing % ap(n_time_levels), &
              forcing % ap_1day_avg(n_time_levels), &
              forcing % nhemi_power(n_time_levels), &
              forcing % nhemi_power_index(n_time_levels), &
              forcing % shemi_power(n_time_levels), &
              forcing % shemi_power_index(n_time_levels), &
              forcing % solarwind_Bt(n_time_levels), &
              forcing % solarwind_angle(n_time_levels), &
              forcing % solarwind_velocity(n_time_levels), &
              forcing % solarwind_Bz(n_time_levels), &
              forcing % solarwind_By(n_time_levels), &
              forcing % solarwind_density(n_time_levels), &
              stat = stat )
    IF ( ipe_alloc_check(stat, msg="Failed to allocate forcing internal arrays", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN


  END SUBROUTINE Allocate_IFP_IPE_Forcing

  SUBROUTINE Trash_IPE_Forcing( forcing, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    INTEGER, OPTIONAL,    INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    DEALLOCATE( forcing % f107, &
                forcing % f107_81day_avg, &
                forcing % kp, &
                forcing % kp_1day_avg, &
                forcing % ap, &
                forcing % ap_1day_avg, &
                forcing % nhemi_power, &
                forcing % nhemi_power_index, &
                forcing % shemi_power, &
                forcing % shemi_power_index, &
                forcing % solarwind_Bt, &
                forcing % solarwind_By, &
                forcing % solarwind_angle, &
                forcing % solarwind_velocity, &
                forcing % solarwind_density, &
                forcing % solarwind_Bz, &
                stat = stat )
    IF ( ipe_dealloc_check( stat, msg="Unable to free up memory", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Trash_IPE_Forcing


  FUNCTION GetAP( forcing ) RESULT( AP )

    CLASS( IPE_Forcing ) :: forcing
    REAL(prec)           :: AP(1:7)

    ! Local
    INTEGER    :: i, im3hr, im6hr, im9hr, im12hr, im36hr

    i       = forcing % current_index
    im3hr   = MAX( 1, i - INT( 3.0_prec*3600.0_prec/forcing % dt ) )
    im6hr   = MAX( 1, i - INT( 6.0_prec*3600.0_prec/forcing % dt ) )
    im9hr   = MAX( 1, i - INT( 9.0_prec*3600.0_prec/forcing % dt ) )
    im12hr  = MAX( 1, i - INT(12.0_prec*3600.0_prec/forcing % dt ) )
    im36hr  = MAX( 1, i - INT(36.0_prec*3600.0_prec/forcing % dt ) )

    AP(1) = forcing % ap_1day_avg(i)
    AP(2) = forcing % ap(i)
    AP(3) = forcing % ap(im3hr)
    AP(4) = forcing % ap(im6hr)
    AP(5) = forcing % ap(im9hr)
    AP(6) = forcing % ap_1day_avg(im12hr)
    AP(7) = forcing % ap_1day_avg(im36hr)

  END FUNCTION GetAP


  FUNCTION GetKP( forcing ) RESULT( KP )

    CLASS( IPE_Forcing ) :: forcing
    REAL(prec)           :: KP(1:7)

    ! Local
    INTEGER    :: i, im3hr, im6hr, im9hr, im12hr, im36hr

    i       = forcing % current_index
    im3hr   = MAX( 1, i - INT( 3.0_prec*3600.0_prec/forcing % dt ) )
    im6hr   = MAX( 1, i - INT( 6.0_prec*3600.0_prec/forcing % dt ) )
    im9hr   = MAX( 1, i - INT( 9.0_prec*3600.0_prec/forcing % dt ) )
    im12hr  = MAX( 1, i - INT(12.0_prec*3600.0_prec/forcing % dt ) )
    im36hr  = MAX( 1, i - INT(36.0_prec*3600.0_prec/forcing % dt ) )

    KP(1) = forcing % kp_1day_avg(i)
    KP(2) = forcing % kp(i)
    KP(3) = forcing % kp(im3hr)
    KP(4) = forcing % kp(im6hr)
    KP(5) = forcing % kp(im9hr)
    KP(6) = forcing % kp_1day_avg(im12hr)
    KP(7) = forcing % kp_1day_avg(im36hr)

  END FUNCTION GetKP


  SUBROUTINE Update_Current_Index( forcing, parameters, mpi_layer, io, deltime, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ),          INTENT(inout) :: forcing
    CLASS( IPE_Model_Parameters ), INTENT(in)    :: parameters
    CLASS( IPE_MPI_Layer ),        INTENT(in)    :: mpi_layer
    CLASS( COMIO_t ),              INTENT(in)    :: io
    REAL(prec),                    INTENT(in)    :: deltime
    INTEGER,                       INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc
    if ( parameters % use_ifp_file .and. parameters % ifp_realtime_interval > 0 .and. &
         mod(INT( deltime ), parameters % ifp_realtime_interval) .eq. 0 ) then
      call forcing % Read_IFP_IPE_Forcing( parameters, &
                                           mpi_layer, &
                                           io, &
                                           localrc )
      IF( ipe_error_check( localrc, msg="call to Read_F107KP_IPE_Forcing failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    end if

    forcing % current_index = INT( deltime / real(forcing % ifp_interval) ) + forcing % ifp_skip + 1

  END SUBROUTINE Update_Current_Index

  SUBROUTINE Check_Write_Lock( forcing, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(in)  :: forcing
    INTEGER,              INTENT(out) :: rc

    character(len=22), parameter :: filename = "input_parameters.wlock"
    logical :: not_ready
    integer :: iostat

    rc = IPE_SUCCESS

    not_ready = .true.
    do while (not_ready)
      inquire(file=filename, exist=not_ready, iostat=iostat)
      if (not_ready) call sleep(1)
    end do

  END SUBROUTINE Check_Write_Lock

  SUBROUTINE Manage_Write_Lock( forcing, me, create, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(in)  :: forcing
    INTEGER,              INTENT(in)  :: me
    LOGICAL,              INTENT(in)  :: create
    INTEGER,              INTENT(out) :: rc

    character(len=26), parameter :: filename = "input_parameters.rlock.ipe"
    character(len=29)            :: lockfile
    integer, parameter           :: unit = 79
    integer                      :: localrc

    rc = IPE_SUCCESS

    write (lockfile, "(A22,I0.3)") filename, me
    open(unit, file=lockfile, status="replace", action="write", iostat=localrc)
      IF( ipe_error_check( localrc, msg="writing lockfile "//trim(lockfile)//" failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    if (create) then
      close(unit, iostat=localrc)
      IF( ipe_error_check( localrc, msg="closing lockfile "//trim(lockfile)//" failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    else
      close(unit, status="delete", iostat=localrc)
      IF( ipe_error_check( localrc, msg="destroying lockfile "//trim(lockfile)//" failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    end if

  END SUBROUTINE Manage_Write_Lock

  SUBROUTINE Read_IFP_IPE_Forcing( forcing, parameters, mpi_layer, io, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ),          INTENT(inout) :: forcing
    CLASS( IPE_Model_Parameters ), INTENT(in)    :: parameters
    CLASS( IPE_MPI_Layer ),        INTENT(in)    :: mpi_layer
    CLASS( COMIO_t ),              INTENT(in)    :: io
    INTEGER,                       INTENT(out)   :: rc

    ! Local
    integer, parameter      :: fmt =  COMIO_FMT_PNETCDF
    integer, pointer        :: dims(:)
    integer                 :: localrc

    rc = IPE_SUCCESS

    if ( parameters % use_ifp_file ) then
      call forcing % check_write_lock(localrc)
      IF( ipe_error_check( localrc, msg="call to Check_Write_Lock failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
      call forcing % manage_write_lock(mpi_layer % rank_id, .true., localrc)
      IF( ipe_error_check( localrc, msg="call to Manage_Write_Lock(create) failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      call io % open(parameters % ifp_file, "r")
      if (io % err % check(msg="Unable to open" // parameters % ifp_file, file=__FILE__,line=__LINE__)) return

      call io % description("skip", forcing % ifp_skip)
      call io % description("ifp_interval", forcing % ifp_interval)
      if (io % err % check(msg="Unable to description", file=__FILE__,line=__LINE__)) return

      call io % domain("f107", dims)
      if (io % err % check(msg="Unable to domain", file=__FILE__,line=__LINE__)) return

      call forcing % allocate_ifp_ipe_forcing(dims(1), localrc)
      IF( ipe_error_check( localrc, msg="call to Allocate_IFP_IPE_Forcing failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      call io % read("f107",  forcing % f107)
      call io % read("f107d", forcing % f107_81day_avg)
      call io % read("kp",    forcing % kp)
      call io % read("kpa",   forcing % kp_1day_avg)
      call io % read("ap",    forcing % ap)
      call io % read("apa",   forcing % ap_1day_avg)
      call io % read("nhp",   forcing % nhemi_power)
      call io % read("nhpi",  forcing % nhemi_power_index)
      call io % read("shp",   forcing % shemi_power)
      call io % read("shpi",  forcing % shemi_power_index)
      call io % read("swden", forcing % solarwind_density)
      call io % read("swang", forcing % solarwind_angle)
      call io % read("swvel", forcing % solarwind_velocity)
      call io % read("swbz",  forcing % solarwind_Bz)
      call io % read("swbt",  forcing % solarwind_Bt)
      call io % close()

      call forcing % manage_write_lock(mpi_layer % rank_id, .false., localrc)
      IF( ipe_error_check( localrc, msg="call to Manage_Write_Lock(destroy) failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    else ! fixed parameters
      call forcing % allocate_ifp_ipe_forcing(int( (parameters % end_time - parameters % start_time)/parameters % time_step ), localrc)
      IF( ipe_error_check( localrc, msg="call to Allocate_IFP_IPE_Forcing failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
      forcing % f107              = parameters % f107
      forcing % f107_81day_avg    = parameters % f107_81day_avg
      forcing % kp                = parameters % kp
      forcing % kp_1day_avg       = parameters % kp_1day_avg
      forcing % ap                = parameters % ap
      forcing % ap_1day_avg       = parameters % ap_1day_avg
      forcing % nhemi_power       = parameters % nhemi_power
      forcing % nhemi_power_index = parameters % nhemi_power_index
      forcing % shemi_power       = parameters % shemi_power
      forcing % shemi_power_index = parameters % shemi_power_index

      forcing % solarwind_angle    = parameters % solarwind_angle
      forcing % solarwind_velocity = parameters % solarwind_velocity
      forcing % solarwind_density  = parameters % solarwind_density
      forcing % solarwind_Bz       = parameters % solarwind_Bz
      forcing % solarwind_Bt       = sqrt( forcing % solarwind_By**2 +&
                                           forcing % solarwind_Bz**2 )
    end if
    ! By is always set by parameter for now
    forcing % solarwind_By       = parameters % solarwind_By
  END SUBROUTINE Read_IFP_IPE_Forcing


  SUBROUTINE Read_Tiros_IPE_Forcing( forcing, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    INTEGER,              INTENT(out)   :: rc

    ! Local
    INTEGER :: fUnit, iBand, irec, iostat
    CHARACTER(LEN=*), PARAMETER :: ion_filename     = "./ionprof"
    CHARACTER(LEN=*), PARAMETER :: spectra_filename = "./tiros_spectra"

    OPEN( UNIT = NewUnit(fUnit), &
          FILE = ion_filename, &
          FORM = 'FORMATTED', &
          STATUS = 'OLD', &
          ACTION = 'READ',&
          IOSTAT = iostat )
    IF ( ipe_iostatus_check( iostat, msg="Error opening tiros file "//ion_filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    READ (fUnit,'(1x,6E13.6)', IOSTAT = iostat) forcing % emaps
    IF ( ipe_iostatus_check( iostat, msg="Error reading emaps from "//ion_filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    READ (fUnit,'(1x,6E13.6)', IOSTAT = iostat) forcing % cmaps
    IF ( ipe_iostatus_check( iostat, msg="Error reading cmaps from "//ion_filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    CLOSE(fUnit, IOSTAT = iostat)
    IF ( ipe_iostatus_check( iostat, msg="Error closing tiros file "//ion_filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    OPEN( UNIT = NewUnit(fUnit), &
          FILE = spectra_filename, &
          FORM = 'FORMATTED', &
          STATUS = 'OLD', &
          ACTION = 'READ',&
          IOSTAT = iostat )
    IF ( ipe_iostatus_check( iostat, msg="Error opening tiros file "//spectra_filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    DO irec = 1, 3
      READ(fUnit,*, IOSTAT = iostat)
      IF ( ipe_iostatus_check( iostat, msg="Error advancing tiros file "//spectra_filename, &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    END DO

    DO iBand = 1, n_bands_ipe

      DO irec = 1, 2
        READ(fUnit,*, IOSTAT = iostat)
        IF ( ipe_iostatus_check( iostat, msg="Error advancing tiros file "//spectra_filename, &
          line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
      END DO

      READ(fUnit,'(1X,5E10.4)', IOSTAT = iostat) forcing % djspectra(1:n_flux_ipe,iBand)
      IF ( ipe_iostatus_check( iostat, msg="Error reading spectral data from file "//spectra_filename, &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      READ(fUnit,*, IOSTAT = iostat)
      IF ( ipe_iostatus_check( iostat, msg="Error advancing tiros file "//spectra_filename, &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ENDDO

    CLOSE(fUnit, IOSTAT = iostat)
    IF ( ipe_iostatus_check( iostat, msg="Error closing tiros file "//spectra_filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Read_Tiros_IPE_Forcing
!
END MODULE IPE_Forcing_Class
