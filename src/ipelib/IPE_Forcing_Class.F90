MODULE IPE_Forcing_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines
USE IPE_Model_Parameters_Class
USE ipe_error_module

IMPLICIT NONE

  TYPE IPE_Forcing

    ! forcing % n_time_levels = f107_kp_data_size    = f107_kp_size
    LOGICAL                 :: coupled

    INTEGER                 :: n_time_levels
    REAL(prec)              :: dt
    REAL(prec), ALLOCATABLE :: time(:)
    REAL(prec)              :: current_time
    INTEGER                 :: current_index
    INTEGER                 :: max_read_index
    !
    REAL(prec), ALLOCATABLE :: f107(:)
    INTEGER, ALLOCATABLE    :: f107_flag(:)
    REAL(prec), ALLOCATABLE :: f107_81day_avg(:)
    REAL(prec), ALLOCATABLE :: kp(:)
    INTEGER, ALLOCATABLE    :: kp_flag(:)
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

      PROCEDURE :: Read_F107KP_IPE_Forcing
      PROCEDURE :: Read_Tiros_IPE_Forcing

      PROCEDURE, PRIVATE :: Estimate_AP_from_KP

  END TYPE IPE_Forcing

CONTAINS

  SUBROUTINE Build_IPE_Forcing( forcing, params, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ),          INTENT(out) :: forcing
    CLASS( IPE_Model_Parameters ), INTENT(in)  :: params
    INTEGER, OPTIONAL,             INTENT(out) :: rc

    ! Local
    INTEGER :: n_time_levels, localrc, stat
    REAL(preC) :: dt

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    forcing % coupled = .false.

    n_time_levels = params % f107_kp_size
    dt            = params % solar_forcing_time_step

    forcing % n_time_levels = params % f107_kp_size + params % f107_kp_read_in_start
    forcing % dt            = params % solar_forcing_time_step
    forcing % current_time  = 0.0_prec
    forcing % current_index = 1

    ALLOCATE( forcing % time(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % f107(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % f107_flag(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % f107_81day_avg(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % kp(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % kp_flag(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % kp_1day_avg(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % ap(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % ap_1day_avg(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % nhemi_power(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % nhemi_power_index(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % shemi_power(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % shemi_power_index(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % solarwind_Bt(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % solarwind_angle(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % solarwind_velocity(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % solarwind_Bz(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % solarwind_By(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % solarwind_density(params % f107_kp_read_in_start+1:forcing % n_time_levels), &
              forcing % emaps(1:maps_ipe_size(1),1:maps_ipe_size(2),1:maps_ipe_size(3)), &
              forcing % cmaps(1:maps_ipe_size(1),1:maps_ipe_size(2),1:maps_ipe_size(3)), &
              forcing % djspectra(1:n_flux_ipe,1:n_bands_ipe), &
              stat = stat )
    IF ( ipe_alloc_check(stat, msg="Failed to allocate forcing internal arrays", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! Default settings
    forcing % f107              = params % f107
    forcing % f107_flag         = params % f107_flag
    forcing % f107_81day_avg    = params % f107_81day_avg
    forcing % kp                = params % kp
    forcing % kp_flag           = params % kp_flag
    forcing % kp_1day_avg       = params % kp_1day_avg
    forcing % ap                = params % ap
    forcing % ap_1day_avg       = params % ap_1day_avg
    forcing % nhemi_power       = params % nhemi_power
    forcing % nhemi_power_index = params % nhemi_power_index
    forcing % shemi_power       = params % shemi_power
    forcing % shemi_power_index = params % shemi_power_index

    forcing % solarwind_angle    = params % solarwind_angle
    forcing % solarwind_velocity = params % solarwind_velocity
    forcing % solarwind_density  = params % solarwind_density
    forcing % solarwind_Bz       = params % solarwind_Bz
    forcing % solarwind_By       = params % solarwind_By
    forcing % solarwind_Bt       = sqrt( forcing % solarwind_By**2 +&
                                         forcing % solarwind_Bz**2 )

    forcing % emaps     = 0.0_prec
    forcing % cmaps     = 0.0_prec
    forcing % djspectra = 0.0_prec
    IF( params % use_f107_kp_file )THEN
      CALL forcing % Read_F107KP_IPE_Forcing( params % f107_kp_file,          &
                                              params % f107_kp_data_size,     &
                                              params % f107_kp_read_in_start+1, &
                                              0,                              &
                                              params % f107_kp_realtime_interval < 0, &
                                              localrc )
      IF( ipe_error_check( localrc, msg="call to Read_F107KP_IPE_Forcing failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    ENDIF

    CALL forcing % Estimate_AP_from_KP( params % f107_kp_read_in_start+1 )

    CALL forcing % Read_Tiros_IPE_Forcing( localrc )
    IF( ipe_error_check( localrc, msg="call to Read_Tiros_IPE_Forcing failed", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Build_IPE_Forcing


  SUBROUTINE Trash_IPE_Forcing( forcing, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    INTEGER, OPTIONAL,    INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    DEALLOCATE( forcing % time, &
                forcing % f107, &
                forcing % f107_flag, &
                forcing % f107_81day_avg, &
                forcing % kp, &
                forcing % kp_flag, &
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


  SUBROUTINE Update_Current_Index( forcing, params, deltime, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ),          INTENT(inout) :: forcing
    CLASS( IPE_Model_Parameters ), INTENT(in)    :: params
    REAL(prec),                    INTENT(in)    :: deltime
    INTEGER,                       INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc

    rc = IPE_SUCCESS

    forcing % current_index = INT( deltime / real(params % f107_kp_interval) ) + &
                                      1 + params % f107_kp_skip_size

    if ( params % use_f107_kp_file .and. forcing % current_index > forcing % max_read_index &
             .and. params % f107_kp_realtime_interval > 0 ) then
      call forcing % read_f107kp_ipe_forcing( params % f107_kp_file,   &
                                              params % f107_kp_realtime_interval, &
                                              forcing % current_index, &
                                              forcing % max_read_index, &
                                              params % f107_kp_realtime_interval < 0, &
                                              localrc )
      IF( ipe_error_check( localrc, msg="call to Read_F107KP_IPE_Forcing failed", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    end if
    CALL forcing % Estimate_AP_from_KP( forcing % current_index )

  END SUBROUTINE Update_Current_index


  SUBROUTINE Read_F107KP_IPE_Forcing( forcing, filename, data_size, read_in_start, read_in_skip, fill, rc )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    CHARACTER(*),         INTENT(in)    :: filename
    INTEGER,              INTENT(in)    :: data_size
    INTEGER,              INTENT(in)    :: read_in_start
    INTEGER,              INTENT(in)    :: read_in_skip
    LOGICAL,              INTENT(in)    :: fill
    INTEGER,              INTENT(out)   :: rc

    ! Local
    INTEGER       :: fUnit
    INTEGER       :: i, iostat, read_in_size
    CHARACTER(20) :: date_work

    rc = IPE_SUCCESS

    OPEN( UNIT   = NewUnit(fUnit), &
          FILE   = TRIM(filename), &
          FORM   = 'FORMATTED',    &
          ACTION = 'READ',         &
          STATUS = 'OLD' ,         &
          IOSTAT = iostat )
    IF ( ipe_iostatus_check( iostat, msg="Error opening forcing file "//filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! Skip over the header information
    DO i = 1, 5

      READ(fUnit, *, IOSTAT = iostat )
      IF ( ipe_iostatus_check( iostat, msg="Error advancing forcing file "//filename, &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    END DO

    DO i = 1, read_in_skip

      READ(fUnit, *, IOSTAT = iostat)
      IF ( ipe_iostatus_check( iostat, msg="Error advancing forcing file "//filename, &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    END DO

    read_in_size = MIN(forcing % n_time_levels, data_size)
    DO i = read_in_start, read_in_size + read_in_start - 1
!      write(6,*) "reading",i
      READ(fUnit, *, IOSTAT = iostat) date_work, &
                                     forcing % f107(i), &
                                     forcing % kp(i), &
                                     forcing % f107_flag(i), &
                                     forcing % kp_flag(i), &
                                     forcing % f107_81day_avg(i), &
                                     forcing % kp_1day_avg(i), &
                                     forcing % nhemi_power(i), &
                                     forcing % nhemi_power_index(i), &
                                     forcing % shemi_power(i), &
                                     forcing % shemi_power_index(i), &
                                     forcing % solarwind_Bt(i), &
                                     forcing % solarwind_angle(i), &
                                     forcing % solarwind_velocity(i), &
                                     forcing % solarwind_Bz(i), &
                                     forcing % solarwind_density(i)
      IF ( ipe_iostatus_check( iostat, msg="Error reading forcing file "//filename, &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    END DO

    CLOSE(fUnit, IOSTAT = iostat)
    IF ( ipe_iostatus_check( iostat, msg="Error closing forcing file "//filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    forcing % max_read_index = read_in_size + read_in_start - 1
    if ( fill ) then
      DO i = read_in_start + read_in_size, forcing % n_time_levels
        forcing % f107(i)               = forcing % f107(read_in_size)
        forcing % f107_81day_avg(i)     = forcing % f107_81day_avg(read_in_size)
        forcing % kp(i)                 = forcing % kp(read_in_size)
        forcing % kp_1day_avg(i)        = forcing % kp_1day_avg(read_in_size)
        forcing % nhemi_power(i)        = forcing % nhemi_power(read_in_size)
        forcing % nhemi_power_index(i)  = forcing % nhemi_power_index(read_in_size)
        forcing % shemi_power(i)        = forcing % shemi_power(read_in_size)
        forcing % shemi_power_index(i)  = forcing % shemi_power_index(read_in_size)
        forcing % solarwind_Bt(i)       = forcing % solarwind_Bt(read_in_size)
        forcing % solarwind_angle(i)    = forcing % solarwind_angle(read_in_size)
        forcing % solarwind_velocity(i) = forcing % solarwind_velocity(read_in_size)
        forcing % solarwind_Bz(i)       = forcing % solarwind_Bz(read_in_size)
        forcing % solarwind_density(i)  = forcing % solarwind_density(read_in_size)

      END DO
    end if

  END SUBROUTINE Read_F107KP_IPE_Forcing


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

  SUBROUTINE Estimate_AP_from_KP( forcing, read_in_start )

    IMPLICIT NONE

    CLASS( IPE_Forcing ), INTENT(inout) :: forcing
    INTEGER,              INTENT(in   ) :: read_in_start

    ! Local
    REAL(prec) :: lookup, remainder
    INTEGER    :: i
    INTEGER, PARAMETER  :: table(1:29) = (/  0,   2,   3, & ! 0-0.67
                                             4,   5,   6, & ! 1-1.67
                                             7,   9,  12, & ! 2-2.67
                                             15,  18,  22, & ! 3-3.67
                                             27,  32,  39, & ! 4-4.67
                                             48,  56,  67, & ! 5-5.67
                                             80,  94, 111, & ! 6-6.67
                                             132, 154, 179, & ! 7-7.67
                                             207, 236, 300, & ! 8-8.67
                                             400, 999/)       ! 9-dumm
    DO i=read_in_start, forcing % n_time_levels

      lookup    = forcing % kp(i)*3.0_prec + 1.0_prec
      remainder = lookup - INT(lookup)
      forcing % ap(i) = (1.0_prec - remainder) * table(INT(lookup)) + remainder *table(INT(lookup)+1)

      lookup    = forcing % kp_1day_avg(i)*3.0_prec + 1.0_prec
      remainder = lookup - INT(lookup)
      forcing % ap_1day_avg(i) = (1.0_prec - remainder) * table(INT(lookup)) + remainder *table(INT(lookup)+1)

    END DO

  END SUBROUTINE Estimate_AP_from_KP
!
END MODULE IPE_Forcing_Class
