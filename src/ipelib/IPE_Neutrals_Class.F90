MODULE IPE_Neutrals_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Grid_Class
  USE IPE_Time_Class
  USE IPE_Forcing_Class
  USE ipe_error_module

  ! MSIS
  USE physics_msis ! gtd7
  USE utils_constants, ONLY : msis_dp => dp

  IMPLICIT NONE

  TYPE IPE_Neutrals

    INTEGER :: nFluxTube, NLP, NMP
    INTEGER :: mp_low, mp_high

    REAL(prec), POINTER :: helium(:,:,:)             ! he
    REAL(prec), POINTER :: oxygen(:,:,:)             ! on
    REAL(prec), POINTER :: molecular_oxygen(:,:,:)   ! o2n
    REAL(prec), POINTER :: molecular_nitrogen(:,:,:) ! n2n
    REAL(prec), POINTER :: nitrogen(:,:,:)           ! n4s
    REAL(prec), POINTER :: hydrogen(:,:,:)           ! hn
    REAL(prec), POINTER :: temperature(:,:,:)
    REAL(prec), POINTER :: temperature_inf(:,:,:)
    REAL(prec), POINTER :: velocity_geographic(:,:,:,:)
    REAL(prec), POINTER :: velocity_apex(:,:,:,:)

    ! Interpolated fields
    REAL(prec), ALLOCATABLE :: geo_helium(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_oxygen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_molecular_oxygen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_molecular_nitrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_nitrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_hydrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_velocity(:,:,:,:)

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Neutrals
      PROCEDURE :: Trash => Trash_IPE_Neutrals

      PROCEDURE :: Update => Update_IPE_Neutrals

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: IPE_Neutrals_Empirical
      PROCEDURE, PRIVATE :: IPE_Neutrals_Extrapolate
      PROCEDURE, PRIVATE :: Geographic_to_Apex_Velocity

  END TYPE IPE_Neutrals


  CHARACTER(250), PARAMETER      :: hwm_path     = './'

  INTEGER,    PARAMETER, PRIVATE :: N_heights    = 72
  INTEGER,    PARAMETER, PRIVATE :: N_Latitudes  = 19
  INTEGER,    PARAMETER, PRIVATE :: N_Longitudes = 36
  REAL(prec), PARAMETER, PRIVATE :: small_power  = -3.0_prec
  REAL(prec), PARAMETER, PRIVATE :: small_number = 1.0e-03_prec
  REAL(prec), PARAMETER, PRIVATE :: min_density  = 1.0e-12_prec


CONTAINS

  SUBROUTINE Build_IPE_Neutrals( neutrals, nFluxTube, NLP, NMP, mp_low, mp_high, rc )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    INTEGER,               INTENT(in)    :: nFluxTube
    INTEGER,               INTENT(in)    :: NLP
    INTEGER,               INTENT(in)    :: NMP
    INTEGER,               INTENT(in)    :: mp_low, mp_high
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    neutrals % nFluxTube = nFluxTube
    neutrals % NLP       = NLP
    neutrals % NMP       = NMP
    neutrals % mp_low    = mp_low
    neutrals % mp_high   = mp_high

    ALLOCATE( neutrals % helium(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % oxygen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % molecular_oxygen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % molecular_nitrogen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % nitrogen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % hydrogen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % temperature(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % temperature_inf(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % velocity_geographic(1:3,1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % velocity_apex(1:3,nFluxTube,1:NLP,mp_low:mp_high), &
              stat = stat )
    IF ( ipe_alloc_check( stat, msg="Failed to allocate neutrals internal arrays", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    neutrals % helium              = 0.0_prec
    neutrals % oxygen              = 0.0_prec
    neutrals % molecular_oxygen    = 0.0_prec
    neutrals % molecular_nitrogen  = 0.0_prec
    neutrals % nitrogen            = 0.0_prec
    neutrals % hydrogen            = 0.0_prec
    neutrals % temperature         = 0.0_prec
    neutrals % temperature_inf     = 0.0_prec
    neutrals % velocity_geographic = 0.0_prec
    neutrals % velocity_apex       = 0.0_prec

    ALLOCATE( neutrals % geo_helium(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_oxygen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_molecular_oxygen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_molecular_nitrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_nitrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo),&
              neutrals % geo_hydrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo),&
              neutrals % geo_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_velocity(1:3,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              stat = stat )
    IF ( ipe_alloc_check( stat, msg="Failed to allocate neutrals internal geo arrays", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    neutrals % geo_helium             = 0.0_prec
    neutrals % geo_oxygen             = 0.0_prec
    neutrals % geo_molecular_oxygen   = 0.0_prec
    neutrals % geo_molecular_nitrogen = 0.0_prec
    neutrals % geo_nitrogen           = 0.0_prec
    neutrals % geo_hydrogen           = 0.0_prec
    neutrals % geo_temperature        = 0.0_prec
    neutrals % geo_velocity           = 0.0_prec

  END SUBROUTINE Build_IPE_Neutrals


  SUBROUTINE Trash_IPE_Neutrals( neutrals, rc )

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    IF ( ASSOCIATED( neutrals % helium              ) ) &
         DEALLOCATE( neutrals % helium              , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % hydrogen            ) ) &
         DEALLOCATE( neutrals % hydrogen            , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % molecular_nitrogen  ) ) &
         DEALLOCATE( neutrals % molecular_nitrogen  , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % molecular_oxygen    ) ) &
         DEALLOCATE( neutrals % molecular_oxygen    , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % nitrogen            ) ) &
         DEALLOCATE( neutrals % nitrogen            , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % oxygen              ) ) &
         DEALLOCATE( neutrals % oxygen              , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % temperature         ) ) &
         DEALLOCATE( neutrals % temperature         , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % temperature_inf     ) ) &
         DEALLOCATE( neutrals % temperature_inf     , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % velocity_apex       ) ) &
         DEALLOCATE( neutrals % velocity_apex       , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % velocity_geographic ) ) &
         DEALLOCATE( neutrals % velocity_geographic , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN


    DEALLOCATE( neutrals % geo_helium, &
                 neutrals % geo_oxygen, &
                 neutrals % geo_molecular_oxygen, &
                 neutrals % geo_molecular_nitrogen, &
                 neutrals % geo_nitrogen, &
                 neutrals % geo_hydrogen, &
                 neutrals % geo_temperature, &
                 neutrals % geo_velocity, &
                 stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Trash_IPE_Neutrals


  SUBROUTINE Update_IPE_Neutrals( neutrals, params, grid, time, forcing, mpi_layer, vertical_wind_limit, rc )

    CLASS( IPE_Neutrals         ), INTENT(inout) :: neutrals
    TYPE ( IPE_Model_Parameters ), INTENT(in   ) :: params
    TYPE ( IPE_Grid             ), INTENT(in   ) :: grid
    TYPE ( IPE_Time             ), INTENT(in   ) :: time
    TYPE ( IPE_Forcing          ), INTENT(in   ) :: forcing
    TYPE ( IPE_MPI_Layer        ), INTENT(in   ) :: mpi_layer
    REAL ( prec ),                 INTENT(in   ) :: vertical_wind_limit

    INTEGER, OPTIONAL,             INTENT(out  ) :: rc

    ! Local
    LOGICAL :: msis_switch
    INTEGER :: localrc

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    msis_switch = mod(time % elapsed_sec,params % msis_time_step) == 0.0

    IF ( msis_switch .and. (time % elapsed_sec > 0._prec .or. .NOT. params % read_apex_neutrals) ) THEN
      IF( mpi_layer % rank_id == 0 ) write(6,*) 'Calling MSIS ', int(time % elapsed_sec / 60), ' Mins UT'
      CALL neutrals % IPE_Neutrals_Empirical( grid, time, forcing, rc=localrc )
      IF ( ipe_error_check( localrc, msg="call to IPE_Neutrals_Empirical failed", rc=rc ) ) RETURN
    ENDIF

    CALL neutrals % IPE_Neutrals_Extrapolate( grid, forcing )

    CALL neutrals % Geographic_to_Apex_Velocity( grid, vertical_wind_limit )

  END SUBROUTINE Update_IPE_Neutrals


  SUBROUTINE IPE_Neutrals_Extrapolate( neutrals, grid, forcing )

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE ( IPE_Grid     ), INTENT(in   ) :: grid
    TYPE ( IPE_Forcing  ), INTENT(in   ) :: forcing

    ! Local
    INTEGER    :: kp, kpp, kpStart, kpStep, kpStop, lp, mp
    REAL(prec) :: r, w
    REAL(prec) :: t( grid % mp_low:grid % mp_high )

    REAL(prec), PARAMETER :: fac = 2.0_prec * G0 / GSCON

    ! -- NOTE: requires correct northern_top_index and southern_top_index arrays

    DO lp = 1, grid % NLP

      ! loop over hemispheres, starting with northern one
      kpStart = grid % northern_top_index(lp)
      kpStop  = grid % flux_tube_midpoint(lp)

      DO kpStep = 1, -1, -2

        DO kp = kpStart + kpStep, kpStop, kpStep

          ! extend temperature
          neutrals % temperature(kp,lp,grid % mp_low:grid % mp_high) = &
            neutrals % temperature(kpStart,lp,grid % mp_low:grid % mp_high)

          ! extend winds
          neutrals % velocity_geographic(:,kp,lp,grid % mp_low:grid % mp_high) = &
            neutrals % velocity_geographic(:,kpStart,lp,grid % mp_low:grid % mp_high)

          r = 0.0_prec
          t = 0.0_prec
          DO kpp = kp - kpStep, kp, kpStep
            w = 1.0_prec + grid % altitude(kpp, lp) / earth_radius
            t = t + w * w * neutrals % temperature(kpp,lp,:)
            r = grid % altitude(kpp, lp) - r
          ENDDO
          t = -fac * r / t

          kpp = kp - kpStep
          ! extrapolating atomic oxygen
          neutrals % oxygen(kp,lp,grid % mp_low:grid % mp_high) = &
            max( min_density, neutrals % oxygen(kpp,lp,grid % mp_low:grid % mp_high) * exp( O_mass * t ) )
          ! extrapolating molecular oxygen
          neutrals % molecular_oxygen(kp,lp,grid % mp_low:grid % mp_high) = &
            max( min_density, neutrals % molecular_oxygen(kpp,lp,grid % mp_low:grid % mp_high) * exp( O2_mass * t ) )
          ! extrapolating molecular nitrogen
          neutrals % molecular_nitrogen(kp,lp,grid % mp_low:grid % mp_high) = &
            max( min_density, neutrals % molecular_nitrogen(kpp,lp,grid % mp_low:grid % mp_high) * exp( N2_mass * t ) )

        ENDDO
        ! preparing for the southern hemisphere
        kpStart = grid % southern_top_index(lp)
        kpStop  = kpStop + 1
      ENDDO
    ENDDO

    ! set exospheric temperature
    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        neutrals % temperature_inf(1:grid % flux_tube_midpoint(lp),lp,mp) = &
          neutrals % temperature(grid % northern_top_index(lp),lp,mp)
        neutrals % temperature_inf(grid % flux_tube_midpoint(lp)+1:grid % flux_tube_max(lp),lp,mp) = &
          neutrals % temperature(grid % southern_top_index(lp),lp,mp)
      ENDDO
    ENDDO

  END SUBROUTINE IPE_Neutrals_Extrapolate


  SUBROUTINE IPE_Neutrals_Empirical( neutrals, grid, time, forcing, rc )
  !
  ! Usage :
  !
  !   CALL neutrals % Update( grid, utime, year, day, f107d, f107a, ap )
  ! ================================================================================================== !

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ),      INTENT(in)    :: grid
    TYPE( IPE_Time ),      INTENT(in)    :: time
    TYPE( IPE_Forcing ),   INTENT(in)    :: forcing
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER    :: localrc
    INTEGER    :: kp, lp, mp, day
    REAL(prec) :: geo_alt, geo_lat, geo_lon, utime
    REAL(prec), DIMENSION(7) :: AP

    INTEGER(4)            :: iyd
    REAL(4)               :: hwm_sec, hwm_f107d, hwm_f107a, hwm_alt, hwm_lat, hwm_lon
    REAL(4), DIMENSION(2) :: hwm_ap, w

    REAL(msis_dp)               :: msis_alt, msis_f107d, msis_f107a, msis_lat, msis_lon, msis_sec, msis_stl
    REAL(msis_dp), DIMENSION(7) :: msis_ap
    REAL(msis_dp), DIMENSION(2) :: temperatures
    REAL(msis_dp), DIMENSION(9) :: densities

    INTEGER, PARAMETER    :: msis_mass = 48


    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    AP = forcing % GetAP( )

    w         = 0.0
    hwm_ap    = REAL(AP(1:2), KIND=4)
    hwm_f107a = REAL(forcing % f107_81day_avg( forcing % current_index ), KIND=4)
    hwm_f107d = REAL(forcing % f107( forcing % current_index ),           KIND=4)

    msis_ap      = REAL(AP,    KIND=msis_dp)
    msis_f107a   = REAL(forcing % f107_81day_avg( forcing % current_index ), KIND=msis_dp)
    msis_f107d   = REAL(forcing % f107( forcing % current_index ),           KIND=msis_dp)
    densities    = 0.0_msis_dp
    temperatures = 0.0_msis_dp

    iyd = 99000 + time % day_of_year    ! Input, year and day as yyddd

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)
          geo_lon = rtd * grid % longitude(kp,lp,mp)
          geo_lat = 90.0 - rtd * grid % colatitude(kp,lp,mp)
          geo_alt = m_to_km * grid % altitude(kp,lp)

          ! -- horizontal & vertical wind
          iF ( .NOT. forcing % coupled ) THEN

            hwm_sec = REAL(time % utime, KIND=4)
            hwm_alt = geo_alt
            hwm_lat = geo_lat
            hwm_lon = geo_lon
            w       = 0.0

            call hwm14( iyd,       &    ! Input, year and day as yyddd
                        hwm_sec,   &    ! Input, universal time ( sec )
                        hwm_alt,   &    ! Input, altitude ( km )
                        hwm_lat,   &    ! Input, geodetic latitude ( degrees )
                        hwm_lon,   &    ! Input, geodetic longitude ( degrees )
                        0.0,       &    ! Input, local apparent solar time ( hrs )[ not used ]
                        hwm_f107a, &    ! Input, 3 month average of f10.7 flux [ not used ]
                        hwm_f107d, &    ! Input, daily average of f10.7 flux for the previous day [ not used ]
                        hwm_ap,    &    ! Input, magnetic index ( daily ), current 3hr ap index
                        hwm_path,  &    ! Input, default datafile path
                        w,         &    ! Ouput, neutral wind velocity meridional-northwards(1) and zonal-eastwards(2) components
                        localrc )       ! Ouput, return code
            IF ( ipe_error_check( localrc, msg="call to hwm14 failed", &
              line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

            neutrals % velocity_geographic(1,kp,lp,mp) = w(2)
            neutrals % velocity_geographic(2,kp,lp,mp) = w(1)
            neutrals % velocity_geographic(3,kp,lp,mp) = 0.0_prec

          ENDIF

          ! -- composition & temperature
          msis_sec     = REAL(time % utime, KIND=msis_dp)
          msis_alt     = geo_alt
          msis_lat     = geo_lat
          msis_lon     = geo_lon
          msis_stl     = REAL(time % utime / 3600.0_prec + geo_lon / 15.0_prec, KIND=msis_dp)
          densities    = 0.0_msis_dp
          temperatures = 0.0_msis_dp

          call gtd7( iyd,         &    ! Input, year and day as yyddd
                     msis_sec,    &    ! Input, universal time ( sec )
                     msis_alt,    &    ! Input, altitude ( km )
                     msis_lat,    &    ! Input, geodetic latitude ( degrees )
                     msis_lon,    &    ! Input, geodetic longitude ( degrees )
                     msis_stl,    &    ! Input, local apparent solar time ( hrs )
                     msis_f107a,  &    ! Input, 3 month average of f10.7 flux
                     msis_f107d,  &    ! Input, daily average of f10.7 flux for the previous day
                     msis_ap,     &    ! Input, magnetic index ( daily ), current, 3,6,9hrs prior 3hr ap index, 12-33 hr prior ap average, 36-57 hr prior ap average
                     msis_mass,   &    ! Mass number ( see src/msis/physics_msis.f90 for more details )
                     densities,   &    ! Ouput, neutral densities in cm-3
                     temperatures )    ! Output, exospheric temperature and temperature at altitude

          IF ( .NOT. forcing % coupled ) THEN

            neutrals % temperature_inf(kp,lp,mp) = temperatures(1)
            neutrals % temperature(kp,lp,mp)     = temperatures(2)

            neutrals % oxygen(kp,lp,mp)             = cm_3_to_m_3 * densities(2)
            neutrals % molecular_nitrogen(kp,lp,mp) = cm_3_to_m_3 * densities(3)
            neutrals % molecular_oxygen(kp,lp,mp)   = cm_3_to_m_3 * densities(4)

          ENDIF

          neutrals % helium(kp,lp,mp)   = cm_3_to_m_3 * densities(1)
          neutrals % hydrogen(kp,lp,mp) = cm_3_to_m_3 * densities(7)
          neutrals % nitrogen(kp,lp,mp) = cm_3_to_m_3 * densities(8)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE IPE_Neutrals_Empirical


  SUBROUTINE Geographic_to_Apex_Velocity( neutrals, grid, vertical_wind_limit )

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ),      INTENT(in)    :: grid
    REAL(prec), INTENT(in) :: vertical_wind_limit

    ! Local
    INTEGER    :: i_D_vec, kp, lp, mp
    REAL(prec) :: dotprod
    REAL(prec) :: neutral_vertical_velocity

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)

          DO i_D_vec = 1, 3 ! D1, D2 and D3

             dotprod = grid % apex_d_vectors(1,i_D_vec,kp,lp,mp)*grid % apex_d_vectors(1,i_D_vec,kp,lp,mp) + &
                       grid % apex_d_vectors(2,i_D_vec,kp,lp,mp)*grid % apex_d_vectors(2,i_D_vec,kp,lp,mp) + &
                       grid % apex_d_vectors(3,i_D_vec,kp,lp,mp)*grid % apex_d_vectors(3,i_D_vec,kp,lp,mp)

             IF ( dotprod > 0.0_prec ) THEN

               if (neutrals % velocity_geographic(3,kp,lp,mp).gt.vertical_wind_limit) then
                 neutral_vertical_velocity = vertical_wind_limit
               else if (neutrals % velocity_geographic(3,kp,lp,mp).lt.0.0 - vertical_wind_limit) then
                 neutral_vertical_velocity = 0.0 - vertical_wind_limit
               else
                 neutral_vertical_velocity = neutrals % velocity_geographic(3,kp,lp,mp)
               endif

               neutrals % velocity_apex(i_D_vec,kp,lp,mp) = &
                 ( grid % apex_d_vectors(1,i_D_vec,kp,lp,mp)*neutrals % velocity_geographic(1,kp,lp,mp) + &
                   grid % apex_d_vectors(2,i_D_vec,kp,lp,mp)*neutrals % velocity_geographic(2,kp,lp,mp) + &
                   grid % apex_d_vectors(3,i_D_vec,kp,lp,mp)*neutral_vertical_velocity ) &
                 / SQRT( dotprod )

             ELSE

                neutrals % velocity_apex(i_D_vec,kp,lp,mp) = 0.0_prec

             END IF

           ENDDO

           ! dbg20110131 : the midpoint values become NaN otherwise because
           ! of inappropriate D1/3 values...
           IF ( lp >= 1 .AND. lp <= 6 )THEN
             IF( kp == grid % flux_tube_midpoint(lp) )THEN
               neutrals % velocity_apex(1:3,kp,lp,mp) = neutrals % velocity_apex(1:3,kp-1,lp,mp)
             ENDIF
           ENDIF

        END DO
      END DO
    END DO

  END SUBROUTINE Geographic_to_Apex_Velocity


END MODULE IPE_Neutrals_Class
