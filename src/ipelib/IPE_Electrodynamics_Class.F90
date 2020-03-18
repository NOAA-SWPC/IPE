#include "IPE_Macros.inc"

MODULE IPE_Electrodynamics_Class


USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines
USE IPE_MPI_Layer_Class
USE IPE_Grid_Class
USE IPE_Forcing_Class
USE IPE_Time_Class
USE IPE_Plasma_Class

USE efield_ipe
USE dynamo_module

#ifdef HAVE_NETCDF
USE netcdf
#endif
!
!! For the TIEGCM Wrapper
!USE module_init_cons
!USE cons_module
!USE magfield_module
!USE dynamo_module
!USE module_sub_dynamo
!USE module_highlat
!USE module_init_heelis




IMPLICIT NONE

  TYPE IPE_Electrodynamics
    INTEGER    :: nFluxTube, NLP, NMP
    INTEGER    :: mp_low, mp_high, mp_halo
    REAL(prec), ALLOCATABLE :: electric_potential(:,:)
    REAL(prec), ALLOCATABLE :: electric_potential2(:,:)
    REAL(prec), ALLOCATABLE :: mhd_electric_potential(:,:)
    REAL(prec), ALLOCATABLE :: electric_field(:,:,:)
    REAL(prec), ALLOCATABLE :: v_ExB_geographic(:,:,:,:)  ! "zonal" and "meridional" direction on the geographic grid
    REAL(prec), ALLOCATABLE :: v_ExB_apex(:,:,:) ! "zonal" and "meridional" direction ( VEXBth, VEXBe ) on the apex grid

    REAL(prec), ALLOCATABLE, PRIVATE :: lat_interp_weights(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude
    REAL(prec), ALLOCATABLE, PRIVATE :: lon_interp_weights(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude
    INTEGER, ALLOCATABLE, PRIVATE    :: lat_interp_index(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude
    INTEGER, ALLOCATABLE, PRIVATE    :: lon_interp_index(:,:) ! Weights for interpolating from magnetic longitude to ipe longitude

    ! Inputs for the potential solver (on the dynamo grid)
    REAL(prec), ALLOCATABLE :: hall_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: pedersen_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: b_parallel_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: neutral_apex_velocity(:,:,:) ! Components are apex directions


    ! Geographic interpolated attributes
    REAL(prec), ALLOCATABLE :: geo_electric_potential(:,:)
    REAL(prec), ALLOCATABLE :: geo_mhd_electric_potential(:,:)
    REAL(prec), ALLOCATABLE :: geo_v_ExB_geographic(:,:,:)    ! ExB transport velocity with geographic components
    REAL(prec), ALLOCATABLE :: geo_hall_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: geo_pedersen_conductivity(:,:)
    REAL(prec), ALLOCATABLE :: geo_b_parallel_conductivity(:,:)


    CONTAINS

      PROCEDURE :: Build => Build_IPE_Electrodynamics
      PROCEDURE :: Trash => Trash_IPE_Electrodynamics

      PROCEDURE :: Update => Update_IPE_Electrodynamics
      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Electrodynamics

      PROCEDURE :: Read_Geospace_Potential

      PROCEDURE :: Read_MHD_Potential
      PROCEDURE :: Write_MHD_Potential
      PROCEDURE :: Interpolate_Geospace_to_MHDpotential

      PROCEDURE, PRIVATE :: Empirical_E_Field_Wrapper
      PROCEDURE, PRIVATE :: Dynamo_Wrapper
      PROCEDURE, PRIVATE :: Regrid_Potential
      PROCEDURE, PRIVATE :: Calculate_Potential_Gradient
      PROCEDURE, PRIVATE :: Calculate_ExB_Velocity

  END TYPE IPE_Electrodynamics

  REAL(prec), PARAMETER, PRIVATE :: fillValue = -999999.9_prec
  LOGICAL, PRIVATE :: dynamo_efield

  REAL(prec), PRIVATE :: theta90_rad(0:nmlat)
  REAL(prec), PRIVATE, ALLOCATABLE :: geospace_latitude(:), geospace_longitude(:), geospace_potential(:,:)

  INTEGER, PRIVATE :: n_lat_geospace, n_lon_geospace

CONTAINS


  SUBROUTINE Build_IPE_Electrodynamics( eldyn, nFluxTube, NLP, NMP, dynamo, mp_low, mp_high, halo )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(out) :: eldyn
    INTEGER, INTENT(in)                       :: nFluxTube
    INTEGER, INTENT(in)                       :: NLP
    INTEGER, INTENT(in)                       :: NMP
    LOGICAL, INTENT(IN)                       :: dynamo
    INTEGER, INTENT(in)                       :: mp_low, mp_high, halo
    ! Local
    INTEGER :: j
    REAL(prec) :: theta130_rad

      eldyn % nFluxTube = nFluxTube
      eldyn % NLP       = NLP
      eldyn % NMP       = NMP
      eldyn % mp_low    = mp_low
      eldyn % mp_high   = mp_high
      eldyn % mp_halo   = halo


      ALLOCATE( eldyn % electric_potential(1:NLP,mp_low-halo:mp_high+halo), &
                eldyn % electric_potential2(1:NLP,mp_low-halo:mp_high+halo), &
                eldyn % mhd_electric_potential(1:NLP,mp_low-halo:mp_high+halo), &
                eldyn % electric_field(1:2,1:NLP,mp_low:mp_high), &
                eldyn % v_ExB_geographic(1:3,1:nFluxTube,1:NLP,mp_low:mp_high), &
                eldyn % v_ExB_apex(1:3,1:NLP,mp_low:mp_high), &
                eldyn % hall_conductivity(1:NLP,mp_low:mp_high), &
                eldyn % pedersen_conductivity(1:NLP,mp_low:mp_high), &
                eldyn % b_parallel_conductivity(1:NLP,mp_low:mp_high), &
                eldyn % neutral_apex_velocity(1:3,1:NLP,mp_low:mp_high), &
                eldyn % lat_interp_weights(1:2,1:NLP), &
                eldyn % lon_interp_weights(1:2,mp_low-halo:mp_high+halo), &
                eldyn % lat_interp_index(1:2,1:NLP), &
                eldyn % lon_interp_index(1:2,mp_low-halo:mp_high+halo) )

      eldyn % electric_potential      = 0.0_prec
      eldyn % electric_potential2     = 0.0_prec
      eldyn % mhd_electric_potential  = 0.0_prec
      eldyn % electric_field          = 0.0_prec
      eldyn % v_ExB_geographic        = 0.0_prec
      eldyn % v_ExB_apex              = 0.0_prec
      eldyn % hall_conductivity       = 0.0_prec
      eldyn % pedersen_conductivity   = 0.0_prec
      eldyn % b_parallel_conductivity = 0.0_prec
      eldyn % neutral_apex_velocity   = 0.0_prec
      eldyn % lat_interp_weights      = 0.0_prec
      eldyn % lon_interp_weights      = 0.0_prec
      eldyn % lat_interp_index        = 0
      eldyn % lon_interp_index        = 0


      ALLOCATE( eldyn % geo_electric_potential(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_mhd_electric_potential(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_v_ExB_geographic(1:3,1:nlon_geo,1:nlat_geo), &
                eldyn % geo_hall_conductivity(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_pedersen_conductivity(1:nlon_geo,1:nlat_geo), &
                eldyn % geo_b_parallel_conductivity(1:nlon_geo,1:nlat_geo) )

      ! When building the Electrodynamics data structure, we set this
      ! module-private switch to true to ensure that the appropriate
      ! initialization is executed for the tiegcm model.
      !tie_gcm_init = .TRUE.

      dynamo_efield = dynamo
      IF ( .not. dynamo ) THEN
        CALL efield_init
        ! Maps latitude from 130km to 90km along flux tube.
        DO j=0,nmlat
          theta130_rad   = ( 180.0_prec - ylatm(j) ) * dtr
          theta90_rad(j) = ASIN(SIN( theta130_rad )*SQRT((earth_radius+90000.0_prec)/(earth_radius+130000.0_prec)))
          IF ( theta130_rad > half_pi ) theta90_rad(j) = pi-theta90_rad(j)
        END DO
      END IF


  END SUBROUTINE Build_IPE_Electrodynamics

  SUBROUTINE Trash_IPE_Electrodynamics( eldyn )
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn

      DEALLOCATE( eldyn % electric_potential, &
                  eldyn % electric_potential2, &
                  eldyn % mhd_electric_potential, &
                  eldyn % electric_field, &
                  eldyn % v_ExB_geographic, &
                  eldyn % v_ExB_apex, &
                  eldyn % hall_conductivity, &
                  eldyn % pedersen_conductivity, &
                  eldyn % b_parallel_conductivity, &
                  eldyn % neutral_apex_velocity, &
                  eldyn % lat_interp_weights, &
                  eldyn % lon_interp_weights, &
                  eldyn % lat_interp_index, &
                  eldyn % lon_interp_index )

      DEALLOCATE( eldyn % geo_electric_potential, &
                  eldyn % geo_v_ExB_geographic, &
                  eldyn % geo_hall_conductivity, &
                  eldyn % geo_pedersen_conductivity, &
                  eldyn % geo_b_parallel_conductivity )

      IF( ALLOCATED( geospace_latitude ) )  DEALLOCATE( geospace_latitude )
      IF( ALLOCATED( geospace_longitude ) ) DEALLOCATE( geospace_longitude )
      IF( ALLOCATED( geospace_potential ) ) DEALLOCATE( geospace_potential )

  END SUBROUTINE Trash_IPE_Electrodynamics


  SUBROUTINE Read_Geospace_Potential( eldyn, filename, error )

    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    CHARACTER(*),                 INTENT(in)    :: filename
    INTEGER,                      INTENT(out)   :: error

    ! Local
    INTEGER :: ncid, ncerr
    INTEGER :: dimid, varid
    INTEGER :: nFluxtube, NLP, NMP
#ifdef HAVE_NETCDF
    CHARACTER(NF90_MAX_NAME) :: nameHolder
#endif

    error = 0

#ifdef HAVE_NETCDF
    _ncCheck( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))

    ! Obtain the dimensions of the Geospace grid
    _ncCheck( nf90_inq_dimid( ncid, "lon", dimid ) )
    _ncCheck( nf90_inquire_dimension( ncid, dimid, nameHolder, n_lon_geospace ) )

    _ncCheck( nf90_inq_dimid( ncid, "lat", dimid ) )
    _ncCheck( nf90_inquire_dimension( ncid, dimid, nameHolder, n_lat_geospace ) )

    IF( .NOT. ALLOCATED( geospace_latitude ) ) ALLOCATE( geospace_latitude(1:n_lat_geospace) )
    IF( .NOT. ALLOCATED( geospace_longitude ) ) ALLOCATE( geospace_longitude(1:n_lat_geospace) )
    IF( .NOT. ALLOCATED( geospace_potential ) ) ALLOCATE( geospace_potential(1:n_lon_geospace,1:n_lat_geospace) )

    _ncCheck( nf90_inq_varid( ncid, "lat", varid ) )
    _ncCheck( nf90_get_var( ncid, varid, geospace_latitude ) )

    _ncCheck( nf90_inq_varid( ncid, "lon", varid ) )
    _ncCheck( nf90_get_var( ncid, varid, geospace_longitude ) )

    _ncCheck( nf90_inq_varid( ncid, "potential", varid ) )
    _ncCheck( nf90_get_var( ncid, varid, geospace_potential ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Read_Geospace_Potential


  SUBROUTINE Interpolate_Geospace_to_MHDpotential( eldyn, grid, time_tracker)

    IMPLICIT NONE

    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)             :: grid
    TYPE( IPE_Time ), INTENT(in)             :: time_tracker

    ! Local
    INTEGER :: j,latidx
    REAL(prec) :: theta110_rad,geospace_latitude_90_rad(1:n_lat_geospace)
    REAL(prec) :: colat_local(1:n_lat_geospace)
    REAL(prec) :: potential_local(1:n_lon_geospace,1:n_lat_geospace)

    DO j = 1, n_lat_geospace
      theta110_rad   = ( 90.0_prec - geospace_latitude(j) ) * dtr
      geospace_latitude_90_rad(j) = ASIN(SIN(theta110_rad)*SQRT((earth_radius+90000.0_prec)/(earth_radius+110000.0_prec)))
      IF ( theta110_rad > half_pi ) geospace_latitude_90_rad(j) = pi - geospace_latitude_90_rad(j)
    ENDDO

    DO j = 1, n_lat_geospace
      latidx = n_lat_geospace-j+1
      colat_local(j)= geospace_latitude_90_rad(latidx)*rtd
      potential_local(:,j)=geospace_potential(:,latidx)
    END DO

!    CALL eldyn % Regrid_Potential( grid, time_tracker, potential_local, geospace_longitude, colat_local, 1, n_lon_geospace, n_lat_geospace )

!    eldyn % mhd_electric_potential= eldyn % electric_potential

  END SUBROUTINE Interpolate_Geospace_to_MHDpotential


  SUBROUTINE Write_MHD_Potential( eldyn, grid, time_tracker, filename, error )

    IMPLICIT NONE

    CLASS( IPE_Electrodynamics ), INTENT(in)  :: eldyn
    TYPE( IPE_Grid ),             INTENT(in)  :: grid
    TYPE( IPE_Time ),             INTENT(in)  :: time_tracker
    CHARACTER(*),                 INTENT(in)  :: filename
    INTEGER,                      INTENT(out) :: error
    ! Local
    REAL(prec) :: time
    INTEGER :: NF90_PREC
    INTEGER :: ncid, ncerr
    INTEGER :: x_dimid, y_dimid, time_dimid, time_varid
    INTEGER :: mhd_phi_varid
    INTEGER :: recStart(1:3), recCount(1:3)


    error = 0

#ifdef HAVE_NETCDF
    recStart = (/ 1, 1, 1 /)
    recCount = (/ grid % NLP, grid % NMP, 1 /)

    time = time_tracker % Calculate_Date_Difference( 2000, 1, 1, 0, 0 )
    IF( prec == sp )THEN
      NF90_PREC = NF90_FLOAT
    ELSE
      NF90_PREC = NF90_DOUBLE
    ENDIF

    _ncCheck( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))

    _ncCheck( nf90_def_dim( ncid, "lp", grid % NLP, x_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "mp", grid % NMP, y_dimid ) )
    _ncCheck( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

    _ncCheck( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "long_name", "minutes since 2000-1-1 00:00 UT" ) )
    _ncCheck( nf90_put_att( ncid, time_varid, "units", "minutes" ) )

    _ncCheck( nf90_def_var( ncid, "mhd_phi", NF90_PREC, (/ x_dimid, y_dimid, time_dimid /) ,mhd_phi_varid ) )
    _ncCheck( nf90_put_att( ncid, mhd_phi_varid, "long_name", "Electric Potential - MHD Component" ) )
    _ncCheck( nf90_put_att( ncid, mhd_phi_varid, "units", "[Unknown]" ) )

    _ncCheck( nf90_enddef(ncid) )

    _ncCheck( nf90_put_var( ncid, time_varid, time ) )
    _ncCheck( nf90_put_var( ncid, mhd_phi_varid, eldyn % mhd_electric_potential ) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Write_MHD_Potential


  SUBROUTINE Read_MHD_Potential( eldyn, filename, error )

    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    CHARACTER(*),                 INTENT(in)    :: filename
    INTEGER,                      INTENT(out)   :: error

    ! Local
    INTEGER :: ncid, ncerr
    INTEGER :: dimid, varid
#ifdef HAVE_NETCDF
    CHARACTER(NF90_MAX_NAME) :: nameHolder
#endif

    error = 0

#ifdef HAVE_NETCDF
    _ncCheck( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))

    _ncCheck( nf90_inq_varid( ncid, "mhd_phi", varid ) )
    _ncCheck( nf90_get_var( ncid, varid, eldyn % mhd_electric_potential) )

    _ncCheck( nf90_close( ncid ) )
#endif

  END SUBROUTINE Read_MHD_Potential


  SUBROUTINE Update_IPE_Electrodynamics( eldyn, grid, forcing, time_tracker, plasma, mpi_layer)
    IMPLICIT NONE
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_Forcing ), INTENT(in)             :: forcing
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    TYPE( IPE_Plasma ), INTENT(in)              :: plasma
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    INTEGER :: lp, mp
    REAL(prec) :: max_v_exb_local
    REAL(prec) :: max_v_exb
#ifdef HAVE_MPI
    INTEGER :: mpiError
#endif

    IF( dynamo_efield ) THEN

      CALL eldyn % Dynamo_Wrapper(grid, forcing, time_tracker, plasma, mpi_layer )
      IF( mpi_layer % rank_id == 0 )THEN
        print *,'Dynamo E field ', int(time_tracker % elapsed_sec / 60), ' mins UT'
      ENDIF

    ELSE

      CALL eldyn % Empirical_E_Field_Wrapper( grid, forcing, time_tracker, mpi_layer)
      IF( mpi_layer % rank_id == 0 )THEN
        print *,'TZU-WEI calling empirical E field'
      ENDIF


!      IF( geospace )THEN

!        CALL  eldyn % Read_Geospace_Potential( filename )
!        CALL  eldyn % Regrid_Geospace(  )
!#ifdef DEBUG
!        CALL eldyn % Write_MHD_Potential( )
!#endif
!        CALL eldyn % Merge_Geospace_Potential

!      ELSEIF( openggcm )THEN

!        CALL  eldyn % Read_OpenGGCM_Potential( filename )
!        CALL  eldyn % Regrid_OpenGGCM(  )
!#ifdef DEBUG
!        CALL eldyn % Write_MHD_Potential( )
!#endif
!        CALL eldyn % Merge_OpenGGCM_Potential

!      ENDIF


      ! Calculate the potential gradient in IPE coordinates.
      CALL eldyn % Calculate_Potential_Gradient( grid )

    ENDIF

    ! Calculate ExB drift velocity
    CALL eldyn % Calculate_ExB_Velocity( grid , mpi_layer, max_v_exb_local )

!   IF( mpi_layer % rank_id == 0 ) print *,'End Update Electrodynamics'

  END SUBROUTINE Update_IPE_Electrodynamics


  SUBROUTINE Empirical_E_Field_Wrapper( eldyn, grid, forcing, time_tracker, mpi_layer)

    IMPLICIT NONE

    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ),             INTENT(in)    :: grid
    TYPE( IPE_Forcing ),          INTENT(in)    :: forcing
    TYPE( IPE_Time ),             INTENT(in)    :: time_tracker
    TYPE( IPE_MPI_Layer ),        INTENT(in)    :: mpi_layer

    ! Local
    INTEGER :: i, j, year
    INTEGER :: lp, mp
    REAL(prec) :: utime, lat
    REAL(prec) :: potent_local(0:nmlon,0:nmlat)
    REAL(prec) :: mlt_local(0:nmlon)
    REAL(prec) :: colat_local(0:nmlat)

    ! iday, iday_m, ut, f107d, bt, angle, v_sw, bz are all variables
    ! declared in efield_ipe.f
    CALL time_tracker % Get_Date( iyear, imo, iday_m )
    iday = time_tracker % day_of_year

    CALL time_tracker % Get_UTime( utime )
    ut=utime/3600.0_prec

    f107d = forcing % f107( forcing % current_index )
    v_sw  = forcing % solarwind_velocity( forcing % current_index )
    bz    = forcing % solarwind_Bz( forcing % current_index )
    by    = forcing % solarwind_By( forcing % current_index )

    iyear = 1999

    call get_efield(mpi_layer % rank_id)

    ! Interpolate the potential to the IPE grid
    potent_local = potent
    mlt_local    = ylonm
    colat_local  = rtd * theta90_rad
    CALL eldyn % Regrid_Potential( grid, mpi_layer, time_tracker, potent_local, mlt_local, colat_local, 0, nmlon, nmlat )

  END SUBROUTINE Empirical_E_Field_Wrapper


  SUBROUTINE Calculate_ExB_Velocity( eldyn, grid, mpi_layer, max_v_exb )
    ! Calculates ExB drift velocity according to Eqs(4.18-4.19) in
    ! A.D. Richmond (1995)
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    INTEGER    :: mp, lp, i_90km
    REAL(prec) :: max_v_exb
    REAL(prec) :: v_boost_factor

    v_boost_factor = 1.0

!   write(1000 + mpi_layer % rank_id, *) 'TIME'

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
         DO i_90km = 1, grid % flux_tube_max(lp)

        eldyn % v_ExB_apex(1,lp,mp) = eldyn % electric_field(2,lp,mp)/grid % apex_be3(lp,mp)
        eldyn % v_ExB_apex(2,lp,mp) = -eldyn % electric_field(1,lp,mp)/grid % apex_be3(lp,mp)
!         eldyn % v_ExB_apex(1,lp,mp) = v_boost_factor * (eldyn % electric_field(2,lp,mp)/grid % apex_be3(lp,mp))
!         eldyn % v_ExB_apex(2,lp,mp) = v_boost_factor * (-eldyn % electric_field(1,lp,mp)/grid % apex_be3(lp,mp))

! TWFANG modify the equation and add i_90km (May 2019)
! geographic eastward....
        eldyn % v_exb_geographic(1,i_90km,lp,mp) = (eldyn % v_ExB_apex(1,lp,mp) * grid % apex_e_vectors(1,1,i_90km,lp,mp)) &
                                                 + (eldyn % v_ExB_apex(2,lp,mp) * grid % apex_e_vectors(1,2,i_90km,lp,mp))
! geographic northward....
        eldyn % v_exb_geographic(2,i_90km,lp,mp) = (eldyn % v_ExB_apex(1,lp,mp) * grid % apex_e_vectors(2,1,i_90km,lp,mp)) &
                                                 + (eldyn % v_ExB_apex(2,lp,mp) * grid % apex_e_vectors(2,2,i_90km,lp,mp))
! geographic upwards....
        eldyn % v_exb_geographic(3,i_90km,lp,mp) = (eldyn % v_ExB_apex(1,lp,mp) * grid % apex_e_vectors(3,1,i_90km,lp,mp)) &
                                                 + (eldyn % v_ExB_apex(2,lp,mp) * grid % apex_e_vectors(3,2,i_90km,lp,mp))
        ENDDO

!       write(1000 + mpi_layer % rank_id, 255) mp,lp,eldyn % v_ExB_apex(1,lp,mp),eldyn % v_ExB_apex(2,lp,mp)
!255    format(2i4,2e12.4)

      ENDDO

    ENDDO
    max_v_exb = maxval(eldyn % v_ExB_apex)

  END SUBROUTINE Calculate_ExB_Velocity


  SUBROUTINE Calculate_Potential_Gradient( eldyn, grid )
    ! Uses 2nd order centered differencing (on the apex grid) to calculate the
    ! electric field components from the potential attribute from the
    ! IPE_Electrodynamics class.
    ! The electric field gradient is calculated using Sections 3 and 4 of A.D.
    ! Richmond (1995)
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ),             INTENT(in)    :: grid

    ! Local
    INTEGER    :: mp, lp
    REAL(prec) :: r, d_lp, d_mp, coslam, sinim


    r = earth_radius + 90000.0_prec

    ! longitude component of e-field ( ed1 ): Eq(4.8)
    ! Note that ed1 is calculated slightly differently here by dividing
    ! through by a factor of sin_Im for accomodating the transport.
    ! The distance conversion from radians to meters involves (here) i
    ! multiplication of the difference in radians (between two grid points)
    ! by sin_Im*r.

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP

        coslam = cos( half_pi - grid % magnetic_colatitude(1,lp) )
        d_mp   = coslam*r*( grid % magnetic_longitude(mp+1) - grid % magnetic_longitude(mp-1) )

        eldyn % electric_field(1,lp,mp) = -( eldyn % electric_potential(lp,mp+1) - &
                                            eldyn % electric_potential(lp,mp-1) )/d_mp

!       if (mp.eq.1.and.lp.eq.40) then
!         write(612,*) 'GHGM POT ', eldyn % electric_potential(lp,mp+1), eldyn % electric_potential(lp,mp-1), d_mp
!       endif
      ENDDO
    ENDDO


    ! latitude component of e-field ( ed2 ): Eq(4.9)
    ! Note that ed1 is calculated slightly differently here by dividing
    ! through by a factor of sin_Im for accomodating the transport.
    ! The distance conversion from radians to meters involves (here) i
    ! multiplication of the difference in radians (between two grid points)
    ! by sin_Im*r.
    DO mp = grid % mp_low, grid % mp_high

      lp = 1
      coslam = cos( half_pi - grid % magnetic_colatitude(1,lp) )
      sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
      d_lp = sinim*r*( grid % magnetic_colatitude(1,lp+1) - grid % magnetic_colatitude(1,lp) )

      eldyn % electric_field(2,lp,mp) = -( eldyn % electric_potential(lp+1,mp) - &
                                          eldyn % electric_potential(lp,mp) )/d_lp

      DO lp = 2, grid % NLP-1

        coslam = cos( half_pi - grid % magnetic_colatitude(1,lp) )
        sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
        d_lp = sinim*r*( grid % magnetic_colatitude(1,lp+1) - grid % magnetic_colatitude(1,lp-1) )

        eldyn % electric_field(2,lp,mp) = -( eldyn % electric_potential(lp+1,mp) - &
                                            eldyn % electric_potential(lp-1,mp) )/d_lp

      ENDDO

      lp = grid % NLP
      coslam = cos( half_pi - grid % magnetic_colatitude(1,lp) )
      sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
      d_lp = sinim*r*( grid % magnetic_colatitude(1,lp) - grid % magnetic_colatitude(1,lp-1) )

      eldyn % electric_field(2,lp,mp) = -( eldyn % electric_potential(lp,mp) - &
                                          eldyn % electric_potential(lp-1,mp) )/d_lp

    ENDDO

  END SUBROUTINE Calculate_Potential_Gradient


  SUBROUTINE Regrid_Potential( eldyn, grid, mpi_layer, time_tracker, potential, mlt, colat, start_index, nlon, nlat )
  ! This subroutine regrids electric potential from a structured grid
  ! on (magnetic local time, magnetic colatitude) to IPE's grid.
  !
  ! Input :
  !
  !   grid
  !
  !   time_tracker
  !
  !   potential
  !
  !   mlt - magnetic local time [ deg ]
  !
  !   colat - magnetic colatitude [ deg ]
  !
  !   start_index - starting index for the potential, mlt, and colat arrays
  !
  !   nlon - last index in the longitude direction
  !
  !   nlat - last index in the latitude direction
  !
    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)                :: grid
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    TYPE( IPE_Time ), INTENT(in)                :: time_tracker
    INTEGER, INTENT(in)                         :: start_index, nlon, nlat
    REAL(prec), INTENT(in)                      :: potential(start_index:nlon,start_index:nlat)
    REAL(prec), INTENT(in)                      :: mlt(start_index:nlon)
    REAL(prec), INTENT(in)                      :: colat(start_index:nlat)

    ! Local
    INTEGER :: mp, lp, i, j, i1, i2, j1, j2, ii
    INTEGER :: ilat(1:grid % NLP), jlon(grid % mp_low - grid % mp_halo:grid % mp_high+grid % mp_halo)
    REAL(prec) :: lat, lon, dlon, difflon, mindiff
    REAL(prec) :: lat_weight(1:2), lon_weight(1:2)
    REAL(prec) :: mlon90_rad(start_index:nlon)
    REAL(prec) :: mlat90_rad(start_index:nlat)

    IF( dynamo_efield ) THEN
      mlon90_rad = dtr * mlt
      mlat90_rad = dtr * colat
    ELSE
      mlon90_rad = MLT_to_MagneticLongitude( mlt, 1999, time_tracker % day_of_year, time_tracker % utime, start_index, nlon )
      mlat90_rad = dtr * colat
    ENDIF


    ! Search for nearest grid points in the magnetic longitude/latitude grid
    DO mp = grid % mp_low - grid % mp_halo, grid % mp_high + grid % mp_halo
      lon = grid % magnetic_longitude(mp)


      ! Note that the mlon90_rad(j) array is not monotonically increasing
      ! with index j. This happens because the empirical model provides data
      ! on the MLT grid so that magnetic longitude=0 may not occur at j=0.
      ! Instead, magnetic longitude =0 may occur at any index in mlon90_rad.
      ! Because of this, we have to perform a full search of the array

      ! Initial condition for the minimum difference is something absurd that
      ! is guaranteed to be larger than the differences calculated.
      mindiff = 1000.0_prec
      DO j = start_index, nlon
        difflon = ABS( mlon90_rad(j) - lon )

        IF( difflon < mindiff )THEN
          jlon(mp) = j
          mindiff = difflon
        ENDIF

      ENDDO

    ENDDO

    DO lp = 1, grid % NLP

      lat = grid % magnetic_colatitude(1,lp)
      ! colatitude decreases with increasing lp
      ilat(lp) = nlat
      DO i = start_index, nlat

        ! Need to pass in colatitude through the call stack (theta90_rad ->
        ! colatitude )
        IF( mlat90_rad(i) < lat )THEN
          ilat(lp) = i

          EXIT
        ENDIF

      ENDDO

    ENDDO

    DO mp = grid % mp_low - grid % mp_halo, grid % mp_high + grid % mp_halo
      DO lp = 1, grid % NLP

        lat = grid % magnetic_colatitude(1,lp)
        lon = grid % magnetic_longitude(mp)

        IF( ilat(lp) == start_index )THEN

          i1 = ilat(lp)
          i2 = ilat(lp)
          lat_weight(1) = 1.0_prec
          lat_weight(2) = 0.0_prec

        ELSE

          i1 = ilat(lp)-1
          i2 = ilat(lp)
          lat_weight(1) =  ( lat - mlat90_rad(i2) )/( mlat90_rad(i1) - mlat90_rad(i2) )
          lat_weight(2) = -( lat - mlat90_rad(i1) )/( mlat90_rad(i1) - mlat90_rad(i2) )

        ENDIF

        IF( jlon(mp) == start_index )THEN

          IF( lon > mlon90_rad(jlon(mp)) ) THEN

            ! In this case, the nearest empricial grid point is the first index
            ! and is to the left of the IPE grid point. The indices bounding
            ! the IPE grid point are jlon(mp) and jlon(mp)+1

            j1 = start_index
            j2 = start_index+1

          ELSE

            ! In this case the first longitude point on the empirical model grid
            ! is greater than the mp-longitude point on the IPE grid. Here, the
            ! bounding points are at empirical longitude points start_index and nlon

            j1 = nlon
            j2 = start_index

          ENDIF


        ELSEIF( jlon(mp) == nlon )THEN

          IF( lon > mlon90_rad(jlon(mp)) ) THEN

            ! In this case the first longitude point on the empirical model grid
            ! is greater than the mp-longitude point on the IPE grid. Here, the
            ! bounding points are at empirical longitude points start_index and nlon

            j1 = nlon
!           j2 = start_index
! GHGM - this special case:
            j2 = start_index + 1
! GHGM

          ELSE

            j1 = nlon-1
            j2 = nlon

          ENDIF

        ELSE

          ! In this case, the IPE longitude point is interior to the empirical
          ! model grid. The point we found (associated with jlon(mp)) is
          ! greater than the IPE longitude point. Because of this, the bounding
          ! points are at [jlon(mp)-1] and at [jlon(mp)].

          IF( lon > mlon90_rad(jlon(mp)) ) THEN

            ! In this case the first longitude point on the empirical model grid
            ! is greater than the mp-longitude point on the IPE grid. Here, the
            ! bounding points are at empirical longitude points start_index and nlon

            j1 = jlon(mp)
            j2 = jlon(mp)+1

          ELSE

            j1 = jlon(mp)-1
            j2 = jlon(mp)

          ENDIF

        ENDIF

        dlon = mlon90_rad(j1) - mlon90_rad(j2)
!       IF( dlon > 0.0_prec )THEN
! GHGM .ge. (I think)
        IF( dlon.ge.0.0_prec )THEN

          dlon = dlon - two_pi

          IF( mp <= 1 )THEN

            lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
            lon_weight(2) = -( lon - (mlon90_rad(j1)-two_pi) )/dlon

          ELSEIF( mp >= grid % NMP )THEN

            lon_weight(1) =  ( lon - (mlon90_rad(j2)+two_pi) )/dlon
            lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

          ENDIF

        ELSE

          lon_weight(1) =  ( lon - mlon90_rad(j2) )/dlon
          lon_weight(2) = -( lon - mlon90_rad(j1) )/dlon

        ENDIF

        eldyn % electric_potential(lp,mp) = ( potential(j1,i1)*lat_weight(1)*lon_weight(1) +&
                                              potential(j2,i1)*lat_weight(1)*lon_weight(2) +&
                                              potential(j1,i2)*lat_weight(2)*lon_weight(1) +&
                                              potential(j2,i2)*lat_weight(2)*lon_weight(2) )

      ENDDO
    ENDDO

  END SUBROUTINE Regrid_Potential


  FUNCTION MLT_to_MagneticLongitude( mlt, year, day_of_year, utime, start_index, nlon ) RESULT( mag_longitude )

    INTEGER    :: start_index, nlon
    REAL(prec) :: mlt(start_index:nlon)
    INTEGER    :: year, day_of_year
    REAL(prec) :: utime
    REAL(prec) :: mag_longitude(start_index:nlon)

    ! Local
    INTEGER    :: i
    REAL       :: sunlons

    ! Map magnetic local time to magnetic longitude
    CALL sunloc( year, day_of_year, utime, sunlons )
    DO i=start_index,nlon

      mag_longitude(i) = dtr * (mlt(i)-180.0_prec) + sunlons
      IF( mag_longitude(i) < 0.0_prec ) mag_longitude(i) = mag_longitude(i) + two_pi
      IF( mag_longitude(i) >= two_pi  ) mag_longitude(i) = mag_longitude(i) - two_pi

    END DO

  END FUNCTION MLT_to_MagneticLongitude


  SUBROUTINE sunloc(iyr,iday,secs,sunlons)

    integer,    intent(in)  :: iyr, iday ! day of year
    REAL(prec), intent(in)  ::  secs    ! ut in seconds
    REAL, INTENT(out) :: sunlons

    integer :: ihr,imn
    real :: sec,date,vp,xmlon ! apex magnetic longitude
    real ::  sbsllat    ! geographic latitude of subsolar point (degrees)
    real ::  sbsllon    ! geographic longitude of subsolar point (degrees)
    real ::  colat      ! Geocentric colatitude of geomagnetic dipole north pole (deg)
    real ::  elon        ! East longitude of geomagnetic dipole north pole (deg)

    ihr = int(secs/3600.)
    imn = int((secs - float(ihr)*3600.)/60.)
    sec = secs - float(ihr)*3600. - float(imn)*60.

    !  calculate subsol point: given universal time
    !          input: iyr,iday,ihr,imn,sec
    !          output: sbsllat,sbsllon
    !
    call subsol_empirical(iyr,iday,ihr,imn,sec ,sbsllat,sbsllon)

    date = float(iyr) + float(iday)/365. + float(ihr)/24./365. + &
           float(imn)/60./24./365.+ sec/60./60./24./365.

    call cofrm_empirical(date)
    call dypol_empirical(colat,elon,vp)

    ! calculate geomagn. diploe longitude
    !        input: aloni,sbsllat,sbsllon,colat,elon
    !        output: xmlon
    call solgmlon_empirical(sbsllat,sbsllon,colat,elon,xmlon)
    sunlons = xmlon*dtr

  END SUBROUTINE sunloc


  SUBROUTINE Interpolate_to_GeographicGrid_IPE_Electrodynamics( eldyn, grid )

    IMPLICIT NONE

    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Grid ), INTENT(in)         :: grid

    CALL grid % Interpolate_2D_to_Geographic_Grid( eldyn % electric_potential, eldyn % geo_electric_potential )
    CALL grid % Interpolate_2D_to_Geographic_Grid( eldyn % mhd_electric_potential, eldyn % geo_mhd_electric_potential )

  END SUBROUTINE Interpolate_to_GeographicGrid_IPE_Electrodynamics


  SUBROUTINE Dynamo_Wrapper(eldyn,grid,forcing,time_tracker, plasma, mpi_layer)

    use module_init_cons!,only:init_cons
    use module_init_heelis!,only:init_heelis
    use module_readin_ascii!,only:readin_ascii
    use module_highlat!,only: highlat
    use module_sub_dynamo!,only: dynamo
    use params_module
    use cons_module
    use dynamo_module
    use heelis_module, only:ctpoten
    use module_magfield
!!     use module_update_fli,ONLY:update_fli

    CLASS( IPE_Electrodynamics ), INTENT(inout) :: eldyn
    TYPE( IPE_Forcing ), INTENT(in)             :: forcing
    TYPE( IPE_Time ),             INTENT(in)    :: time_tracker
    TYPE( IPE_Grid ),             INTENT(in)    :: grid
    TYPE( IPE_Plasma ),           INTENT(in)    :: plasma
    TYPE( IPE_MPI_Layer ),        INTENT(in)    :: mpi_layer

    INTEGER, PARAMETER :: dyn_midpoint=48 
    REAL(prec) :: ed_IPE(2,grid % NLP,grid % NMP,2) !1:N/SH....4:ed1/2
    REAL(prec) :: eldyn_conductivities(6,2,dyn_midpoint,grid % NMP) !1:N/SH....4:ed1/2
    REAL(prec) :: ed_conductivities(6,kmlon+1,kmlat)
    INTEGER    :: year,i,i_dyn,j,k,l,ilat_dyn,ilon_dyn
    INTEGER    :: diffval,imlat_plas,imlat_dyn
    INTEGER    :: tube_need(dyn_midpoint),ihem,lp_dyn
    REAL(prec) :: mlat_plas,mlat_dyn
    REAL(prec) :: fkp
    REAL(prec) :: ed_dyn(170,80,2)
    REAL(prec) :: theta0(kmlat)
    REAL(prec) :: xlonm_deg_map(82)
    REAL(prec) :: ylatm_deg_map(kmlat)
    REAL(prec) :: ed1dy_map(82,kmlat)
    REAL(prec) :: ed2dy_map(82,kmlat)

    eldyn_conductivities(:,:,:,:)=0.
    ed_conductivities(:,:,:)=0.
    CALL init_cons

    sangle= forcing % solarwind_angle ( forcing % current_index )
    bt= forcing % solarwind_Bt (forcing % current_index )
    swvel= forcing % solarwind_velocity (forcing % current_index )
    swden= forcing % solarwind_density (forcing % current_index )
    stilt= get_tilt(time_tracker%year,time_tracker%month,time_tracker%day,time_tracker%utime)

    fkp= forcing % kp ( forcing % current_index )
    ctpoten= 15.+15.*fkp+0.8*fkp**2
    CALL init_heelis

    year=2000

    CALL sunloc( year, time_tracker % day_of_year, time_tracker % utime, sunlons )
!!---
!!! read in integrals
!        call readin_ascii

    DO i = 1, grid % NLP
      mlat_plas  = 90. - grid % magnetic_colatitude(1,i)*rtd
      imlat_plas = INT(mlat_plas*10.)
      DO i_dyn=1,dyn_midpoint !from SH toward eq
        mlat_dyn  = xlatm_deg(i_dyn) ![deg]
        imlat_dyn = INT(mlat_dyn*10.)

        diffval = imlat_plas -abs(imlat_dyn)

        if (abs(diffval) < 1)  then
          tube_need(i_dyn)=i
        endif
      ENDDO
    ENDDO
    DO i=1, grid % NMP
      DO j= 1, 48
        DO k=1,2
         if (tube_need(j).lt.1) then
         eldyn_conductivities(1:6,k,j,i)=0.
         else
         eldyn_conductivities(1:6,k,j,i)=plasma % conductivities(1:6,k,tube_need(j),i)
         endif
! no FLI values on 1,47,48,49(mag. equ.)
! eldyn_conductivities(1:6,1:2,2:47,1:grid % NMP)=plasma % conductivities(1:6,1:2,tube_need,1:grid % NMP)
! print *,'eldyn',eldyn_conductivities(l,k,j,i),'plasma',plasma % conductivities(l,k,tube_need(j),i)
         ENDDO
       ENDDO
    ENDDO
    
    DO j=1,6
      DO lp_dyn=2,46 !from SH toward eq
        DO ihem=1,2
          if (ihem==1) then !NH
            ilat_dyn = kmlat - lp_dyn +1
          else if ( ihem==2 ) then !SH
            ilat_dyn = lp_dyn
          end if
          DO i=1,kmlon
            ilon_dyn = i + (kmlon / 2)
            if ( ilon_dyn > kmlon+1 ) ilon_dyn = ilon_dyn - kmlon

            ed_conductivities(j,ilon_dyn,ilat_dyn) = eldyn_conductivities(j,ihem,lp_dyn,i)

          ENDDO
        ENDDO
      ENDDO

!perdiodic solution: kmlon+1<--mp=1
!dbg20140804         dum(kmlon+1,:)=dum(1,:)
      ed_conductivities(j,1,:)=ed_conductivities(j,kmlon+1,:)!dbg20140804

! magnetic poles: is this correct? all longitudes should be collapsed to one
! same value?
      do i=1,kmlon+1
        ed_conductivities(j,i,1)=ed_conductivities(j,i,2) !South
        ed_conductivities(j,i,kmlat)=ed_conductivities(j,i,kmlat-1) !North
        ed_conductivities(j,i,49)=(ed_conductivities(j,i,46)+ed_conductivities(j,i,52))*0.5 !magnetic eq
        ed_conductivities(j,i,48)=(ed_conductivities(j,i,49)+ed_conductivities(j,i,46))*0.5
        ed_conductivities(j,i,47)=(ed_conductivities(j,i,48)+ed_conductivities(j,i,46))*0.5
        ed_conductivities(j,i,50)=(ed_conductivities(j,i,49)+ed_conductivities(j,i,52))*0.5
        ed_conductivities(j,i,51)=(ed_conductivities(j,i,50)+ed_conductivities(j,i,52))*0.5
      enddo

    ENDDO

    zigm11=ed_conductivities(1,:,:)
    zigm22=ed_conductivities(2,:,:)
    zigm2 =ed_conductivities(3,:,:)
    zigmc =ed_conductivities(4,:,:)
    rim(:,:,1)=ed_conductivities(5,:,:)
    rim(:,:,2)=ed_conductivities(6,:,:)

!UNDERCONSTRUCTION!!!
!dbg20150608: copied from module_readin_ascii.f
! am 10/04 so far no value at the equator therefore set it
! but this shouldn't be in the code
    j = kmlat/2+1
    do i = 1,kmlonp1
      zigm11(i,j)= .125*(zigm11(i,j-1)+ zigm11(i,j+1))
      zigm22(i,j)= .125*(zigm22(i,j-1)+ zigm22(i,j+1))
      zigmc(i,j) = .125*(zigmc(i,j-1) + zigmc(i,j+1))
      zigm2(i,j) = .06 *(zigm2(i,j-1) + zigm2(i,j+1))
      rim(i,j,1) = .06 *(rim(i,j-1,1) + rim(i,j+1,1))
      rim(i,j,2) = .06 *(rim(i,j-1,2) + rim(i,j+1,2))
    enddo ! i = 1,kmlon
!
! am 10/04 change sign of K_(m lam)^D in the SH- that's what TIEGCM dynamo
! expects
    do j = 1,(kmlat+1)/2
!     rim(:,lp,2) = rim(:,lp,2)*1.e8
      rim(:,j,2) = -rim(:,j,2)
    enddo

! // TODO // !
! HEY YOU ! PAY ATTENTION *!
! Can we push the dynamo grid to match the IPE grid (NMP,2*NLP)?
! ************************ !

    call highlat

    call dynamo

    dlonm = two_pi/float(kmlon)

    do i=1,kmlonp1
      xlonm(i) = -pi+float(i-1)*dlonm
      xlonm_deg(i) = xlonm(i)*rtd
    enddo ! i=1,kmlonp1

    xlonm_deg_map(2:41)=xlonm_deg(41:80)
    xlonm_deg_map(42:81)=xlonm_deg(1:40)+360.
    xlonm_deg_map(1)=xlonm_deg(40)
    xlonm_deg_map(82)=xlonm_deg(41)+360.
    ylatm_deg_map=90.-xlatm_deg

    ed1dy_map(2:41,:)=ed1dy(41:80,:)
    ed1dy_map(42:81,:)=ed1dy(1:40,:)
    ed1dy_map(1,:)=ed1dy(40,:)
    ed1dy_map(82,:)=ed1dy(41,:)

    ed2dy_map(2:41,:)=ed2dy(41:80,:)
    ed2dy_map(42:81,:)=ed2dy(1:40,:)
    ed2dy_map(1,:)=ed2dy(40,:)
    ed2dy_map(82,:)=ed2dy(41,:)

    CALL eldyn % Regrid_Potential( grid,mpi_layer, time_tracker,ed1dy_map,xlonm_deg_map,ylatm_deg_map, 1, 82,kmlat )
      eldyn % electric_field(1,:,:) = eldyn % electric_potential

    CALL eldyn % Regrid_Potential( grid,mpi_layer, time_tracker,ed2dy_map,xlonm_deg_map,ylatm_deg_map, 1, 82,kmlat )
      eldyn % electric_field(2,:,:) = eldyn % electric_potential

  END SUBROUTINE Dynamo_Wrapper

        FUNCTION GET_TILT(YEAR,MONTH,DAY,HOUR)  RESULT(get_tilt_angle)
!
!-----------------------------------------------------------------------
!  It has been changed to return the dipole tilt from this function
!  call.      
!         
!      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
!      TRANSFORMATIONS, IDENTIFIED BY K.
!          K=1 TRANSFORMS GSE to GEO
!          K=2     "      GEO to MAG
!          K=3     "      GSE to MAG
!          K=4     "      GSE to GSM
!          K=5     "      GEO to GSM
!          K=6     "      GSM to MAG
!          K=7     "      GSE to GEI
!          K=8     "      GEI to GEO
!          K=9     "      GSM to SM 
!          K=10    "      GEO to SM 
!          K=11    "      MAG to SM 
!
!
!      The formal names of the coordinate systems are:
!       GSE - Geocentric Solar Ecliptic
!       GEO - Geographic
!       MAG - Geomagnetic
!       GSM - Geocentric Solar Magnetospheric
!       SM  - Solar Magnetic
!       
!      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
!      ST(I) AND CT(I) ARE SINES & COSINES.       
!-------------------------------------------------------------------------

       implicit none
!
!-----------------------------Return Value--------------------------
!
        real(prec)  ::  get_tilt_angle,HOUR
!  ------------------------------Arguments--------------------------------
!
!       INTEGER YEAR, MONTH, DAY,JULDAY_WAM
        INTEGER YEAR, MONTH, DAY
!
!-----------------------------Arrays-----------------------------

        real(prec) CX(9),ST(6),CT(6),AM(3,3,11)

!-----------------------------Parameters------------------------------
!
      real(prec) , parameter :: EPOCH=1980.,TH0=11.19,PH0=-70.76,   &
                          DIPOLE=.30574

        INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
        INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
        INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
        INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
        INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

        PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
        PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
        PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
        PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
        PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
        PARAMETER (MAGSM =11,SMMAG =-11)
!
!---------------------------Local variables-----------------------------
!
        integer j, k, jd, iyr, i, mjd

        REAL(prec) UT, T0, GMSTD, GMSTH, ECLIP, MA, LAMD, SUNLON, pi
        real(prec) b32, b33, b3
!
!-----------------------------------------------------------------------


        pi=3.141592653


        IF(YEAR.LT.1900)THEN
          IYR=1900+YEAR
        ELSE
          IYR=YEAR
        ENDIF
        UT=HOUR
        JD=JULDAY_WAM(MONTH,DAY,IYR)
        MJD=JD-2400001
!       T0=(real(MJD,r8)-51544.5)/36525.0
        T0=(float(MJD)-51544.5)/36525.0
        GMSTD=100.4606184 +36000.770*T0 +3.87933E-4*T0*T0 +   &
              15.0410686*UT
        CALL ADJUST(GMSTD)
        GMSTH=GMSTD*24./360.
        ECLIP=23.439 - 0.013*T0
        MA=357.528 + 35999.050*T0 + 0.041066678*UT
        CALL ADJUST(MA)
        LAMD=280.460 + 36000.772*T0 + 0.041068642*UT
        CALL ADJUST(LAMD)
        SUNLON=LAMD + (1.915-0.0048*T0)*SIN(MA*pi/180.) + 0.020*  &
               SIN(2.*MA*pi/180.)
        CALL ADJUST(SUNLON)

               CX(1)= GMSTD
        CX(2) = ECLIP
        CX(3) = SUNLON
        CX(4) = TH0
        CX(5) = PH0
! Derived later:
!       CX(6) = Dipole tilt angle  
!       CX(7) = Angle between sun and magnetic pole
!       CX(8) = Subsolar point latitude
!       CX(9) = Subsolar point longitude

        DO I=1,5
          ST(I) = SIN(CX(I)*pi/180.)
          CT(I) = COS(CX(I)*pi/180.)
        ENDDO
!         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0.
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)
!         
      AM(1,1,GEIGEO) = CT(1)
      AM(1,2,GEIGEO) = ST(1)
      AM(1,3,GEIGEO) = 0.
      AM(2,1,GEIGEO) = -ST(1)
      AM(2,2,GEIGEO) = CT(1)
      AM(2,3,GEIGEO) = 0.
      AM(3,1,GEIGEO) = 0.
      AM(3,2,GEIGEO) = 0.
      AM(3,3,GEIGEO) = 1.
!         
      DO I=1,3
      DO J=1,3
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) +   &
      AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
!         
      AM(1,1,GEOMAG) = CT(4)*CT(5)
      AM(1,2,GEOMAG) = CT(4)*ST(5)
      AM(1,3,GEOMAG) =-ST(4)
      AM(2,1,GEOMAG) =-ST(5)
      AM(2,2,GEOMAG) = CT(5)
      AM(2,3,GEOMAG) = 0.
      AM(3,1,GEOMAG) = ST(4)*CT(5)
      AM(3,2,GEOMAG) = ST(4)*ST(5)
      AM(3,3,GEOMAG) = CT(4)
!         
      DO I=1,3
      DO J=1,3
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) +   &
      AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
!         
      B32 = AM(3,2,GSEMAG)
      B33 = AM(3,3,GSEMAG)
      B3  = SQRT(B32*B32+B33*B33)
      IF (B33.LE.0.) B3 = -B3
!         
      AM(2,2,GSEGSM) = B33/B3
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)
      AM(3,2,GSEGSM) = B32/B3
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)
      AM(1,1,GSEGSM) = 1.
      AM(1,2,GSEGSM) = 0.
      AM(1,3,GSEGSM) = 0.
      AM(2,1,GSEGSM) = 0.
      AM(3,1,GSEGSM) = 0.
!         
      DO I=1,3
      DO J=1,3
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) +    &
      AM(I,2,GSEGSM)*AM(J,2,GSEGEO) +                       &
      AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
!         
      DO I=1,3
      DO J=1,3
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) +    &
      AM(I,2,GEOMAG)*AM(J,2,GEOGSM) +                       &
      AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO
!
        ST(6) = AM(3,1,GSEMAG)
!        CT(6) = SQRT(1.-ST(6)*ST(6))
        CX(6) = ASIN(ST(6)*pi/180.)

             GET_TILT_ANGLE = ST(6)

      END FUNCTION GET_TILT

!================================================================================================

!       INTEGER FUNCTION JULDAY_WAM(MM,ID,IYYY)

       FUNCTION JULDAY_WAM(MM,ID,IYYY)  RESULT(JULDAY)
!
!-----------------------------------------------------------------------
!
!     use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none
!
      integer :: JULDAY

!------------------------------Arguments--------------------------------
!
      integer mm, id, iyyy
!
!-----------------------------Parameters------------------------------
!
      integer igreg
      PARAMETER (IGREG=15+31*(10+12*1582))
!
!---------------------------Local variables-----------------------------
!
      integer ja, jm, jy
!
!-----------------------------------------------------------------------
!
!!!compiler warning      IF (IYYY.EQ.0) PAUSE 'There is no Year Zero.'
      IF (IYYY.EQ.0) STOP 'There is no Year Zero.'
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      END FUNCTION JULDAY_WAM

!================================================================================================

        SUBROUTINE ADJUST(ANGLE)
!
!-----------------------------------------------------------------------
!       ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
!-----------------------------------------------------------------------
!
!
!       use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none
!
!------------------------------Arguments--------------------------------
!
        real(prec) angle
!
!-----------------------------------------------------------------------
!
! 10     CONTINUE
        IF(ANGLE.LT.0.)THEN
          ANGLE=ANGLE+360.
!          GOTO 10
        ENDIF
! 20     CONTINUE
        IF(ANGLE.GE.360.)THEN
          ANGLE=ANGLE-360.
!          GOTO 20
        ENDIF
!        RETURN
        END SUBROUTINE ADJUST
                              


!-----------------------------------------------------------------------
!
END MODULE IPE_Electrodynamics_Class
