MODULE IPE_Plasma_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Grid_Class
  USE IPE_Neutrals_Class
  USE IPE_Forcing_Class
  USE IPE_Time_Class
  USE IPE_MPI_Layer_Class
  USE IPE_Common_Routines

  USE HDF5

  IMPLICIT NONE

  TYPE IPE_Plasma
    INTEGER :: nFluxTube, NLP, NMP
    INTEGER :: mp_low, mp_high, mp_halo

    REAL(prec), ALLOCATABLE :: ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_density2(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: ion_densities_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperature_old(:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density_old(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature_old(:,:,:)

    REAL(prec), ALLOCATABLE :: ionization_rates(:,:,:,:)
    REAL(prec), ALLOCATABLE :: conductivities(:,:,:,:)

    ! Interpolated Fields
    REAL(prec), ALLOCATABLE :: geo_ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ion_velocities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ion_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_velocity(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_electron_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_ionization_rates(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_tec(:,:)
    REAL(prec), ALLOCATABLE :: geo_nmf2(:,:)
    REAL(prec), ALLOCATABLE :: geo_hmf2(:,:)


    CONTAINS

      PROCEDURE :: Build => Build_IPE_Plasma
      PROCEDURE :: Trash => Trash_IPE_Plasma

      PROCEDURE :: Update => Update_IPE_Plasma
      PROCEDURE :: Update_Halos => Update_Halos_IPE_Plasma
      PROCEDURE :: Interpolate_to_GeographicGrid => Interpolate_to_GeographicGrid_IPE_Plasma
      PROCEDURE :: Read_Legacy_Input => Read_Legacy_Input_IPE_Plasma

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: Clean_Data
      PROCEDURE, PRIVATE :: Buffer_Old_State
      PROCEDURE, PRIVATE :: Calculate_Pole_Values
      PROCEDURE, PRIVATE :: Cross_Flux_Tube_Transport
      PROCEDURE, PRIVATE :: Cross_Flux_Tube_Transport2
      PROCEDURE, PRIVATE :: High_Latitude_Flux_Tube_Transport
      PROCEDURE, PRIVATE :: Test_Transport_Time_step
      PROCEDURE, PRIVATE :: Auroral_Precipitation
      PROCEDURE, PRIVATE :: FLIP_Wrapper
      PROCEDURE, PRIVATE :: Write_Electron_Density_to_HDF5
      PROCEDURE :: Calculate_Field_Line_Integrals

  END TYPE IPE_Plasma

  INTEGER, PARAMETER, PRIVATE    :: n_ion_species = 9
  REAL(prec), PARAMETER, PRIVATE :: safe_density_minimum = 1.0e+06_prec
  REAL(prec), PARAMETER, PRIVATE :: safe_temperature_minimum = 100.0_prec
  REAL(prec), PARAMETER, PRIVATE :: colfac = 1.5_prec
  REAL(prec), PARAMETER, PRIVATE :: qeoNao10 = 9.6489E7_prec        !  qe/m_e*1000 [C/g]
  REAL(prec), PARAMETER, PRIVATE :: qeomeo10 = 1.7588028E11_prec    ! qe/m_e*1000 [C/g]
  REAL(prec), PARAMETER, PRIVATE :: rmassinv_nop = 1.0_prec / NO_mass ! inverted rmass
  REAL(prec), PARAMETER, PRIVATE :: rmassinv_o1  = 1.0_prec / O_mass  ! inverted rmass
  REAL(prec), PARAMETER, PRIVATE :: rmassinv_o2  = 1.0_prec / O2_mass ! inverted rmass


  ! ::::::::::::::::: Cross_Flux_Tube_Transport  PARAMETERs ::::::::::::::::: !
  !
  ! transport_min_altitude - Perpendicular transport only occurs
  !                          for altitudes greater than
  !                          "transport_min_altitude"
  !
  !
  ! n_transport_species    - Number of ion species advected.
  !                          Only ion_densitities(1:n_transport_species)
  !                          are advected in the Cross_Flux_Tube_Transport
  !                          subroutine.
  !
  ! transport_highlat_lp   - The lp index corresponding to the last flux tube
  !                          adiabatic compression is not used.
  !
  ! perp_transport_max_lp  - Maximum of the lp-loop in
  !                          Cross_Flux_Tube_Transport. Plasma properties are
  !                          only advected on magnetic flux tubes with lp
  !                          between 1 and perp_transport_max_lp ( inclusive ).
  !
  ! ------------------------------------------------------------------------ !

  INTEGER, PARAMETER, PRIVATE    :: n_transport_species    = 4
  REAL(prec), PARAMETER, PRIVATE :: transport_min_altitude = 150000.0_prec
  INTEGER, PARAMETER , PRIVATE   :: transport_highlat_lp   = 30
  INTEGER, PARAMETER , PRIVATE   :: perp_transport_max_lp  = 151
! GHGM 
!  INTEGER, PARAMETER , PRIVATE   :: perp_transport_max_lp  = 130

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: !

  ! ::::::::::::::::::::::::::: Flip Parameters :::::::::::::::::::::::::::: !
  !
  !
  !
  !
  !
  !
  ! ------------------------------------------------------------------------ !
  REAL(dp), PARAMETER, PRIVATE :: DTMIN    = 1.0D0
  REAL(dp), PARAMETER, PRIVATE :: FPAS     = 0.0D0
  REAL(dp), PARAMETER, PRIVATE :: HEPRAT   = 9.0D-2
  REAL(dp), PARAMETER, PRIVATE :: COLFACX  = 1.7D0
  REAL(dp), PARAMETER, PRIVATE :: HPEQ     = 0.0D0
  ! IHEPLS,INPLS turn on diffusive solutions if > 0. no solution if 0, chemical equilibrium if < 0
  INTEGER, PARAMETER, PRIVATE  :: IHEPLS   = 1
  INTEGER, PARAMETER, PRIVATE  :: INPLS    = 1
  INTEGER, PARAMETER, PRIVATE  :: INNO     = 0
  integer :: istop
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: !


#ifdef HAVE_MPI
  INTEGER, ALLOCATABLE, PRIVATE :: ion_requestHandle(:)
  INTEGER, ALLOCATABLE, PRIVATE :: ion_requestStats(:,:)
#endif

CONTAINS

  SUBROUTINE Build_IPE_Plasma( plasma, nFluxTube, NLP, NMP, mp_low, mp_high, halo )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(out) :: plasma
    INTEGER, INTENT(in)              :: nFluxTube
    INTEGER, INTENT(in)              :: NLP
    INTEGER, INTENT(in)              :: NMP
    INTEGER, INTENT(in)              :: mp_low, mp_high, halo

      plasma % nFluxTube = nFluxTube
      plasma % NLP       = NLP
      plasma % NMP       = NMP
      plasma % mp_low    = mp_low
      plasma % mp_high   = mp_high
      plasma % mp_halo   = halo

      ALLOCATE( plasma % ion_densities(1:n_ion_species,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % ion_velocities(1:n_ion_species,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % ion_temperature(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_density(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_density2(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_velocity(1:3,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_temperature(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % ion_densities_old(1:n_ion_species,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % ion_velocities_old(1:n_ion_species,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % ion_temperature_old(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_density_old(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_velocity_old(1:3,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % electron_temperature_old(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                plasma % conductivities(1:6,1:2,1:NLP,1:NMP), &
                plasma % ionization_rates(1:4,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo) )

      plasma % ion_densities           = safe_density_minimum
      plasma % ion_velocities          = 0.0_prec
      plasma % ion_temperature         = safe_temperature_minimum
      plasma % electron_density        = safe_density_minimum
      plasma % electron_velocity       = 0.0_prec
      plasma % electron_temperature    = safe_temperature_minimum
      plasma % ionization_rates        = 0.0_prec
      plasma % conductivities          = 0.0_prec

      ALLOCATE( plasma % geo_ion_densities(1:n_ion_species,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ion_velocities(1:n_ion_species,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ion_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_density(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_velocity(1:3,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_electron_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_ionization_rates(1:4,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
                plasma % geo_tec(1:nlon_geo,1:nlat_geo), &
                plasma % geo_nmf2(1:nlon_geo,1:nlat_geo), &
                plasma % geo_hmf2(1:nlon_geo,1:nlat_geo) )

      plasma % geo_ion_densities           = 0.0_prec
      plasma % geo_ion_velocities          = 0.0_prec
      plasma % geo_ion_temperature         = 0.0_prec
      plasma % geo_electron_density        = 0.0_prec
      plasma % geo_electron_temperature    = 0.0_prec
      plasma % geo_electron_velocity       = 0.0_prec
      plasma % geo_ionization_rates        = 0.0_prec
      plasma % geo_tec                     = 0.0_prec
      plasma % geo_nmf2                    = 0.0_prec
      plasma % geo_hmf2                    = 0.0_prec

#ifdef HAVE_MPI
      ALLOCATE( ion_requestHandle(1:16), ion_requestStats(MPI_STATUS_SIZE,1:16) )
      ion_requestHandle = 0
      ion_requestStats  = 0
#endif

  END SUBROUTINE Build_IPE_Plasma
!
  SUBROUTINE Trash_IPE_Plasma( plasma )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma


    DEALLOCATE( plasma % ion_densities, &
                plasma % ion_velocities, &
                plasma % ion_temperature, &
                plasma % electron_density, &
                plasma % electron_density2, &
                plasma % electron_velocity, &
                plasma % electron_temperature, &
                plasma % ion_densities_old, &
                plasma % ion_velocities_old, &
                plasma % ion_temperature_old, &
                plasma % electron_density_old, &
                plasma % electron_velocity_old, &
                plasma % electron_temperature_old, &
                plasma % ionization_rates, &
                plasma % geo_ion_velocities, &
                plasma % geo_ion_temperature, &
                plasma % geo_electron_density, &
                plasma % geo_electron_velocity, &
                plasma % geo_electron_temperature, &
                plasma % geo_ionization_rates , &
                plasma % geo_tec, &
                plasma % geo_nmf2, &
                plasma % conductivities, &
                plasma % geo_hmf2 )

#ifdef HAVE_MPI
      DEALLOCATE( ion_requestHandle, ion_requestStats )
#endif

  END SUBROUTINE Trash_IPE_Plasma
!
  SUBROUTINE Update_IPE_Plasma_no_convection( plasma, grid, neutrals, forcing, time_tracker, mpi_layer, v_ExB, time_step, transport_time_step )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP,grid % mp_low:grid % mp_high)
    REAL(prec), INTENT(in)             :: time_step
    REAL(prec), INTENT(in)             :: transport_time_step
    INTEGER  :: nflag_t(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    INTEGER  :: nflag_d(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    ! Local
    INTEGER    :: n_transport_timesteps
    INTEGER    :: i, lp, mp, j
    REAL(prec) :: max_transport_convection_ratio_local
    REAL(prec) :: max_transport_convection_ratio
    REAL(prec) :: transport_time_step2
#ifdef HAVE_MPI
    INTEGER :: mpiError
#endif

        CALL plasma % Buffer_Old_State( grid )
        CALL plasma % Update_Halos( grid, mpi_layer )

#ifdef HAVE_MPI

        CALL MPI_WAITALL( 16, &
                         ion_requestHandle, &
                         ion_requestStats, &
                         mpiError)


#endif

!       CALL plasma % Cross_Flux_Tube_Transport( grid, v_ExB, transport_time_step2, mpi_layer )


      ! GHGM moved precipitation to here
      CALL plasma % FLIP_Wrapper( grid,         &
                                  neutrals,     &
                                  forcing,      &
                                  time_tracker, &
                                  time_step,nflag_t,nflag_d )

      !TWFANG, calculate field line integrals for dynamo solver
      CALL plasma % Calculate_Field_Line_Integrals(grid, neutrals, mpi_layer)

! GHGM
!     CALL plasma % Write_Electron_Density_to_HDF5( grid, mpi_layer, "IPE_Electron_Density."//time_tracker % DateStamp( )//".h5" )


  END SUBROUTINE Update_IPE_Plasma_no_convection


  SUBROUTINE Update_IPE_Plasma( plasma, grid, neutrals, forcing, time_tracker, mpi_layer, v_ExB, time_step )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP,grid % mp_low:grid % mp_high)
    REAL(prec), INTENT(in)             :: time_step
    INTEGER  :: nflag_t(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    INTEGER  :: nflag_d(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    ! Local
    INTEGER    :: n_transport_timesteps
    INTEGER    :: i, lp, mp, j
    REAL(prec) :: max_transport_convection_ratio_local
    REAL(prec) :: max_transport_convection_ratio
    REAL(prec) :: transport_time_step2
#ifdef HAVE_MPI
    INTEGER :: mpiError
#endif

! GHGM no transport for first call .....
      if(time_tracker % utime.gt.0.00001) then
      CALL plasma % Test_Transport_Time_step( grid, v_ExB, time_step, mpi_layer, &
                                              max_transport_convection_ratio_local )
!     write(6,7999) mpi_layer % rank_id , grid % mp_low, grid % mp_high, max_transport_convection_ratio_local, &
!                                          v_ExB(1,5:7,grid % mp_low), v_ExB(2,5:7,grid % mp_high)     
!7999 format('GHGM convect ratio ', 3i4 , f12.1, 6e10.2)

#ifdef HAVE_MPI
      CALL MPI_ALLREDUCE( max_transport_convection_ratio_local, &
                          max_transport_convection_ratio, &
                          1, &
                          mpi_layer % mpi_prec, &
                          MPI_MAX, &
                          mpi_layer % mpi_communicator, &
                          mpiError )

!     write(7000 + mpi_layer % rank_id,*) 'GHGM MAX ', max_transport_convection_ratio

      ! sub-stepping for advection
!     n_transport_timesteps = INT( time_step/transport_time_step )
      n_transport_timesteps = ceiling(max_transport_convection_ratio)

#else

      n_transport_timesteps = ceiling(max_transport_convection_ratio_local)

#endif

! GHGM minimun transport timestep of 1 minute .....
      if (n_transport_timesteps.lt.3) then
        if (mpi_layer % rank_id.eq.0) then
          write(6,*) 'GHGM TRANSPORT TIMESTEP LT 3 ', n_transport_timesteps
        endif
        n_transport_timesteps = 3
      endif
! GHGM minimun transport timestep of 1 minute .....

      transport_time_step2 = time_step / float(n_transport_timesteps)

!     write(6,1000) mpi_layer % rank_id, time_tracker % utime, &
!                                                      n_transport_timesteps, &
!                                                      transport_time_step2

!1000 format('GHGM timestep ', i4,f12.2,i4,f12.4)


      DO i = 1, n_transport_timesteps

!       write(1000 + mpi_layer % rank_id,*) ' GHGM TRANSPORT LOOP ', i , ' OF ', n_transport_timesteps

        CALL plasma % Buffer_Old_State( grid )
        CALL plasma % Update_Halos( grid, mpi_layer )

#ifdef HAVE_MPI

        CALL MPI_WAITALL( 16, &
                         ion_requestHandle, &
                         ion_requestStats, &
                         mpiError)


#endif

        CALL plasma % Cross_Flux_Tube_Transport( grid, v_ExB, transport_time_step2, mpi_layer )

      ENDDO
! GHGM no transport for first call .....
      endif  ! if(time_tracker % utime.gt.0.00001) then

      CALL plasma % Auroral_Precipitation( grid, &
                                           neutrals, &
                                           forcing, &
                                           time_tracker )

      CALL plasma % FLIP_Wrapper( grid,         &
                                  neutrals,     &
                                  forcing,      &
                                  time_tracker, &
                                  time_step,nflag_t,nflag_d )

      !TWFANG, calculate field line integrals for dynamo solver
      CALL plasma % Calculate_Field_Line_Integrals(grid, neutrals, mpi_layer)

! GHGM
!     CALL plasma % Write_Electron_Density_to_HDF5( grid, mpi_layer, "IPE_Electron_Density."//time_tracker % DateStamp( )//".h5" )


  END SUBROUTINE Update_IPE_Plasma
!
  SUBROUTINE Clean_Data( plasma, grid )
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    ! Local
    INTEGER :: i, j, lp, mp

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, plasma % NLP
          DO i = 1, grid % flux_tube_max(lp)

            DO j = 1, n_ion_species

              IF( plasma % ion_densities(j,i,lp,mp) < safe_density_minimum )THEN
                plasma % ion_densities(j,i,lp,mp) = safe_density_minimum
              ENDIF

            ENDDO

            IF( plasma % ion_temperature(i,lp,mp) < safe_temperature_minimum )THEN
              plasma % ion_temperature(i,lp,mp) = safe_temperature_minimum
            ENDIF

            IF( plasma % electron_temperature(i,lp,mp) < safe_temperature_minimum )THEN
              plasma % electron_temperature(i,lp,mp) = safe_temperature_minimum
            ENDIF

          ENDDO
        ENDDO
      ENDDO

  END SUBROUTINE Clean_Data

  SUBROUTINE Buffer_Old_State( plasma, grid )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    ! Local
    INTEGER :: mp, lp, i

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)

            plasma % ion_densities_old(1:n_ion_species,i,lp,mp)  = plasma % ion_densities(1:n_ion_species,i,lp,mp)
            plasma % ion_velocities_old(1:n_ion_species,i,lp,mp) = plasma % ion_velocities(1:n_ion_species,i,lp,mp)

            plasma % ion_temperature_old(i,lp,mp)      = plasma % ion_temperature(i,lp,mp)
            plasma % electron_density_old(i,lp,mp)     = plasma % electron_density(i,lp,mp)
            plasma % electron_temperature_old(i,lp,mp) = plasma % electron_temperature(i,lp,mp)

          ENDDO
        ENDDO
      ENDDO


  END SUBROUTINE Buffer_Old_State
!
  SUBROUTINE Update_Halos_IPE_Plasma( plasma, grid, mpi_layer )
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    INTEGER :: mp, mp_c, lp, i, iError, nhandles

#ifdef HAVE_MPI

    IF( mpi_layer % n_ranks == 1 )THEN

      ! Update the halos on the west end of the domain
      DO mp = plasma % mp_low-plasma % mp_halo, plasma % mp_low-1
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)

            mp_c = plasma % NMP + mp
            plasma % ion_densities_old(1:n_ion_species,i,lp,mp)  = plasma % ion_densities(1:n_ion_species,i,lp,mp_c)
            plasma % ion_velocities_old(1:n_ion_species,i,lp,mp) = plasma % ion_velocities(1:n_ion_species,i,lp,mp_c)

            plasma % ion_temperature_old(i,lp,mp)      = plasma % ion_temperature(i,lp,mp_c)
            plasma % electron_density_old(i,lp,mp)     = plasma % electron_density(i,lp,mp_c)
            plasma % electron_temperature_old(i,lp,mp) = plasma % electron_temperature(i,lp,mp_c)

          ENDDO
        ENDDO
      ENDDO

      ! Update the halos on the east end of the domain
      DO mp = plasma % mp_high+1, plasma % mp_high+plasma % mp_halo
        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)

            mp_c = mp - plasma % NMP + 1
            plasma % ion_densities_old(1:n_ion_species,i,lp,mp)  = plasma % ion_densities(1:n_ion_species,i,lp,mp_c)
            plasma % ion_velocities_old(1:n_ion_species,i,lp,mp) = plasma % ion_velocities(1:n_ion_species,i,lp,mp_c)

            plasma % ion_temperature_old(i,lp,mp)      = plasma % ion_temperature(i,lp,mp_c)
            plasma % electron_density_old(i,lp,mp)     = plasma % electron_density(i,lp,mp_c)
            plasma % electron_temperature_old(i,lp,mp) = plasma % electron_temperature(i,lp,mp_c)

          ENDDO
        ENDDO
      ENDDO

    ELSE  ! mpi_layer % n_ranks == 1

    ! Here, we send data to the rank on the right and receive data from the rank
    ! on the left

      ! Each rank will receive data in the mp_low-halo:mp_low range in the
      ! _old buffer. Data will be sent from the current ion state arrays

      CALL MPI_IRECV( plasma % ion_densities_old(:,:,:,plasma % mp_low-plasma % mp_halo), &
                     n_ion_species*plasma % nFluxTube*plasma % NLP, &
                     mpi_layer % MPI_PREC,   &
                     mpi_layer % neighbor_rank(1), 0,  &
                     mpi_layer % mpi_communicator,   &
                     ion_requestHandle(1), iError )

      CALL MPI_ISEND( plasma % ion_densities(:,:,:,plasma % mp_high), &
                      n_ion_species*plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 0,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(2), iError )


      CALL MPI_IRECV( plasma % ion_velocities_old(:,:,:,plasma % mp_low-plasma % mp_halo), &
                     n_ion_species*plasma % nFluxTube*plasma % NLP, &
                     mpi_layer % MPI_PREC,   &
                     mpi_layer % neighbor_rank(1), 1,  &
                     mpi_layer % mpi_communicator,   &
                     ion_requestHandle(3), iError )

      CALL MPI_ISEND( plasma % ion_velocities(:,:,:,:plasma % mp_high), &
                      n_ion_species*plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 1,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(4), iError )

      CALL MPI_IRECV( plasma % ion_temperature_old(:,:,plasma % mp_low-plasma % mp_halo), &
                     plasma % nFluxTube*plasma % NLP, &
                     mpi_layer % MPI_PREC,   &
                     mpi_layer % neighbor_rank(1), 2,  &
                     mpi_layer % mpi_communicator,   &
                     ion_requestHandle(5), iError )

      CALL MPI_ISEND( plasma % ion_temperature(:,:,plasma % mp_high), &
                      plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 2,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(6), iError )

      CALL MPI_IRECV( plasma % electron_temperature_old(:,:,plasma % mp_low-plasma % mp_halo), &
                     plasma % nFluxTube*plasma % NLP, &
                     mpi_layer % MPI_PREC,   &
                     mpi_layer % neighbor_rank(1), 3,  &
                     mpi_layer % mpi_communicator,   &
                     ion_requestHandle(7), iError )

      CALL MPI_ISEND( plasma % electron_temperature(:,:,plasma % mp_high), &
                      plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 3,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(8), iError )

   ! Here, we send data to the rank on the left and receive data from the rank
   ! on the right

     ! Each rank will receive data in the mp_high:mp_low+plasma % mp_halo range in the
     ! _old buffer. Data will be sent from the current ion state arrays
      CALL MPI_IRECV( plasma % ion_densities_old(:,:,:,plasma % mp_high+plasma % mp_halo), &
                     n_ion_species*plasma % nFluxTube*plasma % NLP, &
                     mpi_layer % MPI_PREC,   &
                     mpi_layer % neighbor_rank(2), 4,  &
                     mpi_layer % mpi_communicator,   &
                     ion_requestHandle(9), iError )

      CALL MPI_ISEND( plasma % ion_densities(:,:,:,plasma % mp_low), &
                      n_ion_species*plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(1), 4,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(10), iError )

      CALL MPI_IRECV( plasma % ion_velocities_old(:,:,:,plasma % mp_high+plasma % mp_halo), &
                      n_ion_species*plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 5,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(11), iError )

      CALL MPI_ISEND( plasma % ion_velocities(:,:,:,plasma % mp_low), &
                      n_ion_species*plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(1), 5,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(12), iError )

      CALL MPI_IRECV( plasma % ion_temperature_old(:,:,plasma % mp_high+plasma % mp_halo), &
                      plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 6,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(13), iError )

      CALL MPI_ISEND( plasma % ion_temperature(:,:,plasma % mp_low), &
                      plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(1), 6,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(14), iError )

      CALL MPI_IRECV( plasma % electron_temperature_old(:,:,plasma % mp_high+plasma % mp_halo), &
                      plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(2), 7,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(15), iError )

      CALL MPI_ISEND( plasma % electron_temperature(:,:,plasma % mp_low), &
                      plasma % nFluxTube*plasma % NLP, &
                      mpi_layer % MPI_PREC,   &
                      mpi_layer % neighbor_rank(1), 7,  &
                      mpi_layer % mpi_communicator,   &
                      ion_requestHandle(16), iError )


    ENDIF ! mpi_layer % n_ranks == 1

#else
    ! In the serial case, only the periodic boundary conditions need to be
    ! updated


      ! Update the halos on the west end of the domain
      mp = 0
      mp_c = plasma % NMP
      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          plasma % ion_densities_old(1:n_ion_species,i,lp,mp)  = plasma % ion_densities(1:n_ion_species,i,lp,mp_c)
          plasma % ion_velocities_old(1:n_ion_species,i,lp,mp) = plasma % ion_velocities(1:n_ion_species,i,lp,mp_c)

          plasma % ion_temperature_old(i,lp,mp)      = plasma % ion_temperature(i,lp,mp_c)
          plasma % electron_density_old(i,lp,mp)     = plasma % electron_density(i,lp,mp_c)
          plasma % electron_temperature_old(i,lp,mp) = plasma % electron_temperature(i,lp,mp_c)

        ENDDO
      ENDDO

      ! Update the halos on the east end of the domain
      mp   = plasma % NMP + 1
      mp_c = 1
      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          plasma % ion_densities_old(1:n_ion_species,i,lp,mp)  = plasma % ion_densities(1:n_ion_species,i,lp,mp_c)
          plasma % ion_velocities_old(1:n_ion_species,i,lp,mp) = plasma % ion_velocities(1:n_ion_species,i,lp,mp_c)

          plasma % ion_temperature_old(i,lp,mp)      = plasma % ion_temperature(i,lp,mp_c)
          plasma % electron_density_old(i,lp,mp)     = plasma % electron_density(i,lp,mp_c)
          plasma % electron_temperature_old(i,lp,mp) = plasma % electron_temperature(i,lp,mp_c)

        ENDDO
      ENDDO

#endif


  END SUBROUTINE Update_Halos_IPE_Plasma
!
  SUBROUTINE Calculate_Pole_Values( plasma, &
                                    grid, &
                                    mpi_layer, &
                                    ion_densities_pole_value, &
                                    ion_temperature_pole_value, &
                                    ion_velocities_pole_value, &
                                    electron_temperature_pole_value )

    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(in)   :: plasma
    TYPE( IPE_Grid ), INTENT(in)      :: grid
    TYPE( IPE_MPI_Layer ), INTENT(in) :: mpi_layer
    integer, PARAMETER :: n_conv_spec = 4   ! number of ions for ExB convection ( = 4)
    REAL(prec), INTENT(out)           :: ion_densities_pole_value(1:n_conv_spec, 1:plasma % nFluxTube)
    REAL(prec), INTENT(out)           :: ion_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec), INTENT(out)           :: ion_velocities_pole_value(1:n_conv_spec, 1:plasma % nFluxTube)
    REAL(prec), INTENT(out)           :: electron_temperature_pole_value(1:plasma % nFluxTube)
    ! Local
    REAL(prec) :: ion_densities_pole_value_loc(1:n_conv_spec, 1:plasma % nFluxTube)
    REAL(prec) :: ion_temperature_pole_value_loc(1:plasma % nFluxTube)
    REAL(prec) :: ion_velocities_pole_value_loc(1:n_conv_spec, 1:plasma % nFluxTube)
    REAL(prec) :: electron_temperature_pole_value_loc(1:plasma % nFluxTube)
    INTEGER    :: i, lp, mp, j
#ifdef HAVE_MPI
    INTEGER    :: error
#endif

      ion_densities_pole_value_loc      = 0.0_prec
      ion_temperature_pole_value_loc    = 0.0_prec
      ion_velocities_pole_value_loc     = 0.0_prec
      electron_temperature_pole_value_loc = 0.0_prec

      DO mp = plasma % mp_low, plasma % mp_high
        DO i = 1, grid % flux_tube_max(1)

          DO j = 1, n_conv_spec

            ion_densities_pole_value_loc(j,i) = ion_densities_pole_value_loc(j,i) + plasma % ion_densities(j,i,1,mp)
            ion_velocities_pole_value_loc(j,i) = ion_velocities_pole_value_loc(j,i) + plasma % ion_velocities(j,i,1,mp)

          ENDDO

          ion_temperature_pole_value_loc(i) = ion_temperature_pole_value_loc(i) + plasma % ion_temperature(i,1,mp)
          electron_temperature_pole_value_loc(i) = electron_temperature_pole_value_loc(i) + plasma % electron_temperature(i,1,mp)

        ENDDO
      ENDDO
#ifdef HAVE_MPI
      CALL MPI_ALLREDUCE( ion_densities_pole_value_loc, &
                          ion_densities_pole_value, &
                          n_conv_spec*grid % flux_tube_max(1), &
                          mpi_layer % mpi_prec, &
                          MPI_SUM, &
                          mpi_layer % mpi_communicator, &
                          error )

      CALL MPI_ALLREDUCE( ion_temperature_pole_value_loc, &
                          ion_temperature_pole_value, &
                          grid % flux_tube_max(1), &
                          mpi_layer % mpi_prec, &
                          MPI_SUM, &
                          mpi_layer % mpi_communicator, &
                          error )

      CALL MPI_ALLREDUCE( ion_velocities_pole_value_loc, &
                          ion_velocities_pole_value, &
                          n_conv_spec*grid % flux_tube_max(1), &
                          mpi_layer % mpi_prec, &
                          MPI_SUM, &
                          mpi_layer % mpi_communicator, &
                          error )

      CALL MPI_ALLREDUCE( electron_temperature_pole_value_loc, &
                          electron_temperature_pole_value, &
                          grid % flux_tube_max(1), &
                          mpi_layer % mpi_prec, &
                          MPI_SUM, &
                          mpi_layer % mpi_communicator, &
                          error )

      ion_densities_pole_value        = ion_densities_pole_value/REAL( plasma % NMP )
      ion_temperature_pole_value      = ion_temperature_pole_value/REAL( plasma % NMP )
      electron_temperature_pole_value = electron_temperature_pole_value/REAL( plasma % NMP )
      ion_velocities_pole_value       = ion_velocities_pole_value/REAL( plasma % NMP )
#else

      ion_densities_pole_value        = ion_densities_pole_value_loc/REAL( plasma % NMP )
      ion_temperature_pole_value      = ion_temperature_pole_value_loc/REAL( plasma % NMP )
      electron_temperature_pole_value = electron_temperature_pole_value_loc/REAL( plasma % NMP )
      ion_velocities_pole_value       = ion_velocities_pole_value_loc/REAL( plasma % NMP )
#endif


  END SUBROUTINE Calculate_Pole_Values


  SUBROUTINE Write_Electron_Density_to_HDF5( plasma, grid, mpi_layer, filename )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)      :: grid
    TYPE( IPE_MPI_Layer ), INTENT(in) :: mpi_layer
    CHARACTER(*), INTENT(in)      :: filename
    ! Local
    CHARACTER(100)   :: groupname
    CHARACTER(13)    :: timeStampString
    CHARACTER(10)    :: zoneID
    INTEGER          :: mp, iEl, N, rank, m_rank, error, istat, nEl, elID
    INTEGER(HSIZE_T) :: dimensions(1:3), global_dimensions(1:3)
    INTEGER(HSIZE_T) :: starts(1:3), counts(1:3), strides(1:3)
    INTEGER(HID_T)   :: file_id, memspace, dataset_id, filespace
    INTEGER(HID_T)   :: model_group_id
    INTEGER(HID_T)   :: plist_id
    REAL(prec)       :: electron_density(1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high)


    IF( mpi_layer % rank_id == 0 )THEN
      PRINT*, '  Writing Electron Density file : '//TRIM(filename)
    ENDIF

#ifdef HAVE_MPI

    CALL MPI_BARRIER( mpi_layer % mpi_communicator, error )
    rank = 3
    ! Local Dimensions
    dimensions = (/ grid % nFluxTube, grid % NLP, 1 /)
    global_dimensions = (/ grid % nFluxTube, grid % NLP, grid % NMP /)


#else

    rank = 3
    ! Local Dimensions
    dimensions = (/ grid % nFluxTube, grid % NLP, grid % NMP /)
    global_dimensions = (/ grid % nFluxTube, grid % NLP, grid % NMP /)

#endif


    CALL h5open_f(error)

#ifdef HAVE_MPI
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, mpi_layer % mpi_communicator, MPI_INFO_NULL, error)

    ! Create a new file using default properties.
    CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
    CALL h5pclose_f(plist_id, error)

    ! Create a dataspace for the global dataset
    CALL h5screate_simple_f(rank, global_dimensions, filespace, error)
    CALL h5screate_simple_f(rank, dimensions, memspace, error)

    ! Set the data creation mode to CHUNK
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    CALL h5pset_chunk_f(plist_id, rank, dimensions, error)

#else

    ! Create a new file using default properties.
    CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_id, error)
    CALL h5screate_simple_f(rank, dimensions, memspace, error)

#endif



    ! Create groups
    groupname = "/apex"
    CALL h5gcreate_f( file_id, TRIM(groupname), model_group_id, error )
    IF( error /= 0 ) STOP


    ! Plasma data
    ! O+
    electron_density(1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) = &
             plasma % ion_densities(1,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(2,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(3,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(4,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(5,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(6,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(7,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(8,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high) + &
             plasma % ion_densities(9,1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high)

    CALL Add_Electron_Variable_to_HDF5( file_id, "/apex/electron_density",&
                               electron_density(1:grid % nFluxTube,1:grid % NLP,mpi_layer % mp_low:mpi_layer % mp_high), &
                               filespace, memspace, plist_id, dimensions,&
                               grid % nFluxTube, grid % NLP, &
                               grid % NMP, &
                               mpi_layer % mp_low, &
                               mpi_layer % mp_high, error )

#ifdef HAVE_MPI
    CALL h5pclose_f( plist_id, error )
    CALL h5sclose_f( filespace, error )
#endif

    CALL h5gclose_f( model_group_id, error )
    CALL h5sclose_f( memspace, error )
    CALL h5fclose_f( file_id, error )

    CALL h5close_f( error )



  END SUBROUTINE Write_Electron_Density_to_HDF5


  SUBROUTINE Add_Electron_Variable_to_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, n_fluxtube, NLP, NMP, mp_low, mp_high, error )
    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: filespace
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HID_T), INTENT(in)   :: plist_id
    INTEGER, INTENT(in)          :: n_fluxtube, NLP, NMP
    INTEGER, INTENT(in)          :: mp_low, mp_high
    REAL(prec), INTENT(in)       :: variable(1:n_fluxtube,1:NLP,mp_low:mp_high)
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:3)
    INTEGER, INTENT(out)         :: error
    ! Local
    INTEGER(HID_T) :: dataset_id
    INTEGER        :: mp
    INTEGER(HSIZE_T)        :: starts(1:3), counts(1:3), strides(1:3)



#ifdef HAVE_MPI

    CALL h5dcreate_f( file_id, TRIM(variable_name), H5T_IEEE_F64LE, filespace, dataset_id, error, plist_id)

    DO mp = mp_low, mp_high
      starts = (/ 0, 0, mp-1 /)
      counts = (/ 1, 1, 1 /)
      strides = (/ 1, 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )

      CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                       variable(1:n_fluxtube,1:NLP,mp), &
                       dimensions, error, memspace, filespace )

    ENDDO

#else

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                       H5T_IEEE_F64LE, memspace, dataset_id, error)
    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)

#endif

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_Electron_Variable_to_HDF5

  SUBROUTINE Test_Transport_Time_step( plasma, grid, v_ExB, time_step, mpi_layer, &
                                       max_transport_convection_ratio )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec), INTENT(in)             :: time_step
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    REAL(prec) :: colat_90km(1:grid % NLP)
    REAL(prec) :: phi_t0 !magnetic longitude,phi[rad] at t0(previous time step)
    REAL(prec) :: theta_t0 !magnetic latitude,theta[rad] at t0
    REAL(prec) :: coslam, sinim
    INTEGER    :: mp, lp, i, lpx, mpx, jth
    REAL(prec), PARAMETER :: rad_to_deg = 57.295779513
    REAL(prec) :: transport_convection_ratio(1:perp_transport_max_lp, grid % mp_low:grid % mp_high)
    REAL(prec) :: max_transport_convection_ratio
    REAL(prec) :: longitude_spacing, r

      colat_90km(1:grid % NLP) = grid % magnetic_colatitude(1,1:grid % NLP)
      r = earth_radius + 90000.0_prec
      longitude_spacing = 360.0_prec / REAL( plasma % NMP )

!     write(7000 + mpi_layer % rank_id, *) 'mp_low mp_high ', plasma % mp_low, plasma % mp_high

      DO 100 mp = plasma % mp_low, plasma % mp_high
        DO 200 lp = 1, perp_transport_max_lp


          transport_convection_ratio(lp,mp) = (abs(v_ExB(1,lp,mp)*time_step/(r*sin( colat_90km(lp)))) * rad_to_deg) / longitude_spacing
!         write(7000 + mpi_layer % rank_id, *) mp , lp , transport_convection_ratio(lp,mp)

 200    CONTINUE
 100  CONTINUE

      max_transport_convection_ratio = maxval(transport_convection_ratio)
!     write(7000 + mpi_layer % rank_id, *) 'Ratio ',maxval(transport_convection_ratio), &
!                                          maxloc(transport_convection_ratio)

  END SUBROUTINE Test_Transport_Time_step


  SUBROUTINE Cross_Flux_Tube_Transport( plasma, grid, v_ExB, time_step, mpi_layer )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec), INTENT(in)             :: time_step
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    INTEGER, PARAMETER :: n_conv_spec = 4   ! number of ions for ExB convection ( = 4)
    REAL(prec) :: ion_densities_pole_value(1:n_conv_spec, 1:plasma % nFluxTube)
    REAL(prec) :: ion_velocities_pole_value(1:n_conv_spec, 1:plasma % nFluxTube)
    REAL(prec) :: ion_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: electron_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: colat_90km(1:grid % NLP)
    REAL(prec) :: phi_t0 !magnetic longitude,phi[rad] at t0(previous time step)
    REAL(prec) :: theta_t0 !magnetic latitude,theta[rad] at t0
    REAL(prec) :: q_value
    REAL(prec) :: phi_i(1:2)
    REAL(prec) :: lp_comp_weight(1:2)
    REAL(prec) :: mp_comp_weight(1:2)
    REAL(prec) :: i_comp_weight(1:2)
    REAL(prec) :: q_int(1:2)
    REAL(prec) :: ion_densities_int(1:n_conv_spec)
    REAL(prec) :: ion_velocities_int(1:n_conv_spec)
    REAL(prec) :: ion_temperature_int
    REAL(prec) :: electron_temperature_int
    REAL(prec) :: r, B_int, max_phi, ksi_fac
    REAL(prec) :: B(1:2,1:2), velocity(1:n_conv_spec,1:2,1:2), density(1:n_conv_spec,1:2,1:2), temperature(1:2,1:2), e_temperature(1:2,1:2)
    REAL(prec) :: coslam, sinim
    INTEGER    :: lp_t0(1:2)
    INTEGER    :: mp_t0(1:2)
    INTEGER    :: lp_min, mp_min, isouth, inorth, ii, ispecial
    INTEGER    :: mp, lp, i, lpx, mpx, jth
    INTEGER    :: i_min(1:2)
    REAL(prec), PARAMETER :: rad_to_deg = 57.295779513

      CALL plasma % Calculate_Pole_Values( grid,                       &
                                           mpi_layer,                  &
                                           ion_densities_pole_value,   &
                                           ion_temperature_pole_value, &
                                           ion_velocities_pole_value,  &
                                           electron_temperature_pole_value )


      colat_90km(1:grid % NLP) = grid % magnetic_colatitude(1,1:grid % NLP)
      r = earth_radius + 90000.0_prec

      DO 100 mp = plasma % mp_low, plasma % mp_high
        DO 200 lp = 1, perp_transport_max_lp

          phi_t0   = grid % magnetic_longitude(mp) - v_ExB(1,lp,mp)*time_step/(r*sin( colat_90km(lp) ) )

          coslam = cos( half_pi - grid % magnetic_colatitude(1,lp) )
          sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
          theta_t0 = colat_90km(lp) - v_ExB(2,lp,mp)*time_step/(r*sinim)

          ! If a Lagrangian trajectory crosses the equator, we clip the colatitude
          ! so that the point resides at the equator.
          IF( theta_t0 > colat_90km( grid % NLP ) )THEN ! NLP ==> Equator
            theta_t0 = colat_90km( grid % NLP )
          ENDIF



          ! lp_min is the nearest point to theta_t0 that has a larger colat value
          lp_min = 0
          IF( lp == 1 )THEN
            ! Check poleward
            IF( theta_t0 < colat_90km(1) )THEN
              lp_min = 1
            ELSE
              lp_min = 2
            ENDIF

          ELSE
            ! Check poleward
            IF( theta_t0 <= colat_90km(lp) .AND. theta_t0 >= colat_90km(lp-1) )THEN
              lp_min = lp
            ! Check equatorward
            ELSEIF( theta_t0 <= colat_90km(lp+1) .AND. theta_t0 >= colat_90km(lp) )THEN
              lp_min = lp+1
            ENDIF
          ENDIF

          IF( lp_min == 0 )THEN
            write(6,*) 'GHGM Need to expand search ', mp , lp
!           STOP 'lp-CFL > 1'
            ! Check poleward
            IF( theta_t0 <= colat_90km(lp-1) .AND. theta_t0 >= colat_90km(lp-2) )THEN
              lp_min = lp - 1
              write(6,*) 'GHGM 2 poleward ', mp , lp
            ! Check equatorward
            ELSEIF( theta_t0 <= colat_90km(lp+2) .AND. theta_t0 >= colat_90km(lp+1) )THEN
              lp_min = lp + 2
              write(6,*) 'GHGM 2 equatorward ', mp , lp
            ENDIF
          ENDIF

          IF( lp_min == 0 )THEN
            write(6,*) 'GHGM Still didnt get it trying again ', mp , lp
            IF( theta_t0 <= colat_90km(lp-2) .AND. theta_t0 >= colat_90km(lp-3) )THEN
              lp_min = lp - 2
              write(6,*) 'GHGM 3 poleward ', mp , lp
            ! Check equatorward
            ELSEIF( theta_t0 <= colat_90km(lp+3) .AND. theta_t0 >= colat_90km(lp+2) )THEN
              lp_min = lp + 3
              write(6,*) 'GHGM 3 equatorward ', mp , lp
            ENDIF
          ENDIF

          IF( lp_min == 0 )THEN
            write(6,*) 'GHGM Still didnt get it trying again AGAIN ', mp , lp
            IF( theta_t0 <= colat_90km(lp-3) .AND. theta_t0 >= colat_90km(lp-4) )THEN
              lp_min = lp - 3
              write(6,*) 'GHGM 4 poleward ', mp , lp
            ! Check equatorward
            ELSEIF( theta_t0 <= colat_90km(lp+4) .AND. theta_t0 >= colat_90km(lp+3) )THEN
              lp_min = lp + 4
              write(6,*) 'GHGM 4 equatorward ', mp , lp
            ENDIF
          ENDIF

          IF( lp_min == 0 )THEN
            write(6,*) 'GHGM OK I GIVE UP ', mp , lp
          ENDIF

          lp_t0(1) = lp_min-1
          lp_t0(2) = lp_min

          mp_min = 0
          IF( phi_t0 <= grid % magnetic_longitude(mp) .AND. phi_t0 >= grid % magnetic_longitude(mp-1) )THEN
            mp_min = mp
          ELSEIF( phi_t0 <= grid % magnetic_longitude(mp+1) .AND. phi_t0 >= grid % magnetic_longitude(mp) )THEN
            mp_min = mp+1
          ENDIF

          IF( mp_min == 0 )THEN
!           STOP 'mp-CFL > 1'
          ENDIF

          mp_t0(1) = mp_min-1
          mp_t0(2) = mp_min

          phi_i(1) = grid % magnetic_longitude(mp_t0(1))
          phi_i(2) = grid % magnetic_longitude(mp_t0(2))
          mp_comp_weight(1) =  ( phi_t0 - phi_i(2) )/( phi_i(1)-phi_i(2) )
          mp_comp_weight(2) = -( phi_t0 - phi_i(1) )/( phi_i(1)-phi_i(2) )

          IF( lp_min == 1 )THEN  ! lp_min == 1

            DO i = 1, grid % flux_tube_max(lp)

              plasma % ion_densities(1:n_conv_spec,i,lp,mp) = ion_densities_pole_value(1:n_conv_spec,i)
              plasma % ion_velocities(1:n_conv_spec,i,lp,mp) = ion_velocities_pole_value(1:n_conv_spec,i)
              plasma % ion_temperature(i,lp,mp)               = ion_temperature_pole_value(i)
              plasma % electron_temperature(i,lp,mp)          = electron_temperature_pole_value(i)

            ENDDO

          ELSE ! lp_min =/= 1 ....

              
              if(lp_t0(2).eq.0) then
                lp_t0(1) = 1
                lp_t0(2) = 2
                write(6,*) 'GHGM LP_T0 ',mp,lp,v_ExB(1,lp,mp), v_ExB(2,lp,mp)
              endif

            lp_comp_weight(1) =  ( theta_t0 - colat_90km(lp_t0(2)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )
            lp_comp_weight(2) = -( theta_t0 - colat_90km(lp_t0(1)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )

            DO 300 i = 1, grid % flux_tube_max(lp)

! only do convection for points over 200km.....

              if(grid % altitude(i,lp).ge.200000.0_prec) then

              ! q interpolation
              q_value = grid % q_factor(i,lp,mp)
              DO mpx = 1, 2
                DO lpx = 1,2

                  B(lpx,mpx) = 0.0_prec
                  density(1:n_conv_spec,lpx,mpx) = 0.0_prec
                  velocity(1:n_conv_spec,lpx,mpx) = 0.0_prec
                  temperature(lpx,mpx) = 0.0_prec
                  e_temperature(lpx,mpx) = 0.0_prec

                  ! We assume by default the last flux tube point should be used
                  ! and we setup weights that prolong the last flux tube value
                  isouth = grid % flux_tube_max(lp_t0(lpx))
                  inorth = isouth
                  i_comp_weight(1) = 1.0_prec
                  i_comp_weight(2) = 0.0_prec
                  ! Search for the nearest q_factor
                  DO ii = 2, grid % flux_tube_max(lp_t0(lpx))
                    IF(  grid % q_factor(ii, lp_t0(lpx), mp_t0(mpx)) < q_value )THEN

                      isouth   = ii
                      inorth   = ii-1
                      q_int(1) = grid % q_factor(isouth, lp_t0(lpx), mp_t0(mpx))
                      q_int(2) = grid % q_factor(inorth, lp_t0(lpx), mp_t0(mpx))

                      i_comp_weight(1) = ( q_value - q_int(2) )/( q_int(1) - q_int(2) )
                      i_comp_weight(2) = -( q_value - q_int(1) )/( q_int(1) - q_int(2) )
                      EXIT

                    ENDIF
                  ENDDO

                  B(lpx,mpx) = grid % magnetic_field_strength(isouth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(1) +&
                               grid % magnetic_field_strength(inorth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(2)

                  density(1:n_conv_spec,lpx,mpx) = plasma % ion_densities_old(1:n_conv_spec,isouth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(1) +&
                                                     plasma % ion_densities_old(1:n_conv_spec,inorth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(2)

                  velocity(1:n_conv_spec,lpx,mpx) = plasma % ion_velocities_old(1:n_conv_spec,isouth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(1) +&
                                                      plasma % ion_velocities_old(1:n_conv_spec,inorth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(2)

                  temperature(lpx,mpx) = plasma % ion_temperature_old(isouth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(1) +&
                                         plasma % ion_temperature_old(inorth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(2)

                  e_temperature(lpx,mpx) = plasma % electron_temperature_old(isouth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(1) +&
                                           plasma % electron_temperature_old(inorth,lp_t0(lpx),mp_t0(mpx))*i_comp_weight(2)




                ENDDO
              ENDDO

              ! mp,lp interpolation
              ! Reduction over mpx, lpx
              ion_densities_int   = 0.0_prec
              ion_velocities_int   = 0.0_prec
              ion_temperature_int = 0.0_prec
              electron_temperature_int = 0.0_prec
              B_int               = 0.0_prec
              DO mpx = 1, 2
                DO lpx = 1, 2
                  B_int = B_int + B(lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)

                  ion_densities_int(1:n_conv_spec) = ion_densities_int(1:n_conv_spec) + density(1:n_conv_spec,lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)
                  ion_velocities_int(1:n_conv_spec) = ion_velocities_int(1:n_conv_spec) + velocity(1:n_conv_spec,lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)

                  ion_temperature_int = ion_temperature_int + temperature(lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)
                  electron_temperature_int = electron_temperature_int + e_temperature(lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)
                ENDDO
              ENDDO

              IF( lp <= transport_highlat_lp )THEN
                ksi_fac = 1.0_prec
              ELSE
                ksi_fac = grid % magnetic_field_strength(i,lp,mp)/B_int
              ENDIF

              plasma % ion_densities(1:n_conv_spec,i,lp,mp) = ion_densities_int(1:n_conv_spec)*( ksi_fac**2 )
              plasma % ion_velocities(1:n_conv_spec,i,lp,mp) = ion_velocities_int(1:n_conv_spec)
!             plasma % ion_temperature(i,lp,mp) = ion_temperature_int
!             plasma % electron_temperature(i,lp,mp) = electron_temperature_int
              plasma % ion_temperature(i,lp,mp) = ion_temperature_int*( ksi_fac**(4.0_prec/3.0_prec) )
              plasma % electron_temperature(i,lp,mp) = electron_temperature_int*( ksi_fac**(4.0_prec/3.0_prec) )



              endif     ! if(grid % altitude(i,lp).ge.200000.0_prec)

 300        CONTINUE  !  i = 1, grid % flux_tube_max(lp)

          ENDIF ! lp_min =/= 1

 200    CONTINUE
 100  CONTINUE

  END SUBROUTINE Cross_Flux_Tube_Transport


  SUBROUTINE High_Latitude_Flux_Tube_Transport( plasma, grid, v_ExB, time_step, mpi_layer, time_tracker, forcing )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec), INTENT(in)             :: time_step
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    REAL(prec) :: v_ExB_local(1:3,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec) :: ion_densities_pole_value(1:n_ion_species, 1:plasma % nFluxTube)
    REAL(prec) :: ion_velocities_pole_value(1:n_ion_species, 1:plasma % nFluxTube)
    REAL(prec) :: ion_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: electron_temperature_pole_value(1:plasma % nFluxTube)
    REAL(prec) :: colat_90km(1:grid % NLP)
    REAL(prec) :: phi_t0 !magnetic longitude,phi[rad] at t0(previous time step)
    REAL(prec) :: theta_t0 !magnetic latitude,theta[rad] at t0
    REAL(prec) :: q_value
    REAL(prec) :: phi_i(1:2)
    REAL(prec) :: lp_comp_weight(1:2)
    REAL(prec) :: mp_comp_weight(1:2)
    REAL(prec) :: i_comp_weight(1:2)
    REAL(prec) :: q_int(1:2)
    REAL(prec) :: ion_densities_int(1:n_ion_species)
    REAL(prec) :: ion_temperature_int
    REAL(prec) :: r, B_int, max_phi, ksi_fac
    REAL(prec) :: coslam, sinim
    REAL(prec) :: q_val(1:plasma % nFluxTube)
    REAL(prec) :: q_val_N_W(1:plasma % nFluxTube)
    REAL(prec) :: q_val_N_E(1:plasma % nFluxTube)
    REAL(prec) :: q_val_S_W(1:plasma % nFluxTube)
    REAL(prec) :: q_val_S_E(1:plasma % nFluxTube)

    REAL(prec) :: ion_den_NW(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_NW(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_NW(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_NW(1:plasma % nFluxTube)

    REAL(prec) :: ion_den_NE(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_NE(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_NE(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_NE(1:plasma % nFluxTube)

    REAL(prec) :: ion_den_SW(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_SW(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_SW(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_SW(1:plasma % nFluxTube)

    REAL(prec) :: ion_den_SE(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_SE(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_SE(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_SE(1:plasma % nFluxTube)

    REAL(prec) :: ion_den_N_W(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_N_W(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_N_W(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_N_W(1:plasma % nFluxTube)
    REAL(prec) :: ion_den_N_E(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_N_E(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_N_E(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_N_E(1:plasma % nFluxTube)
    REAL(prec) :: ion_den_S_W(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_S_W(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_S_W(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_S_W(1:plasma % nFluxTube)
    REAL(prec) :: ion_den_S_E(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_vel_S_E(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temp_S_E(1:plasma % nFluxTube)
    REAL(prec) :: electron_temp_S_E(1:plasma % nFluxTube)

    REAL(prec) :: lon_fac, lat_fac
    REAL(prec) :: ion_densities_N(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_densities_S(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_velocities_N(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_velocities_S(1:4,1:plasma % nFluxTube)
    REAL(prec) :: ion_temperature_N(1:plasma % nFluxTube)
    REAL(prec) :: ion_temperature_S(1:plasma % nFluxTube)
    REAL(prec) :: electron_temperature_N(1:plasma % nFluxTube)
    REAL(prec) :: electron_temperature_S(1:plasma % nFluxTube)

    REAL(prec) :: ion_densities_t0(1:4,1:plasma % nFluxTube,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec) :: ion_velocities_t0(1:4,1:plasma % nFluxTube,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec) :: ion_temperature_t0(1:plasma % nFluxTube,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec) :: electron_temperature_t0(1:plasma % nFluxTube,1:grid % NLP, grid % mp_low:grid % mp_high)

    INTEGER    :: lp_N, lp_S
    INTEGER    :: mp_W, mp_E
    INTEGER    :: lp_min, mp_min, isouth, inorth, ispecial
    INTEGER    :: mp, lp, i, ii, lpx, mpx, jth
    INTEGER    :: i_min(1:2)
    INTEGER    :: high_lat_convection_max_lp
    REAL(prec), PARAMETER :: r2d = 57.295779513
    INTEGER    :: dont_convect(1:grid % NLP, grid % mp_low:grid % mp_high)

          dont_convect(:,:) = 0
          ion_densities_t0(:,:,:,:) = 0.0
          ion_velocities_t0(:,:,:,:) = 0.0
          ion_temperature_t0(:,:,:) = 0.0
          electron_temperature_t0(:,:,:) = 0.0


      CALL plasma % Calculate_Pole_Values( grid,                       &
                                           mpi_layer,                  &
                                           ion_densities_pole_value,   &
                                           ion_temperature_pole_value, &
                                           ion_velocities_pole_value,  &
                                           electron_temperature_pole_value )


      colat_90km(1:grid % NLP) = grid % magnetic_colatitude(1,1:grid % NLP)
      r = earth_radius + 90000.0_prec
      high_lat_convection_max_lp = 30

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, high_lat_convection_max_lp

          v_exb_local(1,lp,mp) = v_exb(1,lp,mp)
          v_exb_local(2,lp,mp) = v_exb(2,lp,mp)
!         v_exb_local(1,lp,mp) = 0.0
!         v_exb_local(2,lp,mp) = 100.0  ! no convection

          ion_densities_N(:,:) = 0.0
          ion_velocities_N(:,:) = 0.0
          ion_temperature_N(:) = 0.0
          electron_temperature_N(:) = 0.0
          ion_densities_S(:,:) = 0.0
          ion_velocities_S(:,:) = 0.0
          ion_temperature_S(:) = 0.0
          electron_temperature_S(:) = 0.0
          ion_den_N_W(:,:) = 0.0
          ion_vel_N_W(:,:) = 0.0
          ion_temp_N_W(:) = 0.0
          electron_temp_N_W(:) = 0.0
          ion_den_N_E(:,:) = 0.0
          ion_vel_N_E(:,:) = 0.0
          ion_temp_N_E(:) = 0.0
          electron_temp_N_E(:) = 0.0
          ion_den_S_W(:,:) = 0.0
          ion_vel_S_W(:,:) = 0.0
          ion_temp_S_W(:) = 0.0
          electron_temp_S_W(:) = 0.0
          ion_den_S_E(:,:) = 0.0
          ion_vel_S_E(:,:) = 0.0
          ion_temp_S_E(:) = 0.0
          electron_temp_S_E(:) = 0.0

          phi_t0   = grid % magnetic_longitude(mp) - v_ExB_local(1,lp,mp)*time_step/(r*sin( colat_90km(lp) ) )

          coslam = cos( half_pi - grid % magnetic_colatitude(1,lp) )
          sinim  = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
          theta_t0 = colat_90km(lp) - v_ExB_local(2,lp,mp)*time_step/(r*sinim)

!          write(3000 + mpi_layer % rank_id, 77) mp , lp , time_step,  &
!                                                grid % magnetic_longitude(mp), &
!                                                colat_90km(lp), &
!                                                phi_t0, &
!                                                theta_t0, &
!                                                v_ExB(1,lp,mp), &
!                                                v_ExB(2,lp,mp)
           write(3000 + mpi_layer % rank_id, 77) mp , lp , &
                                                 grid % magnetic_longitude(mp), &
                                                 colat_90km(lp), &
                                                 grid % longitude(1,lp,mp), &
                                                 grid % colatitude(1,lp,mp), &
                                                 v_ExB_local(1,lp,mp), &
                                                 v_ExB_local(2,lp,mp)
 77        format(2i6,6f12.4)

          IF( phi_t0 < grid % magnetic_longitude(plasma % mp_low - plasma % mp_halo) .OR. &
              phi_t0 > grid % magnetic_longitude(plasma % mp_high +  plasma % mp_halo) )THEN

            write(3000 + mpi_layer % rank_id, *)  'phi_t0 outside of this rank range : [',&
                        grid % magnetic_longitude(plasma % mp_low - plasma % mp_halo),',', &
                        grid % magnetic_longitude(plasma % mp_high +  plasma % mp_halo),']'
            write(3000 + mpi_layer % rank_id, *)  'phi_t0 =', phi_t0
            write(3000 + mpi_layer % rank_id, *)  'mp =', mp
            write(3000 + mpi_layer % rank_id, *)  'longitude =', grid % magnetic_longitude(mp)
            write(3000 + mpi_layer % rank_id, *)  'Advective Distance =', ABS(phi_t0-grid % magnetic_longitude(mp))
!           STOP

          ENDIF

          ! We are trying to find the lp indices for the flux tubes nearest to
          ! theta_t0. Output is lp_t0(1:2)

          lp_min = grid % NLP

          DO lpx = 1, grid % NLP

            IF( theta_t0 <= colat_90km(lpx) )THEN
              lp_min = lpx
              EXIT
            ENDIF

          ENDDO

          lp_N = lp_min-1
          lp_S = lp_min

           mp_min = plasma % mp_high + plasma % mp_halo

          DO mpx = plasma % mp_low - plasma % mp_halo+1, plasma % mp_high + plasma % mp_halo

             IF( phi_t0 <= grid % magnetic_longitude(mpx) )THEN
               mp_min = mpx
               EXIT
             ENDIF

           ENDDO

           mp_W = mp_min-1  !this -1 because of +1 above
           mp_E = mp_min


        if(((mp.ne.mp_W).and.(mp.ne.mp_E)).or.(lp_N.eq.0)) then

        write(3000 + mpi_layer % rank_id, *)  'mp ', mp , mp_W , mp_E, '<<<<'
        write(3000 + mpi_layer % rank_id, *)  'lp ', lp , lp_N , lp_S, '<<<<'
        dont_convect(lp,mp) = 1

        else

        write(3000 + mpi_layer % rank_id, *)  'mp ', mp , mp_W , mp_E
        write(3000 + mpi_layer % rank_id, *)  'lp ', lp , lp_N , lp_S

        q_val(1:grid % flux_tube_max(lp)) = grid % q_factor(1:grid % flux_tube_max(lp),lp,mp)

        q_val_N_W(1:grid % flux_tube_max(lp_N)) = grid % q_factor(1:grid % flux_tube_max(lp_N),lp_N,mp_W)
        q_val_N_E(1:grid % flux_tube_max(lp_N)) = grid % q_factor(1:grid % flux_tube_max(lp_N),lp_N,mp_E)
        q_val_S_W(1:grid % flux_tube_max(lp_S)) = grid % q_factor(1:grid % flux_tube_max(lp_S),lp_S,mp_W)
        q_val_S_E(1:grid % flux_tube_max(lp_S)) = grid % q_factor(1:grid % flux_tube_max(lp_S),lp_S,mp_E)

        ion_den_NW(1:4,1:grid % flux_tube_max(lp_N)) =  plasma %  ion_densities_old(1:4,1:grid % flux_tube_max(lp_N),lp_N,mp_W)
        ion_vel_NW(1:4,1:grid % flux_tube_max(lp_N)) =  plasma %  ion_velocities_old(1:4,1:grid % flux_tube_max(lp_N),lp_N,mp_W)
        ion_temp_NW(1:grid % flux_tube_max(lp_N)) =  plasma %  ion_temperature_old(1:grid % flux_tube_max(lp_N),lp_N,mp_W)
        electron_temp_NW(1:grid % flux_tube_max(lp_N)) =  plasma %  electron_temperature_old(1:grid % flux_tube_max(lp_N),lp_N,mp_W)

        ion_den_NE(1:4,1:grid % flux_tube_max(lp_N)) =  plasma %  ion_densities_old(1:4,1:grid % flux_tube_max(lp_N),lp_N,mp_E)
        ion_vel_NE(1:4,1:grid % flux_tube_max(lp_N)) =  plasma %  ion_velocities_old(1:4,1:grid % flux_tube_max(lp_N),lp_N,mp_E)
        ion_temp_NE(1:grid % flux_tube_max(lp_N)) =  plasma %  ion_temperature_old(1:grid % flux_tube_max(lp_N),lp_N,mp_E)
        electron_temp_NE(1:grid % flux_tube_max(lp_N)) =  plasma %  electron_temperature_old(1:grid % flux_tube_max(lp_N),lp_N,mp_E)

        ion_den_SW(1:4,1:grid % flux_tube_max(lp_S)) =  plasma %  ion_densities_old(1:4,1:grid % flux_tube_max(lp_S),lp_S,mp_W)
        ion_vel_SW(1:4,1:grid % flux_tube_max(lp_S)) =  plasma %  ion_velocities_old(1:4,1:grid % flux_tube_max(lp_S),lp_S,mp_W)
        ion_temp_SW(1:grid % flux_tube_max(lp_S)) =  plasma %  ion_temperature_old(1:grid % flux_tube_max(lp_S),lp_S,mp_W)
        electron_temp_SW(1:grid % flux_tube_max(lp_S)) =  plasma %  electron_temperature_old(1:grid % flux_tube_max(lp_S),lp_S,mp_W)

        ion_den_SE(1:4,1:grid % flux_tube_max(lp_S)) =  plasma %  ion_densities_old(1:4,1:grid % flux_tube_max(lp_S),lp_S,mp_E)
        ion_vel_SE(1:4,1:grid % flux_tube_max(lp_S)) =  plasma %  ion_velocities_old(1:4,1:grid % flux_tube_max(lp_S),lp_S,mp_E)
        ion_temp_SE(1:grid % flux_tube_max(lp_S)) =  plasma %  ion_temperature_old(1:grid % flux_tube_max(lp_S),lp_S,mp_E)
        electron_temp_SE(1:grid % flux_tube_max(lp_S)) =  plasma %  electron_temperature_old(1:grid % flux_tube_max(lp_S),lp_S,mp_E)

        call interpolate_in_q(grid % flux_tube_max(lp), grid % flux_tube_max(lp_N), q_val, q_val_N_W, &
                             ion_den_NW, ion_vel_NW, ion_temp_NW, electron_temp_NW, ion_den_N_W, ion_vel_N_W, ion_temp_N_W, electron_temp_N_W)
        call interpolate_in_q(grid % flux_tube_max(lp), grid % flux_tube_max(lp_N), q_val, q_val_N_E, &
                             ion_den_NE, ion_vel_NE, ion_temp_NE, electron_temp_NE, ion_den_N_E, ion_vel_N_E, ion_temp_N_E, electron_temp_N_E)
        call interpolate_in_q(grid % flux_tube_max(lp), grid % flux_tube_max(lp_S), q_val, q_val_S_W, &
                             ion_den_SW, ion_vel_SW, ion_temp_SW, electron_temp_SW, ion_den_S_W, ion_vel_S_W, ion_temp_S_W, electron_temp_S_W)
        call interpolate_in_q(grid % flux_tube_max(lp), grid % flux_tube_max(lp_S), q_val, q_val_S_E, &
                             ion_den_SE, ion_vel_SE, ion_temp_SE, electron_temp_SE, ion_den_S_E, ion_vel_S_E, ion_temp_S_E, electron_temp_S_E)

        lon_fac = (phi_t0 - grid % magnetic_longitude(mp_W))/(grid % magnetic_longitude(mp_E) - grid % magnetic_longitude(mp_W))
        lat_fac = (theta_t0 - colat_90km(lp_N))/(colat_90km(lp_S) - colat_90km(lp_N))

        do ii = 1 , grid % flux_tube_max(lp)

        ion_densities_N(1:4,ii) = (lon_fac * (ion_den_N_E(1:4,ii) - ion_den_N_W(1:4,ii))) + ion_den_N_W(1:4,ii)
        ion_densities_S(1:4,ii) = (lon_fac * (ion_den_S_E(1:4,ii) - ion_den_S_W(1:4,ii))) + ion_den_S_W(1:4,ii)
        ion_densities_t0(1:4,ii,lp,mp) = (lat_fac * (ion_densities_S(1:4,ii) - ion_densities_N(1:4,ii))) + ion_densities_N(1:4,ii)

        ion_velocities_N(1:4,ii) = (lon_fac * (ion_vel_N_E(1:4,ii) - ion_vel_N_W(1:4,ii))) + ion_vel_N_W(1:4,ii)
        ion_velocities_S(1:4,ii) = (lon_fac * (ion_vel_S_E(1:4,ii) - ion_vel_S_W(1:4,ii))) + ion_vel_S_W(1:4,ii)
        ion_velocities_t0(1:4,ii,lp,mp) = (lat_fac * (ion_velocities_S(1:4,ii) - ion_velocities_N(1:4,ii))) + ion_velocities_N(1:4,ii)

        ion_temperature_N(ii) = (lon_fac * (ion_temp_N_E(ii) - ion_temp_N_W(ii))) + ion_temp_N_W(ii)
        ion_temperature_S(ii) = (lon_fac * (ion_temp_S_E(ii) - ion_temp_S_W(ii))) + ion_temp_S_W(ii)
        ion_temperature_t0(ii,lp,mp) = (lat_fac * (ion_temperature_S(ii) - ion_temperature_N(ii))) + ion_temperature_N(ii)

        electron_temperature_N(ii) = (lon_fac * (electron_temp_N_E(ii) - electron_temp_N_W(ii))) + electron_temp_N_W(ii)
        electron_temperature_S(ii) = (lon_fac * (electron_temp_S_E(ii) - electron_temp_S_W(ii))) + electron_temp_S_W(ii)
        electron_temperature_t0(ii,lp,mp) = (lat_fac * (electron_temperature_S(ii) - electron_temperature_N(ii))) + electron_temperature_N(ii)

        enddo ! ii

!       if ((mp.eq.21).and.(lp.eq.6)) then
        if (lp.eq.6) then

        write(7050 + mpi_layer % rank_id, *) plasma % mp_low , plasma % mp_halo+1, plasma % mp_high , plasma % mp_halo
        write(7050 + mpi_layer % rank_id, *) plasma % mp_low - plasma % mp_halo+1, plasma % mp_high + plasma % mp_halo
        write(7050 + mpi_layer % rank_id, 12)  mp , mp_W , mp_E, phi_t0 * r2d, &
                                               grid % magnetic_longitude(mp_W) * r2d, &
                                               grid % magnetic_longitude(mp_E) * r2d
 12     format('mp ',3i4,3f10.4)
        write(7050 + mpi_layer % rank_id, 13)  v_ExB_local(1,lp,mp)
 13     format('mp vel ',f10.2)
        write(7050 + mpi_layer % rank_id, 14)  lp , lp_N , lp_S, theta_t0  * r2d, &
                                               colat_90km(lp_N) * r2d, &
                                               colat_90km(lp_S) * r2d
 14     format('lp ',3i4,3f10.4)
        write(7050 + mpi_layer % rank_id, 15)  v_ExB_local(2,lp,mp)
 15     format('lp vel ',f10.2)
        write(7050 + mpi_layer % rank_id, *) 'fac_lon , fac_lat ', lon_fac, lat_fac

        do ii = 40 , 50
        write(7050 + mpi_layer % rank_id, 39) ii, q_val(ii),q_val_N_W(ii),q_val_N_E(ii),q_val_S_W(ii),q_val_S_E(ii)
 39     format(i5,5e16.6)
        enddo

        write(7050 + mpi_layer % rank_id, *) '***************NW NE SW SE*(pre q interpolation)******************'
        do ii = 40 , 50
        write(7050 + mpi_layer % rank_id, 33) ii, ion_den_NW(1,ii), ion_den_NE(1,ii), ion_den_SW(1,ii),ion_den_SE(1,ii)
        enddo

        write(7050 + mpi_layer % rank_id, *) '***************N_W N_E S_W S_E*(post q interpolation)******************'
        do ii = 40 , 50
            write(7050 + mpi_layer % rank_id, 33) ii , ion_den_N_W(1,ii), ion_den_N_E(1,ii), ion_den_S_W(1,ii),ion_den_S_E(1,ii)
 33         format(i5,4e16.6)
        enddo
        write(7050 + mpi_layer % rank_id, *) '****************** N S t0 ******************************************'
        do ii = 40 , 50
            write(7050 + mpi_layer % rank_id, 38) ii , ion_densities_N(1,ii),ion_densities_S(1,ii),ion_densities_t0(1,ii,lp,mp)
 38         format(i5,3e16.6)
         enddo

        endif ! mp 21 lp 6

        endif ! outside mp, lp box

        ENDDO
      ENDDO

! finally copy across.....

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, high_lat_convection_max_lp

        if (dont_convect(lp,mp).eq.0) then

        do ii = 1 , grid % flux_tube_max(lp)

        if( grid % altitude(ii,lp) > transport_min_altitude ) then

        plasma % ion_densities(1:4,ii,lp,mp) = ion_densities_t0(1:4,ii,lp,mp)
        plasma % ion_velocities(1:4,ii,lp,mp) = ion_velocities_t0(1:4,ii,lp,mp)
        plasma % ion_temperature(ii,lp,mp) = ion_temperature_t0(ii,lp,mp)
        plasma % electron_temperature(ii,lp,mp) = electron_temperature_t0(ii,lp,mp)

        endif ! 150

        enddo ! ii

        endif ! dont_convect(lp,mp)

        ENDDO
      ENDDO

  END SUBROUTINE High_Latitude_Flux_Tube_Transport


  SUBROUTINE Cross_Flux_Tube_Transport2( plasma, grid, v_ExB, time_step, mpi_layer, time_tracker, forcing)
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec), INTENT(in)             :: time_step
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    ! Local
    INTEGER    :: mp, lp, i
    INTEGER    :: midpoint
    INTEGER    :: day_of_the_year
    INTEGER    :: ions
    REAL(prec), PARAMETER :: rad_to_deg = 57.295779513
    REAL(prec) :: ut_hours
    REAL(prec) :: local_time_hours
    REAL(prec) :: geo_lon_degrees
    REAL(prec) :: v_upwards_at_apex
    REAL(prec) :: f107d
    REAL(prec) :: apex_height_km
    REAL(prec) :: xmlon_300
    REAL(prec) :: xmlat_300
    REAL(prec) :: magnitude_e2_at_300
    REAL(prec) :: magnitude_e2_at_apex
    REAL(prec) :: e2_factor_300km_to_apex
    REAL(prec) :: re_apex(1:grid % NLP)
    REAL(prec) :: q_val(1:grid % nFluxTube,1:grid % NLP)
    REAL(prec) :: ion_den(1:9,1:grid % nFluxTube,1:grid % NLP)
    REAL(prec) :: ion_vel(1:9,1:grid % nFluxTube,1:grid % NLP)
    REAL(prec) :: ion_temp(1:grid % nFluxTube,1:grid % NLP)
    REAL(prec) :: electron_temp(1:grid % nFluxTube,1:grid % NLP)
    REAL(prec) :: ion_den_final(1:9,1:grid % nFluxTube)
    REAL(prec) :: ion_vel_final(1:9,1:grid % nFluxTube)
    REAL(prec) :: ion_temp_final(1:grid % nFluxTube)
    REAL(prec) :: electron_temp_final(1:grid % nFluxTube)

      ut_hours = time_tracker % utime / 3600.0
      day_of_the_year = time_tracker % day_of_year
      F107D = forcing % f107( forcing % current_index )
!     write(8000 + mpi_layer % rank_id, *) day_of_the_year, F107D, ut_hours

      DO mp = plasma % mp_low, plasma % mp_high

        do lp = 1 , grid % NLP
        do i = 1 , grid % flux_tube_max(lp)
        re_apex(lp) = earth_radius + grid % altitude(grid % flux_tube_midpoint(lp),lp)
        q_val(i,lp) = grid % q_factor(i,lp,mp)
        ion_temp(i,lp) = plasma % ion_temperature(i,lp,mp)
        electron_temp(i,lp) = plasma % electron_temperature(i,lp,mp)
        do ions = 1 , 9
        ion_den(ions,i,lp) = plasma % ion_densities(ions,i,lp,mp)
        ion_vel(ions,i,lp) = plasma % ion_velocities(ions,i,lp,mp)
        enddo
        enddo
        enddo

        DO lp = 20, perp_transport_max_lp

           midpoint = grid % flux_tube_midpoint(lp)
           geo_lon_degrees = grid % longitude(midpoint,lp,mp) * rad_to_deg
           local_time_hours = (geo_lon_degrees / 15.0) + ut_hours
           if(local_time_hours > 24.) local_time_hours = local_time_hours - 24.0

! GHGM - now use the mixed Fejer/Richmond version....
!          call fejer_exb_model(local_time_hours,geo_lon_degrees, &
!                               day_of_the_year,F107d,v_upwards_at_apex)

           apex_height_km = grid % altitude(grid % flux_tube_midpoint(lp),lp) / 1000.

           xmlon_300 = grid % magnetic_longitude(mp) * rad_to_deg
! GHGM - this 10 needs sorting - index for 300km...
           xmlat_300 = 90. - (grid % magnetic_colatitude(10,lp) * rad_to_deg)
              !g

! GHGM - this 10 needs sorting - index for 300km...
           magnitude_e2_at_300 = &
           sqrt((grid % apex_e_vectors(1,2,10,lp,mp)*grid % apex_e_vectors(1,2,10,lp,mp)) + &
           (grid % apex_e_vectors(2,2,10,lp,mp)*grid % apex_e_vectors(2,2,10,lp,mp)) + &
           (grid % apex_e_vectors(3,2,10,lp,mp)*grid % apex_e_vectors(3,2,10,lp,mp)))
              !g
           magnitude_e2_at_apex = &
           sqrt((grid % apex_e_vectors(1,2,midpoint,lp,mp)*grid % apex_e_vectors(1,2,midpoint,lp,mp)) + &
           (grid % apex_e_vectors(2,2,midpoint,lp,mp)*grid % apex_e_vectors(2,2,midpoint,lp,mp)) + &
           (grid % apex_e_vectors(3,2,midpoint,lp,mp)*grid % apex_e_vectors(3,2,midpoint,lp,mp)))

           e2_factor_300km_to_apex = magnitude_e2_at_apex / magnitude_e2_at_300

           call get_empirical_vertical_exb(apex_height_km,xmlat_300,xmlon_300, &
                                           local_time_hours,geo_lon_degrees, &
                                           ut_hours,day_of_the_year,f107d,e2_factor_300km_to_apex, &
                                           v_upwards_at_apex)

!          write(8000 + mpi_layer % rank_id, 77) mp , lp , v_upwards_at_apex

!                midpoint, &
!                grid % altitude(midpoint,lp) / 1000., &
!                grid % longitude(midpoint,lp,mp) * rad_to_deg
 77        format(2i6,2f20.3)
!          write(8000 + mpi_layer % rank_id,*) time_tracker % utime, &
!                                              time_tracker % day_of_year
           CALL TUBES_LINEAR_INTERPOLATE(lp,grid % flux_tube_max, &
                                         re_apex,v_upwards_at_apex,time_step, &
                                         q_val, &
                                         ion_den, &
                                         ion_vel, &
                                         ion_temp, &
                                         electron_temp, &
                                         ion_den_final, &
                                         ion_vel_final, &
                                         ion_temp_final, &
                                         electron_temp_final)

        do i = 1 , grid % flux_tube_max(lp)
        ion_temp(i,lp) = ion_temp_final(i)
        electron_temp(i,lp) = electron_temp_final(i)
        do ions = 1 , 9
        ion_den(ions,i,lp) = ion_den_final(ions,i)
        ion_vel(ions,i,lp) = ion_vel_final(ions,i)
        enddo
        enddo


        ENDDO

        DO lp = 45, perp_transport_max_lp
        do i = 1 , grid % flux_tube_max(lp)
        plasma % ion_temperature(i,lp,mp) = ion_temp(i,lp)
        plasma % electron_temperature(i,lp,mp) = electron_temp(i,lp)
        do ions = 1 , 9
          plasma % ion_densities(ions,i,lp,mp) = ion_den(ions,i,lp)
          plasma % ion_velocities(ions,i,lp,mp) = ion_vel(ions,i,lp)
        enddo
        enddo
        enddo

      ENDDO

  END SUBROUTINE Cross_Flux_Tube_Transport2


SUBROUTINE get_empirical_vertical_exb(apex_height_km,xmlat_300,xmlon_300,lt,geolon_apex, &
  ut,iday,f107,e2_factor_300km_to_apex, &
  v_upwards_at_apex)

  implicit none
  REAL(kind=8), INTENT(in)  :: apex_height_km
  REAL(kind=8), INTENT(in)  :: xmlat_300
  REAL(kind=8), INTENT(in)  :: xmlon_300
  REAL(kind=8), INTENT(in)  :: lt
  REAL(kind=8), INTENT(in)  :: geolon_apex
  REAL(kind=8), INTENT(in)  :: ut
  REAL(kind=8), INTENT(in)  :: f107
  REAL(kind=8), INTENT(in)  :: e2_factor_300km_to_apex
  REAL(kind=8) v_upwards_at_apex_fejer
  REAL(kind=8) v_upwards_at_apex_richmond
  REAL(kind=8) ve
  REAL(kind=8) pot
  REAL(kind=8) dayno
  REAL(kind=8) factor_less_than_300
  REAL(kind=8) apex_factor
  REAL(kind=8) v_outwards_at_300km
  INTEGER :: iday
  INTEGER :: iseasav
  INTEGER :: iutav
  REAL(kind=8), INTENT(out)  :: v_upwards_at_apex
!g



!g  This combines the low-latitude electric field models
!g  of Richmond and Fejer.
!g  The Fejer model is used for flux-tubes with apex heights up to 1000km.
!g  The Richmond model is used for flux-tubes with apex heights greater than
!g  2000km.  Between 1000km and 2000km a linear interpolation between the two
!g  is used.
!g
  if(apex_height_km < 2000.) then
  !g
  !g  use the Fejer model....
  !g
      call fejer_exb_model(LT,Geolon_apex,IDAY,F107,v_upwards_at_apex_fejer)
  !g
  !g  if the apex height is less than 300km then tail the upwards drift off
  !g  in height between a full value at 300km and 0.0_prec at 100km....
  !g
      if(apex_height_km < 300.) then
          factor_less_than_300 = (apex_height_km - 100. ) / 200.
          v_upwards_at_apex_fejer = factor_less_than_300 * v_upwards_at_apex_fejer
          if(apex_height_km < 100.) v_upwards_at_apex_fejer = 0.0
      endif
  !g
      v_upwards_at_apex = v_upwards_at_apex_fejer
  !g
  endif

  if(apex_height_km > 1000.) then
  !g
  !g  use the Richmond model....
  !g
      Iseasav = 3
      Iutav = 0
      dayno = float(iday)
      call Richmond_lowlat_Efield_model(XMLat_300,XMLon_300,dayno,UT,ISEasav,IUTav, &
      POT,v_outwards_at_300km,VE)
  !g
      v_upwards_at_apex_richmond = v_outwards_at_300km * e2_factor_300km_to_apex
  !g
      v_upwards_at_apex = v_upwards_at_apex_richmond
  !g
  endif

  if(apex_height_km > 1000. .AND. apex_height_km < 2000.) then
  !g
  !g  interpolate between the 2 models according to the apex height....
  !g
      apex_factor = (apex_height_km - 1000.) / 1000.
      v_upwards_at_apex = apex_factor * (v_upwards_at_apex_richmond - v_upwards_at_apex_fejer) &
      + v_upwards_at_apex_fejer
  endif
!g
  return



end SUBROUTINE get_empirical_vertical_exb


SUBROUTINE fejer_exb_model(SLT,GL,IDAY,F107,VPE)

  IMPLICIT NONE

    REAL(kind=8), INTENT(IN)  :: SLT, GL, F107
    INTEGER,      INTENT(IN)  :: IDAY

    REAL(kind=8), INTENT(OUT) :: VPE

!*******************************************************************************

! ROUTINE TO DETERMINE VERTICAL EXB DRIFT
! (SCHERLIESS, L., AND B.G. FEJER,
! J. GEOPHYS. RES., 104, 6829-6842, 1999)

!*******************************************************************************

! SLT:   SOLAR LOCAL TIME
! GL:    GEOGRAPHIC LONGITUDE (+ EAST)
! IDAY:  DAY OF YEAR
! F107:  F10.7
! VPE:   EQUATORIAL VERTICAL DRIFT

!*******************************************************************************

  INTEGER :: i, kk, ind, k , j
  REAL(kind=8) :: BSPL4, funct, coeff
  REAL(kind=8) :: BSPL4T , BSPL4L

  DIMENSION COEFF(624),FUNCT(6)



  DATA (COEFF(I),I=1,60)/ &
  -10.80592, -9.63722,-11.52666, -0.05716, -0.06288,  0.03564, &
  -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138, &
  2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171, &
  -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656, &
  -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419, &
  -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984, &
  -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062, &
  -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297, &
  1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023, &
  5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166/
  DATA (COEFF(I),I=61,120)/ &
  3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456, &
  7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008, &
  -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507, &
  -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477, &
  -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528, &
  -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075, &
  14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723, &
  12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484, &
  18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249, &
  4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394/
  DATA (COEFF(I),I=121,180)/ &
  14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739, &
  7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626, &
  7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039, &
  10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177, &
  21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715, &
  19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434, &
  26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631, &
  21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903, &
  28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622, &
  22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206/
  DATA (COEFF(I),I=181,240)/ &
  31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322, &
  46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982, &
  13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970, &
  21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589, &
  16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570, &
  18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991, &
  10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775, &
  12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851, &
  -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433, &
  -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385/
  DATA (COEFF(I),I=241,300)/ &
  2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932, &
  3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084, &
  17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516, &
  0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331, &
  15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934, &
  4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308, &
  9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124, &
  13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465, &
  5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222, &
  9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418/
  DATA (COEFF(I),I=301,360)/ &
  9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919, &
  13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801, &
  10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208, &
  10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751, &
  6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705, &
  5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151, &
  29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325, &
  17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398, &
  8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374, &
  -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187/
  DATA (COEFF(I),I=361,420)/ &
  8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156, &
  14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061, &
  14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129, &
  6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993, &
  7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631, &
  -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188, &
  -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214, &
  -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655, &
  2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650, &
  5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376/
  DATA (COEFF(I),I=421,480)/ &
  13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223, &
  -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585, &
  -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463, &
  -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679, &
  -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355, &
  -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593, &
  -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839, &
  -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226, &
  -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004, &
  -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699/
  DATA (COEFF(I),I=481,540)/ &
  -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796, &
  -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959, &
  -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792, &
  -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535, &
  -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143, &
  -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688, &
  -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422, &
  -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978, &
  -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144, &
  -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133/
  DATA (COEFF(I),I=541,600)/ &
  -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949, &
  -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583, &
  -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213, &
  -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053, &
  -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085, &
  -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726, &
  -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529, &
  -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544, &
  -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807, &
  -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214/
  DATA (COEFF(I),I=601,624)/ &
  -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652, &
  -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235, &
  -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214, &
  -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/

  CALL fejer_G(IDAY,F107,FUNCT,GL)


  VPE=0.
  DO I=1,13
      DO J=1,8,1
          KK=8*(I-1)+J
          DO K=1,6
              IND=6*(KK-1)+K

              call FUNCTION_BSPL4T(I,SLT,BSPL4T)
              call FUNCTION_BSPL4L(J,GL,BSPL4L)

              BSPL4=BSPL4T*BSPL4L

              VPE=VPE+BSPL4*FUNCT(K)*COEFF(IND)
          ENDDO
      ENDDO
  ENDDO

  RETURN




end SUBROUTINE fejer_exb_model



SUBROUTINE fejer_G(IDAY,F107,FUNCT,T)

  IMPLICIT NONE

  REAL(kind=8) :: f107, T
  INTEGER      :: iday
  REAL(kind=8) :: FUNCT, flux
  REAL(kind=8) :: t1, sigma, gauss, cflux, a
  INTEGER      :: i, kk

  DIMENSION FUNCT(6)

  IF(F107 < 75.) THEN
      FLUX=75.
  ELSEIF(F107 > 230.) THEN
      FLUX=230.
  ELSE
      FLUX=F107
  ENDIF
  CFLUX=FLUX

  A=0.
  IF(IDAY >= 120 .AND. IDAY <= 240) THEN
      A=170.
      SIGMA=60.
  ENDIF
  IF(IDAY >= 60 .OR. IDAY >= 300.) THEN
      A=170.
      SIGMA=40.
  ENDIF
  IF(FLUX < 95. .AND. A > 0.) THEN
      GAUSS=EXP(-0.5*((T-A)**2)/SIGMA**2)
      CFLUX=GAUSS*95.+(1.-GAUSS)*FLUX
  ENDIF

  DO 11 I=1,6
      FUNCT(I)=0.
  11 ENDDO

  IF(IDAY >= 135 .AND. IDAY <= 230) THEN
      FUNCT(1)=1.
  ENDIF
  IF(IDAY <= 45 .OR. IDAY >= 320) THEN
      FUNCT(2)=1.
  ENDIF
  IF(IDAY >= 75 .AND. IDAY <= 105) THEN
      FUNCT(3)=1.
  ENDIF
  IF(IDAY >= 260 .AND. IDAY <= 290) THEN
      FUNCT(3)=1.
  ENDIF

  IF(IDAY >= 45 .AND. IDAY <= 75) THEN    ! W-E
      FUNCT(2)=1.-(IDAY-45)/30.
      FUNCT(3)=1.-FUNCT(2)
  ENDIF
  IF(IDAY >= 105 .AND. IDAY <= 135) THEN  ! E-S
      FUNCT(3)=1.-(IDAY-105)/30.
      FUNCT(1)=1.-FUNCT(3)
  ENDIF
  IF(IDAY >= 230 .AND. IDAY <= 260) THEN  ! S-E
      FUNCT(1)=1.-(IDAY-230)/30.
      FUNCT(3)=1.-FUNCT(1)
  ENDIF
  IF(IDAY >= 290 .AND. IDAY <= 320) THEN  ! E-W
      FUNCT(3)=1.-(IDAY-290)/30.
      FUNCT(2)=1.-FUNCT(3)
  ENDIF

  FUNCT(4)=(CFLUX-140.)*FUNCT(1)
  FUNCT(5)=(CFLUX-140.)*FUNCT(2)
  FUNCT(6)=(CFLUX-140.)*FUNCT(3)

  RETURN




end SUBROUTINE fejer_G

SUBROUTINE FUNCTION_BSPL4T(I,T1,BSPL4T)

  IMPLICIT NONE

  integer, INTENT(in) :: i
  REAL(kind=8), INTENT(in) :: T1
  REAL(kind=8), INTENT(out) :: bspl4t

  INTEGER :: j,k
  REAL(kind=8) :: T
  REAL(kind=8) :: tt,b

  DIMENSION TT(0:39),B(20,20)

  DATA TT/ 0.00, 2.75, 4.75, 5.50, 6.25, &
  7.25,10.00,14.00,17.25,18.00, &
  18.75,19.75,21.00,24.00,26.75, &
  28.75,29.50,30.25,31.25,34.00, &
  38.00,41.25,42.00,42.75,43.75, &
  45.00,48.00,50.75,52.75,53.50, &
  54.25,55.25,58.00,62.00,65.25, &
  66.00,66.75,67.75,69.00,72.00/


  T=T1
  IF(I >= 0 .AND. T < TT(I)) THEN
      T=T+24.
  ENDIF
  DO J=I,I+4-1
      IF(T >= TT(J) .AND. T < TT(J+1)) THEN
          B(J,1)=1.
      ELSE
          B(J,1)=0.
      ENDIF
  ENDDO
  DO J=2,4
      DO K=I,I+4-J
          B(K,J)=(T-TT(K))/(TT(K+J-1)-TT(K))*B(K,J-1)
          B(K,J)=B(K,J)+(TT(K+J)-T)/(TT(K+J)-TT(K+1))*B(K+1,J-1)
      ENDDO
  ENDDO

  BSPL4T=B(I,4)

  RETURN




end SUBROUTINE FUNCTION_BSPL4T



SUBROUTINE FUNCTION_BSPL4L(I,T1,bspl4l)

  IMPLICIT NONE

  integer, INTENT(in) :: i
  REAL(kind=8), INTENT(in) :: T1
  REAL(kind=8), INTENT(out) :: bspl4l

  INTEGER :: j,k
  REAL(kind=8) :: TL, T, B

  DIMENSION TL(0:24),B(20,20)

  DATA TL/  0, 10,100,190,200,250,280,310, &
  360,370,460,550,560,610,640,670, &
  720,730,820,910,920,970,1000,1030,1080/

  T=T1
  IF(I >= 0 .AND. T < TL(I)) THEN
      T=T+360.
  ENDIF
  DO J=I,I+4-1
      IF(T >= TL(J) .AND. T < TL(J+1)) THEN
          B(J,1)=1.
      ELSE
          B(J,1)=0.
      ENDIF
  ENDDO

  DO J=2,4
      DO K=I,I+4-J
          B(K,J)=(T-TL(K))/(TL(K+J-1)-TL(K))*B(K,J-1)
          B(K,J)=B(K,J)+(TL(K+J)-T)/(TL(K+J)-TL(K+1))*B(K+1,J-1)
      ENDDO
  ENDDO

  BSPL4L=B(I,4)

  RETURN




end SUBROUTINE FUNCTION_BSPL4L


SUBROUTINE Richmond_lowlat_Efield_model(XMLat,XMLon,DAYno,UT,ISEasav,IUTav,POT,VU,VE)
  IMPLICIT NONE

  REAL(kind=8) :: &
  a , ang , cl , cml , ct , cts , dang , DAYno , daynop , fs , &
  ft , fut , hrang , p , pa , pb , POT , q , rad , rb
  REAL(kind=8) :: &
  rbt , rr , rs , rsm , rsrs , sl , sml , sq2 , st , tl , tla , &
  tlp , tua , tup , UT , va , vb , VE , vp , VU
  REAL(kind=8) :: x , xm , XMLat , xmlatp , XMLon , xnms , xns , xz , z , zp
  INTEGER :: i , icpt , imax , ISEasav , IUTav , j , jf , k , kf , kp , &
  l , lend , lf , lp , lst , m , mf , mm
  INTEGER :: mmf , mp , mpf , mpp , n , nf , np
!g
!g  Equatorial Electric Field model.......
!g
! gives quiet-day ionospheric electrostatic pseudo-potential and e x b
! drifts at 300 km for solar minimum conditions.  see richmond et al. (jgr,
! 1980, p. 4658) for definitions of magnetic coordinates, pseudo-potential,
! and drift components.
!***********************************************************************
! 8/2/2 routine has been modified by deleting shortcuts that depend
! on the invalid assumption that variables are saved between
! successive calls to it.  to reactivate the shortcuts it will be
! necessary to use the save statement to retain all needed variables.
! input  PARAMETERs -
! xmlat, xmlon are magnetic latitude and east longitude in degrees.
! dayno is day number of the year from 1. to 365.24, with 1. being jan. 1.
! ut is universal time in hours.
! iseasav is 0 if no seasonal averaging is desired.
! is 1 for average over nov. - feb.
! is 2 for average over may -aug.
! is 3 for average over mar., apr., sept., oct.
! is 4 for average over entire year.
! if iseasav.ne.0, dayno is ignored.
! iutav is 0 if no ut averaging is desired.
! is 1 for average over all ut at the fixed local time given by
! ut + (xmlon - 69.)/15.

! output  PARAMETERs -
! pot is the electrostatic pseudo-potential in volts.
! vu is the drift velocity component perpendicular to the geomagnetic field
! in the upward/poleward direction in the magnetic meridian plane, in m/s.
! ve is the drift velocity component in the magnetic eastward direction, in
! m/s.
! the output values are geophysically meaningful only for latitudes between
! about -65 and +65 degrees.  if iseasav or iutav is out of range pot, vu,
! and ve are set to -1/0.
  DIMENSION kf(128) , lf(128) , mf(128) , nf(128) , jf(128) , q(5) , &
  rs(16) , fut(5) , rr(16) , rsrs(16) , p(16) , vp(16) , &
  pa(5,9) , va(5,9) , fs(3) , ft(3,3) , sml(4) , cml(4) , &
  a(128) , pb(9) , vb(9)
  DATA a/ - 70. , -183. , 31. , -112. , 19. , -39. , -2. , 2. , &
  -33. , 2. , 2. , -111. , 46. , -4. , -5. , 7. , 9. , -17. , &
  2. , 9. , -10. , 2. , -9. , 22. , 145. , -57. , -42. , -6. , &
  6. , -5. , -2. , 20. , 16. , 16. , -77. , -18. , 13. , -8. , &
  16. , -52. , -10. , 7. , 2. , 11. , -28. , 2. , -85. , -82. , &
  3. , -281. , -71. , -25. , -57. , -50. , 21. , -10. , 10. , &
  -81. , 24. , 7. , 5. , 30. , 32. , 5. , -5. , 11. , -31. , &
  8. , 10. , 20. , -15. , -42. , 32. , 7. , -19. , 7. , 34. , &
  -11. , -15. , 26. , 21. , 1. , 22. , 12. , -2. , 275. , &
  777. , -318. , -320. , -208. , 47. , 429. , -523. , 8. , &
  -35. , -224. , -450. , -66. , -7. , -8. , -231. , 55. , 6. , &
  -28. , -51. , -81. , 48. , 9. , 2. , -10. , 54. , 16. , &
  112. , 69. , -33. , 120. , -47. , 5. , -19. , -17. , -23. , &
  -40. , -22. , -21. , -7. , -30. , 15. , 3./
!***********************************************************************
! 8/2/2 deactivate shortcut, since it won't necessarily work right if
! variables are not retained with a save statement.
!***********************************************************************
! (through statement 510) set up constant  PARAMETERs in first call to
! routine.
! set up values of xmlatp, tup, tlp, and daynop which are not equal to xmlat,
! ut, magnetic local time, and dayno, respectively.
  xmlatp = -361.
  IF ( xmlatp == XMLat ) xmlatp = 0.
  tup = -25.
  IF ( tup == UT ) tup = 0.
  tl = UT + (XMLon-69.)/15.
  tlp = -25.
  IF ( tlp == tl ) tlp = 0.
  daynop = -366.
  IF ( daynop == DAYno ) daynop = 0.
! (through statement 100) select only those terms in series for which
! coefficients a are defined to be non-0.0_prec.
  i = 0
  DO 100 kp = 1 , 3
      lst = 4 - kp
      lend = 2 + kp
      DO 50 lp = lst , lend
          l = lp - 3
          IF ( MOD(kp+lp,2) == 0 ) THEN
              DO 10 np = 2 , 8
                  n = np - 1
                  IF ( kp == 3 .OR. n <= 6 ) THEN
                      IF ( IABS(l) /= 2 .OR. n <= 5 ) THEN
                          DO 2 mp = 1 , 9
                              m = mp - 5
                              IF ( IABS(m) <= n ) THEN
                                  IF ( IABS(m) <= 3 .OR. IABS(l) /= 2 ) THEN
                                      IF ( MOD(n-m,2) == 0 ) THEN
                                          i = i + 1
                                          kf(i) = kp
                                          lf(i) = lp
                                          nf(i) = np
                                          mf(i) = mp
                                      ENDIF
                                  ENDIF
                              ENDIF
                          2 ENDDO
                      ENDIF
                  ENDIF
              10 ENDDO
          ENDIF
      50 ENDDO
  100 ENDDO
  imax = i
  ft(1,1) = .75*SQRT(6.E0)/3.1415926535898
  ft(1,2) = 2.E0*ft(1,1)
  ft(1,3) = 1.E0
  ft(2,1) = ft(1,1)
  ft(2,2) = -ft(1,2)
  ft(2,3) = 1.E0
  ft(3,1) = ft(2,2)
  ft(3,2) = 0.
  ft(3,3) = 1.E0
  hrang = 3.1415926535898/12.
  dang = 3.1415926535898/182.62
  sq2 = SQRT(2.E0)
  rad = 180./3.1415926535898
  rb = -6.671E6*5.2E-5
! rb is -(earth radius + 3.e5 m) times dipole magnetic field at pole at 300 km.
  DO 200 i = 1 , imax
      mm = IABS(mf(i)-5)
  ! jf gives appropriate index of legendre polynomials as ordered between
  ! statements 530 and 595.
      jf(i) = (2*(nf(i)+7*mm)-(mm-1)**2+4)/4
  200 ENDDO
! (through statement 500) compute rr (defined as r(n,m)*r(n-1,m)) and rsrs
! (defined as r(n-1,m)**2 + r(n-2,m)**2) needed for legendre polynomial
! generating recursion relations, where r(n,m) is defined as
! sqrt(n**2 - m**2)/sqrt(4*n**2 - 1).  ordering is same as for p(n,m).
  j = 0
  DO 300 mp = 1 , 5
      m = mp - 1
      xm = m
      IF ( m /= 0 ) q(mp) = SQRT((2.*xm+1.)/(2.*xm))
      DO 250 np = mp , 8 , 2
          n = np - 1
          xns = n*n
          xnms = (n-1)**2
          j = j + 1
          rs(j) = (xns-xm*xm)/(4.*xns-1.)
          rsm = (xnms-xm*xm)/(4.*xnms-1.)
          rr(j) = SQRT(rs(j)*rsm)
          IF ( np /= mp ) rsrs(j) = rsm + rs(j-1)
      250 ENDDO
  300 ENDDO
  IF ( IUTav >= 0 .AND. IUTav <= 1 ) THEN
      IF ( IABS(ISEasav-2) <= 2 ) THEN
          icpt = 1
      ! (through statement 530) fs(1), fs(2), fs(3) are factors for amplitude of
      ! semiannual, annual, yearly average components, respectively.
      ! if iseasav = 0, compute fs for given day of the year.
          IF ( ISEasav /= 0 ) THEN
          ! if iseasav = 4, use only yearly average component.
              IF ( ISEasav == 4 ) THEN
                  fs(1) = 0.
                  fs(2) = 0.
                  fs(3) = 1.
              ELSE
              ! if iseasav = 1 - 3, compute fs for appropriate seasonal average.
                  DO 305 k = 1 , 3
                      fs(k) = ft(ISEasav,k)
                  305 ENDDO
              ENDIF
              GOTO 320
          ENDIF
      !***********************************************************************
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (dayno.eq.daynop) go to 530
      !***********************************************************************
          daynop = DAYno
          icpt = 1
          ang = (DAYno+9.)*dang
          fs(1) = sq2*COS(2.*ang)
          fs(2) = sq2*COS(ang)
          fs(3) = 1.
      ! if magnetic latitude is same as in previous call to this routine, skip to
      ! statement 596.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (xmlat.eq.xmlatp) go to 596
      !***********************************************************************
          320 icpt = 1
          xmlatp = XMLat
          ct = SIN(XMLat/rad)
          cts = ct*ct
          st = SQRT(1.-cts)
          rbt = rb*SQRT(.25+.75*cts)
      ! (through statement 595) calculate legendre polynomials p(n,m) as well as vp,
      ! defined as (dp(n,m)/d(colatitude))/(rb*ct).
          j = 0
          DO 340 mp = 1 , 5
          ! for mp=1, p=p(n,m).  for mp.gt.1, p=p(n,m)/st.  this difference is so that
          ! program never divides by st.
              j = j + 1
              xm = mp - 1
              mpp = mp + 2
              IF ( mp > 1 ) THEN
                  p(j) = q(mp)*x
                  x = p(j)*st
                  vp(j) = xm*p(j)/rb
                  xz = 2.*st*st/rb
              ELSE
                  x = 1.
                  p(1) = 1.E0
                  xz = 2.*st/rb
                  vp(1) = 0.
              ENDIF
              DO 330 np = mpp , 8 , 2
                  j = j + 1
                  z = 0.
                  zp = 0.
                  IF ( np /= mpp ) THEN
                      z = rr(j-1)*p(j-2)
                      zp = rr(j-1)*vp(j-2)
                  ENDIF
                  p(j) = ((cts-rsrs(j))*p(j-1)-z)/rr(j)
                  vp(j) = ((cts-rsrs(j))*vp(j-1)-zp-xz*p(j-1))/rr(j)
              330 ENDDO
          340 ENDDO
      ! (through statement 600) calculate arrays of fourier coefficients, with lp
      ! indicating harmonic of ut and mpf indicating harmonic of magnetic local
      ! time.  if arrays are same as in previous call to this routine, skip to
      ! statement 601.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (icpt.eq.0) go to 601
      !***********************************************************************
          DO 360 mpf = 1 , 9
              DO 350 lp = 1 , 5
                  pa(lp,mpf) = 0.
                  va(lp,mpf) = 0.
              350 ENDDO
          360 ENDDO
          DO 380 i = 1 , imax
              lp = lf(i)
              mpf = mf(i)
              j = jf(i)
              x = a(i)*fs(kf(i))
              pa(lp,mpf) = pa(lp,mpf) + x*p(j)
              va(lp,mpf) = va(lp,mpf) + x*vp(j)
          380 ENDDO
      ! (through statement 607) calculate fourier coefficients pb and vb at given ut
      ! for harmonics of tl.
      ! if ut is same as in previous call to this routine, skip to statement 603.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (ut.eq.tup) go to 603
      !***********************************************************************
          tup = UT
          icpt = 1
          tua = UT*hrang
          sl = SIN(tua)
          cl = COS(tua)
          fut(3) = 1.
          fut(2) = sq2*cl
          fut(4) = sq2*sl
          fut(1) = cl*fut(2) - sl*fut(4)
          fut(5) = cl*fut(4) + sl*fut(2)
          IF ( icpt /= 0 ) THEN
              DO 390 mpf = 1 , 9
                  pb(mpf) = 0.
                  vb(mpf) = 0.
                  DO 385 lp = 1 , 5
                      IF ( IUTav == 0 .OR. lp == 3 ) THEN
                          pb(mpf) = pb(mpf) + fut(lp)*pa(lp,mpf)
                          vb(mpf) = vb(mpf) + fut(lp)*va(lp,mpf)
                      ENDIF
                  385 ENDDO
              390 ENDDO
          ENDIF
      ! tl is magnetic local time.
          tl = UT + (XMLon-69.)/15.
      ! if tl is same as in previous call to this routine, skip to statement 630.
      ! 8/2/2 deactivate shortcut, since it won't necessarily work right if
      ! variables are not retained with a save statement.
      ! if (tl.eq.tlp) go to 630
      !***********************************************************************
          tlp = tl
      ! (through statement 610) calculate sines and cosines, times sq2, of harmonics
      ! of magnetic local time.
          tla = tl*hrang
          sl = SIN(tla)
          cl = COS(tla)
          sml(1) = sq2*sl
          cml(1) = sq2*cl
          DO 400 m = 2 , 4
              sml(m) = cl*sml(m-1) + sl*cml(m-1)
              cml(m) = cl*cml(m-1) - sl*sml(m-1)
          400 ENDDO
      ! calculate pot, ve, and vu by summing fourier coefficients multiplied by
      ! appropriate sines and cosines.
          POT = 0.
          VE = 0.
          VU = 0.
          DO 420 m = 1 , 4
              mpf = m + 5
              mmf = 5 - m
              xm = m
              POT = POT + pb(mpf)*sml(m) + pb(mmf)*cml(m)
              VE = VE + vb(mpf)*sml(m) + vb(mmf)*cml(m)
              VU = VU + xm*(pb(mpf)*cml(m)-pb(mmf)*sml(m))
          420 ENDDO
          POT = POT*st + pb(5)
          VE = VE + vb(5)
          VU = VU/rbt
          GOTO 500
      ENDIF
  ENDIF
  POT = 999.
  VU = 0.0
  VE = 0.0
  500 RETURN

end SUBROUTINE Richmond_lowlat_Efield_model


SUBROUTINE TUBES_LINEAR_INTERPOLATE(lp,flux_tube_max,re_apex,v_upwards_at_apex,time_step, &
                                    q_val, &
                                    ion_den, &
                                    ion_vel, &
                                    ion_temp, &
                                    electron_temp, &
                                    ion_den_final, &
                                    ion_vel_final, &
                                    ion_temp_final, &
                                    electron_temp_final)

!*********************************************
! *
! Interpolation routine for incorporating  *
! EXB drift into the 'fixed tubes' version *
! *
!*********************************************

  IMPLICIT NONE
  INTEGER, PARAMETER :: nlp = 170
  INTEGER, PARAMETER :: npts = 1115
  INTEGER :: i , lp , istop

  REAL(kind=8) :: re_apex(nlp)
  REAL(kind=8) :: reback, factor,factor2
  REAL(kind=8) :: lambda_m_inner,lambda_m_outer,lambda_m_stepped_backwards
  REAL(kind=8),  PARAMETER :: R0=6.370E06

  integer :: lp_out,lp_in,lpdum,ip,isouth,inorth, iww,ispecial
  integer :: flux_tube_max(nlp)
  integer :: ions

  REAL(kind=8) :: sqrt_part
  REAL(kind=8) :: v_upwards_at_apex
  REAL(kind=8) :: time_step
  REAL(kind=8) :: q_val(npts,nlp)

  REAL(kind=8) :: ion_den(9,npts,nlp)
  REAL(kind=8) :: ion_vel(9,npts,nlp)
  REAL(kind=8) :: ion_temp(npts,nlp)
  REAL(kind=8) :: electron_temp(npts,nlp)

  REAL(kind=8) :: ion_den_inner(9,npts)
  REAL(kind=8) :: ion_vel_inner(9,npts)
  REAL(kind=8) :: ion_temp_inner(npts)
  REAL(kind=8) :: electron_temp_inner(npts)

  REAL(kind=8) :: ion_den_outer(9,npts)
  REAL(kind=8) :: ion_vel_outer(9,npts)
  REAL(kind=8) :: ion_temp_outer(npts)
  REAL(kind=8) :: electron_temp_outer(npts)

  REAL(kind=8) :: ion_den_final(9,npts)
  REAL(kind=8) :: ion_vel_final(9,npts)
  REAL(kind=8) :: ion_temp_final(npts)
  REAL(kind=8) :: electron_temp_final(npts)


! step backwards to imagined previous flux-tube position....

! stop=1
! f(istop.eq.1) stop
  iww=0
! f(mp.eq.12.and.lp.eq.10) iww=1
!g
  reback=re_apex(lp) - (v_upwards_at_apex * time_step)
!g
  istop=1
! f(istop.eq.1) stop

! we now need  PARAMETERs on this imagined tube (by interpolation).....

  do lpdum=1,NLP
      if(reback > re_apex(lpdum)) then
          lp_out=lpdum-1
          lp_in=lpdum
          if(iww == 1) write(88,*) lp_out,lp_in
          goto 3456
      endif
  enddo
  lp_out=NLP-1
  lp_in=NLP
  3456 continue
  if(lp_out == 0) then
      lp_in = 2
      lp_out = 1
  endif
!g
! factor=(reback-re(mp,lp_in))/(re(mp,lp_out)-re(mp,lp_in))
!g
!g  For the apex code we use the modified apex latitude as our  PARAMETER for
!g  interpolating in the inwards/outwards direction....
!g

  lambda_m_inner = acos(sqrt((r0+90000.)/re_apex(lp_in)))
  lambda_m_outer = acos(sqrt((r0+90000.)/re_apex(lp_out)))

!cg   For the very small flux tubes it is possible that stepping backwards
!cg   leads to an apex height which is lower than the base height of 90km.
!cg   This leads to a NAN because the ACOS won't compute.
!cg   So we do the 'sqrt_part' bit below and check that this is less than 1.0

  sqrt_part = sqrt((r0+90000.)/reback)

  if(sqrt_part.le.0.999999) then
       lambda_m_stepped_backwards = acos(sqrt_part)
  else
       lambda_m_stepped_backwards = 0.0
  endif

  factor=(lambda_m_stepped_backwards - lambda_m_inner) / (lambda_m_outer - lambda_m_inner)
!g
!g  Need to flag when factor gets too big or small.....
!g  Factor of less than 0.0_prec means that the position of
!g  the imagined tubes is inside of the inner most REAL tube.
!g  The interpolation then becomes extrapolation (which isn't
!g  necessarily a problem in itself but can produce nasty
!g  negative densities if the gradients are large enough).
!g  In this case set the factor to 0.0_prec which means that you
!g  end up just using the inner flux_tube values...
!g  (Should be safe)
!g
  if(factor < 0.0) then
  ! write(6,*) 'FACTOR ',factor
  ! write(6,*) lp_in,lp_out
  ! write(6,*) reback,re(mp,lp_in),re(mp,lp_out)
  ! write(6,*) '     '
      factor = 0.0
  endif
!g
!g  The same could be true if factor is greater than 1.0
!g  The extrapolation isn't in itself a problem but large
!g  negative gradients in density (temp, whatever) could
!g  lead to negative values.  For the time being lets just
!g  flag this and see if it is ever a problem....
!g
  if(factor > 1.0) then
      write(6,*) 'FACTOR ',factor
      write(6,*) lp_in,lp_out
      write(6,*) reback,re_apex(lp_in),re_apex(lp_out)
      write(6,*) '     '
  endif
!g
!g first interpolate in q for the inner tube.....
!g
  if(iww == 1) write(88,*) '++++++++++ INNER TUBE ++++++++++++'
  do 7000 ip = 1 , flux_tube_max(lp)
      if(iww == 1) write(88,*) 'point ',ip
      ispecial=0
  !g
  !g loop over the inner tube ......
  !g
      do i = 1 , flux_tube_max(lp_in)
          if(q_val(i,lp_in) < q_val(ip,lp)) then
              isouth=i
              inorth=i-1
              goto 3488
          endif
      enddo
  !g
      ispecial=1
      isouth=flux_tube_max(lp_in)
      inorth=isouth-1

      3488 continue

      if(isouth == 1) then
          isouth= 2
          inorth= 1
          ispecial=2
      endif
      24 format(i4,2x,3f9.5)
  !g
      if(iww == 1) write(88,*) 'ispecial ',ispecial
      if(ispecial == 0) then

          factor2=(q_val(ip,lp)-q_val(isouth,lp_in))/(q_val(inorth,lp_in)-q_val(isouth,lp_in))

          if(iww == 1) write(88,*) 'factor2 inner',factor2

          do ions = 1 , 9
            ion_den_inner(ions,ip)=(factor2*(ion_den(ions,inorth,lp_in) - ion_den(ions,isouth,lp_in))) + ion_den(ions,isouth,lp_in)
            ion_vel_inner(ions,ip)=(factor2*(ion_vel(ions,inorth,lp_in) - ion_vel(ions,isouth,lp_in))) + ion_vel(ions,isouth,lp_in)
          enddo
          ion_temp_inner(ip)=(factor2*(ion_temp(inorth,lp_in) - ion_temp(isouth,lp_in))) + ion_temp(isouth,lp_in)
          electron_temp_inner(ip)=(factor2*(electron_temp(inorth,lp_in) - electron_temp(isouth,lp_in))) + electron_temp(isouth,lp_in)

      elseif(ispecial == 1) then

          do ions = 1 , 9
            ion_den_inner(ions,ip)=ion_den(ions,isouth,lp_in)
            ion_vel_inner(ions,ip)=ion_vel(ions,isouth,lp_in)
          enddo
          ion_temp_inner(ip)=ion_temp(isouth,lp_in)
          electron_temp_inner(ip)=electron_temp(isouth,lp_in)

      elseif(ispecial == 2) then

          do ions = 1 , 9
            ion_den_inner(ions,ip)=ion_den(ions,inorth,lp_in)
            ion_vel_inner(ions,ip)=ion_vel(ions,inorth,lp_in)
          enddo
          ion_temp_inner(ip)=ion_temp(inorth,lp_in)
          electron_temp_inner(ip)=electron_temp(inorth,lp_in)

      endif

      if(iww == 1) write(88,*) factor2,isouth,inorth,ion_den_inner(1,ip)
  7000 ENDDO
!g
!g then interpolate in q for the outer tube.....
!g
  if(iww == 1) write(88,*) '++++++++++ OUTER TUBE ++++++++++++'
  do 8000 ip = 1 , flux_tube_max(lp)
      if(iww == 1) write(88,*) 'point ',ip
      ispecial=0
      do i = 1 , flux_tube_max(lp_out)
          if(q_val(i,lp_out) < q_val(ip,lp)) then
              isouth=i
              inorth=i-1
              goto 3489
          endif
      enddo
  !g
      ispecial=1
      isouth = flux_tube_max(lp_out)
      inorth = isouth - 1
      3489 continue
      if(isouth == 1) then
          isouth = 2
          inorth = 1
          ispecial=2
      endif
  !g
      if(iww == 1) write(88,*) 'ispecial ',ispecial
      if(ispecial == 0) then

          factor2=(q_val(ip,lp)-q_val(isouth,lp_out))/(q_val(inorth,lp_out)-q_val(isouth,lp_out))

          if(iww == 1) write(88,*) 'factor2 outer',factor2

          do ions = 1 , 9
            ion_den_outer(ions,ip)=(factor2*(ion_den(ions,inorth,lp_out) - ion_den(ions,isouth,lp_out))) + ion_den(ions,isouth,lp_out)
            ion_vel_outer(ions,ip)=(factor2*(ion_vel(ions,inorth,lp_out) - ion_vel(ions,isouth,lp_out))) + ion_vel(ions,isouth,lp_out)
          enddo
          ion_temp_outer(ip)=(factor2*(ion_temp(inorth,lp_out) - ion_temp(isouth,lp_out))) + ion_temp(isouth,lp_out)
          electron_temp_outer(ip)=(factor2*(electron_temp(inorth,lp_out) - electron_temp(isouth,lp_out))) + electron_temp(isouth,lp_out)

      elseif(ispecial == 1) then

          do ions = 1 , 9
            ion_den_outer(ions,ip)=ion_den(ions,isouth,lp_out)
            ion_vel_outer(ions,ip)=ion_vel(ions,isouth,lp_out)
          enddo
          ion_temp_outer(ip)=ion_temp(isouth,lp_out)
          electron_temp_outer(ip)=electron_temp(isouth,lp_out)

      elseif(ispecial == 2) then

          do ions = 1 , 9
            ion_den_outer(ions,ip)=ion_den(ions,inorth,lp_out)
            ion_vel_outer(ions,ip)=ion_vel(ions,inorth,lp_out)
          enddo
          ion_temp_outer(ip)=ion_temp(inorth,lp_out)
          electron_temp_outer(ip)=electron_temp(inorth,lp_out)

      endif

!     if(iww == 1) write(88,*) factor2,isouth,inorth,ni1_out(ip)
!     if(iww == 1) write(88,*) ni(inorth,mp,1), &
!     ni(isouth,mp,1),inorth,isouth,mp,lp_out

  8000 ENDDO

! if(istop.eq.1) stop
!g
!g linear interpolation ......
!g
  do 9000 i = 1 , flux_tube_max(lp)

      do ions = 1 , 9
      ion_den_final(ions,i)=((ion_den_outer(ions,i) - ion_den_inner(ions,i))*factor) + ion_den_inner(ions,i)
      ion_vel_final(ions,i)=((ion_vel_outer(ions,i) - ion_vel_inner(ions,i))*factor) + ion_vel_inner(ions,i)
      enddo

      ion_temp_final(i)=((ion_temp_outer(i) - ion_temp_inner(i))*factor) + ion_temp_inner(i)
      electron_temp_final(i)=((electron_temp_outer(i) - electron_temp_inner(i))*factor) + electron_temp_inner(i)

  !g
  !g  If any of these  PARAMETERs have gone -ve we got some
  !g  trub.........
  !g
  9000 ENDDO

  return



end SUBROUTINE TUBES_LINEAR_INTERPOLATE

SUBROUTINE interpolate_in_q(flux_tube_max, &
                            flux_tube_max_target, &
                            q_val, &
                            q_val_target, &
                            ion_den, &
                            ion_vel, &
                            ion_temp, &
                            electron_temp, &
                            ion_den_target, &
                            ion_vel_target, &
                            ion_temp_target, &
                            electron_temp_target)

  IMPLICIT NONE
  INTEGER, PARAMETER :: npts = 1115
  INTEGER :: i , lp , istop

  REAL(kind=8) :: q_factor
  REAL(kind=8) :: lambda_m_target,lambda_m_outer,lambda_m_stepped_backwards

  integer :: mp_t,lp_t,lpdum,ip,isouth,inorth, iww,ispecial
  integer :: flux_tube_max
  integer :: flux_tube_max_target
  REAL(kind=8) :: q_val(npts),q_val_target(npts)
  integer :: ions

  REAL(kind=8) :: ion_den(4,npts)
  REAL(kind=8) :: ion_vel(4,npts)
  REAL(kind=8) :: ion_temp(npts)
  REAL(kind=8) :: electron_temp(npts)

  REAL(kind=8) :: ion_den_target(4,npts)
  REAL(kind=8) :: ion_vel_target(4,npts)
  REAL(kind=8) :: ion_temp_target(npts)
  REAL(kind=8) :: electron_temp_target(npts)

  do 7000 ip = 1 , flux_tube_max ! loop over the points on the current flux-tube

      ispecial=0
  !g
  !g loop over the target tube ......
  !g
      do i = 1 , flux_tube_max_target
          if(q_val_target(i) < q_val(ip)) then
              isouth=i
              inorth=i-1
              goto 3488
          endif
      enddo
  !g
      ispecial=1
      isouth=flux_tube_max_target
      inorth=isouth-1

      3488 continue

      if(isouth == 1) then
          isouth= 2
          inorth= 1
          ispecial=2
      endif
  !g
      if(ispecial == 0) then

          q_factor=(q_val(ip)-q_val_target(isouth))/(q_val_target(inorth)-q_val_target(isouth))

          do ions = 1 , 4
            ion_den_target(ions,ip)=(q_factor*(ion_den(ions,inorth) - ion_den(ions,isouth))) + ion_den(ions,isouth)
            ion_vel_target(ions,ip)=(q_factor*(ion_vel(ions,inorth) - ion_vel(ions,isouth))) + ion_vel(ions,isouth)
          enddo
          ion_temp_target(ip)=(q_factor*(ion_temp(inorth) - ion_temp(isouth))) + ion_temp(isouth)
          electron_temp_target(ip)=(q_factor*(electron_temp(inorth) - electron_temp(isouth))) + electron_temp(isouth)

      elseif(ispecial == 1) then

          do ions = 1 , 4
            ion_den_target(ions,ip)=ion_den(ions,isouth)
            ion_vel_target(ions,ip)=ion_vel(ions,isouth)
          enddo
          ion_temp_target(ip)=ion_temp(isouth)
          electron_temp_target(ip)=electron_temp(isouth)

      elseif(ispecial == 2) then

          do ions = 1 , 4
            ion_den_target(ions,ip)=ion_den(ions,inorth)
            ion_vel_target(ions,ip)=ion_vel(ions,inorth)
          enddo
          ion_temp_target(ip)=ion_temp(inorth)
          electron_temp_target(ip)=electron_temp(inorth)

      endif

  7000 ENDDO
  return



end SUBROUTINE interpolate_in_q



  SUBROUTINE Auroral_Precipitation( plasma, grid, neutrals, forcing, time_tracker )
  ! Previously : tiros_ionize_ipe
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    ! Local
    INTEGER, PARAMETER :: jmaxwell = 6
    INTEGER    :: i, lp, mp, m, j, iband, l
    INTEGER    :: i1, i2, j1, j2, k, kk
    INTEGER    :: tiros_activity_level
    REAL(prec) :: ratio(1:21)
    REAL(prec) :: mlt, dfac
    REAL(prec) :: ch, chi, diff, alpha
    REAL(prec) :: eflux, essa, dl_lower,dl_upper,qiont_lower,qiont_upper
    REAL(prec) :: ri, rj, th, THMagd, GW
    REAL(prec) :: pres(plasma % nFluxTube)
    REAL(prec) :: den(plasma % nFluxTube), grav(plasma % nFluxTube)
    REAL(prec) :: ntot(plasma % nFluxTube), meanmass(plasma % nFluxTube)
    REAL(prec) :: rno, rang, pr, ratioz, rlamz, q
    REAL(prec) :: TE11(21),TE15(21)
    REAL(prec) :: ratio_ch,en_maxwell(jmaxwell),dl(jmaxwell)
    REAL(prec) :: qion_maxwell(plasma % nFluxTube)
    REAL(prec) :: en(1:15), width(1:15), RLAM(1:21), ionchr(1:21)

    REAL(prec), PARAMETER :: fc = 1.6e-06_prec

       en(1:15) = (/ 0.37_prec, 0.6_prec, 0.92_prec, 1.37_prec, 2.01_prec, &
                     2.91_prec, 4.19_prec, 6.0_prec, 8.56_prec, 12.18_prec, &
                     17.3_prec, 24.49_prec, 36.66_prec, 54.77_prec, 81.82_prec /)

       width(1:15) = (/ 0.158_prec, 0.315_prec, 0.315_prec, 0.63_prec, 0.631_prec, &
                        1.261_prec, 1.26_prec, 2.522_prec, 2.522_prec, 5.043_prec, &
                        5.043_prec, 10.0_prec, 14.81_prec, 22.13_prec, 33.06_prec /)

       rlam(1:21) = (/ 1.49_prec, 1.52_prec, 1.51_prec, 1.48_prec, 1.43_prec, &
                                  1.37_prec, 1.30_prec, 1.22_prec, 1.12_prec, 1.01_prec, &
                                  0.895_prec, 0.785_prec, 0.650_prec, 0.540_prec, 0.415_prec, &
                                  0.320_prec, 0.225_prec, 0.14_prec, 0.08_prec, 0.04_prec, 0.0_prec /)

       ionchr(1:21) = (/ 0.378_prec, 0.458_prec, 0.616_prec, 0.773_prec, 0.913_prec, &
                                     1.088_prec, 1.403_prec, 1.718_prec, 2.033_prec, 2.349_prec, &
                                     2.979_prec, 3.610_prec, 4.250_prec, 4.780_prec, 6.130_prec, &
                                     7.392_prec, 8.653_prec, 9.914_prec, 12.436_prec, 14.957_prec, 17.479_prec /)

      tiros_activity_level = forcing % nhemi_power_index( forcing % current_index )
      GW                   = forcing % nhemi_power( forcing % current_index )

      DO m = 1, 21
        ratio(m) = REAL( (m-1), prec )*0.05_prec
      ENDDO

      DO iband=1,21

        te15(iband)=0.0_prec
        te11(iband)=0.0_prec

        DO m = 1, 15
          te15(iband) = te15(iband) + fc * forcing % djspectra(m,iband)*en(m)*width(m)
        ENDDO
        DO m = 1,11
          te11(iband) = te15(iband) + fc * forcing % djspectra(m,iband)*en(m)*width(m)
        ENDDO

      ENDDO

      DO j = 1, jmaxwell
        en_maxwell(j) = REAL( j, prec )*0.05_prec - 0.025_prec
      ENDDO

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, grid % NLP

          DO i=1,grid % flux_tube_max(lp)
            plasma % ionization_rates(1,i,lp,mp) = 0.0_prec
            plasma % ionization_rates(2,i,lp,mp) = 0.0_prec
            plasma % ionization_rates(3,i,lp,mp) = 0.0_prec
            plasma % ionization_rates(4,i,lp,mp) = 0.0_prec
            qion_maxwell(i)                      = 0.0_prec
          ENDDO

          ! convert magnetic latitude from radians to degrees
          ! Use the geomagnetic latitude of the foot point
          thmagd = rtd * ( half_pi - grid % magnetic_colatitude(1,lp) )
          th = abs(thmagd) - 50.0_prec

          IF ( abs(thmagd) > 50.0_prec )  THEN

            mlt = time_tracker % utime/3600.0_prec + &
                  rtd * grid % longitude(1,lp,mp)/15.0_prec

            essa = (mlt + 12.0_prec)*15.0_prec

            IF ( essa >= 360.0_prec ) THEN
              essa = essa - 360.0_prec
            ELSEIF ( essa < 0.0_prec ) THEN
              essa = essa + 360.0_prec
            ENDIF

            dfac = 1.0_prec
            IF( tiros_activity_level > 9 .AND. GW > 96.0_prec ) THEN
              dfac = GW/96.0_prec
            ENDIF

            l = tiros_activity_level - 2
            IF ( l < 1 ) THEN
              l = 1
            ELSEIF ( l > 7 ) THEN
              l = 7
            ENDIF

            ri = essa/18.0_prec + 11.0_prec
            i1 = ri
            ri = ri - i1
            IF ( i1 > 20 ) THEN
              i1 = i1 - 20
            ENDIF

            i2 = i1 + 1
            IF ( i2 > 20 ) THEN
              i2 = i2 - 20
            ENDIF

            rj = th/2.0_prec + 1.0_prec
            j1 = rj
            rj = rj - j1
            j2 = j1 + 1

            eflux = rj*ri*forcing % emaps(j2,i2,l) + &
                    (1.0_prec-rj)*ri*forcing % emaps(j1,i2,l) + &
                    rj*(1.0_prec-ri)*forcing % emaps(j2,i1,l) + &
                    (1.0_prec-rj)*(1.0_prec-ri)*forcing % emaps(j1,i1,l)
            eflux = 10.0_prec**(eflux)/1000.0_prec

            ch = rj*ri*forcing % cmaps(j2,i2,l) + &
                 (1.0_prec-rj)*ri*forcing % cmaps(j1,i2,l) + &
                 rj*(1.0_prec-ri)*forcing % cmaps(j2,i1,l) + &
                 (1.0_prec-rj)*(1.0_prec-ri)*forcing % cmaps(j1,i1,l)


            IF ( ch < 0.378_prec ) THEN
              ch = 0.379_prec
            ENDIF

            DO kk = 2 , 21
              IF ( ch <= ionchr(kk) ) THEN
                k = kk - 1
                EXIT
              ENDIF
            ENDDO

            kk = k+1
            chi = ch - ionchr(k)
            diff = ionchr(kk) - ionchr(k)
            ratio_ch = chi/diff

            DO i = 1 , grid % flux_tube_max(lp)

              IF ( grid % altitude(i,lp) <= 1.0e+06_prec )THEN

                grav(i) = -grid % grx(i,lp,mp)

                ntot(i) = ( neutrals % oxygen(i,lp,mp)  +&
                            neutrals % molecular_oxygen(i,lp,mp) +&
                            neutrals % molecular_nitrogen(i,lp,mp) +&
                            neutrals % hydrogen(i,lp,mp)  +&
                            neutrals % helium(i,lp,mp) )

                pres(i) = ntot(i)*kBoltz*neutrals % temperature(i,lp,mp)

                meanmass(i) = ( neutrals % oxygen(i,lp,mp)*O_mass   +&
                                neutrals % molecular_oxygen(i,lp,mp)*O2_mass +&
                                neutrals % molecular_nitrogen(i,lp,mp)*N2_mass +&
                                neutrals % hydrogen(i,lp,mp)*H_mass   +&
                                neutrals % helium(i,lp,mp)*He_mass )/ntot(i)

                den(i) = pres(i)*meanmass(i)/(gscon*neutrals % temperature(i,lp,mp))

                DO l = 1 , 15

                  dl_lower = forcing % djspectra(l,k)
                  dl_upper = forcing % djspectra(l,kk)
                  rang = 4.57e-05_prec*en(l)**1.75_prec
                  pr = rang*grav(i)
                  ratioz = pres(i)/pr

                  IF ( ratioz > 1.0_prec ) THEN

                    rlamz = 0.0_prec

                  ELSE

                    DO m = 2 , 21 !! Used to loop from 1,21
                      IF ( ratioz <= ratio(m) ) THEN
                        rlamz = rlam(m-1) + (ratioz-ratio(m-1))*&
                                (rlam(m)-rlam(m-1))/(ratio(m)-ratio(m-1))
                        EXIT
                      ENDIF
                    ENDDO

                  ENDIF

                  qiont_lower = den(i)*en(l)*rlamz*dl_lower*width(l)*1.0e+07_prec/rang/e0
                  qiont_upper = den(i)*en(l)*rlamz*dl_upper*width(l)*1.0e+07_prec/rang/e0

                  plasma % ionization_rates(1,i,lp,mp) = plasma % ionization_rates(1,i,lp,mp) + &
                                                         (ratio_ch*qiont_upper+(1.0_prec-ratio_ch)*qiont_lower)*&
                                                         eflux/(ratio_ch*te11(kk)+(1.0_prec-ratio_ch)*te11(k))

                ENDDO

                alpha = ch/2.0_prec
                rno = eflux*6.24e+12_prec/2.0_prec/alpha**3

                DO l = 1 , JMAXWELL

                  dl(l) = rno*en_maxwell(l)*EXP(-en_maxwell(l)/alpha)
                  rang = 4.57e-05_prec*en_maxwell(l)**1.75_prec
                  pr = rang*grav(i)
                  ratioz = pres(i)/pr

                  IF ( ratioz > 1.0_prec ) THEN
                    rlamz = 0.0_prec
                  ELSE

                    DO m = 2 , 21 !! Used to loop from 1, 21
                      IF ( ratioz <= ratio(m) ) THEN
                        rlamz = rlam(m-1) + (ratioz-ratio(m-1))*(rlam(m)-rlam(m-1))/(ratio(m)-ratio(m-1))
                        EXIT
                      ENDIF
                    ENDDO

                  ENDIF

                  qion_maxwell(i) = qion_maxwell(i) + den(i)*en_maxwell(l)*rlamz*dl(l)*width_maxwell/rang/e0

                ENDDO

                plasma % ionization_rates(1,i,lp,mp) = ( plasma % ionization_rates(1,i,lp,mp) + qion_maxwell(i) )*dfac

                q = plasma % ionization_rates(1,i,lp,mp)/( 0.92_prec*neutrals % molecular_nitrogen(i,lp,mp) +&
                                                           1.50_prec*neutrals % molecular_oxygen(i,lp,mp)  +&
                                                           0.56_prec*neutrals % oxygen(i,lp,mp))

                plasma % ionization_rates(2,i,lp,mp) = ( 0.50_prec*neutrals % molecular_oxygen(i,lp,mp) +& ! O+ Rate
                                                         0.56_prec*neutrals % oxygen(i,lp,mp) )*q

                plasma % ionization_rates(3,i,lp,mp) = neutrals % molecular_oxygen(i,lp,mp)*q             ! O2+ Rate

                plasma % ionization_rates(4,i,lp,mp) = 0.92_prec*neutrals % molecular_nitrogen(i,lp,mp)*q ! N2+ Rate

              ENDIF

            ENDDO
          ENDIF

        ENDDO
      ENDDO

  END SUBROUTINE Auroral_Precipitation

  SUBROUTINE FLIP_Wrapper( plasma, grid, neutrals, forcing, time_tracker, flip_time_step, &
    nflag_t,nflag_d )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    REAL(prec), INTENT(in)             :: flip_time_step
    ! Local
    INTEGER  :: i, lp, mp, iprint, ii
    INTEGER  :: JMINX, JMAXX
    INTEGER  :: EFLAG(11,11)
    INTEGER  :: nflag_t(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    INTEGER  :: nflag_d(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    REAL(dp) :: PCO, UTHR, HPEQ_flip
    REAL(dp) :: ZX(1:grid % nFluxTube)
    REAL(dp) :: SLX(1:grid % nFluxTube)
    REAL(dp) :: GLX(1:grid % nFluxTube)
    REAL(dp) :: BMX(1:grid % nFluxTube)
    REAL(dp) :: GRX(1:grid % nFluxTube)
    REAL(dp) :: OX(1:grid % nFluxTube)
    REAL(dp) :: HX(1:grid % nFluxTube)
    REAL(dp) :: N2X(1:grid % nFluxTube)
    REAL(dp) :: O2X(1:grid % nFluxTube)
    REAL(dp) :: HEX(1:grid % nFluxTube)
    REAL(dp) :: N4SX(1:grid % nFluxTube)
    REAL(dp) :: TNX(1:grid % nFluxTube)
    REAL(dp) :: TINFX(1:grid % nFluxTube)
    REAL(dp) :: UNX(1:grid % nFluxTube)
    REAL(dp) :: XIONNX(1:9,1:grid % nFluxTube)
    REAL(dp) :: XIONVX(1:9,1:grid % nFluxTube)
    REAL(dp) :: TE_TIX(1:3,1:grid % nFluxTube)
    REAL(dp) :: EHTX(1:3,1:grid % nFluxTube)
    REAL(dp) :: AUR_PROD(1:3,1:grid % nFluxTube)
    REAL(dp) :: NNOX(1:grid % nFluxTube)
    REAL(dp) :: NHEAT(1:grid % nFluxTube)
    REAL(dp) :: SZA(1:grid % nFluxTube)
    REAL(dp) :: dotprod, sini
    REAL(sp) :: F107D, F107A

      F107D = forcing % f107( forcing % current_index )
      F107A = forcing % f107_81day_avg( forcing % current_index )
      UTHR  = time_tracker % hour + (time_tracker % minute) / 60.0
      HPEQ_flip = HPEQ
      nflag_t = 0
      nflag_d = 0

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, plasma % NLP

          ! Copy over the grid information (for now)
          ZX(1:grid % flux_tube_max(lp))  = grid % altitude(1:grid % flux_tube_max(lp),lp)/1000.0_prec !convert from m to km
          PCO                             = grid % p_value(lp)  !Pvalue is a single value
          SLX(1:grid % flux_tube_max(lp)) = grid % foot_point_distance(1:grid % flux_tube_max(lp),lp,mp)
          GLX(1:grid % flux_tube_max(lp)) = half_pi - grid % magnetic_colatitude(1:grid % flux_tube_max(lp),lp)  ! magnetic latitude [radians]
          BMX(1:grid % flux_tube_max(lp)) = grid % magnetic_field_strength(1:grid % flux_tube_max(lp),lp,mp)   !Tesla
! Need a -SinI factor for the GRX
          DO i=1, grid % flux_tube_max(lp)
            dotprod = SQRT( grid % apex_d_vectors(1,3,i,lp,mp)**2 + &
                        grid % apex_d_vectors(2,3,i,lp,mp)**2 + &
                        grid % apex_d_vectors(3,3,i,lp,mp)**2  )

            IF( Almost_Equal( dotprod, 0.0_prec )) THEN
              SinI = 0.0_prec
            ELSE
              SinI = grid % apex_d_vectors(3,3,i,lp,mp)/dotprod
            ENDIF
            GRX(i) = 0.0_prec - (grid % grx(i,lp,mp) * SinI)
          ENDDO

          ! Copy over neutrals
          OX(1:grid % flux_tube_max(lp))    = neutrals % oxygen(1:grid % flux_tube_max(lp),lp,mp) !(m-3)
          HX(1:grid % flux_tube_max(lp))    = neutrals % hydrogen(1:grid % flux_tube_max(lp),lp,mp)
          N2X(1:grid % flux_tube_max(lp))   = neutrals % molecular_nitrogen(1:grid % flux_tube_max(lp),lp,mp)
          O2X(1:grid % flux_tube_max(lp))   = neutrals % molecular_oxygen(1:grid % flux_tube_max(lp),lp,mp)
          HEX(1:grid % flux_tube_max(lp))   = neutrals % helium(1:grid % flux_tube_max(lp),lp,mp)
          N4SX(1:grid % flux_tube_max(lp))  = neutrals % nitrogen(1:grid % flux_tube_max(lp),lp,mp)
          TNX(1:grid % flux_tube_max(lp))   = neutrals % temperature(1:grid % flux_tube_max(lp),lp,mp)
          TINFX(1:grid % flux_tube_max(lp)) = neutrals % temperature_inf(1:grid % flux_tube_max(lp),lp,mp)
          UNX(1:grid % flux_tube_max(lp))   = -neutrals % velocity_apex(3,1:grid % flux_tube_max(lp),lp,mp)
          SZA(1:grid % flux_tube_max(lp)) = Solar_Zenith_Angle( time_tracker % utime, &
                                                                time_tracker % day_of_year, &
                                                                grid % colatitude(1:grid % flux_tube_max(lp),lp,mp), &
                                                                grid % longitude(1:grid % flux_tube_max(lp),lp,mp), &
                                                                grid % flux_tube_max(lp) )

          AUR_PROD(1:3,1:grid % flux_tube_max(lp)) = plasma % ionization_rates(2:4,1:grid % flux_tube_max(lp),lp,mp)

          EFLAG(1:11,1:11) = 0

!         write(966,*) 'GHGM 966 ', mp,lp, plasma % ion_densities(6,1,lp,mp)

          DO i=1, grid % flux_tube_max(lp)

! GHGM - for some reason the O2Plus at the very bottom point can be negative
! GHGM - sort that out
            do ii = 1,9
            if (plasma % ion_densities(ii,i,lp,mp).lt.0.0) then
!            write(966,*) 'NEGATIVE ',ii, i, mp, lp, ' sorting '
             plasma % ion_densities(ii,i,lp,mp) = 0.1  ! just some default value
            endif
            enddo
! GHGM end of that bodge


            ! Ion Densities
            XIONNX(1:9,i) = plasma % ion_densities(1:9,i,lp,mp)
            ! Along Flux Tube Ion Velocities
            XIONVX(1:9,i) = plasma % ion_velocities(1:9,i,lp,mp)

            ! Ion Temperatures
            TE_TIX(1,i) = plasma % ion_temperature(i,lp,mp)
            TE_TIX(2,i) = plasma % ion_temperature(i,lp,mp)
            ! Electron Temperature
            TE_TIX(3,i) = plasma % electron_temperature(i,lp,mp)

            EHTX(1:3,i) = 0.0_dp
            NNOX(i)     = 0.0_dp
            NHEAT(i)    = 0.0_dp

          ENDDO

          JMINX = 1
          JMAXX = grid % flux_tube_max(lp)

          CALL CTIPINT( JMINX, & !.. index of the first point on the field line
                        JMAXX, & !.. index of the last point on the field line
                        grid % flux_tube_max(lp), & !.. CTIPe array dimension, must equal to FLDIM
                        ZX(1:JMAXX), & !.. array, altitude (km)
                        PCO, & !.. p coordinate (L-shell)
                        SLX(1:JMAXX), & !.. array, distance of point from northern hemisphere (meter)
                        GLX(1:JMAXX), & !.. array, magnetic latitude (radians)
                        BMX(1:JMAXX), & !.. array, magnetic field strength, (Tesla)
                        GRX(1:JMAXX), & !.. array, gravity, m2 s-1
                        OX(1:JMAXX), & !.. array, O density (m-3)
                        HX(1:JMAXX), & !.. array, H density (m-3)
                        N2X(1:JMAXX), & !.. array, N2 density (cm-3)
                        O2X(1:JMAXX), & !.. array, O2 density (cm-3)
                        HEX(1:JMAXX), & !.. array, He density (cm-3)
                        N4SX(1:JMAXX), & !.. array, N(4S) density (cm-3)
                        INNO, &
                        NNOX(1:JMAXX), & !.. array, NO density (cm-3)
                        TNX(1:JMAXX), & !.. array, Neutral temperature (K)
                        TINFX(1:JMAXX), & !.. array, Exospheric Neutral temperature (K)
                        UNX(1:JMAXX), & !.. array, Neutral wind (m/s), field aligned component, positive SOUTHward
                        flip_time_step, & !.. CTIPe time step (secs)
                        DTMIN, & !.. Minimum time step allowed (>=10 secs?)
                        F107D, & !.. Daily F10.7
                        F107A, & !.. 81 day average F10.7
                        SZA(1:JMAXX), & !.. Solar Zenith angle (radians)
                        FPAS, & !.. Pitch angle scattering fraction
                        HPEQ_flip, & !.. Sets initial equatorial H+ density. See declaration below
                        HEPRAT, & !.. Intial He+/H+ ratio (.01 to 1.0)
                        COLFACX, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
                        IHEPLS, &
                        INPLS, &
                        UTHR, &  !.. Universal time in hours
                        EHTX(1:3,1:JMAXX), & !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
                        AUR_PROD(1:3,1:JMAXX), & ! IN 2D array, ionization rates for O+, O2+, and N2+
                        TE_TIX(1:3,1:JMAXX), & !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
                        XIONNX(1:9,1:JMAXX), &
                        XIONVX(1:9,1:JMAXX), & !.. IN/OUT: 2D array, Storage for ion densities and velocities
                        NHEAT(1:JMAXX), & !.. OUT: array, Neutral heating rate (eV/cm^3/s)
                        EFLAG,mp,lp,nflag_t(lp,mp),nflag_d(lp,mp) ) !.. OUT: 2D array, Error Flags


          DO i=1, grid % flux_tube_max(lp)

            ! Ion Densities
            plasma % ion_densities(1:9,i,lp,mp) = XIONNX(1:9,i)
            ! Electron Density
!           plasma % electron_density(i,lp,mp) = XIONNX(1,i) + XIONNX(2,i) + XIONNX(3,i) + XIONNX(4,i) + XIONNX(5,i) + XIONNX(6,i) + XIONNX(7,i) + XIONNX(8,i) + XIONNX(9,i)
! Just make electron density be Oplus + Hplus for now
            plasma % electron_density(i,lp,mp) = XIONNX(1,i) + XIONNX(2,i)
            ! Along Flux Tube Ion Velocities
            plasma % ion_velocities(1:9,i,lp,mp) = XIONVX(1:9,i)

            ! Ion Temperatures
            plasma % ion_temperature(i,lp,mp) = TE_TIX(1,i)
            ! Electron Temperature
            plasma % electron_temperature(i,lp,mp) = TE_TIX(3,i)

          ENDDO


      ENDDO
    ENDDO

  END SUBROUTINE FLIP_Wrapper
!
  FUNCTION Solar_Zenith_Angle( utime, day, colatitude, longitude, flux_tube_max ) RESULT( sza )
    IMPLICIT NONE
    REAL(prec) :: utime
    INTEGER    :: day
    INTEGER    :: flux_tube_max               ! Number of flux tube points on this flux tube
    REAL(prec) :: colatitude(1:flux_tube_max) ! Geographic Colatitude
    REAL(prec) :: longitude(1:flux_tube_max)  ! Geographic longitude
    REAL(dp)   :: sza(1:flux_tube_max)        ! Solar Zenith Angle ( Radians )
    ! Local
    REAL(dp) :: ty, sda, ssa, rlt, rlat, utime12, cos_sza
    INTEGER  :: i

        sza = 0.0_dp
        DO i=1, flux_tube_max

          ty = ( REAL(day, prec) + 15.5_prec )*12.0_prec/365.0_prec
          IF ( ty > 12.0_prec )THEN
            ty = ty - 12.0_prec
          ENDIF

          sda = ATAN(0.434_prec*SIN(pi/6.0_prec*(ty-3.17_prec)))

          IF ( utime >= 43200.0_prec ) THEN
            utime12 = utime - 43200.0_prec
          ELSE
            utime12 = utime + 43200.0_prec
          ENDIF

          ssa = rtd * longitude(i) + utime12/240.0_prec
          rlt = 180.0_prec + ssa
          IF ( rlt > 360.0_prec ) THEN
            rlt = rlt - 360.0_prec
          ENDIF
          rlt = dtr * rlt

          rlat = half_pi - colatitude(i)
          cos_sza = -COS(rlat)*COS(sda)*COS(rlt)+SIN(rlat)*SIN(sda)

          IF( cos_sza > 1.0_prec )THEN
            cos_sza = 1.0_prec
          ELSEIF( cos_sza < -1.0_prec )THEN
            cos_sza = -1.0_prec
          ENDIF

          sza(i)  = ACOS( cos_sza )

       END DO

  END FUNCTION Solar_Zenith_Angle

  SUBROUTINE Interpolate_to_GeographicGrid_IPE_Plasma( plasma, grid )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)         :: grid
    ! Local
    INTEGER :: i

     DO i = 1, n_ion_species
        CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_densities(i,1:grid % nFluxtube,1:grid % NLP,1:grid % NMP), plasma % geo_ion_densities(i,:,:,:) )
     ENDDO

     CALL grid % Interpolate_to_Geographic_Grid( plasma % ion_temperature(1:grid % nFluxtube,1:grid % NLP,1:grid % NMP), plasma % geo_ion_temperature )
     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_temperature(1:grid % nFluxtube,1:grid % NLP,1:grid % NMP), plasma % geo_electron_temperature )

!    CALL Calculate_Tec_Nmf2_Hmf2( plasma % geo_ion_densities , plasma % geo_tec , plasma % geo_nmf2 , plasma % geo_hmf2 , nlon_geo, nlat_geo, nheights_geo)
!    plasma % electron_density(1:grid % nFluxtube,1:grid % NLP,1:grid % NMP) = &
!     plasma % ion_densities(1,1:grid % nFluxtube,1:grid % NLP,1:grid % NMP) + &
!     plasma % ion_densities(2,1:grid % nFluxtube,1:grid % NLP,1:grid % NMP)

     CALL grid % Interpolate_to_Geographic_Grid( plasma % electron_density2(1:grid % nFluxtube,1:grid % NLP,1:grid % NMP), plasma % geo_electron_density )
    ! CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(1,:,:,:), plasma % geo_ionization_rates(1,:,:,:) )
    ! CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(2,:,:,:), plasma % geo_ionization_rates(2,:,:,:) )
    ! CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(3,:,:,:), plasma % geo_ionization_rates(3,:,:,:) )
    ! CALL grid % Interpolate_to_Geographic_Grid( plasma % ionization_rates(4,:,:,:), plasma % geo_ionization_rates(4,:,:,:) )


  END SUBROUTINE Interpolate_to_GeographicGrid_IPE_Plasma


  SUBROUTINE Calculate_Tec_Nmf2_Hmf2( geo_ion_densities , tec , nmf2 , hmf2 , nlon_geo, nlat_geo, nheights_geo)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: nlon_geo , nlat_geo , nheights_geo
    REAL(dp), INTENT(in) :: geo_ion_densities(1:9,nlon_geo, nlat_geo, nheights_geo)
    REAL(dp), INTENT(out) :: tec(nlon_geo,nlat_geo) , nmf2(nlon_geo,nlat_geo) , hmf2(nlon_geo,nlat_geo)
    INTEGER :: i_max_location , i1 , i2 , i3 , iheight , ilon , ilat
    REAL(dp) :: x1 , x2 , x3 , y1 , y2 , y3 , a , b , c
    REAL(dp) :: electron_density_profile(nheights_geo)
    REAL(dp) :: geo_electron_density(nlon_geo, nlat_geo, nheights_geo)

! GHGM make electron density just O+ and H+ for now....
    geo_electron_density(:,:,:) = geo_ion_densities(1,:,:,:) +&
                                geo_ion_densities(2,:,:,:) ! +&
!                               geo_ion_densities(3,:,:,:) +&
!                               geo_ion_densities(4,:,:,:)  ! +&
!                               geo_ion_densities(5,:,:,:) +&
!                               geo_ion_densities(6,:,:,:) +&
!                               geo_ion_densities(7,:,:,:) +&
!                               geo_ion_densities(8,:,:,:) +&
!                               geo_ion_densities(9,:,:,:)


! Total Electron Content....
  tec = 0.0
  do 500 iheight = 1 , nheights_geo

    tec(1:nlon_geo,1:nlat_geo) = tec(1:nlon_geo,1:nlat_geo) + (geo_electron_density(1:nlon_geo,1:nlat_geo,iheight) * 5000.0)
  500 enddo

  tec = 1.0e-16_prec * tec

! NmF2 and hmF2....
  do 660 ilon = 1 , 90
  do 650 ilat = 1 , 91
    electron_density_profile = geo_electron_density(ilon,ilat,1:nheights_geo)

    i_max_location = maxval(maxloc(electron_density_profile))

    i1 = i_max_location - 1
    i2 = i_max_location
    i3 = i_max_location + 1
    x1 = (float(i1 - 1) * 5.0) + 90.0
    x2 = (float(i2 - 1) * 5.0) + 90.0
    x3 = (float(i3 - 1) * 5.0) + 90.0
    y1 = electron_density_profile(i1)
    y2 = electron_density_profile(i2)
    y3 = electron_density_profile(i3)

    c = (x3*y1 - x3*y2 + x1*y2 + x2*y3 - x2*y1 - x1*y3) / (x3*x1*x1 - x3*x2*x2 + x1*x2*x2 - x2*x1*x1 + x2*x3*x3 - x1*x3*x3)
    b = (y2 - (c*x2*x2) + (c*x1*x1) - y1) / (x2 - x1)
    a = y1 - (b*x1) - (c*x1*x1)
    hmf2(ilon,ilat) = (0.0 - b) / (2*c)
    nmf2(ilon,ilat) = a + (b*hmf2(ilon,ilat)) + (c*hmf2(ilon,ilat)*hmf2(ilon,ilat))

  650 enddo
  660 enddo

  END SUBROUTINE Calculate_Tec_Nmf2_Hmf2


  SUBROUTINE Read_Legacy_Input_IPE_Plasma( plasma, grid, filename )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    CHARACTER(100), INTENT(in)         :: filename
    ! Local
    INTEGER               :: ispec, i, lp, mp, ii, fUnit
    REAL(4), ALLOCATABLE :: dumm(:,:,:)

      ALLOCATE( dumm(1:grid % npts2d, 1:grid % NMP, 1:12) )

      dumm = 0.0_sp

      OPEN( UNIT   = NewUnit( fUnit ), &
            FILE   = TRIM( filename ), &
            FORM   = 'UNFORMATTED', &
            STATUS = 'OLD' )


      DO ispec=1,12
        READ( fUnit ) dumm(:,:,ispec)
      ENDDO

      CLOSE( fUnit )

      DO mp = 1, grid % NMP

        ii = 0

        DO lp = 1, grid % NLP
          DO i = 1, grid % flux_tube_max(lp)

            ii = ii + 1
            plasma % ion_densities(1:9,i,lp,mp)    = dumm(ii,mp,1:9)
            plasma % electron_temperature(i,lp,mp) = dumm(ii,mp,10)
            plasma % ion_temperature(i,lp,mp)      = dumm(ii,mp,11)
 !           plasma % ion_velocities(1:3,i,lp,mp)   = dumm(ii,mp,13:15)

          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE( dumm )

  END SUBROUTINE Read_Legacy_Input_IPE_Plasma

  SUBROUTINE Calculate_Field_Line_Integrals( plasma, grid, neutrals, mpi_layer )


    ! this routine takes in the plasma, grid, and neutral fields and returns the six conductivities on IPE
    ! mp,lp grid. the output of this routine will likely need to be interpolated onto the dynamo solver grid
    ! before calling the dynamo solver
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout)  :: plasma
    TYPE( IPE_Grid ), INTENT(in)        :: grid
    TYPE( IPE_Neutrals ), INTENT(in)    :: neutrals
    TYPE( IPE_MPI_Layer ), INTENT(in)   :: mpi_layer
    ! Local
    INTEGER    :: i, ihem, istart, istop, istep, mp, lp
    REAL(prec) :: Ue(1:2),r_factor,qe_fac
    REAL(prec) :: ds, apex_d1d1, apex_d2d2, apex_d1d2, apex_D, apex_BMAG
    REAL(prec) :: effective_temp, o_cm3, o2_cm3, n2_cm3,abs_ds
    REAL(prec) :: ion_neutral_collisionfrequency, ion_mass_amu, rnu_o2p, rnu_op, rnu_nop, rnu_ne
    REAL(prec) ::  integral513
    REAL(prec) ::  integral514
    REAL(prec) ::  integral517
    REAL(prec) ::  integral518
    REAL(prec) ::  integral519
    REAL(prec) ::  integral520
    REAL(prec) ::  electron_density
    REAL(prec) ::  electron_charge_Coulombs
    REAL(prec) ::  sigma_ped,sigma_hall
#ifdef HAVE_MPI
    INTEGER :: mpierror, sendcount
    REAL(prec) :: conductivities(1:6,1:2,1:grid % NLP, grid % mp_low:grid % mp_high)
#endif

! Need to pick up only latitudes match dynamo grid for the calculation!

   DO mp = grid % mp_low , grid % mp_high
     DO lp = 1, grid % NLP

        DO ihem = 1, 2

          if (ihem==2) then ! southern hemisphere
             istart=grid % flux_tube_max(lp)
             istop=grid % flux_tube_midpoint(lp)
             istep=-1
          else if (ihem==1) then ! northern hemisphere
             istart=1
             istop=grid % flux_tube_midpoint(lp) - 1
             istep=1
          endif

          integral513=0.0_prec
          integral514=0.0_prec
          integral517=0.0_prec
          integral518=0.0_prec
          integral519=0.0_prec
          integral520=0.0_prec

          DO i = istart, istop, istep

            IF ( i==grid % flux_tube_max(lp) ) THEN
              ds = grid % foot_point_distance(i,lp,mp) - &
                      grid % foot_point_distance(i-1,lp,mp)
            ELSE
              ds = grid % foot_point_distance(i+1,lp,mp) - &
                      grid % foot_point_distance(i,lp,mp)
            END IF

            Ue(1) =  grid % apex_d_vectors(1,1,i,lp,mp)*neutrals % velocity_geographic(1,i,lp,mp) + &
                     grid % apex_d_vectors(2,1,i,lp,mp)*neutrals % velocity_geographic(2,i,lp,mp) + &
                     grid % apex_d_vectors(3,1,i,lp,mp)*neutrals % velocity_geographic(3,i,lp,mp)

            Ue(2) =  grid % apex_d_vectors(1,2,i,lp,mp)*neutrals % velocity_geographic(1,i,lp,mp) + &
                     grid % apex_d_vectors(2,2,i,lp,mp)*neutrals % velocity_geographic(2,i,lp,mp) + &
                     grid % apex_d_vectors(3,2,i,lp,mp)*neutrals % velocity_geographic(3,i,lp,mp)

            apex_d1d1 = grid % apex_d_vectors(1,1,i,lp,mp)**2 +&
                        grid % apex_d_vectors(2,1,i,lp,mp)**2 +&
                        grid % apex_d_vectors(3,1,i,lp,mp)**2

            apex_d2d2 = grid % apex_d_vectors(1,2,i,lp,mp)**2 +&
                        grid % apex_d_vectors(2,2,i,lp,mp)**2 +&
                        grid % apex_d_vectors(3,2,i,lp,mp)**2

            apex_d1d2 = grid % apex_d_vectors(1,1,i,lp,mp)*grid % apex_d_vectors(1,2,i,lp,mp) +&
                        grid % apex_d_vectors(2,1,i,lp,mp)*grid % apex_d_vectors(2,2,i,lp,mp) +&
                        grid % apex_d_vectors(3,1,i,lp,mp)*grid % apex_d_vectors(3,2,i,lp,mp)

            apex_D = grid % d1xd2_magnitude(i,lp,mp)
            apex_BMAG = grid % magnetic_field_strength(i,lp,mp)
!         if (lp.gt.135.and.lp.le.138.and.mp.eq.78.and.ihem.eq.1) then
!          print *,'p1',lp,Ue(1),Ue(2),apex_d1d1,apex_d2d2,apex_d1d2
!          print *,'p1',apex_D,apex_BMAG
!         endif

      !      ni_oplus_1d(i1d)=plasma_3d(i,lp_plas,mp,1)
      !      ni_hplus_1d(i1d)=plasma_3d(i,lp_plas,mp,2)

      !      n_plus( i1d)=plasma_3d(i,lp_plas,mp,4)
      !      no_plus(i1d)=plasma_3d(i,lp_plas,mp,5)
      !      o2_plus(i1d)=plasma_3d(i,lp_plas,mp,6)
      !      n2_plus(i1d)=plasma_3d(i,lp_plas,mp,7)
      !      ti_oplus_1d(i1d)=plasma_3d(i,lp_plas,mp,11)

      !      tn(i1d)=TN_k(  i,lp_plas,mp) ![K]
      !      o(i1d)=ON_m3( i,lp_plas,mp) !m-3
      !      o2(i1d)=O2N_m3(i,lp_plas,mp)
      !      n2(i1d)=N2N_m3(i,lp_plas,mp)
      !      grid % apex_be3(lp_plas,mp)
             effective_temp = (neutrals % temperature(i,lp,mp)+plasma % ion_temperature(i,lp,mp))/2.0_prec
             IF ( effective_temp < neutrals % temperature(i,lp,mp) ) effective_temp = neutrals % temperature(i,lp,mp)

             CALL IONNEUT_PLAS(neutrals % oxygen(i,lp,mp),&
                               neutrals % molecular_oxygen(i,lp,mp),&
                               neutrals % molecular_nitrogen(i,lp,mp),&
                               plasma % ion_densities(1,i,lp,mp),&
                               plasma % ion_densities(5,i,lp,mp),&
                               plasma % ion_densities(6,i,lp,mp),&
                               effective_temp,&
                               ion_neutral_collisionfrequency,&
                               ion_mass_amu )

              !! calculates collision frequencies see TGCM
              o_cm3  = neutrals % oxygen(i,lp,mp)*1.e-6   ! convert from #/m3 to #/cm3
              o2_cm3 = neutrals % molecular_oxygen(i,lp,mp)*1.e-6
              n2_cm3 = neutrals % molecular_nitrogen(i,lp,mp)*1.e-6
              !!
              CALL calc_collfreq(o_cm3,                                   &
                                 o2_cm3,                                  &
                                 n2_cm3,                                  &
                                 effective_temp,                          &
                                 neutrals % temperature(i,lp,mp),         &
                                 grid % magnetic_field_strength(i,lp,mp), &
                                 rnu_o2p,rnu_op,rnu_nop,rnu_ne)

!! get pedersen & hall conductivities
!! more ion spieces can be added for calculating electron density
               electron_density = plasma % ion_densities(1,i,lp,mp) + plasma % ion_densities(5,i,lp,mp) + plasma % ion_densities(6,i,lp,mp)
!              electron_density = plasma % ion_densities(1,i,lp,mp) + plasma % ion_densities(5,i,lp,mp) + plasma % ion_densities(6,i,lp,mp)+ &
!                                 plasma % ion_densities(2,i,lp,mp) + plasma % ion_densities(3,i,lp,mp) + plasma % ion_densities(4,i,lp,mp)+ &
!                                 plasma % ion_densities(7,i,lp,mp)
               electron_charge_Coulombs=1.6022E-19
!
               if (electron_density.gt.1.e-10) then
!
               r_factor  = ion_mass_amu*ion_neutral_collisionfrequency/electron_charge_Coulombs/apex_BMAG
               qe_fac = electron_charge_Coulombs/apex_BMAG
!
!! densities #/m^3
!! Pedersen conductivity
        sigma_ped = qe_fac* &
    &         ((plasma % ion_densities(1,i,lp,mp)*rnu_op  /(1.+rnu_op **2))+ &
    &          (plasma % ion_densities(6,i,lp,mp) *rnu_o2p/(1.+rnu_o2p**2))+ &
    &          (plasma % ion_densities(5,i,lp,mp)*rnu_nop /(1.+rnu_nop**2))+ &
    &          (electron_density*rnu_ne  /(1.+rnu_ne**2)))

!! Hall conductivity
        sigma_hall = qe_fac* &
    &          (electron_density /(1.+rnu_ne **2)- &
    &           plasma % ion_densities(1,i,lp,mp)/(1.+rnu_op **2)- &
    &           plasma % ion_densities(6,i,lp,mp)/(1.+rnu_o2p**2)- &
    &           plasma % ion_densities(5,i,lp,mp)/(1.+rnu_nop**2))

! get integrals
                abs_ds = ABS(ds)

!         if (lp.gt.134.and.lp.le.139.and.mp.eq.78.and.ihem.eq.1) then
!           print *,i,lp,electron_density,plasma % ion_densities(1,i,lp,mp),plasma % ion_densities(5,i,lp,mp),plasma % ion_densities(6,i,lp,mp)
!           print *,i,lp,sigma_ped,sigma_hall
!         endif
!g
!g  The following integrals all come from page 203 and 204 of the paper.  They are numbered
!g  to match the equations in the paper...
!g  The integral parts are calculated here and then some additional factors are applied below.
!g  Note, however, we ignore the |sin I_m| factor wherever it appears since this is applied
!g  later on within the Dynamo solver...
!g

               integral513 = integral513 + sigma_ped*apex_d1d1*abs_ds/apex_D

               integral514 = integral514 + sigma_ped*apex_d2d2*abs_ds/apex_D

               integral517 = integral517 + sigma_hall*abs_ds

               integral518 = integral518 + sigma_ped*apex_d1d2*abs_ds/apex_D

               integral519 = integral519 +(sigma_ped*apex_d1d1*Ue(2)/apex_D &
     &               + (sigma_hall-sigma_ped*apex_d1d2/apex_D)*Ue(1))*abs_ds

               integral520 = integral520 + ((sigma_hall+sigma_ped*apex_d1d2 &
     &               /apex_D )*Ue(2)- sigma_ped*apex_d2d2*Ue(1)/apex_D)*abs_ds

! inputs to the dynamo solver
!               IF ( sw_3DJ==1 ) THEN
! calculation of Je1 and Je2, will be useful later
!eq (5.7)
!                 ed11 = Ed1_90(1,lp,mp) + Ue(ipts,2) * Apex_BE3
!                 ed21 = Ed2_90(1,lp,mp) - Ue(ipts,1) * Apex_BE3
!                 Je(ipts,1) = sigma_ped(ipts)*Apex_d1d1(ipts) * ed11 &
!                 & +         ( sigma_ped(ipts)*Apex_d1d2(ipts) -
!                               sigma_hall(ipts)*apex_D(ipts) ) * ed21
!eq (5.8)
!                 Je(ipts,2) =( sigma_ped(ipts)*Apex_d1d2(ipts) + sigma_hall(ipts)*apex_D(ipts) ) * ed11 &
!                      & +            sigma_ped(ipts)*Apex_d2d2(ipts) * ed21
!              END IF !( sw_3DJ==1 ) THEN

          end if !(electron_density.gt.1.e-10)

!g
!g  Integrals 5.13 and 5.14 do not !include the |sin I_m| or 1/|sin I_m| factors
!respectively
!g  as these are dealt with in the dynamo module....
!g
          !sigma_phph_dsi(ihem,mp,lp) = integral513   !(5.13) divided by |sin I_m |
          plasma % conductivities(1,ihem,lp,mp) = integral513


          !sigma_lmlm_msi(ihem,mp,lp) = integral514   !(5.14) multiplied by | sin I_m |
          plasma % conductivities(2,ihem,lp,mp) = integral514
!g
!g  Integrals 5.17 and 5.18.....
!g
          !sigma_h(ihem,mp,lp) = integral517       !(5.17)
          plasma % conductivities(3,ihem,lp,mp) = integral517


          !sigma_c(ihem,mp,lp) = integral518       !(5.18)
          plasma % conductivities(4,ihem,lp,mp) = integral518
!g
!g  integral 5.19 is multiplied by BE3.  However we have not multiplied by |sinI_m| because
!g  this is done within the dynamo module itself (I've said this enough yeh ?)
!g
!          Kdmph_dsi(ihem,mp,lp) =  Apex_BE3*integral519  !(5.19) divided by |sin I_m |
          plasma % conductivities(5,ihem,lp,mp) = grid % apex_be3(lp,mp)*integral519

!g
!g  The following is equation 5.20.  The integral is multiplied by BE3.  There is also a minus
!g  sign for the northern hemisphere part (and a plus sign for the southern hemisphere part).
!g  The +- statement is that the upper (lower) sign applies to the northern southern) magnetic
!g  hemisphere.  This statement is written on page 200 of the paper - just below equation 3.22
!g
          if(ihem == 2) then
            !Kdmlm(ihem,mp,lp) = Apex_BE3*integral520  !(5.20) plus for southern hemi
            plasma % conductivities(6,ihem,lp,mp) = grid % apex_be3(lp,mp)*integral520
          else if(ihem == 1) then
            !Kdmlm(ihem,mp,lp) = -Apex_BE3*integral520  !(5.20) minus for northern hemi
            plasma % conductivities(6,ihem,lp,mp) = -grid % apex_be3(lp,mp)*integral520
          end if

          ENDDO
        ENDDO
     ENDDO
   ENDDO

#ifdef HAVE_MPI

   DO mp = grid % mp_low , grid % mp_high
     DO lp = 1, grid % NLP
        DO ihem = 1, 2 
           conductivities(1:6,ihem,lp,mp) = plasma % conductivities(1:6,ihem,lp,mp)
        ENDDO
     ENDDO
   ENDDO

   ! This AllGather requires that NMP is evenly divisible by the number of MPI ranks
   sendcount = 6*2*grid % NLP*( grid % mp_high-grid % mp_low + 1 )
   CALL MPI_Allgather( conductivities(1:6,1:2,1:grid % NLP,grid % mp_low:grid % mp_high), &
                       sendcount,&
                       mpi_layer % mpi_prec,&
                       plasma % conductivities, &
                       sendcount, &
                       mpi_layer % mpi_prec, &
                       mpi_layer % mpi_communicator, &
                       mpierror)
    WRITE( 5000 + mpi_layer % rank_id,* ) plasma % conductivities
#endif

  END SUBROUTINE Calculate_Field_Line_Integrals


  SUBROUTINE calc_collfreq(o1_cm3,o2_cm3,n2_cm3,tnti,tn,apex_Bmag,rnu_o2p,rnu_op,rnu_nop,rnu_ne)
      IMPLICIT NONE
      REAL(prec), INTENT(in)::o1_cm3,o2_cm3,n2_cm3,tnti,tn,apex_Bmag
      REAL(prec), INTENT(out):: rnu_o2p,rnu_op,rnu_nop,rnu_ne

! local
      REAL(prec) :: te, sqrt_te, omega_op, omega_o2p, omega_nop, omega_op_inv, omega_o2p_inv, omega_nop_inv, omega_e, omega_e_inv

      REAL(kind=8) ::                                    &
           rnu_o2p_o2, & ! O2+ ~ O2 collision freq (resonant, temperature dependent)
           rnu_op_o2 , & ! O+  ~ O2 collision freq (non-resonant)
           rnu_nop_o2, & ! NO+ ~ O2 collision freq (non-resonant)
           rnu_o2p_o, &  ! O2+ ~ O  collision freq (non-resonant)
           rnu_op_o , &  ! O+  ~ O  collision freq (resonant, temperature dependent)
           rnu_nop_o, &  ! NO+ ~ O  collision freq (non-resonant)
           rnu_o2p_n2, & ! O2+ ~ N2 collision freq (non-resonant)
           rnu_op_n2 , & ! O+  ~ N2 collision freq (non-resonant)
           rnu_nop_n2    ! NO+ ~ N2 collision freq (non-resonant)

! gyrofrequencies: omega_i = eB/m_i  [1/s]
!                  omega_e = eB/m_e  [1/s]
! with qeoNao10 = e/Na [C/mol g/kg ]
!      qeomeo10 = e/m_e [C/g g/kg ]
! 1000 in qeoNao10 and qeomeo10 for conversion from 1/g to 1/kg
!
            omega_op     = qeoNao10*apex_Bmag*rmassinv_o1
            omega_o2p    = qeoNao10*apex_Bmag*rmassinv_o2
            omega_nop    = qeoNao10*apex_Bmag*rmassinv_nop
            omega_op_inv = 1./omega_op
            omega_o2p_inv= 1./omega_o2p
            omega_nop_inv= 1./omega_nop
            omega_e      = qeomeo10*apex_Bmag
            omega_e_inv  = 1./omega_e
!
! O2 collision frequencies:
           rnu_o2p_o2 = 2.59E-11*sqrt(tnti)*  &! O2+ ~ O2 (resonant)
            (1.-0.073*log10(tnti))**2
           rnu_op_o2  = 6.64E-10                  ! O+  ~ O2
           rnu_nop_o2 = 4.27E-10                  ! NO+ ~ O2
!
! O collision frequencies:
           rnu_o2p_o = 2.31E-10                   ! O2+ ~ O
           rnu_op_o  = 3.67e-11*sqrt(tnti)* &  ! O+  ~ O (resonant)
            (1.-0.064*log10(tnti))**2*colfac
           rnu_nop_o = 2.44E-10                   ! NO+ ~ O
!
! N2 collision frequencies:
           rnu_o2p_n2 = 4.13E-10                  ! O2+ ~ N2
           rnu_op_n2  = 6.82E-10                  ! O+  ~ N2
           rnu_nop_n2 = 4.34E-10                  ! NO+ ~ N2

! collision frequency nu_in for each ion [1/s]
!    by multiplying with neutral number density [1/cm^3] and sum over
!    neutrals
! nu_in is divided by gyrofrequency omega_i
! nu_in/omega_i [-]:
! rnu_o2p = [[o2p~o2]n(o2)+[o2p~o]n(o)+[o2p~n2]n(n2)]/w(o2p)
! rnu_op  = [[op ~o2]n(o2)+[op ~o]n(o)+[op ~n2]n(n2)]/w(op )
! rnu_nop = [[nop~o2]n(o2)+[nop~o]n(o)+[nop~n2]n(n2)]/w(nop)
           rnu_o2p = (rnu_o2p_o2*o2_cm3 + &
                           rnu_o2p_o *o1_cm3 + &
                           rnu_o2p_n2*n2_cm3)*omega_o2p_inv
           rnu_op  = (rnu_op_o2 *o2_cm3 + &
                           rnu_op_o  *o1_cm3 + &
                           rnu_op_n2 *n2_cm3)*omega_op_inv
           rnu_nop = (rnu_nop_o2*o2_cm3 + &
                           rnu_nop_o *o1_cm3 + &
                           rnu_nop_n2*n2_cm3)*omega_nop_inv
!
! neutral~electron collision frequency (from Banks & Kockards) nu_en
! divided by gyrofrequency omega_2:
! nu_en/omega_e [-]
!
            te = tn       ! Te is not available approxilate Te = Tn
            sqrt_te = sqrt(te)

            rnu_ne =  &
             (2.33e-11*n2_cm3*te     *(1.-1.21e-4*te     )+ &
              1.82e-10*o2_cm3*sqrt_te*(1.+3.60e-2*sqrt_te)+ &
              8.90e-11*o1_cm3*sqrt_te*(1.+5.70e-4*te    ))* &
              omega_e_inv
!
! 6/2/06 btf: Multiply rnu_ne by 4, as per Richmond:
!
! The effective electron-neutral collision frequency is increased in
! an an hoc manner by a factor of 4 in order for the model to produce
! electric fields and currents below 105 km that agree better with
! observations, as recommended by Gagnepain et al. (J. Atmos. Terr.
! Phys., 39, 1119-1124, 1977).
!
           rnu_ne = rnu_ne*4.


      END SUBROUTINE calc_collfreq
!
 SUBROUTINE IONNEUT_PLAS(P1,P2,P3,PI1,PI2,PI3,T,VIN,AMIn)
  IMPLICIT NONE
  REAL(prec), INTENT(in) :: P1, P2, P3, PI1, PI2, PI3, T
  REAL(prec), INTENT(out) :: AMIn, VIN
  ! Local
  REAL(prec) :: a(1:3), b(1:3), amu, factor, v1, v2, sumPI, summol, mi1, mi2, mi3

      mi1 = 16.0_prec
      mi2 = 30.0_prec
      mi3 = 32.0_prec
      a(1:3) = (/3.42E-11, 6.66E-10, 6.82E-10/)
      b(1:3) = (/2.44E-10, 4.28E-10, 4.34E-10/)
      amu = 1.66E-27
      factor=1.0

      summol = PI2 + PI3
      sumPI = PI1 + PI2 + PI3
      v2 = b(1)*P1 + b(2)*P2+ b(3)*P3
      v1 = a(3)*P3 + a(2)*P2 + a(1)*P1*factor*SQRT(T)*(1.08-0.139*log10(T)+4.51E-03*log10(T)**2)
      if(summol < 1.d-90) summol=0.0
      if(v1 < 1.d-90) v1=0.0
      if(v2 < 1.d-90) v2=0.0
      VIN = (v1*PI1+v2*summol)*1.E-06/sumPI
      AMIn = (PI1*mi1+PI2*mi2+PI3*mi3)*amu/sumPI

  END SUBROUTINE IONNEUT_PLAS
!
!
!

END MODULE IPE_Plasma_Class
