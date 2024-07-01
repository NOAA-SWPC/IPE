MODULE IPE_Plasma_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Grid_Classw
  USE IPE_Neutrals_Class
  USE IPE_Forcing_Class
  USE IPE_Time_Class
  USE IPE_MPI_Layer_Class
  USE IPE_Common_Routines
  USE ipe_error_module

  USE HDF5

  IMPLICIT NONE

  TYPE IPE_Plasma
    INTEGER :: nFluxTube, NLP, NMP
    INTEGER :: mp_low, mp_high, mp_halo

    REAL(prec), POINTER     :: ion_densities(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities(:,:,:,:)
    REAL(prec), POINTER     :: ion_temperature(:,:,:)

    REAL(prec), POINTER     :: electron_density(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_density2(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity(:,:,:,:)
    REAL(prec), POINTER     :: electron_temperature(:,:,:)

    REAL(prec), ALLOCATABLE :: ion_densities_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_velocities_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: ion_temperature_old(:,:,:)

    REAL(prec), ALLOCATABLE :: electron_density_old(:,:,:)
    REAL(prec), ALLOCATABLE :: electron_velocity_old(:,:,:,:)
    REAL(prec), ALLOCATABLE :: electron_temperature_old(:,:,:)

    REAL(prec), ALLOCATABLE :: ionization_rates(:,:,:,:)
    REAL(prec), ALLOCATABLE :: conductivities(:,:,:,:)


    CONTAINS

      PROCEDURE :: Build => Build_IPE_Plasma
      PROCEDURE :: Trash => Trash_IPE_Plasma

      PROCEDURE :: Update => Update_IPE_Plasma
      PROCEDURE :: Update_Halos => Update_Halos_IPE_Plasma

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: Clean_Data
      PROCEDURE, PRIVATE :: Buffer_Old_State
      PROCEDURE, PRIVATE :: Calculate_Pole_Values
      PROCEDURE, PRIVATE :: Cross_Flux_Tube_Transport
      PROCEDURE, PRIVATE :: Test_Transport_Time_step
      PROCEDURE, PRIVATE :: Auroral_Precipitation
      PROCEDURE, PRIVATE :: FLIP_Wrapper
      PROCEDURE :: Calculate_Field_Line_Integrals

  END TYPE IPE_Plasma

  INTEGER, PARAMETER, PRIVATE    :: n_ion_species = 9
  REAL(prec), PARAMETER, PRIVATE :: safe_density_minimum = 1.0e+06_prec
  REAL(prec), PARAMETER, PRIVATE :: safe_temperature_minimum = 100.0_prec
!  REAL(prec), PARAMETER, PRIVATE :: colfac = 1.3_prec
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
!  INTEGER, PARAMETER , PRIVATE   :: transport_highlat_lp   = 30
!  INTEGER, PARAMETER , PRIVATE   :: perp_transport_max_lp  = 151

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
!  REAL(dp), PARAMETER, PRIVATE :: HPEQ     = 0.0D0
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

  SUBROUTINE Build_IPE_Plasma( plasma, nFluxTube, NLP, NMP, mp_low, mp_high, halo, rc )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(out) :: plasma
    INTEGER,             INTENT(in)  :: nFluxTube
    INTEGER,             INTENT(in)  :: NLP
    INTEGER,             INTENT(in)  :: NMP
    INTEGER,             INTENT(in)  :: mp_low, mp_high, halo
    INTEGER, OPTIONAL,   INTENT(out) :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

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
                plasma % ionization_rates(1:4,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
                stat = stat )
      IF ( ipe_alloc_check( stat, msg="Failed to allocate plasma internal arrays", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      plasma % ion_densities           = safe_density_minimum
      plasma % ion_velocities          = 0.0_prec
      plasma % ion_temperature         = safe_temperature_minimum
      plasma % electron_density        = safe_density_minimum
      plasma % electron_velocity       = 0.0_prec
      plasma % electron_temperature    = safe_temperature_minimum
      plasma % ionization_rates        = 0.0_prec
      plasma % conductivities          = 0.0_prec

#ifdef HAVE_MPI
      ALLOCATE( ion_requestHandle(1:16), ion_requestStats(MPI_STATUS_SIZE,1:16) )
      ion_requestHandle = 0
      ion_requestStats  = 0
#endif

  END SUBROUTINE Build_IPE_Plasma
!
  SUBROUTINE Trash_IPE_Plasma( plasma, rc )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    INTEGER, OPTIONAL,   INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS


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
                plasma % conductivities, &
                stat = stat )
    IF ( ipe_dealloc_check( stat, msg="Unable to free up memory", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

#ifdef HAVE_MPI
    DEALLOCATE( ion_requestHandle, ion_requestStats, stat=stat )
    IF ( ipe_dealloc_check( stat, msg="Unable to free up memory", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
#endif

  END SUBROUTINE Trash_IPE_Plasma
 

  SUBROUTINE Update_IPE_Plasma( plasma, grid, neutrals, forcing, time_tracker, mpi_layer, v_ExB, time_step, colfac, hpeq, &
                                transport_highlat_lp, perp_transport_max_lp, rc )
    IMPLICIT NONE
    CLASS( IPE_Plasma ),   INTENT(inout) :: plasma
    TYPE( IPE_Grid ),      INTENT(in)    :: grid
    TYPE( IPE_Neutrals ),  INTENT(in)    :: neutrals
    TYPE( IPE_Forcing ),   INTENT(in)    :: forcing
    TYPE( IPE_Time ),      INTENT(in)    :: time_tracker
    TYPE( IPE_MPI_Layer ), INTENT(in)    :: mpi_layer
    REAL(prec),            INTENT(in)    :: v_ExB(1:3,1:grid % NLP,grid % mp_low:grid % mp_high,2) ! am2023.04 add hemisphere
    REAL(prec),            INTENT(in)    :: time_step
    REAL(prec),            INTENT(in)    :: colfac
    REAL(prec),            INTENT(in)    :: hpeq
    INTEGER,               INTENT(in)    :: transport_highlat_lp
    INTEGER,               INTENT(in)    :: perp_transport_max_lp
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER    :: n_transport_timesteps
    INTEGER    :: i, lp, mp, j, localrc
    INTEGER    :: nflag_t(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    INTEGER    :: nflag_d(1 : plasma % NLP,plasma % mp_low : plasma % mp_high)
    REAL(prec) :: max_transport_convection_ratio_local
    REAL(prec) :: max_transport_convection_ratio
    REAL(prec) :: transport_time_step2
#ifdef HAVE_MPI
    INTEGER :: mpiError
#endif

      IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

      CALL plasma % Test_Transport_Time_step( grid, v_ExB, time_step, mpi_layer, &
                                              max_transport_convection_ratio_local, &
                                              perp_transport_max_lp )
					          
!7999 format('GHGM convect ratio ', f7.1, 3i4 , f12.1, 6e10.2)
!       write(6,*) 'transport_ratio ', mpi_layer % rank_id,max_transport_convection_ratio_local

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

!        write(1000 + mpi_layer % rank_id,*) ' GHGM TRANSPORT LOOP ', i , ' OF ', n_transport_timesteps
!         write(6,*) ' Inter-LP GHGM TRANSPORT LOOP ', i , ' OF ', n_transport_timesteps," ", mpi_layer % rank_id

        CALL plasma % Buffer_Old_State( grid )
        CALL plasma % Update_Halos( grid, mpi_layer )

#ifdef HAVE_MPI

        CALL MPI_WAITALL( 16, &
                         ion_requestHandle, &
                         ion_requestStats, &
                         mpiError)


#endif
! am2023.07 added the time_tracker for debugging
        CALL plasma % Cross_Flux_Tube_Transport( grid, v_ExB, transport_time_step2, time_tracker,&
                                                 transport_highlat_lp,perp_transport_max_lp, &
                                                 mpi_layer, localrc )
        IF ( ipe_error_check( localrc, msg="call to Cross_Flux_Tube_Transport failed", &
          line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      ENDDO

      CALL plasma % Auroral_Precipitation( grid, &
                                           neutrals, &
                                           forcing, &
                                           time_tracker )

       if (mpi_layer % rank_id.eq.0) then
       write(6,899) time_tracker % year, time_tracker % month, time_tracker % day, &       
                    time_tracker % hour, time_tracker % minute
 899   format('Calling Plasma         ', i4,x,i2.2,x,i2.2,2x,i2.2,':'i2.2)
       endif

      CALL plasma % FLIP_Wrapper( grid,         &
                                  neutrals,     &
                                  forcing,      &
                                  time_tracker, &
                                  time_step, colfac, hpeq, nflag_t,nflag_d )
       
!       write(6,890) time_tracker % year, time_tracker % month, time_tracker % day, &       
!                    time_tracker % hour, time_tracker % minute,mpi_layer % rank_id
! 890   format('Calling Field_line         ', i4,x,i2.2,x,i2.2,2x,i2.2,':'i2.2,x,i3)
 
      !TWFANG, calculate field line integrals for dynamo solver
      CALL plasma % Calculate_Field_Line_Integrals(grid, neutrals, colfac, mpi_layer)
      
!       write(6,880) time_tracker % year, time_tracker % month, time_tracker % day, &       
!                    time_tracker % hour, time_tracker % minute,mpi_layer % rank_id
! 880   format('After Field_line         ', i4,x,i2.2,x,i2.2,2x,i2.2,':'i2.22,x,i3)

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


  SUBROUTINE Test_Transport_Time_step( plasma, grid, v_ExB, time_step, mpi_layer, &
                                       max_transport_convection_ratio, perp_transport_max_lp )          
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    REAL(prec), INTENT(in)             :: v_ExB(1:3,1:grid % NLP, grid % mp_low:grid % mp_high,2) ! am2023.04 add hemisphere
    REAL(prec), INTENT(in)             :: time_step
    REAL(prec), INTENT(inout)          :: max_transport_convection_ratio                          ! am2023.06 was not defined as inout before
    TYPE( IPE_MPI_Layer ), INTENT(in)  :: mpi_layer
    INTEGER, INTENT(in)                :: perp_transport_max_lp
    ! Local
    
    REAL(prec) :: colat_90km(1:grid % NLP)
    REAL(prec) :: phi_t0 !magnetic longitude,phi[rad] at t0(previous time step)
    REAL(prec) :: theta_t0 !magnetic latitude,theta[rad] at t0
    REAL(prec) :: coslam, sinim
    INTEGER    :: mp, lp, i, lpx, mpx, jth, ih                       ! am2023.04 add hemisphere
    REAL(prec), PARAMETER :: rad_to_deg = 57.295779513
    REAL(prec) :: transport_convection_ratio(1:grid % NLP, grid % mp_low:grid % mp_high)
    REAL(prec) :: longitude_spacing, r, v_exb_max, exb_avg            ! am2023.06 exb_avg only for debugging
  
    transport_convection_ratio      = 0.0_prec
    max_transport_convection_ratio  = 0.0_prec

      colat_90km(1:grid % NLP) = grid % magnetic_colatitude(1,1:grid % NLP)
      r = earth_radius + 90000.0_prec
      longitude_spacing = 360.0_prec / REAL( plasma % NMP )
     
      DO 100 mp = plasma % mp_low, plasma % mp_high
        DO 200 lp = 1, perp_transport_max_lp
	
          v_exb_max = MAX(abs(v_ExB(1,lp,mp,1)),abs(v_ExB(1,lp,mp,2)))       ! am2023.06 take the max of the two hemispheres; test was done with avg ExB

          transport_convection_ratio(lp,mp) = (abs(v_ExB_max*time_step/(r*sin( colat_90km(lp)))) * rad_to_deg) / longitude_spacing ! am2023.06 changes for two hemispheres
  
 200    CONTINUE
 100  CONTINUE

      max_transport_convection_ratio = maxval(transport_convection_ratio)
      ! write(6,*) 'Ratio ', mpi_layer%rank_id,maxval(transport_convection_ratio),max_transport_convection_ratio

  END SUBROUTINE Test_Transport_Time_step


  SUBROUTINE Cross_Flux_Tube_Transport( plasma, grid, v_ExB, time_step, time_tracker, &
                                        transport_highlat_lp,perp_transport_max_lp, mpi_layer, rc )  
    IMPLICIT NONE
    CLASS( IPE_Plasma ),   INTENT(inout) :: plasma
    TYPE( IPE_Grid ),      INTENT(in)    :: grid
    TYPE( IPE_Time ),      INTENT(in)    :: time_tracker
    REAL(prec),            INTENT(in)    :: v_ExB(1:3,1:grid % NLP, grid % mp_low:grid % mp_high,2)  ! add hemisphere am2023.04
    REAL(prec),            INTENT(in)    :: time_step
    TYPE( IPE_MPI_Layer ), INTENT(in)    :: mpi_layer
    INTEGER,               INTENT(in)    :: transport_highlat_lp
    INTEGER,               INTENT(in)    :: perp_transport_max_lp
    INTEGER,               INTENT(out)   :: rc
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
    REAL(prec) :: z , z_factor,isign    ! am_2023.06  to add the other hemisphere
    INTEGER    :: lp_t0(1:2)
    INTEGER    :: mp_t0(1:2)
    INTEGER    :: lp_min, mp_min, isouth, inorth, ii, ispecial, ih, ibottom_flux_tube,iiq_start,iiq_end ! add hemisphere am2023.04
    INTEGER    :: mp, lp, i, lpx, mpx, jth
    real(prec) :: v1_exb_ns, v2_exb_ns   ! use the average of NH & SH ExB drift at the footpoints am2023.04
    INTEGER    :: i_min(1:2)
    LOGICAL    :: i_convection_too_far_in_lp
    REAL(prec), PARAMETER :: rad_to_deg = 57.295779513
    CHARACTER(len=128) :: errmsg

      rc = IPE_SUCCESS
 
      CALL plasma % Calculate_Pole_Values( grid,                       &
                                           mpi_layer,                  &
                                           ion_densities_pole_value,   &
                                           ion_temperature_pole_value, &
                                           ion_velocities_pole_value,  &
                                           electron_temperature_pole_value )


      colat_90km(1:grid % NLP) = grid % magnetic_colatitude(1,1:grid % NLP)
      r = earth_radius + 90000.0_prec

      DO 100 mp = plasma % mp_low, plasma % mp_high       ! magnetic longitudes
        DO 200 lp = 1, perp_transport_max_lp              ! flux fubes
	  do 400 ih=1,2                                   ! am_2023.06 since the dritf will be different in the two hemispheres loop over both NH, SH
	    ! am_2023.06.20 specify what part of the flux tube is convected 
	    
	    if(ih.eq.1) then    ! NH
	      ibottom_flux_tube = 1                         ! NH goes from footpoint to midpoint, increasing order
	    else
	      ibottom_flux_tube = grid % flux_tube_max(lp)  ! SH goes from footpoint to midpoint, decreasing order
	    endif

          i_convection_too_far_in_lp = .FALSE.
	  isign = 1.  
	  if(ih.eq.2.) isign = -1   ! am_2023.06  ih=1 NH and ih=2 SH (see subroutine Calculate_ExB_Velocity)
          ! 
          phi_t0   = grid % magnetic_longitude(mp) - v_ExB(1,lp,mp,ih)*time_step/(r*sin( colat_90km(lp) ) )  ! am_2023.06 - add hemisphere to ExB; sinus is sym about eq- no change necessary
          !v1_ExB_ns = 0.5*(v_ExB(1,lp,mp,1)+v_ExB(1,lp,mp,2))                                               ! am_2023.07 only for testing average exb drift
	  !phi_t0   = grid % magnetic_longitude(mp) - v1_ExB_ns*time_step/(r*sin( colat_90km(lp) ) )         ! am_2023.07 test 

          coslam   = cos( isign*(half_pi - grid % magnetic_colatitude(1,lp)) )                               ! am_2023.06 cos(latitude) -> sym about equator; therefore isign cosmetic
          sinim    = 2.0_prec*sqrt( 1.0_prec - coslam*coslam )/sqrt( 4.0_prec - 3.0_prec*coslam*coslam )
          
	  !
	  ! v2_ExB_ns = 0.5*(v_ExB(2,lp,mp,1)+v_ExB(2,lp,mp,2))         				     ! am_2023.07 only for testing average exb drift
	  ! theta_t0 = colat_90km(lp) - v2_ExB_ns*time_step/(r*sinim)   				     ! am_2023.07 test
           theta_t0 = colat_90km(lp) - v_ExB(2,lp,mp,ih)*time_step/(r*sinim)   				     ! am_2023.06 ve2 is down/equator. in both hemisphere; colat it will be different in 2 hemis which is ok
         								       				     ! both hemisphere can be treated the same wrt to finding the nearest points
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

            if((lp.ge.3).and.(lp.le.grid % NLP - 2)) then ! Make sure lp is within bounds                   

            ! Check poleward
              IF( theta_t0 <= colat_90km(lp-1) .AND. theta_t0 >= colat_90km(lp-2) )THEN
                lp_min = lp - 1
              ! Check equatorward
              ELSEIF( theta_t0 <= colat_90km(lp+2) .AND. theta_t0 >= colat_90km(lp+1) )THEN
                lp_min = lp + 2
              ENDIF

            endif ! Make sure lp is within bounds

          ENDIF

          IF( lp_min == 0 )THEN

            if((lp.ge.4).and.(lp.le.grid % NLP - 3)) then ! Make sure lp is within bounds                   

            IF( theta_t0 <= colat_90km(lp-2) .AND. theta_t0 >= colat_90km(lp-3) )THEN
              lp_min = lp - 2
            ! Check equatorward
            ELSEIF( theta_t0 <= colat_90km(lp+3) .AND. theta_t0 >= colat_90km(lp+2) )THEN
              lp_min = lp + 3
            ENDIF

            endif ! Make sure lp is within bounds

          ENDIF

          IF( lp_min == 0 )THEN

            if((lp.ge.5).and.(lp.le.grid % NLP - 4)) then ! Make sure lp is within bounds                   

            IF( theta_t0 <= colat_90km(lp-3) .AND. theta_t0 >= colat_90km(lp-4) )THEN
              lp_min = lp - 3
            ! Check equatorward
            ELSEIF( theta_t0 <= colat_90km(lp+4) .AND. theta_t0 >= colat_90km(lp+3) )THEN
              lp_min = lp + 4
            ENDIF

            endif ! Make sure lp is within bounds

          ENDIF

          IF( lp_min == 0 )THEN
            i_convection_too_far_in_lp = .TRUE.
            write(errmsg,*) 'GHGM convection too far ', mp , lp
            CALL ipe_warning_log( msg=errmsg, line=__LINE__, file=__FILE__ )
          ENDIF

! GHGM - check that lp_min is not greater than NLP
          
!         IF( lp_min.gt. grid % NLP )THEN
!           write(6,*) 'GHGM LP_MIN is greater than NLP - out of bounds ', mp , lp
!           lp_min = grid % NLP
!         ENDIF

! GHGM - check that lp_min is not greater than NLP
!         
!         IF( lp_min.lt. 2 )THEN
!           write(6,*) 'GHGM LP_MIN is less than 2 - out of bounds ', mp , lp
!           lp_min = 2
!         ENDIF

          if(i_convection_too_far_in_lp) then
! Very rare problem where convection in lp is greater than +-4 - for this case
! we set both lp indexes to just lp which means there will be no lp convection
! for this tube
            lp_t0(1) = lp
            lp_t0(2) = lp
          else
! The normal situation...
            lp_t0(1) = lp_min-1
            lp_t0(2) = lp_min
          endif

          IF( phi_t0 <= grid % magnetic_longitude(mp) .AND. phi_t0 >= grid % magnetic_longitude(mp-1) )THEN
            mp_min = mp
          ELSEIF( phi_t0 <= grid % magnetic_longitude(mp+1) .AND. phi_t0 >= grid % magnetic_longitude(mp) )THEN
            mp_min = mp+1
          ELSE
	     write(6,'("convect_loop mp_min not found",1(x,i4),4(x,f15.7))') mp,phi_t0,grid % magnetic_longitude(mp-1),grid % magnetic_longitude(mp),grid % magnetic_longitude(mp+1)
          ENDIF

          mp_t0(1) = mp_min-1
          mp_t0(2) = mp_min

          phi_i(1) = grid % magnetic_longitude(mp_t0(1))
          phi_i(2) = grid % magnetic_longitude(mp_t0(2))
          mp_comp_weight(1) =  ( phi_t0 - phi_i(2) )/( phi_i(1)-phi_i(2) )
          mp_comp_weight(2) = -( phi_t0 - phi_i(1) )/( phi_i(1)-phi_i(2) )

          IF( lp_min == 1 )THEN  ! lp_min == 1

            DO i =ibottom_flux_tube,grid % flux_tube_midpoint(lp),isign          ! am_2023.06.20 convect NH & SH flux tube part separately

              plasma % ion_densities(1:n_conv_spec,i,lp,mp)   = ion_densities_pole_value(1:n_conv_spec,i)
              plasma % ion_velocities(1:n_conv_spec,i,lp,mp)  = ion_velocities_pole_value(1:n_conv_spec,i)
              plasma % ion_temperature(i,lp,mp)               = ion_temperature_pole_value(i)
              plasma % electron_temperature(i,lp,mp)          = electron_temperature_pole_value(i)

            ENDDO

          ELSE ! lp_min =/= 1 ....
         
              if(lp_t0(2).eq.0) then
                lp_t0(1) = 1
                lp_t0(2) = 2
                !write(6,*) 'GHGM LP_T0 ',mp,lp,v1_ExB_ns, v2_ExB_ns
              endif

            lp_comp_weight(1) =  ( theta_t0 - colat_90km(lp_t0(2)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )
            lp_comp_weight(2) = -( theta_t0 - colat_90km(lp_t0(1)) )/( colat_90km(lp_t0(1))-colat_90km(lp_t0(2)) )
	 
            DO 300 i =ibottom_flux_tube,grid % flux_tube_midpoint(lp),isign   ! am_2023.06.20 convect NH & SH flux tube part separately

              ! q interpolation
              q_value = grid % q_factor(i,lp,mp) 
	      !if(mpi_layer % rank_id.eq.0) write(99,'("Inter-LPMP",4(x,i4),2(x,f15.7),3(x,e15.7))')lp,mp,ih,i, grid % magnetic_longitude(mp),grid % magnetic_colatitude(i,lp), &
	      !    grid % altitude(i,lp),plasma % ion_velocities_old(1,i,lp,mp),plasma % ion_velocities_old(2,i,lp,mp)
	      
              DO mpx = 1, 2
                DO lpx = 1,2
	         if(ih.eq.1) then       ! am_2023.06 for the ii loop which has to work for both hemisphere NH q goes from 1 (footpoint) to 0(apex), SH  q goes from 0 (apex) to -1 (footpoint)
	           iiq_start = ibottom_flux_tube+isign
	           iiq_end   = grid % flux_tube_midpoint(lp_t0(lpx))+1
	         else
	           iiq_start = grid % flux_tube_midpoint(lp_t0(lpx))-1
	           iiq_end   = grid % flux_tube_max(lp_t0(lpx))  ! ibottom_flux_tube+isign
	         end if	

                  B(lpx,mpx) = 0.0_prec
                  density(1:n_conv_spec,lpx,mpx)  = 0.0_prec
                  velocity(1:n_conv_spec,lpx,mpx) = 0.0_prec
                  temperature(lpx,mpx)   = 0.0_prec
                  e_temperature(lpx,mpx) = 0.0_prec

                  ! We assume by default the last flux tube point should be used
                  ! and we setup weights that prolong the last flux tube value
                  isouth = grid % flux_tube_max(lp_t0(lpx))
                  inorth = isouth
                  i_comp_weight(1) = 1.0_prec
                  i_comp_weight(2) = 0.0_prec
                  ! Search for the nearest q_factor
                  DO ii = iiq_start,iiq_end                                            ! am_2023.06.20 q is decreasing along the flux-tube NH (>0) to SH (<0), only for high lat flux tube from ~1 to -1
	                
		    IF(  grid % q_factor(ii, lp_t0(lpx), mp_t0(mpx)) <= q_value )THEN   ! this should also be fine for the SH... maybe I would not need the change of the loop above but do not need values of whole flux tube

                      isouth   = ii
                      inorth   = ii-1
                      q_int(1) = grid % q_factor(isouth, lp_t0(lpx), mp_t0(mpx))
                      q_int(2) = grid % q_factor(inorth, lp_t0(lpx), mp_t0(mpx))

                      i_comp_weight(1) =  ( q_value - q_int(2) )/( q_int(1) - q_int(2) )
                      i_comp_weight(2) = -( q_value - q_int(1) )/( q_int(1) - q_int(2) ) 
		              
                      EXIT  ! found a value
                    ENDIF
                  ENDDO
!
! GHGM - Below the bottom of the tube at the Northern end
! the Q interpolation factors are greater than 1 (incorrect).
! The following resets these values so the interpolation comes from the point
! at 90km
!
                  if((abs(i_comp_weight(1)).gt.1).and.(abs(i_comp_weight(2)).gt.1)) then
                      i_comp_weight(1) = 0.0
                      i_comp_weight(2) = 1.0
                  endif

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

                ENDDO    ! lpx loop
              ENDDO      ! mpx loop

              ! mp,lp interpolation
              ! Reduction over mpx, lpx
              ion_densities_int   = 0.0_prec
              ion_velocities_int  = 0.0_prec
              ion_temperature_int = 0.0_prec
              electron_temperature_int = 0.0_prec
              B_int               = 0.0_prec
              DO mpx = 1, 2
                DO lpx = 1, 2
                  B_int = B_int + B(lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)

                  ion_densities_int(1:n_conv_spec)  = ion_densities_int(1:n_conv_spec) + density(1:n_conv_spec,lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)
                  ion_velocities_int(1:n_conv_spec) = ion_velocities_int(1:n_conv_spec) + velocity(1:n_conv_spec,lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)

                  ion_temperature_int      = ion_temperature_int + temperature(lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)
                  electron_temperature_int = electron_temperature_int + e_temperature(lpx,mpx)*lp_comp_weight(lpx)*mp_comp_weight(mpx)
                ENDDO
              ENDDO

              IF( lp <= transport_highlat_lp )THEN
                ksi_fac = 1.0_prec
              ELSE
                ksi_fac = grid % magnetic_field_strength(i,lp,mp)/B_int
              ENDIF

              z = grid % altitude(i,lp)/1000.0_prec
              if ((z.lt.200.0).and.(z.gt.100.0)) then
              z_factor = (z - 100.0) / 100.0
              ion_densities_int(1:n_conv_spec) = (z_factor * (ion_densities_int(1:n_conv_spec) - plasma % ion_densities_old(1:n_conv_spec,i,lp,mp))) &
                                               + plasma % ion_densities_old(1:n_conv_spec,i,lp,mp)
              ion_velocities_int(1:n_conv_spec) = (z_factor * (ion_velocities_int(1:n_conv_spec) - plasma % ion_velocities_old(1:n_conv_spec,i,lp,mp))) &
                                               + plasma % ion_velocities_old(1:n_conv_spec,i,lp,mp)
              ion_temperature_int = (z_factor * (ion_temperature_int - plasma % ion_temperature_old(i,lp,mp))) &
                                               + plasma % ion_temperature_old(i,lp,mp)
              electron_temperature_int = (z_factor * (electron_temperature_int - plasma % electron_temperature_old(i,lp,mp))) &
                                               + plasma % electron_temperature_old(i,lp,mp)
              endif
              if (z.lt.100.0) then
                 ion_densities_int(1:n_conv_spec) = plasma % ion_densities_old(1:n_conv_spec,i,lp,mp)
                 ion_velocities_int(1:n_conv_spec) = plasma % ion_velocities_old(1:n_conv_spec,i,lp,mp)
                 ion_temperature_int = plasma % ion_temperature_old(i,lp,mp)
                 electron_temperature_int = plasma % electron_temperature_old(i,lp,mp)
              endif

              plasma % ion_densities(1:n_conv_spec,i,lp,mp)  = ion_densities_int(1:n_conv_spec)*( ksi_fac**2 )
              plasma % ion_velocities(1:n_conv_spec,i,lp,mp) = ion_velocities_int(1:n_conv_spec)
!             plasma % ion_temperature(i,lp,mp) = ion_temperature_int
!             plasma % electron_temperature(i,lp,mp) = electron_temperature_int
              plasma % ion_temperature(i,lp,mp) = ion_temperature_int*( ksi_fac**(4.0_prec/3.0_prec) )
              plasma % electron_temperature(i,lp,mp) = electron_temperature_int*( ksi_fac**(4.0_prec/3.0_prec) )
	      
 300        CONTINUE  !  i = ibottom_flux_tube,grid % flux_tube_midpoint(lp),isign

          ENDIF ! lp_min =/= 1

 400     CONTINUE  ! am_2023.06.23 loop over both hemispheres
 200    CONTINUE   ! lp = 1, perp_transport_max_lp 
 100  CONTINUE     ! mp = plasma % mp_low, plasma % mp_high   ! magnetic longitudes

  END SUBROUTINE Cross_Flux_Tube_Transport


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


  SUBROUTINE FLIP_Wrapper( plasma, grid, neutrals, forcing, time_tracker, flip_time_step, colfac, hpeq, nflag_t, nflag_d )
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout) :: plasma
    TYPE( IPE_Grid ), INTENT(in)       :: grid
    TYPE( IPE_Neutrals ), INTENT(in)   :: neutrals
    TYPE( IPE_Forcing ), INTENT(in)    :: forcing
    TYPE( IPE_Time ), INTENT(in)       :: time_tracker
    REAL(prec), INTENT(in)             :: flip_time_step
    REAL(prec), INTENT(in)             :: colfac
    REAL(prec), INTENT(in)             :: hpeq
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
    CHARACTER(len=95) :: ERRMSG
    CHARACTER(len=10) :: mp_lp_string
    LOGICAL, EXTERNAL :: CTIP_CHECK_EFLAG
    logical, parameter :: debug=.false.  ! debug one timestep am2023.05

      F107D = forcing % f107( forcing % current_index )
      F107A = forcing % f107_81day_avg( forcing % current_index )
      UTHR  = time_tracker % hour + (time_tracker % minute) / 60.0
      HPEQ_flip = HPEQ
      nflag_t = 0
      nflag_d = 0

      DO mp = plasma % mp_low, plasma % mp_high
        DO lp = 1, plasma % NLP
	  if(debug) write(6,'("In FLIP_Wrapper1st,2nd loop",5(x,i4))') mp,lp, plasma % mp_low, plasma % mp_high,plasma % NLP

          ! Copy over the grid information (for now)
          ZX(1:grid % flux_tube_max(lp))  = grid % altitude(1:grid % flux_tube_max(lp),lp)/1000.0_prec !convert from m to km
          PCO                             = grid % p_value(lp)  !Pvalue is a single value
          SLX(1:grid % flux_tube_max(lp)) = grid % foot_point_distance(1:grid % flux_tube_max(lp),lp,mp)
          GLX(1:grid % flux_tube_max(lp)) = half_pi - grid % magnetic_colatitude(1:grid % flux_tube_max(lp),lp)  ! magnetic latitude [radians]
          BMX(1:grid % flux_tube_max(lp)) = grid % magnetic_field_strength(1:grid % flux_tube_max(lp),lp,mp)   !Tesla
! Need a -SinI factor for the GRX
          DO i=1, grid % flux_tube_max(lp)
	    if(debug) write(6,*) 'In FLIP_Wrapper 3rdloop', mp,lp,i
	    
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
	    if(debug) write(6,'("In FLIP_Wrapper 4rdloop",6(x,i4))') mp,lp,i, plasma % mp_low, plasma % mp_high,plasma % NLP,grid % flux_tube_max(lp)

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
	    
	    if(debug) write(6,'("In FLIP_Wrapper before CTIPINIT ",6(x,i4))') mp,lp,i, plasma % mp_low, plasma % mp_high,plasma % NLP,grid % flux_tube_max(lp)

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
                        COLFAC, & !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
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
	 

          IF ( CTIP_CHECK_EFLAG( ERRMSG, EFLAG ) ) THEN
            write(mp_lp_string,"(2i4)") mp,lp
            CALL ipe_warning_log( msg=trim(ERRMSG)//trim(mp_lp_string), line=__LINE__, file=__FILE__ )
          ENDIF

          DO i=1, grid % flux_tube_max(lp)

	    if(debug) write(6,*) 'In FLIP_Wrapper 5th loop', mp,lp,i
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


  SUBROUTINE Calculate_Field_Line_Integrals( plasma, grid, neutrals, colfac, mpi_layer )


    ! this routine takes in the plasma, grid, and neutral fields and returns the six conductivities on IPE
    ! mp,lp grid. the output of this routine will likely need to be interpolated onto the dynamo solver grid
    ! before calling the dynamo solver
    IMPLICIT NONE
    CLASS( IPE_Plasma ), INTENT(inout)  :: plasma
    TYPE( IPE_Grid ), INTENT(in)        :: grid
    TYPE( IPE_Neutrals ), INTENT(in)    :: neutrals
    REAL(prec),            INTENT(in)   :: colfac
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
                                 rnu_o2p,rnu_op,rnu_nop,rnu_ne,colfac)

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
#endif

  END SUBROUTINE Calculate_Field_Line_Integrals


  SUBROUTINE calc_collfreq(o1_cm3,o2_cm3,n2_cm3,tnti,tn,apex_Bmag,rnu_o2p,rnu_op,rnu_nop,rnu_ne,colfac)
      IMPLICIT NONE
      REAL(prec), INTENT(in)::o1_cm3,o2_cm3,n2_cm3,tnti,tn,apex_Bmag,colfac
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
