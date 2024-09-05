MODULE IPE_MPI_Layer_Class

#ifdef HAVE_MPI
  USE mpi
#endif

  USE IPE_Precision

  IMPLICIT NONE

  ! The Cross_Flux_Tube_Transport is the only routine that will introduce data
  ! dependency between different ranks. If we parallelize in the MP-direction,
  ! then only "magnetic zonal" advection results in data necessarily accessed
  ! from another rank. Data that will need to be accessed from another rank
  ! includes
  !   * plasma % ion_densities_old
  !   * plasma % ion_temperature_old
  !   * grid % magnetic_field_strength
  !   * grid % q_factor
  !
  ! Each rank will have arrays (mp_low-mp_halo_size:mp_high+mp_halo_size)
  !
  TYPE IPE_MPI_Layer
    LOGICAL :: enabled
    LOGICAL :: initialized
    INTEGER :: global_nmp
    INTEGER :: nmp
    INTEGER :: lp_low, lp_high ! This processes nlp lower and upper bounds
    INTEGER :: mp_low, mp_high ! This processes nmp lower and upper bounds
    INTEGER :: mp_halo_size       ! how many mp grid points to the right and left should I keep updated from my neighbors
    integer :: mpi_prec
    integer :: mpi_communicator
    integer :: mpi_info
    INTEGER :: n_ranks
    INTEGER :: rank_id
    INTEGER :: neighbor_rank(1:2)

    CONTAINS

    PROCEDURE :: Initialize => Initialize_MPI_Layer
    PROCEDURE :: Finalize   => Finalize_MPI_Layer
    PROCEDURE :: Set_Domain => Set_Domain_On_MPI_Layer

  END TYPE IPE_MPI_Layer


CONTAINS

  SUBROUTINE Initialize_MPI_Layer( mpi_layer, comm )

    CLASS( IPE_MPI_Layer ), INTENT(inout) :: mpi_layer
    integer, OPTIONAL,      INTENT(in)    :: comm

    ! Local
    INTEGER :: mpiErr

    ! If we are coupled to another model, we are expecting that the MPI
    ! communicator has been assigned via the mediator, and we should
    ! not call the MPI_INIT routine. For IPE standalone, MPI_INIT needs
    ! to be called. The call to MPI_INIT occurs when the CPPFLAG
    ! -DCOUPLED is passed at compile time
    !
    ! Note that the mpi_communicator attribute should be filled in before
    ! calling this routine when running in COUPLED mode. For this reason,
    ! the "mpi_layer" structure has intent(inout)

    mpi_layer % enabled          = .false.
    mpi_layer % initialized      = .false.
    mpi_layer % global_nmp       = 0
    mpi_layer % nmp              = 0
    mpi_layer % lp_low           = 0
    mpi_layer % lp_high          = 0
    mpi_layer % mp_low           = 0
    mpi_layer % mp_high          = 0
    mpi_layer % mp_halo_size     = 0
!    mpi_layer % mpi_prec         = 0
!    mpi_layer % mpi_communicator = 0
!    mpi_layer % mpi_info         = 0
    mpi_layer % n_ranks          = 0
    mpi_layer % rank_id          = 0
    mpi_layer % neighbor_rank    = 0

#ifdef HAVE_MPI
    IF( prec == dp )THEN
      mpi_layer % mpi_prec = MPI_DOUBLE
    ELSE
      mpi_layer % mpi_prec = MPI_FLOAT
    ENDIF

    ! Check to see if MPI has been initialized
    CALL MPI_INITIALIZED( mpi_layer % initialized, mpiErr )

    IF( mpi_layer % initialized )THEN
      IF (PRESENT(comm)) mpi_layer % mpi_communicator = comm
    ELSE
      CALL MPI_INIT( mpiErr )
      mpi_layer % mpi_communicator = MPI_COMM_WORLD
    ENDIF

    CALL MPI_COMM_GET_INFO(mpi_layer % mpi_communicator, mpi_layer % mpi_info, mpiErr)

    CALL MPI_COMM_SIZE( mpi_layer % mpi_communicator, &
                        mpi_layer % n_ranks, &
                        mpiErr )

    CALL MPI_COMM_RANK( mpi_layer % mpi_communicator, &
                        mpi_layer % rank_id, &
                        mpiErr )

    mpi_layer % enabled = .true.
#else
    mpi_layer % rank_id = 0
    mpi_layer % n_ranks = 1
#endif

  END SUBROUTINE Initialize_MPI_Layer


  SUBROUTINE Set_Domain_On_MPI_Layer( mpi_layer, NLP, NMP, error )

    CLASS( IPE_MPI_Layer ), INTENT(inout) :: mpi_layer
    INTEGER,                INTENT(in)    :: NLP
    INTEGER,                INTENT(in)    :: NMP
    INTEGER,                INTENT(out)   :: error

    ! Local
    INTEGER :: remainder
    INTEGER :: mpiErr
    LOGICAL :: initialized

    error = 0

    mpi_layer % lp_low       = 1
    mpi_layer % lp_high      = NLP
    mpi_layer % global_nmp   = NMP
    mpi_layer % mp_halo_size = 1

    ! Calculate how many slices each rank will have
    mpi_layer % nmp = mpi_layer % global_nmp/(mpi_layer % n_ranks)
    remainder = MOD( mpi_layer % global_nmp, mpi_layer % n_ranks )

    IF( remainder > 0 )THEN
      PRINT *, 'NMP(80) must be evenly divisible by number of MPI ranks'
      FLUSH(6)
      error = -1
      RETURN
    ENDIF

    IF( mpi_layer % rank_id < remainder )THEN

      mpi_layer % nmp       = mpi_layer % nmp + 1
      mpi_layer % mp_low    = mpi_layer % rank_id*(mpi_layer % nmp) + 1
      mpi_layer % mp_high   = mpi_layer % mp_low + mpi_layer % nmp - 1

    ELSE

      mpi_layer % mp_low  = mpi_layer % rank_id*(mpi_layer % nmp) + 1 + remainder
      mpi_layer % mp_high = mpi_layer % mp_low + mpi_layer % nmp - 1

    ENDIF

    mpi_layer % neighbor_rank(1) = mpi_layer % rank_id - 1 ! Rank to the "left"
    mpi_layer % neighbor_rank(2) = mpi_layer % rank_id + 1 ! Rank to the "right"

    ! Set up the neighbor point-to-point communication patterns
    IF( mpi_layer % neighbor_rank(1) == -1 )THEN
      mpi_layer % neighbor_rank(1) = mpi_layer % n_ranks-1
    ENDIF

    IF( mpi_layer % neighbor_rank(2) == mpi_layer % n_ranks )THEN
      mpi_layer % neighbor_rank(2) = 0
    ENDIF

  END SUBROUTINE Set_Domain_On_MPI_Layer


  SUBROUTINE Finalize_MPI_Layer( mpi_layer )

    CLASS( IPE_MPI_Layer ), INTENT(inout) :: mpi_layer

    INTEGER :: mpiErr
#ifdef HAVE_MPI
    IF ( .NOT. mpi_layer % initialized ) CALL MPI_FINALIZE( mpiErr )
#endif

  END SUBROUTINE Finalize_MPI_Layer

END MODULE IPE_MPI_Layer_Class
