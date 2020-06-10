MODULE IPE_Grid_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Common_Routines
  USE IPE_Model_Parameters_Class
  USE IPE_MPI_Layer_Class
  USE ipe_error_module

  USE COMIO

  IMPLICIT NONE

  TYPE IPE_Grid
    INTEGER :: nFluxTube, NLP, NMP
    INTEGER :: mp_low, mp_high, mp_halo
    INTEGER :: nheights_geo, nlat_geo, nlon_geo

    REAL(prec), ALLOCATABLE :: altitude(:,:)
    REAL(prec), ALLOCATABLE :: colatitude(:,:,:)
    REAL(prec), ALLOCATABLE :: longitude(:,:,:)
    REAL(prec), ALLOCATABLE :: grx(:,:,:)    ! GR index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: foot_point_distance(:,:,:)     ! ISL index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: magnetic_field_strength(:,:,:) ! IBM index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: magnetic_colatitude(:,:) ! plasma_grid_GL
    REAL(prec), ALLOCATABLE :: magnetic_longitude(:)
    REAL(prec), ALLOCATABLE :: r_meter(:,:) ! rmeter2D

    REAL(prec), ALLOCATABLE :: p_value(:) ! Pvalue
    REAL(prec), ALLOCATABLE :: q_factor(:,:,:) ! Q index of plasma_grid_3d
    REAL(prec), ALLOCATABLE :: l_magnitude(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: apex_d_vectors(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: d1xd2_magnitude(:,:,:)
    REAL(prec), ALLOCATABLE :: apex_e_vectors(:,:,:,:,:)
    REAL(prec), ALLOCATABLE :: apex_be3(:,:)

    INTEGER, ALLOCATABLE :: flux_tube_midpoint(:)
    INTEGER, ALLOCATABLE :: flux_tube_max(:)
    INTEGER, ALLOCATABLE :: southern_top_index(:)
    INTEGER, ALLOCATABLE :: northern_top_index(:)

    ! Geographic Interpolation attributes
    REAL(prec), ALLOCATABLE :: latitude_geo(:)
    REAL(prec), ALLOCATABLE :: longitude_geo(:)
    REAL(prec), ALLOCATABLE :: altitude_geo(:)
    REAL(prec), ALLOCATABLE :: facfac_interface(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dd_interface(:,:,:,:)
    INTEGER, ALLOCATABLE    :: ii1_interface(:,:,:,:)
    INTEGER, ALLOCATABLE    :: ii2_interface(:,:,:,:)
    INTEGER, ALLOCATABLE    :: ii3_interface(:,:,:,:)
    INTEGER, ALLOCATABLE    :: ii4_interface(:,:,:,:)

    CONTAINS

      PROCEDURE :: Build  => Build_IPE_Grid
      PROCEDURE :: Trash  => Trash_IPE_Grid

      PROCEDURE :: WriteFile => ipe_grid_write_file

      ! Functions
      PROCEDURE :: SinI

  END TYPE IPE_Grid

REAL(prec), PARAMETER, PRIVATE :: dlonm90km = 4.5_prec

CONTAINS

  SUBROUTINE Build_IPE_Grid( grid, io, mpl, params, filename, rc )
    CLASS(IPE_Grid)                              :: grid
    CLASS(COMIO_T)                               :: io
    CLASS( IPE_MPI_Layer ),        INTENT(inout) :: mpl
    TYPE ( IPE_Model_Parameters ), INTENT(in)    :: params
    CHARACTER(len=*),              INTENT(in)    :: filename
    INTEGER, OPTIONAL,             INTENT(out)   :: rc

    INTEGER :: localrc

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    CALL ipe_grid_read_file(grid, io, mpl, filename, localrc)
    IF ( ipe_error_check( localrc, &
      msg="Unable to read grid from file "//filename, &
      file=__FILE__,line=__LINE__)) RETURN

    CALL Initialize_IPE_Grid( grid, params, localrc )
    IF ( ipe_error_check( localrc, &
      msg="Unable to initialize grid", &
      file=__FILE__,line=__LINE__)) RETURN

  END SUBROUTINE Build_IPE_Grid


  SUBROUTINE Allocate_IPE_Grid( grid, rc )

    IMPLICIT NONE

    CLASS( IPE_Grid ), INTENT(inout) :: grid
    INTEGER,           INTENT(out)   :: rc

    ! Local
    INTEGER :: i, stat
    INTEGER :: nFluxTube, NLP, NMP
    INTEGER :: mp_low, mp_high, halo

    ! Begin
    rc = IPE_SUCCESS

    nFluxTube = grid % nFluxTube
    NLP       = grid % NLP
    NMP       = grid % NMP

    mp_low  = grid % mp_low
    mp_high = grid % mp_high
    halo    = grid % mp_halo

    grid % nheights_geo = nheights_geo
    grid % nlat_geo     = nlat_geo
    grid % nlon_geo     = nlon_geo

    ALLOCATE( grid % altitude(1:nFluxTube,1:NLP), &
              grid % colatitude(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % longitude(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % grx(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % foot_point_distance(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % magnetic_field_strength(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % magnetic_colatitude(1:nFluxTube,1:NLP), &
              grid % magnetic_longitude(mp_low-halo:mp_high+halo), &
              grid % r_meter(1:nFluxTube,1:NLP), &
              grid % p_value(1:NLP), &
              grid % q_factor(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % l_magnitude(1:3,1:2,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % apex_e_vectors(1:3,1:2,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % apex_d_vectors(1:3,1:3,1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % d1xd2_magnitude(1:nFluxTube,1:NLP,mp_low-halo:mp_high+halo), &
              grid % apex_be3(1:NLP,mp_low-halo:mp_high+halo), &
              grid % flux_tube_midpoint(1:NLP), &
              grid % flux_tube_max(1:NLP), &
              grid % southern_top_index(1:NLP), &
              grid % northern_top_index(1:NLP), &
              grid % facfac_interface(1:3,1:nheights_geo,1:nlat_geo,1:nlon_geo) ,&
              grid % dd_interface(1:3,1:nheights_geo,1:nlat_geo,1:nlon_geo) ,&
              grid % ii1_interface(1:3,1:nheights_geo,1:nlat_geo,1:nlon_geo) ,&
              grid % ii2_interface(1:3,1:nheights_geo,1:nlat_geo,1:nlon_geo) ,&
              grid % ii3_interface(1:3,1:nheights_geo,1:nlat_geo,1:nlon_geo) ,&
              grid % ii4_interface(1:3,1:nheights_geo,1:nlat_geo,1:nlon_geo), &
              grid % latitude_geo(1:nlat_geo), &
              grid % longitude_geo(1:nlon_geo), &
              grid % altitude_geo(1:nheights_geo), stat=stat )
    if (ipe_alloc_check(stat, msg="Failed to allocate grid internal arrays", &
      line=__LINE__, file=__FILE__, rc=rc)) return

    grid % altitude                = 0.0_prec
    grid % colatitude              = 0.0_prec
    grid % longitude               = 0.0_prec
    grid % grx                     = 0.0_prec
    grid % foot_point_distance     = 0.0_prec
    grid % magnetic_field_strength = 0.0_prec
    grid % magnetic_colatitude     = 0.0_prec
    grid % magnetic_longitude      = 0.0_prec
    grid % r_meter                 = 0.0_prec
    grid % p_value                 = 0.0_prec
    grid % q_factor                = 0.0_prec
    grid % l_magnitude             = 0.0_prec
    grid % apex_e_vectors          = 0.0_prec
    grid % apex_d_vectors          = 0.0_prec
    grid % d1xd2_magnitude         = 0.0_prec
    grid % apex_be3                = 0.0_prec

    grid % flux_tube_midpoint = 0
    grid % flux_tube_max      = 0
    grid % southern_top_index = 0
    grid % northern_top_index = 0

    grid % facfac_interface = 0.0_prec
    grid % dd_interface     = 0.0_prec
    grid % ii1_interface    = 0
    grid % ii2_interface    = 0
    grid % ii3_interface    = 0
    grid % ii4_interface    = 0


    DO i = mp_low-halo, mp_high+halo
      grid % magnetic_longitude(i) = REAL( (i-1),prec )*dlonm90km*dtr
    ENDDO

    DO i = 1, nlat_geo
      grid % latitude_geo(i) = -90.0_prec + REAL(i-1,prec)*180.0_prec/REAL( nlat_geo-1,prec )
    ENDDO

    DO i= 1, nlon_geo
      grid % longitude_geo(i) = REAL(i-1,prec)*360.0_prec/REAL( nlon_geo,prec )
    ENDDO

    DO i = 1, nheights_geo
      grid % altitude_geo(i) = REAL( (i-1)*5, prec ) + 90.0_prec
    ENDDO

  END SUBROUTINE Allocate_IPE_Grid


  SUBROUTINE Trash_IPE_Grid( grid, rc )

    IMPLICIT NONE

    CLASS( IPE_Grid ), INTENT(inout) :: grid
    INTEGER, OPTIONAL, INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    DEALLOCATE( grid % altitude, &
                grid % colatitude, &
                grid % longitude, &
                grid % grx, &
                grid % foot_point_distance, &
                grid % magnetic_field_strength, &
                grid % magnetic_colatitude, &
                grid % magnetic_longitude, &
                grid % r_meter, &
                grid % p_value, &
                grid % q_factor, &
                grid % l_magnitude, &
                grid % apex_e_vectors, &
                grid % apex_d_vectors, &
                grid % d1xd2_magnitude, &
                grid % apex_be3, &
                grid % flux_tube_midpoint, &
                grid % flux_tube_max, &
                grid % southern_top_index, &
                grid % northern_top_index, &
                grid % facfac_interface, &
                grid % dd_interface, &
                grid % ii1_interface, &
                grid % ii2_interface, &
                grid % ii3_interface, &
                grid % ii4_interface, &
                grid % longitude_geo, &
                grid % latitude_geo, &
                grid % altitude_geo, &
                stat=stat )
    IF ( ipe_dealloc_check(stat, msg="Failed to free up memory", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Trash_IPE_Grid


  SUBROUTINE Initialize_IPE_Grid( grid, params, rc )

    IMPLICIT NONE

    CLASS( IPE_Grid ),             INTENT(inout) :: grid
    TYPE ( IPE_Model_Parameters ), INTENT(in)    :: params
    INTEGER,                       INTENT(out)   :: rc

    INTEGER :: im, in, is, lp, jtop

    rc = IPE_SUCCESS

    ! -- compute northern, southern, and midpoint (apex) index along magnetic field line
    ! -- NOTE: field lines are assumed to be symmetrical around midpoint

    grid % northern_top_index = 0
    grid % southern_top_index = 0
    grid % flux_tube_midpoint = 0

    in = 1
    DO lp = 1, grid % NLP

      is = grid % flux_tube_max(lp)
      im = (is + in) / 2
      grid % flux_tube_midpoint(lp) = im

      IF ( grid % altitude(im,lp) > params % mesh_height_max ) THEN
        jtop = MAXLOC( grid % altitude(in:im,lp), DIM = 1, MASK = grid % altitude(in:im,lp) <= params % mesh_height_max )
        IF ( jtop > 0 ) THEN
          grid % northern_top_index(lp) = jtop + in - 1
          grid % southern_top_index(lp) = is - jtop + 1
        ELSE
          CALL ipe_error_set( msg="Unable to find grid point corresponding to mesh_height_max", &
            line=__LINE__, file=__FILE__, rc=rc )
          RETURN
        ENDIF
      ELSE
        grid % northern_top_index(lp) = im
        grid % southern_top_index(lp) = im + 1
      ENDIF

    ENDDO

  END SUBROUTINE Initialize_IPE_Grid


  subroutine ipe_grid_read_file(grid, io, mpl, filename, rc)
    class(IPE_Grid)                   :: grid
    class(COMIO_T)                    :: io
    class(IPE_MPI_Layer)              :: mpl
    character(len=*),     intent(in)  :: filename
    integer,              intent(out) :: rc

    ! -- local
    integer :: i, j, is, ie, js, je
    integer :: mp_beg, mp_end, n, rank

    integer, dimension(3) :: gdims, mstart, mcount
    integer, dimension(:), pointer :: fdims => null()

    character(len=34) :: dset_name

    ! -- begin
    rc = IPE_FAILURE

    ! -- open grid file for reading
    call io % open(filename, "r")
    if (io % err % check(msg="Unable to open file "//filename, &
      file=__FILE__,line=__LINE__)) return

    ! -- retrieve domain dimensions from longitude variable
    nullify(fdims)
    dset_name = "/apex_grid/longitude"
    call io % domain(dset_name, fdims)
    if (io % err % check(msg="Unable to retrieve domain dimensions from file "//filename, &
      file=__FILE__,line=__LINE__)) return

    ! -- domaain must be 3D
    if (io % err % check(size(fdims) /= 3, msg="Input grid domain must be 3D", &
      file=__FILE__,line=__LINE__)) return

    ! -- setup data decomposition arrays
    grid % nFluxTube = fdims(1)
    grid % NLP       = fdims(2)
    grid % NMP       = fdims(3)

    ! -- setup domain decomposition across MPI tasks
    call mpl % Set_Domain( grid % NLP, grid % NMP, rc )
    if (io % err % check(rc /= 0, &
      msg="Unable to setup domain decomposition over MPI", &
      file=__FILE__,line=__LINE__)) return

    grid % mp_low  = mpl % mp_low
    grid % mp_high = mpl % mp_high
    grid % mp_halo = mpl % mp_halo_size

    ! -- allocate internal arrays
    call Allocate_IPE_Grid( grid, rc )
    if (ipe_error_check(rc, msg="Unable to allocate mmeory for grid", &
      file=__FILE__,line=__LINE__)) return

    ! -- setup local domain
    mp_beg = max( 1,          grid % mp_low  - grid % mp_halo )
    mp_end = min( grid % NMP, grid % mp_high + grid % mp_halo )

    mcount(1:2) = fdims(1:2)
    mcount(3)   = mp_end - mp_beg + 1

    mstart(1:2) = 1
    mstart(3)   = mp_beg

    call io % domain(fdims, mstart, mcount)
    if (io % err % check(msg="Unable to setup I/O data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- read apex 3D variables

    ! -- geographic colatitude
    dset_name = "/apex_grid/colatitude"
    call io % read(dset_name, grid % colatitude(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- geographic longitude
    dset_name = "/apex_grid/longitude"
    call io % read(dset_name, grid % longitude(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- gravity
    dset_name = "/apex_grid/gravity"
    call io % read(dset_name, grid % grx(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- foot point distance
    dset_name = "/apex_grid/foot_point_distance"
    call io % read(dset_name, grid % foot_point_distance(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- magnetic field strength
    dset_name = "/apex_grid/magnetic_field_strength"
    call io % read(dset_name, grid % magnetic_field_strength(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- Q value
    dset_name = "/apex_grid/q"
    call io % read(dset_name, grid % q_factor(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- L magnitude
    dset_name = ""
    is = lbound(grid % l_magnitude, dim=1)
    ie = ubound(grid % l_magnitude, dim=1)
    js = lbound(grid % l_magnitude, dim=2)
    je = ubound(grid % l_magnitude, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/l_magnitude_",2i0)') i, j
        call io % read(dset_name, grid % l_magnitude(i,j,:,:,mp_beg:mp_end))
        if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex e
    dset_name = ""
    is = lbound(grid % apex_e_vectors, dim=1)
    ie = ubound(grid % apex_e_vectors, dim=1)
    js = lbound(grid % apex_e_vectors, dim=2)
    je = ubound(grid % apex_e_vectors, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/apex_e_",2i0)') i, j
        call io % read(dset_name, grid % apex_e_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex d
    dset_name = ""
    is = lbound(grid % apex_d_vectors, dim=1)
    ie = ubound(grid % apex_d_vectors, dim=1)
    js = lbound(grid % apex_d_vectors, dim=2)
    je = ubound(grid % apex_d_vectors, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/apex_d_",2i0)') i, j
        call io % read(dset_name, grid % apex_d_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex d1xd2 mag
    dset_name = "/apex_grid/d1xd2_mag"
    call io % read(dset_name, grid % d1xd2_magnitude(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- apex 2d variables
    ! -- reset domain
    call io % domain(fdims(2:3), mstart(2:3), mcount(2:3))
    if (io % err % check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/apex_be3"
    call io % read(dset_name, grid % apex_be3(:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- reset domain
    call io % domain(fdims(1:2), mstart(1:2), mcount(1:2))
    if (io % err % check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/altitude"
    call io % read(dset_name, grid % altitude)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/m_colat"
    call io % read(dset_name, grid % magnetic_colatitude)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/r_meter"
    call io % read(dset_name, grid % r_meter)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return


    ! -- apex 1d variables
    ! -- reset domain
    call io % domain(fdims(2:2), mstart(2:2), mcount(2:2))
    if (io % err % check(msg="Unable to setup 1D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/p_value"
    call io % read(dset_name, grid % p_value)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_midpoint"
    call io % read(dset_name, grid % flux_tube_midpoint)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_max"
    call io % read(dset_name, grid % flux_tube_max)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/southern_top"
    call io % read(dset_name, grid % southern_top_index)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/northern_top"
    call io % read(dset_name, grid % northern_top_index)
    if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- periodic boundary conditions for B and q_factor if needed
    if (grid % mp_halo > 0) then
      ! -- reset domain to halo slice
      mcount(3)   = grid % mp_halo

      ! -- left boundary
      if (mp_beg == 1) then
        ! -- reset domain
        mstart(3) = grid % NMP - grid % mp_halo + 1
        call io % domain(fdims, mstart, mcount)
        if (io % err % check(msg="Unable to apex I/O halo decomposition",file=__FILE__,line=__LINE__)) return

        ! -- magnetic field strength
        dset_name = "/apex_grid/magnetic_field_strength"
        call io % read(dset_name, grid % magnetic_field_strength(:,:,mp_beg-grid % mp_halo:mp_beg-1))
        if (io % err % check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

        ! -- Q value
        dset_name = "/apex_grid/q"
        call io % read(dset_name, grid % q_factor(:,:,mp_beg-grid % mp_halo:mp_beg-1))
        if (io % err % check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

      end if

      ! -- right boundary
      if (mp_end == grid % NMP) then
        ! -- reset domain
        mstart(3) = 1
        call io % domain(fdims, mstart, mcount)
        if (io % err % check(msg="Unable to apex I/O halo decomposition",file=__FILE__,line=__LINE__)) return

        ! -- magnetic field strength
        dset_name = "/apex_grid/magnetic_field_strength"
        call io % read(dset_name, grid % magnetic_field_strength(:,:,mp_end+1:mp_end+grid % mp_halo))
        if (io % err % check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

        ! -- Q value
        dset_name = "/apex_grid/q"
        call io % read(dset_name, grid % q_factor(:,:,mp_end+1:mp_end+grid % mp_halo))
        if (io % err % check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

      end if

    end if

    ! -- read geographic 3D variables
    fdims  = (/ grid % nheights_geo, grid % nlat_geo, grid % nlon_geo /)
    mcount = fdims
    mstart = 1

    call io % domain(fdims, mstart, mcount)
    if (io % err % check(msg="Unable to setup I/O geographic data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- fac
    dset_name = ""
    is = lbound(grid % facfac_interface, dim=1)
    ie = ubound(grid % facfac_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/fac_",i0)') i
      call io % read(dset_name, grid % facfac_interface(i,:,:,:))
      if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- dd
    dset_name = ""
    is = lbound(grid % dd_interface, dim=1)
    ie = ubound(grid % dd_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/dd_",i0)') i
      call io % read(dset_name, grid % dd_interface(i,:,:,:))
      if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii1
    dset_name = ""
    is = lbound(grid % ii1_interface, dim=1)
    ie = ubound(grid % ii1_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii1_",i0)') i
      call io % read(dset_name, grid % ii1_interface(i,:,:,:))
      if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii2
    dset_name = ""
    is = lbound(grid % ii2_interface, dim=1)
    ie = ubound(grid % ii2_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii2_",i0)') i
      call io % read(dset_name, grid % ii2_interface(i,:,:,:))
      if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii3
    dset_name = ""
    is = lbound(grid % ii3_interface, dim=1)
    ie = ubound(grid % ii3_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii3_",i0)') i
      call io % read(dset_name, grid % ii3_interface(i,:,:,:))
      if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii4
    dset_name = ""
    is = lbound(grid % ii4_interface, dim=1)
    ie = ubound(grid % ii4_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii4_",i0)') i
      call io % read(dset_name, grid % ii4_interface(i,:,:,:))
      if (io % err % check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- close grid file
    call io % close()
    if (io % err % check(msg="Unable to close grid file "//filename, file=__FILE__,line=__LINE__)) return

    rc = IPE_SUCCESS

  end subroutine ipe_grid_read_file


  subroutine ipe_grid_write_file(grid, io, filename, rc)
    class(IPE_Grid)               :: grid
    class(COMIO_T)                :: io
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: rc

    ! -- local
    integer :: i, j, is, ie, js, je
    integer :: mp_beg, mp_end

    integer, dimension(3) :: fdims, gdims, mstart, mcount

    character(len=34) :: dset_name

    ! -- begin
    rc = IPE_FAILURE

    ! -- setup data decomposition arrays
    fdims(1) = grid % nFluxTube
    fdims(2) = grid % NLP
    fdims(3) = grid % NMP


    mp_beg = max( 1,          grid % mp_low  - grid % mp_halo )
    mp_end = min( grid % NMP, grid % mp_high + grid % mp_halo )

    mcount(1:2)  = fdims(1:2)
    mcount(3)    = mp_end - mp_beg + 1

    mstart(1:2) = 1
    mstart(3)   = mp_beg

    ! -- assume io has been initialized

    ! -- use same domain as parent
    call io % domain(fdims, mstart, mcount)
    if (io % err % check(msg="Unable to setup I/O data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- open grid file for writing
    call io % open(filename, "c")
    if (io % err % check(msg="Unable to create file "//trim(filename), &
      file=__FILE__,line=__LINE__)) return

    ! -- write apex 3D variables

    ! -- geographic colatitude
    dset_name = "/apex_grid/colatitude"
    call io % write(dset_name, grid % colatitude(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name,file=__FILE__,line=__LINE__)) return

    ! -- geographic longitude
    dset_name = "/apex_grid/longitude"
    call io % write(dset_name, grid % longitude(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name,file=__FILE__,line=__LINE__)) return

    ! -- gravity
    dset_name = "/apex_grid/gravity"
    call io % write(dset_name, grid % grx(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- foot point distance
    dset_name = "/apex_grid/foot_point_distance"
    call io % write(dset_name, grid % foot_point_distance(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- magnetic field strength
    dset_name = "/apex_grid/magnetic_field_strength"
    call io % write(dset_name, grid % magnetic_field_strength(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- Q value
    dset_name = "/apex_grid/q"
    call io % write(dset_name, grid % q_factor(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- L magnitude
    dset_name = ""
    is = lbound(grid % l_magnitude, dim=1)
    ie = ubound(grid % l_magnitude, dim=1)
    js = lbound(grid % l_magnitude, dim=2)
    je = ubound(grid % l_magnitude, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/l_magnitude_",2i0)') i, j
        call io % write(dset_name, grid % l_magnitude(i,j,:,:,mp_beg:mp_end))
        if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex e
    dset_name = ""
    is = lbound(grid % apex_e_vectors, dim=1)
    ie = ubound(grid % apex_e_vectors, dim=1)
    js = lbound(grid % apex_e_vectors, dim=2)
    je = ubound(grid % apex_e_vectors, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/apex_e_",2i0)') i, j
        call io % write(dset_name, grid % apex_e_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex d
    dset_name = ""
    is = lbound(grid % apex_d_vectors, dim=1)
    ie = ubound(grid % apex_d_vectors, dim=1)
    js = lbound(grid % apex_d_vectors, dim=2)
    je = ubound(grid % apex_d_vectors, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/apex_d_",2i0)') i, j
        call io % write(dset_name, grid % apex_d_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex d1xd2 mag
    dset_name = "/apex_grid/d1xd2_mag"
    call io % write(dset_name, grid % d1xd2_magnitude(:,:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- apex 2d variables
    ! -- reset domain
    call io % domain(fdims(2:3), mstart(2:3), mcount(2:3))
    if (io % err % check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/apex_be3"
    call io % write(dset_name, grid % apex_be3(:,mp_beg:mp_end))
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- the following datasets are available on all processses, so can be written just by one process
    ! -- to increase performance, pause output from all other processes
    call io % pause(mp_beg /= 1)

    ! -- apex 2d variables
    ! -- reset domain
    call io % domain(fdims(1:2), mstart(1:2), mcount(1:2))
    if (io % err % check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/altitude"
    call io % write(dset_name, grid % altitude)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/m_colat"
    call io % write(dset_name, grid % magnetic_colatitude)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/r_meter"
    call io % write(dset_name, grid % r_meter)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- apex 1d variables
    ! -- reset domain
    call io % domain(fdims(2:2), mstart(2:2), mcount(2:2))
    if (io % err % check(msg="Unable to setup 1D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/p_value"
    call io % write(dset_name, grid % p_value)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_midpoint"
    call io % write(dset_name, grid % flux_tube_midpoint)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_max"
    call io % write(dset_name, grid % flux_tube_max)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/southern_top"
    call io % write(dset_name, grid % southern_top_index)
    if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/northern_top"
    call io % write(dset_name, grid % northern_top_index)
    if (io % err % check(msg="Unable to write dataset", file=__FILE__,line=__LINE__)) return

    ! -- geographic 3D variables
    fdims  = (/ grid % nheights_geo, grid % nlat_geo, grid % nlon_geo /)
    mcount = fdims
    mstart = 1

    call io % domain(fdims, mstart, mcount)
    if (io % err % check(msg="Unable to setup I/O geographic data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- fac
    dset_name = ""
    is = lbound(grid % facfac_interface, dim=1)
    ie = ubound(grid % facfac_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/fac_",i0)') i
      call io % write(dset_name, grid % facfac_interface(i,:,:,:))
      if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- dd
    dset_name = ""
    is = lbound(grid % dd_interface, dim=1)
    ie = ubound(grid % dd_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/dd_",i0)') i
      call io % write(dset_name, grid % dd_interface(i,:,:,:))
      if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii1
    dset_name = ""
    is = lbound(grid % ii1_interface, dim=1)
    ie = ubound(grid % ii1_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii1_",i0)') i
      call io % write(dset_name, grid % ii1_interface(i,:,:,:))
      if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii2
    dset_name = ""
    is = lbound(grid % ii2_interface, dim=1)
    ie = ubound(grid % ii2_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii2_",i0)') i
      call io % write(dset_name, grid % ii2_interface(i,:,:,:))
      if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii3
    dset_name = ""
    is = lbound(grid % ii3_interface, dim=1)
    ie = ubound(grid % ii3_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii3_",i0)') i
      call io % write(dset_name, grid % ii3_interface(i,:,:,:))
      if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii4
    dset_name = ""
    is = lbound(grid % ii4_interface, dim=1)
    ie = ubound(grid % ii4_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii4_",i0)') i
      call io % write(dset_name, grid % ii4_interface(i,:,:,:))
      if (io % err % check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- restore output from all processes
    call io % pause(.false.)

    ! -- close grid file
    call io % close()
    if (io % err % check(msg="Unable to close grid file "//filename, file=__FILE__,line=__LINE__)) return

    rc = IPE_SUCCESS

  end subroutine ipe_grid_write_file


  REAL(prec) FUNCTION SinI( grid, i, lp, mp )
    IMPLICIT NONE
    CLASS( IPE_Grid ) :: grid
    INTEGER           :: i, lp, mp
    ! Local
    REAL(prec) :: dotprod

    dotprod = SQRT( grid % apex_d_vectors(1,3,i,lp,mp) * grid % apex_d_vectors(1,3,i,lp,mp) + &
                    grid % apex_d_vectors(2,3,i,lp,mp) * grid % apex_d_vectors(2,3,i,lp,mp) + &
                    grid % apex_d_vectors(3,3,i,lp,mp) * grid % apex_d_vectors(3,3,i,lp,mp) )

    IF( Almost_Equal( dotprod, 0.0_prec )) THEN

      SinI = 0.0_prec

    ELSE

      SinI = grid % apex_d_vectors(3,3,i,lp,mp)/dotprod

    ENDIF

  END FUNCTION SinI


END MODULE IPE_Grid_Class
