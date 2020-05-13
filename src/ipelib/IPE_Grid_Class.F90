#include "IPE_Macros.inc"

MODULE IPE_Grid_Class

USE IPE_Precision
USE IPE_Constants_Dictionary
USE IPE_Common_Routines
USE IPE_Model_Parameters_Class
USE IPE_MPI_Layer_Class
USE IPE_IO_Class

USE HDF5

IMPLICIT NONE

  TYPE IPE_Grid
    INTEGER :: nFluxTube, NLP, NMP, NPTS2D
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
      PROCEDURE :: Create => ipe_grid_build
      PROCEDURE :: Trash  => Trash_IPE_Grid

      PROCEDURE :: Initialize => Initialize_IPE_Grid

      ! -- Legacy I/O -- !
      PROCEDURE :: Read_IPE_Grid
      PROCEDURE :: Calculate_Grid_Attributes_From_Legacy_Input
      ! ---------------- !

      PROCEDURE :: ReadFile  => ipe_grid_read_file
      PROCEDURE :: WriteFile => ipe_grid_write_file

      PROCEDURE :: Interpolate_to_Geographic_Grid
      PROCEDURE :: Interpolate_2D_to_Geographic_Grid

      ! Functions
      PROCEDURE :: SinI

  END TYPE IPE_Grid

REAL(prec), PARAMETER, PRIVATE :: dlonm90km = 4.5_prec

CONTAINS

  subroutine ipe_grid_build( grid, mpl, params )
    class( IPE_Grid )                          :: grid
    class( IPE_MPI_Layer ),        intent(in)  :: mpl
    type ( IPE_Model_Parameters ), intent(in)  :: params

    ! -- begin

    call Build_IPE_Grid( grid,               &
                         params % nFluxTube, &
                         params % NLP,       &
                         params % NMP,       &
                         params % NPTS2D,    &
                         mpl % mp_low,       &
                         mpl % mp_high,      &
                         mpl % mp_halo_size )

  end subroutine ipe_grid_build

  SUBROUTINE Build_IPE_Grid( grid, nFluxTube, NLP, NMP, NPTS2D, mp_low, mp_high, halo )

    IMPLICIT NONE

    CLASS( IPE_Grid ), INTENT(out) :: grid
    INTEGER, INTENT(in)            :: nFluxTube
    INTEGER, INTENT(in)            :: NLP
    INTEGER, INTENT(in)            :: NMP
    INTEGER, INTENT(in)            :: NPTS2D
    INTEGER, INTENT(in)            :: mp_low, mp_high, halo

    ! Local
    INTEGER :: i

    grid % nFluxTube = nFluxTube
    grid % NLP       = NLP
    grid % NMP       = NMP
    grid % NPTS2D    = NPTS2D

    grid % mp_low    = mp_low
    grid % mp_high   = mp_high
    grid % mp_halo   = halo

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
              grid % altitude_geo(1:nheights_geo) )

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

  END SUBROUTINE Build_IPE_Grid


  SUBROUTINE Trash_IPE_Grid( grid )

    IMPLICIT NONE

    CLASS( IPE_Grid ), INTENT(inout) :: grid

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
                grid % altitude_geo )

  END SUBROUTINE Trash_IPE_Grid


  SUBROUTINE Initialize_IPE_Grid( grid, params, error )

    IMPLICIT NONE

    CLASS( IPE_Grid ),             INTENT(inout) :: grid
    TYPE ( IPE_Model_Parameters ), INTENT(in)    :: params
    INTEGER,                       INTENT(out)   :: error

    INTEGER :: im, in, is, lp, jtop

    error = 0

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
          error = -1
          RETURN
        ENDIF
      ELSE
        grid % northern_top_index(lp) = im
        grid % southern_top_index(lp) = im + 1
      ENDIF

    ENDDO

    grid % npts2d = params % npts2d

  END SUBROUTINE Initialize_IPE_Grid


  SUBROUTINE Read_IPE_Grid( grid, filename )

    IMPLICIT NONE

    CLASS( IPE_Grid ), INTENT(inout) :: grid
    CHARACTER(*), INTENT(in)         :: filename

    ! Local
    INTEGER    :: i, lp, mp , istop
    INTEGER    :: fUnit, ii, NMP, NLP, NPTS2D, MaxFluxTube
    INTEGER    :: jmin(1:grid % NMP,1:grid % NLP)
    INTEGER    :: jmax(1:grid % NMP,1:grid % NLP)
    REAL(prec) ::  Be3_all1(1:grid % NMP, 1:grid % NLP)
    REAL(prec) ::  Be3_all2(1:grid % NMP, 1:grid % NLP)
    REAL(prec), ALLOCATABLE ::  dum0(:,:)
    REAL(prec), ALLOCATABLE ::  dum1(:,:)
    REAL(prec), ALLOCATABLE ::  dum2(:,:)
    REAL(prec), ALLOCATABLE ::  dum3(:,:)
    REAL(prec), ALLOCATABLE ::  dum4(:,:,:)
    REAL(prec), ALLOCATABLE ::  dum5(:,:,:)
    REAL(prec), ALLOCATABLE ::  dum6(:,:,:)


    ! Place all of these dummies on the heap so that we don't wind up
    ! with a stack overflow. For typical values of NPTS2D, NLP, NMP,
    ! the "dum*" arrays will easily take up 100's of MB.

    ALLOCATE( dum0(1:grid % NPTS2D,1:grid % NMP), &
              dum1(1:grid % NPTS2D,1:grid % NMP), &
              dum2(1:grid % NPTS2D,1:grid % NMP), &
              dum3(1:grid % NPTS2D,1:grid % NMP), &
              dum4(1:3,1:grid % NPTS2D,1:grid % NMP), &
              dum5(1:3,1:grid % NPTS2D,1:grid % NMP), &
              dum6(1:3,1:grid % NPTS2D,1:grid % NMP) )

    OPEN( UNIT   = NewUnit(fUnit), &
          FILE   = TRIM(filename), &
          STATUS = 'OLD', &
          FORM   = 'FORMATTED', &
          ACTION = 'READ' )

    READ( fUnit, * ) jmin, jmax

    DO lp = 1, grid % NLP

      grid % flux_tube_max(lp) = jmax(1,lp) - jmin(1,lp) + 1

    ENDDO

    MaxFluxTube = maxval(grid % flux_tube_max)

    ! Just in case we set up the grid initially with the improper number of
    ! flux tube points, we'll reset the grid here.
    IF( MaxFluxTube /= grid % nFluxTube )THEN

      NLP = grid % NLP
      NMP = grid % NMP
      NPTS2D = grid % NPTS2D

      CALL grid % Trash( )
      CALL grid % Build( MaxFluxTube, NLP, NMP, NPTS2D, 1, NMP, 1 )

    ENDIF

   ! grid % flux_tube_max = grid % flux_tube_max - jmin(1,:) + 1
    DO lp = 1, grid % NLP

      grid % flux_tube_midpoint(lp) = 1 + ( grid % flux_tube_max(lp) - 1)/2

    ENDDO

    READ( fUnit, * ) dum0, dum1, dum2, dum3

      DO i = 1, grid % flux_tube_max(lp)

        ii = jmin(1,lp) + (i-1)
        grid % r_meter(i,lp) =  dum0(ii,1)
        grid % altitude(i,lp) = grid % r_meter(i,lp) - earth_radius

      ENDDO

    DO mp = 1, grid % NMP
      DO lp = 1, grid % NLP

        DO i = 1, grid % flux_tube_max(lp)

          ii = jmin(1,lp) + (i-1)
          grid % colatitude(i,lp,mp) = dum1(ii,mp)
          grid % longitude(i,lp,mp)  = dum2(ii,mp)
          grid % q_factor(i,lp,mp)   = dum3(ii,mp)

        ENDDO

      ENDDO
    ENDDO


    READ( fUnit, * ) dum0

    DO lp = 1, grid % NLP
      DO i = 1, grid % flux_tube_max(lp)

        ii = jmin(1,lp) + (i-1)
        grid % magnetic_colatitude(i,lp) = dum0(ii,1)

      ENDDO
    ENDDO

    READ( fUnit, * ) dum0, dum1

    DO mp = 1, grid % NMP
      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          ii = jmin(1,lp) + (i-1)
          grid % foot_point_distance(i,lp,mp)     = dum0(ii,mp)
          grid % magnetic_field_strength(i,lp,mp) = dum1(ii,mp)

        ENDDO
      ENDDO
    ENDDO


    READ( fUnit, * ) dum4, dum5, dum6

    DO mp = 1, grid % NMP
      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          ii = jmin(1,lp) + (i-1)
          grid % apex_d_vectors(1,1,i,lp,mp) = dum4(1,ii,mp) !D1
          grid % apex_d_vectors(2,1,i,lp,mp) = dum4(2,ii,mp)
          grid % apex_d_vectors(3,1,i,lp,mp) = dum4(3,ii,mp)

          grid % apex_d_vectors(1,2,i,lp,mp) = dum5(1,ii,mp) !D2
          grid % apex_d_vectors(2,2,i,lp,mp) = dum5(2,ii,mp)
          grid % apex_d_vectors(3,2,i,lp,mp) = dum5(3,ii,mp)

          grid % apex_d_vectors(1,3,i,lp,mp) = dum6(1,ii,mp) !D3
          grid % apex_d_vectors(2,3,i,lp,mp) = dum6(2,ii,mp)
          grid % apex_d_vectors(3,3,i,lp,mp) = dum6(3,ii,mp)

        ENDDO
      ENDDO
    ENDDO


    READ( fUnit, * ) dum4, dum5

    DO mp = 1, grid % NMP
      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          ii = jmin(1,lp) + (i-1)
          grid % apex_e_vectors(1,1,i,lp,mp) = dum4(1,ii,mp)
          grid % apex_e_vectors(2,1,i,lp,mp) = dum4(2,ii,mp)
          grid % apex_e_vectors(3,1,i,lp,mp) = dum4(3,ii,mp)

          grid % apex_e_vectors(1,2,i,lp,mp) = dum5(1,ii,mp)
          grid % apex_e_vectors(2,2,i,lp,mp) = dum5(2,ii,mp)
          grid % apex_e_vectors(3,2,i,lp,mp) = dum5(3,ii,mp)

        ENDDO
      ENDDO
    ENDDO

    READ( fUnit, * ) Be3_all1, Be3_all2

    DO mp=1, grid % NMP
      DO lp=1, grid % NLP

        grid % apex_be3(lp,mp) = Be3_all1(mp,lp)

      ENDDO
    ENDDO

    CLOSE( fUnit )

    DEALLOCATE( dum0,&
                dum1,&
                dum2,&
                dum3,&
                dum4,&
                dum5,&
                dum6 )

    CALL grid % Calculate_Grid_Attributes_From_Legacy_Input( )

  END SUBROUTINE Read_IPE_Grid


  subroutine ipe_grid_read_file(grid, io, filename, rc)
    class(IPE_Grid)               :: grid
    class(IPE_IO)                 :: io
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: rc

    ! -- local
    integer :: i, j, is, ie, js, je
    integer :: mp_beg, mp_end

    integer, dimension(3) :: fdims, gdims, moffset, mcount

    character(len=34) :: dset_name

    ! -- begin
    rc = -1

    ! -- setup data decomposition arrays
    fdims(1) = grid % nFluxTube
    fdims(2) = grid % NLP
    fdims(3) = grid % NMP


    mp_beg = max( 1,          grid % mp_low  - grid % mp_halo )
    mp_end = min( grid % NMP, grid % mp_high + grid % mp_halo )

    mcount(1:2)  = fdims(1:2)
    mcount(3)    = mp_end - mp_beg + 1

    moffset(1:2) = 0
    moffset(3)   = mp_beg - 1

    ! -- assume io has been initialized

    ! -- use same domain as parent
    call io % domain_set(fdims, moffset, mcount)
    if (io % error_check(msg="Unable to setup I/O data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- open grid file for reading
    call io % file_open(filename, "r")
    if (io % error_check(msg="Unable to open file "//filename, &
      file=__FILE__,line=__LINE__)) return

    ! -- read apex 3D variables

    ! -- geographic colatitude
    dset_name = "/apex_grid/colatitude"
    call io % dset_read(dset_name, grid % colatitude(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- geographic longitude
    dset_name = "/apex_grid/longitude"
    call io % dset_read(dset_name, grid % longitude(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- gravity
    dset_name = "/apex_grid/gravity"
    call io % dset_read(dset_name, grid % grx(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- foot point distance
    dset_name = "/apex_grid/foot_point_distance"
    call io % dset_read(dset_name, grid % foot_point_distance(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- magnetic field strength
    dset_name = "/apex_grid/magnetic_field_strength"
    call io % dset_read(dset_name, grid % magnetic_field_strength(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- Q value
    dset_name = "/apex_grid/q"
    call io % dset_read(dset_name, grid % q_factor(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- L magnitude
    dset_name = ""
    is = lbound(grid % l_magnitude, dim=1)
    ie = ubound(grid % l_magnitude, dim=1)
    js = lbound(grid % l_magnitude, dim=2)
    je = ubound(grid % l_magnitude, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/l_magnitude_",2i0)') i, j
        call io % dset_read(dset_name, grid % l_magnitude(i,j,:,:,mp_beg:mp_end))
        if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
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
        call io % dset_read(dset_name, grid % apex_e_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
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
        call io % dset_read(dset_name, grid % apex_d_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex d1xd2 mag
    dset_name = "/apex_grid/d1xd2_mag"
    call io % dset_read(dset_name, grid % d1xd2_magnitude(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- apex 2d variables
    ! -- reset domain
    call io % domain_set(fdims(2:3), moffset(2:3), mcount(2:3))
    if (io % error_check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/apex_be3"
    call io % dset_read(dset_name, grid % apex_be3(:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- reset domain
    call io % domain_set(fdims(1:2), moffset(1:2), mcount(1:2))
    if (io % error_check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/altitude"
    call io % dset_read(dset_name, grid % altitude)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/m_colat"
    call io % dset_read(dset_name, grid % magnetic_colatitude)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/r_meter"
    call io % dset_read(dset_name, grid % r_meter)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return


    ! -- apex 1d variables
    ! -- reset domain
    call io % domain_set(fdims(2:2), moffset(2:2), mcount(2:2))
    if (io % error_check(msg="Unable to setup 1D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/p_value"
    call io % dset_read(dset_name, grid % p_value)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_midpoint"
    call io % dset_read(dset_name, grid % flux_tube_midpoint)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_max"
    call io % dset_read(dset_name, grid % flux_tube_max)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/southern_top"
    call io % dset_read(dset_name, grid % southern_top_index)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/northern_top"
    call io % dset_read(dset_name, grid % northern_top_index)
    if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- periodic boundary conditions for B and q_factor if needed
    if (grid % mp_halo > 0) then
      ! -- reset domain to halo slice
      mcount(3)   = grid % mp_halo

      ! -- left boundary
      if (mp_beg == 1) then
        ! -- reset domain
        moffset(3) = grid % NMP - grid % mp_halo
        call io % domain_set(fdims, moffset, mcount)
        if (io % error_check(msg="Unable to apex I/O halo decomposition",file=__FILE__,line=__LINE__)) return

        ! -- magnetic field strength
        dset_name = "/apex_grid/magnetic_field_strength"
        call io % dset_read(dset_name, grid % magnetic_field_strength(:,:,mp_beg-grid % mp_halo:mp_beg-1))
        if (io % error_check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

        ! -- Q value
        dset_name = "/apex_grid/q"
        call io % dset_read(dset_name, grid % q_factor(:,:,mp_beg-grid % mp_halo:mp_beg-1))
        if (io % error_check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

      end if

      ! -- right boundary
      if (mp_end == grid % NMP) then
        ! -- reset domain
        moffset(3) = 0
        call io % domain_set(fdims, moffset, mcount)
        if (io % error_check(msg="Unable to apex I/O halo decomposition",file=__FILE__,line=__LINE__)) return

        ! -- magnetic field strength
        dset_name = "/apex_grid/magnetic_field_strength"
        call io % dset_read(dset_name, grid % magnetic_field_strength(:,:,mp_end+1:mp_end+grid % mp_halo))
        if (io % error_check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

        ! -- Q value
        dset_name = "/apex_grid/q"
        call io % dset_read(dset_name, grid % q_factor(:,:,mp_end+1:mp_end+grid % mp_halo))
        if (io % error_check(msg="Unable to read dataset", file=__FILE__,line=__LINE__)) return

      end if

    end if

    ! -- read geographic 3D variables
    fdims   = (/ grid % nheights_geo, grid % nlat_geo, grid % nlon_geo /)
    mcount  = fdims
    moffset = 0

    call io % domain_set(fdims, moffset, mcount)
    if (io % error_check(msg="Unable to setup I/O geographic data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- fac
    dset_name = ""
    is = lbound(grid % facfac_interface, dim=1)
    ie = ubound(grid % facfac_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/fac_",i0)') i
      call io % dset_read(dset_name, grid % facfac_interface(i,:,:,:))
      if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- dd
    dset_name = ""
    is = lbound(grid % dd_interface, dim=1)
    ie = ubound(grid % dd_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/dd_",i0)') i
      call io % dset_read(dset_name, grid % dd_interface(i,:,:,:))
      if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii1
    dset_name = ""
    is = lbound(grid % ii1_interface, dim=1)
    ie = ubound(grid % ii1_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii1_",i0)') i
      call io % dset_read(dset_name, grid % ii1_interface(i,:,:,:))
      if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii2
    dset_name = ""
    is = lbound(grid % ii2_interface, dim=1)
    ie = ubound(grid % ii2_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii2_",i0)') i
      call io % dset_read(dset_name, grid % ii2_interface(i,:,:,:))
      if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii3
    dset_name = ""
    is = lbound(grid % ii3_interface, dim=1)
    ie = ubound(grid % ii3_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii3_",i0)') i
      call io % dset_read(dset_name, grid % ii3_interface(i,:,:,:))
      if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii4
    dset_name = ""
    is = lbound(grid % ii4_interface, dim=1)
    ie = ubound(grid % ii4_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii4_",i0)') i
      call io % dset_read(dset_name, grid % ii4_interface(i,:,:,:))
      if (io % error_check(msg="Unable to read dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- close grid file
    call io % file_close()
    if (io % error_check(msg="Unable to close grid file "//filename, file=__FILE__,line=__LINE__)) return

    rc = 0

  end subroutine ipe_grid_read_file


  subroutine ipe_grid_write_file(grid, io, filename, rc)
    class(IPE_Grid)               :: grid
    class(IPE_IO)                 :: io
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: rc

    ! -- local
    integer :: i, j, is, ie, js, je
    integer :: mp_beg, mp_end

    integer, dimension(3) :: fdims, gdims, moffset, mcount

    character(len=34) :: dset_name

    ! -- begin
    rc = -1

    ! -- setup data decomposition arrays
    fdims(1) = grid % nFluxTube
    fdims(2) = grid % NLP
    fdims(3) = grid % NMP


    mp_beg = max( 1,          grid % mp_low  - grid % mp_halo )
    mp_end = min( grid % NMP, grid % mp_high + grid % mp_halo )

    mcount(1:2)  = fdims(1:2)
    mcount(3)    = mp_end - mp_beg + 1

    moffset(1:2) = 0
    moffset(3)   = mp_beg - 1

    ! -- assume io has been initialized

    ! -- use same domain as parent
    call io % domain_set(fdims, moffset, mcount)
    if (io % error_check(msg="Unable to setup I/O data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- open grid file for writing
    call io % file_open(filename, "c")
    if (io % error_check(msg="Unable to create file "//trim(filename), &
      file=__FILE__,line=__LINE__)) return

    ! -- create group
    call io % grp_build("apex_grid")
    if (io % error_check(msg="Error creating group /apex_grid",file=__FILE__,line=__LINE__)) return

    ! -- write apex 3D variables

    ! -- geographic colatitude
    dset_name = "/apex_grid/colatitude"
    call io % dset_write(dset_name, grid % colatitude(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name,file=__FILE__,line=__LINE__)) return

    ! -- geographic longitude
    dset_name = "/apex_grid/longitude"
    call io % dset_write(dset_name, grid % longitude(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name,file=__FILE__,line=__LINE__)) return

    ! -- gravity
    dset_name = "/apex_grid/gravity"
    call io % dset_write(dset_name, grid % grx(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- foot point distance
    dset_name = "/apex_grid/foot_point_distance"
    call io % dset_write(dset_name, grid % foot_point_distance(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- magnetic field strength
    dset_name = "/apex_grid/magnetic_field_strength"
    call io % dset_write(dset_name, grid % magnetic_field_strength(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- Q value
    dset_name = "/apex_grid/q"
    call io % dset_write(dset_name, grid % q_factor(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- L magnitude
    dset_name = ""
    is = lbound(grid % l_magnitude, dim=1)
    ie = ubound(grid % l_magnitude, dim=1)
    js = lbound(grid % l_magnitude, dim=2)
    je = ubound(grid % l_magnitude, dim=2)
    do j = js, je
      do i = is, ie
        write(dset_name, '("/apex_grid/l_magnitude_",2i0)') i, j
        call io % dset_write(dset_name, grid % l_magnitude(i,j,:,:,mp_beg:mp_end))
        if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
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
        call io % dset_write(dset_name, grid % apex_e_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
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
        call io % dset_write(dset_name, grid % apex_d_vectors(i,j,:,:,mp_beg:mp_end))
        if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
      end do
    end do

    ! -- apex d1xd2 mag
    dset_name = "/apex_grid/d1xd2_mag"
    call io % dset_write(dset_name, grid % d1xd2_magnitude(:,:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- apex 2d variables
    ! -- reset domain
    call io % domain_set(fdims(2:3), moffset(2:3), mcount(2:3))
    if (io % error_check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/apex_be3"
    call io % dset_write(dset_name, grid % apex_be3(:,mp_beg:mp_end))
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- the following datasets are available on all processses, so can be written just by one process
    ! -- to increase performance, pause output from all other processes
    call io % pause(mp_beg /= 1)

    ! -- apex 2d variables
    ! -- reset domain
    call io % domain_set(fdims(1:2), moffset(1:2), mcount(1:2))
    if (io % error_check(msg="Unable to setup 2D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/altitude"
    call io % dset_write(dset_name, grid % altitude)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/m_colat"
    call io % dset_write(dset_name, grid % magnetic_colatitude)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/r_meter"
    call io % dset_write(dset_name, grid % r_meter)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    ! -- apex 1d variables
    ! -- reset domain
    call io % domain_set(fdims(2:2), moffset(2:2), mcount(2:2))
    if (io % error_check(msg="Unable to setup 1D apex I/O decomposition",file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/p_value"
    call io % dset_write(dset_name, grid % p_value)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_midpoint"
    call io % dset_write(dset_name, grid % flux_tube_midpoint)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/tube_max"
    call io % dset_write(dset_name, grid % flux_tube_max)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/southern_top"
    call io % dset_write(dset_name, grid % southern_top_index)
    if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return

    dset_name = "/apex_grid/northern_top"
    call io % dset_write(dset_name, grid % northern_top_index)
    if (io % error_check(msg="Unable to write dataset", file=__FILE__,line=__LINE__)) return

    ! -- geographic 3D variables
    fdims   = (/ grid % nheights_geo, grid % nlat_geo, grid % nlon_geo /)
    mcount  = fdims
    moffset = 0

    call io % domain_set(fdims, moffset, mcount)
    if (io % error_check(msg="Unable to setup I/O geographic data decomposition",file=__FILE__,line=__LINE__)) return

    ! -- fac
    dset_name = ""
    is = lbound(grid % facfac_interface, dim=1)
    ie = ubound(grid % facfac_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/fac_",i0)') i
      call io % dset_write(dset_name, grid % facfac_interface(i,:,:,:))
      if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- dd
    dset_name = ""
    is = lbound(grid % dd_interface, dim=1)
    ie = ubound(grid % dd_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/dd_",i0)') i
      call io % dset_write(dset_name, grid % dd_interface(i,:,:,:))
      if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii1
    dset_name = ""
    is = lbound(grid % ii1_interface, dim=1)
    ie = ubound(grid % ii1_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii1_",i0)') i
      call io % dset_write(dset_name, grid % ii1_interface(i,:,:,:))
      if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii2
    dset_name = ""
    is = lbound(grid % ii2_interface, dim=1)
    ie = ubound(grid % ii2_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii2_",i0)') i
      call io % dset_write(dset_name, grid % ii2_interface(i,:,:,:))
      if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii3
    dset_name = ""
    is = lbound(grid % ii3_interface, dim=1)
    ie = ubound(grid % ii3_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii3_",i0)') i
      call io % dset_write(dset_name, grid % ii3_interface(i,:,:,:))
      if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- ii4
    dset_name = ""
    is = lbound(grid % ii4_interface, dim=1)
    ie = ubound(grid % ii4_interface, dim=1)
    do i = is, ie
      write(dset_name, '("/apex_grid/ii4_",i0)') i
      call io % dset_write(dset_name, grid % ii4_interface(i,:,:,:))
      if (io % error_check(msg="Unable to write dataset "//dset_name, file=__FILE__,line=__LINE__)) return
    end do

    ! -- restore output from all processes
    call io % pause(.false.)

    ! -- close grid file
    call io % file_close()
    if (io % error_check(msg="Unable to close grid file "//filename, file=__FILE__,line=__LINE__)) return

    rc = 0

  end subroutine ipe_grid_write_file


  SUBROUTINE Calculate_Grid_Attributes_From_Legacy_Input( grid )

    IMPLICIT NONE

    CLASS( IPE_Grid ) :: grid

    ! Local
    INTEGER    :: i, lp, mp, ii
    REAL(prec) :: a(1:3), b(1:3), c(1:3)
    REAL(prec) :: bhat(1:3), ufac

    DO lp = 1, grid % NLP

      ! "PValue" calculation
      grid % p_value(lp) = grid % r_meter(1,lp)/( earth_radius*(SIN(grid % magnetic_colatitude(1,lp) ))**2 )

      ! Calculate the northern top index
      DO i= 1, grid % flux_tube_midpoint(lp)

        grid % northern_top_index(lp) = i

        IF ( grid % altitude(i,lp) <= mesh_height_max .AND. grid % altitude(i+1,lp) > mesh_height_max )THEN
          EXIT
        ENDIF

      ENDDO

      ! Calculate the southern top index
      DO i = grid % flux_tube_max(lp), grid % flux_tube_midpoint(lp), -1

        grid % southern_top_index(lp) = i

        IF ( grid % altitude(i,lp) <= mesh_height_max .AND. grid % altitude(i-1,lp) > mesh_height_max )THEN
          EXIT
        ENDIF

      ENDDO

    ENDDO

    DO mp = 1, grid % NMP


      DO lp = 1, grid % NLP
        DO i = 1, grid % flux_tube_max(lp)

          IF( Almost_Equal( grid % apex_d_vectors(1,1,i,lp,mp), 0.0_prec ) .AND. &
              Almost_Equal( grid % apex_d_vectors(1,2,i,lp,mp), 0.0_prec ) )THEN

             ii = i-1

             IF( i /= grid % flux_tube_midpoint(lp) ) then

               PRINT*, ' '
               PRINT*, '  Module IPE_Grid_Class : S/R Calculate_Grid_Attributes_From_Legacy_Input'
               PRINT*, '    WARNING : Invalid apex_d_vectors'
               PRINT*, ' '

             ENDIF

          ELSE

            ii = i

          ENDIF

          ! calculate D from eq 3.15: | d1 X d2 |
          a(1) = grid % apex_d_vectors(1,1,ii,lp,mp)
          a(2) = grid % apex_d_vectors(2,1,ii,lp,mp)
          a(3) = grid % apex_d_vectors(3,1,ii,lp,mp)

          b(1) = grid % apex_d_vectors(1,2,ii,lp,mp)
          b(2) = grid % apex_d_vectors(2,2,ii,lp,mp)
          b(3) = grid % apex_d_vectors(3,2,ii,lp,mp)

          c(1) = a(2)*b(3) - a(3)*b(2)
          c(2) = a(3)*b(1) - a(1)*b(3)
          c(3) = a(1)*b(2) - a(2)*b(1)

          grid % d1xd2_magnitude(i,lp,mp) = SQRT ( c(1)**2 + c(2)**2 + c(3)**2 )

          bhat(1) = grid % apex_d_vectors(1,3,ii,lp,mp)*grid % d1xd2_magnitude(i,lp,mp)
          bhat(2) = grid % apex_d_vectors(2,3,ii,lp,mp)*grid % d1xd2_magnitude(i,lp,mp)
          bhat(3) = grid % apex_d_vectors(3,3,ii,lp,mp)*grid % d1xd2_magnitude(i,lp,mp)

          ufac  = SQRT(bhat(1)**2 + bhat(2)**2)

          !(1) magnetic eastward exactly horizontal
          !    l_e = bhat x k /|bhat x k|
          IF( ufac > 0.0_prec ) THEN
            grid % l_magnitude(1,1,i,lp,mp) =  bhat(2)/ufac
            grid % l_magnitude(2,1,i,lp,mp) = -bhat(1)/ufac
            grid % l_magnitude(3,1,i,lp,mp) =  0.0_prec
          ENDIF

          ! l_u = l_e x bhat
          grid % l_magnitude(1,2,i,lp,mp) = grid % l_magnitude(2,1,i,lp,mp)*bhat(3) - grid % l_magnitude(3,1,i,lp,mp)*bhat(2)
          grid % l_magnitude(2,2,i,lp,mp) = grid % l_magnitude(3,1,i,lp,mp)*bhat(1) - grid % l_magnitude(1,1,i,lp,mp)*bhat(3)
          grid % l_magnitude(3,2,i,lp,mp) = grid % l_magnitude(1,1,i,lp,mp)*bhat(2) - grid % l_magnitude(2,1,i,lp,mp)*bhat(1)

          grid % grx(i,lp,mp) = ( G0*earth_radius**2 )/( grid % r_meter(i,lp)**2 )!*grid % sinI(i,lp,mp)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Calculate_Grid_Attributes_From_Legacy_Input


  REAL(prec) FUNCTION SinI( grid, i, lp, mp )
    IMPLICIT NONE
    CLASS( IPE_Grid ) :: grid
    INTEGER           :: i, lp, mp
    ! Local
    REAL(prec) :: dotprod

    dotprod = SQRT( grid % apex_d_vectors(1,3,i,lp,mp)**2 + &
                    grid % apex_d_vectors(2,3,i,lp,mp)**2 + &
                    grid % apex_d_vectors(3,3,i,lp,mp)**2  )

    IF( Almost_Equal( dotprod, 0.0_prec )) THEN

      SinI = 0.0_prec

    ELSE

      SinI = grid % apex_d_vectors(3,3,i,lp,mp)/dotprod

    ENDIF

  END FUNCTION SinI


END MODULE IPE_Grid_Class
