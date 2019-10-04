#include "IPE_Macros.inc"

MODULE IPE_HDF5

USE IPE_Precision
USE HDF5

IMPLICIT NONE

  INTERFACE Get_Variable_from_HDF5
    MODULE PROCEDURE MPI_Get_3D_Variable_from_HDF5, MPI_Get_2D_Variable_from_HDF5, Get_3D_Grid_Variable_from_HDF5, Get_2D_Grid_Variable_from_HDF5, &
                     Get_1D_Grid_Variable_from_HDF5, Get_Integer_3D_Grid_Variable_from_HDF5, Get_Integer_1D_Grid_Variable_from_HDF5
  END INTERFACE Get_Variable_from_HDF5
  INTERFACE Add_Variable_to_HDF5
    MODULE PROCEDURE MPI_Add_3D_Grid_Variable_to_HDF5, MPI_Add_2D_Grid_Variable_to_HDF5, Add_3D_Grid_Variable_to_HDF5, Add_2D_Grid_Variable_to_HDF5, &
                     Add_1D_Grid_Variable_To_HDF5, Add_Integer_3D_Grid_Variable_to_HDF5, Add_Integer_1D_Grid_Variable_to_HDF5
  END INTERFACE Add_Variable_to_HDF5
contains
  SUBROUTINE MPI_Get_3D_Variable_from_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, n_fluxtube, NLP, NMP, mp_low, mp_high, error )

    INTEGER(HID_T), INTENT(in)    :: file_id
    CHARACTER(*), INTENT(in)      :: variable_name
    INTEGER(HID_T), INTENT(in)    :: plist_id
    INTEGER(HID_T), INTENT(out)   :: filespace
    INTEGER(HID_T), INTENT(in)    :: memspace
    INTEGER, INTENT(in)           :: n_fluxtube, NLP, NMP
    INTEGER, INTENT(in)           :: mp_low, mp_high
    REAL(prec), INTENT(out)       :: variable(1:n_fluxtube,1:NLP,mp_low:mp_high)
    INTEGER(HSIZE_T), INTENT(in)  :: dimensions(1:3)
    INTEGER, INTENT(out)          :: error

    ! Local
    INTEGER(HID_T) :: dataset_id
    INTEGER        :: mp
    INTEGER(HSIZE_T)        :: starts(1:3), counts(1:3), strides(1:3)

    CALL h5dopen_f( file_id, TRIM(variable_name), dataset_id, error )
    IF ( error /= 0 ) RETURN
#ifdef HAVE_MPI

! -- PARALLEL -- !
    CALL h5dget_space_f( dataset_id, filespace, error )
    IF ( error /= 0 ) RETURN

    DO mp = mp_low, mp_high
      starts = (/ 0, 0, mp-1 /)
      counts = (/ 1, 1, 1 /)
      strides = (/ 1, 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )
      IF ( error /= 0 ) RETURN

      CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                      variable(1:n_fluxtube,1:NLP,mp), &
                      dimensions, error, memspace, filespace )
      IF ( error /= 0 ) RETURN

    ENDDO

#else

! -- SERIAL -- !

    CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                    variable, dimensions, error)
    IF ( error /= 0 ) RETURN

#endif
    CALL h5dclose_f( dataset_id, error)
    IF ( error /= 0 ) RETURN

  END SUBROUTINE MPI_Get_3D_Variable_from_HDF5


  SUBROUTINE MPI_Get_2D_Variable_from_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, NLP, NMP, mp_low, mp_high, error )

    INTEGER(HID_T), INTENT(in)    :: file_id
    CHARACTER(*), INTENT(in)      :: variable_name
    INTEGER(HID_T), INTENT(in)    :: plist_id
    INTEGER(HID_T), INTENT(out)   :: filespace
    INTEGER(HID_T), INTENT(in)    :: memspace
    INTEGER, INTENT(in)           :: NLP, NMP
    INTEGER, INTENT(in)           :: mp_low, mp_high
    REAL(prec), INTENT(out)       :: variable(1:NLP,mp_low:mp_high)
    INTEGER(HSIZE_T), INTENT(in)  :: dimensions(1:2)
    INTEGER, INTENT(out)          :: error

    ! Local
    INTEGER(HID_T) :: dataset_id
    INTEGER        :: mp
    INTEGER(HSIZE_T)        :: starts(1:2), counts(1:2), strides(1:2)

    CALL h5dopen_f( file_id, TRIM(variable_name), dataset_id, error )
    IF ( error /= 0 ) RETURN
#ifdef HAVE_MPI

! -- PARALLEL -- !
    CALL h5dget_space_f( dataset_id, filespace, error )
    IF ( error /= 0 ) RETURN

    DO mp = mp_low, mp_high
      starts = (/ 0, mp-1 /)
      counts = (/ 1, 1 /)
      strides = (/ 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )
      IF ( error /= 0 ) RETURN

      CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                      variable(1:NLP,mp), &
                      dimensions, error, memspace, filespace )
      IF ( error /= 0 ) RETURN

    ENDDO

#else

! -- SERIAL -- !

    CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                    variable, dimensions, error)
    IF ( error /= 0 ) RETURN

#endif
    CALL h5dclose_f( dataset_id, error)
    IF ( error /= 0 ) RETURN

  END SUBROUTINE MPI_Get_2D_Variable_from_HDF5


  SUBROUTINE Get_3D_Grid_Variable_from_HDF5( file_id, variable_name, variable, memspace, dimensions, x, y, z, error )

    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:3)
    INTEGER, INTENT(in)          :: x, y, z
    REAL(prec), INTENT(out)      :: variable(1:x,1:y,1:z)
    INTEGER, INTENT(out)         :: error

    ! Local
    INTEGER(HID_T) :: dataset_id

    CALL h5dopen_f( file_id, TRIM(variable_name), &
                     dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Get_3D_Grid_Variable_from_HDF5

  SUBROUTINE Get_2D_Grid_Variable_from_HDF5( file_id, variable_name, variable, memspace, dimensions, x, y, error )

    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:2)
    INTEGER, INTENT(in)          :: x, y
    REAL(prec), INTENT(out)      :: variable(1:x,1:y)
    INTEGER, INTENT(out)         :: error

    ! Local
    INTEGER(HID_T) :: dataset_id

    CALL h5dopen_f( file_id, TRIM(variable_name), &
                     dataset_id, error)
    IF ( error /= 0 ) RETURN
    CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Get_2D_Grid_Variable_from_HDF5

  SUBROUTINE Get_1D_Grid_Variable_from_HDF5( file_id, variable_name, variable, memspace, dimensions, x, error )

    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1)
    INTEGER, INTENT(in)          :: x
    REAL(prec), INTENT(out)      :: variable(1:x)
    INTEGER, INTENT(out)         :: error

    ! Local
    INTEGER(HID_T) :: dataset_id

    CALL h5dopen_f( file_id, TRIM(variable_name), &
                     dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dread_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Get_1D_Grid_Variable_from_HDF5

  SUBROUTINE Get_Integer_3D_Grid_Variable_from_HDF5( file_id, variable_name, variable, memspace, dimensions, x, y, z, error )

    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:3)
    INTEGER, INTENT(in)          :: x, y, z
    INTEGER, INTENT(out)          :: variable(1:x,1:y,1:z)
    INTEGER, INTENT(out)         :: error

    ! Local
    INTEGER(HID_T) :: dataset_id

    CALL h5dopen_f( file_id, TRIM(variable_name), &
                     dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Get_Integer_3D_Grid_Variable_from_HDF5

  SUBROUTINE Get_Integer_1D_Grid_Variable_from_HDF5( file_id, variable_name, variable, memspace, dimensions, x, error )

    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1)
    INTEGER, INTENT(in)          :: x
    INTEGER, INTENT(out)         :: variable(1:x)
    INTEGER, INTENT(out)         :: error

    ! Local
    INTEGER(HID_T) :: dataset_id

    CALL h5dopen_f( file_id, TRIM(variable_name), &
                     dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dread_f( dataset_id, H5T_STD_I32LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Get_Integer_1D_Grid_Variable_from_HDF5

  SUBROUTINE MPI_Add_3D_Grid_Variable_to_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, n_fluxtube, NLP, NMP, mp_low, mp_high, error )

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
    IF ( error /= 0 ) RETURN

    DO mp = mp_low, mp_high
      starts = (/ 0, 0, mp-1 /)
      counts = (/ 1, 1, 1 /)
      strides = (/ 1, 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )
      IF ( error /= 0 ) RETURN

      CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                       variable(1:n_fluxtube,1:NLP,mp), &
                       dimensions, error, memspace, filespace )
      IF ( error /= 0 ) RETURN

    ENDDO
#else

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_IEEE_F64LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

#endif

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE MPI_Add_3D_Grid_Variable_to_HDF5


  SUBROUTINE MPI_Add_2D_Grid_Variable_to_HDF5( file_id, variable_name, variable, filespace, memspace, plist_id, dimensions, NLP, NMP, mp_low, mp_high, error )

    INTEGER(HID_T), INTENT(in)   :: file_id
    CHARACTER(*), INTENT(in)     :: variable_name
    INTEGER(HID_T), INTENT(in)   :: filespace
    INTEGER(HID_T), INTENT(in)   :: memspace
    INTEGER(HID_T), INTENT(in)   :: plist_id
    INTEGER, INTENT(in)          :: NLP, NMP
    INTEGER, INTENT(in)          :: mp_low, mp_high
    REAL(prec), INTENT(in)       :: variable(1:NLP,mp_low:mp_high)
    INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:2)
    INTEGER, INTENT(out)         :: error

    ! Local
    INTEGER(HID_T) :: dataset_id
    INTEGER        :: mp
    INTEGER(HSIZE_T)        :: starts(1:2), counts(1:2), strides(1:2)


#ifdef HAVE_MPI
    CALL h5dcreate_f( file_id, TRIM(variable_name), H5T_IEEE_F64LE, filespace, dataset_id, error, plist_id)
    IF ( error /= 0 ) RETURN

    DO mp = mp_low, mp_high
      starts = (/ 0, mp-1 /)
      counts = (/ 1, 1 /)
      strides = (/ 1, 1 /)

      CALL h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, starts, counts, error, strides, dimensions )
      IF ( error /= 0 ) RETURN

      CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                       variable(1:NLP,mp), &
                       dimensions, error, memspace, filespace )
      IF ( error /= 0 ) RETURN

    ENDDO

#else

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_IEEE_F64LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

#endif

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE MPI_Add_2D_Grid_Variable_to_HDF5


  SUBROUTINE Add_3D_Grid_Variable_to_HDF5( file_id, variable_name, variable, memspace, dimensions, x, y, z, error )

  INTEGER(HID_T), INTENT(in)   :: file_id
  CHARACTER(*), INTENT(in)     :: variable_name
  INTEGER(HID_T), INTENT(in)   :: memspace
  INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:3)
  INTEGER, INTENT(in)          :: x, y, z
  REAL(prec), INTENT(in)       :: variable(1:x,1:y,1:z)
  INTEGER, INTENT(out)         :: error

  ! Local
  INTEGER(HID_T) :: dataset_id

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_IEEE_F64LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_3D_Grid_Variable_to_HDF5

  SUBROUTINE Add_Integer_3D_Grid_Variable_to_HDF5( file_id, variable_name, variable, memspace, dimensions, x, y, z, error )

  INTEGER(HID_T), INTENT(in)   :: file_id
  CHARACTER(*), INTENT(in)     :: variable_name
  INTEGER(HID_T), INTENT(in)   :: memspace
  INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:3)
  INTEGER, INTENT(in)          :: x, y, z
  INTEGER, INTENT(in)          :: variable(1:x,1:y,1:z)
  INTEGER, INTENT(out)         :: error

  ! Local
  INTEGER(HID_T) :: dataset_id

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_STD_I32LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_STD_I32LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_Integer_3D_Grid_Variable_to_HDF5


  SUBROUTINE Add_2D_Grid_Variable_to_HDF5( file_id, variable_name, variable, memspace, dimensions, x, y, error )

  INTEGER(HID_T), INTENT(in)   :: file_id
  CHARACTER(*), INTENT(in)     :: variable_name
  INTEGER(HID_T), INTENT(in)   :: memspace
  INTEGER(HSIZE_T), INTENT(in) :: dimensions(1:2)
  INTEGER, INTENT(in)          :: x, y
  REAL(prec), INTENT(in)       :: variable(1:x,1:y)
  INTEGER, INTENT(out)         :: error

  ! Local
  INTEGER(HID_T) :: dataset_id

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_IEEE_F64LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_2D_Grid_Variable_to_HDF5

  SUBROUTINE Add_1D_Grid_Variable_to_HDF5( file_id, variable_name, variable, memspace, dimensions, x, error )

  INTEGER(HID_T), INTENT(in)   :: file_id
  CHARACTER(*), INTENT(in)     :: variable_name
  INTEGER(HID_T), INTENT(in)   :: memspace
  INTEGER(HSIZE_T), INTENT(in) :: dimensions(1)
  INTEGER, INTENT(in)          :: x
  REAL(prec), INTENT(in)       :: variable(1:x)
  INTEGER, INTENT(out)         :: error

  ! Local
  INTEGER(HID_T) :: dataset_id

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_IEEE_F64LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_IEEE_F64LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_1D_Grid_Variable_to_HDF5

  SUBROUTINE Add_Integer_1D_Grid_Variable_to_HDF5( file_id, variable_name, variable, memspace, dimensions, x, error )
! 1D
  INTEGER(HID_T), INTENT(in)   :: file_id
  CHARACTER(*), INTENT(in)     :: variable_name
  INTEGER(HID_T), INTENT(in)   :: memspace
  INTEGER(HSIZE_T), INTENT(in) :: dimensions(1)
  INTEGER, INTENT(in)          :: x
  INTEGER, INTENT(in)          :: variable(1:x)
  INTEGER, INTENT(out)         :: error

  ! Local
  INTEGER(HID_T) :: dataset_id

    CALL h5dcreate_f( file_id, TRIM(variable_name), &
                      H5T_STD_I32LE, memspace, dataset_id, error)
    IF ( error /= 0 ) RETURN

    CALL h5dwrite_f( dataset_id, H5T_STD_I32LE, &
                     variable, dimensions, error)
    IF ( error /= 0 ) RETURN

    CALL h5dclose_f( dataset_id, error)

  END SUBROUTINE Add_Integer_1D_Grid_Variable_to_HDF5


END MODULE IPE_HDF5


