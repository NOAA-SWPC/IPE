MODULE IPE_Common_Routines

USE IPE_Precision

#ifdef HAVE_NETCDF
USE netcdf
#endif

IMPLICIT NONE

CONTAINS

  LOGICAL FUNCTION Almost_Equal( a, b )
    IMPLICIT NONE
    REAL(prec) :: a, b


      IF( a == 0.0_prec .OR. b == 0.0_prec )THEN

        IF( ABS(a-b) <= EPSILON(1.0_prec) )THEN

          Almost_Equal = .TRUE.

        ELSE

          Almost_Equal = .FALSE.

        ENDIF

      ELSE

        IF( (ABS(a-b) <= EPSILON(1.0_prec)*ABS(a)) .OR. (ABS(a-b) <= EPSILON(1.0_prec)*ABS(b)) )THEN

          Almost_Equal = .TRUE.

        ELSE

          Almost_Equal = .FALSE.

        ENDIF

      ENDIF

  END FUNCTION Almost_Equal

  INTEGER FUNCTION NewUnit(thisunit)
    IMPLICIT NONE
    INTEGER, INTENT(out), optional :: thisunit
    ! Local
    INTEGER, PARAMETER :: unitMin=100, unitMax=1000
    LOGICAL :: isopened
    INTEGER :: iUnit

      newunit=-1

      DO iUnit=unitMin, unitMax

         INQUIRE(UNIT=iUnit,opened=isopened)

         if( .not. isopened )then
            NewUnit = iUnit
            EXIT
         ENDif

      ENDDO

      IF( PRESENT(thisunit) ) thisunit = NewUnit

  END FUNCTION NewUnit

#ifdef HAVE_NETCDF
  SUBROUTINE Check(status)
    IMPLICIT NONE
    INTEGER, INTENT (in) :: status

     IF(status /= nf90_noerr) THEN
       PRINT *, trim(nf90_strerror(status))
       STOP "NetCDF Error, Stopped"
     ENDIF
  END SUBROUTINE Check
#endif

  SUBROUTINE Hat_Weights( x_grid, x_interp, weights, indices, N )
    INTEGER, INTENT(in)     :: N
    REAL(prec), INTENT(in)  :: x_grid(0:N)
    REAL(prec), INTENT(in)  :: x_interp
    REAL(prec), INTENT(out) :: weights(1:2)
    INTEGER, INTENT(out)    :: indices(1:2)
    ! Local
    INTEGER :: i, j


      weights = 0.0_prec
      indices = 0
      j = 0


      ! ---- Check Endpoints ---- !
      IF( x_interp < x_grid(0) .AND. x_interp < x_grid(1) )THEN

        j = j+1
        weights(j) = 1.0_prec
        indices(j) = 0

        j = j+1
        weights(j) = 0.0_prec
        indices(j) = 1
        RETURN

      ENDIF

      IF( x_interp >= x_grid(N) )THEN
        j = j+1
        weights(j) = 1.0_prec
        indices(j) = N

        j = j+1
        weights(j) = 0.0_prec
        indices(j) = N-1

        RETURN

      ENDIF
      ! ------------------------------------- !

      IF( x_interp >= x_grid(0) .AND. x_interp < x_grid(1) )THEN

        j = j + 1
        weights(j) = ( x_interp - x_grid(1) )/( x_grid(0) - x_grid(1) )
        indices(j) = 0

      ENDIF


      DO i = 1, N-1

        IF( x_interp >= x_grid(i-1) .AND. x_interp < x_grid(i) )THEN

          j = j+1
          weights(j) = ( x_interp - x_grid(i-1) )/( x_grid(i) - x_grid(i-1) )
          indices(j) = i

        ELSEIF( x_interp >= x_grid(i) .AND. x_interp < x_grid(i+1) )THEN

          j = j+1
          weights(j) = ( x_interp - x_grid(i+1) )/( x_grid(i) - x_grid(i+1) )
          indices(j) = i

        ENDIF

        IF( j == 2 )THEN
          RETURN
        ENDIF

      ENDDO

      IF( x_interp >= x_grid(N-1) .AND. x_interp < x_grid(N) )THEN

        j = j+1
        weights(j) = ( x_interp - x_grid(N-1) )/( x_grid(N) - x_grid(N-1) )
        indices(j) = N

      ENDIF

  END SUBROUTINE Hat_Weights

END MODULE IPE_Common_Routines
