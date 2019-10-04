PROGRAM Eldyn_Converter

USE IPE_Precision
USE IPE_Model_Class

IMPLICIT NONE


  TYPE( IPE_Model ) :: ipe
  LOGICAL           :: init_success
  REAL(prec)        :: t0_ipe
  INTEGER           :: i, rc
  CHARACTER(50)     :: filename

    filename='potential_SMG_20150316_200100.nc'
    CALL ipe % Build( init_success )

    IF( init_success )THEN

      DO i = 1, ipe % parameters % n_model_updates

        t0_ipe = ipe % parameters % start_time + REAL(i-1,prec)*ipe % parameters % file_output_frequency

        CALL ipe % time_tracker % Update( t0_ipe )

        ! Here is where the electrodynamics "Update" routine is called.
        ! You can feel free to comment it out, and replace it with calls
        ! to your "Read", "Regrid", and "Merge" routines.

        CALL ipe % eldyn % Read_Geospace_Potential( filename, rc )

        print *,'reading Geospace'
        CALL ipe % eldyn % Interpolate_Geospace_to_MHDpotential ( ipe % grid,ipe % time_tracker)

        ! This syntax will also work
        !CALL Read_Geospace_Potential( ipe % eldyn, filename, rc )


        CALL ipe % eldyn % Write_MHD_Potential( ipe % grid, &
                                                ipe % time_tracker, &
                                               "MHD_Potential"//ipe % time_tracker % DateStamp( )//".nc", &
                                                rc )
        print *,'writing geospace'

      ENDDO

    ENDIF



CONTAINS


! In this CONTAINS region ( before the END PROGRAM statement and after
! "CONTAINS" ), you can add whatever subroutines/functions you want and can call
! them within this program. Think of this area as your scratch space for quickly
! prototyping routines that accomplish what  you want. Once you're happy with
! the routine, then we can work on pushing it into the
! IPE_Electrodynamics_Class.F90 module.

END PROGRAM Eldyn_Converter
