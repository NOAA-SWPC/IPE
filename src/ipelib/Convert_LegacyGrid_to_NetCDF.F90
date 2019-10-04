PROGRAM Convert_LegacyGrid_to_NetCDF

USE IPE_Grid_Class
USE IPE_Model_Parameters_Class
USE IPE_MPI_Layer_Class


IMPLICIT NONE

  TYPE( IPE_Grid )             :: grid
  TYPE( IPE_Model_Parameters ) :: params
  TYPE( IPE_MPI_Layer )        :: mpi_layer

  CHARACTER(500) :: gridFile, interpFile, ncFile
  LOGICAL :: initSuccess, readSuccess
  INTEGER :: rc


    CALL InitializeFromCommandLine( gridFile, interpFile, ncFile, initSuccess )

    IF( initSuccess )THEN

      CALL params % Build( mpi_layer, readSuccess )

      IF( readSuccess )THEN

        CALL grid % Build( params % nFluxTube, params % NLP, params % NMP, params % NPTS2D, 1, params % NMP, 1 )

        CALL grid % Read_IPE_Grid( gridFile )

        IF( TRIM(interpFile) /= '' )THEN

          CALL grid % Read_Legacy_Geographic_Weights( interpFile )

        ENDIF

        CALL grid % Write_IPE_Grid_NetCDF( ncFile, rc )

        CALL grid % Trash( )

      ELSE

        PRINT*, ' '
        PRINT*, '  Re-run Legacy2NetCDF with the generated IPE.inp file'
        PRINT*, ' '

      ENDIF

    ENDIF

CONTAINS

 SUBROUTINE InitializeFromCommandLine( gridFile, weightFile, outFile, initSuccess )
   IMPLICIT NONE
   CHARACTER(*), INTENT(out) :: gridFile, weightFile, outFile
   LOGICAL, INTENT(out)      :: initSuccess
   ! Local
   LOGICAL        :: interpSuccess, gridSuccess
   INTEGER        :: nArg, argID
   CHARACTER(500) :: argname
   LOGICAL        :: fileExists
   LOGICAL        :: gridGiven, gridFileGiven, helpRequested
   LOGICAL        :: weightsGiven, weightFileGiven
   LOGICAL        :: outGiven, outFileGiven

     gridFileGiven   = .FALSE.
     gridGiven       = .FALSE.
     weightsGiven    = .FALSE.
     weightFileGiven = .FALSE.
     outGiven        = .FALSE.
     outFileGiven    = .FALSE.
     helpRequested   = .FALSE.
     initSuccess     = .FALSE.

     ! Default grid file
     gridFile      = './ipe_grid'
     outFile       = './IPE_Grid.nc'
     weightFile    = ''

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN

       DO argID = 1, nArg

         CALL get_command_argument( argID, argName )

         SELECT CASE( TRIM(argName) )

           CASE("--legacy-grid")

             gridGiven = .TRUE.

           CASE("--interp-file")

             weightsGiven = .TRUE.

           CASE("--nc-file")

             outGiven = .TRUE.

           CASE("--help")

             helpRequested = .TRUE.

           CASE DEFAULT

             IF( gridGiven )THEN

               gridFileGiven = .TRUE.
               gridFile  = TRIM(argName)
               gridGiven = .FALSE.

             ENDIF

             IF( weightsGiven )THEN

               weightFileGiven = .TRUE.
               weightFile      = TRIM(argName)
               weightsGiven    = .FALSE.

             ENDIF

             IF( outGiven )THEN

               outFileGiven = .TRUE.
               outFile      = TRIM(argName)
               outGiven     = .FALSE.

             ENDIF

         END SELECT
       ENDDO

     ENDIF

     IF( helpRequested )THEN


       PRINT*, '  grid2netcdf '
       PRINT*, ' '
       PRINT*, '  Usage : grid2netcdf [OPTIONS] '
       PRINT*, ' '
       PRINT*, '    Options '
       PRINT*, ' '
       PRINT*, '      --help '
       PRINT*, '          Displays this message '
       PRINT*, ' '
       PRINT*, '      --legacy-grid <grid-file>'
       PRINT*, '          Specifies to use <grid-file> as the legacy grid file '
       PRINT*, '          to be converted to NetCDF. If this option is not  '
       PRINT*, '          provided, ./ipe_grid is assumed.  '
       PRINT*, ' '
       PRINT*, '      --interp-file <interpolation-weights-file>'
       PRINT*, '          Specifies to use <interpolation-weights-file> as the '
       PRINT*, '          legacy file containing weights to map from apex grid '
       PRINT*, '          to the geographic grid. '
       PRINT*, ' '
       PRINT*, '      --nc-file <output-nc-file>'
       PRINT*, '          Specifies the name of the output netcdf file to use.'
       PRINT*, '          If this option is not provided, IPE_Grid.nc is '
       PRINT*, '          assumed. '
       PRINT*, ' '

     ENDIF

     IF( .NOT. helpRequested )THEN

       INQUIRE( FILE = TRIM(gridFile), &
                EXIST = fileExists )

       IF( .NOT.(fileExists) )THEN
         PRINT*, '  Grid file :'//TRIM(gridFile)//': not found.'
         PRINT*, '  Setup Failed.'
         gridSuccess = .FALSE.
       ELSE
         PRINT*, ' '
         PRINT*, ' Using Legacy File :', TRIM(gridFile)
         PRINT*, ' '
         gridSuccess = .TRUE.
       ENDIF

       IF( weightFileGiven )THEN

         INQUIRE( FILE = TRIM(weightFile), &
                  EXIST = fileExists )

         IF( .NOT.(fileExists) )THEN
           PRINT*, '  Interp file :'//TRIM(weightFile)//': not found.'
           PRINT*, '  Setup Failed.'
           interpSuccess = .FALSE.
         ELSE
           PRINT*, ' '
           PRINT*, ' Using Legacy File :', TRIM(gridFile)
           PRINT*, ' '
           interpSuccess = .TRUE.
         ENDIF

       ELSE

         interpSuccess = .TRUE.

       ENDIF

       initSuccess = (gridSuccess .AND. interpSuccess)

     ENDIF

 END SUBROUTINE InitializeFromCommandLine

END PROGRAM Convert_LegacyGrid_to_NetCDF
