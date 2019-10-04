PROGRAM Convert_LegacyPlasma_to_NetCDF

USE IPE_Model_Class


IMPLICIT NONE

  TYPE( IPE_Model ) :: ipe
  
  CHARACTER(500) :: plasmaFile
  LOGICAL :: initSuccess, readSuccess


    CALL InitializeFromCommandLine( plasmaFile, initSuccess )

    IF( initSuccess )THEN

      CALL ipe % Build( readSuccess )

      IF( readSuccess )THEN

        CALL ipe % plasma % Read_Legacy_Input( ipe % grid, plasmaFile )
        
        CALL ipe % Write_NetCDF_IPE( "IPE_State.apex."//ipe % time_tracker % DateStamp( )//".nc" ) 

        CALL ipe % Geographic_Interpolation( ) 
        CALL ipe % Write_Geographic_NetCDF_IPE( "IPE_State.geo."//ipe % time_tracker % DateStamp( )//".nc" ) 

        CALL ipe % Trash( )

      ELSE

        PRINT*, ' '
        PRINT*, '  Re-run Legacy2NetCDF with the generated IPE.inp file' 
        PRINT*, ' '

      ENDIF
 
    ENDIF

CONTAINS

 SUBROUTINE InitializeFromCommandLine( plasmaFile, initSuccess )
   IMPLICIT NONE
   CHARACTER(*), INTENT(out) :: plasmaFile
   LOGICAL, INTENT(out)      :: initSuccess
   ! Local
   INTEGER        :: nArg, argID
   CHARACTER(500) :: argname
   LOGICAL        :: fileExists
   LOGICAL        :: helpRequested
   LOGICAL        :: plasmaFileGiven, fileGiven

     helpRequested   = .TRUE.
     initSuccess     = .FALSE.
     fileGiven       = .FALSE.

     ! Default grid file

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
       DO argID = 1, nArg
  
         CALL get_command_argument( argID, argName )

         SELECT CASE( TRIM(argName) )

           CASE("--plasma-file")

             fileGiven = .TRUE.
             helpRequested = .FALSE.

           CASE("--help")

             helpRequested = .TRUE.

           CASE DEFAULT

             IF( fileGiven )THEN

               plasmaFileGiven = .TRUE.
               plasmaFile      = TRIM(argName)
               fileGiven       = .FALSE.

             ENDIF

         END SELECT
       ENDDO

     ENDIF

     IF( helpRequested )THEN


       PRINT*, '  legacy2netcdf ' 
       PRINT*, ' '
       PRINT*, '  Usage : legacy2netcdf [OPTIONS] '
       PRINT*, ' '
       PRINT*, '    Options '
       PRINT*, ' '
       PRINT*, '      --help '
       PRINT*, '          Displays this message '
       PRINT*, ' '
       PRINT*, '      --plasma-file <file>'
       PRINT*, '          Specifies to read in <file> as the legacy IPE plasma file '
       PRINT*, '          to be converted to NetCDF. If this option is not  '
       PRINT*, '          provided, this program will fail  '
       PRINT*, ' '

     ENDIF

     IF( .NOT. helpRequested )THEN 

       INQUIRE( FILE = TRIM(plasmaFile), &
                EXIST = fileExists )

       IF( .NOT.(fileExists) )THEN
         PRINT*, '  Plasma file :'//TRIM(plasmaFile)//': not found.'
         PRINT*, '  Setup Failed.'
         initSuccess = .FALSE.
       ELSE
         PRINT*, ' '
         PRINT*, ' Using Plasma File :', TRIM(plasmaFile)
         PRINT*, ' '
         initSuccess = .TRUE.
       ENDIF


     ENDIF

 END SUBROUTINE InitializeFromCommandLine

END PROGRAM Convert_LegacyPlasma_to_NetCDF
