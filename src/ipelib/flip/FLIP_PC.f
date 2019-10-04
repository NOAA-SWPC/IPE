C::::::::::::::::::::::: PC_CLOSE :::::::::::::::::
C... This routine for PC version of FLIP. 
C... CLOSE file on VAX to avoid file locking but REWIND on PC. 
      SUBROUTINE PC_CLOSE(UNIT)
      INTEGER UNIT
      REWIND(UNIT)
      RETURN
      END
