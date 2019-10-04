C:::::::::::::::::::::GET_CONTR_PARAMS::::::::::::::::::::::::::::::::::::
C.... This routine is used to extract control parameter values from
C.... control file(eg.FYYYYDDD.RUN or FLIPRINT.RUN) and write them into files
C.... instead of input from keyboard. 
C.... Use &=START=& and &=END=& to identify the range of data.
C.... Use ! to identify the data in a line
C.... This routine is called in OPEN_FILE() routine.

      SUBROUTINE GET_CONTR_PARAMS(INUNIT,OUTUNIT)
	INTEGER INUNIT                 ! Input unit 
	INTEGER OUTUNIT                ! Output unit 
         
      CHARACTER*200 LINE             ! Line buffer
      
	!... read each line from input file until find out flag
  50	READ(INUNIT,'(A200)',END=100) LINE
      DO I=1,200
	  IF(INDEX(LINE,'&=START=&').NE.0) GO TO 60 
	ENDDO
	GO TO 50
	!... parse the line and pick out the data
  60  READ(INUNIT,'(A200)',END=100) LINE
      
	!... if meet end flag,return
      IF(INDEX(LINE,'&=END=&').NE.0) GO TO 100
	!... use "!" as a flag to identify data
      I1=INDEX(LINE,'!')
	!... output data into output file
	IF(I1.NE.0) THEN
	   WRITE(OUTUNIT,*) LINE(1:I1-1)
	ENDIF    
	GO TO 60  
  100 RETURN 
      END