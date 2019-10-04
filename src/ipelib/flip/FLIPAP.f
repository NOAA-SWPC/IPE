C............................. FLIPAP.FOR .......................
C.... .. Major rewrite by P. Richards March 2009 to simplify the logic
C....... Routine GET_F107AP replaces FINDAP
C....... This subroutine reads a data file and provides the F10.7 and 7 Ap 
C....... indices for the MSIS model given the time history of the Kp at 
C....... 3-hour intervals.
C::::::::::::::::::::::::: GET_F107AP :::::::::::::::::::::::::::
C....  This program reads the F10KP.DDD data file and returns Ap, F10.7 and F10.7A
C.... Written by P. Richards November 2009
      SUBROUTINE GET_F107AP(FSTDAY   !.. Start day for run (YYYYddd)
     >                      ,UTHRS   !.. Current time on run (hrs)
     >                     ,RUNHRS   !.. Length of run (hrs)     
     >                         ,AP   !.. Output: Ap values
     >                         ,SW   !.. Output: MSIS switch values
     >                       ,F107   !.. Output: F10.7
     >                     ,F107A    !.. Output: F10.7A
     >                     ,KPREAL)   !.. Output: Kp for current time
      IMPLICIT NONE
	INTEGER APDIM
      PARAMETER (APDIM=33001)
      INTEGER I,IND         !.. loop control variables
	INTEGER FSTDAY,STDAY
      INTEGER AP3HR(APDIM),T107(APDIM),T107A(APDIM),KP3HR(APDIM)
	INTEGER NUMAP3,IDAY,INDSAV,IOPEN
	REAL SW(25),UTHRS,RUNHRS,F107,F107A,AP(7),KPREAL
	LOGICAL APREAD
	DATA APREAD/.FALSE./,INDSAV/0/

      !..On first call, get the Ap,F10.7, and F10.7A into temporary arrays
      IF (APREAD.EQV..FALSE.) THEN
	   APREAD=.TRUE.
	   STDAY=FSTDAY  !.. save the first day
         CALL TFILE(20,IOPEN)  !... Test to see if the data file exists
         IF(IOPEN.EQ.0) THEN
           WRITE(6,*) '*** File with 3 hour Ap and F10.7 is missing ***'
           CALL RUN_ERROR    !.. print ERROR warning in output files
           STOP
         ENDIF
         CALL READ_F10KP(APDIM,STDAY,(RUNHRS+UTHRS),NUMAP3,T107,T107A,
     >      AP3HR,KP3HR)
         CLOSE(20)
         !.. call TSELEC and GSELEC to allow MSIS and HWM to use the 3 hour AP
         SW(9)= -1.0
         CALL TSELEC(SW) 
         CALL GSELEC(SW) 
      ENDIF

      !.. IND - index of the 3 hour AP for the current time - UTHRS
	!.. The 25 takes care of the 3 days before the current date
      IND = 25 + INT(UTHRS/3)

      !.... If the index has not changed do not recalculate APs
      IF (IND .EQ. INDSAV) THEN
         RETURN
      ELSE
         INDSAV = IND
      ENDIF

      F107 = REAL(T107(IND))   !.. Transfer F10.7 to MSIS parameters
      F107A= REAL(T107A(IND))

      !... Calculate the daily Ap (AP(1)) by first determining the day which
      !... corresponds to the current UTHRS, then calculate the index that
      !... corresponds to the first 3 hour Ap of that day. Then sum 8 3 hour
      !... APs for that day
      IDAY = 25 + INT(UTHRS/24) * 8
      AP(1)=0.0
      DO I = IDAY, (IDAY+7)
         AP(1) = AP(1)+REAL(AP3HR(I))
      ENDDO
      AP(1)=AP(1)/8.0

      !... Fill the next 4 elements of AP with the corresponding 3 hour values
      DO I = 2, 5
        AP(I) = REAL(AP3HR(IND-I+2))
      ENDDO

      !.... AP(6) is the average of eight 3 hour indices from 12 to 33 hours
      !.... prior to current time
      AP(6) = 0.0
      DO I = 4, 11
        AP(6) = AP(6) + REAL(AP3HR(IND-I))
      ENDDO
      AP(6) = AP(6)/8.0

      !.... AP(7) is the average of eight 3 hour indices from 36 to 60 hours
      !.... prior to current time
      AP(7) = 0.0
      DO I = 12, 19
        AP(7) = AP(7) + REAL(AP3HR(IND-I))
      ENDDO
      AP(7) = AP(7)/8.0

      KPREAL=REAL(KP3HR(IND))/10.0

 90   CONTINUE
      RETURN
      END
C::::::::::::::::::::::::: SUBROUTINE CONVT_DATE :::::::::::::::::::::::::::
C---- Scot A Braze,  8/7/95
C---- This subroutine converts the date from YYYY MM DD to YYYY DDD. 
      SUBROUTINE CONVT_DATE(MM,DD,DDD,YYYY)
      IMPLICIT NONE
      INTEGER DDD,MM,DD,NCONV(12),LYCONV(12),YYYY

      DATA NCONV/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA LYCONV/0,31,60,91,121,152,182,213,244,274,305,335/

      IF(MOD(YYYY,4).EQ.0) THEN
          DDD=LYCONV(MM)+DD
      ELSE
          DDD=NCONV(MM)+DD
      ENDIF
      END
C::::::::::::::::::::::::: SUBROUTINE KP_START_DATE :::::::::::::::::::::::::::
C---- This subroutine takes the start date and then goes back 3 days to get Kp history. 
C---- Written by P. Richards March 2009.
      SUBROUTINE KP_START_DATE(STDAY   !.. Start date for the FLIP run
     >                     ,YYYYDDD)   !.. Start date for the Kp/Ap values
      IMPLICIT NONE
      INTEGER STDAY,YYYYDDD,SYYYY,SDDD,FYYYY,FDDD
      !.. Get the day the FLIP run begins.  FYYYY and FDDD are
      !.. the flip starting year and day for the run.  
      FYYYY=STDAY/1000
      FDDD=MOD(STDAY,1000)
      !.. Check to see if 3 days prior to run is in the previous year.
      IF(FDDD.LT.4) THEN
          IF(MOD(FYYYY-1,4).EQ.0) THEN
              SDDD=366-(3-FDDD)
              SYYYY=FYYYY-1
          ELSE
              SDDD=365-(3-FDDD)
              SYYYY=FYYYY-1
          ENDIF
      ELSE
          SDDD=FDDD-3
          SYYYY=FYYYY
      ENDIF
      YYYYDDD=SYYYY*1000+SDDD !.. return the day 3 days before start date
      RETURN
      END
C::::::::::::::::::::::::::READ_F10KP ::::::::::::::::::::::::::
C..... This subroutine reads in the F10.7 and Kp indices from the master file,
C..... converts Kp to Ap and returns AP and F10.7  in arrays. 
C..... Written by P. Richards March 2009.
      SUBROUTINE READ_F10KP(APDIM      !.. Dimension of arrays
     >                     ,STDAY      !.. Start day for run (YYYYddd)
     >                    ,RUNHRS      !.. Length of run (hrs)     
     >                    ,NUMAP3      !.. Output: # of 3r Ap values
     >                      ,T107      !.. Output: F10.7 array
     >                     ,T107A      !.. Output: F10.7A array
     >                     ,AP3HR      !.. Output: 3 hour Ap array
     >                    ,KP3HR)      !.. Output: 3 hour Kp array
      IMPLICIT NONE
	INTEGER APDIM
      INTEGER I,J,K,IND         !.. loop control variables
	INTEGER YYYY,MM,DD,DDD,YYYYDDD,KP(8),F107,F107A,YSTART,STDAY
      INTEGER AP3HR(APDIM),T107(APDIM),T107A(APDIM),NDAYS,KDAY
      INTEGER NUMAP3,KPTOAP(0:90),KP3HR(APDIM)
	REAL UTAP,RUNHRS
      !.. Conversion factors from Kp = 1..9. AP = KPTOAP(KP)
      !.. Kp= 0    1    2   3    4    5    6     7     8    9
      !.. AP= 0,   3,   7,  15,  27,  48,  80,  140,  240, 400
      !.. 90 Ap values for easy conversion to Ap
      DATA KPTOAP/3*0,4*2,3*3,3*4,4*5,3*6,3*7,4*9,3*12,3*15,4*18,3*22
     >   ,3*27,4*32,3*39,3*49,4*56,3*67,3*80,4*94,3*111,3*132,4*154
     >   ,3*179,3*207,4*236,3*300,400/

      !.. read headers off the input file
 10   READ(20,*,ERR=10,END=99) YYYY,MM,DD,(KP(J),J=1,8),F107,F107A

      !.. Gets the day 3 days before current day for AP history
      CALL KP_START_DATE(STDAY,YSTART)

      !.. Convert year, month, day to YYYYddd format
      CALL CONVT_DATE(MM,DD,DDD,YYYY)
      YYYYDDD=YYYY*1000+DDD 
      IF(YSTART.LT.YYYYDDD) THEN
        WRITE(6,*) 
     >  ' ** Start date is before first date in F10Kp file **'
        CALL RUN_ERROR    !.. print ERROR warning in output files
        STOP
      ENDIF

      NDAYS = NINT((RUNHRS)/24.0)+3  !.. total number of days in the run
      K=0
      UTAP=0.0  !.. only used for debugging
      !.. Read in the F107 and Kp data
      DO I=1,APDIM-1
        !.. Convert year, month, day to YYYYddd format
        CALL CONVT_DATE(MM,DD,DDD,YYYY)
        YYYYDDD=YYYY*1000+DDD 
        IF(KDAY.GT.NDAYS) GOTO 91
        !.. Now, if file date exceeds start date write the Kp at 3 hour intervals 
        IF(YYYYDDD.GE.YSTART) THEN
          KDAY=KDAY+1
          DO J=1,8
            K=K+1
            KP3HR(K)=KP(J)  !.. Note Kp mult. by 10 in data file
            T107(K)=F107
            T107A(K)=F107A
            IF(YYYYDDD.GE.STDAY) UTAP=UTAP+3.0
          ENDDO
        ENDIF
        !.. Get the next day of data from file
	  READ(20,*,ERR=99,END=91) YYYY,MM,DD,(KP(J),J=1,8),F107,F107A
      ENDDO
 91   CONTINUE
      NUMAP3=K-1

      !.. Test F107 for bad data.
      DO I=1,NUMAP3
        IF(T107(I).LT.50.0.OR.T107(I).GT.400.OR.
     >    T107A(I).LT.50.0.OR.T107A(I).GT.300) THEN
          WRITE(6,96) 
 96       FORMAT(/'  ERROR!, Invalid F10.7 or F107A value'
     >    ,1X,'in F10_Kp data file'/)
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
        ENDIF
      ENDDO

      !.. Check the Kp for bad data and convert to Ap. Note Kp mult. by 10 in 
      !.. data file
      DO I=1,NUMAP3
         IF(KP3HR(I).LT.0.OR.KP3HR(I).GT.90) THEN
            WRITE(6,*) ' **** Kp index in file must be 0<= Kp <=9  ****'
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ELSE 
            AP3HR(I) = KPTOAP(KP3HR(I))  !.. Convert Kp to Ap
         ENDIF
      ENDDO
      RETURN

 99   CONTINUE
      WRITE(6,*) ' ** Not enough days in the F10Kp file for FLIP run **'
      CALL RUN_ERROR    !.. print ERROR warning in output files
      STOP
      RETURN
      END
C:::::::::::::::::::::: GSELEC ::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GSELEC(SV)
C        SET SWITCHES
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
      DIMENSION SV(1),SAV(25),SVV(1)
      COMMON/PR_CSW/SW(25),ISW,SWC(25)
      DO 100 I = 1,25
        SAV(I)=SV(I)
        SW(I)=AMOD(SV(I),2.)
        IF(ABS(SV(I)).EQ.1.OR.ABS(SV(I)).EQ.2.) THEN
          SWC(I)=1.
        ELSE
          SWC(I)=0.
        ENDIF
  100 CONTINUE
      ISW=64999
      RETURN
      ENTRY TRETRW(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
      