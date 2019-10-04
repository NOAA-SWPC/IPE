C........................ FLIPWIND.FOR .................
C.. This file contains subroutines associated with neutral winds
C.. Created in July 1998
C.. It also has the routines for the modified MSIS Tn
       !.. evaluate the magnetic field line cmpt. of neutral wind UN(J)
       !.. from Hedin's model (BCOMPU=component along magnetic meridian)
       !.. evaluate the magnetic field line cmpt. of neutral wind UN(J)
       !.. from Hedin's model.
       SUBROUTINE HEDWIN(JDAY,SX,ALT,GLATJ,GLON,SAT,F107A,F107
     >         ,AP,BCOMPU)
       IMPLICIT NONE
       INTEGER JDAY
       REAL SX,ALT,GLATJ,GLON,SAT,F107A,F107,AP(7),UHED(2),BCOMPU
     >   ,DECMAG

c       CALL GWS4(JDAY,SX,ALT,GLATJ,GLON,SAT,F107A,F107
c     >         ,AP,UHED)
       CALL hwm14(JDAY,SX,ALT,GLATJ,GLON,SAT,F107A,F107
     >         ,AP,UHED)
       CALL MAGDEC(GLATJ,GLON,DECMAG) !-- magnetic declination for winds
       DECMAG=DECMAG/57.296
       BCOMPU = (UHED(1)*COS(DECMAG)+UHED(2)*SIN(DECMAG))
       RETURN
       END
C::::::::::::::::::::::::::: ALTWIN ::::::::::::::::::::::::
C...... This subroutine provides alternative winds to Hedin's winds
C...... Winds can be read directly from a file or calculated from hmF2
C...... Either from the file data or the IRI model.
C...... Written by P. Richards in April 1991
      SUBROUTINE ALTWIN(DT,Z,N,GL,ISOL,IWIND,TIMSAV,TSTOP,IWTIM,IWW
     > ,JMAX,COSDIP,UNN,UNC,HMF2N,HMF2S,NMF2N,NMF2S,UN,DNMF2N
     > ,DNMF2S,DHMF2N,DHMF2S,RATIN,RATIS,TOPTEN,TOPTES,ZTOPTE,APFILE)
C
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      INTEGER NWTIMN,NWTIMS,NCOLN,NCOLS,UTOLTN,UTOLTS,HMOWNN,HMOWNS
     > ,TPRIN,JDAY,RCLPD,TMODON,APFILE
      REAL SEC,F107A,SATN,SATF,GLONN,GLONF,GLATN,GLATF,RHMF2,RNMF2
     > ,DEC,BLON,AP,D(19),T(2),ETRAN,F107,CHIN,CHIF,FOE,ZTOPTE,SX
     > ,BCOMPUN,BCOMPUS,NDFNOR,NDFSOU,TINFN,TINFS,DELTIN
     > ,DELTIS,TMSN,TMSS,EXFACN,EXFACS
      DIMENSION COSDIP(401),UN(401),DATINN(6),DATINS(6),Z(401),N(4,401)
     >  ,GL(401)
      COMMON/WPAR/NWTIMN,NWTIMS,NCOLN(9),NCOLS(9),UTOLTN,UTOLTS
     >  ,HMOWNN,HMOWNS
      COMMON/AWIND/UNNSAV,UNCSAV
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/NDFAC/TPRIN,NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS
     >   ,TMSN,TMSS,EXFACN,EXFACS,TMODON
            
      DATA IDELAY/1/, JWN,JWS/1,1/, HMOWNN,HMOWNS/2*0.0/,TMODON/1/
     >  ,NDFNOR/1/,NDFSOU/1/, TINFN,TINFS/0.0,0.0/, NCOLN/9*0/
     >  ,NCOLS/9*0/, IVTN/0/,IVTS/0/, EXFACN,EXFACS/1.0,1.0/

C
C----- Find altitude points for setting lat, long etc
      IF(JWN.EQ.1) THEN
         DO 20 IL=1,JMAX/2
            IF(Z(IL).LT.300) JWN=IL
 20      CONTINUE
         JWS=JMAX-JWN+1
      ENDIF

      CALL ACTUAL_DAY(IDAY,SEC,JDAY,SX) 
C
C........ Call AMBS to get time, latitude and longitude near the F2 peak
      CALL AMBS(1,Z(JWN),GL(JWN),D,T,CHIN,SATN,GLONN,GLATN)
      CALL AMBS(JMAX,Z(JWS),GL(JWS),D,T,CHIF,SATF,GLONF,GLATF)
C
C........ variables used in calculating winds from HmF2
      IF(ISOL.EQ.0) IDELAY=0
      IF(SEC-TIMSAV.GE.200) THEN
         TIMSAV=SEC
         IDELAY=1
         UNNSAV=UNN/100.0
         UNCSAV=UNC/100.0
      ENDIF

      DNMF2N=0.0
      TOPTEN=0.0
      DNMF2S=0.0
      TOPTES=0.0
      !... TINFS=0.0;  TINFN=0.0  !.. Cannot reset
      !.. neutral wind by reading data from a file and interpolating
      !.. NCOLN(1)= # of columns read from hmF2 file. If NCOL(i)=1
      !.. that column of data will be used by FLIP
      IF(IWIND.GT.0) THEN
        CALL NTHWIN(NWTIMN,SEC/3600.,SATN,TSTOP,COSDIP(JWN)
     >    ,HMF2N,UNNSAV,DATINN,NCOLN,UTOLTN,HMOWNN,ZTOPTE)
        UNN=DATINN(2)
        IF(HMOWNN.EQ.1.AND.NCOLN(2).EQ.1) DHMF2N=DATINN(1)
        IF(NCOLN(3).NE.0) DNMF2N=DATINN(3)
        IF(NCOLN(4).NE.0) TOPTEN=DATINN(4)
        IF(NCOLN(6).NE.0) NDFNOR=DATINN(6)
        IF(NCOLN(5).NE.0) TINFN=DATINN(5)
        CALL HEDWIN(JDAY,SX,250.0,GLATN,GLONN,SATN,F107A,F107
     >         ,AP,BCOMPUN)
        IF(NCOLN(2).EQ.-2) UNN=100.0*BCOMPUN+DATINN(2)

        !.... Southern hemisphere
        CALL STHWIN(NWTIMS,SEC/3600.,SATF,TSTOP,COSDIP(JWN)
     >    ,HMF2S,UNCSAV,DATINS,NCOLS,UTOLTS,HMOWNS,ZTOPTE)
        UNC=DATINS(2)
        IF(HMOWNS.EQ.1.AND.NCOLS(2).EQ.1) DHMF2S=DATINS(1)
        IF(NCOLS(3).NE.0) DNMF2S=DATINS(3)
        IF(NCOLS(4).NE.0) TOPTES=DATINS(4)
        IF(NCOLS(6).NE.0) NDFSOU=DATINS(6)
        IF(NCOLS(5).NE.0) TINFS=DATINS(5)
        CALL HEDWIN(JDAY,SX,250.0,GLATF,GLONF,SATF,F107A,F107
     >         ,AP,BCOMPUS)
        IF(NCOLS(2).EQ.-2) UNC=100.0*BCOMPUS+DATINS(2)

        !..... UN=0 for ISOL=0 for getting winds from hmF2
        IF(HMOWNN.EQ.1.AND.IDELAY.EQ.0) UNN=0.0
        IF(HMOWNS.EQ.1.AND.IDELAY.EQ.0) UNC=0.0

      ENDIF


C........ wind from IRI hmF2. UN=0 at start for ISOL=0 for getting 
C........ winds from hmF2. Call IRI model to get hmF2 and nmF2 in North
         IF(NCOLN(1).LT.2.OR.IWIND.LT.0.OR.DHMF2N.GT.550) THEN
            CALL IRIHMF(GLATN,GLONN,F107A,SNGL(HMF2N),JDAY,SATN
     >       ,RHMF2,RCLPD,RNMF2,FOE)
            DHMF2N = RHMF2
            DNMF2N = RNMF2/1.0E6
            CALL HMF2W(SATN,COSDIP(JWN),DHMF2N,HMF2N,UNNSAV,UNN)
            UNN=UNN*100
            IF(IDELAY.EQ.0) UNN=0.0
         ENDIF
C........ IRI in Southern hemisphere
         IF(NCOLS(1).LT.2.OR.IWIND.LT.0.OR.DHMF2S.GT.550) THEN
            CALL IRIHMF(GLATF,GLONF,F107A,SNGL(HMF2S),JDAY,SATF
     >       ,RHMF2,RCLPD,RNMF2,FOE)
            DHMF2S = RHMF2
            DNMF2S = RNMF2/1.0E6
            CALL HMF2W(SATF,COSDIP(JWN),DHMF2S,HMF2S,UNCSAV,UNC)
            UNC=UNC*100
            IF(IDELAY.EQ.0) UNC=0.0
         ENDIF

      !... Adjusting hmF2 to match NmF2 in North
      IF(NCOLN(9).EQ.-2) CALL HM_ADJUST(1,SATN,COSDIP(JWN)
     > ,HMF2N,DHMF2N,NMF2N,DNMF2N,UNNSAV,UNN)
      !... Adjusting hmF2 to match NmF2 in South
      IF(NCOLS(9).EQ.-2) CALL HM_ADJUST(2,SATF,COSDIP(JWN)
     > ,HMF2S,DHMF2S,NMF2S,DNMF2S,UNCSAV,UNC)

      !-- Make sure wind is not too large for printing
      CALL VALIDU(HMF2N,DHMF2N,UNNSAV,UNN)
      CALL VALIDU(HMF2S,DHMF2S,UNCSAV,UNC)

      !.. Set up wind along  grid; sign change for poleward wind. If hmF2 >
      !.. equatorial altitude set wind to 0.0
      DO 15 J=1,JMAX/2+1
         UN(J)= -COSDIP(J)*UNN
         IF(HMOWNN.EQ.1.AND.DHMF2N.GE.Z(JMAX/2)) UN(J)=0.0
 15   CONTINUE
      DO 151 J=JMAX/2,JMAX
         UN(J)= -COSDIP(J)*UNC
         IF(HMOWNS.EQ.1.AND.DHMF2S.GE.Z(JMAX/2)) UN(J)=0.0
 151  CONTINUE
C
C------ calculate nmF2 normalizing factors
      RATIN=-1.0
      RATIS=-1.0
      !.. avoids problem near equator
      IF(DHMF2N.GE.Z(JMAX/2).OR.DHMF2S.GE.Z(JMAX/2)) RETURN
      IF(ISOL.EQ.0) RETURN

      !.. using MSIS modification alg. Richards et al. 1998
      IF(NCOLN(9).EQ.-1.OR.NCOLS(9).EQ.-1.OR.IWIND.EQ.-3) GO TO 178 
      !... Using Richards et al. 1995 density mod.
      IF(NCOLN(3).EQ.1.OR.IWIND.EQ.-2)   RATIN=DNMF2N/NMF2N
      IF(NCOLS(3).EQ.1.OR.IWIND.EQ.-2)   RATIS=DNMF2S/NMF2S
      RETURN
C
C------ Calculate factors for normalizing the MSIS exospheric Tn and [O]
C----- TINFN, TINFS are Tn in north and south, NDFNOR, NDFSOU are also 
C----- used in AMBS to normalize the [O] 
 178  CONTINUE

      !.. Turn on the MSIS modification but not the first time
      IF(TMODON.GE.0) THEN
        TMODON=-1
        RETURN
	ENDIF

      IF(NCOLN(5).NE.0.OR.NCOLS(5).NE.0) THEN
        WRITE(6,767) 
 767    FORMAT(/' *** CANNOT use both Ti/Tn from hmF2 file and'
     >  ,/1X,'Richards et al.[1998] MSIS modification alg. at the'
     >  ,1X,'same time ***'/)
        CALL RUN_ERROR    !.. print ERROR warning in output files
        STOP
      ENDIF

      DELTIN=0.0
      DELTIS=0.0
      IF(SATN.GT.-1.OR.ABS(CHIN-1.5708).GT.0.15) THEN
      IF(NCOLN(3).EQ.1.OR.IWIND.EQ.-3) THEN
         WRITE(6,98) SEC/3600,NMF2N/DNMF2N,NDFNOR,TINFN-TMSN
     >     ,TINFN,TMSN
         IF(IVTN.EQ.0) WRITE(24,86)
         IVTN=1
         WRITE(24,85) (SEC-DT)/3600,DHMF2N,DNMF2N,TMSN,TMSN-TINFN
     >    ,NDFNOR,TINFN,UNN/100,BCOMPUN,UNN/100-BCOMPUN,EXFACN
         IF(TINFN.GT.200) 
     >     CALL GETINF(DT,NMF2N,DNMF2N,NDFNOR,TINFN,DELTIN,TMSN)
      ENDIF
      ENDIF

      IF(SATF.GT.-1.OR.ABS(CHIF-1.5708).GT.0.15) THEN
      IF(NCOLS(3).EQ.1.OR.IWIND.EQ.-3) THEN
         WRITE(6,99) SEC/3600,NMF2S/DNMF2S,NDFSOU,TINFS-TMSS
     >     ,TINFS,TMSS
         IF(IVTS.EQ.0) WRITE(25,87)
         IVTS=1
         WRITE(25,85) (SEC-DT)/3600,DHMF2S,DNMF2S,TMSS,TMSS-TINFS
     >    ,NDFSOU,TINFS,UNC/100,BCOMPUS,UNC/100-BCOMPUS,EXFACS
         IF(TINFS.GT.200) 
     >     CALL GETINF(DT,NMF2S,DNMF2S,NDFSOU,TINFS,DELTIS,TMSS)
      ENDIF
      ENDIF

 85   FORMAT(F10.2,F8.2,1P,E10.2,0P,11F9.2)
 86   FORMAT(1X,'!.. Northern Hemisphere variable Tn run.'
     > ,/1X,'!.. hmF2 & NmF2 are from the current input file'
     > ,/7X,'UT    hmF2     NmF2    MSIS_Tinf DelTn '
     > ,2X,'O_ratio  FLIP_Tinf   Wind    HWM14   DELW  TexFac')
 87   FORMAT(1X,'!.. Southern Hemisphere variable Tn run.'
     > ,/1X,'!.. hmF2 & NmF2 are from the current input file'
     > ,/7X,'UT    hmF2     NmF2    MSIS_Tinf DelTn '
     > ,2X,'O_ratio  FLIP_Tinf   Wind    HWM14   DELW  TexFac')
 95   FORMAT(/'  *** ERROR: Cannot use IRI hmF2 to get winds when AP>0'
     > ,/3X,'in the FLIP run file because IRI uses actual geophysical'
     > ,/3X,'indices for that day.'
     > ,/3X,'Put IWIND=0 in FLIPRUN.DDD to use Hedin''s HWM model or'
     > ,/3X,'read both North and South winds/hmF2 from a file instead'
     > ,1X,'(IWIND=1)')
 98   FORMAT(' North:- UT=',F7.2,' R_NmF2=',F5.2,' OFAC=',F5.2,
     > ' DTn=',F5.0,' TFLIP=',F6.0,' TMSIS=',F6.0)
 99   FORMAT(' South:- UT=',F7.2,' R_NmF2=',F5.2,' OFAC=',F5.2,
     > ' DTn=',F5.0,' TFLIP=',F6.0,' TMSIS=',F6.0)
      RETURN
      END
C::::::::::::::::::::::::::::::: NTHWIN :::::::::::::::::::::::::::::
C....... This subroutine provides the northern hemisphere winds returns
C....... data from a file for the FLIP model. RESULT(2) is wind, NmF2
C....... is RESULT(3), Te is RESULT(4, Ti/Tn is RESULT(5). RESULT(1)
C....... contains interpolated hmF2 if the data file contains hmF2 
C....... rather than WIND.
C....... UTHRS is the universal time the wind will be determined.
C....... RUNHRS is the length of the run (hours). 
C....... SAT is the local times in the north hemisphere.
C....... NWTIM is the number of wind times on the data file. It should
C....... be < 0 on the first call to force the data file to be read.
C....... COSDIP=cos(dip angle of matgnetic field) should be provided by 
C....... main program. HPREV=hmF2 at previous time step and UPREV=wind
C....... at previous time step in north hemisphere. Both them should 
C....... also provided by main program.
C....... The data file is read in a two dimensional array MEDATA(N,M).
C....... The data file contains winds or hmF2 (read in MEDATA(N,2) )
C....... in the northern hemisphere at specied UT or LT (read in 
C....... MEDATA(N,1) ) determined by two sentinels HMORWN and UTORLT 
C....... NCOLS(1)= the number of columns in the data file. The data
C....... file may also contains densities (read in MEDATA(N,3)) and
C....... temperatures (read in MEDATA(N,4)).In July 1998 an additional
C....... row of integers were added to indicate which columns will 
C....... actually be used by FLIP. These are stored in NCOLS(2:9).
C....... NCOLS(1:9) must be 0 on first call. 
C....... If the data file contains winds rather than hmF2, they are
C....... are calculated at UTHRS or SAT, by linear interpolation.
C....... If the data file contain hmF2, it calls HWINPT to calculate the
C....... winds based on FLIP model. JUNIT is the unit number of the data
C....... file which will be fixed as 21 for data file in northern 
C....... hemisphere. 
C....... Programmer - D. A. Burns, June  1988. Modifications by 
C....... Rama Lavu and P. Richards  26 JAN 89, Kevin Zou July 4, 93
C....... Additional modifications in July 1998 to specify colums (by PR)

      SUBROUTINE NTHWIN(NWTIM,UTHRS,SAT,RUNHRS,COSDIP,HPREV,UPREV,
     >               RESULT,NCOLS,UTORLT,HMORWN,ZTOPTE)
C
      PARAMETER (N=99999, M=6, JUNIT=21)
      DOUBLE PRECISION RESULT(M),COSDIP,HPREV,UPREV,RUNHRS
      INTEGER NWTIM, NCOLS(9),UTORLT,HMORWN
      REAL MEDATA(N,M), SAT,UTHRS,ZTOPTE
      SAVE MEDATA
C
C.....If this is the first call to READWD , call RDWIND to read
C..... the MEDATA and NWTIM from the data file
      IF (NWTIM .LT. 0)
     >  CALL RDWIND(N,M,JUNIT,UTHRS,RUNHRS,MEDATA,NWTIM,NCOLS,
     >            UTORLT,HMORWN,ZTOPTE)
C
C..... If column number is 0, no data file exist return to main program.
      IF(NCOLS(1).LT.2) NWTIM=0
      IF(NCOLS(1) .EQ. 0) RETURN
C
C.....  If HMORWN is 1, data file contains hmF2 instead of winds,
C.....  call HWINPT to calculate wind.
      IF (HMORWN .EQ. 1) THEN
          CALL HWINPT(N,M,NWTIM,NCOLS,UTHRS,SAT,RUNHRS,COSDIP,HPREV,
     >                 UPREV,MEDATA,UTORLT,RESULT)
      RETURN
C
C.... If HMORWN is 0, data file contains wind, just call INTWIN do
C.... linear interpolation about the data wind.
      ELSE IF (HMORWN .EQ. 0) THEN      
C....... If UTORLT =1, UT is used in the hemisphere else if UTORLT=0
C....... calculate the winds at the specified local time SAT.
         IF (UTORLT .EQ. 1) THEN
              CALL INTWIN(N,M,NWTIM,NCOLS,UTHRS,MEDATA,RESULT)
         ELSE IF (UTORLT .EQ. 0) THEN
              CALL INTWIN(N,M,NWTIM,NCOLS,SAT,MEDATA,RESULT)
         ELSE
	      WRITE(6,*) 'Problem with UT/LT switch value in data file'
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ENDIF
C
      ELSE
         WRITE(6,*) 'Problem with HMF2/WIND switch value in data file'
         CALL RUN_ERROR    !.. print ERROR warning in output files
         STOP
      ENDIF
C.....convert WIND from m/sec to cm/sec.
      RESULT(2)=100.0*RESULT(2)
C
      RETURN
      END
C
C
C::::::::::::::::::::::::::::::: STHWIN :::::::::::::::::::::::::::::
C....... This subroutine provides the southern hemisphere winds returns
C....... data from a file for the FLIP model. RESULT(2) is wind, NmF2
C....... is RESULT(3), Te is RESULT(4, Ti/Tn is RESULT(5). RESULT(1)
C....... contains interpolated hmF2 if the data file contains hmF2 
C....... rather than WIND.
C....... UTHRS is the universal time the wind will be determined.
C....... RUNHRS is the length of the run (hours). 
C....... SAT is the local times in the south hemisphere.
C....... NWTIM is the number of wind times on the data file. It should
C....... be < 0 on the first call to force the data file to be read.
C....... COSDIP=cos(dip angle of matgnetic field) should be provided by 
C....... main program. HPREV=hmF2 at previous time step and UPREV=wind
C....... at previous time step in south hemisphere. Both them should 
C....... also provided by main program.
C....... The data file is read in a two dimensional array MEDATA(N,M).
C....... The data file contains winds or hmF2 (read in MEDATA(N,2) )
C....... in the southern hemisphere at specied UT or LT (read in 
C....... MEDATA(N,1) ) determined by two sentinels HMORWN and UTORLT 
C....... NCOLS(1)= the number of columns in the data file. The data
C....... file may also contains densities (read in MEDATA(N,3)) and
C....... temperatures (read in MEDATA(N,4)).In July 1998 an additional
C....... row of integers were added to indicate which columns will 
C....... actually be used by FLIP. These are stored in NCOLS(2:9).
C....... NCOLS(1:9) must be 0 on first call. 
C....... If the data file contains winds rather than hmF2, they are
C....... are calculated at UTHRS or SAT, by linear interpolation.
C....... If the data file contain hmF2, it calls HWINPT to calculate the
C....... winds based on FLIP model. JUNIT is the unit number of the data
C....... file which will be fixed as 22 for data file in southern 
C....... hemisphere. 
C....... Programmer - D. A. Burns, June  1988. Modifications by 
C....... Rama Lavu and P. Richards  26 JAN 89, Kevin Zou July 4, 93
C....... Additional modifications in July 1998 to specify colums (by PR)
C
      SUBROUTINE STHWIN(NWTIM,UTHRS,SAT,RUNHRS,COSDIP,HPREV,UPREV,
     >                 RESULT,NCOLS,UTORLT,HMORWN,ZTOPTE) 
C
      PARAMETER (N=99999, M=6, JUNIT=22)
      DOUBLE PRECISION RESULT(M),COSDIP,HPREV,UPREV,RUNHRS
      INTEGER NWTIM, NCOLS(9),UTORLT,HMORWN
      REAL MEDATA(N,M), SAT, UTHRS,ZTOPTE
      SAVE MEDATA
C
C.....If this is the first call to READWD , call RDWIND to read
C..... the MEDATA and NWTIM from the data file
      IF (NWTIM .LT. 0)
     > CALL RDWIND(N,M,JUNIT,UTHRS, RUNHRS, MEDATA,NWTIM,NCOLS,
     >             UTORLT,HMORWN,ZTOPTE)
C
C..... If column number in data file is 0 return to main program.
      IF(NCOLS(1).LT.2) NWTIM=0
      IF(NCOLS(1).EQ.0)  RETURN  
C
C.....  If HMORWN is 1, file contains hmF2 instead of winds, call 
C.....  HWINPT to calculate wind.
      IF (HMORWN .EQ. 1) THEN
        CALL HWINPT(N,M,NWTIM,NCOLS,UTHRS,SAT,RUNHRS,COSDIP,HPREV,UPREV,
     >               MEDATA,UTORLT,RESULT)
      RETURN
C
C.... If HMORWN is 0, data file contains wind, just call INTWIN do
C.... linear interpolation about the data wind.
      ELSE IF (HMORWN .EQ. 0) THEN      
C....... If UTORLT =1, UT is used in the hemisphere else if UTORLT=0
C....... calculate the winds at the specified local time SAT.
         IF (UTORLT .EQ. 1) THEN
              CALL INTWIN(N,M,NWTIM,NCOLS,UTHRS,MEDATA,RESULT)
         ELSE IF (UTORLT .EQ. 0) THEN
              CALL INTWIN(N,M,NWTIM,NCOLS,SAT,MEDATA,RESULT)
         ELSE
            WRITE(6,*) 'Problem with UT/LT switch value in data file'
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ENDIF
C
      ELSE
         WRITE(6,*) 'Problem with HMF2/WIND switch value in data file'
         CALL RUN_ERROR    !.. print ERROR warning in output files
         STOP
      ENDIF
C
C.....convert WIND from m/sec to cm/sec.
      RESULT(2)=100.0*RESULT(2)
C
      RETURN
      END
C
C
C::::::::::::::::::::::::::: INTWIN :::::::::::::::::::::::::::::::::::
C.... This subroutine calculates the values of RESULT (including wind,
C.... density, temparture etc..) at required time TIME by using the 
C.... subroutine FTW. It calls BISPLT to locate the data used for 
C.... interpolation. MEDATA(PDT,NCOLS) and MEDATA(PDT+1, NCOLS) will
C.... be used for interpolation. 
C.... Dr. Phlip Richards, Kevin Zou  at July 8,1993
C
      SUBROUTINE INTWIN(N,M,NWTIM,NCOLS,TIME,MEDATA,RESULT)
C
      INTEGER PDT,NCOLS(9),NWTIM,BEGIN
      DOUBLE PRECISION RESULT(M)
      REAL MEDATA(N,M),TIME
C
C.....Use subroutine BISPLW to locate the data used for interpolation.
      BEGIN=1
      CALL BISPLW(N,M,BEGIN,NWTIM,TIME,MEDATA,PDT)
C----  WRITE(6,*) 'TIME=',TIME,'   ','PDT=',PDT
C....... Use subroutine FTW to interpolate at TIME for RESULT
      CALL FTW(N,M,NWTIM,NCOLS(1),TIME,PDT,MEDATA,RESULT)
      RETURN
      END
C
C 
C:::::::::::::::::::::::::::::: FTW ::::::::::::::::::::::::::::::::::
C...... Simple linear interpolating program. If you need a more precisive 
C...... interpolating for the data, you only need modification this 
C...... suroutine which in fact is the statement FACTOR and RESULT in FTW.
C...... FACTOR is the adjustment coefficient.  
C...... Dr. Phlip Richards, Kevin Zou    July 18, 1993
      SUBROUTINE FTW(N,M,NWTIM,NCOLS,TIME,PDT,MEDATA,RESULT)
C
      INTEGER PDT,NCOLS,NWTIM
      DOUBLE PRECISION FACTOR, RESULT(M)      
      REAL MEDATA(N,M),TIME
C
C
      FACTOR=(TIME-MEDATA(PDT,1))/(MEDATA(PDT+1,1)-MEDATA(PDT,1))
      DO 1 I=2,NCOLS 
      RESULT(I)=FACTOR*(MEDATA(PDT+1,I)-MEDATA(PDT,I))+MEDATA(PDT,I)
 1    CONTINUE
      RETURN
      END
C
C 
C::::::::::::::::::::::: RDWIND :::::::::::::::::::::::::::::
C.... This subroutine reads various data (including universial time(UT)/
C.... local time(LT) in MEDATA(N,1), the wind/hmf2 in MEDATA(N,2), and
C.... if they exist, the densities in MEDATA(N,3), the temparture in 
C.... MEDATA(N,4) etc.. 
C.... Three sentinel number give to: 
C....  1. UTORLT which is used to decide UT/LT;
C....  2. HMORWN which is used to decide hmf2/wind. 
C....  2. NCOLS which tell us the coloum number in the data file. 
C.... A count of the winds read in is made NWTIM. JUNIT is the data file's   
C.... unit number which should be passed from the calling subroutine.
C
      SUBROUTINE RDWIND(N,M,JUNIT,UTHRS,RUNHRS,MEDATA,
     >   NWTIM,NCOLS,UTORLT,HMORWN,ZTOPTE)
C
      REAL MEDATA(N,M),HEADR,HMFWND,UTHRS,ZTOPTE
      INTEGER NCOLS(9),UTORLT,HMORWN,NWTIM,JUNIT,NUMCOL,ICOL
      DOUBLE PRECISION RUNHRS
      CHARACTER*99 CHEADER
      LOGICAL DEBUG
      DATA  DEBUG / .FALSE. /

      DO 3 I=1,9
 3    NCOLS(I)=0

      !...Test to see if the data file exists
      CALL TFILE(JUNIT,IOPEN)
      IF(IOPEN.EQ.0) THEN
        RETURN
      ENDIF

      !.. Loop through input file headers until sentinel value.
      DO 10 I = 1,99
        READ(JUNIT, 87, ERR=10) HEADR 
        UTORLT=NINT(HEADR) 
 87     FORMAT(F8.0)
        GO TO 20
 10   CONTINUE
 20   CONTINUE
      READ(JUNIT,*) HMFWND
      HMORWN=NINT(HMFWND)
      READ(JUNIT,*) NUMSAV
      READ(JUNIT,*) (NCOLS(ICOL),ICOL=1,6)
      DO ICOL=1,6
         IF(NCOLS(ICOL).NE.0) NUMCOL=ICOL
      ENDDO
      READ(JUNIT,*) ZTOPTE

      !.. check for using Te at hmF2 in column 6
      IF(NCOLS(4).EQ.-1) ZTOPTE=1        !.. switches to hmF2
      NCOLS(4)=ABS(NCOLS(4))
      IF(NCOLS(4).NE.1) ZTOPTE=-99.0

      READ(JUNIT,'(80A)') CHEADER   !... junk column header

      !.. use MSIS mod alg.
      IF(NCOLS(3).EQ.-1.AND.NCOLS(2).EQ.1) NCOLS(9)=-1 
      IF(NCOLS(3).EQ.-2.AND.NCOLS(2).EQ.1) NCOLS(9)=-2
      NCOLS(3)=ABS(NCOLS(3))
      IF(NCOLS(3).NE.0) NCOLS(5)=0       !.. can't do both
      IF(NCOLS(3).NE.0) NCOLS(6)=0       !.. can't do both
      NCOLS(1)=NUMCOL

C........ Read in the WTIME, and corresponding MEDATA
      DO  25 I = 1, N 
        READ(JUNIT, *,END=27,ERR=26)(MEDATA(I,J),J=1,NCOLS(1)) 	  
        IF (DEBUG) WRITE (6,*) (MEDATA(I,J), J=1,NCOLS(1)) 	

        !.. Test the HMF2 to be sure it is within the correct range.
        IF(HMORWN .EQ. 1) THEN
          IF (MEDATA(I,2) .GT. 550.0) THEN
            WRITE(6,100) MEDATA(I,2),I
 100        FORMAT(/1X,'WARNING! in hmF2/wind file, HMF2 = ',F12.2,
     >      ' at column 2, row',I7,', IRI hmF2 will be substituted' )
            MEDATA(I,2)=9999
          ELSE IF (MEDATA(I,2) .LT. 180.0) THEN
            WRITE(6,101) MEDATA(I,2),I
 101        FORMAT(/1X,'WARNING! in hmF2/wind file, HMF2 = ',F12.2,
     >      ' at column 2, row',I7,', IRI hmF2 will be substituted' )
            MEDATA(I,2)=9999
          ENDIF

        !.. Test the WIND to be sure it is within the correct range.
        ELSE IF(HMORWN .EQ. 0) THEN
          IF (ABS(MEDATA(I,2)) .GT. 2000.0) THEN
            WRITE(6,110) MEDATA(I,2),I
 110        FORMAT(/1X,'** ERROR!  WIND SPEED is ',F10.3,2X,'TOO BIG!!' 
     >/1X,'** Check the wind data file at column 2  and row ',I3)
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
          ENDIF
        ELSE 
          WRITE(6,*)' ERROR! Problem with HMF2/WIND switch in data file'
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
        ENDIF
C
C.... Check the density, temperature etc., if they are exist
C.... in data file, to be sure they are within correct range.
      IF (NCOLS(3) .EQ. 1) THEN
        IF ((MEDATA(I,3).LE.0.0) .OR. (MEDATA(I,3).GE.1.0E7)) THEN
          WRITE(6,150) MEDATA(I,3),I
 150      FORMAT(/1X,'ERROR!  DENSITY is ',E10.3,' OUT OF RANGE!'
     >     ,/1X,'** Check the hmF2 data file at column 2, row ',I3)
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
        ENDIF
      ENDIF

      IF (NCOLS(4).EQ.1) THEN
        IF ((MEDATA(I,4).LE.400.0) .OR. (MEDATA(I,4).GE.1.0E4)) THEN
          WRITE(6,151) MEDATA(I,4),I
 151      FORMAT(/1X,'ERROR! TEMPERATURE is ',E10.3,' OUT OF RANGE!'
     >    ,/1X,'** Check the hmF2 data file at column 4, row ',I3)
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
        ENDIF
      ENDIF

      IF (NCOLS(5).EQ.1) THEN
        IF ((MEDATA(I,5).LE.400.0) .OR. (MEDATA(I,5).GE.4000)) THEN
          WRITE(6,152) MEDATA(I,5),I
 152      FORMAT(/1X,'ERROR! TEMPERATURE is ',E10.3,' OUT OF RANGE!'
     >    ,/1X,'** Check the hmF2 data file at column 5, row ',I3)
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
        ENDIF
      ENDIF

      IF (NCOLS(6).EQ.1) THEN
        IF ((MEDATA(I,6).LE.0.1) .OR. (MEDATA(I,6).GE.5.0)) THEN
          WRITE(6,153) MEDATA(I,6),I
 153      FORMAT(/1X,'ERROR! OXYGEN factor is ',E10.3,' OUT OF RANGE!'
     >    ,/1X,'** Check the hmF2 data file at column 6, row ',I3)
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
        ENDIF
      ENDIF

25    CONTINUE
26    CONTINUE      !... an error was encountered
      IF(I.LT.N) THEN
        WRITE(6,*) ' **** Bad data in Wind/hmF2 file at line = ',I
      ELSE
        WRITE(6,*) ' **** Too many lines of data in Wind/hmF2 file = '
      ENDIF
      CALL RUN_ERROR    !.. print ERROR warning in output files
      STOP
27    CONTINUE      

      NWTIM = I - 1
      IF (DEBUG) WRITE(6,*)'NWTIM = ',NWTIM
C......... Test for the times to be in ascending order
      DO 30 I = 2,NWTIM
      IF (DEBUG) WRITE(6,*)'LOOP 3 ',I
      IF (MEDATA(I,1) .LT. MEDATA(I-1,1)) THEN
          IF(JUNIT.EQ.21)WRITE(6,781)(MEDATA(I,J),J=1,NCOLS(1))
          IF(JUNIT.EQ.22)WRITE(6,782)(MEDATA(I,J),J=1,NCOLS(1))
          WRITE(6,*) '  '
          CALL RUN_ERROR    !.. print ERROR warning in output files
          STOP
      ENDIF 
C
 781   FORMAT(//3X,'ERROR! TIME NOT IN ASCENDING ORDER in NORTHern'
     >  ,/3X,'hemisphere wind/hmF2 file near the following line'
     >  ,/3X,F7.2,F8.2,1P,E10.2,0P,2F10.2)
 782   FORMAT(//3X,'ERROR! TIME NOT IN ASCENDING ORDER in SOUTHern'
     >  ,/3X,'hemisphere wind/hmF2 file near the following line'
     >  ,/3X,F7.2,F8.2,1P,E10.2,0p,2F10.2)
C
30    CONTINUE
C....... Test for UTHRS larger than the first WTIME
      IF ((MEDATA(1,1).GT.UTHRS).AND.(IABS(UTORLT).EQ.1)) THEN
         WRITE (6,391) UTHRS,MEDATA(1,1) 
 391   FORMAT(2X,'ERROR! Need winds at earlier times in wind data'
     > ,1X,'file.',/5X
     > ,'Run starts at ',F8.2,' hrs but first data time is ',F8.2)
         CALL RUN_ERROR    !.. print ERROR warning in output files
         STOP
      ENDIF
C....... Test for the total hours of run (UTHRS + RUNHRS) is less
C....... than the last WTIME
      IF ((MEDATA(NWTIM,1).LE.(UTHRS+RUNHRS)) .AND.
     >        (IABS(UTORLT).EQ.1)) THEN
         RUNHRS=MEDATA(NWTIM,1)-UTHRS-0.02
         IF(JUNIT.EQ.21) WRITE (6,392) 'NORTH',RUNHRS
         IF(JUNIT.EQ.22) WRITE (6,392) 'SOUTH',RUNHRS
 392     FORMAT(/2X,'***WARNING***  run time exceeds the wind times in '
     >   ,/2X,A5,' wind file. TSTOP has been reset to ',F10.2,' hours')
         IF(RUNHRS.LE.0.0) THEN
            WRITE(6,*) ' NOT ENOUGH VALID DATA'
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ENDIF
      ENDIF
C
C....... Make sure local times are OK
      IF (IABS(UTORLT) .EQ. 0) THEN
         CALL RUN_ERROR    !.. print ERROR warning in output files
         IF (ABS(MEDATA(1,1)) .GT. 0.001)  WRITE(6,393)
 393     FORMAT(2X,'ERROR! First local time must be zero')
         IF (ABS(MEDATA(NWTIM,1)-24.0) .GT. 0.001) WRITE(6,394)
 394     FORMAT(2X,'ERROR! last local time must be 24.0')
         IF (ABS(MEDATA(1,1)).GT.0.001 .OR.
     >      ABS(MEDATA(NWTIM,1)-24.0).GT.0.001) STOP
         MEDATA(1,1) = 0.0
         MEDATA(NWTIM,1) = 24.0
      ENDIF
      CLOSE (JUNIT)
      RETURN
      END
C
C 
C::::::::::::::::::::::::::::::: HWINPT :::::::::::::::::::::::::::::
C....... This subroutine provides the northern or southern hemisphere
C....... winds RESULT(2) from measured hmf2 (MEDATA(NWTIM,2)) for 
C....... the FLIP model given UTHRS is the universal time at which the
C....... wind will be determined, RUNHRS is the length of the run
C....... (hours). 
C....... NWTIM is the number of hmf2 times on the data file. SATK are
C....... the local times in the north or south hemispheres. The data
C....... file contains hmf2 in the northern or southern hemispheres.
C....... To calculate the HMF2 at the input UTHRS or SAT it does a 
C....... simple linear interpolation. Then it call HMFW2 to evaluate the 
C....... wind that corresponds to the HMF2. This subroutine may also provide 
C....... the density in RESULT(2), the temperature in RESULT(3), etc.,
C....... if there such data in the data file. COSDIP=cos(dip angle of the
C....... magnetic field), HPREV=hmf2 at previous time step and
C....... UPREV=wind at previous time step are provide by calling
C....... program.
C....... Programmer - D. A. Burns, June  1988.
C....... Modifications by Rama Lavu and P. Richards  26 JAN 89
C....... Dr. Phlip Richards, Kevin Zou    July 9, 93
C
      SUBROUTINE HWINPT(N,M,NWTIM,NCOLS,UTHRS,SATK,RUNHRS,COSDIP,
     >     HPREV,UPREV,MEDATA,UTORLT,RESULT)
      INTEGER  UTORLT,NWTIM,NCOLS(9)
      DOUBLE PRECISION WINDX,COSDIP,HPREV,UPREV,RESULT(M),ALPHA(25)
     >  ,RUNHRS
      REAL MEDATA(N,M), SAT,SATK,UTHRS,PHASE

      !.. DH/DU for January 1986
      DATA ALPHA/.25,.26,.27,.28,.29,.29,.30,.37,.48,.55,.58,.58,.58,.
     >  59,.57,.53,.48,.42,.38,.37,.33,.29,.27,.26,.25/

      !.. The FLIP local time is adjusted to take into account the lag 
      !.. in ionospheric response. This is related to the ALPHA factor
	PHASE= 0.00000* 0.1*ALPHA(NINT(SATK)+1)/ALPHA(1)	!.. Phase for reaction time
	!..WRITE(6,'(3F8.2)') UTHRS,SATK,PHASE

      SAT = SATK+PHASE
      IF(SAT .GE. 24.0) THEN SAT = SAT-24.0  !.. correct to 0 - 24

      !.. If UTORLT =1 UT is used in the hmF2 file, else LT
      IF (UTORLT .EQ. 1) THEN
        CALL INTWIN(N,M,NWTIM,NCOLS,UTHRS+PHASE,MEDATA,RESULT)
      ELSE IF(UTORLT .EQ. 0) THEN
        CALL INTWIN(N,M,NWTIM,NCOLS,SAT,MEDATA,RESULT)
      ELSE
         WRITE(6,*) ' Problem with UT/LT switch in data file'
         CALL RUN_ERROR    !.. print ERROR warning in output files
         STOP
      ENDIF

C..... Store the interpolation value HMF2 at required time
C..... in RESULT(1)
      RESULT(1)=RESULT(2)

C.... Now evaluate the wind that corresponds to the hmF2
      CALL HMF2W(SAT,COSDIP,RESULT(2),HPREV,UPREV,WINDX)
C
C.....Convert wind result from m/sec to cm/sec and store the wind 
C.....result in RESULT(2) and 
      RESULT(2)=100.0*WINDX
C
      RETURN
      END
C
C
C::::::::::::::::::::::::::: HMF2W :::::::::::::::::::::::::::
C..... This subroutine adjusts the neutral wind to try to reproduce
C..... the stipulated HMF2. It uses the hmF2 (HPREV) and the wind
C..... (UPREV) from the previous time step and assumes a linear
C..... relationship between the wind and hmF2. Written by P. Richards
C..... December 1990. See Miller et al. references (WITS paper 1990)
C
      SUBROUTINE HMF2W(TIME,COSDIP,HMF2,HPREV,UPREV,WIND)

      REAL TIME,DHDU
      DOUBLE PRECISION COSDIP,HMF2,HPREV,UPREV,WIND,ALPHA(25),COLFAC
     > ,CFAC,GL,DHDU_FAC
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/CFACT/DHDU_FAC,COLFAC,CFAC(9),ISCH(9)
	DATA STPWRITE/-1/

      !.. DH/DU for January 1986
      DATA ALPHA/.25,.26,.27,.28,.29,.29,.30,.37,.48,.55,.58,.58,.58,.
     >  59,.57,.53,.48,.42,.38,.37,.33,.29,.27,.26,.25/

      !.. The 1.7 factor enables collision factor scaling. 3.21 dip effect
      DHDU=-(ALPHA(NINT(TIME)+1)) * (1.7/COLFAC)
     > *(3.21*COSDIP*SIN(ACOS(COSDIP)))
       !... * ((F107A+F107)/75.0/2.0)

	DHDU=-0.1    !.. used to be DHDU_FAC
      WIND = (HMF2-HPREV)/DHDU + UPREV
      IF(STPWRITE.LT.-999)
     >  WRITE(6,90) TIME,COLFAC,F107,COSDIP,HMF2,HPREV,UPREV,DHDU,WIND
	STPWRITE=1.0

 90   FORMAT(1X,'DHDU',F9.2,1P,9E10.2)
      RETURN
      END
C::::::::::::::::::::::::::: BISPLW ::::::::::::::::::::::::::::::::::
C.... Use bisection to find the nearest time MD(PDT,1) to desired
C.....time (STA or UTHRS). Return the position value PDT. Notice that the 
C.....input data of time MD(NWTIM,1) is in ascending order.
C.....Dr. Phlip Richards, Kevin Zou at July 9, 1993
C
      SUBROUTINE BISPLW(N,M,BEGIN,NWTIM,T,MD,PDT)
C
      INTEGER PDT,FIRST,LAST,MID,BEGIN,NWTIM
      REAL MD(N,M),T
C
      FIRST=BEGIN
      LAST=NWTIM
C    
 10   CONTINUE
           MID=(FIRST+LAST)/2+1
           IF(T.LE.MD(MID,1)) THEN
              LAST=MID-1
           ELSE
              FIRST=MID
           ENDIF  
      IF(FIRST.LT.LAST) GO TO 10
      IF(T.LT.MD(FIRST,1)) THEN
	   PTD=FIRST-1
      ELSE
	   PDT=FIRST
      ENDIF
      IF(T.LE.MD(BEGIN+1,1)) PDT=BEGIN
C
      RETURN
      END
C
C::::::::::::::::::::::::::: PWINDS ::::::::::::::::::::::::
C...... This subroutine prints data winds and ExB drifts read from data files
C...... Written by P. Richards in August 18 1993, modified September 1998
      SUBROUTINE PWINDS(IWIND,SATN,SATF,SEC,TSTOP,COSDIP,ZTOPTE
     >     ,EXB_TYPE,IODDN,IVIBN2)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,SATN,SATF,ZTOPTE
      DIMENSION DATINN(6),DATINS(6)
      INTEGER NWTIMN,NWTIMS,NCOLN,NCOLS,UTOLTN,UTOLTS,HMOWNN,HMOWNS
     > ,EXB_TYPE,IWIND
      COMMON/WPAR/NWTIMN,NWTIMS,NCOLN(9),NCOLS(9),UTOLTN,UTOLTS
     >  ,HMOWNN,HMOWNS

       !.. Print information about winds only if reading file
       IF(IWIND.GT.0) THEN
          NWTIMN=-1
          NWTIMS=-1
          CALL NTHWIN(NWTIMN,SEC/3600.,SATN,TSTOP,COSDIP,3.0D2
     >               ,0.0D0,DATINN,NCOLN,UTOLTN,HMOWNN,ZTOPTE)
          CALL PRNWIN(NCOLN,UTOLTN,HMOWNN,'NORTHern',ZTOPTE)
          CALL STHWIN(NWTIMS,SEC/3600.,SATF,TSTOP,COSDIP,3.0D2
     >               ,0.0D2,DATINS,NCOLS,UTOLTS,HMOWNS,ZTOPTE)
          CALL PRNWIN(NCOLS,UTOLTS,HMOWNS,'SOUTHern',ZTOPTE)
        ENDIF
        
        !.... info on ExB drifts
        WRITE(6,*) ' '
        IF(EXB_TYPE.EQ.2) THEN
           WRITE(6,*) '   EExBB drift:: Velocities are read from file'
        ELSE IF(EXB_TYPE.EQ.-1) THEN
           WRITE(6,*) '   EExBB drift from Khazanov JGR, 1994 page 5923'
        ELSE
           WRITE(6,*) '   EExBB drift:: NONE on this run'
        ENDIF
        RETURN
        END
C
C:::::::::::::::::::::::::::: PRNWIN :::::::::::::::::::::::::::::::::
C..... This subroutine writes the information about the data in the 
C..... wind/hmf2 data file.
C..... NCOLS= number of columns in the file
C..... UTORLT indicates universial (UTORLT=1) or local (UTORLT=0) time
C..... HMORWN indicates hmF2 (HMORWN=1) or wind (HMORWN=0) in column 2
C..... HEMIS is used to indentify "southern" or "northern" hemisphere.
C..... Dr. Philip Richards, Kevin Zou   at July 29, 1993
C
      SUBROUTINE PRNWIN(NCOLS,UTORLT,HMORWN,HEMIS,ZTOPTE)
      CHARACTER*8 HEMIS
      INTEGER NCOLS(9),UTORLT,HMORWN
      REAL ZTOPTE
C
       !... Indicate if the wind data file exists. NCOLS(1)=0 if not
      IF(NCOLS(1) .LT. 2) THEN
         WRITE(6,11) HEMIS
 11      FORMAT(/1X,'-- No wind/hmF2 data file in the ',A9,
     >  ' hemisphere: Using IRI hmF2 to get equivalent wind instead --')
        RETURN
      ENDIF
C
      WRITE(6,45) HEMIS
 45   FORMAT(/1X,'-- Characteristics of the data in the '
     >  ,A,' hemisphere wind/hmF2 file --')
      !...Show whether UT or LT in data file.
      IF(UTORLT .EQ. 1) THEN
         WRITE(6,*) ' 1st column in data file has UT'
      ELSE
         WRITE(6,*) ' 1st column in data file has LT'
      ENDIF
C
      !...Show the whether hmF2 or wind in data file.
      IF(HMORWN .EQ. 1) THEN
         WRITE(6,*) ' Using HMF2 in 2nd column to get equivalent WIND '
      ELSE
         WRITE(6,*) ' Using WIND interpolated from 2nd column data'
      ENDIF
C
      !...... File contains density
      IF(NCOLS(3).EQ.1)  
     >  WRITE(6,*) ' Using NmF2 from 3rd column to normalize FLIP NmF2'
C
      !...... File contains topside electron temperature
      IF(NCOLS(4).EQ.1.AND.ZTOPTE.GT.200.0) WRITE(6,90) NINT(ZTOPTE)
 90   FORMAT('  Adjusting plasmasphere heating to match Te at'
     > ,I4,' km')
      IF(NCOLS(4).EQ.1.AND.ZTOPTE.LT.200.0) WRITE(6,*) 
     >  ' Adjusting plasmasphere heating to match Te at hmF2'

      !...... File contains ion temperature
      IF(NCOLS(5).EQ.1)
     >    WRITE(6,*) ' Using Tn/Ti from 5th column as MSIS T(exos)'
      !... using Te at hmF2 from column 6

      IF(NCOLS(6).EQ.1) WRITE(6,*)
     >  ' Multiplying MSIS O density by factor in column 6'
      RETURN
      END
C:::::::::::::::::::::::::GETINF ::::::::::::::::::::::::::::
C..... This subroutine determines the change in  MSIS Tinf that is needed
C..... to bring the FLIP NmF2 into agreement with the measurement
C..... See Richards et al [1997a,b]
      SUBROUTINE GETINF(DT,     ! time step in seconds
     >                  NMF2N,  ! Model NmF2
     >                  DNMF2N, ! Measured NmF2
     >                  NDFNOR, ! Output: Ratio of measured and modeled NmF2
     >                  TINFN,  ! Output: Adjusted Tinf
     >                  DELTIN, ! Output: Increment in Tinf
     >                  TMSIS)  ! Output: Saved value of Tinf increment 
      DOUBLE PRECISION NMF2N,DNMF2N,DT
      REAL NDFNOR,TINFN,DELTIN,TMSIS,DENRAT
         !.. The previous value is included to avoid sharp changes
         !.. NDFNOR, DELTIN = factors to adjust MSIS O, Tn
         DENRAT=NMF2N/DNMF2N
         IF(DENRAT.GE.1.0) THEN
            DELTIN=0.5*(DELTIN+1000.0*ALOG(SQRT(DENRAT))*DENRAT**2)
            NDFNOR=0.5*(NDFNOR+(TMSIS/TINFN))  ! factor to adjust [O]
         ELSE
            DELTIN=0.5*(DELTIN+1000.0*ALOG(SQRT(DENRAT))/DENRAT**2)
            NDFNOR=0.5*(NDFNOR+(TMSIS/TINFN)**1.00000)  ! factor to adjust [O]
         ENDIF
         IF(DELTIN.GT.999) DELTIN=999
         IF(DELTIN.LT.-999) DELTIN=-999
      RETURN
      END
C:::::::::::::::::::::: FLIPTN :::::::::::::::::::::::::::
C----- This subroutine allows FLIP to adjust the MSIS exospheric 
C----- temperature (TINF) based on either measured Ti or to allow FLIP to bring 
C----- NmF2 into agreement with observations
C----- Added by P. Richards 8/13/97
      SUBROUTINE FLIPTN(GLAT,TINF)
      REAL NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS
     > ,TMSN,TMSS,UPBND,LOBND,MINTN
      INTEGER NWTIMN,NWTIMS,NCOLN,NCOLS,UTOLTN,UTOLTS,HMOWNN,HMOWNS
     >   ,TPRIN,TMODON
      COMMON/NDFAC/TPRIN,NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS
     >  ,TMSN,TMSS,EXFACN,EXFACS,TMODON
      COMMON/WPAR/NWTIMN,NWTIMS,NCOLN(9),NCOLS(9),UTOLTN,UTOLTS
     >  ,HMOWNN,HMOWNS
      DATA UPBND,LOBND,MINTN /800.0,200,550/
      
      LOBND=0.25*TINF  !... max difference between MSIS when FLIP lower
      UPBND=0.5*TINF   !... max difference between MSIS when FLIP higher
      !... For printing phase with modified temps from file. Use stored values
      IF(TPRIN.EQ.-1) RETURN
      IF(TPRIN.EQ.1.AND.GLAT.GE.0.AND.TINFN.GT.400) TINF=TINFN
      IF(TPRIN.EQ.1.AND.GLAT.LT.0.AND.TINFS.GT.400) TINF=TINFS
      IF(TPRIN.EQ.1) RETURN

      !.. Test to see if MSIS needs modifying (NCOLN(5)=1 use meas Ti)
      IF(TMODON.GE.0.AND.NCOLN(5)+NCOLS(5).EQ.0) RETURN
      !.. Adjust northern or southern hemisphere based on latitude sign
      IF(GLAT.GE.0) THEN
        TMSN=TINF     ! save MSIS value for testing FLIP value
        !... reading Tn from file if NCOLN(5) .NE. 0 
        IF(NCOLN(5).EQ.1.AND.TINFN.GE.MINTN) THEN
          TINF=TINFN   
        ELSEIF(NCOLN(5).EQ.-1.AND.TINFN.LT.MINTN) THEN
          TINF=TMSN-TINFN   !.. using Delta Tn
        ELSE
          TINF=TMSN+DELTIN
          IF(TINF.LT.TMSN-LOBND) TINF=TMSN-LOBND
          IF(TINF.GT.TMSN+UPBND) TINF=TMSN+UPBND
          IF(TINF.LT.MINTN) TINF=MINTN
          TINFN=TINF
        ENDIF
      ENDIF
      !... southern hemisphere
      IF(GLAT.LT.0) THEN
        TMSS=TINF
        IF(NCOLS(5).EQ.1.AND.TINFS.GE.MINTN) THEN
          TINF=TINFS
        ELSEIF(NCOLS(5).EQ.-1.AND.TINFS.LT.MINTN) THEN
          TINF=TMSS-TINFS
        ELSE
          TINF=TMSS+DELTIS
          IF(TINF.LT.TMSS-LOBND) TINF=TMSS-LOBND
          IF(TINF.GT.TMSS+UPBND) TINF=TMSS+UPBND
          IF(TINF.LT.MINTN) TINF=MINTN
          TINFS=TINF
        ENDIF
      ENDIF
      RETURN
      END

C:::::::::::::::::::::::::::::: HM_ADJUST :::::::::::::::::::::::
C.. Adjusting hmF2 to match NmF2 
       SUBROUTINE HM_ADJUST(IW,SATN,COSDIP
     > ,HMF2N,DHMF2N,NMF2N,DNMF2N,UNNSAV,UNN)
      IMPLICIT NONE
      INTEGER IW
      REAL SATN
      DOUBLE PRECISION COSDIP,HMF2N,DHMF2N,NMF2N,DNMF2N
     > ,UNNSAV,UNN,HM_NEW,DEL_HM
      DEL_HM=((DNMF2N-NMF2N)/DNMF2N)*100.00
      IF(DEL_HM.GT.25) DEL_HM=25
      IF(DEL_HM.LT.-25) DEL_HM=-25
      HM_NEW=DHMF2N+DEL_HM
      IF(HM_NEW.LT.200) HM_NEW=200
      IF(HM_NEW.GT.500) HM_NEW=500
      CALL HMF2W(SATN,COSDIP,HM_NEW,HMF2N,UNNSAV,UNN)
      UNN=UNN*100
      IF(IW.EQ.1) WRITE(6,91) HMF2N,DHMF2N,DEL_HM,NMF2N/DNMF2N
      IF(IW.EQ.2) WRITE(6,92) HMF2N,DHMF2N,DEL_HM,NMF2N/DNMF2N
 91   FORMAT('  Varying hmF2 in North:: h_FLIP=',F5.1
     >  ,1X,'h_data=',F5.1,' dh=',F6.2,' N_FLIP/N_data=',F5.2) 
 92   FORMAT('  Varying hmF2 in South-- h_FLIP=',F5.1
     >  ,1X,'h_data=',F5.1,' dh=',F6.2,' N_FLIP/N_data=',F5.2) 
      RETURN
      END
C:::::::::::::::::::::::::::: NEWTINF :::::::::::::::::::
C.... This function returns the MSIS exospheric temperature and
C.... factors for multiplying OX and N2. Tinf is increased by
C.... using the difference between the actual Tinf (TNORM) and 
C.... the Tinf for low Ap (TLOW). The ratio of Tinf to Tlow is
C.... then used as multiplication factors for O and N2
C.... P. Richards August 2001
      SUBROUTINE NEWTINF(IYD,ALT,GLAT,SEC24,TLOW,TNORM,TINF,FACOX
     >   ,FACN2)
      IMPLICIT NONE
      INTEGER IYD,SECSAVN,SECSAVS,SECS,TPRIN,TMODON,IDAY
      REAL ALT,GLAT,SEC24,TLOW,TNORM,TINF,FACOX,FACN2,NDFNOR,NDFSOU
     >  ,TINFN,TINFS,DELTIN,DELTIS,TMSN,TMSS,EXFACN,EXFACS
     >  ,AP,DEC,ETRAN,BLON,F107,F107A,SEC
	DOUBLE PRECISION GLx
      COMMON/NDFAC/TPRIN,NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS
     >   ,TMSN,TMSS,EXFACN,EXFACS,TMODON
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      DATA SECSAVN,SECSAVS/-100,-100/

      !.. adjust Tinf - $$FLAG
      IF(TNORM.GT.TLOW) THEN
         TINF=3*TNORM-2*TLOW      !... Double the Tinf difference
         FACN2=TINF/TLOW
         FACOX=TLOW/TINF
      ELSE
         TINF=TNORM
         FACN2=1.0
         FACOX=1.0
      ENDIF
      
      !... Record the O factor for printing later
      IF(ALT.LT.300.AND.GLAT.GE.0) EXFACN=FACOX
      IF(ALT.LT.300.AND.GLAT.LT.0) EXFACS=FACOX

      !... print once every UT
      IF(GLAT.GT.0.AND.SEC.GT.SECSAVN+1) THEN
        SECSAVN=SEC
        WRITE(6,91) SEC/3600.,NINT(TLOW),NINT(TNORM),NINT(TINF),FACOX
      ENDIF

      IF(GLAT.LE.0.AND.SEC.GT.SECSAVS+1) THEN
        SECSAVS=SEC
        WRITE(6,92) SEC/3600.,NINT(TLOW),NINT(TNORM),NINT(TINF),FACOX
      ENDIF
  91  FORMAT('  MOD MSIS, NORTH:  UT=',F8.2,'  TLOW=',I5,'  TNORM='
     >    ,I5,' TNEW=',I5,'  OFAC=',F7.3)
  92  FORMAT('  MOD MSIS, SOUTH:  UT=',F8.2,'  TLOW=',I5,'  TNORM='
     >    ,I5,' TNEW=',I5,'  OFAC=',F7.3)
      RETURN
      END
C:::::::::::::::::::::::::::: VALIDU :::::::::::::::::::
C... check that wind (UN) is within range
      SUBROUTINE VALIDU(HMF2,DHMF2,UNSAV,UN)
	IMPLICIT NONE
	DOUBLE PRECISION HMF2,DHMF2,UN,UNSAV
      !IF(HMF2.LT.200.OR.DHMF2.LT.200) UN=0 !.. now use HWM in calling routine 
	!.. set a maximum or minimum wind
      IF(UN.LT.-39900)  UN=-39900
      IF(UN.GT.39900) UN=39900

      RETURN
	END
C:::::::::::::::::::::::::::: HM_OPLUSN :::::::::::::::::::
C.. This routine finds the North O+ peak height by searching down the
C.. density profile. Not in use 2003-07-28.
      SUBROUTINE HM_OPLUSN(JSTART,  !.. starting index
     >                      JSTOP,   !.. ending index 
     >                          N,    !.. Ion density array
     >                          Z,    !.. Altitude array
     >                       HMOP)   !.. Output: height of North O+ peak
      IMPLICIT NONE
      INTEGER J,JSTART,JSTOP
      DOUBLE PRECISION N(4,401),Z(401),HMOP

	HMOP=500
      DO J=JSTART-1,JSTOP+1,-1
      IF(Z(J).GT.210.AND.N(1,J).GE.N(1,J-1).AND.N(1,J).GT.N(1,J+1)) THEN 
         HMOP=Z(J)
         !..WRITE(6,90) JSTART,JSTOP,J,NINT(HMOP),N(1,J)
	ENDIF
      ENDDO
      RETURN
 90   FORMAT('  North ',4I6,1PE10.2)
      END
C:::::::::::::::::::::::::::: HM_OPLUSS :::::::::::::::::::
C.. This routine finds the O+ peak height by searching down the
C.. density profile. Not in use 2003-07-28.
      SUBROUTINE HM_OPLUSS(JSTART,  !.. starting index
     >                      JSTOP,   !.. ending index 
     >                          N,    !.. Ion density array
     >                          Z,    !.. Altitude array
     >                       HMOP)   !.. Output: height of O+ peak
      IMPLICIT NONE
      INTEGER J,JSTART,JSTOP
      DOUBLE PRECISION N(4,401),Z(401),HMOP

	HMOP=500
      DO J=JSTART+1,JSTOP-1
      IF(Z(J).GT.210.AND.N(1,J).GE.N(1,J-1).AND.N(1,J).GT.N(1,J+1)) THEN 
         HMOP=Z(J)
         !..WRITE(6,90) JSTART,JSTOP,J,NINT(HMOP),N(1,J)
	ENDIF
      ENDDO
      RETURN
 90   FORMAT('  South ',4I6,1PE10.2)
      END
