C:::::::::::::::::::::: IRI95.FOR - IRIHMF ::::::::::::::::::::::
C.... This is the interface routine for the IRI95 model
C.... Programmed by P. Richards, May 2000. This program has a
C.... test driver located on directory IRI95_FOR
      SUBROUTINE IRIHMF(DLATI,DLONGI,COV,CCIR_URSI,INDATE,
     >  DHOUR,RHMF2,RCLPD,RNMF2,FOE)
      IMPLICIT NONE
      INTEGER INDATE,RCLPD,JMAG,IYYYY,MMDD,I
      REAL DLATI,DLONGI,COV,CCIR_URSI,DHOUR,RHMF2,RNMF2,FOE
      !.. for altitude profiles. Not used here
      REAL HEIBEG,HEIEND,HEISTP,AP,DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      REAL  OUTF(11,100),OARR(38)   !... output
      LOGICAL JF(20)  !... switches for vaious features
      DOUBLE PRECISION GL
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC,GL(401)
      DATA JMAG,IYYYY,MMDD,I,HEIBEG,HEIEND,HEISTP
     >    /  0 , 1990, -075, 0, 300. , 300. , 1.0/
      DATA JF/20*.TRUE./
      !... use URSI instead of CCIR coeffs
      IF(CCIR_URSI.GE.0) JF(5)=.FALSE.

      !.. This section is for picking the IRI geophysical index. If Ap 
      !.. and F10.7 constant, use F10.7A index and convert to sunspot # 
      !.. for the IRI model. Ensures IRI index is compatible with MSIS
      IF(AP(7).LT.0) JF(17)=.FALSE.
      !.. If JF(17) false, IRI will use RZ in OARR(33), otherwise, it reads 
      !.. RZ from its own file 
      OARR(33)=(-0.728+SQRT(0.52998-3.56E-3*(63.75-F107A)))/1.78E-3
      IF(ABS(DLATI).GT.90) THEN
        IF(JF(17).EQV..FALSE.) 
     >   WRITE(6,*) ' *** Using IRI95 model with 3-month ave F10.7 ***'
        IF(JF(17).EQV..TRUE.) 
     >   WRITE(6,*) ' *** Using IRI95 model with 12-month ave F10.7 ***'
        RETURN
      ENDIF

      !... find current day =ddd if negative
      MMDD=-MOD(INDATE,1000)   
      IYYYY=(INDATE+MMDD)/1000   !... year 
 
      !... call the IRI95 model
      CALL IRIS13(JF,JMAG,DLATI,DLONGI,IYYYY,MMDD,DHOUR,
     >    HEIBEG,HEIEND,HEISTP,OUTF,OARR)

      RHMF2=OARR(2)     !... hmF2 (real)
      !... check hmF2 reasonable for FLIP winds
      IF(RHMF2.LT.180) THEN
         RHMF2=180.
         RCLPD=1
      ELSE
         RCLPD=0
      ENDIF

      RNMF2=OARR(1)    !... NmF2 (real)
      FOE=SQRT(OARR(5)/1.24E10)
      RETURN
      END

C:::::::::::::::::::::: IRIS13 (IRI95) ::::::::::::::::::::::
C IRIS13.FOR -------------------------------------------- May 1998
C 
C includes subroutines IRIS13 to compute IRI parameters for specified
C location, date, time, and altitude range and subroutine IRI_WEB to 
c computes IRI parameters for specified location, date, time and
c variable range; variable can be altitude, latitude, longitude, year, 
c month, day of month, day of year, or hour (UT or LT).
C 
C*****************************************************************
C CHANGES FROM  IRIS11.FOR  TO   IRIS12.FOR:
C    - CIRA-1986 INSTEAD OF CIRA-1972 FOR NEUTRAL TEMPERATURE
C    - 10/30/91 VNER FOR NIGHTTIME LAY-VERSION:  ABS(..)
C    - 10/30/91 XNE(..) IN CASE OF LAY-VERSION
C    - 10/30/91 CHANGE SSIN=F/T TO IIQU=0,1,2
C    - 10/30/91 Te > Ti > Tn ENFORCED IN FINAL PROFILE
C    - 10/30/91 SUB ALL NAMES WITH 6 OR MORE CHARACTERS
C    - 10/31/91 CORRECTED HF1 IN HST SEARCH:  NE(HF1)>NME
C    - 11/14/91 C1=0 IF NO F1-REGION
C    - 11/14/91 CORRECTED HHMIN AND HZ FOR LIN. APP.
C    -  1/28/92 RZ12=0 included
C    -  1/29/92 NEQV instead of NE between URSIF2 and URSIFO
C    -  5/ 1/92 CCIR and URSI input as in IRID12
C    -  9/ 2/94 Decimal month (ZMONTH) for IONCOM
C    -  9/ 2/94 Replace B0POL with B0_TAB; better annually
C    -  1/ 4/95 DY for h>hmF2
C    -  2/ 2/95 IG for foF2, topside; RZ for hmF2, B0_TAB, foF1, NmD
C    -  2/ 2/95 winter no longer exclusive for F1 occurrrence
C    -  2/ 2/95 RZ and IG included as DATA statement; smooth annual var.
C CHANGES FROM  IRIS12.FOR  TO   IRIS13.FOR:
C    - 10/26/95 include year as input and corrected MODA; nrm for zmonth
C    - 10/26/95 use TCON and month-month interpolation in foF2, hmF2
C    - 10/26/95 TCON only if date changes
C    - 11/25/95 take out logicals TOPSI, BOTTO, and BELOWE
C    - 12/ 1/95 UT_LT for (date-)correct UT<->LT conversion
C    - 12/22/95 Change ZETA cov term to cov < 180; use cov inst covsat
C    -  2/23/96 take covmax(R<150) for topside; lyear,.. for lt
C    -  3/26/96 topside: 94.5/BETA inst 94.45/..; cov -> covsat(<=188)
C    -  5/01/96 No longer DY for h>hmF2 (because of discontinuity)
C    - 12/01/96 IRIV13: HOUR for IVAR=1 (height)
C    -  4/25/97 D-region: XKK le 10 with D1 calc accordingly.
C    -  1/12/97 DS model for lower ion compoistion DY model
C    -  5/19/98 seamon=zmonth if lati>0; zmonth= ...(1.0*iday)/..
C    -  5/19/98 DY ion composition model below 300 km now DS model
C    -  5/19/98 DS model includes N+, Cl down to 75 km HNIA changed
C    -  5/28/98 User input for Rz12, foF1/NmF1, hmF1, foE/NmE, hmE; jf(17)
C    -  9/ 2/98 1 instead of 0 in MODA after UT_LT call
C    -  4/30/99 program constants moved from DATA statement into program
C    -  4/30/99 changed konsol-unit to 13 (12 is for IG_RZ).
C    -  5/29/99 the limit for IG computed from Rz12-input is 174 not 274
C    - 11/08/99 jf(18)=.true. simple UT to LT conversion, otherwise UT_LT
C    - 11/09/99 added COMMON/const1/humr,dumr also for CIRA86
C
C*****************************************************************
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C*****************************************************************
C**************** ALL-IN-ONE SUBROUTINE  *************************
C*****************************************************************
C
C
       SUBROUTINE IRIS13(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
C-----------------------------------------------------------------
C
C INPUT:  JF(1:20)      true/false switches for several options
C         JMAG          =0 geographic   = 1 geomagnetic coordinates
C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         IYYYY         Year as YYYY, e.g. 1985
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
C                          HOURS
C         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
C          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100
C
C         JF switches to turn off/on (true/false) several options
C
C       jf(i)=          .true.                  .false.
C      ---------------------------------------------------------------------
C        i=1         Ne computed            Ne not computed
C          2         Te, Ti computed        Te, Ti not computed
C          3         Ni computed            Ni not computed
C          4         B0 - Table option      B0 - Gulyaeva (1987)
C          5         foF2 - CCIR            foF2 - URSI
C          6         Ni - Standard          Ni - Danilov-Yaichnikov-Smironova
C          7         Ne - Standard Topside  Ne - IRI-79 Topside
C          8         foF2 from model        foF2 or NmF2 - user input
C          9         hmF2 from model        hmF2 or M3000F2 - user input
C         10         Te - Standard          Te - Using Te/Ne correlation
C         11         Ne - Standard Profile  Ne - Lay-function formalism
C         12         Messages are written to unit
C                             6                          12
C                    {you can also turn off all messages by setting
C                     KONSOL=1 internally}
C         13         foF1 from model        foF1 or NmF1 - user input
C         14         hmF1 from model        hmF1 - user input
C         15         foE  from model        foE or NmE - user input
C         16         hmE  from model        hmE - user input
C         17         Rz12 from file         Rz12 - user input
C         18         simple ut<->lt         using ut_lt subroutine
C         19    free
C         20    free
C      ---------------------------------------------------------------------
C
C  Depending on the jf() settings additional INPUT parameters may be required:
C
C          Setting              INPUT parameter
C       ------------------------------------------------------
C       jf(8)=.false.      OARR(1)=user input for foF2/MHz or NmF2/m-3
C       jf(9)=.false.      OARR(2)=user input for hmF2/km or M(3000)F2
C       jf(10)=.false.     OARR(15),OARR(16),OARR(17)=user input for 
C                            Ne(300km),Ne(400km),Ne(600km)/m-3; use 
C                            OARR()=-1 if one of these values is not
C                            available.
C       jf(13)=.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3 
C       jf(14)=.false.     OARR(4)=user input for hmF1/km
C       jf(15)=.false.     OARR(5)=user input for foE/MHz or NmE/m-3 
C       jf(16)=.false.     OARR(6)=user input for hmE/km
C       jf(17)=.false.     OARR(33)=user input for Rz12
C
C  OUTPUT:  OUTF(1:11,1:100)
C               OUTF(1,*)  ELECTRON DENSITY/M-3
C               OUTF(2,*)  NEUTRAL TEMPERATURE/K
C               OUTF(3,*)  ION TEMPERATURE/K
C               OUTF(4,*)  ELECTRON TEMPERATURE/K
C               OUTF(5,*)  O+ ION DENSITY/M-3
C               OUTF(6,*)  H+ ION DENSITY/M-3
C               OUTF(7,*)  HE+ ION DENSITY/M-3
C               OUTF(8,*)  O2+ ION DENSITY/M-3
C               OUTF(9,*)  NO+ ION DENSITY/M-3
C                 AND, IF JF(6)=.FALSE.:
C               OUTF(10,*)  PERCENTAGE OF CLUSTER IONS IN %
C               OUTF(11,*)  PERCENTAGE OF N+ IONS IN %
C
C            OARR(1:35)   ADDITIONAL OUTPUT PARAMETERS         
C
C              #OARR(1) = NMF2/M-3           #OARR(2) = HMF2/KM
C              #OARR(3) = NMF1/M-3           #OARR(4) = HMF1/KM
C              #OARR(5) = NME/M-3            #OARR(6) = HME/KM
C               OARR(7) = NMD/M-3             OARR(8) = HMD/KM
C               OARR(9) = HHALF/KM            OARR(10) = B0/KM
C               OARR(11) =VALLEY-BASE/M-3     OARR(12) = VALLEY-TOP/KM
C               OARR(13) = TE-PEAK/K          OARR(14) = TE-PEAK HEIGHT/KM
C              #OARR(15) = TE-MOD(300KM)     #OARR(16) = TE-MOD(400KM)/K
C              #OARR(17) = TE-MOD(600KM)      OARR(18) = TE-MOD(1400KM)/K
C               OARR(19) = TE-MOD(3000KM)     OARR(20) = TE(120KM)=TN=TI/K
C               OARR(21) = TI-MOD(430KM)      OARR(22) = X/KM, WHERE TE=TI
C               OARR(23) = SOL ZENITH ANG/DEG OARR(24) = SUN DECLINATION/DEG
C               OARR(25) = DIP/deg            OARR(26) = DIP LATITUDE/deg
C               OARR(27) = MODIFIED DIP LAT.  OARR(28) = DELA
C               OARR(29) = sunrise/dec. hours OARR(30) = sunset/dec. hours
C               OARR(31) = ISEASON (1=spring) OARR(32) = NSEASON (northern)
C              #OARR(33) = Rz12               OARR(34) = Covington Index
C               OARR(35) = B1                 oarr(36) = M(3000)F2
C     place_holders only:
C               oarr(36) = TEC/m-2            oarr(38) = TEC_top/TEC*100.
C
C      # INPUT as well as OUTPUT parameter
C      lower case:      special for IRIWeb
C-------------------------------------------------------------------
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
C***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
C***     TEMPERATURES              120 KM        3000 KM       ***
C***     ION DENSITIES             100 KM        1000 KM       ***
C*****************************************************************
C*****************************************************************
C*********            INTERNALLY                    **************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C*****************************************************************
C*****************************************************************
      INTEGER           DAYNR,DDO,DO2,SEASON,SEADAY
      REAL              LATI,LONGI,MO2,MO,MODIP,NMF2,MAGBR
      REAL              NMF1,NME,NMD,MM,MLAT,MLONG,NOBO2
c      CHARACTER FILNAM*12
      CHARACTER FILNAM*53
      DIMENSION  F(3),RIF(4),E(4),XDELS(4),DNDS(4)
      DIMENSION  FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2)
      DIMENSION  FF0N(988),XM0N(441),F2N(13,76,2),FM3N(9,49,2)
      DIMENSION  AMP(4),HXL(4),SCL(4),B0B1(5)
      DIMENSION  XSM(4),MM(5),DTI(4)
      DIMENSION  AHH(7),STTE(6),DTE(5),ATE(7),TEA(6),HOA(3),XNAR(3)
      DIMENSION  PG1O(80),PG2O(32),PG3O(80),PF1O(12),PF2O(4),PF3O(12)
      DIMENSION  HO(4),MO(5),DDO(4),HO2(2),MO2(3),DO2(2),DION(7)
      DIMENSION  OUTF(11,100),OARR(38),ARIG(3),RZAR(3)
      LOGICAL           EXT,SCHALT,NIGHT,TECON(3),sam_mon,sam_yea
      LOGICAL           F1REG,FOF2IN,HMF2IN,URSIF2,LAYVER,DY,GULB0
      LOGICAL           FOF1IN,HMF1IN,FOEIN,HMEIN,RZIN
      LOGICAL           NODEN,NOTEM,NOION,TENEOP,sam_doy
      LOGICAL           OLD79,JF(20),URSIFO
      COMMON    /BLOCK1/HMF2,NMF2,HMF1  /CONST/UMR      /const1/humr,dumr
     &          /BLOCK2/B0,B1,C1      /BLOCK3/HZ,T,HST,STR
     &          /BLOCK4/HME,NME,HEF   /BLOCK5/NIGHT,E
     &          /BLOCK6/HMD,NMD,HDX   /BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2
     &          /BLOCK8/HS,TNHS,XSM,MM,DTI,MXSM         
     &          /BLOTN/XSM1,TEXOS,TLBDH,SIGMA /BLOTE/AHH,ATE1,STTE,DTE
     &          /BLO10/BETA,ETA,DELTA,ZETA       /ARGEXP/ARGMAX
      EXTERNAL          XE1,XE2,XE3,XE4,XE5,XE6,TEDER
      DATA      B0B1  /.755566,.778596,.797332,.812928,.826146/

CTEST   save    icalls,rssn,rsat,cov,ttt
        save

        DO 7397 KI=1,11
        do 7397 kk=1,100
7397    OUTF(KI,kk)=-1.

        do 8398 kind=7,14,1
8398    oarr(kind)=-1.
        do 8378 kind=18,32,1
8378    oarr(kind)=-1.
        oarr(34)=-1.
        oarr(35)=-1.
        oarr(36)=-1.
        oarr(37)=-1.
        oarr(38)=-1.

C
C PROGRAM CONSTANTS
C
        icalls=icalls+1
        ARGMAX=88.0
        pi=ATAN(1.0)*4.
        UMR=pi/180.
        humr=pi/12.
        dumr=pi/182.5
        ALOG2=ALOG(2.)
        ALG100=ALOG(100.)
        numhei=int(abs(heiend-heibeg)/abs(heistp))+1
        if(numhei.gt.100) numhei=100
C
C Code inserted to aleviate block data problem for PC version.
C Thus avoiding DATA statement with parameters from COMMON block.
C
        HOA(1)=300.
        HOA(2)=400.
        HOA(3)=600.
        XDELS(1)=5.
        XDELS(2)=5.
        XDELS(3)=5.
        XDELS(4)=10.
        DNDS(1)=.016
        DNDS(2)=.01
        DNDS(3)=.016
        DNDS(4)=.016
        DDO(1)=9
        DDO(2)=5
        DDO(3)=5
        DDO(4)=25
        DO2(1)=5
        DO2(2)=5
        XNAR(1)=0.0
        XNAR(2)=0.0
        XNAR(3)=0.0
        AHH(1)=120.
        AHH(2)=0.
        AHH(3)=300.
        AHH(4)=400.
        AHH(5)=600.
        AHH(6)=1400.
        AHH(7)=3000.
        DTE(1)=5.
        DTE(2)=5.
        DTE(3)=10.
        DTE(4)=20.
        DTE(5)=20.
        DTI(1)=10.
        DTI(2)=10.
        DTI(3)=20.
        DTI(4)=20.
C
C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
C AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
C IUURSI=UNIT NUMBER FOR CCIR COEFFICIENTS ........................
c set konsol=1 if you do not want the konsol information
c
      MONITO=6
      IUCCIR=10
      IUURSI=11     !... PGR unit added by P. Richards in May 2000
c web ---- special for web version
c web      KONSOL=1
      KONSOL=1
        if(.not.jf(12)) konsol=13
c
c selection of density and ion composition options ..................
c

      NODEN=(.not.jf(1))
      NOTEM=(.not.jf(2))
      NOION=(.not.jf(3))
      DY=(.not.jf(6))
      LAYVER=(.not.jf(11))
      OLD79=(.not.jf(7))
      GULB0=(.not.jf(4))
c
c rz12 input option ....................................................
c
      RZIN=(.not.jf(17))
       IF(RZIN) THEN
          ARZIN=OARR(33)
        else
          oarr(33)=-1.
          ENDIF
c
c F2 peak density ....................................................
c
      FOF2IN=(.not.jf(8))
       IF(FOF2IN) THEN
          AFOF2=OARR(1)
          IF(AFOF2.GT.100.) AFOF2=SQRT(AFOF2/1.24E10)
        else
          oarr(1)=-1.
          ENDIF
      URSIF2=(.not.jf(5))
c
c F2 peak altitude ..................................................
c
      HMF2IN=(.not.jf(9))
       IF(HMF2IN) then
                AHMF2=OARR(2)
        else
                oarr(2)=-1.
        endif
c
c F1 peak density ....................................................
c
      FOF1IN=(.not.jf(13))
       IF(FOF1IN) THEN
          AFOF1=OARR(3)
          IF(AFOF1.GT.100.) AFOF1=SQRT(AFOF1/1.24E10)
        else
          oarr(3)=-1.
          ENDIF
c
c F1 peak altitude ..................................................
c
      HMF1IN=(.not.jf(14))
       IF(HMF1IN) then
                AHMF1=OARR(4)
                if(.not.layver.and.(konsol.gt.1)) write(konsol,1939)
1939  format(' *Ne* User input of hmF1 is only possible for the LAY-',
     &          'version')
        else
                oarr(4)=-1.
        endif
c
c E peak density ....................................................
c
      FOEIN=(.not.jf(15))
       IF(FOEIN) THEN
          AFOE=OARR(5)
          IF(AFOE.GT.100.) AFOE=SQRT(AFOE/1.24E10)
        else
          oarr(5)=-1.
          ENDIF
c
c E peak altitude ..................................................
c
      HMEIN=(.not.jf(16))
       IF(HMEIN) then
                AHME=OARR(6)
        else
                oarr(6)=-1.
        endif
c
C TE-NE MODEL OPTION ..............................................
C
      TENEOP=(.not.jf(10))
        IF(TENEOP) THEN
           DO 8154 JXNAR=1,3
              XNAR(JXNAR)=OARR(JXNAR+14)
              TECON(JXNAR)=.FALSE. 
8154          IF(XNAR(JXNAR).GT.0.) TECON(JXNAR)=.TRUE. 
        else
           oarr(3)=-1.
           oarr(4)=-1.
           oarr(5)=-1.
           ENDIF
c
c lists the selected options before starting the table
c

      if(icalls.gt.1.or.konsol.eq.1) goto 8201
        write(*,*) '*** IRI parameters are being calculated ***'
      if(NODEN) goto 2889
        if(LAYVER) write(*,*) 'Ne, E-F: The LAY-Version is ',
     &          'prelimenary. Erroneous profile features can occur.'
        if(GULB0) write(*,*) 'Ne, B0: Bottomside thickness is ',
     &          'obtained with Gulyaeva-1987 model.'
        if(OLD79) write(*,*) 'Ne: Using IRI-79. Correction',
     &          ' of equatorial topside is not included.'
        if(FOF2IN) then
                write(*,*) 'Ne, foF2/NmF2: provided by user.'
                goto 2889
                endif
        if(URSIF2) then
                write(*,*) 'Ne, foF2: URSI model is used.'
        else
                write(*,*) 'Ne, foF2: CCIR model is used.'
        endif
2889    continue
        if(HMF2IN) write(*,*) 'Ne, hmF2/M3000F2: provided by user.'
        if(fof1in) write(*,*) 'Ne, foF1/NmF1: provided by user.'
        if(HMF1IN) write(*,*) 'Ne, hmF1: prvided by user.'
        if(foein) write(*,*) 'Ne, foE/NmE: provided by user.'
        if(HMEIN) write(*,*) 'Ne, hmE: prvided by user.'
        if((.not.NOION).and.(DY)) 
     &    write(*,*) 'Ion Com.: Using Danilov et al. 1985/95.'
        if((.not.NOTEM).and.(TENEOP))
     &    write(*,*) 'Te: Temperature-density correlation is used.'
8201    continue

C
C CALCULATION OF GEOG. OR GEOM. COORDINATES IN DEG....................
C CALCULATION OF MAGNETIC INCLINATION (DIP), DECLINATION (DEC)........
C   DIP LATITUDE (MAGBR) AND MODIFIED DIP (MODIP). ALL IN DEGREE......
C
        IF(JMAG.GT.0) THEN
           MLAT=ALATI
           MLONG=ALONG
        ELSE
           LATI=ALATI
           LONGI=ALONG
        ENDIF
        CALL GGM(JMAG,LONGI,LATI,MLONG,MLAT)
        ABSLAT=ABS(LATI)
        CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP)
        MAGBR=ATAN(0.5*TAN(DIP*UMR))/UMR
        ABSMLT=ABS(MLAT)
        ABSMDP=ABS(MODIP)
        ABSMBR=ABS(MAGBR)
C
C CALCULATION OF DAY OF YEAR AND SUN DECLINATION......................
C CALCULATION OF UT/LT AND RELATED YEAR, MONTH, DAYNRs ...............
C CALCULATION OF (UT-)SEASON (SUMMER=2, WINTER=4).....................
C
        iyear=iyyyy
        if(iyear.lt.100) iyear=iyear+1900
        if(MMDD.lt.0) then
                DAYNR=-MMDD
                call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
        else
                MONTH=MMDD/100
                IDAY=MMDD-MONTH*100
                call MODA(0,iyear,MONTH,IDAY,DAYNR,nrdaym)
        endif
c
c lyear,lmonth,lday,ldaynre,lnrday related to LT
c
        lyear=iyear
        lmonth=month
        lday=iday
        ldaynr=daynr
        lnrday=nrdaym

      IF(DHOUR.le.24.0) goto 2619
        UT=DHOUR-25.
        iytmp=iyear
        idtmp=daynr
        if(jf(18)) then
                hour=ut+longi/15.
                if(hour.gt.24.) hour=hour-24.
        else
                call ut_lt(0,ut,hour,longi,iytmp,idtmp)
                if(idtmp.ne.ldaynr) then
                        lyear=iytmp
                        ldaynr=idtmp
                        call MODA(1,lyear,LMONTH,LDAY,LDAYNR,lnrday)
                        endif
        endif
        goto 2629

2619  HOUR=DHOUR
        iytmp=lyear
        idtmp=ldaynr
        if(jf(18)) then
                ut=hour-longi/15.
                if(ut.lt.0) ut=ut+24.
        else
                call ut_lt(1,ut,hour,longi,iytmp,idtmp)
                if(idtmp.ne.daynr) then
                        iyear=iytmp
                        daynr=idtmp
                        call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
                        endif
        endif
2629  zmonth = lmonth + (lday*1.)/lnrday

      SEASON=INT((DAYNR+45.0)/92.0)
      IF(SEASON.LT.1) SEASON=4
      NSEASN=SEASON
      seaday=daynr
      seamon=zmonth
      IF(LATI.GT.0.0) GOTO 5592
        SEASON=SEASON-2
        IF(SEASON.LT.1) SEASON=SEASON+4
        seamon=zmonth+6.
        if(seamon.ge.13.0) seamon=seamon-12.
        seaday=daynr+183
        if(seaday.gt.366) seaday=seaday-366
C
C CALCULATION OF MEAN F10.7CM SOLAR RADIO FLUX (COV)................
C CALCULATION OF RESTRICTED SOLAR ACTIVITIES (RSAT,COVSAT)..............
C

5592    sam_mon=(month.eq.montho)
        sam_yea=(iyear.eq.iyearo)
        sam_doy=(daynr.eq.idayno)
        if(sam_yea.and.sam_doy) goto 2910
                call tcon(iyear,month,iday,daynr,rzar,arig,ttt,nmonth)
                if(nmonth.lt.0) goto 3330
                if(RZIN) then
                  rrr = arzin
                  rzar(1) = rrr
                  rzar(2) = rrr
                  rzar(3) = rrr
                  zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
                  if(zi.gt.174.0) zi=174.0
                  arig(1) = zi
                  arig(2) = zi
                  arig(3) = zi
                  endif
        rssn=rzar(3)
        rsat=arig(3)
        COV=63.75+RSSN*(0.728+RSSN*0.00089)
        rlimit=rsat
        COVSAT=63.75+rlimit*(0.728+rlimit*0.00089)
        if(covsat.gt.188.) covsat=188

C
C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG).........................
C NOON VALUE (XHINON).................................................
C

2910    CALL SOCO(ldaynr,HOUR,LATI,LONGI,SUNDEC,XHI,SAX,SUX)
        CALL SOCO(ldaynr,12.0,LATI,LONGI,SUNDE1,XHINON,SAXNON,SUXNON)

        NIGHT=.FALSE.
        if(abs(sax).gt.25.0) then
                if(sax.lt.0.0) NIGHT=.TRUE.
                goto 1334
                endif
        if(SAX.le.SUX) goto 1386
                if((hour.gt.sux).and.(hour.lt.sax)) night=.true.
        goto 1334
1386            IF((HOUR.GT.SUX).OR.(HOUR.LT.SAX)) NIGHT=.TRUE.

C
C CALCULATION OF ELECTRON DENSITY PARAMETERS................
C

1334  HNEA=65.
      IF(NIGHT) HNEA=80.
      HNEE=2000.
      IF(NODEN) GOTO 4933
      DELA=4.32
      IF(ABSMDP.GE.18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
      DELL=1+EXP(-(ABSLAT-20.)/10.)
C!!!!!!! F-REGION PARAMETERS AND E-PEAK !!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(FOEIN) THEN
          FOE=AFOE
        ELSE
          FOE=FOEEDI(COV,XHI,XHINON,ABSLAT)
        ENDIF
        NME=1.24E10*FOE*FOE

        IF(HMEIN) THEN
          HME=AHMF2
        ELSE
          HME=110.0
        ENDIF
C
C READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH ............
C
      IF((FOF2IN).AND.(HMF2IN)) GOTO 501
      IF(URSIF2.NEQV.URSIFO) GOTO 7797
      IF(sam_mon.AND.(nmonth.EQ.nmono).and.sam_yea) GOTO 4292
      IF(sam_mon) GOTO 4293
c
c the program expects the coefficients files in ASCII format; if you
C want to use the binary version of the coefficients, please use the
C the statements that are commented-out below and comment-out the
C ASCII-related statements.
c
7797    URSIFO=URSIF2
c
        !... changed by P. Richards in May 2000 so that all
        !... the data can be read from one file.
        !... read down to correct month
        DO IMON=1,MONTH
           READ(IUCCIR,4689) F2,FM3
        ENDDO
4689    FORMAT(1X,4E15.8)
        CALL PC_CLOSE(IUCCIR)
C
C then URSI if chosen ....................................
C
        if(URSIF2) then
           !... changed by P. Richards in May 2000 so that all
           !... the data can be read from one file.
           !... read down to correct month
           DO IMON=1,MONTH
             READ(IUURSI,4689) F2
           ENDDO
           CALL PC_CLOSE(IUURSI)
        endif

C
C READ CCIR AND URSI COEFFICIENT SET FOR NMONTH, i.e. previous 
c month if day is less than 15 and following month otherwise 
C

4293    continue

c
c first CCIR ..............................................
c
        !... changed by P. Richards in May 2000 so that all
        !... the data can be read from one file.
        !... read down to correct month
        DO IMON=1,NMONTH
           READ(IUCCIR,4689) F2N,FM3N
        ENDDO
        CALL PC_CLOSE(IUCCIR)

C
C then URSI if chosen .....................................
C
        if(URSIF2) then
           !... changed by P. Richards in May 2000 so that all
           !... the data can be read from one file.
           !... read down to correct month
           DO IMON=1,NMONTH
             READ(IUURSI,4689) F2N
           ENDDO
           CALL PC_CLOSE(IUURSI)
        endif

        nmono=nmonth
        MONTHO=MONTH
        iyearo=iyear
        idayno=daynr
        GOTO 4291
        
8448    WRITE(MONITO,8449) FILNAM
8449    FORMAT(1X////,
     &    ' The file ',A30,'is not in your directory.')
        GOTO 3330
C
C LINEAR INTERPOLATION IN SOLAR ACTIVITY. RSAT used for foF2
C

4291    RR2=ARIG(1)/100.
        RR2N=ARIG(2)/100.
        RR1=1.-RR2
        RR1N=1.-RR2N
        DO 20 I=1,76
        DO 20 J=1,13
        K=J+13*(I-1)
        FF0N(K)=F2N(J,I,1)*RR1N+F2N(J,I,2)*RR2N
20      FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2

        RR2=RZAR(1)/100.
        RR2N=RZAR(2)/100.
        RR1=1.-RR2
        RR1N=1.-RR2N
        DO 30 I=1,49
        DO 30 J=1,9
        K=J+9*(I-1)
        XM0N(K)=FM3N(J,I,1)*RR1N+FM3N(J,I,2)*RR2N
30      XM0(K)=FM3(J,I,1)*RR1+FM3(J,I,2)*RR2

4292  zfof2  =  FOUT(MODIP,LATI,LONGI,UT,FF0)
      fof2n  =  FOUT(MODIP,LATI,LONGI,UT,FF0N)
      zm3000 = XMOUT(MODIP,LATI,LONGI,UT,XM0)
      xm300n = XMOUT(MODIP,LATI,LONGI,UT,XM0N)

        midm=15
        if(month.eq.2) midm=14
        if (iday.lt.midm) then
                yfof2 = fof2n + ttt * (zfof2-fof2n)
                xm3000= xm300n+ ttt * (zm3000-xm300n)
        else
                yfof2 = zfof2 + ttt * (fof2n-zfof2)
                xm3000= zm3000+ ttt * (xm300n-zm3000)
        endif
501     IF(FOF2IN) THEN
          FOF2=AFOF2
        ELSE
          FOF2=YFOF2
        ENDIF
        NMF2=1.24E10*FOF2*FOF2

        IF(HMF2IN) THEN
          HMF2=AHMF2
          IF(AHMF2.LT.50.0) HMF2=HMF2ED(MAGBR,RSSN,FOF2/FOE,AHMF2)
        ELSE
          HMF2=HMF2ED(MAGBR,RSSN,FOF2/FOE,XM3000)
        ENDIF
c
c topside profile parameters .............................
c
      COS2=COS(MLAT*UMR)
      COS2=COS2*COS2
        FLU=(COVSAT-40.0)/30.0
corr    FLU=(COV-40.0)/30.0
corr    flueta=188.
corr    flumax=(flueta-40.)/30.0
      IF(OLD79) then
        ETA1=-0.0070305*COS2
      else
        EX=EXP(-MLAT/15.)
        EX1=EX+1
        EPIN=4.*EX/(EX1*EX1)
        ETA1=-0.02*EPIN
      endif
      ETA=0.058798+ETA1+FLU*(-0.014065+0.0069724*COS2)+
     &(0.0024287+0.0042810*COS2-0.00015280*FOF2)*FOF2
        fluu=flu
corr    if(fluu.gt.flumax) fluu=flumax
      ZETA=0.078922-0.0046702*COS2-0.019132*FLUU+0.0076545*FLU*COS2+
     &(0.0032513+0.0060290*COS2-0.00020872*FOF2)*FOF2
      BETA=-128.03+20.253*COS2+FLU*(-8.0755-0.65896*COS2)+(0.44041
     &+0.71458*COS2-0.042966*FOF2)*FOF2
        Z=EXP(94.5/BETA)
corr    Z=EXP(94.45/BETA)
      Z1=Z+1
      Z2=Z/(BETA*Z1*Z1)
      DELTA=(ETA/Z1-ZETA/2.0)/(ETA*Z2+ZETA/400.0)
c
c bottomside profile parameters .............................
C
1501    HMF1=HMF2
        HZ=HMF2
        HEF=HME
        B1=3.0
C!!!!!!! INTERPOLATION FOR B0 OUT OF ARRAY B0F !!!!!!!!!!!!!!!!!!!!!
        if(GULB0) then
          call ROGUL(SEADAY,XHI,SEAX,GRAT)
ctest
          if(NIGHT) GRAT=0.91-HMF2/4000.
          B0CNEW=HMF2*(1.-GRAT)
          B0=B0CNEW/B0B1(1)
        else
          B0 = B0_TAB(HOUR,SAX,SUX,NSEASN,RSSN,MODIP)
        endif
C!!!!!!! F1-REGION PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      F1REG=.FALSE.
      HMF1=0.
      PNMF1=0.
      C1=0.  
      IF(NIGHT) GOTO 150
        IF(FOF1IN) THEN
          FOF1=AFOF1
        ELSE
          FOF1=FOF1ED(ABSMBR,RSSN,XHI)
        ENDIF
        IF(FOF1.LT.1.E-3) GOTO 150
          F1REG=.TRUE.
          C1=.09+.11/DELA
          PNMF1=1.24E10*FOF1*FOF1
150   NMF1=PNMF1
C!!!!!!! PARAMETER FOR E AND VALLEY-REGION !!!!!!!!!!!!!!!!!!!!!
      XDEL=XDELS(SEASON)/DELA
      DNDHBR=DNDS(SEASON)/DELA
      HDEEP=HPOL(HOUR,10.5/DELA,28.,SAX,SUX,1.,1.)
      WIDTH=HPOL(HOUR,17.8/DELA,45.+22./DELA,SAX,SUX,1.,1.)
      DEPTH=HPOL(HOUR,XDEL,81.,SAX,SUX,1.,1.)
      DLNDH=HPOL(HOUR,DNDHBR,.06,SAX,SUX,1.,1.)
      IF(DEPTH.LT.1.0) GOTO 600
        IF(NIGHT) DEPTH=-DEPTH
        CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
        IF(.NOT.EXT) GOTO 667
        if(konsol.gt.1) WRITE(KONSOL,650)
650     FORMAT(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
600   WIDTH=.0
667   HEF=HME+WIDTH
      VNER = (1. - ABS(DEPTH) / 100.) * NME
c
c Parameters below E  .............................
c

2727  NMD=XMDED(XHI,RSSN,4.0E8)
      HMD=HPOL(HOUR,81.0,88.0,SAX,SUX,1.,1.)
      F(1)=HPOL(HOUR,0.02+0.03/DELA,0.05,SAX,SUX,1.,1.)
      F(2)=HPOL(HOUR,4.6,4.5,SAX,SUX,1.,1.)
      F(3)=HPOL(HOUR,-11.5,-4.0,SAX,SUX,1.,1.)
      FP1=F(1)
      FP2=-FP1*FP1/2.0
      FP30=(-F(2)*FP2-FP1+1.0/F(2))/(F(2)*F(2))
      FP3U=(-F(3)*FP2-FP1-1.0/F(3))/(F(3)*F(3))
c indermediate region between D and E region; parameters xkk
c and d1 are found such that the function reaches hdx/xdx/dxdh
      HDX=HMD+F(2)
      X=HDX-HMD
      XDX=NMD*EXP(X*(FP1+X*(FP2+X*FP30)))
      DXDX=XDX*(FP1+X*(2.0*FP2+X*3.0*FP30))
      X=HME-HDX
      XKK=-DXDX*X/(XDX*ALOG(XDX/NME))
c if exponent xkk is larger than xkkmax, then xkk will be set to xkkmax and
c d1 will be determined such that the point hdx/xdx is reached; derivative
c is no longer continuous.
        xkkmax=5.
        if(xkk.gt.xkkmax) then
                xkk=xkkmax
                d1=-alog(xdx/nme)/(x**xkk)
        else
                D1=DXDX/(XDX*XKK*X**(XKK-1.0))
        endif
C
C SEARCH FOR HMF1 ..................................................
C

2726  if(LAYVER) goto 6153

924   IF(.not.F1REG) GOTO 380
        XE2H=XE2(HEF)
        CALL REGFA1(HEF,HMF2,XE2H,NMF2,0.001,NMF1,XE2,SCHALT,HMF1)
        IF(.not.SCHALT) GOTO 380
        if(konsol.gt.1) WRITE(KONSOL,11)
11      FORMAT(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2')
        IREGFA=1
c
c change B1 and try again ..........................................
c

9244    IF(B1.GT.4.5) GOTO (7398,8922) IREGFA
                B1=B1+0.5
                if(konsol.gt.1) WRITE(KONSOL,902) B1-0.5,B1
902     FORMAT(6X,'CORR.: B1(OLD)=',F4.1,' B1(NEW)=',F4.1)
                IF(GULB0) then
                        ib1=int(b1*2.-5.)
                        B0=B0CNEW/b0b1(ib1)
                        endif
                GOTO 924
c
c omit F1 feature ....................................................
c

7398  if(konsol.gt.1) WRITE(KONSOL,9269)
9269  FORMAT(1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
        HMF1=0.
        NMF1=0.
        C1=0.0
        B1=3.
        F1REG=.FALSE.

C
C SEARCH FOR HST [NE3(HST)=NME] ..........................................
C

380     RRRR=0.5
        IF(F1REG) then
                hf1=hmf1
                xf1=nmf1
                GOTO 3972
                ENDIF
        RATHH=0.5
3973    hf1=hef+(hmf2-hef)*RATHH
        xf1=xe3(hf1)
        IF(XF1.LT.NME) THEN
                RATHH=RATHH+.1
                GOTO 3973
                ENDIF
3972    h=hf1
        deh=10.
        XXMIN=XF1
        HHMIN=HF1
3895    h=h-deh
        if(h.lt.HEF) then
          h=h+2*deh
          deh=deh/10.
          if(deh.lt.1.) goto 3885
          endif
        XE3H=XE3(h)
        IF(XE3H.LT.XXMIN) then
          XXMIN=XE3H
          HHMIN=h
          endif
        if(XE3H.gt.NME) goto 3895
      CALL REGFA1(h,HF1,XE3H,XF1,0.001,NME,XE3,SCHALT,HST)
        STR=HST
        IF(.not.SCHALT) GOTO 360
3885  if(konsol.gt.1) WRITE(KONSOL,100)
100   FORMAT(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
        IREGFA=2
        IF(XXMIN/NME.LT.1.3) GOTO 9244
c
c assume linear interpolation between HZ and HEF ..................
c

8922    HZ=HHMIN+(HF1-HHMIN)*RRRR
        XNEHZ=XE3(HZ)
        if(xnehz-nme.lt.0.001) then
          RRRR=RRRR+.1
          GOTO 8922
          endif
        if(konsol.gt.1) WRITE(KONSOL,901) HZ,HEF
901     FORMAT(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,
     &          ' AND HEF=',F5.1)
        T=(XNEHZ-NME)/(HZ-HEF)
        HST=-333.
        GOTO 4933
c
c calculate HZ, D and T ............................................
c

360     HZ=(HST+HF1)/2.0
        D=HZ-HST
        T=D*D/(HZ-HEF-D)
        GOTO 4933
C
C LAY-functions for middle ionosphere
C

6153    IF(HMF1IN) THEN
          HMF1M=AHMF1
        ELSE
          HMF1M=165.+0.6428*XHI
        ENDIF
        HHALF = GRAT * HMF2
        HV1R = HME + WIDTH
        HV2R = HME + HDEEP
        HHMF2 = HMF2
        CALL INILAY(NIGHT,NMF2,NMF1,NME,VNER,HHMF2,HMF1M,HME,
     &                  HV1R,HV2R,HHALF,HXL,SCL,AMP,IIQU)
        IF((IIQU.EQ.1).and.(konsol.gt.1)) WRITE(KONSOL,7733)
7733    FORMAT('*NE* LAY amplitudes found with 2nd choice of HXL(1).')
        IF((IIQU.EQ.2).and.(konsol.gt.1)) WRITE(KONSOL,7722)
7722    FORMAT('*NE* LAY amplitudes could not be found.')

C---------- CALCULATION OF NEUTRAL TEMPERATURE PARAMETER-------

4933  HTA=120.0
      HTE=3000.0
        IF(NOTEM) GOTO 240
      SEC=UT*3600.
      CALL CIRA86(DAYNR,SEC,LATI,LONGI,HOUR,COV,TEXOS,TN120,SIGMA)
        IF(HOUR.NE.0.0) THEN
                iyz=lyear
                idz=ldaynr
                if(jf(18)) then
                        secni=(24.-longi/15)*3600.
                else
                        call ut_lt(1,utni,0.0,longi,iyz,idz)
                        SECNI=utni*3600.
                endif
      CALL CIRA86(DAYNR,SECNI,LATI,LONGI,0.,COV,TEXNI,TN1NI,SIGNI)
        ELSE
      TEXNI=TEXOS
      TN1NI=TN120
      SIGNI=SIGMA
        ENDIF
      TLBDH=TEXOS-TN120
      TLBDN=TEXNI-TN1NI
C
C--------- CALCULATION OF ELECTRON TEMPERATURE PARAMETER--------
C
881   CONTINUE

C !!!!!!!!!! TE(120KM)=TN(120KM) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ATE(1)=TN120

C !!!!!!!!!! TE-MAXIMUM (JICAMARCA,ARECIBO) !!!!!!!!!!!!!!!!!!!!
      HMAXD=60.*EXP(-(MLAT/22.41)**2)+210.
      HMAXN=150.
      AHH(2)=HPOL(HOUR,HMAXD,HMAXN,SAX,SUX,1.,1.)
      TMAXD=800.*EXP(-(MLAT/33.)**2)+1500.
      TMAXN=TN(HMAXN,TEXNI,TLBDN,SIGNI)+20
      ATE(2)=HPOL(HOUR,TMAXD,TMAXN,SAX,SUX,1.,1.)

C !!!!!!!!!! TE(300,400KM)=TE-AE-C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! TE(1400,3000KM)=TE-ISIS !!!!!!!!!!!!!!!!!!!!!!!!!!!
        DIPLAT=MAGBR
      CALL TEBA(DIPLAT,HOUR,NSEASN,TEA)
      ATE(3)=TEA(1)
      ATE(4)=TEA(2)
      ATE(6)=TEA(3)
      ATE(7)=TEA(4)

C !!!!!!!!!! TE(600KM)=TE-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ETT=EXP(-MLAT/11.35)
      TET=2900.-5600.*ETT/((ETT+1)**2.)
      TEN=839.+1161./(1.+EXP(-(ABSMLT-45.)/5.))
      ATE(5)=HPOL(HOUR,TET,TEN,SAX,SUX,1.5,1.5)

C !!!!!!!!!! OPTION TO USE TE-NE-RELATION !!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! AT 300, 400 OR 600 KM  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(TENEOP) THEN
        DO 3395 I=1,3
3395      IF(TECON(I)) ATE(I+2)=TEDE(HOA(I),XNAR(I),-COV)
        ENDIF
C !!!!!!!!!! TE'S ARE CORRECTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! ALSO TE > TN ENFORCED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TNAHH2=TN(AHH(2),TEXOS,TLBDH,SIGMA)
      IF(ATE(2).LT.TNAHH2) ATE(2)=TNAHH2
      STTE1=(ATE(2)-ATE(1))/(AHH(2)-AHH(1))
      DO 1901 I=2,6
       TNAHHI=TN(AHH(I+1),TEXOS,TLBDH,SIGMA)
       IF(ATE(I+1).LT.TNAHHI) ATE(I+1)=TNAHHI
       STTE2=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
       ATE(I)=ATE(I)-(STTE2-STTE1)*DTE(I-1)*ALOG2
1901  STTE1=STTE2
C !!!!!!!!!! GRADIENTS ARE CALCULATED WITH !!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!! CORRECTED REGION BOUNDARIES !!!!!!!!!!!!!!!!!!!!!!
      DO 1902 I=1,6
1902  STTE(I)=(ATE(I+1)-ATE(I))/(AHH(I+1)-AHH(I))
      ATE1=ATE(1)
887   CONTINUE
C
C------------ CALCULATION OF ION TEMPERATURE PARAMETERS--------
C
C !!!!!!!!!! TI(430KM,DAY)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!!!
      XSM1=430.0
      XSM(1)=XSM1
      Z1=EXP(-0.09*MLAT)
      Z2=Z1+1.
      TID1 = 1240.0 - 1400.0 * Z1 / ( Z2 * Z2 )
      MM(2)=HPOL(HOUR,3.0,0.0,SAX,SUX,1.,1.)
C !!!!!!!!!!  TI < TE   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TED1=TEA(6)+30.
        IF(TID1.GT.TED1) TID1=TED1

C !!!!!!!!!! TI(430KM,NIGHT)=TI-AEROS !!!!!!!!!!!!!!!!!!!!!!!!!
      Z1=ABSMLT
      Z2=Z1*(0.47+Z1*0.024)*UMR
      Z3=COS(Z2)
      TIN1=1200.0-300.0*SIGN(1.0,Z3)*SQRT(ABS(Z3))
C !!!!!!!!!! TN < TI < TE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TEN1=TEA(5)
        TNN1=TN(XSM1,TEXNI,TLBDN,SIGNI)
        IF(TEN1.LT.TNN1) TEN1=TNN1
        IF(TIN1.GT.TEN1) TIN1=TEN1
        IF(TIN1.LT.TNN1) TIN1=TNN1

C !!!!!!!!!! TI(430KM,LT) FROM STEP FUNCTION !!!!!!!!!!!!!!!!!!
        TI1=TIN1  
        IF(TID1.GT.TIN1) TI1=HPOL(HOUR,TID1,TIN1,SAX,SUX,1.,1.)

C !!!!!!!!!! TANGENT ON TN DETERMINES HS !!!!!!!!!!!!!!!!!!!!!!
        TI13=TEDER(130.)
        TI50=TEDER(500.)
      CALL REGFA1(130.0,500.0,TI13,TI50,0.01,TI1,TEDER,SCHALT,HS)
      IF(SCHALT) HS=200.
      TNHS=TN(HS,TEXOS,TLBDH,SIGMA)
      MM(1)=DTNDH(HS,TEXOS,TLBDH,SIGMA)
      IF(SCHALT) MM(1)=(TI1-TNHS)/(XSM1-HS)
      MXSM=2

C !!!!!!!!!! XTETI ALTITTUDE WHERE TE=TI !!!!!!!!!!!!!!!!!!!!!!
2391    XTTS=500.
        X=500.
2390    X=X+XTTS
        IF(X.GE.AHH(7)) GOTO 240
        TEX=ELTE(X)
        TIX=TI(X)
        IF(TIX.LT.TEX) GOTO 2390
        X=X-XTTS
        XTTS=XTTS/10.
        IF(XTTS.GT.0.1) GOTO 2390
        XTETI=X+XTTS*5.

C !!!!!!!!!! TI=TE ABOVE XTETI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MXSM=3
        MM(3)=STTE(6)
        XSM(2)=XTETI
        IF(XTETI.GT.AHH(6)) GOTO 240
        MXSM=4
        MM(3)=STTE(5)
        MM(4)=STTE(6)
        XSM(3)=AHH(6)
        IF(XTETI.GT.AHH(5)) GOTO 240
        MXSM=5
        DTI(1)=5.
        DTI(2)=5.
        MM(3)=STTE(4)
        MM(4)=STTE(5)
        MM(5)=STTE(6)
        XSM(3)=AHH(5)
        XSM(4)=AHH(6)
C
C CALCULATION OF ION DENSITY PARAMETER..................
C
240   IF(NOION) GOTO 141
      HNIA=100.
        if(DY) hnia=75
      HNIE=2000.
C
C INPUT OF THE ION DENSITY PARAMETER ARRAYS PF1O,PF2O AND PF3O......
C
      RIF(1)=2.
      IF(ABSLAT.LT.30.0) RIF(1)=1.
      RIF(2)=2.
      IF(COV.LT.100.0) RIF(2)=1.
      RIF(3)=SEASON
      IF(SEASON.EQ.1) RIF(3)=3.
      RIF(4)=1.
      IF(NIGHT) RIF(4)=2.
      CALL KOEFP1(PG1O)
      CALL KOEFP2(PG2O)
      CALL KOEFP3(PG3O)
      CALL SUFE(PG1O,RIF,12,PF1O)
      CALL SUFE(PG2O,RIF, 4,PF2O)
      CALL SUFE(PG3O,RIF,12,PF3O)
c
c calculate O+ profile parameters
c
      IF(ABS(XHI).LE.90.0) THEN
        ZZZ1=COS(XHI*UMR)
      ELSE
        ZZZ1=0.0
      ENDIF
      msumo=4
      RDOMAX=100.0
      MO(1)=EPSTEP(PF1O(1),PF1O(2),PF1O(3),PF1O(4),ZZZ1)
      MO(2)=EPSTEP(PF1O(5),PF1O(6),PF1O(7),PF1O(8),ZZZ1)
      MO(3)=0.0
      HO(1)=EPSTEP(PF1O(9),PF1O(10),PF1O(11),PF1O(12),ZZZ1)
      HO(2)=290.0
      IF((RIF(2).EQ.2.).AND.(RIF(3).EQ.2.)) HO(2)=237.0
      HO(4)=PF2O(1)
        ho05=pf2o(4)
      MO(4)=PF2O(2)
      MO(5)=PF2O(3)
c
c adjust gradient MO(4) of O+ profile segment above F peak
c
7100    HO(3)=(ALG100-MO(5)*(HO(4)-ho05))/MO(4)+HO(4)
        IF(HO(3).LE.HO(2)+20.) THEN
                MO(4)=MO(4)-0.001
                GOTO 7100
                endif
        hfixo=(ho(2)+ho(3))/2.
c
c find height H0O of maximum O+ relative density
c
      DELX=5.0
      X=HO(2)
      YMAXX=0.0
7102  X=X+DELX
      Y=RPID(X,HFIXO,RDOMAX,msumo,MO,DDO,HO)
      IF(Y.LE.YMAXX) then
        if(delx.le.0.1) GOTO 7104
        x=x-delx
        delx=delx/5.
      ELSE
        YMAXX=Y
      ENDIF
      GOTO 7102 
7104  H0O=X-DELX/2.
7101    if(y.lt.100.0) goto 7103
          rdomax=rdomax-0.01
        y=rpid(h0o,hfixo,rdomax,msumo,mo,ddo,ho)
        goto 7101
7103    yo2h0o=100.-y
        yoh0o=y
c
c calculate parameters for O2+ profile
c
        hfixo2  = pf3o(1)
        rdo2mx = pf3o(2) 
      DO 7105 L=1,2
                I = L * 2
                HO2(L)=PF3O(1+I)+PF3O(2+I)*ZZZ1
7105            MO2(L+1)=PF3O(7+I)+PF3O(8+I)*ZZZ1
      MO2(1)=PF3O(7)+PF3O(8)*ZZZ1
        if(hfixo2.gt.ho2(1)) then
           ymo2z=mo2(2)
        else
           ymo2z=mo2(1)
        endif
        aldo21=alog(rdo2mx)+ymo2z*(ho2(1)-hfixo2)
        hfixo2=(ho2(2)+ho2(1))/2.
        rdo2mx=exp(aldo21+mo2(2)*(hfixo2-ho2(1)))
c
c make sure that rd(O2+) is less or equal 100-rd(O+) at O+ maximum
c
7106  Y=RPID(H0O,hfixo2,rdo2mx,2,MO2,DO2,HO2)
      IF(Y.GT.yo2h0o) then
        MO2(3)=MO2(3)-0.02
        GOTO 7106
        endif
c
C use ratio of NO+ to O2+ density at O+ maximum to calculate
c NO+ density above the O+ maximum (H0O)
c
      IF(y.LT.1.) then
        NOBO2=0.0
      ELSE
        NOBO2= (yo2h0o-y)/y
      ENDIF
C
C CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
C
141   IF(.NOT.F1REG) HMF1=HZ

        height=heibeg
        kk=1

300   IF(NODEN) GOTO 330
      IF((HEIGHT.GT.HNEE).OR.(HEIGHT.LT.HNEA)) GOTO 330
        IF(LAYVER) THEN
          ELEDE=-9.
          IF(IIQU.LT.2) ELEDE=XEN(HEIGHT,HMF2,NMF2,HME,4,HXL,SCL,AMP)
        ELSE
          ELEDE=XE(HEIGHT)
        ENDIF
        OUTF(1,kk)=ELEDE
330   IF(NOTEM) GOTO 7108
      IF((HEIGHT.GT.HTE).OR.(HEIGHT.LT.HTA)) GOTO 7108
        TNH=TN(HEIGHT,TEXOS,TLBDH,SIGMA)
        TIH=TNH
        IF(HEIGHT.GE.HS) TIH=TI(HEIGHT)
        TEH=ELTE(HEIGHT)
        IF(TIH.LT.TNH) TIH=TNH
        IF(TEH.LT.TIH) TEH=TIH
        OUTF(2,kk)=TNH
        OUTF(3,kk)=TIH
        OUTF(4,kk)=TEH
7108  IF(NOION) GOTO 7118
      IF((HEIGHT.GT.HNIE).OR.(HEIGHT.LT.HNIA)) GOTO 7118
        if(DY) then
c      call IONCOM(HEIGHT,XHI,LATI,COV,seamon,DION)
      call ioncom_new(height,xhi,lati,cov,seamon,dion)
      ROX=DION(1)
      RHX=DION(2)
      RNX=DION(3)
      RHEX=DION(4)
      RNOX=DION(5)
      RO2X=DION(6)
      RCLUST=DION(7)
        else
      ROX=RPID(HEIGHT,HFIXO,RDOMAX,msumo,MO,DDO,HO)
      RO2X=RPID(HEIGHT,HFIXO2,rdo2mx,2,MO2,DO2,HO2)
      CALL RDHHE(HEIGHT,H0O,ROX,RO2X,NOBO2,10.,RHX,RHEX)
      RNOX=RDNO(HEIGHT,H0O,RO2X,ROX,NOBO2)
      RNX=-1.
      RCLUST=-1.
        endif
c
c ion densities are given in percent of total electron density;
c to get ion densities in cm-3 use the following statement
c        xnorm=elede/100.
        xnorm=1.
      OUTF(5,kk)=ROX*xnorm
      OUTF(6,kk)=RHX*xnorm
      OUTF(7,kk)=RHEX*xnorm
      OUTF(8,kk)=RO2X*xnorm
      OUTF(9,kk)=RNOX*xnorm
      OUTF(10,kk)=RCLUST*xnorm
      OUTF(11,kk)=RNX*xnorm

7118    height=height+heistp
        kk=kk+1
        if(kk.le.numhei) goto 300
C
C ADDITIONAL PARAMETER FIELD OARR
C
        IF(NODEN) GOTO 6192
      OARR(1)=NMF2
      OARR(2)=HMF2
      OARR(3)=NMF1
      OARR(4)=HMF1
      OARR(5)=NME
      OARR(6)=HME
      OARR(7)=NMD
      OARR(8)=HMD
      OARR(9)=HHALF
      OARR(10)=B0
      OARR(11)=VNER
      OARR(12)=HEF
6192    IF(NOTEM) GOTO 6092
      OARR(13)=ATE(2)
      OARR(14)=AHH(2)
      OARR(15)=ATE(3)
      OARR(16)=ATE(4)
      OARR(17)=ATE(5)
      OARR(18)=ATE(6)
      OARR(19)=ATE(7)
      OARR(20)=ATE(1)
      OARR(21)=TI1
      OARR(22)=XTETI
6092  OARR(23)=XHI
      OARR(24)=SUNDEC
      OARR(25)=DIP
      OARR(26)=MAGBR
      OARR(27)=MODIP
      OARR(28)=DELA
      OARR(29)=SAX
      OARR(30)=SUX
      OARR(31)=SEASON
      OARR(32)=NSEASN
      OARR(33)=rssn
      OARR(34)=COV
      OARR(35)=B1
      OARR(36)=xm3000
3330  CONTINUE
      RETURN
      END
c
        subroutine iri_web(jmag,jf,alati,along,iyyyy,mmdd,iut,dhour,
     &          height,h_tec_max,ivar,vbeg,vend,vstp,a,b)
c----------------------------------------------------------------
c changes:
c       11/16/99 jf(20) instead of jf(17)
c
c----------------------------------------------------------------
c input:        jmag,alati,along,iyyyy,mmdd,dhour  see IRIS12
c               height  height in km
c               h_tec_max  =0 no TEC otherwise upper boundary for integral
c               iut     =1 for UT       =0 for LT
c               ivar    =1      altitude
c                       =2,3    latitude,longitude
c                       =4,5,6  year,month,day
c                       =7      day of year
c                       =8      hour (UT or LT)
c               vbeg,vend,vstp  variable range (begin,end,step)
c output:       a       similar to outf in IRIS12
c               b       similar to oarr in IRIS12
c
c
c               numstp  number of steps; maximal 100
c----------------------------------------------------------------
        dimension outf(11,100),oar(38),oarr(38),a(11,100),b(38,100)
        dimension       xvar(8)
        logical         jf(20)


        numstp=int((vend-vbeg)/vstp)+1
        if(numstp.gt.100) numstp=100

            do 6249 i=1,38
6249            oar(i)=b(i,1) 

        if(ivar.eq.1) then
            do 1249 i=1,38
1249            oarr(i)=oar(i) 
            xhour=dhour+iut*25.
            call IRIS13(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,XHOUR,
     &                  VBEG,VEND,VSTP,a,OARR)
            if(h_tec_max.gt.50.) then 
                call iri_tec (50.,h_tec_max,2,tec,tect,tecb)
                oarr(37)=tec
                oarr(38)=tect
                endif
            do 1111 i=1,38
1111            b(i,1)=oarr(i)
            return
            endif

        if(height.le.0.0) height=100
        xvar(2)=alati
        xvar(3)=along
        xvar(4)=iyyyy
        xvar(5)=mmdd/100
        xvar(6)=mmdd-xvar(5)*100
        xvar(7)=abs(mmdd*1.)
        xvar(8)=dhour

        xvar(ivar)=vbeg

        alati=xvar(2)
        along=xvar(3)
        iyyyy=int(xvar(4))
        if(ivar.eq.7) then
                mmdd=-int(vbeg)
        else
                mmdd=int(xvar(5)*100+xvar(6))
        endif
        dhour=xvar(8)+iut*25.

        do 1 i=1,numstp
                do 1349 iii=1,38
1349                    oarr(iii)=oar(iii) 
                call IRIS13(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &                  height,height,1.,OUTF,OARR)
                if(h_tec_max.gt.50.) then
                        call iri_tec (50.,h_tec_max,2,tec,tect,tecb)
                        oarr(37)=tec
                        oarr(38)=tect
                        endif
                do 2 ii=1,11
2                       a(ii,i)=outf(ii,1)
                do 2222 ii=1,38
2222                    b(ii,i)=oarr(ii)

                xvar(ivar)=xvar(ivar)+vstp

                alati=xvar(2)
                along=xvar(3)
                iyyyy=int(xvar(4))
                if(ivar.eq.7) then
                        mmdd=-xvar(7)
                else
                        mmdd=int(xvar(5)*100+xvar(6))
                endif
                dhour=xvar(8)+iut*25.
1       continue
        return
        end
C IRIT13.FOR ------------------------------------ Oct 20, 1995 
C
C contains IRIT13, IONCORR, IRI_TEC
C
C ----------------------------------------------------------------
C Corrections
C
C  3/25/96 jmag in IRIT13 as input
C  8/31/97 hu=hr(i+1) i=6 out of bounds condition corrected
C  9/16/98 JF(17) added to input parameters; OUTF(11,50->100)
C 11/15/99 JF(20) instead of JF(17)
C
C ----------------------------------------------------------------
C
C
      subroutine IRIT13(ALATI,ALONG,jmag,jf,iy,md,hour,hbeg,hend,
     &                          tec,tecb,tect)
C-----------------------------------------------------------------
c Program for numerical integration of IRI-94 profiles from h=100km
C to h=alth. 
C       
C  INPUT:  ALATI,ALONG  LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C          jmag         =0 geographic   =1 geomagnetic coordinates
C          jf(1:20)     =.true./.false. flags; explained in IRIS13.FOR
C          iy,md        date as yyyy and mmdd (or -ddd)
C          hour         decimal hours LT (or UT+25)
c          hbeg,hend    upper and lower integration limits in km
C 
C  OUTPUT: TEC          Total Electron Content in M-2
C          tecb,tect    percentage of bottomside and topside content
C------------------------------------------------------------------
C  Changes:
C       11/08/99 jf(17) to jf(20)
C  
C------------------------------------------------------------------

        dimension       outf(11,100),oarr(35)
        logical         jf(20)
c       data            jf/20*.true./

ctest
c        save

c
c  select various options and choices for IRI-94
c

c       jf(2)=.false.           ! no temperatures
c       jf(3)=.false.           ! no ion composition
Ctest   jf(4)=.false.           ! Gulyaeva-B0
c       jf(5)=.false.           ! URSI-88 for foF2
C       jf(12)=.false.          ! konsol output to file (unit=12)

        tec = -111.
        tect= -111.
        tecb= -111.

c
c initialize IRI parameter in COMMON blocks
c
        abeg=hbeg
        aend=hend
        astp=hend-hbeg
        call IRIS13(JF,JMAG,ALATI,ALONG,IY,MD,HOUR,
     &          abeg,aend,astp,OUTF,OARR)

c
c  calculate total electron content (TEC) in m-2 using the
c  stepsize selection 2 (highest accuracy)
c

        call iri_tec (hbeg,hend,2,tec,tect,tecb)
ctest   call iri_tec (hbeg,hend,0,tec,tect,tecb)
c
c  to get ionospheric correction for frequency f/Hz use the following
C  formula   ioncorr/m = 40.3 * TEC / (f*f)
c
        return
        end
c
c
        real function ioncorr(tec,f)
c-------------------------------------------------------------------
c computes ionospheric correction IONCORR (in m) for given vertical
c ionospheric electron content TEC (in m-2) and frequency f (in Hz)
c-------------------------------------------------------------------
        ioncorr = 40.3 * tec / (f*f)
        return
        end
c
c
        subroutine iri_tec (hstart,hend,istep,tectot,tectop,tecbot)
C-------------------------------------------------------------------
C subroutine to compute the total ionospheric content
C   INPUT:      
C       hstart  altitude (in km) where integration should start
C       hend    altitude (in km) where integration should end
C       istep   =0 [fast, but higher uncertainty <5%]
C               =1 [standard, recommended]
C               =2 [stepsize of 1km throughout; best TEC, but needs time]
C   OUTPUT:
C       tectot  total ionospheric content in tec-units (10^16 m^-2)
C       tectop  topside content (in %)
C       tecbot  bottomside content (in %)
C
C   The different stepsizes for the numerical integration are defined as
C   follows (h1=100km, h2=hmF2-10km, h3=hmF2+10km, h4=hmF2+150km, h5=
C   hmF2+250km):
C       istep   h1-h2   h2-h3   h3-h4   h4-h5   h5-hend
C       0       2.0km   1.0km   2.5km   exponential approximation
C       1       2.0km   1.0km   2.5km   10.0km  30.0km
C       2       1.0km   0.5km   1.0km   1.0km   1.0km   
C
C--------------------------------------------------------------------

        logical         expo
        dimension       step(5),hr(6)
        common  /block1/hmf2,xnmf2,hmf1

ctest   
        save

        expo = .false.
        numstep = 5
        xnorm = xnmf2/1000.

        hr(1) = 100.
        hr(2) = hmf2-10.
        hr(3) = hmf2+10.
        hr(4) = hmf2+150.
        hr(5) = hmf2+250.
        hr(6) = hend

        if (istep.eq.0) then 
                step(1)=2.0
                step(2)=1.0
                step(3)=2.5
                step(4)=5.
                if (hend.gt.hr(5)) expo=.true.
                endif

        if (istep.eq.1) then
                step(1)=2.0
                step(2)=1.0
                step(3)=2.5
                step(4)=10.0
                step(5)=30.0
                endif

        if (istep.eq.2) then
                step(1)=1.0
                step(2)=0.5
                step(3)=1.0
                step(4)=1.0
                step(5)=1.0
                endif

        sumtop = 0.0
        sumbot = 0.0
C
C find the starting point for the integration
C

        i=0
        ia=1
3       i=i+1
        h=hr(i)
        if(hstart.gt.h) then
                hr(i)=hstart
                ia=i
                goto 3
                endif
C
C start the numerical integration
C
        i=ia
        h=hr(i)
        hu=hr(i+1)
        delx = step(i)
1       h = h + delx
        hx = h - delx/2.
        hh = h
        if (h.ge.hu) then
                delx = hu - h + delx
                hx = hu - delx/2.
                YNE = XE(hx)
                if((hx.gt.hmf2).and.(yne.gt.xnmf2)) yne=xnmf2
                yyy = yne * delx / xnorm
                i=i+1
                h = hr(i)
                hu = hend
                if(i.lt.5) hu = hr(i+1)
                delx = step(i)
        else
                YNE = XE(hx)
                yyy = yne * delx / xnorm
        endif
        if (hx.le.hmf2) then
                sumbot = sumbot + yyy
        else
                sumtop = sumtop + yyy
        endif
        if (expo.and.(hh.ge.hr(4))) goto 5
        if (hh.lt.hend) goto 1

        zzz = sumtop + sumbot
        tectop = sumtop / zzz * 100.
        tecbot = sumbot / zzz * 100.
        tectot = zzz * xnmf2    
        return

5       num_step = 3
        hei_top = hr(4)
        hei_end = hend
        top_end = hei_end - hei_top
        del_hei = top_end / num_step
        xntop = xe(hei_end)/xnmf2

        if(xntop.gt.0.9999) then
                ss_t = top_end  
                goto 2345
                endif

        hei_2 = hei_top
        hei_3 = hei_2 + del_hei
        hei_4 = hei_3 + del_hei
        hei_5 = hei_end

        hss = top_end / 4.
C       hss = 360.
        xkk = exp ( - top_end / hss ) - 1.
        x_2 = hei_2
        x_3 =hei_top-hss*alog(xkk*(hei_3 - hei_top)/top_end + 1.) 
        x_4 =hei_top-hss*alog(xkk*(hei_4 - hei_top)/top_end + 1.)
        x_5 = hei_end

        ed_2 = xe(x_2)/xnmf2
          if(ed_2.gt.1.) ed_2=1.
        ed_3 = xe(x_3)/xnmf2
          if(ed_3.gt.1.) ed_3=1.
        ed_4 = xe(x_4)/xnmf2
          if(ed_4.gt.1.) ed_4=1.
        ed_5 = xntop
        if(ed_3.eq.ed_2) then
         ss_2 = ed_3 * (x_3 - x_2)
        else
         ss_2=( ed_3 - ed_2 ) * ( x_3 - x_2 ) / alog ( ed_3 / ed_2 )
        endif
        if(ed_4.eq.ed_3) then
         ss_3 = ed_4 * (x_4 - x_3)
        else
         ss_3=( ed_4 - ed_3 ) * ( x_4 - x_3 ) / alog ( ed_4 / ed_3 )
        endif
        if(ed_5.eq.ed_4) then
         ss_4 = ed_5 * (x_5 - x_4)
        else
         ss_4=( ed_5 - ed_4 ) * ( x_5 - x_4 ) / alog ( ed_5 / ed_4 )
        endif

        ss_t = ss_2 + ss_3 + ss_4 

2345    sumtop = sumtop + ss_t * 1000.
        
        zzz = sumtop + sumbot
        tectop = sumtop / zzz * 100.
        tecbot = sumbot / zzz * 100.
        tectot = zzz * xnmf2

      RETURN
      END

C IRIF13.FOR 
C**************************************************************  
c changes from IRIFU9 to IRIF10:
c       SOCO for solar zenith angle 
c       ACOS and ASIN argument forced to be within -1 / +1
c       EPSTEIN functions corrected for large arguments
C**************************************************************  
c changes from IRIF10 to IRIF11: 
c       LAY subroutines introduced
c       TEBA corrected for 1400 km
C**************************************************************  
c changes from IRIF11 to IRIF12:
C       Neutral temperature subroutines now in CIRA86.FOR 
C       TEDER changed
C       All names with 6 or more characters replaced 
C       10/29/91 XEN: 10^ in loop, instead of at the end
C       1/21/93 B0_TAB instead of B0POL
C       9/22/94 Alleviate underflow condition in IONCOM exp()
C**************************************************************  
c changes from IRIF12 to IRIF13:
C        9/18/95 MODA: add leap year and number of days in month
C        9/29/95 replace F2out with FOUT and XMOUT.
C       10/ 5/95 add TN and DTNDH; earlier in CIRA86.FOR
C       10/ 6/95 add TCON for reading indices
C       10/20/95 MODA: IN=1 MONTH=IMO
C       10/20/95 TCON: now includes RZ interpolation
C       11/05/95 IONCOM->IONCO1, added IONCOM_new, IONCO2
C       11/05/95 LSTID added for strom-time updating
C       11/06/95 ROGUL: transition 20. instead of 15.
C       12/01/95 add UT_LT for (date-)correct UT<->LT conversion
C       01/16/96 TCON: add IMST to SAVE statement
C       02/02/96 ROGUL: 15. reinstated
C       02/07/96 UT_LT: ddd, dddend integer, no leap year 2000
C       03/15/96 ZERO: finding delta for topside
C       03/18/96 UT_LT: mode=1, change of year
C       12/09/96 since 2000 is leap, delete y/100*100 condition
C       04/25/97 XMDED: minimal value also daytime
C       05/18/98 TCON: changes to IG_RZ (update date); -R = Cov
C       05/19/98 Replaced IONCO2&APROK; HEI,XHI in IONCOM_NEW
C       10/01/98 added INITIALIZE
C       04/30/99 MODA: reset bb(2)=28
C       11/08/99 avoid negative x value in function XE2. Set x=0.
C       11/09/99 added COMMON/const1/humr,dumr also for CIRA86
C       02/01/99 EXIT statement in APROK replaced with GOTO
C
C**************************************************************  
C********** INTERNATIONAL REFERENCE IONOSPHERE ****************  
C**************************************************************  
C****************  FUNCTIONS,SUBROUTINES  *********************
C**************************************************************
C** IMPORTANT!! INITIALIZE (needs to be called before using 
C**                 subroutines or functions)
C** NE:         XE1,ZERO,DXE1N,XE2,XE3,XE4,XE5,XE6,XE
C** TE/TI:      TEBA,SPHARM,ELTE,TEDE,TI,TEDER,TN,DTNDH
C** NI:         RPID,RDHHE,RDNO,KOEFP1,KOEFP2,KOEFP3,SUFE
C** PEAKS:      FOUT,XMOUT,HMF2ED,FOF1ED,FOEEDI,XMDED,GAMMA1
C** MAG. FIELD: GGM,FIELDG
C** FUNCTIONS:  REGFA1,TAL
C** TIME:       SOCO,HPOL,MODA,UT_LT
C** INTERPOL.:  B0POL,B0_TAB
C** EPSTEIN:    RLAY,D1LAY,D2LAY,EPTR,EPST,EPSTEP,EPLA
C** LAY:        XE2TO5,XEN,ROGUL,VALGUL,LNGLSN,LSKNM,INILAY
C** NI-new:     IONCOM,IONCOM_NEW,IONCO1,IONCO2,APROK
C** INDICES:    TCON
C** Updating:   LSTID
C**************************************************************  
C  
C**************************************************************  
C***  -------------------ADDRESSES------------------------  ***
C***  I  PROF. K. RAWER             DR. D. BILITZA       I  ***
C***  I  HERRENSTR. 43              GSFC CODE 933        I  ***
C***  I  7801 MARCH 1               GREENBELT MD 20771   I  ***
C***  I  F.R.G.                     USA                  I  ***
C***  ----------------------------------------------------  ***
C**************************************************************  
C**************************************************************  
C
C
        subroutine initialize
      COMMON    /CONST/UMR      /ARGEXP/ARGMAX  /const1/humr,dumr
        ARGMAX=88.0
        pi=atan(1.0)*4.
        UMR=pi/180.
        humr=pi/12.
        dumr = pi / 182.5
        return 
        end
c
C        
C*************************************************************   
C*************** ELECTRON DENSITY ****************************   
C*************************************************************   
C
C
      FUNCTION XE1(H)    
c----------------------------------------------------------------
C REPRESENTING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE   
C (H=HMF2....1000 KM) BY HARMONIZED BENT-MODEL ADMITTING 
C VARIABILITY OFGLOBAL PARAMETER ETA,ZETA,BETA,DELTA WITH        
C GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY     
C (SEE MAIN PROGRAM).    
C [REF.:K.RAWER,S.RAMAKRISHNAN,1978]     
c----------------------------------------------------------------
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1
     &          /BLO10/         BETA,ETA,DELTA,ZETA
     &          /ARGEXP/        ARGMAX
                    
      DXDH = (1000.-HMF2)/700.
        x0 = 300. - delta
      xmx0 = (H-HMF2)/DXDH
         x = xmx0 + x0
        eptr1 = eptr(x,beta,394.5) - eptr(x0,beta,394.5)
        eptr2 = eptr(x,100.,300.0) - eptr(x0,100.,300.0) 

      y = BETA * ETA * eptr1 + ZETA * (100. * eptr2 - xmx0)

        Y = y * dxdh
        if(abs(Y).gt.argmax) Y = sign(argmax,Y)
      XE1 = XNMF2 * EXP(-Y)                             
      RETURN          
      END             
C 
C
        REAL FUNCTION ZERO(DELTA)
C FOR A PEAK AT X0 THE FUNCTION ZERO HAS TO BE EQUAL TO 0.
        COMMON  /BLO10/         BETA,ETA,DEL,ZETA
     &          /ARGEXP/        ARGMAX
        arg1=delta/100.
        if (abs(arg1).lt.argmax) then
                z1=1./(1.+exp(arg1))
        else if (arg1.lt.0) then
                z1=1.
        else
                z1=0.
        endif
        arg1=(delta+94.5)/beta
        if (abs(arg1).lt.argmax) then
                z2=1./(1.+exp(arg1))
        else if (arg1.lt.0) then
                z2=1.
        else
                z2=0.
        endif
        zero=zeta*(1.-z1) - eta*z2
        return
        end
C
C
      FUNCTION DXE1N(H)                            
C LOGARITHMIC DERIVATIVE OF FUNCTION XE1 (KM-1).   
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1
     &          /BLO10/         BETA,ETA,DELTA,ZETA                    

        x0 = 300. - delta
      X=(H-HMF2)/(1000.0-HMF2)*700.0 + x0
        epst2 = epst(x,100.0,300.0)
        epst1 = epst(x,beta ,394.5)
      DXE1N = - ETA * epst1 + ZETA * (1. - epst2)             
      RETURN          
      END             
C
C
      REAL FUNCTION XE2(H)                         
C ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).                   
      COMMON    /BLOCK1/HMF2,XNMF2,HMF1
     &          /BLOCK2/B0,B1,C1        /ARGEXP/ARGMAX
      X=(HMF2-H)/B0
        if(x.le.0.0) x=0.0
        z=x**b1
        if(z.gt.argmax) z=argmax
      XE2=XNMF2*EXP(-z)/COSH(X)                 
      RETURN          
      END             
C
C
      REAL FUNCTION XE3(H)                         
C ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1). 
      COMMON    /BLOCK1/        HMF2,XNMF2,HMF1
     &          /BLOCK2/        B0,B1,C1
      XE3=XE2(H)+XNMF2*C1*SQRT(ABS(HMF1-H)/B0)      
      RETURN          
      END             
C
C
      REAL FUNCTION XE4(H)                         
C ELECTRON DENSITY FOR THE INDERMEDIUM REGION (HEF..HZ).                        
      COMMON    /BLOCK3/        HZ,T,HST,STR
     &          /BLOCK4/        HME,XNME,HEF
      IF(HST.LT.0.) GOTO 100                       
      XE4=XE3(HZ+T/2.0-SIGN(1.0,T)*SQRT(T*(HZ-H+T/4.0)))                        
      RETURN          
100   XE4=XNME+T*(H-HEF)                            
      RETURN          
      END             
C
C
      REAL FUNCTION XE5(H)                         
C ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).   
      LOGICAL NIGHT   
      COMMON    /BLOCK4/        HME,XNME,HEF
     &          /BLOCK5/        NIGHT,E(4)                    
      T3=H-HME        
      T1=T3*T3*(E(1)+T3*(E(2)+T3*(E(3)+T3*E(4))))  
      IF(NIGHT) GOTO 100                           
      XE5=XNME*(1+T1)  
      RETURN          
100   XE5=XNME*EXP(T1)                              
      RETURN          
      END             
C
C
      REAL FUNCTION XE6(H)                         
C ELECTRON DENSITY FOR THE D REGION (HA...HME).    
      COMMON    /BLOCK4/        HME,XNME,HEF
     &          /BLOCK6/        HMD,XNMD,HDX
     &          /BLOCK7/        D1,XKK,FP30,FP3U,FP1,FP2    
      IF(H.GT.HDX) GOTO 100                        
      Z=H-HMD         
      FP3=FP3U        
      IF(Z.GT.0.0) FP3=FP30                        
      XE6=XNMD*EXP(Z*(FP1+Z*(FP2+Z*FP3)))           
      RETURN          
100   Z=HME-H         
      XE6=XNME*EXP(-D1*Z**XKK)
      RETURN          
      END             
C
C
      REAL FUNCTION XE(H)                          
C ELECTRON DENSITY BEETWEEN HA(KM) AND 1000 KM     
C SUMMARIZING PROCEDURES  NE1....6;                
      COMMON    /BLOCK1/HMF2,XNMF2,HMF1         
     &          /BLOCK3/HZ,T,HST,STR
     &          /BLOCK4/HME,XNME,HEF
      IF(H.LT.HMF2) GOTO 100                       
      XE=XE1(H)     
      RETURN          
100   IF(H.LT.HMF1) GOTO 300                       
      XE=XE2(H)       
      RETURN          
300   IF(H.LT.HZ) GOTO 400                         
      XE=XE3(H)       
      RETURN          
400   IF(H.LT.HEF) GOTO 500                        
      XE=XE4(H)       
      RETURN          
500   IF(H.LT.HME) GOTO 600                        
      XE=XE5(H)       
      RETURN          
600   XE=XE6(H)       
      RETURN          
      END             
C                     
C**********************************************************                     
C***************** ELECTRON TEMPERATURE ********************                    
C**********************************************************                     
C
      SUBROUTINE TEBA(DIPL,SLT,NS,TE) 
C CALCULATES ELECTRON TEMPERATURES TE(1) TO TE(4) AT ALTITUDES
C 300, 400, 1400 AND 3000 KM FOR DIP-LATITUDE DIPL/DEG AND 
C LOCAL SOLAR TIME SLT/H USING THE BRACE-THEIS-MODELS (J. ATMOS.
C TERR. PHYS. 43, 1317, 1981); NS IS SEASON IN NORTHERN
C HEMISOHERE: IS=1 SPRING, IS=2 SUMMER ....
C ALSO CALCULATED ARE THE TEMPERATURES AT 400 KM ALTITUDE FOR
C MIDNIGHT (TE(5)) AND NOON (TE(6)).   
      DIMENSION C(4,2,81),A(82),TE(6)
      COMMON    /CONST/UMR      /const1/humr,dumr
      DATA (C(1,1,J),J=1,81)/                      
     &.3100E1,-.3215E-2,.2440E+0,-.4613E-3,-.1711E-1,.2605E-1,                  
     &-.9546E-1,.1794E-1,.1270E-1,.2791E-1,.1536E-1,-.6629E-2,                  
     &-.3616E-2,.1229E-1,.4147E-3,.1447E-2,-.4453E-3,-.1853,                    
     &-.1245E-1,-.3675E-1,.4965E-2,.5460E-2,.8117E-2,-.1002E-1,                 
     &.5466E-3,-.3087E-1,-.3435E-2,-.1107E-3,.2199E-2,.4115E-3,                 
     &.6061E-3,.2916E-3,-.6584E-1,.4729E-2,-.1523E-2,.6689E-3,                  
     &.1031E-2,.5398E-3,-.1924E-2,-.4565E-1,.7244E-2,-.8543E-4,                 
     &.1052E-2,-.6696E-3,-.7492E-3,.4405E-1,.3047E-2,.2858E-2,                  
     &-.1465E-3,.1195E-2,-.1024E-3,.4582E-1,.8749E-3,.3011E-3,                  
     &.4473E-3,-.2782E-3,.4911E-1,-.1016E-1,.27E-2,-.9304E-3,                   
     &-.1202E-2,.2210E-1,.2566E-2,-.122E-3,.3987E-3,-.5744E-1,                  
     &.4408E-2,-.3497E-2,.83E-3,-.3536E-1,-.8813E-2,.2423E-2,                   
     &-.2994E-1,-.1929E-2,-.5268E-3,-.2228E-1,.3385E-2,                         
     &.413E-1,.4876E-2,.2692E-1,.1684E-2/          
      DATA (C(1,2,J),J=1,81)/.313654E1,.6796E-2,.181413,.8564E-1,               
     &-.32856E-1,-.3508E-2,-.1438E-1,-.2454E-1,.2745E-2,.5284E-1,               
     &.1136E-1,-.1956E-1,-.5805E-2,.2801E-2,-.1211E-2,.4127E-2,                 
     &.2909E-2,-.25751,-.37915E-2,-.136E-1,-.13225E-1,.1202E-1,                 
     &.1256E-1,-.12165E-1,.1326E-1,-.7123E-1,.5793E-3,.1537E-2,                 
     &.6914E-2,-.4173E-2,.1052E-3,-.5765E-3,-.4041E-1,-.1752E-2,                
     &-.542E-2,-.684E-2,.8921E-3,-.2228E-2,.1428E-2,.6635E-2,-.48045E-2,        
     &-.1659E-2,-.9341E-3,.223E-3,-.9995E-3,.4285E-1,-.5211E-3,                 
     &-.3293E-2,.179E-2,.6435E-3,-.1891E-3,.3844E-1,.359E-2,-.8139E-3,          
     &-.1996E-2,.2398E-3,.2938E-1,.761E-2,.347655E-2,.1707E-2,.2769E-3,         
     &-.157E-1,.983E-3,-.6532E-3,.929E-4,-.2506E-1,.4681E-2,.1461E-2,           
     &-.3757E-5,-.9728E-2,.2315E-2,.6377E-3,-.1705E-1,.2767E-2,                 
     &-.6992E-3,-.115E-1,-.1644E-2,.3355E-2,-.4326E-2,.2035E-1,.2985E-1/        
      DATA (C(2,1,J),J=1,81)/.3136E1,.6498E-2,.2289,.1859E-1,-.3328E-1,         
     &-.4889E-2,-.3054E-1,-.1773E-1,-.1728E-1,.6555E-1,.1775E-1,                
     &-.2488E-1,-.9498E-2,.1493E-1,.281E-2,.2406E-2,.5436E-2,-.2115,            
     &.7007E-2,-.5129E-1,-.7327E-2,.2402E-1,.4772E-2,-.7374E-2,                 
     &-.3835E-3,-.5013E-1,.2866E-2,.2216E-2,.2412E-3,.2094E-2,.122E-2           
     &,-.1703E-3,-.1082,-.4992E-2,-.4065E-2,.3615E-2,-.2738E-2,                 
     &-.7177E-3,.2173E-3,-.4373E-1,-.375E-2,.5507E-2,-.1567E-2,                 
     &-.1458E-2,-.7397E-3,.7903E-1,.4131E-2,.3714E-2,.1073E-2,                  
     &-.8991E-3,.2976E-3,.2623E-1,.2344E-2,.5608E-3,.4124E-3,.1509E-3,          
     &.5103E-1,.345E-2,.1283E-2,.7238E-3,-.3464E-4,.1663E-1,-.1644E-2,          
     &-.71E-3,.5281E-3,-.2729E-1,.3556E-2,-.3391E-2,-.1787E-3,.2154E-2,         
     &.6476E-2,-.8282E-3,-.2361E-1,.9557E-3,.3205E-3,-.2301E-1,                 
     &-.854E-3,-.1126E-1,-.2323E-2,-.8582E-2,.2683E-1/                          
      DATA (C(2,2,J),J=1,81)/.3144E1,.8571E-2,.2539,.6937E-1,-.1667E-1,         
     &.2249E-1,-.4162E-1,.1201E-1,.2435E-1,.5232E-1,.2521E-1,-.199E-1,          
     &-.7671E-2,.1264E-1,-.1551E-2,-.1928E-2,.3652E-2,-.2019,.5697E-2,          
     &-.3159E-1,-.1451E-1,.2868E-1,.1377E-1,-.4383E-2,.1172E-1,                 
     &-.5683E-1,.3593E-2,.3571E-2,.3282E-2,.1732E-2,-.4921E-3,-.1165E-2         
     &,-.1066,-.1892E-1,.357E-2,-.8631E-3,-.1876E-2,-.8414E-4,.2356E-2,         
     &-.4259E-1,-.322E-2,.4641E-2,.6223E-3,-.168E-2,-.1243E-3,.7393E-1,         
     &-.3143E-2,-.2362E-2,.1235E-2,-.1551E-2,.2099E-3,.2299E-1,.5301E-2         
     &,-.4306E-2,-.1303E-2,.7687E-5,.5305E-1,.6642E-2,-.1686E-2,                
     &.1048E-2,.5958E-3,.4341E-1,-.8819E-4,-.333E-3,-.2158E-3,-.4106E-1         
     &,.4191E-2,.2045E-2,-.1437E-3,-.1803E-1,-.8072E-3,-.424E-3,                
     &-.26E-1,-.2329E-2,.5949E-3,-.1371E-1,-.2188E-2,.1788E-1,                  
     &.6405E-3,.5977E-2,.1333E-1/                  
      DATA (C(3,1,J),J=1,81)/.3372E1,.1006E-1,.1436,.2023E-2,-.5166E-1,         
     &.9606E-2,-.5596E-1,.4914E-3,-.3124E-2,-.4713E-1,-.7371E-2,                
     &-.4823E-2,-.2213E-2,.6569E-2,-.1962E-3,.3309E-3,-.3908E-3,                
     &-.2836,.7829E-2,.1175E-1,.9919E-3,.6589E-2,.2045E-2,-.7346E-2             
     &,-.89E-3,-.347E-1,-.4977E-2,.147E-2,-.2823E-5,.6465E-3,                   
     &-.1448E-3,.1401E-2,-.8988E-1,-.3293E-4,-.1848E-2,.4439E-3,                
     &-.1263E-2,.317E-3,-.6227E-3,.1721E-1,-.199E-2,-.4627E-3,                  
     &.2897E-5,-.5454E-3,.3385E-3,.8432E-1,-.1951E-2,.1487E-2,                  
     &.1042E-2,-.4788E-3,-.1276E-3,.2373E-1,.2409E-2,.5263E-3,                  
     &.1301E-2,-.4177E-3,.3974E-1,.1418E-3,-.1048E-2,-.2982E-3,                 
     &-.3396E-4,.131E-1,.1413E-2,-.1373E-3,.2638E-3,-.4171E-1,                  
     &-.5932E-3,-.7523E-3,-.6883E-3,-.2355E-1,.5695E-3,-.2219E-4,               
     &-.2301E-1,-.9962E-4,-.6761E-3,.204E-2,-.5479E-3,.2591E-1,                 
     &-.2425E-2,.1583E-1,.9577E-2/                 
      DATA (C(3,2,J),J=1,81)/.3367E1,.1038E-1,.1407,.3622E-1,-.3144E-1,         
     &.112E-1,-.5674E-1,.3219E-1,.1288E-2,-.5799E-1,-.4609E-2,                  
     &.3252E-2,-.2859E-3,.1226E-1,-.4539E-2,.1310E-2,-.5603E-3,                 
     &-.311,-.1268E-2,.1539E-1,.3146E-2,.7787E-2,-.143E-2,-.482E-2              
     &,.2924E-2,-.9981E-1,-.7838E-2,-.1663E-3,.4769E-3,.4148E-2,                
     &-.1008E-2,-.979E-3,-.9049E-1,-.2994E-2,-.6748E-2,-.9889E-3,               
     &.1488E-2,-.1154E-2,-.8412E-4,-.1302E-1,-.4859E-2,-.7172E-3,               
     &-.9401E-3,.9101E-3,-.1735E-3,.7055E-1,.6398E-2,-.3103E-2,                 
     &-.938E-3,-.4E-3,-.1165E-2,.2713E-1,-.1654E-2,.2781E-2,                    
     &-.5215E-5,.2258E-3,.5022E-1,.95E-2,.4147E-3,.3499E-3,                     
     &-.6097E-3,.4118E-1,.6556E-2,.3793E-2,-.1226E-3,-.2517E-1,                 
     &.1491E-3,.1075E-2,.4531E-3,-.9012E-2,.3343E-2,.3431E-2,                   
     &-.2519E-1,.3793E-4,.5973E-3,-.1423E-1,-.132E-2,-.6048E-2,                 
     &-.5005E-2,-.115E-1,.2574E-1/                 
      DATA (C(4,1,J),J=1,81)/.3574E1,.0,.7537E-1,.0,-.8459E-1,                  
     &0.,-.294E-1,0.,.4547E-1,-.5321E-1,0.,.4328E-2,0.,.6022E-2,                
     &.0,-.9168E-3,.0,-.1768,.0,.294E-1,.0,.5902E-3,.0,-.9047E-2,               
     &.0,-.6555E-1,.0,-.1033E-2,.0,.1674E-2,.0,.2802E-3,-.6786E-1               
     &,.0,.4193E-2,.0,-.6448E-3,.0,.9277E-3,-.1634E-1,.0,-.2531E-2              
     &,.0,.193E-4,.0,.528E-1,.0,.2438E-2,.0,-.5292E-3,.0,.1555E-1               
     &,.0,-.3259E-2,.0,-.5998E-3,.3168E-1,.0,.2382E-2,.0,-.4078E-3              
     &,.2312E-1,.0,.1481E-3,.0,-.1885E-1,.0,.1144E-2,.0,-.9952E-2               
     &,.0,-.551E-3,-.202E-1,.0,-.7283E-4,-.1272E-1,.0,.2224E-2,                 
     &.0,-.251E-2,.2434E-1/                        
      DATA (C(4,2,J),J=1,81)/.3574E1,-.5639E-2,.7094E-1,                        
     &-.3347E-1,-.861E-1,-.2877E-1,-.3154E-1,-.2847E-2,.1235E-1,                
     &-.5966E-1,-.3236E-2,.3795E-3,-.8634E-3,.3377E-2,-.1071E-3,                
     &-.2151E-2,-.4057E-3,-.1783,.126E-1,.2835E-1,-.242E-2,                     
     &.3002E-2,-.4684E-2,-.6756E-2,-.7493E-3,-.6147E-1,-.5636E-2                
     &,-.1234E-2,-.1613E-2,-.6353E-4,-.2503E-3,-.1729E-3,-.7148E-1              
     &,.5326E-2,.4006E-2,.6484E-3,-.1046E-3,-.6034E-3,-.9435E-3,                
     &-.2385E-2,.6853E-2,.151E-2,.1319E-2,.9049E-4,-.1999E-3,                   
     &.3976E-1,.2802E-2,-.103E-2,.5599E-3,-.4791E-3,-.846E-4,                   
     &.2683E-1,.427E-2,.5911E-3,.2987E-3,-.208E-3,.1396E-1,                     
     &-.1922E-2,-.1063E-2,.3803E-3,.1343E-3,.1771E-1,-.1038E-2,                 
     &-.4645E-3,-.2481E-3,-.2251E-1,-.29E-2,-.3977E-3,-.516E-3,                 
     &-.8079E-2,-.1528E-2,.306E-3,-.1582E-1,-.8536E-3,.1565E-3,                 
     &-.1252E-1,.2319E-3,.4311E-2,.1024E-2,.1296E-5,.179E-1/                    
        IF(NS.LT.3) THEN
           IS=NS
        ELSE IF(NS.GT.3) THEN
           IS=2
           DIPL=-DIPL
        ELSE
           IS=1
        ENDIF
      COLAT=UMR*(90.-DIPL)                    
      AZ=humr*SLT    
      CALL SPHARM(A,8,8,COLAT,AZ)
        IF(IS.EQ.2) THEN
           KEND=3
        ELSE
           KEND=4
        ENDIF                  
      DO 2 K=1,KEND      
      STE=0.          
      DO 1 I=1,81     
1       STE=STE+A(I)*C(K,IS,I)                       
2     TE(K)=10.**STE
        IF(IS.EQ.2) THEN
           DIPL=-DIPL
           COLAT=UMR*(90.-DIPL)                    
           CALL SPHARM(A,8,8,COLAT,AZ)
           STE=0.          
           DO 11 I=1,81     
11            STE=STE+A(I)*C(4,2,I)                       
           TE(4)=10.**STE
        ENDIF

C---------- TEMPERATURE AT 400KM AT MIDNIGHT AND NOON
      DO 4 J=1,2      
        STE=0.          
        AZ=humr*(J-1)*12.                           
        CALL SPHARM(A,8,8,COLAT,AZ)                  
        DO 3 I=1,81     
3         STE=STE+A(I)*C(2,IS,I)                       
4       TE(J+4)=10.**STE                             
      RETURN          
      END             
C
      SUBROUTINE SPHARM(C,L,M,COLAT,AZ)            
C CALCULATES THE COEFFICIENTS OF THE SPHERICAL HARMONIC                         
C EXPANSION THAT WAS USED FOR THE BRACE-THEIS-MODELS.                           
      DIMENSION C(82)                              
      C(1)=1.         
      K=2             
      X=COS(COLAT)    
      C(K)=X          
      K=K+1           
      DO 10 I=2,L     
      C(K)=((2*I-1)*X*C(K-1)-(I-1)*C(K-2))/I       
10    K=K+1           
      Y=SIN(COLAT)    
      DO 20 MT=1,M    
      CAZ=COS(MT*AZ)  
      SAZ=SIN(MT*AZ)  
      C(K)=Y**MT      
      K=K+1           
      IF(MT.EQ.L) GOTO 16                          
      C(K)=C(K-1)*X*(2*MT+1)                       
      K=K+1           
      IF((MT+1).EQ.L) GOTO 16                      
      DO 15 I=2+MT,L  
      C(K)=((2*I-1)*X*C(K-1)-(I+MT-1)*C(K-2))/(I-MT)                            
15    K=K+1           
16    N=L-MT+1        
      DO 18 I=1,N     
      C(K)=C(K-N)*CAZ                              
      C(K-N)=C(K-N)*SAZ                            
18    K=K+1           
20    CONTINUE        
      RETURN          
      END             
C
C
      REAL FUNCTION ELTE(H)
c----------------------------------------------------------------
C ELECTRON TEMPERATURE PROFILE BASED ON THE TEMPERATURES AT 120                 
C HMAX,300,400,600,1400,3000 KM ALTITUDE. INBETWEEN CONSTANT                    
C GRADIENT IS ASSUMED. ARGMAX IS MAXIMUM ARGUMENT ALLOWED FOR
C EXP-FUNCTION.
c----------------------------------------------------------------
      COMMON /BLOTE/AH(7),ATE1,ST(6),D(5)
C
      SUM=ATE1+ST(1)*(H-AH(1))                     
      DO 1 I=1,5
        aa = eptr(h    ,d(i),ah(i+1))
        bb = eptr(ah(1),d(i),ah(i+1))
1     SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*D(I)                
      ELTE=SUM        
      RETURN          
      END             
C
C
      FUNCTION TEDE(H,DEN,COV)                     
C ELECTRON TEMEPERATURE MODEL AFTER BRACE,THEIS .  
C FOR NEG. COV THE MEAN COV-INDEX (3 SOLAR ROT.) IS EXPECTED.                   
C DEN IS THE ELECTRON DENSITY IN M-3.              
      Y=1051.+(17.01*H-2746.)*                     
     &EXP(-5.122E-4*H+(6.094E-12-3.353E-14*H)*DEN) 
      ACOV=ABS(COV)   
      YC=1.+(.117+2.02E-3*ACOV)/(1.+EXP(-(ACOV-102.5)/5.))                      
      IF(COV.LT.0.)   
     &YC=1.+(.123+1.69E-3*ACOV)/(1.+EXP(-(ACOV-115.)/10.))                      
      TEDE=Y*YC       
      RETURN          
      END             
C
C                     
C*************************************************************                  
C**************** ION TEMPERATURE ****************************
C*************************************************************                  
C
C
      REAL FUNCTION TI(H)
c----------------------------------------------------------------
C ION TEMPERATURE FOR HEIGHTS NOT GREATER 1000 KM AND NOT LESS HS               
C EXPLANATION SEE FUNCTION RPID.                   
c----------------------------------------------------------------
      REAL              MM
      COMMON  /BLOCK8/  HS,TNHS,XSM(4),MM(5),G(4),M

      SUM=MM(1)*(H-HS)+TNHS                        
      DO 100 I=1,M-1  
        aa = eptr(h ,g(i),xsm(i))
        bb = eptr(hs,g(i),xsm(i))
100     SUM=SUM+(MM(I+1)-MM(I))*(AA-BB)*G(I)                
      TI=SUM          
      RETURN          
      END             
C
C
      REAL FUNCTION TEDER(H)                       
C THIS FUNCTION ALONG WITH PROCEDURE REGFA1 ALLOWS TO FIND                      
C THE  HEIGHT ABOVE WHICH TN BEGINS TO BE DIFFERENT FROM TI                     
      COMMON    /BLOTN/XSM1,TEX,TLBD,SIG
      TNH = TN(H,TEX,TLBD,SIG)                        
      DTDX = DTNDH(H,TEX,TLBD,SIG)                        
      TEDER = DTDX * ( XSM1 - H ) + TNH                    
      RETURN          
      END             
C
C
      FUNCTION TN(H,TINF,TLBD,S)
C--------------------------------------------------------------------
C       Calculate Temperature for MSIS/CIRA-86 model
C--------------------------------------------------------------------
      ZG2 = ( H - 120. ) * 6476.77 / ( 6356.77 + H )
      TN = TINF - TLBD * EXP ( - S * ZG2 )
      RETURN
      END
C
C
      FUNCTION DTNDH(H,TINF,TLBD,S)
C---------------------------------------------------------------------
      ZG1 = 6356.77 + H
      ZG2 = 6476.77 / ZG1
      ZG3 = ( H - 120. ) * ZG2
      DTNDH = - TLBD * EXP ( - S * ZG3 ) * ( S / ZG1 * ( ZG3 - ZG2 ) )
      RETURN
      END
C
C                     
C*************************************************************                  
C************* ION RELATIVE PRECENTAGE DENSITY *****************                
C*************************************************************                  
C
C
      REAL FUNCTION RPID (H, H0, N0, M, ST, ID, XS)
c------------------------------------------------------------------
C D.BILITZA,1977,THIS ANALYTIC FUNCTION IS USED TO REPRESENT THE                
C RELATIVE PRECENTAGE DENSITY OF ATOMAR AND MOLECULAR OXYGEN IONS.              
C THE M+1 HEIGHT GRADIENTS ST(M+1) ARE CONNECTED WITH EPSTEIN-                  
C STEP-FUNCTIONS AT THE STEP HEIGHTS XS(M) WITH TRANSITION                      
C THICKNESSES ID(M). RPID(H0,H0,N0,....)=N0.       
C ARGMAX is the highest allowed argument for EXP in your system.
c------------------------------------------------------------------
      REAL              N0         
      DIMENSION         ID(4), ST(5), XS(4)                
      COMMON  /ARGEXP/  ARGMAX

      SUM=(H-H0)*ST(1)                             
      DO 100  I=1,M   
              XI=ID(I)
                aa = eptr(h ,xi,xs(i))
                bb = eptr(h0,xi,xs(i))
100           SUM=SUM+(ST(I+1)-ST(I))*(AA-BB)*XI 
      IF(ABS(SUM).LT.ARGMAX) then
        SM=EXP(SUM)
      else IF(SUM.Gt.0.0) then
        SM=EXP(ARGMAX)
      else
        SM=0.0
      endif
      RPID= n0 * SM        
      RETURN          
      END             
C
c
      SUBROUTINE RDHHE (H,HB,RDOH,RDO2H,RNO,PEHE,RDH,RDHE)                      
C BILITZA,FEB.82,H+ AND HE+ RELATIVE PERECENTAGE DENSITY BELOW                  
C 1000 KM. THE O+ AND O2+ REL. PER. DENSITIES SHOULD BE GIVEN                   
C (RDOH,RDO2H). HB IS THE ALTITUDE OF MAXIMAL O+ DENSITY. PEHE                  
C IS THE PRECENTAGE OF HE+ IONS COMPARED TO ALL LIGHT IONS.                     
C RNO IS THE RATIO OF NO+ TO O2+DENSITY AT H=HB.   
      RDHE=0.0        
      RDH=0.0         
      IF(H.LE.HB) GOTO 100                         
      REST=100.0-RDOH-RDO2H-RNO*RDO2H              
      RDH=REST*(1.-PEHE/100.)                      
      RDHE=REST*PEHE/100.                          
100   RETURN          
      END             
C
C
      REAL FUNCTION RDNO(H,HB,RDO2H,RDOH,RNO)      
C D.BILITZA, 1978. NO+ RELATIVE PERCENTAGE DENSITY ABOVE 100KM.                 
C FOR MORE INFORMATION SEE SUBROUTINE RDHHE.       
      IF (H.GT.HB) GOTO 200                        
      RDNO=100.0-RDO2H-RDOH                        
      RETURN          
200   RDNO=RNO*RDO2H  
      RETURN          
      END
C
C
      SUBROUTINE  KOEFP1(PG1O)                     
C THIEMANN,1979,COEFFICIENTS PG1O FOR CALCULATING  O+ PROFILES                  
C BELOW THE F2-MAXIMUM. CHOSEN TO APPROACH DANILOV-                             
C SEMENOV'S COMPILATION.                           
      DIMENSION PG1O(80)                           
      REAL FELD (80)  
      DATA FELD/-11.0,-11.0,4.0,-11.0,0.08018,     
     &0.13027,0.04216,0.25  ,-0.00686,0.00999,     
     &5.113,0.1 ,170.0,180.0,0.1175,0.15,-11.0,    
     &1.0 ,2.0,-11.0,0.069,0.161,0.254,0.18,0.0161,                             
     &0.0216,0.03014,0.1,152.0,167.0,0.04916,      
     &0.17,-11.0,2.0,2.0,-11.0,0.072,0.092,0.014,0.21,                          
     &0.01389,0.03863,0.05762,0.12,165.0,168.0,0.008,                           
     &0.258,-11.0,1.0,3.0,-11.0,0.091,0.088,       
     &0.008,0.34,0.0067,0.0195,0.04,0.1,158.0,172.0,                            
     &0.01,0.24,-11.0,2.0,3.0, -11.0,0.083,0.102,  
     &0.045,0.03,0.00127,0.01,0.05,0.09,167.0,185.0,                            
     &0.015,0.18/     
      K=0             
      DO 10 I=1,80    
      K=K+1           
10    PG1O(K)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE KOEFP2(PG2O)                      
C THIEMANN,1979,COEFFICIENTS FOR CALCULATION OF O+ PROFILES                     
C ABOVE THE F2-MAXIMUM (DUMBS,SPENNER:AEROS-COMPILATION)                        
      DIMENSION PG2O(32)                           
      REAL FELD(32)   
      DATA FELD/1.0,-11.0,-11.0,1.0,695.0,-.000781,                             
     &-.00264,2177.0,1.0,-11.0,-11.0,2.0,570.0,    
     &-.002,-.0052,1040.0,2.0,-11.0,-11.0,1.0,695.0,                            
     &-.000786,-.00165,3367.0,2.0,-11.0,-11.0,2.0, 
     &575.0,-.00126,-.00524,1380.0/                
      K=0             
      DO 10 I=1,32    
      K=K+1           
10    PG2O(K)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE  KOEFP3(PG3O)                     
C THIEMANN,1979,COEFFICIENTS FOR CALCULATING O2+ PROFILES.                      
C CHOSEN AS TO APPROACH DANILOV-SEMENOV'S COMPILATION.                          
      DIMENSION PG3O(80)                           
      REAL FELD(80)   
      DATA FELD/-11.0,1.0,2.0,-11.0,160.0,31.0,130.0,                           
     &-10.0,198.0,0.0,0.05922,-0.07983,            
     &-0.00397,0.00085,-0.00313,0.0,-11.0,2.0,2.0,-11.0,                        
     &140.0,30.0,130.0,-10.0,                      
     &190.0,0.0,0.05107,-0.07964,0.00097,-0.01118,-0.02614,                     
     &-0.09537,       
     &-11.0,1.0,3.0,-11.0,140.0,37.0,125.0,0.0,182.0,                           
     &0.0,0.0307,-0.04968,-0.00248,                
     &-0.02451,-0.00313,0.0,-11.0,2.0,3.0,-11.0,   
     &140.0,37.0,125.0,0.0,170.0,0.0,              
     &0.02806,-0.04716,0.00066,-0.02763,-0.02247,-0.01919,                      
     &-11.0,-11.0,4.0,-11.0,140.0,45.0,136.0,-9.0, 
     &181.0,-26.0,0.02994,-0.04879,                
     &-0.01396,0.00089,-0.09929,0.05589/           
      K=0             
      DO 10 I=1,80    
      K=K+1           
10    PG3O(K)=FELD(I)                              
      RETURN          
      END             
C
C
      SUBROUTINE SUFE (FIELD,RFE,M,FE)             
C SELECTS THE REQUIRED ION DENSITY PARAMETER SET.
C THE INPUT FIELD INCLUDES DIFFERENT SETS OF DIMENSION M EACH                
C CARACTERISED BY 4 HEADER NUMBERS. RFE(4) SHOULD CONTAIN THE                   
C CHOSEN HEADER NUMBERS.FE(M) IS THE CORRESPONDING SET.                         
      DIMENSION RFE(4),FE(12),FIELD(80),EFE(4)     
      K=0             
100   DO 101 I=1,4    
      K=K+1           
101   EFE(I)=FIELD(K)                              
      DO 111 I=1,M    
      K=K+1           
111   FE(I)=FIELD(K)  
      DO 120 I=1,4    
      IF((EFE(I).GT.-10.0).AND.(RFE(I).NE.EFE(I))) GOTO 100                     
120   CONTINUE        
      RETURN          
      END             
C
C                     
C*************************************************************                  
C************* PEAK VALUES ELECTRON DENSITY ******************                  
C*************************************************************                  
C
C
      real function FOUT(XMODIP,XLATI,XLONGI,UT,FF0)
C CALCULATES CRITICAL FREQUENCY FOF2/MHZ USING SUBROUTINE GAMMA1.      
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME 
C (DEC. HOURS), FF0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
      DIMENSION FF0(988)
      INTEGER QF(9)
      DATA QF/11,11,8,4,1,0,0,0,0/
      FOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,6,QF,9,76,13,988,FF0)
      RETURN
      END
C
C
      real function XMOUT(XMODIP,XLATI,XLONGI,UT,XM0)
C CALCULATES PROPAGATION FACTOR M3000 USING THE SUBROUTINE GAMMA1.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME 
C (DEC. HOURS), XM0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
      DIMENSION XM0(441)
      INTEGER QM(7)
      DATA QM/6,7,5,2,1,0,0/
      XMOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,4,QM,7,49,9,441,XM0)
      RETURN
      END
C
C
      REAL FUNCTION HMF2ED(XMAGBR,R,X,XM3)         
C CALCULATES THE PEAK HEIGHT HMF2/KM FOR THE MAGNETIC                           
C LATITUDE XMAGBR/DEG. AND THE SMOOTHED ZUERICH SUNSPOT                         
C NUMBER R USING CCIR-M3000 XM3 AND THE RATIO X=FOF2/FOE.                       
C [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]                       
C D.BILITZA,1980.     
      F1=(2.32E-3)*R+0.222                         
      F2=1.2-(1.16E-2)*EXP((2.39E-2)*R)            
      F3=0.096*(R-25.0)/150.0                      
      DELM=F1*(1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0))/(X-F2)+F3                
      HMF2ED=1490.0/(XM3+DELM)-176.0               
      RETURN          
      END             
C
C
      REAL FUNCTION FOF1ED(YLATI,R,CHI)
c--------------------------------------------------------------
C CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
C FOR   DIP-LATITUDE (YLATI/DEGREE)
c       SMOOTHED ZURICH SUNSPOT NUMBER (R)
c       SOLAR ZENITH ANGLE (CHI/DEGREE)
C REFERENCE: 
c       E.D.DUCHARME ET AL., RADIO SCIENCE 6, 369-378, 1971
C                                      AND 8, 837-839, 1973
c       HOWEVER WITH MAGNETIC DIP LATITUDE INSTEAD OF GEOMAGNETIC
c       DIPOLE LATITUDE, EYFRIG, 1979                    
C--------------------------------------------- D. BILITZA, 1988.   
        COMMON/CONST/UMR
        FOF1 = 0.0
        DLA =  YLATI
                CHI0 = 49.84733 + 0.349504 * DLA
                CHI100 = 38.96113 + 0.509932 * DLA
                CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
                IF(CHI.GT.CHIM) GOTO 1 
        F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA ) 
        F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
        FS = F0 + ( F100 - F0 ) * R / 100.0
        XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
        FOF1 = FS * COS( CHI * UMR ) ** XMUE
1       FOF1ED = FOF1     
        RETURN
        END             
C
C
      REAL FUNCTION FOEEDI(COV,XHI,XHIM,XLATI)
C-------------------------------------------------------
C CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.      
C INPUT: MEAN 10.7CM SOLAR RADIO FLUX (COV), GEOGRAPHIC
C LATITUDE (XLATI/DEG), SOLAR ZENITH ANGLE (XHI/DEG AND 
C XHIM/DEG AT NOON).
C REFERENCE: 
C       KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
C       TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
C               to improve the nighttime varition)
C D.BILITZA--------------------------------- AUGUST 1986.    
      COMMON/CONST/UMR
C variation with solar activity (factor A) ...............
      A=1.0+0.0094*(COV-66.0)                      
C variation with noon solar zenith angle (B) and with latitude (C)
      SL=COS(XLATI*UMR)
        IF(XLATI.LT.32.0) THEN
                SM=-1.93+1.92*SL                             
                C=23.0+116.0*SL                              
        ELSE
                SM=0.11-0.49*SL                              
                C=92.0+35.0*SL  
        ENDIF
        if(XHIM.ge.90.) XHIM=89.999
        B = COS(XHIM*UMR) ** SM
C variation with solar zenith angle (D) ..........................        
        IF(XLATI.GT.12.0) THEN
                SP=1.2
        ELSE
                SP=1.31         
        ENDIF
C adjusted solar zenith angle during nighttime (XHIC) .............
      XHIC=XHI-3.*ALOG(1.+EXP((XHI-89.98)/3.))   
      D=COS(XHIC*UMR)**SP       
C determine foE**4 ................................................
      R4FOE=A*B*C*D     
C minimum allowable foE (sqrt[SMIN])...............................
      SMIN=0.121+0.0015*(COV-60.)
      SMIN=SMIN*SMIN
      IF(R4FOE.LT.SMIN) R4FOE=SMIN                     
      FOEEDI=R4FOE**0.25                           
      RETURN          
      END   
C
C
      REAL FUNCTION XMDED(XHI,R,YW)                
C D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.                   
C XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER              
C AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.     
C [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,BOULDER,1981]
C corrected 4/25/97 - D. Bilitza
c
      COMMON/CONST/UMR
c
        if(xhi.ge.90) goto 100
      Y = 6.05E8 + 0.088E8 * R
        yy = cos ( xhi * umr )
        ymd = y * exp( -0.1 / ( yy**2.7 ) )
        if (ymd.lt.yw) ymd = yw
        xmded=ymd
      RETURN          
100   XMDED=YW        
      RETURN          
      END
C
C
      REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,IHARM,NQ,
     &                          K1,M,MM,M3,SFE)      
C CALCULATES GAMMA1=FOF2 OR M3000 USING CCIR NUMERICAL MAP                      
C COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
C GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)  
C AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS).
C NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN 
C LATITUDE FOR EACH LONGITUDE HARMONIC.                  
C M=1+NQ1+2(NQ2+1)+2(NQ3+1)+... .                  
C SHEIKH,4.3.77.      
      REAL*8 C(12),S(12),COEF(100),SUM             
      DIMENSION NQ(K1),XSINX(13),SFE(M3)           
      COMMON/CONST/UMR
      HOU=(15.0*HOUR-180.0)*UMR                    
      S(1)=SIN(HOU)   
      C(1)=COS(HOU)   
      DO 250 I=2,IHARM                             
      C(I)=C(1)*C(I-1)-S(1)*S(I-1)                 
      S(I)=C(1)*S(I-1)+S(1)*C(I-1)                 
250   CONTINUE        
      DO 300 I=1,M    
      MI=(I-1)*MM     
      COEF(I)=SFE(MI+1)                            
      DO 300 J=1,IHARM                             
      COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)                       
300   CONTINUE        
      SUM=COEF(1)     
      SS=SIN(SMODIP*UMR)                           
      S3=SS           
      XSINX(1)=1.0    
      INDEX=NQ(1)     
      DO 350 J=1,INDEX                             
      SUM=SUM+COEF(1+J)*SS                         
      XSINX(J+1)=SS   
      SS=SS*S3        
350   CONTINUE        
      XSINX(NQ(1)+2)=SS                            
      NP=NQ(1)+1      
      SS=COS(SLAT*UMR)                             
      S3=SS           
      DO 400 J=2,K1   
      S0=SLONG*(J-1.)*UMR                          
      S1=COS(S0)      
      S2=SIN(S0)      
      INDEX=NQ(J)+1   
      DO 450 L=1,INDEX                             
      NP=NP+1         
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S1              
      NP=NP+1         
      SUM=SUM+COEF(NP)*XSINX(L)*SS*S2              
450   CONTINUE        
      SS=SS*S3        
400   CONTINUE        
      GAMMA1=SUM      
      RETURN          
      END             
C
C                     
C************************************************************                   
C*************** EARTH MAGNETIC FIELD ***********************                   
C**************************************************************                 
C
C
      SUBROUTINE GGM(ART,LONG,LATI,MLONG,MLAT)            
C CALCULATES GEOMAGNETIC LONGITUDE (MLONG) AND LATITUDE (MLAT) 
C FROM GEOGRAFIC LONGITUDE (LONG) AND LATITUDE (LATI) FOR ART=0
C AND REVERSE FOR ART=1. ALL ANGLES IN DEGREE.
C LATITUDE:-90 TO 90. LONGITUDE:0 TO 360 EAST.         
      INTEGER ART     
      REAL MLONG,MLAT,LONG,LATI
      COMMON/CONST/FAKTOR
      ZPI=FAKTOR*360.                              
      CBG=11.4*FAKTOR                              
      CI=COS(CBG)     
      SI=SIN(CBG)
      IF(ART.EQ.0) GOTO 10                         
      CBM=COS(MLAT*FAKTOR)                           
      SBM=SIN(MLAT*FAKTOR)                           
      CLM=COS(MLONG*FAKTOR)                          
      SLM=SIN(MLONG*FAKTOR)
      SBG=SBM*CI-CBM*CLM*SI
        IF(ABS(SBG).GT.1.) SBG=SIGN(1.,SBG)
      LATI=ASIN(SBG)
      CBG=COS(LATI)     
      SLG=(CBM*SLM)/CBG  
      CLG=(SBM*SI+CBM*CLM*CI)/CBG
        IF(ABS(CLG).GT.1.) CLG=SIGN(1.,CLG)                  
      LONG=ACOS(CLG)  
      IF(SLG.LT.0.0) LONG=ZPI-LONG
      LATI=LATI/FAKTOR    
      LONG=LONG/FAKTOR  
      LONG=LONG-69.8    
      IF(LONG.LT.0.0) LONG=LONG+360.0                 
      RETURN          
10    YLG=LONG+69.8    
      CBG=COS(LATI*FAKTOR)                           
      SBG=SIN(LATI*FAKTOR)                           
      CLG=COS(YLG*FAKTOR)                          
      SLG=SIN(YLG*FAKTOR)                          
      SBM=SBG*CI+CBG*CLG*SI                        
        IF(ABS(SBM).GT.1.) SBM=SIGN(1.,SBM)
      MLAT=ASIN(SBM)   
      CBM=COS(MLAT)     
      SLM=(CBG*SLG)/CBM                            
      CLM=(-SBG*SI+CBG*CLG*CI)/CBM
        IF(ABS(CLM).GT.1.) CLM=SIGN(1.,CLM) 
      MLONG=ACOS(CLM)
      IF(SLM.LT..0) MLONG=ZPI-MLONG
      MLAT=MLAT/FAKTOR    
      MLONG=MLONG/FAKTOR  
      RETURN          
      END             
C
C
      SUBROUTINE FIELDG(DLAT,DLONG,ALT,X,Y,Z,F,DIP,DEC,SMODIP)                  
C THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD                    
C LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.                  
C INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),                 
C        ALT=ALTITUDE/KM.                          
C OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT                  
C        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).          
C        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.             
C SHEIK,1977.         
      DIMENSION H(144),XI(3),G(144),FEL1(72),FEL2(72)
      COMMON/CONST/UMR                           
      DATA FEL1/0.0, 0.1506723,0.0101742, -0.0286519, 0.0092606,                
     & -0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,                
     & 0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,                
     & 0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,                
     & 0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,                 
     & 0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,               
     &0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,                     
     &-0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,                    
     & -0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,               
     & -0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,                
     &  0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,                 
     & 0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,                
     & -0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,                 
     & -0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,               
     &0.0083033,  0.0077187/                       
      DATA FEL2/0.0586055,0.0102236,-0.0396107,    
     & -0.0167860, -0.2019911, -0.5810815,0.0379916,  3.7508268,                
     & 1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,                
     & -0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,                  
     &  0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,                 
     & -0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,               
     & 0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,                 
     &  0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,               
     &0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,                  
     &-6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,               
     &37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,               
     & -0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,                   
     &  0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,                
     & 0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,                 
     & 215.5209961, -27.1419220,22.3405762,1108.6394043/                        
      K=0             
      DO 10 I=1,72    
      K=K+1           
      G(K)=FEL1(I)    
10    G(72+K)=FEL2(I)                              
      RLAT=DLAT*UMR   
      CT=SIN(RLAT)    
      ST=COS(RLAT)    
      NMAX=11         
      D=SQRT(40680925.0-272336.0*CT*CT)            
      RLONG=DLONG*UMR                              
      CP=COS(RLONG)   
      SP=SIN(RLONG)   
      ZZZ=(ALT+40408589.0/D)*CT/6371.2             
      RHO=(ALT+40680925.0/D)*ST/6371.2             
      XXX=RHO*CP      
      YYY=RHO*SP      
      RQ=1.0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ)             
      XI(1)=XXX*RQ    
      XI(2)=YYY*RQ    
      XI(3)=ZZZ*RQ    
      IHMAX=NMAX*NMAX+1                            
      LAST=IHMAX+NMAX+NMAX                         
      IMAX=NMAX+NMAX-1                             
      DO 100 I=IHMAX,LAST                          
100   H(I)=G(I)       
      DO 200 K=1,3,2  
      I=IMAX          
      IH=IHMAX        
300   IL=IH-I         
      F1=2./(I-K+2.)  
      X1=XI(1)*F1     
      Y1=XI(2)*F1     
      Z1=XI(3)*(F1+F1)                             
      I=I-2           
      IF((I-1).LT.0) GOTO 400                      
      IF((I-1).EQ.0) GOTO 500                      
      DO 600 M=3,I,2  
      H(IL+M+1)=G(IL+M+1)+Z1*H(IH+M+1)+X1*(H(IH+M+3)-H(IH+M-1))-                
     &Y1*(H(IH+M+2)+H(IH+M-2))                     
      H(IL+M)=G(IL+M)+Z1*H(IH+M)+X1*(H(IH+M+2)-H(IH+M-2))+                      
     &Y1*(H(IH+M+3)+H(IH+M-1))                     
600   CONTINUE        
500   H(IL+2)=G(IL+2)+Z1*H(IH+2)+X1*H(IH+4)-Y1*(H(IH+3)+H(IH))                  
      H(IL+1)=G(IL+1)+Z1*H(IH+1)+Y1*H(IH+4)+X1*(H(IH+3)-H(IH))                  
400   H(IL)=G(IL)+Z1*H(IH)+2.0*(X1*H(IH+1)+Y1*H(IH+2))                          
700   IH=IL           
      IF(I.GE.K) GOTO 300                          
200   CONTINUE        
      S=0.5*H(1)+2.0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))                         
      XT=(RQ+RQ)*SQRT(RQ)                          
      X=XT*(H(3)-S*XXX)                            
      Y=XT*(H(4)-S*YYY)                            
      Z=XT*(H(2)-S*ZZZ)                            
      F=SQRT(X*X+Y*Y+Z*Z)                          
      BRH0=Y*SP+X*CP  
      Y=Y*CP-X*SP     
      X=Z*ST-BRH0*CT  
      Z=-Z*CT-BRH0*ST 
        zdivf=z/f
        IF(ABS(zdivf).GT.1.) zdivf=SIGN(1.,zdivf)
      DIP=ASIN(zdivf)
        ydivs=y/sqrt(x*x+y*y)  
        IF(ABS(ydivs).GT.1.) ydivs=SIGN(1.,ydivs)
      DEC=ASIN(ydivs)
        dipdiv=DIP/SQRT(DIP*DIP+ST)
        IF(ABS(dipdiv).GT.1.) dipdiv=SIGN(1.,dipdiv)
      SMODIP=ASIN(dipdiv)
      DIP=DIP/UMR     
      DEC=DEC/UMR     
      SMODIP=SMODIP/UMR                            
      RETURN          
      END             
C
C
C************************************************************                   
C*********** INTERPOLATION AND REST ***************************                 
C**************************************************************                 
C
C
      SUBROUTINE REGFA1(X11,X22,FX11,FX22,EPS,FW,F,SCHALT,X) 
C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE                
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL                      
C HAS BECOME LESS THAN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)                  
C THEN SCHALT=.TRUE.  
      LOGICAL L1,LINKS,K,SCHALT                    
      SCHALT=.FALSE.
      EP=EPS  
      X1=X11          
      X2=X22          
      F1=FX11-FW     
      F2=FX22-FW     
      K=.FALSE.       
      NG=2       
      LFD=0     
      IF(F1*F2.LE.0.0) GOTO 200
        X=0.0           
        SCHALT=.TRUE.   
        RETURN
200   X=(X1*F2-X2*F1)/(F2-F1)                      
      GOTO 400        
300     L1=LINKS        
        DX=(X2-X1)/NG
        IF(.NOT.LINKS) DX=DX*(NG-1)
        X=X1+DX
400   FX=F(X)-FW
      LFD=LFD+1
      IF(LFD.GT.20) THEN
        EP=EP*10.
        LFD=0
      ENDIF 
      LINKS=(F1*FX.GT.0.0)
      K=.NOT.K        
      IF(LINKS) THEN
        X1=X            
        F1=FX           
      ELSE
        X2=X 
        F2=FX 
      ENDIF   
      IF(ABS(X2-X1).LE.EP) GOTO 800               
      IF(K) GOTO 300  
      IF((LINKS.AND.(.NOT.L1)).OR.(.NOT.LINKS.AND.L1)) NG=2*NG                  
      GOTO 200        
800   RETURN          
      END             
C
C
      SUBROUTINE TAL(SHABR,SDELTA,SHBR,SDTDH0,AUS6,SPT)                         
C CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
C Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5               
C TO FIT THE VALLEY IN Y, REPRESENTED BY:                
C Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),                    
C THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE                       
C DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).                        
C IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY                     
C REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..     
C FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION                        
C Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).           
      DIMENSION SPT(4)                             
      LOGICAL AUS6    
      Z1=-SDELTA/(100.0*SHABR*SHABR)               
      IF(SDELTA.GT.0.) GOTO 500                    
      SDELTA=-SDELTA  
      Z1=ALOG(1.-SDELTA/100.)/(SHABR*SHABR)        
500   Z3=SDTDH0/(2.*SHBR)                          
      Z4=SHABR-SHBR   
      SPT(4)=2.0*(Z1*(SHBR-2.0*SHABR)*SHBR+Z3*Z4*SHABR)/                        
     &  (SHABR*SHBR*Z4*Z4*Z4)                        
      SPT(3)=Z1*(2.0*SHBR-3.0*SHABR)/(SHABR*Z4*Z4)-
     &  (2.*SHABR+SHBR)*SPT(4)          
      SPT(2)=-2.0*Z1/SHABR-2.0*SHABR*SPT(3)-3.0*SHABR*SHABR*SPT(4)              
      SPT(1)=Z1-SHABR*(SPT(2)+SHABR*(SPT(3)+SHABR*SPT(4)))                      
      AUS6=.FALSE.    
      B=4.*SPT(3)/(5.*SPT(4))+SHABR                
      C=-2.*SPT(1)/(5*SPT(4)*SHABR)                
      Z2=B*B/4.-C     
      IF(Z2.LT.0.0) GOTO 300                       
      Z3=SQRT(Z2)     
      Z1=B/2.         
      Z2=-Z1+Z3       
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
      IF (ABS(Z3).GT.1.E-15) GOTO 400              
      Z2=C/Z2         
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
      RETURN          
400   Z2=-Z1-Z3       
      IF(Z2.GT.0.0.AND.Z2.LT.SHBR) AUS6=.TRUE.     
300   RETURN          
      END             
C
C
C******************************************************************
C********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
C******************************************************************
C
C
        subroutine soco (ld,t,flat,Elon,
     &          DECLIN, ZENITH, SUNRSE, SUNSET)
c--------------------------------------------------------------------
c       s/r to calculate the solar declination, zenith angle, and
c       sunrise & sunset times  - based on Newbern Smith's algorithm
c       [leo mcnamara, 1-sep-86, last modified 16-jun-87]
c       {dieter bilitza, 30-oct-89, modified for IRI application}
c
c in:   ld      local day of year
c       t       local hour (decimal)
c       flat    northern latitude in degrees
c       elon    east longitude in degrees
c
c out:  declin      declination of the sun in degrees
c       zenith      zenith angle of the sun in degrees
c       sunrse      local time of sunrise in hours 
c       sunset      local time of sunset in hours 
c-------------------------------------------------------------------
c
        common/const/   dtr     /const1/humr,dumr
c amplitudes of Fourier coefficients  --  1955 epoch.................
        data    p1,p2,p3,p4,p6 /
     &  0.017203534,0.034407068,0.051610602,0.068814136,0.103221204 /
c
c s/r is formulated in terms of WEST longitude.......................
        wlon = 360. - Elon
c
c time of equinox for 1980...........................................
        td = ld + (t + Wlon/15.) / 24.
        te = td + 0.9369
c
c declination of the sun..............................................
        dcl = 23.256 * sin(p1*(te-82.242)) + 0.381 * sin(p2*(te-44.855))
     &      + 0.167 * sin(p3*(te-23.355)) - 0.013 * sin(p4*(te+11.97))
     &      + 0.011 * sin(p6*(te-10.41)) + 0.339137
        DECLIN = dcl
        dc = dcl * dtr
c
c the equation of time................................................
        tf = te - 0.5
        eqt = -7.38*sin(p1*(tf-4.)) - 9.87*sin(p2*(tf+9.))
     &      + 0.27*sin(p3*(tf-53.)) - 0.2*cos(p4*(tf-17.))
        et = eqt * dtr / 4.
c
        fa = flat * dtr
        phi = humr * ( t - 12.) + et
c
        a = sin(fa) * sin(dc)
        b = cos(fa) * cos(dc)
        cosx = a + b * cos(phi)
        if(abs(cosx).gt.1.) cosx=sign(1.,cosx)
        zenith = acos(cosx) / dtr
c
c calculate sunrise and sunset times --  at the ground...........
c see Explanatory Supplement to the Ephemeris (1961) pg 401......
c sunrise at height h metres is at...............................
c       chi(h) = 90.83 + 0.0347 * sqrt(h)........................
c this includes corrections for horizontal refraction and........
c semi-diameter of the solar disk................................
        ch = cos(90.83 * dtr)
        cosphi = (ch -a ) / b
c if abs(secphi) > 1., sun does not rise/set.....................
c allow for sun never setting - high latitude summer.............
        secphi = 999999.
        if(cosphi.ne.0.) secphi = 1./cosphi
        sunset = 99.
        sunrse = 99.
        if(secphi.gt.-1.0.and.secphi.le.0.) return
c allow for sun never rising - high latitude winter..............
        sunset = -99.
        sunrse = -99.
        if(secphi.gt.0.0.and.secphi.lt.1.) return
c
        if(cosphi.gt.1.) cosphi=sign(1.,cosphi)
        phi = acos(cosphi)
        et = et / humr
        phi = phi / humr
        sunrse = 12. - phi - et
        sunset = 12. + phi - et
        if(sunrse.lt.0.) sunrse = sunrse + 24.
        if(sunset.ge.24.) sunset = sunset - 24.
c
        return
        end
c
C
      FUNCTION HPOL(HOUR,TW,XNW,SA,SU,DSA,DSU)            
C-------------------------------------------------------
C PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN  
C STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE 
C STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
C TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO 
C BE INTERPOLATED. SA AND SU ARE TIME OF SUNRIES AND 
C SUNSET IN DECIMAL HOURS.
C BILITZA----------------------------------------- 1979.
        IF(ABS(SU).GT.25.) THEN
                IF(SU.GT.0.0) THEN
                        HPOL=TW
                ELSE
                        HPOL=XNW
                ENDIF
                RETURN
        ENDIF
      HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+
     &  (XNW-TW)*EPST(HOUR,DSU,SU) 
      RETURN          
      END       
C      
C
        SUBROUTINE MODA(IN,IYEAR,MONTH,IDAY,IDOY,NRDAYMO)
C-------------------------------------------------------------------
C CALCULATES DAY OF YEAR (IDOY, ddd) FROM YEAR (IYEAR, yy or yyyy), 
C MONTH (MONTH, mm) AND DAY OF MONTH (IDAY, dd) IF IN=0, OR MONTH 
C AND DAY FROM YEAR AND DAY OF YEAR IF IN=1. NRDAYMO is an output 
C parameter providing the number of days in the specific month.
C-------------------------------------------------------------------
        DIMENSION       MM(12)
        DATA            MM/31,28,31,30,31,30,31,31,30,31,30,31/

        IMO=0
        MOBE=0
c
c  leap year rule: years evenly divisible by 4 are leap years, except
c  years also evenly divisible by 100 are not leap years, except years also
c  evenly divisible by 400 are leap years. The year 2000 therefore is a
c  leap year. The 100 and 400 year exception rule
c       if((iyear/4*4.eq.iyear).and.(iyear/100*100.ne.iyear)) mm(2)=29
c  will become important again in the year 2100 which is not a leap year.
c
        mm(2)=28
        if(iyear/4*4.eq.iyear) mm(2)=29

        IF(IN.GT.0) GOTO 5
                mosum=0
                if(month.gt.1) then
                        do 1234 i=1,month-1 
1234                            mosum=mosum+mm(i)
                        endif
                idoy=mosum+iday
                nrdaymo=mm(month)
                RETURN

5       IMO=IMO+1
                IF(IMO.GT.12) GOTO 55
                MOOLD=MOBE
                nrdaymo=mm(imo)
                MOBE=MOBE+nrdaymo
                IF(MOBE.LT.IDOY) GOTO 5
55              MONTH=IMO
                IDAY=IDOY-MOOLD
        RETURN
        END             
c
c
        subroutine ut_lt(mode,ut,slt,glong,iyyy,ddd)
c -----------------------------------------------------------------
c Converts Universal Time UT (decimal hours) into Solar Local Time
c SLT (decimal hours) for given date (iyyy is year, e.g. 1995; ddd
c is day of year, e.g. 1 for Jan 1) and geodatic longitude in degrees.
C For mode=0 UT->LT and for mode=1 LT->UT
c Please NOTE that iyyy and ddd are input as well as output parameters
c since the determined LT may be for a day before or after the UT day.
c ------------------------------------------------- bilitza nov 95
        integer         ddd,dddend

        xlong=glong
        if(glong.gt.180) xlong=glong-360
        if(mode.ne.0) goto 1
c
c UT ---> LT
c
        SLT=UT+xlong/15.
        if((SLT.ge.0.).and.(SLT.le.24.)) goto 2
        if(SLT.gt.24.) goto 3
                SLT=SLT+24.
                ddd=ddd-1
                if(ddd.lt.1.) then
                        iyyy=iyyy-1
                        ddd=365
c
c leap year if evenly divisible by 4 and not by 100, except if evenly
c divisible by 400. Thus 2000 will be a leap year.
c
                        if(iyyy/4*4.eq.iyyy) ddd=366
                        endif
                goto 2
3               SLT=SLT-24.
                ddd=ddd+1
                dddend=365
                if(iyyy/4*4.eq.iyyy) dddend=366
                if(ddd.gt.dddend) then
                        iyyy=iyyy+1
                        ddd=1
                        endif
                goto 2
c
c LT ---> UT
c
1       UT=SLT-xlong/15.
        if((UT.ge.0.).and.(UT.le.24.)) goto 2
        if(UT.gt.24.) goto 5
                UT=UT+24.
                ddd=ddd-1
                if(ddd.lt.1.) then
                        iyyy=iyyy-1
                        ddd=365
                        if(iyyy/4*4.eq.iyyy) ddd=366
                        endif
                goto 2
5               UT=UT-24.
                ddd=ddd+1
                dddend=365
                if(iyyy/4*4.eq.iyyy) dddend=366
                if(ddd.gt.dddend) then
                        iyyy=iyyy+1
                        ddd=1
                        endif
2       return
        end
c
C
        REAL FUNCTION B0_TAB ( HOUR, SAX, SUX, NSEASN, R, ZMODIP)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and low and middle modified dip latitudes 18 and 45
C degress (ILATI=1,2). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C corrected to include a smooth transition at the modip equator
C and no discontinuity at the equatorial change in season.
C JAN 1993 ---------------------------------------- Dieter Bilitza
C
      REAL      NITVAL
      DIMENSION B0F(2,4,2,2),bfr(2,2,2),bfd(2,2),zx(5),g(6),dd(5)
      DATA      B0F/114.,64.0,134.,77.0,128.,66.0,75.,73.0,
     &              113.,115.,150.,116.,138.,123.,94.,132.,
     &              72.0,84.0,83.0,89.0,75.0,85.0,57.,76.0,
     &              102.,100.,120.,110.,107.,103.,76.,86.0/
        data    zx/45.,70.,90.,110.,135./,dd/2.5,2.,2.,2.,2.5/

C jseasn is southern hemisphere season
        jseasn=nseasn+2
        if(jseasn.gt.4) jseasn=jseasn-4

        zz = zmodip + 90.
        zz0 = 0.

C Interpolation in Rz12: linear from 10 to 100
        DO 7035 ISL=1,2
          DO 7034 ISD=1,2
            bfr(isd,1,isl) = b0f(isd,nseasn,1,isl) +
     &      (b0f(isd,nseasn,2,isl) - b0f(isd,nseasn,1,isl))/90.*(R-10.)
            bfr(isd,2,isl) = b0f(isd,jseasn,1,isl) +
     &      (b0f(isd,jseasn,2,isl) - b0f(isd,jseasn,1,isl))/90.*(R-10.)
7034      continue

C Interpolation day/night with transitions at SAX (sunrise) and SUX (sunset)
          do 7033 iss=1,2
                DAYVAL = BFR(1,ISS,ISL)
                NITVAL = BFR(2,ISS,ISL)
                BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
7033      continue
7035    continue

C Interpolation with epstein-transitions in modified dip latitude.
C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
C
C g(1:5) are the latitudinal slopes; g(1) is for the region from -90
C to -45 degrees, g(2) for -45/-20, g(3) for -20/0, g(4) for 0/20,
C g(5) for 20/45, and g(6) for 45/90. B0=bfd(2,2) at modip = -90, 
C bfd(2,2) at modip = -45, bfd(2,1) at modip = -20, bfd(2,1)+delta at 
C modip = -10 and 0, bfd(1,1) at modip = 20, bfd(1,2) at modip = 45 and 90.

        g(1) = 0.
        g(2) = ( bfd(2,1) - bfd(2,2) ) / 25.
        g(5) = ( bfd(1,2) - bfd(1,1) ) / 25.
        g(6) = 0.
        if(bfd(2,1).gt.bfd(1,1)) then
                g(3) = g(2) / 4.
                yb4 = bfd(2,1) + 20. * g(3) 
                g(4) = ( bfd(1,1) - yb4 ) / 20.
        else
                g(4) = g(5) / 4.
                yb5 = bfd(1,1) - 20. * g(4)
                g(3) = ( yb5 - bfd(2,1) ) / 20.
        endif
        bb0 = bfd(2,2)
      SUM = bb0                     
      DO 1 I=1,5
        aa = eptr(zz ,dd(i),zx(i))
        bb = eptr(zz0,dd(i),zx(i))
        DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)                
        SUM = SUM + DSUM
1       continue
      B0_TAB = SUM        

        RETURN
        END
c
c
        REAL FUNCTION B0_NEW ( HOUR, SAX, SUX, NSEASN, R, ZMODIP)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and modified dip latitudes of 0, 18 and 45
C degress (ILATI=1,2,3). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C corrected to include a smooth transition at the modip equator
C and no discontinuity at the equatorial change in season.
C JAN 1993 ---------------------------------------- Dieter Bilitza
C
      REAL      NITVAL
      DIMENSION B0F(2,4,2,3),bfr(2,2,3),bfd(2,3),zx(5),g(6),dd(5)
      DATA      B0F/185.,82.,185.,82.,185.,82.,185.,82.,
     &              255.,70.,255.,70.,255.,70.,255.,70.,
     &              114.,64.0,134.,77.0,128.,66.0,75.,73.0,
     &              113.,115.,150.,116.,138.,123.,94.,132.,
     &              72.0,84.0,83.0,89.0,75.0,85.0,57.,76.0,
     &              102.,100.,120.,110.,107.,103.,76.,86.0/
        data    zx/45.,70.,90.,110.,135./,dd/2.5,2.,2.,2.,2.5/

        num_lat=3

C jseasn is southern hemisphere season
        jseasn=nseasn+2
        if(jseasn.gt.4) jseasn=jseasn-4

        zz = zmodip + 90.
        zz0 = 0.

C Interpolation in Rz12: linear from 10 to 100
        DO 7035 ISL=1,num_lat
          DO 7034 ISD=1,2
            bfr(isd,1,isl) = b0f(isd,nseasn,1,isl) +
     &      (b0f(isd,nseasn,2,isl) - b0f(isd,nseasn,1,isl))/90.*(R-10.)
            bfr(isd,2,isl) = b0f(isd,jseasn,1,isl) +
     &      (b0f(isd,jseasn,2,isl) - b0f(isd,jseasn,1,isl))/90.*(R-10.)
7034      continue

C Interpolation day/night with transitions at SAX (sunrise) and SUX (sunset)
          do 7033 iss=1,2
                DAYVAL = BFR(1,ISS,ISL)
                NITVAL = BFR(2,ISS,ISL)
                BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
7033      continue
7035    continue

C Interpolation with epstein-transitions in modified dip latitude.
C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
C
C g(1:5) are the latitudinal slopes; g(1) is for the region from -90
C to -45 degrees, g(2) for -45/-20, g(3) for -20/0, g(4) for 0/20,
C g(5) for 20/45, and g(6) for 45/90. B0=bfd(2,2) at modip = -90, 
C bfd(2,2) at modip = -45, bfd(2,1) at modip = -20, bfd(2,1)+delta at 
C modip = -10 and 0, bfd(1,1) at modip = 20, bfd(1,2) at modip = 45 and 90.

        g(1) = 0.
        g(2) = ( bfd(2,2) - bfd(2,3) ) / 25.
        g(3) = ( bfd(2,1) - bfd(2,2) ) / 20.
        g(4) = ( bfd(1,2) - bfd(1,1) ) / 20.
        g(5) = ( bfd(1,3) - bfd(1,2) ) / 25.
        g(6) = 0.
c       if(bfd(2,1).gt.bfd(1,1)) then
c               g(3) = g(2) / 4.
c               yb4 = bfd(2,1) + 20. * g(3) 
c               g(4) = ( bfd(1,1) - yb4 ) / 20.
c       else
c               g(4) = g(5) / 4.
c               yb5 = bfd(1,1) - 20. * g(4)
c               g(3) = ( yb5 - bfd(2,1) ) / 20.
c       endif
        bb0 = bfd(2,3)
      SUM = bb0                     
      DO 1 I=1,5
        aa = eptr(zz ,dd(i),zx(i))
        bb = eptr(zz0,dd(i),zx(i))
        DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)                
        SUM = SUM + DSUM
1       continue
      B0_NEW = SUM        

        RETURN
        END
c
C
        REAL FUNCTION B0POL ( HOUR, SAX, SUX, ISEASON, R, DELA)
C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON=1 spring), low and high
C solar activity (IR=1,2), and low and middle modified dip
C latitudes (ILATI=1,2). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
      REAL              NITVAL
      DIMENSION         B0F(2,4,2,2),SIPH(2),SIPL(2)
      DATA      B0F/114.,64.0,134.,77.0,128.,66.0,75.,73.0,
     &              113.,115.,150.,116.,138.,123.,94.,132.,
     &              72.0,84.0,83.0,89.0,75.0,85.0,57.,76.0,
     &              102.,100.,120.,110.,107.,103.,76.,86.0/

        DO 7033 ISR=1,2
                DO 7034 ISL=1,2
                        DAYVAL   = B0F(1,ISEASON,ISR,ISL)
                        NITVAL = B0F(2,ISEASON,ISR,ISL)

C Interpolation day/night with transitions at SAX (sunrise) and SUX (sunset)
7034                    SIPH(ISL) = HPOL(HOUR,DAYVAL,NITVAL,
     &                          SAX,SUX,1.,1.)

C Interpolation low/middle modip with transition at 30 degrees modip
7033            SIPL(ISR) = SIPH(1) + (SIPH(2) - SIPH(1)) / DELA

C Interpolation low/high Rz12: linear from 10 to 100
        B0POL=SIPL(1)+(SIPL(2)-SIPL(1))/90.*(R-10.)
        RETURN
        END
c
C
C *********************************************************************
C ************************ EPSTEIN FUNCTIONS **************************
C *********************************************************************
C REF:  H. G. BOOKER, J. ATMOS. TERR. PHYS. 39, 619-623, 1977
C       K. RAWER, ADV. SPACE RES. 4, #1, 11-15, 1984
C *********************************************************************
C
C
        REAL FUNCTION  RLAY ( X, XM, SC, HX )
C -------------------------------------------------------- RAWER  LAYER
        Y1  = EPTR ( X , SC, HX )
        Y1M = EPTR ( XM, SC, HX )
        Y2M = EPST ( XM, SC, HX )
        RLAY = Y1 - Y1M - ( X - XM ) * Y2M / SC
        RETURN
        END
C
C
        REAL FUNCTION D1LAY ( X, XM, SC, HX )
C ------------------------------------------------------------ dLAY/dX
        D1LAY = ( EPST(X,SC,HX) - EPST(XM,SC,HX) ) /  SC
        RETURN
        END
C
C
        REAL FUNCTION D2LAY ( X, XM, SC, HX )
C ---------------------------------------------------------- d2LAY/dX2
        D2LAY = EPLA(X,SC,HX) /  (SC * SC)
        RETURN
        END
C
C
        REAL FUNCTION EPTR ( X, SC, HX )
C ------------------------------------------------------------ TRANSITION
        COMMON/ARGEXP/ARGMAX
        D1 = ( X - HX ) / SC
        IF (ABS(D1).LT.ARGMAX) GOTO 1
        IF (D1.GT.0.0) THEN
          EPTR = D1
        ELSE
          EPTR = 0.0
        ENDIF
        RETURN
1       EPTR = ALOG ( 1. + EXP( D1 ))
        RETURN
        END
C
C
        REAL FUNCTION EPST ( X, SC, HX )
C -------------------------------------------------------------- STEP
        COMMON/ARGEXP/ARGMAX
        D1 = ( X - HX ) / SC
        IF (ABS(D1).LT.ARGMAX) GOTO 1
        IF (D1.GT.0.0) THEN
          EPST = 1.
        ELSE
          EPST = 0.
        ENDIF
        RETURN
1       EPST = 1. / ( 1. + EXP( -D1 ))
        RETURN
        END
C
C
        REAL FUNCTION EPSTEP ( Y2, Y1, SC, HX, X)
C---------------------------------------------- STEP FROM Y1 TO Y2      
        EPSTEP = Y1 + ( Y2 - Y1 ) * EPST ( X, SC, HX)
        RETURN
        END
C
C
        REAL FUNCTION EPLA ( X, SC, HX )
C ------------------------------------------------------------ PEAK 
        COMMON/ARGEXP/ARGMAX
        D1 = ( X - HX ) / SC
        IF (ABS(D1).LT.ARGMAX) GOTO 1
                EPLA = 0
                RETURN  
1       D0 = EXP ( D1 )
        D2 = 1. + D0
        EPLA = D0 / ( D2 * D2 )
        RETURN
        END
c
c
        FUNCTION XE2TO5(H,HMF2,NL,HX,SC,AMP)
C----------------------------------------------------------------------
C NORMALIZED ELECTRON DENSITY (N/NMF2) FOR THE MIDDLE IONOSPHERE FROM 
C HME TO HMF2 USING LAY-FUNCTIONS.
C----------------------------------------------------------------------
        DIMENSION       HX(NL),SC(NL),AMP(NL)
        SUM = 1.0
        DO 1 I=1,NL
           YLAY = AMP(I) * RLAY( H, HMF2, SC(I), HX(I) )
           zlay=10.**ylay
1          sum=sum*zlay
        XE2TO5 = sum
        RETURN
        END
C
C
        REAL FUNCTION XEN(H,HMF2,XNMF2,HME,NL,HX,SC,AMP)
C----------------------------------------------------------------------
C ELECTRON DENSITY WITH NEW MIDDLE IONOSPHERE
C----------------------------------------------------------------------
        DIMENSION       HX(NL),SC(NL),AMP(NL)
C
        IF(H.LT.HMF2) GOTO 100
                XEN = XE1(H)
                RETURN
100     IF(H.LT.HME) GOTO 200
                XEN = XNMF2 * XE2TO5(H,HMF2,NL,HX,SC,AMP)
                RETURN
200     XEN = XE6(H)
        RETURN
        END
C
C
        SUBROUTINE VALGUL(XHI,HVB,VWU,VWA,VDP)
C --------------------------------------------------------------------- 
C   CALCULATES E-F VALLEY PARAMETERS; T.L. GULYAEVA, ADVANCES IN
C   SPACE RESEARCH 7, #6, 39-48, 1987.
C
C       INPUT:  XHI     SOLAR ZENITH ANGLE [DEGREE]
C       
C       OUTPUT: VDP     VALLEY DEPTH  (NVB/NME)
C               VWU     VALLEY WIDTH  [KM]
C               VWA     VALLEY WIDTH  (SMALLER, CORRECTED BY RAWER)
C               HVB     HEIGHT OF VALLEY BASE [KM]
C -----------------------------------------------------------------------
C
        COMMON  /CONST/UMR
C
        CS = 0.1 + COS(UMR*XHI)
        ABC = ABS(CS)
        VDP = 0.45 * CS / (0.1 + ABC ) + 0.55
        ARL = ( 0.1 + ABC + CS ) / ( 0.1 + ABC - CS)
        ZZZ = ALOG( ARL )
        VWU = 45. - 10. * ZZZ
        VWA = 45. -  5. * ZZZ
        HVB = 1000. / ( 7.024 + 0.224 * CS + 0.966 * ABC )
        RETURN
        END
C
C
        SUBROUTINE ROGUL(IDAY,XHI,SX,GRO)
C --------------------------------------------------------------------- 
C   CALCULATES RATIO H0.5/HMF2 FOR HALF-DENSITY POINT (NE(H0.5)=0.5*NMF2)
C   T.L. GULYAEVA, ADVANCES IN SPACE RESEARCH 7, #6, 39-48, 1987.
C
C       INPUT:  IDAY    DAY OF YEAR
C               XHI     SOLAR ZENITH ANGLE [DEGREE]
C       
C       OUTPUT: GRO     RATIO OF HALF DENSITY HEIGHT TO F PEAK HEIGHT
C               SX      SMOOTHLY VARYING SEASON PARAMTER (SX=1 FOR 
C                       DAY=1; SX=3 FOR DAY=180; SX=2 FOR EQUINOX)
C -----------------------------------------------------------------------
C
        common  /const1/humr,dumr
        SX = 2. - COS ( IDAY * dumr )
        XS = ( XHI - 20. * SX) / 15.
        GRO = 0.8 - 0.2 / ( 1. + EXP(XS) )
c same as gro=0.6+0.2/(1+exp(-xs))
        RETURN
        END
C
C
        SUBROUTINE LNGLSN ( N, A, B, AUS)
C --------------------------------------------------------------------
C SOLVES QUADRATIC SYSTEM OF LINEAR EQUATIONS:
C
C       INPUT:  N       NUMBER OF EQUATIONS (= NUMBER OF UNKNOWNS)
C               A(N,N)  MATRIX (LEFT SIDE OF SYSTEM OF EQUATIONS)
C               B(N)    VECTOR (RIGHT SIDE OF SYSTEM)
C
C       OUTPUT: AUS     =.TRUE.   NO SOLUTION FOUND
C                       =.FALSE.  SOLUTION IS IN  A(N,J) FOR J=1,N
C --------------------------------------------------------------------
C
        DIMENSION       A(5,5), B(5), AZV(10)
        LOGICAL         AUS
C
        NN = N - 1
        AUS = .FALSE.
        DO 1 K=1,N-1
                IMAX = K
                L    = K
                IZG  = 0
                AMAX  = ABS( A(K,K) )
110             L = L + 1
                IF (L.GT.N) GOTO 111
                HSP = ABS( A(L,K) )
                IF (HSP.LT.1.E-8) IZG = IZG + 1
                IF (HSP.LE.AMAX) GOTO 110
111             IF (ABS(AMAX).GE.1.E-10) GOTO 133
                        AUS = .TRUE.
                        RETURN
133             IF (IMAX.EQ.K) GOTO 112
                DO 2 L=K,N
                        AZV(L+1)  = A(IMAX,L)
                        A(IMAX,L) = A(K,L)
2                       A(K,L)    = AZV(L+1)
                AZV(1)  = B(IMAX)
                B(IMAX) = B(K)
                B(K)    = AZV(1)
112             IF (IZG.EQ.(N-K)) GOTO 1
                AMAX = 1. / A(K,K)
                AZV(1) = B(K) * AMAX
                DO 3 M=K+1,N
3                       AZV(M+1) = A(K,M) * AMAX
                DO 4 L=K+1,N
                        AMAX = A(L,K)
                        IF (ABS(AMAX).LT.1.E-8) GOTO 4
                        A(L,K) = 0.0
                        B(L) = B(L) - AZV(1) * AMAX
                        DO 5 M=K+1,N
5                               A(L,M) = A(L,M) - AMAX * AZV(M+1)
4               CONTINUE
1       CONTINUE
        DO 6 K=N,1,-1
                AMAX = 0.0
                IF (K.LT.N) THEN
                        DO 7 L=K+1,N
7                               AMAX = AMAX + A(K,L) * A(N,L)
                        ENDIF
                IF (ABS(A(K,K)).LT.1.E-6) THEN
                        A(N,K) = 0.0
                ELSE
                        A(N,K) = ( B(K) - AMAX ) / A(K,K)
                ENDIF
6       CONTINUE
        RETURN
        END
C
C
        SUBROUTINE LSKNM ( N, M, M0, M1, HM, SC, HX, W, X, Y, VAR, SING)
C --------------------------------------------------------------------
C   DETERMINES LAY-FUNCTIONS AMPLITUDES FOR A NUMBER OF CONSTRAINTS:
C
C       INPUT:  N       NUMBER OF AMPLITUDES ( LAY-FUNCTIONS)
C               M       NUMBER OF CONSTRAINTS
C               M0      NUMBER OF POINT CONSTRAINTS
C               M1      NUMBER OF FIRST DERIVATIVE CONSTRAINTS
C               HM      F PEAK ALTITUDE  [KM]
C               SC(N)   SCALE PARAMETERS FOR LAY-FUNCTIONS  [KM]
C               HX(N)   HEIGHT PARAMETERS FOR LAY-FUNCTIONS  [KM]
C               W(M)    WEIGHT OF CONSTRAINTS
C               X(M)    ALTITUDES FOR CONSTRAINTS  [KM]
C               Y(M)    LOG(DENSITY/NMF2) FOR CONSTRAINTS
C
C       OUTPUT: VAR(M)  AMPLITUDES
C               SING    =.TRUE.   NO SOLUTION
C ------------------------------------------------------------------------
C
        LOGICAL         SING
        DIMENSION       VAR(N), HX(N), SC(N), W(M), X(M), Y(M),
     &                  BLI(5), ALI(5,5), XLI(5,10)
C
        M01=M0+M1
        SCM=0
        DO 1 J=1,5
                BLI(J) = 0.
                DO 1 I=1,5
1                       ALI(J,I) = 0. 
        DO 2 I=1,N
                DO 3 K=1,M0
3                       XLI(I,K) = RLAY( X(K), HM, SC(I), HX(I) )
                DO 4 K=M0+1,M01
4                       XLI(I,K) = D1LAY( X(K), HM, SC(I), HX(I) )
                DO 5 K=M01+1,M
5                       XLI(I,K) = D2LAY( X(K), HM, SC(I), HX(I) )
2       CONTINUE
                DO 7 J=1,N
                DO 6 K=1,M
                        BLI(J) = BLI(J) + W(K) * Y(K) * XLI(J,K)
                        DO 6 I=1,N
6                               ALI(J,I) = ALI(J,I) + W(K) * XLI(I,K) 
     &                                  * XLI(J,K)
7       CONTINUE
        CALL LNGLSN( N, ALI, BLI, SING )
        IF (.NOT.SING) THEN
                DO 8 I=1,N
8                       VAR(I) = ALI(N,I)
                ENDIF
        RETURN
        END
C
C
        SUBROUTINE INILAY(NIGHT,XNMF2,XNMF1,XNME,VNE,HMF2,HMF1, 
     &                          HME,HV1,HV2,HHALF,HXL,SCL,AMP,IQUAL)
C-------------------------------------------------------------------
C CALCULATES AMPLITUDES FOR LAY FUNCTIONS
C D. BILITZA, DECEMBER 1988
C
C INPUT:        NIGHT   LOGICAL VARIABLE FOR DAY/NIGHT DISTINCTION
C               XNMF2   F2 PEAK ELECTRON DENSITY [M-3]
C               XNMF1   F1 PEAK ELECTRON DENSITY [M-3]
C               XNME    E  PEAK ELECTRON DENSITY [M-3]
C               VNE     ELECTRON DENSITY AT VALLEY BASE [M-3]
C               HMF2    F2 PEAK ALTITUDE [KM]
C               HMF1    F1 PEAK ALTITUDE [KM]
C               HME     E  PEAK ALTITUDE [KM]
C               HV1     ALTITUDE OF VALLEY TOP [KM]
C               HV2     ALTITUDE OF VALLEY BASE [KM]
C               HHALF   ALTITUDE OF HALF-F2-PEAK-DENSITY [KM]
C
C OUTPUT:       HXL(4)  HEIGHT PARAMETERS FOR LAY FUNCTIONS [KM] 
C               SCL(4)  SCALE PARAMETERS FOR LAY FUNCTIONS [KM]
C               AMP(4)  AMPLITUDES FOR LAY FUNCTIONS
C               IQUAL   =0 ok, =1 ok using second choice for HXL(1)
C                       =2 NO SOLUTION
C---------------------------------------------------------------  
        DIMENSION       XX(8),YY(8),WW(8),AMP(4),HXL(4),SCL(4)
        LOGICAL         SSIN,NIGHT
c
c constants --------------------------------------------------------
                NUMLAY=4
                NC1 = 2
                ALG102=ALOG10(2.)
c
c constraints: xx == height     yy == log(Ne/NmF2)    ww == weights
c -----------------------------------------------------------------
                ALOGF = ALOG10(XNMF2)
                ALOGEF = ALOG10(XNME) - ALOGF
                XHALF=XNMF2/2.
                XX(1) = HHALF
                XX(2) = HV1
                XX(3) = HV2
                XX(4) = HME
                XX(5) = HME - ( HV2 - HME )
                YY(1) = -ALG102
                YY(2) = ALOGEF
                YY(3) = ALOG10(VNE) - ALOGF
                YY(4) = ALOGEF
                YY(5) = YY(3)
                YY(7) = 0.0
                WW(2) = 1.
                WW(3) = 2.
                WW(4) = 5.
c
c geometric paramters for LAY -------------------------------------
c difference to earlier version:  HXL(3) = HV2 + SCL(3)
c
                SCL0 = 0.7 * ( 0.216 * ( HMF2 - HHALF ) + 56.8 )
                SCL(1) = 0.8 * SCL0
                SCL(2) = 10.
                SCL(3) = 9.
                SCL(4) = 6.
                HXL(3) = HV2
c
C DAY CONDITION--------------------------------------------------
c earlier tested:       HXL(2) = HMF1 + SCL(2)
c 
            IF(NIGHT) GOTO 7711
                NUMCON = 8
                HXL(1) = 0.9 * HMF2
                  HXL1T  = HHALF
                HXL(2) = HMF1
                HXL(4) = HME - SCL(4)
                XX(6) = HMF1
                XX(7) = HV2
                XX(8) = HME
                YY(8) = 0.0
                WW(5) = 1.
                WW(7) = 50.
                WW(8) = 500.
c without F-region ----------------------------------------------
                IF(XNMF1.GT.0) GOTO 100
                        HXL(2)=(HMF2+HHALF)/2.
                        YY(6) = 0.
                        WW(6) = 0.
                        WW(1) = 1.
                        GOTO 7722
c with F-region --------------------------------------------
100             YY(6) = ALOG10(XNMF1) - ALOGF
                WW(6) = 3.
                IF((XNMF1-XHALF)*(HMF1-HHALF).LT.0.0) THEN
                  WW(1)=0.5
                ELSE
                  ZET = YY(1) - YY(6)
                  WW(1) = EPST( ZET, 0.1, 0.15)
                ENDIF
                IF(HHALF.GT.HMF1) THEN
                  HFFF=HMF1
                  XFFF=XNMF1
                ELSE
                  HFFF=HHALF
                  XFFF=XHALF
                ENDIF
                GOTO 7722
c
C NIGHT CONDITION---------------------------------------------------
c different HXL,SCL values were tested including: 
c       SCL(1) = HMF2 * 0.15 - 27.1     HXL(2) = 200.   
c       HXL(2) = HMF1 + SCL(2)          HXL(3) = 140.
c       SCL(3) = 5.                     HXL(4) = HME + SCL(4)
c       HXL(4) = 105.                   
c
7711            NUMCON = 7
                HXL(1) = HHALF
                  HXL1T  = 0.4 * HMF2 + 30.
                HXL(2) = ( HMF2 + HV1 ) / 2.
                HXL(4) = HME
                XX(6) = HV2
                XX(7) = HME
                YY(6) = 0.0
                WW(1) = 1.
                WW(3) = 3.
                WW(5) = 0.5
                WW(6) = 50.
                WW(7) = 500.
                HFFF=HHALF
                XFFF=XHALF
c
C are valley-top and bottomside point compatible ? -------------
C
7722    IF((HV1-HFFF)*(XNME-XFFF).LT.0.0) WW(2)=0.5
        IF(HV1.LE.HV2+5.0) WW(2)=0.5
c
C DETERMINE AMPLITUDES-----------------------------------------
C
            NC0=NUMCON-NC1
            IQUAL=0
2299        CALL LSKNM(NUMLAY,NUMCON,NC0,NC1,HMF2,SCL,HXL,WW,XX,YY,
     &          AMP,SSIN)
                IF(IQUAL.gt.0) GOTO 1937
            IF((ABS(AMP(1)).GT.10.0).OR.(SSIN)) THEN
                IQUAL=1
                HXL(1)=HXL1T
                GOTO 2299
                ENDIF
1937        IF(SSIN) IQUAL=2
            RETURN
            END
c
c
        subroutine ioncom(h,zd,fd,fs,t,cn)
c---------------------------------------------------------------
c ion composition model
c A.D. Danilov and A.P. Yaichnikov, A New Model of the Ion
c   Composition at 75 to 1000 km for IRI, Adv. Space Res. 5, #7,
c   75-79, 107-108, 1985
c
c       h       altitude in km
c       zd      solar zenith angle in degrees
c       fd      latitude in degrees
c       fs      10.7cm solar radio flux
c       t       season (decimal month)
c       cn(1)   O+  relative density in percent
c       cn(2)   H+  relative density in percent
c       cn(3)   N+  relative density in percent
c       cn(4)   He+ relative density in percent
c       cn(5)   NO+ relative density in percent
c       cn(6)   O2+ relative density in percent
c       cn(7)   cluster ions  relative density in percent
c---------------------------------------------------------------
c
        dimension       cn(7),cm(7),hm(7),alh(7),all(7),beth(7),
     &                  betl(7),p(5,6,7),var(6),po(5,6),ph(5,6),
     &                  pn(5,6),phe(5,6),pno(5,6),po2(5,6),pcl(5,6)

        common  /argexp/argmax
        common  /const/ umr

        data po/4*0.,98.5,4*0.,320.,4*0.,-2.59E-4,2.79E-4,-3.33E-3,
     &          -3.52E-3,-5.16E-3,-2.47E-2,4*0.,-2.5E-6,1.04E-3,
     &          -1.79E-4,-4.29E-5,1.01E-5,-1.27E-3/
        data ph/-4.97E-7,-1.21E-1,-1.31E-1,0.,98.1,355.,-191.,
     &          -127.,0.,2040.,4*0.,-4.79E-6,-2.E-4,5.67E-4,
     &          2.6E-4,0.,-5.08E-3,10*0./
        data pn/7.6E-1,-5.62,-4.99,0.,5.79,83.,-369.,-324.,0.,593.,
     &          4*0.,-6.3E-5,-6.74E-3,-7.93E-3,-4.65E-3,0.,-3.26E-3,
     &          4*0.,-1.17E-5,4.88E-3,-1.31E-3,-7.03E-4,0.,-2.38E-3/
        data phe/-8.95E-1,6.1,5.39,0.,8.01,4*0.,1200.,4*0.,-1.04E-5,
     &          1.9E-3,9.53E-4,1.06E-3,0.,-3.44E-3,10*0./
        data pno/-22.4,17.7,-13.4,-4.88,62.3,32.7,0.,19.8,2.07,115.,
     &          5*0.,3.94E-3,0.,2.48E-3,2.15E-4,6.67E-3,5*0.,
     &          -8.4E-3,0.,-3.64E-3,2.E-3,-2.59E-2/
        data po2/8.,-12.2,9.9,5.8,53.4,-25.2,0.,-28.5,-6.72,120.,
     &          5*0.,-1.4E-2,0.,-9.3E-3,3.3E-3,2.8E-2,5*0.,4.25E-3,
     &          0.,-6.04E-3,3.85E-3,-3.64E-2/
        data pcl/4*0.,100.,4*0.,75.,10*0.,4*0.,-9.04E-3,-7.28E-3,
     &          2*0.,3.46E-3,-2.11E-2/


        z=zd*umr
        f=fd*umr
        DO 8 I=1,5
        DO 8 J=1,6
                p(i,j,1)=po(i,j)
                p(i,j,2)=ph(i,j)
                p(i,j,3)=pn(i,j)
                p(i,j,4)=phe(i,j)
                p(i,j,5)=pno(i,j)
                p(i,j,6)=po2(i,j)
                p(i,j,7)=pcl(i,j)
8       continue

        s=0.
        do 5 i=1,7
          do 7 j=1,6
                var(j) = p(1,j,i)*cos(z) + p(2,j,i)*cos(f) +
     &                   p(3,j,i)*cos(0.013*(300.-fs)) +
     &                   p(4,j,i)*cos(0.52*(t-6.)) + p(5,j,i)
7         continue
          cm(i)  = var(1)
          hm(i)  = var(2)
          all(i) = var(3)
          betl(i)= var(4)
          alh(i) = var(5)
          beth(i)= var(6)
          hx=h-hm(i)
          if(hx) 1,2,3
1               arg = hx * (hx * all(i) + betl(i))
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
                goto 4
2               cn(i) = cm(i)
                goto 4
3               arg = hx * (hx * alh(i) + beth(i))
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
4         continue
          if(cn(i).LT.0.005*cm(i)) cn(i)=0.
          if(cn(i).GT.cm(i)) cn(i)=cm(i)
          s=s+cn(i)
5       continue
        do 6 i=1,7
6               cn(i)=cn(i)/s*100.
        return
        end
c
c

      Subroutine ionco2(hei,xhi,it,F,R1,R2,R3,R4)
*----------------------------------------------------------------
*     INPUT DATA :
*      hei  -  altitude in km
*      xhi  -  solar zenith angle in degree
*      it   -  season (month)
*      F    -  10.7cm solar radio flux
*     OUTPUT DATA :
*     R1 -  NO+ concentration (in percent)
*     R2 -  O2+ concentration (in percent) 
*     R3 -  Cb+ concentration (in percent) 
*     R4 -  O+  concentration (in percent) 
*
*  A.D. Danilov and N.V. Smirnova, Improving the 75 to 300 km ion 
*  composition model of the IRI, Adv. Space Res. 15, #2, 171-177, 1995.
*
*-----------------------------------------------------------------
      dimension j1ms70(7),j2ms70(7),h1s70(13,7),h2s70(13,7),
     *       R1ms70(13,7),R2ms70(13,7),rk1ms70(13,7),rk2ms70(13,7),
     *       j1ms140(7),j2ms140(7),h1s140(13,7),h2s140(13,7), 
     *       R1ms140(13,7),R2ms140(13,7),rk1ms140(13,7),rk2ms140(13,7),
     *       j1mw70(7),j2mw70(7),h1w70(13,7),h2w70(13,7),
     *       R1mw70(13,7),R2mw70(13,7),rk1mw70(13,7),rk2mw70(13,7),
     *       j1mw140(7),j2mw140(7),h1w140(13,7),h2w140(13,7), 
     *       R1mw140(13,7),R2mw140(13,7),rk1mw140(13,7),rk2mw140(13,7),
     *       j1mr70(7),j2mr70(7),h1r70(13,7),h2r70(13,7),
     *       R1mr70(13,7),R2mr70(13,7),rk1mr70(13,7),rk2mr70(13,7),
     *       j1mr140(7),j2mr140(7),h1r140(13,7),h2r140(13,7), 
     *       R1mr140(13,7),R2mr140(13,7),rk1mr140(13,7),rk2mr140(13,7)
      data j1ms70/11,11,10,10,11,9,11/ 
      data j2ms70/13,11,10,11,11,9,11/
      data h1s70/75,85,90,95,100,120,130,200,220,250,270,0,0,
     *        75,85,90,95,100,120,130,200,220,250,270,0,0, 
     *        75,85,90,95,100,115,200,220,250,270,0,0,0,
     *        75,80,95,100,120,140,200,220,250,270,0,0,0,
     *        75,80,95,100,120,150,170,200,220,250,270,0,0,
     *        75,80,95,100,140,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,110,145,200,220,250,270,0,0/
      data h2s70/75,80,90,95,100,120,130,140,150,200,220,250,270,
     *        75,80,90,95,100,120,130,200,220,250,270,0,0, 
     *        75,80,90,95,100,115,200,220,250,270,0,0,0,
     *        75,80,95,100,120,140,150,200,220,250,270,0,0,
     *        75,80,95,100,120,150,170,200,220,250,270,0,0,
     *        75,80,95,100,140,200,220,250,270,0,0,0,0,
     *        75,80,90,95,100,110,145,200,220,250,270,0,0/
      data R1ms70/6,30,60,63,59,59,66,52,20,4,2,0,0,
     *         6,30,60,63,69,62,66,52,20,4,2,0,0, 
     *         6,30,60,63,80,68,53,20,4,2,0,0,0,
     *         4,10,60,85,65,65,52,25,12,4,0,0,0, 
     *         4,10,60,89,72,60,60,52,30,20,10,0,0, 
     *         4,10,60,92,68,54,40,25,13,0,0,0,0, 
     *         1,8,20,60,95,93,69,65,45,30,20,0,0/ 
      data R2ms70/4,10,30,32,41,41,32,29,34,28,15,3,1,
     *         4,10,30,32,31,38,32,28,15,3,1,0,0,
     *         4,10,30,32,20,32,28,15,3,1,0,0,0,
     *         2,6,30,15,35,30,34,26,19,8,3,0,0,
     *         2,6,30,11,28,38,29,29,25,12,5,0,0,
     *         2,6,30,8,32,30,20,14,8,0,0,0,0,
     *         1,2,10,20,5,7,31,23,18,15,10,0,0/ 
      data rk1ms70/2.4,6.,.6,-.8,0,.7,-.2,-1.6,-.533,-.1,-.067,0,0,
     *         2.4,6.,.6,1.2,-.35,.4,-.2,-1.6,-.533,-.1,-.067,0,0, 
     *         2.4,6.,.6,3.4,-.8,-.176,-1.65,-.533,-.1,-.067,0,0,0,
     *         1.2,3.333,5.,-1.,0,-.216,-1.35,-.433,-.4,-.1,0,0,0,
     *         1.2,3.333,5.8,-.85,-.4,0,-.267,-1.1,-.333,-.4,-.2,0,0, 
     *         1.2,3.333,6.4,-.6,-.233,-.7,-.5,-.6,-.267,0,0,0,0, 
     *         1.4,2.4,4.,7.,-.2,-.686,-.072,-1.,-.5,-.5,-.5,0,0/
      data rk2ms70/1.2,2.,.4,1.8,0,-.9,-.3,.5,-.12,-.65,-.4,-.1,-.033,
     *         1.2,2.,.4,-.2,.35,-.6,-.057,-.65,-.4,-.1,-.033,0,0,
     *         1.2,2.,.4,-2.4,.8,-.047,-.65,-.4,-.1,-.033,0,0,0,
     *         .8,1.6,-3.,1.,-.25,.4,-.16,-.35,-.367,-.25,-.1,0,0,
     *         .8,1.6,-3.8,.85,.333,-.45,0,-.2,-.433,-.35,-.1,0,0,
     *         .8,1.6,-4.4,.6,-.033,-.5,-.2,-.3,-.2,0,0,0,0,
     *         .2,.8,2.,-3.,.2,.686,-.145,-.25,-.1,-.25,-.2,0,0/
      data j1ms140/11,11,10,10,9,9,12/ 
      data j2ms140/11,11,10,9,10,10,12/
      data h1s140/75,85,90,95,100,120,130,140,200,220,250,0,0,
     *        75,85,90,95,100,120,130,140,200,220,250,0,0,
     *        75,85,90,95,100,120,140,200,220,250,0,0,0,
     *        75,80,95,100,120,140,200,220,250,270,0,0,0,
     *        75,80,95,100,120,200,220,250,270,0,0,0,0,
     *        75,80,95,100,130,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,110,140,180,200,220,250,270,0/
      data h2s140/75,80,90,95,100,120,130,155,200,220,250,0,0,
     *        75,80,90,95,100,120,130,160,200,220,250,0,0,
     *        75,80,90,95,100,120,165,200,220,250,0,0,0,
     *        75,80,95,100,120,180,200,250,270,0,0,0,0,
     *        75,80,95,100,120,160,200,220,250,270,0,0,0,
     *        75,80,95,100,130,160,200,220,250,270,0,0,0,
     *        75,80,90,95,100,110,140,180,200,220,250,270,0/
      data R1ms140/6,30,60,63,59,59,66,66,38,14,1,0,0,
     *         6,30,60,63,69,62,66,66,38,14,1,0,0,
     *         6,30,60,63,80,65,65,38,14,1,0,0,0,
     *         4,10,60,85,66,66,38,22,9,1,0,0,0,
     *         4,10,60,89,71,42,26,17,10,0,0,0,0,
     *         4,10,60,93,71,48,35,22,10,0,0,0,0,
     *         1,8,20,60,95,93,72,60,58,40,26,13,0/ 
      data R2ms140/4,10,30,32,41,41,30,30,10,6,1,0,0,
     *         4,10,30,32,31,38,31,29,9,6,1,0,0,
     *         4,10,30,32,20,35,26,9,6,1,0,0,0,
     *         2,6,30,15,34,24,10,5,1,0,0,0,0,
     *         2,6,30,11,28,37,21,14,8,5,0,0,0,
     *         2,6,30,7,29,36,29,20,13,5,0,0,0,
     *         1,2,10,20,5,7,28,32,28,20,14,7,0/ 
      data rk1ms140/2.4,6.,.6,-.8,0,.7,0,-.467,-1.2,-.433,0,0,0,
     *         2.4,6.,.6,1.2,-.35,.4,0,-.467,-1.2,-.433,0,0,0,    
     *         2.4,6.,.6,3.4,-.75,0,-.45,-1.2,-.433,0,0,0,0,
     *         1.2,3.333,5.,-.95,0,-.467,-.8,-.433,-.4,0,0,0,0,
     *         1.2,3.333,5.8,-.9,-.363,-.8,-.3,-.35,-.3,0,0,0,0,
     *         1.2,3.333,6.6,-.733,-.329,-.65,-.433,-.6,-.267,0,0,0,0,
     *         1.4,2.4,4.,7.,-.2,-.7,-.3,-.1,-.9,-.467,-.65,-.333,0/
      data rk2ms140/1.2,2.,.4,1.8,0,-1.1,0,-.444,-.2,-.166,0,0,0,
     *         1.2,2.,.4,-.2,.35,-.7,-.067,-.5,-.15,-.166,0,0,0,
     *         1.2,2.,.4,-2.4,.75,-.2,-.486,-.15,-.166,0,0,0,0,
     *         .8,1.6,-3.,.95,-.167,-.7,-.1,-.2,0,0,0,0,0,
     *         .8,1.6,-3.8,.85,.225,-.4,-.35,-.2,-.15,-.133,0,0,0,
     *         .8,1.6,-4.6,.733,.233,-.175,-.45,-.233,-.4,-.1,0,0,0, 
     *         .2,.8,2.,-3.,.2,.7,.1,-.2,-.4,-.2,-.35,-.167,0/
      data j1mr70/12,12,12,9,10,11,13/ 
      data j2mr70/9,9,10,13,12,11,11/
      data h1r70/75,80,90,95,100,120,140,180,200,220,250,270,0,
     *        75,80,90,95,100,120,145,180,200,220,250,270,0, 
     *        75,80,90,95,100,120,145,180,200,220,250,270,0,  
     *        75,95,100,110,140,180,200,250,270,0,0,0,0,
     *        75,95,125,150,185,195,200,220,250,270,0,0,0,
     *        75,95,100,150,160,170,190,200,220,250,270,0,0,
     *        75,80,85,95,100,140,160,170,190,200,220,250,270/
      data h2r70/75,95,100,120,180,200,220,250,270,0,0,0,0,
     *        75,95,100,120,180,200,220,250,270,0,0,0,0, 
     *        75,95,100,120,130,190,200,220,250,270,0,0,0, 
     *        75,80,85,95,100,110,130,180,190,200,220,250,270,
     *        75,80,85,95,100,125,150,190,200,220,250,270,0,
     *        75,80,85,95,100,150,190,200,220,250,270,0,0, 
     *        75,85,95,100,140,180,190,200,220,250,270,0,0/
      data R1mr70/13,17,57,57,30,53,58,38,33,14,6,2,0,
     *         13,17,57,57,37,56,56,38,33,14,6,2,0, 
     *         13,17,57,57,47,58,55,37,33,14,6,2,0, 
     *         5,65,54,58,58,38,33,9,1,0,0,0,0, 
     *         5,65,65,54,40,40,45,26,17,10,0,0,0,    
     *         5,65,76,56,57,48,44,51,35,22,10,0,0, 
     *         3,11,35,75,90,65,63,54,54,50,40,26,13/ 
      data R2mr70/7,43,70,47,15,17,10,4,0,0,0,0,0,
     *         7,43,63,44,17,17,10,4,0,0,0,0,0, 
     *         7,43,53,42,42,13,17,10,4,0,0,0,0,
     *         3,5,26,34,46,42,41,23,16,16,10,1,0,
     *         3,5,26,34,35,35,42,25,22,14,8,5,0, 
     *         3,5,26,34,24,41,31,26,20,13,5,0,0,
     *         3,15,15,10,35,35,30,34,20,14,7,0,0/ 
      data rk1mr70/.8,4.,0,-5.4,1.15,.25,-.5,-.25,-.95,-.267,-.2,
     *             -.067,0,
     *         .8,4.,0,-4.,.95,0,-.514,-.25,-.95,-.267,-.2,-.067,0, 
     *         .8,4.,0,-2.,.55,-.12,-.514,-.2,-.95,-.267,-.2,-.067,0, 
     *         3.,-2.2,.4,0,-.5,-.25,-.48,-.4,-.033,0,0,0,0,   
     *         3.,0,-.44,-.466,0,1.0,-.95,-.3,-.35,-.3,0,0,0, 
     *         3.,2.2,-.4,0.1,-.9,-.2,.7,-.8,-.433,-.6,-.267,0,0, 
     *         1.6,4.8,4.,3.,-.625,-.1,-.9,0,-.4,-.5,-.467,-.65,-.3/
      data rk2mr70/1.8,5.4,-1.15,-.533,.1,-.35,-.2,-.2,0,0,0,0,0,
     *         1.8,4.,-.95,-.45,0,-.35,-.2,-.2,0,0,0,0,0,
     *         1.8,2.,-.55,0,-.483,.4,-.35,-.2,-.2,0,0,0,0,    
     *         .4,4.2,.8,2.4,-.4,-.05,-.36,-.7,0,-.3,-.3,-.05,0,
     *         .4,4.2,.8,.2,0,.28,-.425,-.3,-.4,-.2,-.15,-.133,0,  
     *         .4,4.2,.8,-2.,.34,-.25,-.5,-.3,-.233,-.4,-.1,0,0,  
     *         1.2,0,-1.,.625,0,-.5,.4,-.7,-.2,-.35,-.167,0,0/
      data j1mr140/12,12,11,12,9,9,13/ 
      data j2mr140/10,9,10,12,13,13,12/
      data h1r140/75,80,90,95,100,115,130,145,200,220,250,270,0,
     *        75,80,90,95,100,110,120,145,200,220,250,270,0, 
     *        75,80,90,95,100,115,150,200,220,250,270,0,0,
     *        75,95,100,120,130,140,150,190,200,220,250,270,0,
     *        75,95,120,150,190,200,220,250,270,0,0,0,0,
     *        75,95,100,145,190,200,220,250,270,0,0,0,0,
     *        75,80,85,95,100,120,160,170,190,200,220,250,270/
      data h2r140/75,95,100,115,130,175,200,220,250,270,0,0,0,
     *        75,95,100,110,175,200,220,250,270,0,0,0,0, 
     *        75,95,100,115,130,180,200,220,250,270,0,0,0, 
     *        75,80,85,95,100,120,130,190,200,220,250,270,0,
     *        75,80,85,95,100,120,140,160,190,200,220,250,270,
     *        75,80,85,95,100,145,165,180,190,200,220,250,270,
     *        75,85,95,100,120,145,170,190,200,220,250,270,0/
      data R1mr140/13,17,57,57,28,51,56,56,12,8,1,0,0,
     *         13,17,57,57,36,46,55,56,10,8,1,0,0, 
     *         13,17,57,57,46,56,55,12,8,1,0,0,0,
     *         5,65,54,59,56,56,53,23,16,13,3,1,0,
     *         5,65,65,54,29,16,16,10,2,0,0,0,0,
     *         5,65,76,58,36,25,20,12,7,0,0,0,0,
     *         3,11,35,75,91,76,58,49,45,32,28,20,12/ 
      data R2mr140/7,43,72,49,44,14,7,4,1,0,0,0,0,
     *         7,43,64,51,14,7,4,1,0,0,0,0,0, 
     *         7,43,54,44,44,13,7,4,1,0,0,0,0,
     *         3,5,26,34,46,41,44,9,11,7,2,1,0,
     *         3,5,26,34,35,35,40,40,16,14,9,5,2, 
     *         3,5,26,34,24,40,40,32,19,20,10,7,3,
     *         3,15,15,9,24,35,40,28,28,20,10,8,0/ 
      data rk1mr140/.8,4.,0,-5.8,1.533,.333,0,-.8,-.2,-.233,-.05,0,0,
     *         .8,4.,0,-4.2,1.3,.6,.04,-.836,-.1,-.233,-.05,0,0,
     *         .8,4.,0,-2.2,.667,-.029,-.86,-.2,-.233,-.05,0,0,0,         
     *         3.,-2.2,.25,-.3,0,-.3,-.75,-.7,-.15,-.333,-.1,-.033,0,
     *         3.,0,-.367,-.625,-1.3,0,-.2,-.4,-.067,0,0,0,0, 
     *         3.,2.2,-.4,-.489,-1.1,-.25,-.267,-.25,-.2,0,0,0,0,
     *         1.6,4.8,4.,3.2,-.75,-.45,-.9,-.2,-1.3,-.2,-.267,-.4,-.3/
      data rk2mr140/1.8,5.8,-1.533,-.333,-.667,-.28,-.15,-.1,-.05,
     *              0,0,0,0,
     *         1.8,4.2,-1.3,-.569,-.28,-.15,-.1,-.05,0,0,0,0,0,
     *         1.8,2.2,-.667,0,-.62,-.3,-.15,-.1,-.05,0,0,0,0,    
     *         .4,4.2,.8,2.4,-.25,.3,-.583,.2,-.2,-.167,-.05,-.033,0,
     *         .4,4.2,.8,.02,0,.25,0,-.6,-.2,-.25,-.133,-.15,-.067,  
     *         .4,4.2,.8,-2.,.356,0,-.533,-1.3,.1,-.5,-.1,-.2,-.1,  
     *         1.2,0,-1.2,.75,.44,.2,-.6,0,-.4,-.333,-.1,-.2,0/
      data j1mw70/13,13,13,13,9,8,9/ 
      data j2mw70/10,10,11,11,9,8,11/
      data h1w70/75,80,85,95,100,110,125,145,180,200,220,250,270,
     *        75,80,85,95,100,110,120,150,180,200,220,250,270,
     *        75,80,85,95,100,110,120,155,180,200,220,250,270,
     *        75,80,90,100,110,120,140,160,190,200,220,250,270, 
     *        75,80,90,110,150,200,220,250,270,0,0,0,0,
     *        75,80,90,100,150,200,250,270,0,0,0,0,0,
     *        75,80,90,100,120,130,140,200,270,0,0,0,0/
      data h2w70/75,90,95,100,110,125,190,200,250,270,0,0,0,
     *        75,90,95,100,110,125,190,200,250,270,0,0,0, 
     *        75,90,95,100,110,120,145,190,200,250,270,0,0,
     *        75,80,95,100,110,120,150,200,220,250,270,0,0,
     *        75,80,90,95,110,145,200,250,270,0,0,0,0,
     *        75,80,90,100,140,150,200,250,0,0,0,0,0,
     *        75,80,85,90,100,120,130,140,160,200,270,0,0/ 
      data R1mw70/28,35,65,65,28,44,46,50,25,25,10,5,0,
     *         28,35,65,65,36,49,47,47,25,25,10,5,0,
     *         28,35,65,65,48,54,51,43,25,25,10,5,0,
     *         16,24,66,54,58,50,50,38,25,25,10,5,0, 
     *         16,24,66,66,46,30,20,6,3,0,0,0,0,  
     *         16,24,66,76,49,32,12,7,0,0,0,0,0,      
     *         6,19,67,91,64,68,60,40,12,0,0,0,0/ 
      data R2mw70/5,35,35,72,56,54,12,12,2,0,0,0,0,
     *         5,35,35,64,51,53,12,12,2,0,0,0,0, 
     *         5,35,35,52,46,49,41,12,12,2,0,0,0,
     *         4,10,40,46,42,50,41,12,7,2,0,0,0,
     *         4,10,30,34,34,51,14,4,2,0,0,0,0,
     *         4,10,30,24,45,48,20,5,0,0,0,0,0,
     *         2,6,17,23,9,36,32,40,40,20,6,0,0/ 
      data rk1mw70/1.4,6.,0,-7.4,1.6,.133,.2,-.714,0,-.75,-.167,-.25,0,
     *         1.4,6.,0,-5.8,1.3,-.2,0,-.733,0,-.75,-.167,-.25,0,
     *         1.4,6.,0,-3.4,.6,-.3,-.229,-.72,0,-.75,-.167,-.25,0,
     *         1.6,4.2,-1.2,.4,-.8,0,-.6,-.433,0,-.75,-.167,-.25,0,   
     *         1.6,4.2,0,-.5,-.32,-.5,-.467,-.15,-.1,0,0,0,0,
     *         1.6,4.2,1.,-.54,-.34,-.4,-.25,-.2,0,0,0,0,0,      
     *         2.6,4.8,2.4,-1.35,.4,-.8,-.333,-.4,-.3,0,0,0,0/
      data rk2mw70/2.,0,7.4,-1.6,-.133,-.646,0,-.2,-.1,0,0,0,0,
     *         2.,0,5.8,-1.3,.133,-.631,0,-.2,-.1,0,0,0,0,    
     *         2.,0,3.4,-.6,.3,-.32,-.644,0,-.2,-.1,0,0,0,
     *         1.2,2.,1.2,-.4,.8,-.3,-.58,-.25,-.167,-.1,0,0,0,
     *         1.2,2.,.8,0,.486,-.673,-.2,-.1,-.066,0,0,0,0,  
     *         1.2,2.,-.6,.525,.3,-.56,-.3,-.1,0,0,0,0,0,  
     *         .8,2.2,1.2,-1.4,1.35,-.4,.8,0,-.5,-.2,-.167,0,0/
      data j1mw140/12,11,11,11,11,10,12/ 
      data j2mw140/10,11,11,11,11,10,12/
      data h1w140/75,80,85,95,100,110,125,145,190,200,220,250,0,
     *        75,80,85,95,100,110,120,150,190,220,250,0,0,
     *        75,80,85,95,100,110,120,155,190,220,250,0,0,
     *        75,80,90,100,110,120,140,160,190,220,250,0,0,
     *        75,80,90,110,150,160,190,200,220,250,270,0,0,
     *        75,80,90,100,150,160,190,200,250,270,0,0,0,
     *        75,80,90,100,120,130,140,160,190,200,250,270,0/
      data h2w140/75,90,95,100,110,125,190,200,220,250,0,0,0,
     *        75,90,95,100,110,120,125,190,200,220,250,0,0,
     *        75,90,95,100,110,120,145,190,200,220,250,0,0,  
     *        75,80,95,100,110,120,150,190,200,220,250,0,0,
     *        75,80,90,95,110,145,190,200,220,250,270,0,0,
     *        75,80,90,100,140,150,200,220,250,270,0,0,0,
     *        75,80,85,90,100,120,130,140,160,180,200,220,0/ 
      data R1mw140/28,35,65,65,28,44,46,50,9,6,2,0,0,
     *         28,35,65,65,36,49,47,47,8,2,0,0,0,
     *         28,35,65,65,48,54,51,43,8,2,0,0,0,
     *         16,24,66,54,58,50,50,42,8,2,0,0,0, 
     *         16,24,66,66,46,49,9,10,7,2,0,0,0,
     *         16,24,66,76,49,54,10,14,4,1,0,0,0,
     *         6,19,67,91,64,68,60,58,11,20,5,2,0/ 
      data R2mw140/5,35,35,72,56,54,5,5,1,0,0,0,0,
     *         5,35,35,64,51,53,53,5,5,1,0,0,0,
     *         5,35,35,52,46,49,41,5,5,1,0,0,0,    
     *         4,10,40,46,42,50,41,5,5,1,0,0,0,   
     *         4,10,30,34,34,51,10,5,3,1,0,0,0,  
     *         4,10,30,24,45,48,4,2,1,0,0,0,0,
     *         2,6,17,23,9,36,32,40,39,29,1,0,0/ 
      data rk1mw140/1.4,6.,0,-7.4,1.6,.133,.2,-.911,-.3,-.2,-.066,0,0,
     *         1.4,6.,0,-5.8,1.3,-.2,0,-.975,-.2,-.066,0,0,0,
     *         1.4,6.,0,-3.4,.6,-.3,-.229,-1.,-.2,-.066,0,0,0,
     *         1.6,4.2,-1.2,.4,-.8,0,-.4,-1.133,-.2,-.066,0,0,0,
     *         1.6,4.2,0,-.5,.3,-1.133,.1,-.15,-.166,-.1,0,0,0,
     *         1.6,4.2,1.,-.54,.5,-1.466,.4,-.2,-.15,-.0333,0,0,0,
     *         2.6,4.8,2.4,-1.35,.4,-.8,-.1,-1.566,.9,-.3,-.15,-.05,0/
      data rk2mw140/2.,0,7.4,-1.6,-.133,-.754,0,-.2,-.033,0,0,0,0,
     *         2.,0,5.8,-1.3,.2,0,-.738,0,-.2,-.033,0,0,0,
     *         2.,0,3.4,-.6,.3,-.32,-.8,0,-.2,-.033,0,0,0,
     *         1.2,2.,1.2,-.4,.8,-.3,-.9,0,-.2,-.033,0,0,0,
     *         1.2,2.,.8,0,.486,-.911,-.5,-.1,-.066,-.05,0,0,0,
     *         1.2,2.,-.6,.525,.3,-.88,-.1,-.033,-.05,0,0,0,0,
     *         .8,2.2,1.2,-1.4,1.35,-.4,.8,-.05,-.5,-1.4,-.05,0,0/

        h = hei
        z = xhi

         if(z.lt.20)z=20
         if(z.gt.90)z=90
        if((it.eq.1).or.(it.eq.2).or.(it.eq.11).or.(it.eq.12))then
         if(f.lt.140)then
           Call aprok(j1mw70,j2mw70,h1w70,h2w70,R1mw70,R2mw70,        
     *                rk1mw70,rk2mw70,h,z,R1,R2) 
           R170=R1
           R270=R2
         endif
         if(f.gt.70)then
           Call aprok(j1mw140,j2mw140,h1w140,h2w140,R1mw140,R2mw140,        
     *                rk1mw140,rk2mw140,h,z,R1,R2) 
           R1140=R1
           R2140=R2
         endif
         if((f.gt.70).and.(f.lt.140))then
           R1=R170+(R1140-R170)*(f-70)/70
           R2=R270+(R2140-R270)*(f-70)/70
         endif
        endif
        if((it.eq.5).or.(it.eq.6).or.(it.eq.7).or.(it.eq.8))then
         if(f.lt.140)then
           Call aprok(j1ms70,j2ms70,h1s70,h2s70,R1ms70,R2ms70,        
     *                rk1ms70,rk2ms70,h,z,R1,R2) 
           R170=R1
           R270=R2
         endif
         if(f.gt.70)then
           Call aprok(j1ms140,j2ms140,h1s140,h2s140,R1ms140,R2ms140,        
     *                rk1ms140,rk2ms140,h,z,R1,R2) 
           R1140=R1
           R2140=R2
         endif
         if((f.gt.70).and.(f.lt.140))then
           R1=R170+(R1140-R170)*(f-70)/70
           R2=R270+(R2140-R270)*(f-70)/70
         endif
        endif
        if((it.eq.3).or.(it.eq.4).or.(it.eq.9).or.(it.eq.10))then
         if(f.lt.140)then
           Call aprok(j1mr70,j2mr70,h1r70,h2r70,R1mr70,R2mr70,        
     *                rk1mr70,rk2mr70,h,z,R1,R2) 
           R170=R1
           R270=R2
         endif
         if(f.gt.70)then
           Call aprok(j1mr140,j2mr140,h1r140,h2r140,R1mr140,R2mr140,
     *                rk1mr140,rk2mr140,h,z,R1,R2) 
           R1140=R1
           R2140=R2
         endif
         if((f.gt.70).and.(f.lt.140))then
           R1=R170+(R1140-R170)*(f-70)/70
           R2=R270+(R2140-R270)*(f-70)/70
         endif
        endif
        R3=0
        R4=0
        if (h.lt.100) R3=100-(R1+R2)
        if (h.ge.100) R4=100-(R1+R2)
         if(R3.lt.0) R3=0
         if(R4.lt.0) R4=0
        R1=ANINT(R1)
        R2=ANINT(R2)
        R3=ANINT(R3)
        R4=ANINT(R4)
 300   continue
        end
c
c
      Subroutine aprok(j1m,j2m,h1,h2,R1m,R2m,rk1m,rk2m,hei,xhi,R1,R2)
c----------------------------------------------------------------- 
      dimension   zm(7),j1m(7),j2m(7),h1(13,7),h2(13,7),R1m(13,7),
     *            R2m(13,7),rk1m(13,7),rk2m(13,7)
      data        zm/20,40,60,70,80,85,90/
      
        h=hei
        z=xhi

         j1=1
         j2=1
         i1=1
       do 1 i=1,7
         i1=i
        if(z.eq.zm(i)) j1=0
        if(z.le.zm(i)) goto 11
 1     continue
 11    continue
          i2=1
         do 2 i=2,j1m(i1)
          i2=i-1
          if(h.lt.h1(i,i1)) goto 22
          i2=j1m(i1)
 2       continue
 22      continue
          i3=1
         do 3 i=2,j2m(i1)
          i3=i-1
          if(h.lt.h2(i,i1)) goto 33
          i3=j2m(i1)
 3       continue
 33      continue
        R01=R1m(i2,i1)
        R02=R2m(i3,i1)
        rk1=rk1m(i2,i1)
        rk2=rk2m(i3,i1)
        h01=h1(i2,i1)
        h02=h2(i3,i1)
        R1=R01+rk1*(h-h01)
        R2=R02+rk2*(h-h02)
        if(j1.eq.1)then
          j1=0
          j2=0
          i1=i1-1
          R11=R1
          R12=R2
          goto 11
        endif
        if(j2.eq.0)then
          rk=(z-zm(i1))/(zm(i1+1)-zm(i1)) 
          R1=R1+(R11-R1)*rk
          R2=R2+(R12-R2)*rk
        endif
       end
c
c
        subroutine ioncom_new(hei,xhia,slati,covi,zmos,dion)
c-------------------------------------------------------
c       see IONCO1 for explanation of i/o parameters
c       NOTE: xhi,xlati are in DEGREES !!!
c       NOTE: zmosea is the seasonal northern month, so 
c             for the southern hemisphere zmosea=month+6
c-------------------------------------------------------
        dimension       dion(7),dup(4),diont(7)
        common  /const/ umr

        do 1122 i=1,7
1122                diont(i)=0.

        h = hei
        xhi = xhia
        xlati = slati
        cov = covi
        zmosea = zmos
        if (h.gt.300.) then
                call ionco1(h,xhi,xlati,cov,zmosea,dup)
                diont(1)=dup(1)
                diont(2)=dup(2)
                diont(3)=dup(3)
                diont(4)=dup(4)
        else
                monsea=int(zmosea)
                call ionco2(h,xhi,monsea,cov,rno,ro2,rcl,ro)
                diont(5)=rno
                diont(6)=ro2
                diont(7)=rcl
                diont(1)=ro
        endif
        do 1 i=1,7
                dion(i)=diont(i)
1               continue
        return
        end
c
c
        subroutine ionco1(h,zd,fd,fs,t,cn)
c---------------------------------------------------------------
c ion composition model
c   A.D. Danilov and A.P. Yaichnikov, A New Model of the Ion
c   Composition at 75 to 1000 km for IRI, Adv. Space Res. 5, #7,
c   75-79, 107-108, 1985
c
c       h       altitude in km
c       zd      solar zenith angle in degrees
c       fd      latitude in degrees
c       fs      10.7cm solar radio flux
c       t       season (decimal month) 
c       cn(1)   O+  relative density in percent
c       cn(2)   H+  relative density in percent
c       cn(3)   N+  relative density in percent
c       cn(4)   He+ relative density in percent
c Please note: molecular ions are now computed in IONCO2
c       [cn(5)   NO+ relative density in percent
c       [cn(6)   O2+ relative density in percent
c       [cn(7)   cluster ions  relative density in percent
c---------------------------------------------------------------
c
c        dimension       cn(7),cm(7),hm(7),alh(7),all(7),beth(7),
c     &                  betl(7),p(5,6,7),var(6),po(5,6),ph(5,6),
c     &                  pn(5,6),phe(5,6),pno(5,6),po2(5,6),pcl(5,6)
        dimension       cn(4),cm(4),hm(4),alh(4),all(4),beth(4),
     &                  betl(4),p(5,6,4),var(6),po(5,6),ph(5,6),
     &                  pn(5,6),phe(5,6)

        common  /argexp/argmax
        common  /const/ umr
        data po/4*0.,98.5,4*0.,320.,4*0.,-2.59E-4,2.79E-4,-3.33E-3,
     &          -3.52E-3,-5.16E-3,-2.47E-2,4*0.,-2.5E-6,1.04E-3,
     &          -1.79E-4,-4.29E-5,1.01E-5,-1.27E-3/
        data ph/-4.97E-7,-1.21E-1,-1.31E-1,0.,98.1,355.,-191.,
     &          -127.,0.,2040.,4*0.,-4.79E-6,-2.E-4,5.67E-4,
     &          2.6E-4,0.,-5.08E-3,10*0./
        data pn/7.6E-1,-5.62,-4.99,0.,5.79,83.,-369.,-324.,0.,593.,
     &          4*0.,-6.3E-5,-6.74E-3,-7.93E-3,-4.65E-3,0.,-3.26E-3,
     &          4*0.,-1.17E-5,4.88E-3,-1.31E-3,-7.03E-4,0.,-2.38E-3/
        data phe/-8.95E-1,6.1,5.39,0.,8.01,4*0.,1200.,4*0.,-1.04E-5,
     &          1.9E-3,9.53E-4,1.06E-3,0.,-3.44E-3,10*0./ 
c       data pno/-22.4,17.7,-13.4,-4.88,62.3,32.7,0.,19.8,2.07,115.,
c    &          5*0.,3.94E-3,0.,2.48E-3,2.15E-4,6.67E-3,5*0.,
c    &          -8.4E-3,0.,-3.64E-3,2.E-3,-2.59E-2/
c       data po2/8.,-12.2,9.9,5.8,53.4,-25.2,0.,-28.5,-6.72,120.,
c    &          5*0.,-1.4E-2,0.,-9.3E-3,3.3E-3,2.8E-2,5*0.,4.25E-3,
c    &          0.,-6.04E-3,3.85E-3,-3.64E-2/
c       data pcl/4*0.,100.,4*0.,75.,10*0.,4*0.,-9.04E-3,-7.28E-3,
c    &          2*0.,3.46E-3,-2.11E-2/

        z=zd*umr
        f=fd*umr

        DO 8 I=1,5
        DO 8 J=1,6
                p(i,j,1)=po(i,j)
                p(i,j,2)=ph(i,j)
                p(i,j,3)=pn(i,j)
                p(i,j,4)=phe(i,j)
c               p(i,j,5)=pno(i,j)
c               p(i,j,6)=po2(i,j)
c               p(i,j,7)=pcl(i,j)
8       continue

        s=0.
c       do 5 i=1,7
        do 5 i=1,4
          do 7 j=1,6
                var(j) = p(1,j,i)*cos(z) + p(2,j,i)*cos(f) +
     &                   p(3,j,i)*cos(0.013*(300.-fs)) +
     &                   p(4,j,i)*cos(0.52*(t-6.)) + p(5,j,i)
7         continue
          cm(i)  = var(1)
          hm(i)  = var(2)
          all(i) = var(3)
          betl(i)= var(4)
          alh(i) = var(5)
          beth(i)= var(6)
          hx=h-hm(i)
          if(hx) 1,2,3
1               arg = hx * (hx * all(i) + betl(i)) 
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
                goto 4
2               cn(i) = cm(i)
                goto 4
3               arg = hx * (hx * alh(i) + beth(i)) 
                cn(i) = 0.
                if(arg.gt.-argmax) cn(i) = cm(i) * exp( arg )
4         continue
          if(cn(i).LT.0.005*cm(i)) cn(i)=0.
          if(cn(i).GT.cm(i)) cn(i)=cm(i)
          s=s+cn(i)
5       continue
c       do 6 i=1,7
        do 6 i=1,4
6               cn(i)=cn(i)/s*100.
        return
        end
c
c
           subroutine tcon(yr,mm,day,idn,rz,ig,rsn,nmonth)
c----------------------------------------------------------------
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c output:       rz(3)           12-month-smoothed solar sunspot number
c               ig(3)           12-month-smoothed IG index
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c rz(1), ig(1) contain the indices for the month mm and rz(2), ig(2)
c for the previous month (if day less than 15) or for the following
c month (otherwise). These indices are for the mid of the month. The
c indices for the given day are obtained by linear interpolation and
c are stored in rz(3) and ig(3).
c
           integer      yr, mm, day, iflag, iyst, iyend,iymst
           integer      imst,iymend
           real         ionoindx(999),indrz(992)
           real         ig(3),rz(3)
           save         ionoindx,indrz,iflag,iyst,iymst,
     &                  iymend,imst
c
c Rz12 and IG are determined from the file IG_RZ.DDD which has the
c following structure: 
c day, month, year of the last update of this file,
c start month, start year, end month, end year,
c the 12 IG indices (13-months running mean) for the first year, 
c the 12 IG indices for the second year and so on until the end year,
c the 12 Rz indices (13-months running mean) for the first year,
c the 12 Rz indices for the second year and so on until the end year.
c The inteporlation procedure also requires the IG and Rz values for
c the month preceeding the start month and the IG and Rz values for the
c month following the end month. These values are also included in IG_RZ.
c 
c A negative Rz index means that the given index is the 13-months-running
c mean of the solar radio flux (F10.7). The close correlation between (Rz)12
c and (F10.7)12 is used to derive the (Rz)12 indices.
c
c An IG index of -111 indicates that no IG values are available for the
c time period. In this case a correlation function between (IG)12 and (Rz)12 
c is used to obtain (IG)12.
c
c The computation of the 13-month-running mean for month M requires the indices
c for the six months preceeding M and the six months following M (month: M-6, 
c ..., M+6). To calculate the current running mean one therefore requires
c predictions of the indix for the next six months. Starting from six months
c before the UPDATE DATE (listed at the top of the file) and onward the
c indices are therefore based on indices predictions.
c
      IRZ=19  !... PGR unit changed by P. Richards in May 2000
      if(iflag.eq.0) then
         !open(unit=IRZ,file='IG_RZ.DDD',status='old')
c-web- special for web version
c-web-        open(unit=12,file=
c-web-     *'/usr/local/etc/httpd/cgi-bin/models/IRI/IG_RZ.DDD',
c-web-     *status='old')
c Read the update date, the start date and the end date (mm,yyyy), and
c get number of data points to read.
          read(IRZ,*) iupd,iupm,iupy
          read(IRZ,*) imst,iyst, imend, iyend
          iymst=iyst*100+imst
          iymend=iyend*100+imend
c inum_vals= 12-imst+1+(iyend-iyst-1)*12 +imend + 2
c            1st year \ full years       \last y\ before & after
          inum_vals= 3-imst+(iyend-iyst)*12 +imend
c Read all the ionoindx and indrz values
          read(IRZ,*) (ionoindx(i),i=1,inum_vals)
          read(IRZ,*) (indrz(i),i=1,inum_vals)
          do 1 jj=1,inum_vals
                rrr=indrz(jj)
                if(rrr.lt.0.0) then
                        covr=abs(rrr)
                        rrr=33.52*sqrt(covr+85.12)-408.99
                        if(rrr.lt.0.0) rrr=0.0
                        indrz(jj)=rrr
                        endif
                if(ionoindx(jj).gt.-90.) goto 1
                  zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
                  if(zi.gt.274.0) zi=274.0
                  ionoindx(jj)=zi
1               continue
          CALL PC_CLOSE(IRZ)
          iflag = 1
        endif

        iytmp=yr*100+mm
        if (iytmp .lt. iymst .or. iytmp .gt. iymend) then
               write(6,8000) iytmp,iymst,iymend
 8000          format(1x,I10,'** OUT OF RANGE **'/,5x,
     &  'The file IG_RZ.DDD which contains the indices Rz12',
     &  ' and IG12'/5x,'currently only covers the time period',
     &  ' (yymm) : ',I6,'-',I6)
               nmonth=-1
               return
               endif

c       num=12-imst+1+(yr-iyst-1)*12+mm+1
        num=2-imst+(yr-iyst)*12+mm

        rz(1)=indrz(num)
        ig(1)=ionoindx(num)
        midm=15
        if(mm.eq.2) midm=14
        call MODA(0,yr,mm,midm,idd1,nrdaym)
        if(day.lt.midm) goto 1926
                imm2=mm+1
                if(imm2.gt.12) then
                        imm2=1
                        iyy2=yr+1
                        idd2=380
c               if((yr/4*4.eq.yr).and.(yr/100*100.ne.yr)) idd2=381
                        if(yr/4*4.eq.yr) idd2=381

                        if(yr/4*4.eq.yr) idd2=381
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num+1)
                ig(2)=ionoindx(num+1)
                rsn=(idn-idd1)*1./(idd2-idd1)
                rz(3)=rz(1)+(rz(2)-rz(1))*rsn
                ig(3)=ig(1)+(ig(2)-ig(1))*rsn
                goto 1927
1926            imm2=mm-1
                if(imm2.lt.1) then
                        imm2=12
                        idd2=-16
                        iyy2=yr-1
                else
                        iyy2=yr
                        midm=15
                        if(imm2.eq.2) midm=14
                        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
                endif
                rz(2)=indrz(num-1)
                ig(2)=ionoindx(num-1)
                rsn=(idn-idd2)*1./(idd1-idd2)
                rz(3)=rz(2)+(rz(1)-rz(2))*rsn
                ig(3)=ig(2)+(ig(1)-ig(2))*rsn

1927    nmonth=imm2
            return
            end
c
c
        subroutine LSTID(FI,ICEZ,R,AE,TM,SAX,SUX,TS70,DF0F2,DHF2)
C*****************************************************************

C   COMPUTER PROGRAM FOR UPDATING FOF2 AND HMF2 FOR EFFECTS OF 
C   THE LARGE SCALE SUBSTORM.
C
C   P.V.Kishcha, V.M.Shashunkina, E.E.Goncharova, Modelling of the
C   ionospheric effects of isolated and consecutive substorms on 
C   the basis of routine magnetic data, Geomagn. and Aeronomy v.32,
C   N.3, 172-175, 1992.
C
C   P.V.Kishcha et al. Updating the IRI ionospheric model for    
C   effects of substorms, Adv. Space Res.(in press) 1992.
C
C   Address: Dr. Pavel V. Kishcha, 
C            Institute of Terrestrial Magnetism,Ionosphere and Radio
C            Wave Propagation, Russian Academy of Sciences, 
C            142092, Troitsk, Moscow Region, Russia
C
C***       INPUT PARAMETERS:
C       FI ------ GEOMAGNETIC LATITUDE,
C       ICEZ ---- INDEX OF SEASON(1-WINTER AND EQUINOX,2-SUMMER),
C       R ------- ZURICH SUNSPOT NUMBER,
C       AE ------ MAXIMUM AE-INDEX REACHED DURING SUBSTORM,
C       TM ------ LOCAL TIME,
C       SAX,SUX - TIME OF SUNSET AND SUNRISE,
C       TS70 ---- ONSET TIME (LT) OF SUBSTORMS ONSET 
C                        STARTING ON FI=70 DEGR.
C***      OUTPUT PARAMETERS:
C       DF0F2,DHF2- CORRECTIONS TO foF2 AND hmF2 FROM IRI OR 
C                   OBSERVATIONAL MEDIAN  OF THOSE VALUES.
C*****************************************************************


        INTEGER ICEZ
        REAL A(7,2,3,2),B(7,2,3,2),C(7,2,3,2),D(7,2,3,2),A1(7,2,2),
     *B1(7,2,2),Y1(84),Y2(84),Y3(84),Y4(84),Y5(28),Y6(28)
        DATA Y1/
     *150.,250.,207.8,140.7,158.3,87.2,158.,
     *150.,250.,207.8,140.7,158.3,87.2,158.,
     *115.,115.,183.5,144.2,161.4,151.9,272.4,
     *115.,115.,183.5,144.2,161.4,151.9,272.4,
     *64.,320.,170.6,122.3,139.,79.6,180.6,
     *64.,320.,170.6,122.3,139.,79.6,180.6,
     *72.,84.,381.9,20.1,75.1,151.2,349.5,
     *120.,252.,311.2,241.,187.4,230.1,168.7,
     *245.,220.,294.7,181.2,135.5,237.7,322.,
     *170.,110.,150.2,136.3,137.4,177.,114.,
     *170.,314.,337.8,155.5,157.4,196.7,161.8,
     *100.,177.,159.8,165.6,137.5,132.2,94.3/
        DATA Y2/
     *2.5,2.0,1.57,2.02,2.12,1.46,2.46,
     *2.5,2.0,1.57,2.02,2.12,1.46,2.46,
     *2.3,1.6,1.68,1.65,2.09,2.25,2.82,
     *2.3,1.6,1.68,1.65,2.09,2.25,2.82,
     *0.8,2.0,1.41,1.57,1.51,1.46,2.2,
     *0.8,2.0,1.41,1.57,1.51,1.46,2.2,
     *3.7,1.8,3.21,3.31,2.61,2.82,2.34,
     *2.8,3.2,3.32,3.33,2.96,3.43,2.44,
     *3.5,2.8,2.37,2.79,2.26,3.4,2.28,
     *3.9,2.,2.22,1.98,2.33,3.07,1.56,
     *3.7,3.,3.3,2.99,3.57,2.98,3.02,
     *2.6,2.8,1.66,2.04,1.91,1.49,0.43/
        DATA Y3/
     *-1.8,-1.9,-1.42,-1.51,-1.53,-1.05,-1.66,
     *-1.8,-1.9,-1.42,-1.51,-1.53,-1.05,-1.66,
     *-1.5,-1.3,-1.46,-1.39,-1.53,-1.59,-1.9,
     *-1.5,-1.3,-1.46,-1.39,-1.53,-1.59,-1.9,
     *-0.7,-2.,-1.41,-1.09,-1.22,-0.84,-1.32,
     *-0.7,-2.,-1.41,-1.09,-1.22,-0.84,-1.32,
     *-1.7,-1.0,-2.08,-1.80,-1.35,-1.55,-1.79,
     *-1.5,-2.,-2.08,-2.16,-1.86,-2.19,-1.70,
     *-2.2,-1.7,-1.57,-1.62,-1.19,-1.89,-1.47,
     *-1.9,-1.5,-1.26,-1.23,-1.52,-1.89,-1.02,
     *-1.7,-1.7,-1.76,-1.43,-1.66,-1.54,-1.24,
     *-1.1,-1.5,-1.09,-1.23,-1.11,-1.14,-0.4/
        DATA Y4/
     *-2.,-5.,-5.,0.,0.,0.,2.,
     *-2.,-5.,-5.,0.,0.,0.,2.,
     *-5.,-5.,6.,0.,1.,5.,2.,
     *-5.,-5.,6.,0.,1.,5.,2.,
     *0.,-7.,-3.,-6.,2.,2.,3.,
     *0.,-7.,-3.,-6.,2.,2.,3.,
     *-5.,-1.,-11.,-6.,0.,-5.,-6.,
     *-5.,-10.,1.,4.,-6.,-2.,1.,
     *2.,-13.,-10.,0.,-8.,10.,-16.,
     *0.,-3.,-7.,-2.,-2.,4.,2.,
     *-11.,-12.,-13.,0.,0.,7.,0.,
     *-8.,6.,-1.,-5.,-7.,4.,-4./
        DATA Y5/
     *0.,0.,-0.1,-0.19,-0.19,-0.25,-0.06,
     *0.,0.,-0.31,-0.28,-0.27,-0.06,0.02,
     *0.,0.,0.18,-0.07,-0.2,-0.1,0.3,
     *0.,0.,-0.24,-0.5,-0.4,-0.27,-0.48/
        DATA Y6/
     *0.,0.,-0.00035,-0.00028,-0.00033,-0.00023,-0.0007,
     *0.,0.,-0.0003,-0.00025,-0.0003,-0.0006,-0.00073,
     *0.,0.,-0.0011,-0.0006,-0.0003,-0.0005,-0.0015,
     *0.,0.,-0.0008,-0.003,-0.0002,-0.0005,-0.0003/

        INN=0
        IF(TS70.GT.12. .AND. TM.LT.SAX)INN=1
        IF(FI.LT.0.)FI=ABS(FI)

        N=0
        DO 001 M=1,2
        DO 001 K=1,3
        DO 001 J=1,2
        DO 001 I=1,7
        N=N+1
        A(I,J,K,M)=Y1(N)
        B(I,J,K,M)=Y2(N)
        C(I,J,K,M)=Y3(N)
 001    D(I,J,K,M)=Y4(N)
        N1=0
        DO 002 M=1,2
        DO 002 J=1,2
        DO 002 I=1,7
        N1=N1+1
        A1(I,J,M)=Y5(N1)
 002    B1(I,J,M)=Y6(N1)
        IF(FI.GT.65..OR.AE.LT.500.)THEN
        WRITE(*,*)'LSTID are for AE>500. and ABS(FI)<65.'
        GOTO 004
        ENDIF
        TS=TS70+(-1.5571*FI+109.)/60.
        IF(TS.LT.SUX.AND.TS.GT.SAX)THEN
        WRITE(*,*)' LSTID are only at night'
        GOTO 004
        ENDIF
        IF(INN.EQ.1)TM=TM+24.
        
        IF(TS.GE.TM.OR.TS.LT.TM-5.)THEN
C        WRITE(*,*)'LSTID are onli if  TM-5.<TS<TM ;Here TS=',TS,'TM=',TM     
        GOTO 004
        ENDIF
        DO 007 I=1,7
        IF(FI.GE.-5.+10.*(I-1) .AND. FI.LT.5.+10.*(I-1))GOTO 008
 007    CONTINUE
 008    J=ICEZ
        IF(AE.GE.500. .AND. AE.LE.755.)K=1
        IF(AE.GT.755. .AND. AE.LT.1000.)K=2
        IF(AE.GE.1000.)K=3
        M=-1
        IF(R.LE.20.)M=1
        IF(R.GE.120.)M=2
        T=TM-TS
        IF(M.LT.0)GOTO 003
C        WRITE(*,*)'A1=',A1(I,J,M),' B1=',B1(I,J,M)
C        WRITE(*,*)'A=',A(I,J,K,M),' B=',B(I,J,K,M),' C=',C(I,J,K,M),
C     *'D=',D(I,J,K,M)
        DF0F2=A1(I,J,M)+B1(I,J,M)*AE
        DHF2=A(I,J,K,M)*(T**B(I,J,K,M))*EXP(C(I,J,K,M)*T)+D(I,J,K,M)
        GOTO 005
 003    DF1=A1(I,J,1)+B1(I,J,1)*AE
        DF2=A1(I,J,2)+B1(I,J,2)*AE
        DF0F2=DF1+(DF2-DF1)*(R-20.)/100.
        DH1=A(I,J,K,1)*(T**B(I,J,K,1))*EXP(C(I,J,K,1)*T)+D(I,J,K,1)
        DH2=A(I,J,K,2)*(T**B(I,J,K,2))*EXP(C(I,J,K,2)*T)+D(I,J,K,2)
        DHF2=DH1+(DH2-DH1)*(R-20.)/100.
        GOTO 005
 004    DHF2=0.
        DF0F2=0.
 005    CONTINUE
        IF(INN.EQ.1)TM=TM-24.
        RETURN
        END
      SUBROUTINE CIRA86(IDAY,SEC,GLAT,GLONG,STL,F107A,TINF,TLB,SIGMA)
C*******************************************************************
C   Calculates neutral temperature parameters for IRI using the
C   MSIS-86/CIRA 1986 Neutral Thermosphere Model. The subroutines
C   GTS5, GLOBE5 and GLOBL5 developed by A.E. Hedin (2/26/87) were
C   modified for use in IRI --------- D. Bilitza -------- March 1991
C
C   CHANGES:
C     11/09/99 always calculated Legendre; 'if glat' and 'if stl' taken out;
C     11/09/99 use UMR, dumr and humr from COMMON 
C
C     INPUT:
C        IDAY - DAY OF YEAR 
C        SEC - UT(SEC)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C
C     OUTPUT: 
C        TINF - EXOSPHERIC TEMPERATURE (K)
C        TLB - TEMPERATURE AT LOWER BOUNDARY (K)
C        SIGMA - SHAPE PARAMETER FOR TEMPERATURE PROFILE
C
C **********************************************************************
      DIMENSION         PLG(9,4)
        common  /const/umr      /const1/hr,dr
      DATA      XL/1000./,TLL/1000./

c    data umr/1.74E-2/,hr/0.2618/,dr/1.74e-2
cDR,DR2/1.72142E-2,0.0344284/,
cSR/7.2722E-5/,
c,HR/.2618/
c,DGTR/1.74533E-2/
c       dr = hr * 24. / 365.
c
        dr2 = dr * 2.
        sr = hr / 3600.

C
C CALCULATE LEGENDRE POLYNOMIALS
C
      IF(XL.EQ.GLAT)   GO TO 15
      C = SIN(GLAT*umr)
      S = COS(GLAT*umr)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5*(3.*C2 -1.)
      PLG(4,1) = 0.5*(5.*C*C2-3.*C)
      PLG(5,1) = (35.*C4 - 30.*C2 + 3.)/8.
      PLG(6,1) = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.
      PLG(2,2) = S
      PLG(3,2) = 3.*C*S
      PLG(4,2) = 1.5*(5.*C2-1.)*S
      PLG(5,2) = 2.5*(7.*C2*C-3.*C)*S
      PLG(6,2) = 1.875*(21.*C4 - 14.*C2 +1.)*S
      PLG(7,2) = (11.*C*PLG(6,2)-6.*PLG(5,2))/5.
      PLG(3,3) = 3.*S2
      PLG(4,3) = 15.*S2*C
      PLG(5,3) = 7.5*(7.*C2 -1.)*S2
      PLG(6,3) = 3.*C*PLG(5,3)-2.*PLG(4,3)
      PLG(4,4) = 15.*S2*S
      PLG(5,4) = 105.*S2*S*C 
      PLG(6,4)=(9.*C*PLG(5,4)-7.*PLG(4,4))/2.
      PLG(7,4)=(11.*C*PLG(6,4)-8.*PLG(5,4))/3.
      XL=GLAT
   15 CONTINUE
      IF(TLL.EQ.STL)   GO TO 16
      STLOC = SIN(HR*STL)
      CTLOC = COS(HR*STL)
      S2TLOC = SIN(2.*HR*STL)
      C2TLOC = COS(2.*HR*STL)
      S3TLOC = SIN(3.*HR*STL)
      C3TLOC = COS(3.*HR*STL)
      TLL = STL
   16 CONTINUE
C
      DFA=F107A-150.
C
C EXOSPHERIC TEMPERATURE
C
C         F10.7 EFFECT
      T1 =  ( 3.11701E-3 - 0.64111E-5 * DFA ) * DFA
        F1 = 1. + 0.426385E-2 * DFA 
        F2 = 1. + 0.511819E-2 * DFA
        F3 = 1. + 0.292246E-2 * DFA
C        TIME INDEPENDENT
      T2 = 0.385528E-1 * PLG(3,1) + 0.303445E-2 * PLG(5,1)
C        SYMMETRICAL ANNUAL AND SEMIANNUAL
        CD14 = COS( DR  * (IDAY+8.45398) )
        CD18 = COS( DR2 * (IDAY-125.818) )
        CD32 = COS( DR  * (IDAY-30.0150) )
        CD39 = COS( DR2 * (IDAY-2.75905) )
      T3 = 0.805486E-2 * CD32 + 0.14237E-1 * CD18
C        ASYMMETRICAL ANNUAL AND SEMIANNUAL
      T5 = F1 * (-0.127371 * PLG(2,1) - 0.302449E-1 * PLG(4,1) ) * CD14
     &       - 0.192645E-1 * PLG(2,1) * CD39
C        DIURNAL
      T71 =  0.123512E-1 * PLG(3,2) * CD14
      T72 = -0.526277E-2 * PLG(3,2) * CD14
      T7 = ( -0.105531 *PLG(2,2) - 0.607134E-2 *PLG(4,2) + T71 ) *CTLOC
     4   + ( -0.115622 *PLG(2,2) + 0.202240E-2 *PLG(4,2) + T72 ) *STLOC
C        SEMIDIURNAL
      T81 = 0.386578E-2 * PLG(4,3) * CD14
      T82 = 0.389146E-2 * PLG(4,3) * CD14
      T8= (-0.516278E-3 *PLG(3,3) - 0.117388E-2 *PLG(5,3) +T81)*C2TLOC
     3   +( 0.990156E-2 *PLG(3,3) - 0.354589E-3 *PLG(5,3) +T82)*S2TLOC
C        TERDIURNAL
      Z1 =  PLG(5,4) * CD14
      Z2 =  PLG(7,4) * CD14
      T14=(0.147284E-2*PLG(4,4)-0.173933E-3*Z1+0.365016E-4*Z2)*S3TLOC
     2   +(0.341345E-3*PLG(4,4)-0.153218E-3*Z1+0.115102E-3*Z2)*C3TLOC 
      T7814 = F2 * ( T7 + T8 + T14 )
C        LONGITUDINAL
      T11= F3 * (( 0.562606E-2 * PLG(3,2) + 0.594053E-2 * PLG(5,2) + 
     $       0.109358E-2 * PLG(7,2) - 0.301801E-2 * PLG(2,2) - 
     $       0.423564E-2 * PLG(4,2) - 0.248289E-2 * PLG(6,2) + 
     $      (0.189689E-2 * PLG(2,2) + 0.415654E-2 * PLG(4,2)) * CD14
     $     ) * COS(umr*GLONG) +
     $     ( -0.11654E-1 * PLG(3,2) - 0.449173E-2 * PLG(5,2) - 
     $       0.353189E-3 * PLG(7,2) + 0.919286E-3 * PLG(2,2) + 
     $       0.216372E-2 * PLG(4,2) + 0.863968E-3 * PLG(6,2) +
     $      (0.118068E-1 * PLG(2,2) + 0.331190E-2 * PLG(4,2)) * CD14
     $     ) * SIN(umr*GLONG) )
C        UT AND MIXED UT,LONGITUDE
      T12 = ( 1. - 0.565411 * PLG(2,1) ) * COS( SR*(SEC-31137.0) ) *
     $ (-0.13341E-1*PLG(2,1)-0.243409E-1*PLG(4,1)-0.135688E-1*PLG(6,1))
     $ +    ( 0.845583E-3 * PLG(4,3) + 0.538706E-3  * PLG(6,3) ) *
     $      COS( SR * (SEC-247.956) + 2.*umr*GLONG )
C  Exospheric temperature TINF/K  [Eq. A7]
      TINF = 1041.3 * ( 1. + T1+T2+T3+T5+T7814+T11+T12 ) * 0.99604
C
C TEMPERATURE DERIVATIVE AT LOWER BOUNDARY
C
C         F10.7 EFFECT
      T1 =  0.252317E-2 * DFA 
C        TIME INDEPENDENT
      T2 = -0.467542E-1 * PLG(3,1) + 0.12026 * PLG(5,1) 
C        ASYMMETRICAL ANNUAL
        CD14 = COS( DR  * (IDAY+8.45398) )
      T5 = -0.13324 * PLG(2,1)  * CD14
C        SEMIDIURNAL
      ZZ = PLG(4,3) * CD14
      T81 = -0.973404E-2 * ZZ
      T82 = -0.718482E-3 * ZZ
      T8 =(0.191357E-1 *PLG(3,3) + 0.787683E-2 *PLG(5,3) + T81) *C2TLOC
     3  + (0.125429E-2 *PLG(3,3) - 0.233698E-2 *PLG(5,3) + T82) *S2TLOC
C  dTn/dh at lower boundary  [Eq. A6]
      G0 = 0.166728E2 * ( 1. + T1+T2+T5+T8 ) * 0.951363
C
C NEUTRAL TEMPERATURE AT LOWER BOUNDARY 120KM
C
        CD9  = COS( DR2 * (IDAY-89.3820) )
        CD11 = COS( DR  * (IDAY+8.45398) )
      T1 = 0.568478E-3 * DFA
      T4 = 0.107674E-1 * CD9
      T5 =-0.192414E-1 * PLG(2,1) * CD11
      T7 = -0.2002E-1 *PLG(2,2) *CTLOC - 0.195833E-2 *PLG(2,2) *STLOC
      T8 = (-0.938391E-2 * PLG(3,3) - 0.260147E-2 * PLG(5,3)
     $       + 0.511651E-4 * PLG(6,3) * CD11 ) * C2TLOC
     $   + ( 0.131480E-1 * PLG(3,3) - 0.808556E-3 * PLG(5,3)
     $       + 0.255717E-2 * PLG(6,3) * CD11 ) * S2TLOC
C  Tn at lower boundary 120km   [Eq. A8]
      TLB = 386.0 * ( 1. + T1+T4+T5+T7+T8 ) * 0.976619
C  Sigma      [Eq. A5]
      SIGMA = G0 / ( TINF - TLB )
      RETURN
      END
