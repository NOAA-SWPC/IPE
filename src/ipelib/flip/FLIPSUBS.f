C................................ FLIPSUBS.FOR :::::::::::::::::::::
C......... This subroutine is used to initialize the arrays in the FLIP
C......... model. This is necessary for the IBM machine
C::::::::::::::::::::::::::::::: FLIPIN :::::::::::::::::::::::::::::
      SUBROUTINE FLIPIN(S)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      INTEGER ION
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,EUV,UVFAC
      DIMENSION S(50000)
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/ND/ON(401),HN(401),N2N(401),O2N(401),PHION(401),TN(401)
      COMMON/SAV/NSAVE(2,401),TISAV(3,401),FY(2,401),UN(401),EHT(3,401)
      COMMON/STAW/STR(12,401),RCON(401)
      COMMON/CFACT/DHDU_FAC,COLFAC,CFAC(9),ISCH(9)

      IVERT=0
      DO 828 IIJJ=1,50000
 828  S(IIJJ)=0.0

      SL(1)=0.0
      JMAX=401

      DO 15 I = 1, 3
      DO 15 J = 1, JMAX
        EHT(I,J) = 0.0D0
        TI(I,J) = 0.0D0
 15     TISAV(I,J) = 0.0D0

      DO 20 I = 1, 6
      DO 20 J = 1, JMAX
          STR (I,J) = 0.0D0
          STR (I+6, J) = 0.0D0
 20   CONTINUE

      DO 30 I = 1, 2
      DO 30 J = 1, JMAX
         N(I,J) = 0.0D0
         N(I+1,J) = 0.0D0
         NSAVE(I,J) = 0.0D0
  30     FY(I,J) = 0.0D0

      DO 35 J = 1, JMAX
         GR(J) = 0.0D0
         GL(J) = 0.0D0
         BM(J) = 0.0D0
         SL(J) = 0.0D0
         ON(J) = 0.0D0
         HN(J) = 0.0D0
         N2N(J) = 0.0D0
         O2N(J) = 0.0D0
         PHION(J) = 0.0D0
         TN(J) = 0.0D0
         Z(J) = 0.0D0
         UN(J) = 0.0D0
         SZA(J) = 0.0D0
 35      RCON(J) = 0.0D0

      DO 45 I = 1,20
 45     F(I) = 0.0D0

      DO I=1,9
         CFAC(I)=0.0
         ISCH(I)=0
      ENDDO
      DO 65 I=1,ISPEC
      DO 65 J=1,JMAX
         XIONN(I,J)=0.0D0
         XIONV(I,J)=0.0D0
 65   CONTINUE
      RETURN
      END
C:::::::::::::::::::::: DIAGPR :::::::::::::::::::::::::::
C... This routine prints FLIP diagnostics
      SUBROUTINE DIAGPR(ISOL,JB,JK,JC,JQ,JNQ,JFQ,DT,HMF2N,HMF2S
     > ,NMF2N,NMF2S,SATN,SATF,CHIN,CHIF,UNN,UNC,VTOT,PCO,VPERPN,COSDIP)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,D,T,CHI,SATN,SATF,SAT
     > ,GLONN,GLONS,CHIN,CHIF,GLATN,GLATS,UTHRS,RE,PLAT,PLON,PI,ANGVEL
     > ,GLON_EQ,GLAT_EQ,SZA_EQ,SAT_EQ
      DIMENSION D(19),T(2),NT(4)
      COMMON/PALT/ZLO,ZHI,JWR,JTI
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/ND/ON(401),HN(401),N2N(401),O2N(401),PHION(401),TN(401)
      COMMON/SAV/NSAVE(2,401),TISAV(3,401),FY(2,401),UN(401),EHT(3,401)
      COMMON/STAW/STR(12,401),RCON(401)
      COMMON/GPAR/RE,PLAT,PLON,PI,ANGVEL
      COMMON/CFACT/DHDU_FAC,COLFAC,CFAC(9),ISCH(9)
      DATA JJJ/0/
      !......... printing densities and temperatures in UNIT=8 ....
 453  UTHRS=SEC/3600.0

      !.... determine the total flux tube contents (NT)
      DO 5 I=1,4
 5    NT(I)=0.0
      DO 29 J=JNQ,JMAX-JNQ
         NT(1)=NT(1)+BM(JB)*N(1,J)*(SL(J)-SL(J-1))/BM(J)
         NT(2)=NT(2)+BM(JB)*N(2,J)*(SL(J)-SL(J-1))/BM(J)
         NT(4)=NT(4)+BM(JB)*XIONN(3,J)*(SL(J)-SL(J-1))/BM(J)
 29   CONTINUE

      !.. determine the indices for printing 
      JZLON=JQ
      J300=JQ
      DO JS=1,JQ
         J=JQ-JS+1                      !.. coming down in altitude
         IF(Z(J).GE.ZLO) JZLON=J        !.. index for specified altitude
         IF(Z(J).GT.300.0) J300=J       !.. index of 300 km altitude
         IF(Z(J).GT.180.0) J180K=J      !.. index of 180 km altitude
      ENDDO

      !.. Make sure closest altitude taken for ZLO alt
      IF(ABS(Z(JZLON-1)-ZLO).LT.ABS(ZLO-Z(JZLON))) JZLON=JZLON-1
      JZLOC=JMAX+1-JZLON        !.. conjugate point to ZLO
      J180C=JMAX+1-J180K        !.. conjugate 180 km altitude


      !.. determine the index (JK) of maximum [e] in NORTHern hemisphere.
      !.. [e] = [O+] + [min+].
      JK=JQ
      DO JS=2,JQ
         J=JQ-JS+1   !.. come down in altitude
         !.. Check for peak density altitude
         IF(Z(J).GT.888) JK=J    
         IF(N(1,J)+N(3,J).GT.N(1,J-1)+N(3,J-1).AND. 
     >      N(1,J)+N(3,J).GT.N(1,J+1)+N(3,J+1)) THEN
               !.. If there are multiple peaks take the largest one
               IF(N(1,J)+N(3,J).GT.N(1,JK)+N(3,JK)) JK=J
         ENDIF
         IF(Z(J).LT.180.0) GO TO 130   !.. hmF2 not found
      ENDDO
 130  CONTINUE

      !... Determine HmF2 and NmF2 by fitting a parabola: If no [e] max above
      !... 180 km just use 180 km.
      IF(ISOL.NE.0) THEN
        IF(N(1,JK)+N(3,JK).GT.N(1,JK-1)+N(3,JK-1). AND. 
     >    N(1,JK)+N(3,JK).GT.N(1,JK+1)+N(3,JK+1) ) THEN
            CALL F2PEAK(N(1,JK-1)+N(3,JK-1),N(1,JK)+N(3,JK)
     >        ,N(1,JK+1)+N(3,JK+1),Z(JK-1),Z(JK),Z(JK+1),NMF2N,HMF2N)
        !.. maximum is at equator.
        ELSEIF(JK.EQ.JQ) THEN
          HMF2N=Z(JQ)
          NMF2N=N(1,JQ)+N(3,JQ)
        !.. No maximum above 180 km
        ELSE
          WRITE(6,970)
          HMF2N=Z(J180K)
          NMF2N=N(1,J180K)+N(3,J180K)
        ENDIF
      ENDIF
 970  FORMAT( ' *** WARNING, North hmF2 not found, using ~180 km')

      !.. determine the index (JC) of maximum [e] in SOUTHern hemisphere.
      !.. [e] = [O+] + [min+].
      JC=JQ
      DO J=JQ+1,JMAX-1
         !.. Check for peak density altitude
         IF(Z(J).GT.888) JC=J    
         IF(N(1,J)+N(3,J).GT.N(1,J-1)+N(3,J-1).AND. 
     >      N(1,J)+N(3,J).GT.N(1,J+1)+N(3,J+1)) THEN
               !.. If there are multiple peaks take the largest one
               IF(N(1,J)+N(3,J).GT.N(1,JC)+N(3,JC)) JC=J
         ENDIF
         IF(Z(J).LT.180.0) GO TO 230   !.. hmF2 not found
      ENDDO
 230  CONTINUE

      !... Determine HmF2 and NmF2 by fitting a parabola: If no [e] max above
      !... 180 km just use 180 km.
      IF(ISOL.NE.0) THEN
        IF(N(1,JC)+N(3,JC).GT.N(1,JC-1)+N(3,JC-1). AND. 
     >    N(1,JC)+N(3,JC).GT.N(1,JC+1)+N(3,JC+1) ) THEN
            CALL F2PEAK(N(1,JC-1)+N(3,JC-1),N(1,JC)+N(3,JC)
     >        ,N(1,JC+1)+N(3,JC+1),Z(JC-1),Z(JC),Z(JC+1),NMF2S,HMF2S)
        !.. maximum is at equator.
        ELSEIF(JC.EQ.JQ) THEN
          HMF2N=Z(JQ)
          NMF2N=N(1,JQ)+N(3,JQ)
        !.. No maximum above 180 km
        ELSE
          WRITE(6,971)
          HMF2S=Z(J180C)
          NMF2S=N(1,J180C)+N(3,J180C)
        ENDIF
      ENDIF
 971  FORMAT( ' *** WARNING, South hmF2 not found, using ~180 km')

      !... if this is first time through, no F2 peak
      IF(ISOL.EQ.0) JK=J180K
      IF(ISOL.EQ.0) JC=J180K
      JFQ=JMAX+1-JNQ

      !.. Determine the fluxes for printing
      FYNOP=N(1,JNQ)*XIONV(1,JNQ)*BM(JB)/BM(JNQ)
      FYNHP=N(2,JNQ)*XIONV(2,JNQ)*BM(JB)/BM(JNQ)
      FYNHEP=XIONN(3,JNQ)*XIONV(3,JNQ)*BM(JB)/BM(JNQ)
      FYNNP=XIONN(4,JNQ)*XIONV(4,JNQ)*BM(JB)/BM(JNQ)
      FYFOP=N(1,JFQ)*XIONV(1,JFQ)*BM(JB)/BM(JFQ)
      FYFHP=N(2,JFQ)*XIONV(2,JFQ)*BM(JB)/BM(JFQ)
      FYFHEP=XIONN(3,JFQ)*XIONV(3,JFQ)*BM(JB)/BM(JFQ)
      FYFNP=XIONN(4,JFQ)*XIONV(4,JFQ)*BM(JB)/BM(JFQ)

      !.. write hemisphere header information
      JJJ=JJJ+1
      IF(JJJ.EQ.1) WRITE(3,351) NINT(ZLO),NINT(ZHI)
      IF(JJJ.EQ.1) WRITE(9,352) NINT(ZLO),NINT(ZHI)
      IF(JJJ.EQ.1) WRITE(4,353)
 351  FORMAT(/'  Important parameters in the NORTHern hemisphere:-'
     > ,1X,' Ion and Neutrals near alt _Z = ',I7,' km; fluxes (Phi)'
     > ,1X,'near',I7,' km.',/2X,'Note - For low L values these print'
     > ,1X,'altitudes may be set to the field line apex altitude.')
 352  FORMAT(/5X,'Important parameters in the SOUTHern hemisphere'
     > ,1X,' Ion and Neutrals near alt _Z = ',I7,' km; fluxes (Phi)'
     > ,1X,'near',I7,' km.',/2X,'Note - For low L values these print'
     > ,1X,'altitudes may be set to the field line apex altitude.')
 353  FORMAT(/5X,'Important parameters in the EQUATorial plane. If ExB'
     > ,1X,'drift is on, the equatorial apex altitude will vary.')
       
      !.. write warning about first 12 hours
      IF(JJJ.EQ.1.AND.ISOL.EQ.0) WRITE(3,951)
      IF(JJJ.EQ.1.AND.ISOL.EQ.0) WRITE(4,951)
      IF(JJJ.EQ.1.AND.ISOL.EQ.0) WRITE(9,951)
 951  FORMAT(/1X,9('*'),'Do not use first 10-12 hours - too close to'
     > ,1X,'initial conditions',9('*'))

      !.. write column header information
      IF(JJJ.EQ.1) WRITE(3,451)
      IF(JJJ.EQ.1) WRITE(4,534)
      IF(JJJ.EQ.1) WRITE(9,451)
 451  FORMAT(/5X,'UT    LT  SZA  Te  HmF2  Wind   NmF2     min+'
     >  ,5X,'PhiH+     PhiO+    PhiHe+    OX_Z     O2_Z     N2_Z'
     >  ,5X,'H_Z     He_Z     Tn   PhiN+   TI_Z   TE_Z      O+_Z'
     >  ,5X,'H+_Z      min+_Z     RON2   TBrace_Z     TBF_Z')

 534  FORMAT(/4X,'UT',4X,'LT',4X,'Ti',4X,'Te',3X,'[O+]',4X
     > ,'[H+]',4X,'[He+]',3X,'PhiH+',5X,'PhiO+',5X,'PhiHe+',4X
     > ,'NT(O+)',3X,'NT(H+)',3X,'NT(He+)',5X,'[N+]',4X,'PhiN+')

      !-------------- Northern Hemisphere neutrals ----------------
      CALL AMBS(1,ZLO,GL(JZLON),D,T,CHI,SAT,GLONN,GLATN)

      NE=N(1,JZLON)+N(2,JZLON)+N(3,JZLON)
      CALL BRACE(Z(JZLON),NE,TI(3,JZLON),TB,TEOTB)

      !... Print north params but get neutral densities first
       WRITE(3,355) UTHRS,SATN,NINT(57.3*CHIN),NINT(TI(3,JK)),
     > HMF2N,NINT(UNN/100.),NMF2N,N(3,JK),FYNHP,FYNOP,FYNHEP
     > ,D(2),D(4),D(3),D(7),D(1),NINT(T(2)),FYNNP,TI(2,JZLON)
     > ,TI(3,JZLON),N(1,JZLON),N(2,JZLON),N(3,JZLON),D(2)/D(3),
     >  TB,TEOTB

      !-------------- Southern Hemisphere neutrals ----------------
      CALL AMBS(JMAX,ZLO,GL(JZLOC),D,T,CHI,SAT,GLONS,GLATS)

      NE=N(1,JZLOC)+N(2,JZLOC)+N(3,JZLOC)
      CALL BRACE(Z(JZLOC),NE,TI(3,JZLOC),TB,TEOTB)

      WRITE(9,355) UTHRS,SATF,NINT(57.3*CHIF),NINT(TI(3,JC))
     > ,HMF2S,NINT(UNC/100.),NMF2S,N(3,JC),-FYFHP,-FYFOP,-FYFHEP
     > ,D(2),D(4),D(3),D(7),D(1),NINT(T(2)),-FYFNP,TI(2,JC)
     > ,TI(3,JZLOC),N(1,JZLOC),N(2,JZLOC),N(3,JZLOC),D(2)/D(3),
     >  TB,TEOTB

      !. Get LT at equator
      CALL AMBS(1,Z(JQ),GL(JQ),D,T,SZA_EQ,SAT_EQ,GLON_EQ,GLAT_EQ)

      WRITE(4,357) UTHRS,SAT_EQ,NINT(TI(2,JQ))
     > ,NINT(TI(3,JQ)),N(1,JQ),N(2,JQ),XIONN(3,JQ)
     > ,N(2,JQ)*XIONV(2,JQ)*BM(JB)/BM(JQ),
     >  N(1,JQ)*XIONV(1,JQ)*BM(JB)/BM(JQ)
     > ,XIONN(3,JQ)*XIONV(3,JQ)*BM(JB)/BM(JQ),NT(1),NT(2),NT(4)
     > ,XIONN(4,JQ),XIONN(4,JQ)*XIONV(4,JQ)*BM(JB)/BM(JQ),PCO

      !... Write parameters that change for E x B drifts
      IF(JJJ.EQ.1) WRITE(66,354) NINT(Z(JNQ))
 354  FORMAT(/5X,' *** Parameters that change during ExB drift ***'
     > ,3X,' ion fluxes at ',I7,' km'
     > ,/4X,'UT   LT_EQ L-val  Vperp   Ilat   latN  longN   latS  longS'
     > ,1X,'tube_vol alt1    netHflx   netOflx  H+_tot   O+_tot'
     > ,3X,'H+_eq    O+_eq    WindN  WindS')

      IF(JJJ.EQ.1) NETHFLX=0.0 
      IF(JJJ.EQ.1) NETOFLX=0.0 
      NETHFLX=NETHFLX+(FYNHP-FYFHP)*DT
      NETOFLX=NETOFLX+(FYNOP-FYFOP)*DT

      VPERPI=VPERPN/(PCO*SQRT(4*PCO-3))     !..E in the ionosphere

      !..evaluate the actual wind in north and south
      SINDIP=SIN(ACOS(COSDIP))
      EQW_N=(UNN-VPERPI/SINDIP)/100
      EQW_S=(UNC-VPERPI/SINDIP)/100

      IF(JJJ.NE.1) WRITE(66,358) UTHRS,SAT_EQ,PCO,VPERPI/100.0
     >  ,GL(J300)*57.29578,GLATN,GLONN,GLATS,GLONS,VTOT,NINT(Z(1))
     >  ,NETHFLX,NETOFLX,NT(2),NT(1),N(2,JQ),N(1,JQ),UNN/100.0,UNC/100.0

 355  FORMAT(F8.2,F6.2,I4,I5,F6.1,I5,1P,E9.2,1P,E9.2,1P,2E10.2,1P,E9.1
     >  ,5E9.2,I5,E9.1,0P,2F7.1,1P,11E10.2)
 356  FORMAT(F8.2,F6.2,I4,I5,I4,I5,1P,E9.2,1P,E9.2,1P,2E10.2,1P,E9.1
     >  ,5E8.1,I5,5E9.1)
 357  FORMAT(F8.2,F6.2,I6,I6,1P,3E8.1,1P,3E10.2,1P,5E9.1,1P,9E10.2)
 358  FORMAT(F8.2,F6.2,F6.3,6F7.1,1P,E9.2,I6,1P,2E10.2,4E9.2,0P,3F7.1)
      RETURN
      END
C::::::::::::::::::::::: MAGDEC :::::::::::::::::::::::::::
C...... This routine uses a table to calculate the magnetic declination
C...... (DECPT) for use with a neutral wind model. The table was constructed 
C...... using " The Earth's magnetic Field" by Robert Merrill and
C...... Michael McElhinny, Academic Press p 18. (UAH library # QC816.M47) 
C...... The declination at the specified location (degrees) is obtained
C...... by bilinear interpolation of the table values
C...... Written by Phil Richards, August 1995
C...... This code is not accurate near the magnetic poles
C...... GLAT, GLONG, and DECPT are in degrees
       SUBROUTINE MAGDEC(GLAT,GLONG,DECPT)
       IMPLICIT NONE
      !..--- NLONG=# longs, NLAT= # lats. I,J,K are array indices
       INTEGER NLONG,NLAT,I,J,K
       PARAMETER (NLONG=13,NLAT=11) 
       REAL GLAT,GLONG,DECPT
       REAL ALONG(NLONG),ALAT(NLAT),DECLS(13,11),RLAT,RLONG
      !... Fill up the longitude, latitude, and magnetic declination arrays
       DATA ALONG/0,30,60,90,120,150,180,210,240,270,300,330,360/
       DATA ALAT/-90,-80,-60,-40,-20,0,20,40,60,80,90/
      !..-- The declinations are loaded in latitude slices. The first NLONG
      !..-- values are for the first latitude in the ALAT array
       DATA DECLS / -20,-45,-67,-95,-125,170,110,80,60,40,20,0,-20
     > ,-20,-40,-63,-82,-120,140,82,65,52,35,16,-2,-20, -20,-35,-55,-70
     > ,-60,35,45,40,38,30,10,-7.5,-20, -27,-28,-42,-35,-5,13,20,21,21
     > ,20,0,-21,-27, -21,-15,-22,-10,2,8,13,13,12,11,-7,-25,-21, -10
     > ,-2,-4,-3.5,2,6,12,9,8,6,-10,-20,-10, -5,1,0,-2,-2,1,10,12,11,5
     > ,-14,-16,-5, -5,3,5,2,-7,-5,7,17,17,2,-20,-16,-5, -7,7,16,10,-12
     > ,-11,6,25,30,-7.5,-36,-15,-7, -10,10,25,20,-11,-11,7,30,42,-40
     > ,-52,-33,-10, -12,10,27,25,-7,-7,7,32,50,-70,-60,-35,-12/
      !.. Find the index in the longitude array
       DO 20 I=1,NLONG
         IF(GLONG.LT.ALONG(I)) GO TO 30
  20   CONTINUE
  30   CONTINUE
       J=I-1
      !.. Find the index in the latitude array
       DO 40 I=1,NLAT
         IF(GLAT.LT.ALAT(I)) GO TO 50
 40    CONTINUE
 50    CONTINUE
       K=I-1
      !..-- Now do the bilinear fit (Press et al. Num. Rec. page 96)
       RLAT=(GLAT-ALAT(K))/(ALAT(K+1)-ALAT(K))
       RLONG=(GLONG-ALONG(J))/(ALONG(J+1)-ALONG(J))
       DECPT=(1-RLAT)*(1-RLONG)*DECLS(J,K)+RLONG*(1-RLAT)*DECLS(J+1,K)
     >   +RLAT*RLONG*DECLS(J+1,K+1)+(1-RLONG)*RLAT*DECLS(J,K+1)
       END
C:::::::::::::::::::::::::::: RUNPRN ::::::::::::::::::::::::::::::::::::::
      !..- This subroutine prints out the characteristics of the FLIP run
      !..- Written by Phil Richards in February 1993
      SUBROUTINE RUNPRN(IFILE,ISKP,ZLBDY,HPEQ,FPAS,HFAC,UVFAC,
     >  AUR_FILE,PCO)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      INTEGER RCLPD,AUR_FILE,FEXISTS
      REAL UVFAC(59),RHMF2N,RNMF2N,FOE
      REAL AMP_RING,TMAX,SIGT,LMAX,SIGL,RC_TIME_CON,PROP
      COMMON/CONFAC/IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP
      COMMON/RC/AMP_RING,TMAX,SIGT,LMAX,SIGL,RC_TIME_CON,PROP
C
      WRITE(IFILE,5)
 5    FORMAT(/1X,9('*'),' Interhemispheric solutions along magnetic'
     > ,1X,'flux tube ',9('*'))

      IF(IHPOP.EQ.1) WRITE(IFILE,*) '  ON: O+ & H+ is calculated'
      IF(IHPOP.NE.1) WRITE(IFILE,*)
     >  ' OFF: O+ & H+ from initial conditions only'

      IF(ITEMP.EQ.1) WRITE(IFILE,*) '  ON: Te & Ti is calculated'
      IF(ITEMP.NE.1) WRITE(IFILE,*) 
     >  ' OFF: Te & Ti from initial conditions only'


      IF(ITEMP.EQ.0.AND.IHPOP.EQ.0) WRITE(IFILE,*)
     >   '  WARNING !! No Te, Ti, [H+], [O+]'

      IF(INPLUS.EQ.1) WRITE(IFILE,*) '  ON: N+ is calculated'
      IF(INPLUS.NE.1) WRITE(IFILE,*) ' OFF: No N+ calculation'

      IF(IHEP.EQ.1) WRITE(IFILE,*) '  ON: He+ is calculated'
      IF(IHEP.NE.1) WRITE(IFILE,*) ' OFF: No He+ calculation'

      WRITE(IFILE,38)
 38   FORMAT(/1X,9('*'),' Solutions on vertical grid in Northern and'
     >  ,1X,'Southern hemisphere ',9('*'))

      IF(IODDN.EQ.1) WRITE(IFILE,*) 
     >  '  ON: Solve N(2D), N(4S), & NO Continuity & momentum equations'
      IF(IODDN.NE.1) WRITE(IFILE,*) 
     >  ' OFF: N(2D), N(4S), and NO from chemical equilibrium densities'

      IF(IVIBN2.EQ.1) WRITE(IFILE,*) 
     >  '  ON: Continuity & momentum for vibrationally excited N2'
      IF(IVIBN2.EQ.0) WRITE(IFILE,*) 
     >  ' OFF: No vibrationally excited N2 densities available'

      IF(IWIND.EQ.0) WRITE(IFILE,55)
      IF(IWIND.EQ.-1) WRITE(IFILE,60)
      IF(IWIND.EQ.-2) WRITE(IFILE,62)
      IF(IWIND.EQ.-3) WRITE(IFILE,63)

      !... print IRI version
      IF(IWIND.LT.0) CALL IRIHMF(999.0,289.0,70.0
     >   ,300.0,1,12.0,RHMF2N,RCLPD,RNMF2N,FOE)
      IF(IWIND.NE.0.AND.PCO.LT.1.1) 
     >  WRITE(IFILE,*) ' ** CAUTION, hmF2 alg. may be unreliable L<1.1'

 55   FORMAT(/1X,'WINDS: Using Hedin''s HWM14 wind model.')
 60     FORMAT(/1X,'WINDS: Using IRI hmF2 to derive winds with'
     >   ,1X,'Richards method, JGR 1991 pp 17839.')
 62     FORMAT(/1X,'NMF2: Using IRI hmF2 and NmF2 to normalize'
     >   ,1X,' FLIP density using Richards GRL 1995 algorithm.')
 63     FORMAT(/1X,'NMF2: Using IRI hmF2 and NmF2 to normalize'
     >   ,1X,' FLIP density with Richards et al. JGR,1998 MSIS'
     >   ,1X,'modification algorithm.')

      WRITE(IFILE,67) ZLBDY
 67   FORMAT(/1X,'Z0: The absolute lower boundary =',F6.1,' km')
      WRITE(IFILE,70) FPAS
 70   FORMAT(1X,'FPAS: Plasmaspheric pitch angle trapping factor ='
     >  ,F5.2)

      IF(HFAC.GT.0) WRITE(IFILE,75) HFAC
 75   FORMAT(/3X,'HFAC: PLASMAspheric heating multiplied by = ',F5.2)
      IF(HFAC.LT.0) WRITE(IFILE,76) -HFAC
 76   FORMAT(3X,'HFAC: IONOspheric heating multiplied by = ',F5.2)

      IF(AUR_FILE.GT.0) WRITE(IFILE,*)
     >  ' AURORA: There is an auroral parameter file on this run'

      IF(AMP_RING.GT.0) WRITE(IFILE,86)
 86   FORMAT(/3X,' There is RING CURRENT heating on this run')

      !.. see if EUV data file exists 
      CALL TFILE(42,FEXISTS)
      IF(FEXISTS.NE.0) WRITE(IFILE,87) 
 87   FORMAT(/3X,'+++ EUV fluxes read from file +++')

      WRITE(IFILE,*) ' '  !.. adds a blank line at the end

      RETURN
      END
C::::::::::::::::::::::::::::: IRIPRN :::::::::::::::::::::::::::::
      !..-- Output for comparison with IRI model 20 APR 1992
      SUBROUTINE IRIPRN(IRIHED,SEC,DT,SATN,CHIN,NMF2N,DNMF2N,HMF2N
     >   ,DHMF2N,SATF,CHIF,NMF2S,DNMF2S,HMF2S,DHMF2S,IDAY
     >   ,GLATN,GLONN,GLATF,GLONF,F107A,N170KMN,N170KMS)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,SATN,SATF,CHIN,CHIF,RDAY,SX
      REAL GLATN,GLONN,GLATF,GLONF,F107A
      REAL RHMF2N,RNMF2N,FOE,RHMF2S,RNMF2S
      INTEGER RCLPD
      IF(IRIHED.NE.1) WRITE(12,517)
      IRIHED=1
      CALL ACTUAL_DAY(IDAY,SEC,JDAY,SX) 
      CALL IRIHMF(GLATN,GLONN,F107A,SNGL(HMF2N),JDAY,SATN
     >       ,RHMF2N,RCLPD,RNMF2N,FOE)
      CALL IRIHMF(GLATF,GLONF,F107A,SNGL(HMF2S),JDAY,SATF
     >       ,RHMF2S,RCLPD,RNMF2S,FOE)
      RDAY=REAL(MOD(IDAY,1000))+SEC/86400.
      WRITE(12,518) SEC/3600.,SATN,57.3*CHIN,NMF2N,DNMF2N
     > ,RNMF2N/1.0E6,HMF2N,DHMF2N,RHMF2N,SATF,57.3*CHIF,NMF2S,DNMF2S
     > ,RNMF2S/1.0E6,HMF2S,DHMF2S,RHMF2S,JDAY
     > ,NMF2N/N170KMN,NMF2S/N170KMS

 517  FORMAT(/5X,' Comparison of FLIP NmF2 and HmF2  with the IRI model'
     >  ,//30X,'Northern hemisphere',50X,'Southern hemisphere'
     >  ,/5X,'UT    LTn    SZAn    NmF2n     Ndatn    N_IRIn   HmF2n'
     >  ,2X,'Hdatn H_IRIn     LTs    SZAs    NmF2s     Ndats     N_IRIs'
     >  ,2X,'HmF2s hdats H_IRIs   DAY    RF2F1_NTH RF2F1_STH')

 518  FORMAT(F8.2,2F7.2,1P,3E10.2,0P,3F7.2,2X,2F7.2,1P,3E10.2,0P,3F7.2
     >  ,I8,1P,2E10.2)
      RETURN
      END
C:::::::::::::::::::::: PRCAUT ::::::::::::::::::::::::::::::::
C------ This subroutine writes a warning about first 12 hours when
C------ creating a new field line
      SUBROUTINE PRCAUT(ICAUTN)
      ICAUTN=1
      WRITE(3,90)
      WRITE(4,90)
      WRITE(6,90)
      WRITE(9,90)
      WRITE(12,90)
      WRITE(56,90)
      WRITE(66,90)
 90   FORMAT(/2X,9('*'),' DISCARD data above this line - too close'
     > ,1X,'to initial conditions ',9('*'))
      END
C::::::::::::::::::::::::::::: EUVPRN :::::::::::::::::::::::::::::
C------ Output for comparison with IRI model 20 APR 1992
      SUBROUTINE EUVPRN(JK,JC,JFREQN,JFREQS,SEC,SATN,CHIN,SATF,CHIF
     >  ,IDAY,F107,F107A,TI,Z)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      IMPLICIT NONE
      INTEGER IDAY,JDAY,JK,JC,JFREQN,JFREQS,ISELECTN,ISELECTS
      REAL SEC,SATN,SATF,CHIN,CHIF,RDAY,SX,F107,F107A
      REAL CE_NTH,CE_STH,LOSS_NTH,LOSS_STH
      REAL ATTEN_NTH,ATTEN_STH,PLTIME,SDAYN,SDAYS
      DOUBLE PRECISION ON,HN,N2N,O2N,TN,TI(3,401),RTS(99),VC,
     >   Z(401),PHIONX
      COMMON/ND/ON(401),HN(401),N2N(401),O2N(401),PHIONX(401),TN(401)

      DATA RDAY/-1.0/,SDAYN,SDAYS/0.0,0.0/
      IF(RDAY.LT.0.0) WRITE(13,517)
 517  FORMAT(/5X,' Diagnostics for O+ chemistry at hmF2'
     >  ,1X,'in the Northern hemisphere'
     >  ,/1X,'  UT      LT     DAY   Jt_hm    Atten   N2_lossf O2_lossf'
     >  ,2X,'Ch_eq    O_N2   O_O2        OX      N2       O2      TINF'
     >  ,4X,'Jt_INF SLT   PHION')

      IF(RDAY.LT.0.0) WRITE(16,518)
 518  FORMAT(/5X,' Diagnostics for O+ chemistry at hmF2'
     >  ,1X,'in the Southern hemisphere'
     >  ,/1X,'  UT      LT     DAY   Jt_hm    Atten   N2_lossf O2_lossf'
     >  ,2X,'Ch_eq    O_N2   O_O2        OX      N2       O2      TINF'
     >  ,4X,'Jt_INF SLT   PHION')

      CALL ACTUAL_DAY(IDAY,SEC,JDAY,SX) 
      RDAY=REAL(MOD(IDAY,1000))+SEC/86400.

      !.... Northern hemisphere chemistry
      CALL RATS(JK,TI(3,JK),TI(2,JK),TN(JK),RTS)
      CE_NTH=PHION(JK)/(N2N(JK)*RTS(3)+O2N(JK)*RTS(4))
      LOSS_NTH=N2N(JK)*RTS(3)+O2N(JK)*RTS(4)
      ATTEN_NTH=EUVION(1,1,JK)*ON(JFREQN)/ON(JK)/
     >  (EUVION(1,1,JFREQN)+1.0E-22)

      !.. Identify a LT of interest in North 
      PLTIME=12.0
      ISELECTN=0  !.. identifies a special time in printout
      IF(INT(RDAY).NE.INT(SDAYN).AND.SATN.GE.PLTIME) THEN
         ISELECTN=1
         SDAYN=RDAY
      ENDIF
      IF(RDAY.GT.-99) WRITE(13,519) SEC/3600.,SATN,JDAY
     > ,PHION(JK)/ON(JK),ATTEN_NTH,N2N(JK)*RTS(3)
     > ,O2N(JK)*RTS(4),CE_NTH,ON(JK)/N2N(JK),ON(JK)/O2N(JK)
     > ,ON(JK),N2N(JK),O2N(JK),TN(JFREQN),PHION(JFREQN)/ON(JFREQN)
     > ,ISELECTN,PHION(JK)
      !... Southern hemisphere - chemistry
      CALL RATS(JK,TI(3,JC),TI(2,JC),TN(JC),RTS)
      CE_STH=PHION(JC)/(N2N(JC)*RTS(3)+O2N(JC)*RTS(4))
      LOSS_STH=N2N(JC)*RTS(3)+O2N(JC)*RTS(4)
      ATTEN_STH=EUVION(1,1,JC)*ON(JFREQS)/ON(JC)/
     >  (EUVION(1,1,JFREQS)+1.0E-22)

      !.. Identify a LT of interest in South 
      PLTIME=12.0
      ISELECTS=0  !.. identifies a special time in printout
      IF(INT(RDAY).NE.INT(SDAYS).AND.SATF.GE.PLTIME) THEN
         ISELECTS=1
         SDAYS=RDAY
      ENDIF
      IF(RDAY.GT.-99) WRITE(16,519) SEC/3600.,SATF,JDAY
     > ,PHION(JC)/ON(JC),ATTEN_STH,N2N(JC)*RTS(3)
     > ,O2N(JC)*RTS(4),CE_STH,ON(JC)/N2N(JC),ON(JC)/O2N(JC)
     > ,ON(JC),N2N(JC),O2N(JC),TN(JFREQS),PHION(JFREQS)/ON(JFREQS)
     > ,ISELECTN,PHION(JC)

 519  FORMAT(F8.2,F6.2,I8,1P,5E9.2,0P,2F7.2,2X,1P,3E9.2
     >  ,0P,F9.2,1P,E9.2,I3,7E9.2,0P,3F8.1,1P,E9.2)
      RETURN
      END
C:::::::::::::::::::::::: MPOLE_LOC ::::::::::::::::::::::::::::::::
C... This function calls the routines that determine the latitude and 
C... longitude of the pseudo north magnetic pole. These are different
C... for the two hemispheres
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_LOC(GLONG,     !... Geographic long. in degrees
     >                     GLAT,      !... Geographic lat. in degrees
     >                     IDAY,      !... year
     >                     PLAT,PLON) !... pole position in radians
      IMPLICIT NONE
      INTEGER IDAY
      REAL GLAT,GLONG,PLAT,PLON
      IF(GLAT.GE.0) CALL MPOLE_NTH(GLONG,IDAY,PLAT,PLON)
      IF(GLAT.LT.0) CALL MPOLE_STH(GLONG,IDAY,PLAT,PLON)
      RETURN
      END
C:::::::::::::::::::::::: MPOLE_NTH ::::::::::::::::::::::::::::::::
C... This function returns the latitude and longitude of the
C... pseudo north magnetic pole for a tilted, centred, dipole. The 
C... pole position varies with geographic location on the Earth so 
C... as to provide fair agreement with L values from the IGRF model.
C... This routine is for the northern latitudes. There is a separate
C... routine for southern latitudes. It is most accurate for L=3 and 
C... not very accurate at the Atlantic longitudes.
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_NTH(GLONG,      !... Geographic long. in degrees
     >                     IDAY,      !... year
     >                     PLAT,PLON)  !... pole position in radians
      IMPLICIT NONE
      INTEGER I,DIM,IDAY
      PARAMETER (DIM=19)
      REAL GLONG,PLAT,PLON,DLAT,DLON,GEOLON(DIM),MAGLON(DIM)
      REAL MAGLAT(DIM),GRAD,XCUT

      DATA GEOLON/0,20,40,60,80,100,120,140,160,180,200,220,240
     > ,260,280,300,320,340,360/

      DATA MAGLAT/78.0,78.0,80.0,81.0,83.0,83.5,83.5,82.0,81.0
     > ,80.0,80.0,79.0,78.0,77.0,77.5,79.0,80.0,79.0,78.0/

      DATA MAGLON/270.0,280.0,284.0,295.0,300.0,295.0,295.0,295.0,295.0
     > ,295.0,295.0,293.0,290.0,290.0,290.0,280.0,267.0,265.0,270.0/


      IF(GLONG.LT.0) GLONG=GLONG+360
      IF(GLONG.GT.360) GLONG=GLONG-360
      PLAT=MAGLAT(1)
      PLON=MAGLON(1)

      !..Loop backwards to locate the GLONG within the GLONG array
      I=DIM+1
  10  I=I-1
      IF(GLONG.LT.GEOLON(I)) GO TO 10

      !... Now find the longitude by linear interpolation
      GRAD=(MAGLON(I+1)-MAGLON(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLON(I)-GEOLON(I)*GRAD
      PLON=GLONG*GRAD+XCUT

      !... now find the latitude by linear interpolation
      GRAD=(MAGLAT(I+1)-MAGLAT(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLAT(I)-GEOLON(I)*GRAD
      PLAT=GLONG*GRAD+XCUT

      CALL MPOLE_INC(IDAY,DLAT,DLON)   !..find change since 1980

      PLON=(DLON+PLON)/57.2958   !.. convert to radians
      PLAT=(DLAT+PLAT)/57.2958

      RETURN
      END
C:::::::::::::::::::::::: MPOLE_STH ::::::::::::::::::::::::::::::::
C... This function returns the latitude and longitude of the
C... pseudo north magnetic pole for a tilted, centred, dipole. The 
C... pole position varies with geographic location on the Earth so 
C... as to provide fair agreement with L values from the IGRF model.
C... This routine is for the southern latitudes. There is a separate
C... routine for northern latitudes. It is most accurate for L=3 and 
C... not very accurate at the Atlantic longitudes.
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_STH(GLONG,      !... Geographic long. in degrees
     >                     IDAY,      !... year
     >                     PLAT,PLON)  !... pole position in radians
      IMPLICIT NONE
      INTEGER I,DIM,IDAY
      PARAMETER (DIM=19)
      REAL GLONG,PLAT,PLON,DLAT,DLON,GEOLON(DIM),MAGLON(DIM)
      REAL MAGLAT(DIM),GRAD,XCUT

      DATA GEOLON/0,20,40,60,80,100,120,140,160,180,200,220,240
     > ,260,280,300,320,340,360/

      DATA MAGLON/285.0,280.0,280.0,280.0,280.0,290.0,300.0,310.0,310.0
     > ,305.0,300.0,297.0,290.0,290.0,290.0,285.0,285.0,285.0,285.0/

      DATA MAGLAT/78.0,78.0,79.0,79.0,78.50,77.0,77.0,78.0,79.0,80.0
     > ,80.0,79.0,78.0,76.0,75.0,75.0,78.0,78.0,78.0/


      IF(GLONG.LT.0) GLONG=GLONG+360
      IF(GLONG.GT.360) GLONG=GLONG-360
      PLAT=MAGLAT(1)
      PLON=MAGLON(1)

      !.. Loop backwards to locate the GLONG within the GLONG array
      I=DIM+1
  10  I=I-1
      IF(GLONG.LT.GEOLON(I)) GO TO 10

      !... Now find the longitude by linear interpolation
      GRAD=(MAGLON(I+1)-MAGLON(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLON(I)-GEOLON(I)*GRAD
      PLON=GLONG*GRAD+XCUT

      !... now find the latitude by linear interpolation
      GRAD=(MAGLAT(I+1)-MAGLAT(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLAT(I)-GEOLON(I)*GRAD
      PLAT=GLONG*GRAD+XCUT

      !... Year correction not needed in South
      PLON=PLON/57.2958   !.. convert to radians
      PLAT=PLAT/57.2958

      RETURN
      END
C:::::::::::::::::::: MPOLE_INC ::::::::::::::::::::::::::::::
C... This function finds returns the change in latitude and longitude 
C... of the north magnetic pole referred to its location in 1980.
C... This is from the NSSDC web site
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_INC(IDAY,DLAT,DLON)
      IMPLICIT NONE
      INTEGER I,DIM,IDAY
      PARAMETER (DIM=12)
      REAL YEAR,PLAT,PLON,MAGYR(DIM),MAGLON(DIM),MAGLAT(DIM),GRAD,XCUT
      REAL DLAT,DLON
      !... Years that data are location specified
      DATA MAGYR/1950,1955,1960,1965,1970,1975,1980,1985,1990,1995,
     >           2000,2005/
      !... Magnetic longitude of magnetic pole
      DATA MAGLON/281.05,280.87,280.66,280.41,280.14,279.84,279.44
     >            ,279.11,278.76,278.35,278.06,277.89/
      !... Magnetic latitude of magnetic pole
      DATA MAGLAT/79.85,79.95,80.05,80.13,80.22,80.34,80.55,80.80
     >             ,81.04,81.28,81.64,81.99/

      YEAR=IDAY/1000
      !... make sure the year is within the legal range
      IF(YEAR.LT.MAGYR(1)) THEN
        WRITE(6,90) NINT(YEAR),NINT(MAGYR(1)),NINT(MAGYR(1))
        PLAT=MAGLAT(1)
        PLON=MAGLON(1)
        RETURN
      ELSE IF(YEAR.GE.MAGYR(DIM)) THEN
        WRITE(6,91) NINT(YEAR),NINT(MAGYR(DIM)),NINT(MAGYR(DIM))
        PLAT=MAGLAT(DIM)
        PLON=MAGLON(DIM)
        RETURN
      ENDIF

 90   FORMAT(1X,' *** Requested year',I5,' is less than',I5,
     >  '; the magnetic pole location has been set to the',I5,
     >  ' value')
 91   FORMAT(1X,' *** Requested year',I5,' is greater than',I5,
     >  '; the magnetic pole location has been set to the',I5,
     >  ' value')

      !.. Loop backwards to locate the year within the year array
      I=DIM+1
  10  I=I-1
      IF(YEAR.LT.MAGYR(I)) GO TO 10

      !... Now find the longitude by linear interpolation
      GRAD=(MAGLON(I+1)-MAGLON(I))/(MAGYR(I+1)-MAGYR(I))
      XCUT=MAGLON(I)-MAGYR(I)*GRAD
      PLON=YEAR*GRAD+XCUT

      !... now find the latitude by linear interpolation
      GRAD=(MAGLAT(I+1)-MAGLAT(I))/(MAGYR(I+1)-MAGYR(I))
      XCUT=MAGLAT(I)-MAGYR(I)*GRAD
      PLAT=YEAR*GRAD+XCUT

      !.. determine the change in pole position since 1980 in radians
      DLON=(PLON-MAGLON(7))
      DLAT=(PLAT-MAGLAT(7))
      RETURN
      END
C::::::::::::::::::::::::::::::::: RUN_ERROR :::::::::::::::::::::::::::::::::
C..   This routine is used to put a generic warning in all the FLIP run output 
C..   files. Written by P. Richards, April 2008
      SUBROUTINE RUN_ERROR

      IMPLICIT NONE
      INTEGER UNITNUM(5)     !.. Unit numbers for writing error message
      INTEGER K

      DATA UNITNUM/3,4,8,9,17/ !.. Unit numbers for writing error message

      DO K=1,5
        WRITE(UNITNUM(K),199)
      ENDDO

 199   FORMAT(//3X,'**** ERROR **** ERROR **** ERROR **** '
     >  ,'Check the'//3X,' TOP and BOTTOM of the output file:- '
     >  ' *LOG*.txt (Unit FOR006) for detailed messages ***')
       WRITE(6,299)
 299  FORMAT(//'  --HELP--HELP--HELP--HELP--HELP--HELP--'
     > ,/2X,'If the FLIP model stops prematurely because of'
     > ,/2X,'Non-Convergence or other errors try modifying the control'
     > ,/2X,'parameters in FLIPRUN.DDD. First try increasing or'
     > ,/2X,'decreasing DTMAX. You can change DTSET, TWIDT, JTIMAX,'
     > ,/2X,'and Z0. All these should have only a small effect on the'
     > ,/2X,'accuracy.'
     > ,//5X,'If this does not clear the problem, you can also'
     > ,/2X,'try changing IWIND, HPEQ, HEPRAT, FPAS one at a time.'
     > ,/2X,'These will have a greater effect on the model values but'
     > ,/2X,'may still be acceptable. For example, the standard is to'
     > ,/2X,'use IRI hmF2 for winds, but the HWM model winds may be OK.'
     > ,/2X,'You can also try switching off various solutions such as'
     > ,/2X,'IODDN and IVIBN2 if these are not important for the'
     > ,/2X,'calculations.'
     > ,//5X,'If you are reading data from a file for winds, hmF2,'
     > ,/2X,'ExB drifts, etc., it may be possible to overcome the'
     > ,/2X,'problem by modifying some of the values in the data file.')
      RETURN
      END
