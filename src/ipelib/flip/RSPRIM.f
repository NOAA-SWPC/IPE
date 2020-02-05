C.................... RSPRIM.FOR ..................................
C------ This routine evaluates the ionization rates for EUV impact
C------ Modified extensively in February 1993 to include nighttime
C------ light sources instead of using subroutine NYTE.
      SUBROUTINE PRIMPR(F107,F107A,UVFAC,COLUM,EUVFLUX)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION    !.. EUV, photoelectron, and auroral production, PHION
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE MINORNEUT !.. N4S N2D NNO N2P N2A O1D O1S
      USE PHOTOIONIZATION_DATA  !.. NPOT LMAX PROB ZLAM SIGION SIGABS TPOT 
      IMPLICIT NONE
      INTEGER IJ,IS,IK,I,L,K1,K
      DOUBLE PRECISION XN(3),COLUMN(3),TAU,FLUX,HEPLS,FBSBN
     > ,DISN,TAUGAM,FLUXG,ALTG,XNSIGF,DSPECT,CLNITE(3)
      DOUBLE PRECISION COLUM(3,FLDIM)  !.. Only used in PEPRIM
      REAL EUVFLUX(37)
      REAL XSNPLS(37),FNITE(37),UVFAC(59)
      REAL F107,F107A,FNFAC,O2LYXS,O2SRXS,FREQSR,FREQLY,TAUN
     > ,FLUXN
      !.. Fluxes for nighttime ion production in the 37 wavelength bins of
      !.. Torr et al GRL 1979. The fluxes are set to reproduce the production
      !.. rates in Strobel et al. PSS, p1027, 1980. Note that most bins are 
      !.. set to zero and that the Strobel production rates are scaled by 
      !.. FNFAC to stabilize the O+ solution below 200 km. Note also that
      !.. the wavelengths in FNITE go from largest (#3=HI) to smallest.
      DATA FNITE/9E5,0.0,9E5,2*0.0,9E6,13*0.0,3E5,8*0.0,3E5,8*0.0/
      DATA FNFAC/1.0/

      !----------- Loop through the altitudes ----------------
      DO IJ=JMIN,JMAX
        !... initialization of production rates. 1.0E-15 stabilizes 
        !... e density evaluation at low altitudes in CMINOR
        DO IS=1,3
          DO IK=1,12
            EUVION(IS,IK,IJ)=1.0E-15
          ENDDO
        ENDDO

        DISN=0.0
        DO I=1,6
          OTHPR2(I,IJ)=1.0E-15
          OTHPR1(I,IJ)=1.0E-15
        ENDDO

        !.. Nighttime He+ production is calculated and stored. Attenuated to
        !.. avoid excess production at low altitudes
        OTHPR1(2,IJ)= 8E-11* EXP(-1.0E-11*N2N(IJ)) *HE(IJ)
        DO I=1,3
          COLUM(I,IJ)=1.0E+25
        ENDDO

        !.. temporary variables for neutral densities
        XN(1)=ON(IJ)
        XN(2)=O2N(IJ)
        XN(3)=N2N(IJ)

        !.. determine if sun is below the horizon ...
        !.. Must now do calculation for night production - Feb 93
        ALTG=(6371.0+Z(IJ))*SIN(3.1416-SZA(IJ))-6371.0
        !..      IF(SZA(IJ).GT.1.57.AND.ALTG.LT.85.) RETURN
        IF(Z(IJ).GT.1500) GOTO 777

        !.. get column densities for scattered light at night
        !CALL SCOLUM(IJ,0.0D0,Z(IJ),GL(IJ),TN(IJ),XN,CLNITE)
        CALL NCOLUM(IJ,0.0D0,Z(IJ),GL(IJ),TN(IJ),TINF(IJ),XN,CLNITE)

        !.. evaluate the neutral column density 
        !CALL SCOLUM(IJ,SZA(IJ),Z(IJ),GL(IJ),TN(IJ),XN,COLUMN)
        CALL NCOLUM(IJ,SZA(IJ),Z(IJ),GL(IJ),TN(IJ),TINF(IJ),XN,COLUMN)
        !.. Store the column densities for the 2-Stream program
        COLUM(1,IJ)=COLUMN(1)
        COLUM(2,IJ)=COLUMN(2)
        COLUM(3,IJ)=COLUMN(3)

        !.. O2 dissociation by Schumann-Runge UV.
        !.. OTHPR1(3,IJ)= dissociation rate. OTHPR1(5,IJ)= Energy
        CALL SCHUMN(IJ,Z(IJ),O2N(IJ),COLUMN,OTHPR1(3,IJ),OTHPR1(5,IJ),
     >    UVFAC)

        !.. Calculate hv + NO ion. freq. from Lyman-a (Brasseur & Solomon)
        !.. OTHPR2(2,IJ) is photodissociation of NO in the SR bands. 
        !.. A small night production from scattered light is included. FREQLY
        !.. varies with solar activity using Richards et al. 1994 page 8981
        !.. LY_a=2.5E11 (Lean), sigi(NO)=2.0E-18 (Brasseur & Solomon page 329)
        DATA O2LYXS,O2SRXS,FREQSR /1.0E-20,1.0E-21,5.0E-6/
        FREQLY=5.0E-7*(1+4.0E-3*(0.5*(F107+F107A)-80.0))
        OTHPR2(1,IJ)=FREQLY*(EXP(-O2LYXS*COLUMN(2))
     >    +0.001*EXP(-O2LYXS*CLNITE(2)))
        OTHPR2(2,IJ)=FREQSR*(EXP(-O2SRXS*COLUMN(2))
     >    +0.001*EXP(-O2SRXS*CLNITE(2)))

        !..  wavelength loop begins here  ----------
        !..  TAU, TAUN = optical depth for day, night 
        HEPLS=0.0
        DO L=1,LMAX
          TAU=0.
          TAUN=0.0
          DO I=1,3
            TAUN=TAUN+SIGABS(I,L)*CLNITE(I)
            TAU=TAU+SIGABS(I,L)*COLUMN(I)
          ENDDO
          IF(TAU.GT.70.0) TAU=70.0
          IF(TAUN.GT.70.0) TAUN=70.0
          !.. evaluate nighttime flux and daytime flux
!nm20120304: For the case that doesn't converge, try increasing values of nighttime production.
!Nighttime production is not well known. So long as the density in the
!E-region around 110 km does not exceed about 1.0E4 /cc it should be OK.
!If the densities were larger than 1.0E4 /cc at 110 km, they would be
!observed by ionosondes, which they are not. In any case, I don't think
!anyone knows what the E-region densities are at night. 
          !FNFAC=FNFAC_flip
          FLUXN=FNFAC*(F107/75.0)*FNITE(L)*EXP(-TAUN)
          FLUX=EUVFLUX(L)*EXP(-TAU) + FLUXN

          !.... he+ production. He+ X-S  = 0.25 N2  X-S. HEPRDN = nite He+
          IF(ZLAM(L).LT.500.) HEPLS=HEPLS+HE(IJ)*0.25*SIGION(3,L)*FLUX

          !.. hv + N -> N+ + e. ion. freq. Richards et al. JGR 1994 page 8989
          DATA XSNPLS/6*0.0,.211,10.294,11.171,10.961,11.244,11.323,
     >    12.098,13.265,12.423,11.951,11.212,11.798,11.758,11.778,11.772
     >    ,11.503,11.016,10.578,9.556,8.15,8.302,7.298,6.413,6.399,5.192
     >    ,5.725,4.787,3.778,2.3,.878,.286/
          OTHPR2(3,IJ)=OTHPR2(3,IJ)+XSNPLS(L)*1.0E-18*FLUX*N4S(IJ)

          IF(ZLAM(L).GE.600.0) THEN
            !.. calculation of total euv absorption-ionization.....
            FBSBN=FLUX*(SIGABS(3,L)-SIGION(3,L))*XN(3)
            !.. Save energy absorbed in the photodissociative process
            OTHPR1(4,IJ)=OTHPR1(4,IJ)+1.24E+4*FBSBN/ZLAM(L)
            !.. production on atomic nitrogen by dissociation
            DISN=DISN+FBSBN
            !.. take into account the large n2 absorption of lyman gamma(972.54)
            IF(NINT(ZLAM(L)).EQ.975) THEN
              TAUGAM=370E-18*COLUMN(3)
              IF(TAUGAM.GT.70.0)TAUGAM=70.0
              FLUXG=UVFAC(34) *0.82E+9 *EXP(-TAUGAM)
              DISN=DISN+FLUXG*370E-18*XN(3)
            ENDIF
          ENDIF

          !***** species loop begins here *****
          DO I=1,3
            !.. calculation of total euv absorption.....
            !..      PEXCIT(2,8,IJ)=PEXCIT(2,8,IJ)
            !..     >     +1.24E+4*FLUX*SIGABS(I,L)*XN(I)/ZLAM(L)
            XNSIGF=XN(I)*SIGION(I,L)*FLUX
            K1=NPOT(I)

            !.. dspect=# ions formed by w-l l by ionization of k state of species i
            DO K=1,K1
              DSPECT=XNSIGF*PROB(I,K,L)
              !... store ion production rates .....
              EUVION(I,K,IJ)=EUVION(I,K,IJ)+DSPECT

              !.. calculation of ion heating rate......
              EUVION(1,10,IJ)=EUVION(1,10,IJ)+DSPECT*TPOT(I,K)
            ENDDO
          ENDDO
        ENDDO    !..---   wavelength loop ends here   -----------

        !...Store UV disoc of N2 2 atoms produced for every dissociation
        OTHPR1(1,IJ)=2.0*DISN
        !.. Transfer He+ production to storage
        OTHPR1(2,IJ)=OTHPR1(2,IJ)+HEPLS
 777  CONTINUE
      ENDDO  !... End of altitude loop
      END
C::::::::::::::::::::::::: SCHUMN ::::::::::::::::::::::::::::::::::::::::::::::
C... production of o(1d) by schumann-runge bands
C... The fluxes are from the SEE web site data. 
C... P. Richards March 2011
      SUBROUTINE SCHUMN (J,Z,O2N,COLUMN,SCHUPR,SCHUHT,UVFAC)
      IMPLICIT NONE
      INTEGER J,LMAX,LSR
      DOUBLE PRECISION Z,O2N,SCHUPR,SCHUHT,COLUMN(3),HSRX,FLD,SRXSCT
      REAL SRFLUX(8),SRXS(8),SRLAM(8),UVFAC(59)
      !.. Schumann-Runge fluxes * 1.0E-11
      DATA SRFLUX/3.3074,1.9699,1.0785,.7083,.5008,.2890,.1651,.1269/
      DATA SRXS/.5,1.5,3.4,6,10,13,15,11/
      DATA SRLAM/1725,1675,1625,1575,1525,1475,1425,1375/

      !.. lmax=# of lambdas in sub. primpr: schuht=heating: schupr=o(1d) prod
      LMAX=37

      !.. Loop over the SR wavelengths
      DO LSR=1,8
        SRXSCT=1.0E-18*SRXS(LSR)  !.. photoabsorption cross section
        HSRX=SRXSCT*COLUMN(2)
        IF(HSRX.GT.70)HSRX=70
        !.. attentuated solar flux
        FLD=UVFAC(LMAX+LSR)*1.E+11*SRFLUX(LSR)*EXP(-HSRX)
        !.. neutral heating SCHUHT and photodissociation rate SCHUPR
        SCHUHT=SCHUHT+1.24E+4*(FLD*SRXSCT)*O2N/SRLAM(LSR)
        SCHUPR=SCHUPR+FLD*SRXSCT
        !WRITE(6,90) LSR,SRXSCT,FLD,SCHUPR,UVFAC(LMAX+LSR)*SRFLUX(LSR)
      ENDDO

      SCHUPR=O2N*SCHUPR
 90   FORMAT(2X,I5,1P,9E10.2)
      RETURN
      END
C:::::::::::::::::::::::::: FACSR ::::::::::::::::::::::::::::::::
C.... The Schumann-Runge factors are scaled according to F10.7
C.... Using SEE fluxes
C.... P. Richards March 2011
      SUBROUTINE FACSR(UVFAC,F107,F107A)
      DIMENSION UVFAC(59),SRFLUX(8),SEEFAC(8)
      !.. ratio of fluxes for day 2002092 (P=196) and day 2007284 (P=68)
      DATA SEEFAC/1.0844,1.0959,1.1176,1.1364,1.1837,1.1418,1.2151,
     >           1.2826/

      F107AV=(F107+F107A)*0.5
      DO LSR=38,45
        I=LSR-37
        UVFAC(I)=1.0
        A=(SEEFAC(I)-1)/128.0
        B=1-A*68.0
          UVFAC(LSR)=A*F107AV+B
          IF(UVFAC(LSR).LT.0.8) UVFAC(LSR)=0.8
          !WRITE(6,'(I6,9F11.3)') LSR,UVFAC(LSR)
      ENDDO
      DO I=46,50
        UVFAC(I)=1.0
      ENDDO
      RETURN
      END
C:::::::::::::::::::::::::: SUMPRD :::::::::::::::::::::::::
C..... Sum the EUV, PE, and AUR production
      SUBROUTINE SUMPRD(JMIN,JMAX
     &,qiont !units number/m3/s
     &)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      INTEGER,INTENT(IN) :: JMIN,JMAX
      REAL*8,dimension(3,JMIN:JMAX),INTENT(IN) :: qiont !1:O,2:O2,3:N2 !units number/m3/s
      INTEGER :: J,IS,IK
      !--- Branching ratios for ion states were updated Sep 91 by P. Richards
      !--- N2+ from Doering and Goembel, JGR 1991, page 16025. fraction of
      !--- N+ from Richards and Torr JGR 1985, page 9917. fraction of O+ from
      !--- O2 dissociation is from Rapp and Englander Golden J. Chem Phys, 1965
      !--- We still need to find cross sections for the higher levels of O+
      REAL ASPRD(3,6)
      DATA ASPRD/.4,.47,.42, .4,.28,.34, .2,.15,.08,0.0,.05,0.0,0.0,
     &   .05,0.0,0.0,0.0,0.16/

!...
      DO 30 J=JMIN,JMAX

      !nm20151031... Calculate aurora ionization rate.
        DO IS=1,3 !o,o2,n2
         DO IK=1,6
          PAUION(IS,IK,J)=PAUION(IS,IK,J)+qiont(IS,J)*ASPRD(IS,IK)*1.E-6 !convert from m-3 to cm-3
         END DO                  !IK
        END DO                    !IS
!d      if (j==16) then
!d      print *,'PAUION1',pauion(1,1:6,j)
!d      print *,'PAUION2',pauion(2,1:6,j)
!d      print *,'PAUION3',pauion(3,1:6,j)
!d      endif

      !..    add contributions from 2 highest O+ metastables to 3 lowest
         EUVION(1,7,J) = EUVION(1,1,J) + EUVION(1,4,J)
         EUVION(1,8,J) = EUVION(1,2,J) + EUVION(1,5,J)/1.3
         EUVION(1,9,J) = EUVION(1,3,J) + EUVION(1,5,J)/4.3

         PEPION(1,7,J) = PEPION(1,1,J) + PEPION(1,4,J)
         PEPION(1,8,J) = PEPION(1,2,J) + PEPION(1,5,J)/1.3
         PEPION(1,9,J) = PEPION(1,3,J) + PEPION(1,5,J)/4.3

         PAUION(1,7,J) = PAUION(1,1,J) + PAUION(1,4,J)
         PAUION(1,8,J) = PAUION(1,2,J) + PAUION(1,5,J)/1.3
         PAUION(1,9,J) = PAUION(1,3,J) + PAUION(1,5,J)/4.3
 
      !..- SUM non-diss. states to get total O2+ and N2+ production
         EUVION(2,7,J)= EUVION(2,1,J)+EUVION(2,2,J)+EUVION(2,3,J)
         PEPION(2,7,J)= PEPION(2,1,J)+PEPION(2,2,J)+PEPION(2,3,J)
         PAUION(2,7,J)= PAUION(2,1,J)+PAUION(2,2,J)+PAUION(2,3,J)

         EUVION(3,7,J)= EUVION(3,1,J)+EUVION(3,2,J)+EUVION(3,3,J)
         PEPION(3,7,J)= PEPION(3,1,J)+PEPION(3,2,J)+PEPION(3,3,J)
         PAUION(3,7,J)= PAUION(3,1,J)+PAUION(3,2,J)+PAUION(3,3,J)
 30   CONTINUE

      DO 50 I=1,3
      DO 50 K=1,12
      DO 50 J=JMIN,JMAX
        SUMION(I,K,J)=EUVION(I,K,J)+PEPION(I,K,J)+PAUION(I,K,J)
        SUMEXC(I,K,J)=PEXCIT(I,K,J)+PAUEXC(I,K,J)
 50   CONTINUE

      !.. total N+ production for MINORA routine
      DO J=JMIN,JMAX
        NPLSPRD(J)=SUMION(3,4,J)+SUMION(3,5,J)+SUMION(3,6,J)+OTHPR2(3,J)
      ENDDO

      RETURN
      END
C:::::::::::::::::::::::::::::::: FACEUV :::::::::::::::::::::::
C----- This routine uses the EUV scaling from Richards et al.[1994]
C----- The EUVAC flux model is based on the F74113 solar reference
C----- spectrum and Hinteregger's scaling factors. This subroutine
C----- just provides the scaling factors as a function of the proxy
C----- (F107+F107A)/2. Modified in March 2009 to return fluxes as 
C-----  well as factors.
      SUBROUTINE FACEUV(F107,F107A,UVFAC,EUVFLUX)
      IMPLICIT NONE
      INTEGER I,L,LMAX
      REAL F107,F107A,F107AV,A,B
      REAL UVFAC(59),HFG200(37),EUVFLUX(37),ZFX(37)
	!..Ratio of EUVAC fluxes for P=200 and P=80 (day 1974113)
      DATA HFG200/2.202,1.855,2.605,3.334,1.333,17.522,4.176,4.0
     >  ,1.4,3.694,1.791,5.385,1.889,1.899,3.427,2.051,1.392,1.619
     >  ,1.439,2.941,1.399,2.416,1.512,1.365,1.570,1.462,2.537,1.393
     >  ,1.572,1.578,1.681,1.598,1.473,1.530,1.622,1.634,1.525/
      DATA LMAX/37/
      !--  Test to see if need to scale - see DATRD2 subroutine      
         !........... EUV scaling
         F107AV=(F107+F107A)*0.5
         DO I=1,LMAX
           A=(HFG200(I)-1)/120.0
           B=1-A*80.0
           UVFAC(I)=A*F107AV+B
           IF(UVFAC(I).LT.0.8) UVFAC(I)=0.8
         ENDDO

      !.. These fluxes, which are divided by 1.0E+9, are for day 1974113.
      !.. Note!!!! fluxes doubled below 250A, while the shortest wavelength 
      !.. fluxes have been tripled
      DATA ZFX/2.4665,2.1,3.5,1.4746,4.4,3.,3.537,1.625,.758,.702,
     > .26,.17,.141,.36,.23,.342,1.59,.53,.357,1.27,.72,.452,.285,
     > .29,.383,.314,.65,.965,6.9,.8,1.679,.21,.46,3.1,4.8,.45,1.2/

      !.. Determine the actual EUV fluxes. Note that for historical reasons, 
	!.. the UV factors and fluxes run in opposite directions
      DO L=1,LMAX
         EUVFLUX(L)=ZFX(L)*1.0E+9*UVFAC(LMAX+1-L)
      ENDDO
      RETURN
      END
C:::::::::::::::::::::::::::::::: FACSEE :::::::::::::::::::::::
C----- This routine produces scaling factors for the SEE solar EUV 
C----- fluxes. It uses P=(F107+F107A)/2. The scaling is based on a solar  
C----- minimum spectrum for day 2002092 (P=196) and day 2007284 (P=68)
C----- Written by P. Richards January 2009
      SUBROUTINE FACSEE(F107,F107A,UVFAC,EUVFLUX)
      IMPLICIT NONE
      INTEGER I,L,LMAX
      REAL F107,F107A,F107AV,A,B
      REAL UVFAC(59),SEEFAC(37),EUVFLUX(37),ZFX(37)
      !.. ratio of fluxes for day 2002092 (P=196) and day 2007284 (P=68)
      DATA SEEFAC/1.82,1.50,1.51,1.78,1.57,14.39,1.84,3.71,1.43,8.02,
     >  1.64,2.94,1.96,1.73,1.90,2.19,1.40,2.15,1.03,3.02,1.51,2.15,
     >  1.65,1.30,1.78,1.27,2.24,1.35,1.29,1.72,1.93,1.85,1.38,2.48,
     >  1.50,1.53,2.58/
      DATA LMAX/37/
      !--  Test to see if need to scale - see DATRD2 subroutine      
         !........... EUV scaling
         F107AV=(F107+F107A)*0.5
         DO I=1,LMAX
           A=(SEEFAC(I)-1)/128.0
           B=1-A*68.0
           UVFAC(I)=A*F107AV+B
           IF(UVFAC(I).LT.0.8) UVFAC(I)=0.8
         ENDDO

      !.. These are SEE fluxes for day 2007284 (P=68), divided by 1.0E+9
      DATA ZFX/3.106,2.1,3.5,1.888,4.4,4.433,4.734,2.42,0.731,0.7,0.26,
     >  0.1832,0.353,0.36,0.5454,0.712,0.795,0.265,0.349,0.791,0.5,
     >  0.775,0.817,0.29,0.649,1.051,0.65,0.4824,6.3736,0.739,4.51,
     >  0.21,0.46,3.558,6.26,0.28,0.4721/

      !.. Determine the actual EUV fluxes. Note that for historical reasons, 
	!.. the UV factors and fluxes run in opposite directions
      DO L=1,LMAX
         EUVFLUX(L)=ZFX(L)*1.0E+9*UVFAC(LMAX+1-L)
      ENDDO
      RETURN
      END
C::::::::::::::::::::::::: NCOLUM :::::::::::::::::::::::::::::::::::::
C++++ this routine evaluates the neutral column density
C++++ see smith & smith jgr 1972 p 3592, subr ambs is the msis
C++++ neutral atmosphere model for z<120 km, chi=solar zenith
C++++ angle, re & ge radius and grav con for earth
C.... Modified by P. Richards September 2010 to eliminate need for calling
C.... the MSIS model at grazing incidence
      SUBROUTINE NCOLUM(J,CHI,ZKM,MLAT,TNIN,TINF,XN,COLUMN)
      IMPLICIT NONE
      INTEGER I,J 
      REAL CHIX, SAT,GLON,GLATD,SN(5),M(3),DG(19),T(2)
      REAL XI,GTN,GN(3)
      DOUBLE PRECISION ZG,CHI,Z,ZKM,TNJ,XN(3),COLUMN(3),ALTG,GE,GR,RE,
     >  RP,SH,XP,Y,ERFY2,CHAPFN,RG,HG,XG,EM,F,G,A,B,C,D
      DOUBLE PRECISION MLAT,TNIN,TINF
      DATA A,B,C,D,F,G/1.0606963,0.55643831,1.0619896,1.7245609
     >  ,0.56498823,0.06651874/
      DATA SN/0.0,0.0,0.0,0.0,0.0/
        DATA    EM   ,   M(1) , M(2) , M(3) ,  RE   , GE
     1    / 1.662E-24 ,   16. ,  32. ,  28. ,6.357E8, 980/
      DATA T,ALTG,ERFY2/0.0,0.0,0.0D0,0.0D0/
      DATA DG/19*0.0/

      TNJ=TNIN  !.. Avoids modifying input TN
      Z=ZKM*1.0E+5  !.. altitude in cm
      DO I=1,3
	  SN(I)=0.0
        COLUMN(I)=1.E+25
	ENDDO

      IF(Z.LT.70*1.0E5) RETURN   !.. below lower boundary of FLIP

      IF(CHI.LT.1.5708) GO  TO 2938      !.. is sza>90.0 degrees

      !.. Calculate grazing incidence parameters. 
      !.. Grazing incidence densities are obtained by extrapolating
      !.. down to grazing incidence altitude
      ALTG=(6371.0E5+Z)*SIN(3.1416-CHI)-6371.0E5
      IF(ALTG.GE.85*1.0E5) THEN
        ZG=ALTG*1.E-5

        !.. Bates-Walker temperature
        XI=(ZG-120.0)*(6357.0+120.0)/(6357.0+ZG)
        !.. GTN is the temperature at grazing incidence
!nm20110126: modified due to ibm xlf compiling error
        GTN=MAX(TINF-(TINF-300.0)*EXP(-0.025*XI),180.0D0)

        !.. Neutral densities are extrapolated from altitude to grazing
        !.. using hydrostatic equilibrium
        !.. altitude. Weighted average Tn and GTn is used 
        GR=GE*(RE/(RE+Z))**2   !.. gravity
        DO I=1,3
          GN(I)=XN(I)*EXP((Z-ALTG)/
     >      ((1.38E-16*(TNIN+GTN*2)/3)/(EM*M(I)*GR)))
        ENDDO
        IF(CHI.GE.1.75.AND.CHI.LT.-1.8)
     .  WRITE(88,'(7F8.2,1P,22E10.2)') Z/1.0E5,ZG,CHI*180/3.1416,
     >      TNIN,TINF,GTN,TNJ
     >,(XN(I),GN(I),I=1,3)  !nm110210
        !.. Make sure that grazing incidence density is not lower than MSIS
        !.. This is for using non MSIS densities like CTIPe
        TNJ=GTN
        DO I=1,3
          SN(I)=GN(I)
          IF(SN(I).LT.XN(I)) SN(I)=XN(I)
        ENDDO
      ELSE
        RETURN
      ENDIF

      !.. sn(1)=o , sn(2)=o2 , sn(3)=n2 , tnj=tn,  gr=gravity, rp=distance
      !.. to pt p, sh=scale height, rg=distance to pt g, hg=scale height at g

 2938 CONTINUE
      GR=GE*(RE/(RE+Z))**2
      RP=RE+Z
      DO 10 I=1,3
      SH=(1.38E-16*TNJ)/(EM*M(I)*GR)
      XP=RP/SH
      Y=SQRT(0.5*XP)*ABS(COS(CHI))
      IF(Y.GT.8) ERFY2=F/(G+Y)
      IF(Y.LT.8) ERFY2=(A+B*Y)/(C+D*Y+Y*Y)
    4 IF(CHI.GT.1.5708)GO  TO 2
      CHAPFN=SQRT(0.5*3.1416*XP)*ERFY2

      COLUMN(I)=XN(I)*SH*CHAPFN
        GO TO 10
    2 RG=RP*SIN(3.1416-CHI)
      HG=1.38E-16*TNJ/
     1    (EM*M(I)*980.*(6371.E5/(6371.E5+ALTG))**2)
      XG=RG/HG
      COLUMN(I)=SQRT(0.5*3.1416*XG)*HG*(2.0*SN(I)-XN(I)*ERFY2)
10       CONTINUE
      RETURN
      END
