C.................... RSDENA_EXB.FOR ..................................
C----------------------------------------------------------------
C   this program was written by Eugene Young in 1978 but has been substantially
C   modified since then by P. Richards. It is responsible for setting up the
C   continuity equations for minimization. it calls routines
C   VEL to get velocities and fluxes,
C   CHEMO for interpolated chemical sources and sinks, and 
C   DAVE for interpolated dn/dt.
C   It sets up the integrated continuity equations as the function F.
C.. Reference  Torr, et al., J. Geophys. Res., 95, 21,147-21,168, 1990.
      SUBROUTINE DFIJ(J,    !.. point on the field line
     >              JSJ,    !.. used to print diagnostics if JSJ NE 0 or 1
     >               DT,    !.. Time step in seconds
     >                N,    !.. O+, H+, He+, minor ion densities array
     >               TI,    !.. Ion and electron temperatures
     >                F,    !.. OUTPUT: The integrated continuity equation
     >            NSAVE)    !.. ion densities at previous time step
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      IMPLICIT NONE
      INTEGER ION      !.. number of ions to solve = 2 (O+,H+)
      INTEGER I,IUN    !.. I= loop control, IUN= unit to write
      INTEGER J,JSJ    !.. See I/O parameters
      DOUBLE PRECISION DT,N(4,FLDIM),TI(3,FLDIM),F(20),NSAVE(2,FLDIM)
      DOUBLE PRECISION VL(2),VU(2)     !.. Velocities at upper and lower points
      DOUBLE PRECISION Q(2),L(2)       !.. Integrated Sources and Sinks
      DOUBLE PRECISION FLU(2),FLL(2),FGR      !.. Flux quantities
      DOUBLE PRECISION TINCR(2),ANM(2),PM(2)  !.. dn/dt quantities
      DOUBLE PRECISION BU,BL        !.. B field at upper and lower 1/2 points
      DATA FLU(1)/0.0/, FLU(2)/0.0/, ION/2/

      !.. Call VEL for the velocities (VL,VU) and fluxes (FLL,FLU) 
      CALL VEL(J-1,ION,VL,FLL,N,TI,JSJ)       !.. Lower 1/2 pt
      CALL VEL(J,ION,VU,FLU,N,TI,JSJ)         !.. Upper 1/2 pt

      !..  Call CHEMO to evaluate the source (Q) & sink (L) terms
      CALL CHEMO(JSJ,ION,J,Q,N,NSAVE,TI,L)

      !..  Call DAVE for dn/dt. PM(AM)=future(present) cpt of dn/dt
      CALL DAVE(ION,J,ANM,PM,N,NSAVE)

      !.. Magnetic field at Upper and Lower 1/2 points
      BU=(BM(J)+BM(J+1))*0.5
      BL=(BM(J)+BM(J-1))*0.5

      !.. this section sets up the continuity eqns F(i) for the 2 ions
      DO I=1,ION
        TINCR(I)=(PM(I)-ANM(I))/DT      !.. dn/dt
        FGR=(FLU(I)/BU-FLL(I)/BL)       !.. flux

!nm20150516(3): production=0
!      q(i)=0.0
!nm20150516(4): loss=0
!      L(I)=0.0
!nm20150516(5): flux=0
!       FGR=0.0
!nm20150516(5): flux=0
!       TINCR(I)=0.0

        F(I)=Q(I)-L(I)-TINCR(I)-FGR     !.. continuity equation

        !.. Store velocity for use in other routines on call from DLOOPS 
        !.. (JSJ=0) but not on call from DMATRX
        IF(JSJ.EQ.0) XIONV(I,J)=VU(I)

      ENDDO
      RETURN
      END
C:::::::::::::::::::::::::::::::: VEL ::::::::::::::::::::::::::::::::
C-------------------------------------------------------------------
C   This function was written by Eugene Young (1978). it calculates the
C   velocities at the midpoint between j and j+1 from the ion momentum eqtns
C   for O+ and H+. Parameters from the grid point z(j) are interpolated to
C   get their values  at the half point. The parameters are from 
C   Refer to St. Maurice and Schunk Planet. Space Sci. [1977]
C.. Refer to Torr, et al., J. Geophys. Res., 95, 21,147-21,168, [1990].
C-------------------------------------------------------------------
      SUBROUTINE VEL(J,    !.. point on the field line
     >             ION,    !.. # of ions (2)
     >               V,    !.. Ion velocities
     >            FLUX,    !.. Ion Fluxes 
     >               ion_density,    !.. O+, H+, He+, minor ion densities array
     >              TI,    !.. Ion and electron temperatures
     >             JSJ)    !.. used to print diagnostics if JSJ NE 0 or 1
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc.  
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      IMPLICIT NONE
      INTEGER ION        !.. number of ions to solve = 2 (O+,H+)
      INTEGER I,K,IUN    !.. I, K= loop control, IUN= unit to write
      INTEGER J,JSJ      !.. See I/O parameters
      DOUBLE PRECISION ion_density(4,FLDIM),TI(3,FLDIM),V(2),FLUX(2)
      DOUBLE PRECISION MASS(4),A(4)      !.. Mass and atomic # of O+,H+, He+
      DOUBLE PRECISION PP(4),NEU,NEL,NE  !.. Interpolated ions and electrons
      !.. Collision stuff from routine JP      
      DOUBLE PRECISION D(4),GAMMA(4),ALPHA(4,4),ALPHAS(4,4),B(4,4)
      !.. These are used for the terms in the momentum equation
      DOUBLE PRECISION GRADNE,ALPP,GRA,DLNI,QION(2),Q1(2),Q2,Q3
      DOUBLE PRECISION DET,F(20)
      DOUBLE PRECISION QSIGN(2)
      DATA QSIGN/-1.,1./
      DATA MASS/26.7616E-24,1.6726E-24,6.6904E-24,0.0/,A/16,1,4,0/

      !.. Exponential interpolation for densities
      PP(1)=DSQRT(ion_density(1,J+1)*ion_density(1,J))          !.. Midpoint O+ density
      PP(2)=DSQRT(ion_density(2,J+1)*ion_density(2,J))          !.. Midpoint H+ density
      PP(3)=DSQRT(XIONN(3,J+1)*XIONN(3,J))  !.. Midpoint He+ density

      NEU=ion_density(1,J+1)+ion_density(2,J+1)+ion_density(3,J+1)+
     >    XIONN(3,J+1)  !.. upper e density
      NEL=ion_density(1,J)+ion_density(2,J)+ion_density(3,J)+
     >    XIONN(3,J)          !.. lower e density
      !.. Midpoint electron densities
      NE=DSQRT(ion_density(1,J+1)*ion_density(1,J))+
     >   DSQRT(ion_density(2,J+1)*ion_density(2,J))+
     >   DSQRT(ion_density(3,J+1)*ion_density(3,J))+
     >   DSQRT(XIONN(3,J+1)*XIONN(3,J))

      !.. call JP to obtain D(i),GAMMA,ALPHA,ALPHAS,and B for V(i)
      CALL MAURICE_SCHUNK_1977(J,PP,D,GAMMA,ALPHA,ALPHAS,B,JSJ,MASS,A)

      !-- calculate altitude derivatives not dependent on ion density
      GRADNE=TEJ(J)*(NEU-NEL)/NE
      ALPP=(ALPHA(1,2)-ALPHAS(1,2))*GRADTI(J)/(PP(1)+PP(2))

      !.. Terms in the momentum equation. Consult St. Maurice 
      !.. and Schunk [1977] equation (19)
      DO I=1,ION
        DLNI=(ion_density(I,J+1)-ion_density(I,J))/DS(J)/PP(I)  !.. ion density gradient
        GRA=-MASS(I)*GRAV(J)                !.. gravity
        K=3-I
        !ALPP=(ALPHA(I,K)-ALPHAS(I,K))*GRADTI(J)/(PP(1)+PP(2))
        QION(I)=QSIGN(I)*(GAMMA(I)*GRADTE(J)-PP(K)*ALPP)
        Q1(I)=DLNI+GRA+GRADTI(J)+GRADNE+GRADTE(J)+QION(I)
        Q2=UNJ(J)*B(I,K)          !.. neutral wind term
!nm20110815: should do the test for collision term with He+: 
!2nd step: take out "*0.000" to have nonzero value or Q3
        Q3=(XIONV(3,J+1)+XIONV(3,J))*B(I,3)*00.000   !.. He+ term
        F(I)=-D(I)*Q1(I)+Q2+Q3       !.. Momentum equation
      ENDDO

      !.. Add He+ collisions to neutral wind collisions
!nm20110815: should do the test for collision term with He+: 
!1st step: UNcomment the following two lines:
c      B(1,2)=B(1,2)+B(1,3)
c      B(2,1)=B(2,1)+B(2,3)
      !-- now calculate the velocities and fluxes 
      DET=B(1,2)+B(2,1)+B(1,2)*B(2,1)
      V(1)=(F(1)+F(2)+B(2,1)*F(1))/DET
      V(2)=(F(1)+F(2)+B(1,2)*F(2))/DET
      FLUX(1)=V(1)*PP(1)
      FLUX(2)=V(2)*PP(2)

      !.. Print the momentum equation terms ..........
      IF(JSJ.EQ.-4) THEN
        !.. write header, but only once
        IF(IUN.NE.172) THEN
          IUN=172
          WRITE(IUN,665) 
        ENDIF
        WRITE(IUN,605) J,NINT(Z(J)),DLNI,GRADTI(J),GRADNE
     >     ,GRADTE(J),GRA,QION(1),Q1(1),F(1),Q2,NE
     >     ,NINT(TI(3,J)),NINT(TI(2,J)),GAMMA(1),GAMMA(2)
        IF(I.EQ.2) WRITE(IUN,*) '   ' !.. blank line
        WRITE(IUN,606) J,Z(J),V(1),V(2),DET,F(1),F(2)
      ENDIF

 605   FORMAT(1X,'VEL',I4,I7,1P,10E10.2,2I6,9E10.2)
 606  FORMAT(1X,'J=',I4,3X,'ALT=',F9.0,2X,'V1=',1P,E13.6,2X,'V2=',E13.6
     >  ,2X,'DET=',E13.6,2X,'F1=',E13.6,2X,'F2=',E13.6)
 665  FORMAT(6X,'J    ALT    DLNI      GRADTI    GRADNE'
     > ,3X,'GRADTE     GRAV      QION       Q1         F'
     > ,8X,'Q2        NE      Te    Ti  GAMMA(1)  GAMMA(2)')
      RETURN
      END
C:::::::::::::::::::::::::::: JP :::::::::::::::::::::::::::::::::
C--- This routine provides the parameters on p910 of St. Maurice and 
C--- Schunk Planet. Space Sci. [1977] equations 20,21,22,23,24,25,26
C--- The ratio of neutral-ion collision frequencies is also calculated 
C--- and returned in array B.
C--- Note this version includes collisions with He+
C-------------------------------------------------------------------
      SUBROUTINE MAURICE_SCHUNK_1977(J,    !.. point on the field line
     >             PP,    !.. interpolated ion and e densities
     >              D,    !.. Diffusion coefficient
     >          GAMMA,    !.. See St. Maurice ref
     >          ALPHA,    !.. See St. Maurice ref
     >         ALPHAS,    !.. See St. Maurice ref
     >              B,    !.. ratio of neutral-ion collision frequencies
     >            JSJ,    !.. used to print diagnostics if JSJ NE 0 or 1
     >           MASS,    !.. Ion masses
     >              A)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values  
      IMPLICIT  NONE
      INTEGER J,JSJ    !.. See I/O parameters
      INTEGER I,K,L        !.. loop control variables, IUN unit #
      INTEGER IUN          !.. IUN unit #
      !.. Collision related stuff (see I/O comments)      
      DOUBLE PRECISION PP(4),D(4),GAMMA(4),ALPHA(4,4),ALPHAS(4,4),
     >   B(4,4),MASS(4),A(4)
      DOUBLE PRECISION AMR(4,4)           !.. Reduced mass M1*M2/(M1+M2)
      DOUBLE PRECISION NU(4,4),NUP(4,4)   !.. Ion-ion collision frequencies
      DOUBLE PRECISION D1(4,4),D4(4,4)    !.. St-M & Sch page 908
      DOUBLE PRECISION BK,X,DELTA(4,4)    !.. Boltzmann, St-M & Sch page 910
      DATA BK/1.3807E-16/                 !.. Boltzmann
      DATA AMR/8.0,0.9412,3.20,7.4667,0.9412,0.50,0.80,0.9333,3.20,
     >  0.80,2.00,3.1111,7.4667,0.9333,3.1111,7.00/
      DATA D1/0.725,-0.1612,0.008,0.6213,2.6623,0.725,1.928,2.6187,
     >  1.928,0.008,0.725,1.8222,0.8347,-0.1547,0.0444,0.725/
      DATA D4/-0.075,-0.0789,-0.192,-0.112,0.9799,-0.075,0.528,
     >  0.952,0.528,-0.192,-0.075,0.4667,-0.032,-0.088,-0.20,-0.075/

      !.. ion-ion collisions from St.Maurice and Schunk  p920
      DO I=1,2
        DO K=1,3
          NU(I,K)=1.2726*PP(K)*DSQRT(AMR(I,K))/(A(I)*TIJ(J)**1.5)
        ENDDO
      ENDDO

      !.. NUP=NU' (from St. Maurice and Schunk p920). NOTE: Equations 
      !.. simplified because T(O+)=T(H+)
      DO I=1,2
        NUP(I,3)=1.25*NU(I,3)*(D4(I,3)+1.5*AMR(I,3)/A(I))  !.. for He+
        K=3-I
        DO L=1,2
          IF(I.NE.L) NUP(I,L)=1.25*NU(I,L)*(D4(I,L)+1.5*AMR(I,L)/A(I))
          IF(I.EQ.L) NUP(I,L)=NU(I,L)+
     >     1.25*NU(I,K)*(D1(I,K)+1.5*AMR(I,K)/A(I))
        ENDDO
      ENDDO

      !.. Factors for St. Maurice and Schunk [1977] equations 24, 25, 26
      !.. St. Maurice and Schunk [1977] equations 24, 25, 26
      !.. Note; for Ti=Tj, the alpha terms in eqns 18 and 19 are the same
      DO I=1,2
        K=3-I
        ALPHA(I,K)=1.875*((PP(I)+PP(K))/PP(K))*(AMR(I,K)/A(I))
     >    *(NU(I,K)*(NUP(K,K)-NUP(K,I)))
     >    /(NUP(I,I)*NUP(K,K)-NUP(I,K)*NUP(K,I))     
        ALPHAS(I,K)=ALPHA(I,K)*(MASS(I)/MASS(K))**2*
     >    (NUP(I,I)-NUP(I,K))/(NUP(K,K)-NUP(K,I))
        DELTA(I,K)=0.6*(PP(K)/(PP(I)+PP(K)))*(AMR(I,K)/A(I))
     >    *(ALPHA(I,K)+ALPHAS(I,K)*A(I)*PP(I)/A(K)/PP(K))
      ENDDO

      !.. Calculate diffusion coefficients and neutral-ion collisions
      DO I=1,2
        GAMMA(I)=0.0   !... gamma is small, set to zero
        K=3-I
        !--- diffusion coeff D(i) St. Maurice and Schunk p910
        D(I)=BK*TIJ(J)/MASS(I)/NU(I,K)/(1.0-DELTA(I,K))
        !.. Ratio of neutral - ion collision freqs.
        B(I,K)=NUX(I,J)/NU(I,K)/(1.0-DELTA(I,K))
        !.. Ratio of He+ - ion collision freqs.
        B(I,3)=NU(I,3)/NU(I,K)/(1.0-DELTA(I,K))
      ENDDO
     
      !..... output diagnostics
      IF(JSJ.EQ.-4) THEN
        IF(IUN.NE.172) THEN
          IUN=172
          WRITE(IUN,665) 
        ENDIF
        WRITE(IUN,99) J,PP(1),PP(2),NU(1,2),NU(2,1),D(1)
     >  ,D(2),B(1,2),B(2,1),DELTA(1,2),DELTA(2,1),ALPHAS(1,2),
     >   ALPHAS(2,1),ALPHA(1,2),ALPHA(2,1),NU(1,3),NU(2,3)
        IF(I.EQ.2) WRITE(IUN,*) '   '
      ENDIF
 665  FORMAT(11X,'PP(1)       PP(2)     NU(1,2)    NU(2,1)     D(1)'
     >  ,7X,'D(2)       B(1,2)     B(2,1)    DELTA1     DELTA2'
     >  ,5X,'ALS12      ALS21       AL12     AL21     NU(1,3)'
     >  ,5X,'NU(2,3)')
 99   FORMAT(1X,'JP',I4,1P,22E11.3)
      RETURN
      END
C::::::::::::::::::::::::::::: AVDEN :::::::::::::::::::::::::::::::::::
C.. This function evaluates the interpolated densities at the midpoints
C.. for those quantities that do not change during the solution procedure. 
C.. It also evaluates the ion-neutral collision frequencies NUX(i,j)  
C.. and the O+ and H+ production and loss rates
C.. .....................................................................
      SUBROUTINE AVDEN(TI,    !.. Ion and electron temperature
     >                JSJ,    !.. used to print diagnostics if JSJ NE 0 or 1
     >                ZLO,    !.. Lower boundary for printing
     >                ZHI)    !.. Upper boundary for printing
      USE FIELD_LINE_GRID   !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc.
      IMPLICIT NONE
      INTEGER J,JSJ,IUN    !.. loop control variables, IUN unit #
      DOUBLE PRECISION TI(3,FLDIM),ZLO,ZHI,RTS(99)
      DOUBLE PRECISION BK
      !.. Values only used here to calculate collision frequencies
      DOUBLE PRECISION OA,HA,N2A,O2A,TNJ,TR
      !.. Ion-neutral coll frequency parameters
      DOUBLE PRECISION RHPH,ROPO,RHPO,ROPH,CHPN,COPN
      DATA BK/1.3807E-16/

      DO J=JMIN,JMAX-1
        !.. Values stored in AVE_PARAMS
        DS(J)=SL(J+1)-SL(J)
        TIJ(J)=0.5*(TI(1,J)+TI(1,J+1))
        TEJ(J)=0.5*(TI(3,J)+TI(3,J+1))/DS(J)/TIJ(J)
        GRAV(J)=0.5*(GR(J)+GR(J+1))/BK/TIJ(J)
        GRADTI(J)=(TI(1,J+1)-TI(1,J))/DS(J)/TIJ(J)
        GRADTE(J)=(TI(3,J+1)-TI(3,J))/DS(J)/TIJ(J)
        UNJ(J)=0.5*(UN(J+1)+UN(J))

        !.. Values only used here to calculate collision frequencies
        OA=SQRT(ON(J)*ON(J+1))
        HA=SQRT(HN(J)*HN(J+1))
        N2A=SQRT(N2N(J)*N2N(J+1))
        O2A=SQRT(O2N(J)*O2N(J+1))
        TNJ=0.5*(TN(J)+TN(J+1))
        TR=(TIJ(J)+TNJ)*0.5

        !.. Ion-neutral coll freqs taken from Schunk and Nagy Rev.
        !.. Geophys., v18, p813, 1980. 
	  !.. First, resonant charge exchange from table 5 p823
        !.. H-H+ resonant charge exchange
        RHPH=2.65E-10*HA*SQRT(TR)*(1-0.083*DLOG10(TR))**2
        !.. O-O+ resonant charge exchange. COLFAC is a scaling factor 
        !.. for O+ - O collision frequency  
        ROPO= COLFAC * 3.67E-11*OA*SQRT(TR)*(1-0.064*DLOG10(TR))**2
        !.. H-O+ resonant charge exchange
        RHPO=6.61E-11*OA*SQRT(TIJ(J))*(1-0.047*DLOG10(TIJ(J)))**2
        !.. O+-H resonant charge exchange
        ROPH=4.63E-12*HA*SQRT(TNJ+TIJ(J)/16.0)

        !.. H+-O2 non-resonant ion neutral interactions table 6
        CHPN=(33.6*N2A+32*O2A)*1.0E-10
        !.. O+-O2 non-resonant ion neutral interactions table 6
        COPN=(6.28*N2A+6.64*O2A)*1.0E-10

        !.. Ion-neutral collision frequencies for O+ and H+
        NUX(1,J)=ROPO+ROPH+COPN
        NUX(2,J)=RHPH+RHPO+CHPN

        !.. get reaction rates for production and loss 
        CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
        !.. Production and loss rates for O+ and H+
        OLOSS(J)=RTS(2)*HN(J)+RTS(3)*N2N(J)+RTS(4)*O2N(J)
        HLOSS(J)=RTS(1)*ON(J)
        HPROD(J)=RTS(2)*HN(J)

        !.. Printing diagnostics
        IF(JSJ.EQ.-4.AND.Z(J).LT.ZLO.AND.Z(J).GT.ZHI) THEN
          IUN=34
          WRITE(IUN,99) J,Z(J),OA,HA,N2A,O2A,TIJ(J),TNJ,TR,RHPH,ROPO,
     >      RHPO
          WRITE(IUN,99) J,Z(J),ROPH,CHPN,COPN,NUX(1,J),NUX(2,J),OLOSS(J)
     >     ,HLOSS(J),HPROD(J)
        ENDIF
      ENDDO
 99   FORMAT(1X,'AVDEN',I4,F7.0,1P,22E10.3)
      RETURN
      END
C:::::::::::::::::::::::::::: TERLIN :::::::::::::::::::::::::::::::::::::::::::::::
C.. This routine integrates the quantity Q along the field line using
C.. linear interpolation between adjacent grid points
C.. This subroutine calculates the integral of Q represented by
C.. the 3 values QJM1, QJ, QJP1 from SL=(S(J)+S(J-1))/2 to
C.. SU=(S(J+1)+S(J))/2. The integral is done in 2 parts,
C.. from SL to SJ and from SJ to SU.
C.. The integration is done from the lower mid point to the upper midpoint.
C.. The lower mid point is halfway between points j and j-1.
C.. The upper mid point is halfway between points j and j+1. 
C.. The dist between the j and j-1 = SL and that between j and j+1 = SU. 
C.. The integration is done from SL/2 to SU/2. Hence the factor 0.25.
      SUBROUTINE TERLIN(QJM1, !.. Input value at point J-1
     >                  QJ,   !.. Input value at point J
     >                QJP1,   !.. Input value at point J+1
     >                  SL,   !.. Distance between points J and J-1
     >                  SU,   !.. Distance between points J and J+1
     >                  QM)   !.. OUTPUT: Intergrated value over the interval
      IMPLICIT NONE
      INTEGER J
      DOUBLE PRECISION QJM1,QJ,QJP1,SL,SU,QL,QM,QU
        QL=.5*(QJM1+QJ)   !.. Average value on the lower interval
        QU=.5*(QJ+QJP1)   !.. Average value on the Upper interval
        !.. Integrated value
        QM=.25*((QL+QJ)*SL+(QU+QJ)*SU)
      RETURN
      END
C::::::::::::::::::::::::::::::::: TEREXP ::::::::::::::::::::::::::::::::::::::::::
C.. This routine integrates the quantity Q along the field line using
C.. exponential interpolation between adjacent grid points
C.. This subroutine calculates the integral of Q represented by
C.. the 3 values QJM1, QJ, QJP1 from SL=(S(J)+S(J-1))/2 to
C.. SU=(S(J+1)+S(J))/2. The integral is done in 2 parts,
C.. from SL to SJ and from SJ to SU.
C.. For exponential interpolation, the mean value of A and B is SQRT(A*B)
C.. The integration is done from the lower mid point to the upper midpoint.
C.. The lower mid point is halfway between points j and j-1.
C.. The upper mid point is halfway between points j and j+1. 
C.. The dist between the j and j-1 = SL and that between j and j+1 = SU. 
C.. The integration is done from SL/2 to SU/2. Hence the factor 0.25.
      SUBROUTINE 
     >  TEREXP(QJM1,   !.. Input value at point J-1
     >           QJ,   !.. Input value at point J
     >         QJP1,   !.. Input value at point J+1
     >           SL,   !.. Distance between points J and J-1
     >           SU,   !.. Distance between points J and J+1
     >           QM)   !.. OUTPUT: Intergrated value over the interval
      IMPLICIT NONE
      INTEGER J
      DOUBLE PRECISION QJM1,QJ,QJP1,SL,SU,QL,QM,QU
        QL=DSQRT(QJM1*QJ)   !.. Average value on the lower interval
        QU=DSQRT(QJ*QJP1)   !.. Average value on the Upper interval
        !.. Integrated value
        QM=0.5*(SU*(QJ-QU)/DLOG(QJ/QU)+SL*(QL-QJ)/DLOG(QL/QJ))      
      RETURN
      END
C:::::::::::::::::::::::::::::: DAVE :::::::::::::::::::::::::::::::::::::::::::::
C.. This routine calculates the integrated components of the integral of dn/dt
C.. It returns the ante and post values of half interval.
      SUBROUTINE DAVE(ION,    !.. # of ions (2)
     >                  J,    !.. point on the field line
     >                ANM,    !.. interpolated value at time t
     >                 PM,    !.. interpolated value at time t+DT
     >                  N,    !.. O+, H+, He+, minor ion densities array
     >              NSAVE)    !.. ion densities at previous time step
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc. 
      IMPLICIT NONE
      INTEGER ION,J,I
      DOUBLE PRECISION ANM(2),PM(2),N(4,FLDIM),NSAVE(2,FLDIM)
      !.. values returned by the interpolation
      DOUBLE PRECISION A,B,C,QM

      !.. Perform the integration for O+ and H+ time derivative     
      DO I=1,ION
        !.. Ante values at time t
        B=NSAVE(I,J-1)/BM(J-1)
        C=NSAVE(I,J)/BM(J)
        A=NSAVE(I,J+1)/BM(J+1)
        CALL TEREXP(B,C,A,DS(J-1),DS(J),QM)
        ANM(I)=QM
        !.. Post values at time t+DT
        B=N(I,J-1)/BM(J-1)
        C=N(I,J)/BM(J)
        A=N(I,J+1)/BM(J+1) 
        CALL TEREXP(B,C,A,DS(J-1),DS(J),QM)
        PM(I)=QM
      ENDDO
      RETURN
      END
C:::::::::::::::::::::::::::::::: CHEMO :::::::::::::::::::::::::::::::::::::::::::
C.. this program determines the integrated production and loss
C.. processes. it calls rates to get the rate constants and terd
C.. to do the interpolation
      SUBROUTINE CHEMO(JSJ,  !.. used to print diagnostics if JSJ NE 0 or 1
     >                 ION,  !.. # of ions (2)
     >                   J,  !.. point on the field line
     >              SOURCE,  !.. OUTPUT: O+, H+ production rate
     >                   N,  !.. O+, H+, He+, minor ion densities array
     >               NSAVE,  !.. ion densities at previous time step
     >                  TI,  !.. Ion and electron temperature
     >                SINK)  !.. OUTPUT: O+, H+ loss rate
      USE FIELD_LINE_GRID  !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE     !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production, PHION
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc.
      IMPLICIT NONE
      INTEGER I,K,IUN    !.. loop controls and print unit
      INTEGER JSJ,ION,J
      DOUBLE PRECISION SOURCE(2),N(4,FLDIM),NSAVE(2,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION SINK(2),Q(2,3),L(2,3)

      !-- Set up the sources (Q) and sinks (L) divided by B field at points 
      !-- j-1, j, and j+1
      DO I=1,3
        K=I+J-2
        Q(1,I)=(HLOSS(K)*N(2,K)+PHION(K))/BM(K)  !.. O+ production
        Q(2,I)=HPROD(K)*N(1,K)/BM(K)             !.. H+ production
        L(1,I)=OLOSS(K)*N(1,K)/BM(K)             !.. O+ loss
        L(2,I)=HLOSS(K)*N(2,K)/BM(K)             !.. O+ loss
      ENDDO

      !.. Interpolate the sources and sinks
      DO I=1,ION
        CALL TERLIN(Q(I,1),Q(I,2),Q(I,3),DS(J-1),DS(J),SOURCE(I))
        CALL TERLIN(L(I,1),L(I,2),L(I,3),DS(J-1),DS(J),SINK(I))
      ENDDO

      !.. Print diagnostics
      IF(JSJ.EQ.-4) THEN
        IUN=35
        WRITE(IUN,99) J,Z(J),Q(2,1),Q(2,2),Q(2,3),SOURCE(1)
     >   ,SOURCE(2),SINK(1),SINK(2),SL(J),DS(J)
      ENDIF
      RETURN
 99   FORMAT(1X,'CHEMO',I4,F7.0,1P,22E10.3)
      END
C:::::::::::::::::::::::::::::::: GET_GAMMA :::::::::::::::::::::::::::::::::::::::::::
C... This routine provides gamma from St. Maurice and Schunk Planet. 
C... Gamma is small for O+ and H+ ions
C... Call from VEL   CALL GET_GAMMA(PP(1),PP(2),TEJ,NE,GAMMA)
      SUBROUTINE GET_GAMMA(N1,  !.. ion density 1
     >                     N2,  !.. ion density 2
     >                     TE,  !.. electron temperature
     >                     NE,  !.. electron density
     >                  GAMMA)  !.. OUTPUT: Gamma
      IMPLICIT NONE
      DOUBLE PRECISION N1,N2,TE,NE,GSIGN,GAMMA(4)
      !.. Collision frequencies
      DOUBLE PRECISION NU1,NU2,NUE,NU12

      NU1=54.5*N1/TE**1.5
      NU2=54.5*N2/TE**1.5
      NUE=38.537*NE/TE**1.5
      NU12=NU1+NU2
	GAMMA(1)=15.0*(NE*NU1/N1-NU12)/(13*NU12+8*NUE)
	GAMMA(2)=15.0*(NU12-NE*NU2/N2)/(13*NU12+8*NUE)
      RETURN
      END
