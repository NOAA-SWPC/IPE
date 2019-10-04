C.................... CMINOR.FOR;16        26-APR-1994 14:53:26.59 
C....... This subroutine evaluates the minor ion and neutral densities
C....... from chemical equilibrium and is also used for printing the
C....... production and loss rates for the FLIP model. P. Richards
C....... Sep 1989
      SUBROUTINE CMINOR(I,J,JP,IHEPLS,INPLS,INNO,FD,ID,N,TI,Z,EFLAG,
     >                  mp,lp)                       
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE MINORNEUT !.. N4S N2D NNO N2P N2A O1D O1S
      IMPLICIT NONE
      integer mp,lp
      INTEGER I,II,J,K,ITS,ID,JITER,JP
      INTEGER INNO  !.. switch to turn on FLIP NO calculation if <0
      INTEGER IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on if > 0
      INTEGER EFLAG(11,11)         !.. error flags
      !.. These variables used in Newton solver for [e]
      DOUBLE PRECISION DNOP,DO2P,ZNE,ZNESAV,B,C,H,DEX,DX
      DOUBLE PRECISION DISN2D,DISN4S,DISNP,UVDISN,DISN2S
      DOUBLE PRECISION HEPLUS,PRHEP
      DOUBLE PRECISION TOP2DI,TOP2PI,TOTO2I
      DOUBLE PRECISION PDISOP,PLYNOP,PN2PEL,PNOP,PO2P
      DOUBLE PRECISION PO1DEL,PO1DSR,PO1SEL,PSRNO
      DOUBLE PRECISION PROD(6),LOSS(6),RTS(99),FEX(2),FD(9)
      DOUBLE PRECISION NOPLUS,N2PLUS,O2PLUS,NPLUS,OP2D,OP2P
      DOUBLE PRECISION N(4,IDIM),TI(3,IDIM),F(20),Z(IDIM)
      DOUBLE PRECISION BOD3,O2SS,O2B1,O2DEL,O2ADEL,O2ASIG
      DOUBLE PRECISION EUVP4S,PEOP4S,PEOP2P

      DO II=1,9
        FD(II)=0.0D0
      ENDDO

      !.. Don't calculate minor species at high altitudes
      IF(Z(J).LT.80.0.OR.Z(J).GT.700) RETURN

      !.. Initialize ion densities with stored values
      NOPLUS=XIONN(5,J)
      O2PLUS=XIONN(6,J)
      N2PLUS=XIONN(7,J)
      OP2D=XIONN(8,J)
      OP2P=XIONN(9,J)
      HEPLUS=XIONN(3,J)
      NPLUS=XIONN(4,J)

      !.. Production rates, consult subroutine prdint
      TOP2PI=SUMION(1,9,J)
      TOP2DI=SUMION(1,8,J)
      TOTO2I=SUMION(2,7,J)
      DISNP=SUMION(3,4,J)+SUMION(3,5,J)+SUMION(3,6,J)
      DISN4S=SUMEXC(3,5,J)+DISNP*0.5
      DISN2D=SUMEXC(3,6,J)+DISNP*0.5
      PDISOP=SUMION(2,4,J)+SUMION(2,5,J)
      UVDISN=OTHPR1(1,J)
      PRHEP=OTHPR1(2,J)
      PO1DSR=OTHPR1(3,J)
      PO1SEL=SUMEXC(1,2,J)
      PN2PEL=SUMEXC(3,7,J)
      PO1DEL=SUMEXC(1,1,J)
      !..  Calculate hv + NO ionization frequency from Lyman-a
      PLYNOP=OTHPR2(1,J)
      !..  Calculate diss and ion rates for NOP - see also CODDN
      !..  Calculate hv + NO dissociation frequency from hv + O2 dissoc rate
      !..  Lyman-a NO+ freq used for NO dissoc by Lyman-a
      PSRNO=OTHPR2(2,J)

      !.. - Odd nitrogen zeroed to prevent convergence problems at low altitude
      !.. IF(Z(J).LT.100) NNO(J)=0.0
      !.. IF(Z(J).LT.100) N4S(J)=0.0

      !.. C*****  obtain reaction rates from subr rats
      CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)

	!.. RTS(99) is used in the RVN2PB file for N2+(v) calculation
      RTS(99)=RTS(42)*RTS(10)

      !.. ... Only recalculate minor ions if not already stored
      IF(I.EQ.1.OR.I.EQ.-1.OR.I.GT.2) GO TO 12

      !.. . First guess for [e] using previous stored values
! GHGM just O+ and H+ here also.....
      ZNE=XIONN(1,J)+XIONN(2,J)+XIONN(3,J)+XIONN(4,J)+XIONN(5,J)+
     >  XIONN(6,J)+XIONN(7,J)
! GHGM
!     ZNE=XIONN(1,J)+XIONN(2,J)
      ZNESAV=ZNE
      CALL CO2P(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PO2P,DO2P
     > ,TOTO2I,N(1,J),OP2D,N2PLUS,NPLUS,N4S(J),NNO(J),OP2P,mp,lp)

      !.. .... no+
      CALL CNOP(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PNOP
     > ,DNOP,N(1,J),N2PLUS,O2PLUS,N4S(J),NNO(J),NPLUS,N2P(J)
     >  ,PLYNOP,N2D(J),OP2D,mp,lp)

      !.. calculation of [e] analytically for first guess in Newton
      B=N(1,J)+N(2,J)+N2PLUS
      C=PNOP/RTS(5)+PO2P/RTS(6)
      ZNE=(B+DSQRT(B**2+4.0*C))/2.0

      !.. Loop through chemistry twice to evaluate dF/dh for Newton
      JITER=0
 9    DO 10 ITS=1,2
        IF(JITER.EQ.0) H=ZNE*0.001
        IF(ITS.EQ.2) ZNE=ZNE+H
      !....... o+(2p)
      CALL COP2P(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,OP2P
     > ,TOP2PI,0.0D0,HEPLUS,N4S(J),NNO(J))
      !....... o+(2d)
      CALL COP2D(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,OP2D
     > ,EUVION(1,8,J),PEPION(1,8,J),OP2P,HEPLUS,N4S(J),NNO(J))

      !...... n2+
      CALL CN2PLS(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,N2PLUS,EUVION(3,1,J),EUVION(3,2,J),EUVION(3,3,J)
     > ,PEPION(3,1,J)+PAUION(3,1,J),PEPION(3,2,J)+PAUION(3,2,J)
     > ,PEPION(3,3,J)+PAUION(3,3,J),OP2D,OP2P,HEPLUS
     > ,NPLUS,NNO(J),N4S(J))

      !...... o2+
      CALL CO2P(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PO2P
     > ,O2PLUS,TOTO2I,N(1,J),OP2D,N2PLUS,NPLUS,N4S(J)
     > ,NNO(J),OP2P,mp,lp)

      !...... no+
      CALL CNOP(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PNOP
     > ,NOPLUS,N(1,J),N2PLUS,O2PLUS,N4S(J),NNO(J),NPLUS,N2P(J)
     > ,PLYNOP,N2D(J),OP2D,mp,lp)

      FEX(ITS)=ZNE-N(1,J)-N(2,J)-NOPLUS-O2PLUS-N2PLUS-
     >  OP2D-OP2P-NPLUS-HEPLUS
 10   CONTINUE

      JITER=JITER+1
      DEX=(FEX(2)-FEX(1))/H
      ZNE=ZNE-H-FEX(1)/DEX
! GHGM - stop ZNE from going to zero....
      if(zne.lt.0.000001) zne = zne + 0.000001
! GHGM
      IF(JITER.GT.9) GO TO 14
      IF(ABS(FEX(1)/DEX/ZNE).GT.1.0E-2) GO TO 9
      DX=ABS(FEX(1)/DEX)

      EFLAG(5,1)=0
 14   IF(JITER.GT.9) THEN
         IF(EFLAG(11,11).EQ.1) WRITE(6,57) J, Z(J),ZNE,ZNESAV,N(1,J)
         EFLAG(5,1)=-1
         RETURN
      ENDIF
 12   CONTINUE

 57   FORMAT(' ** NEWTON non-convergence in CMINOR at PT. J=',I3
     >  ,' ALT= ',F8.2,' NE= ',1P,E11.2,' NESAV= ',1P,E11.2,' O+ ='
     >  ,1P,E11.2)

      !...n2(a3sigma)
      CALL CN2A(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,N2A(J),SUMEXC(3,1,J),SUMEXC(3,2,J),SUMEXC(3,3,J),SUMEXC(3,4,J))
      !.. densities for O2 atmospheric bands
        CALL CO2(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,BOD3,O1D(J),O2SS,TN(J),O2B1,O2DEL,O2ADEL,O2ASIG)
      !..O(1S)
      CALL CO1S(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,O1S(J),O2PLUS,N4S(J),PO1SEL,NPLUS,N2A(J),PO1DSR,O2ADEL,O2DEL)

      !.... local equilibrium densities for metastables and odd nitrogen 

      !.. N(2P) density
      CALL CN2P(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PROD(5)
     > ,LOSS(5),N2P(J),PN2PEL,UVDISN,O2PLUS,NNO(J),N2PLUS)
      IF(LOSS(5).GT.0) N2P(J)=PROD(5)/LOSS(5)

      !.. N(2D) density
      CALL CN2D(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),NOPLUS,ZNE,PROD(1)
     > ,LOSS(1),N2PLUS,DISN2D,UVDISN,NPLUS,N2P(J)
     > ,N2D(J),N(1,J),NNO(J),N2A(J))
      IF(LOSS(1).GT.0) N2D(J)=PROD(1)/LOSS(1)

      !.... NO density
      CALL CNO(J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PROD(3),LOSS(3)
     >  ,N2D(J),N4S(J),N2P(J),NNO(J),O2PLUS,N(1,J),PSRNO,PLYNOP,N2A(J),
     >  NPLUS)
      !.. Calculate chemical equilibrium NO density if INNO < 0
      IF(INNO.LT.0.AND.LOSS(3).GT.0) THEN
        NNO(J)=PROD(3)/LOSS(3)
        !.. Set a floor on NO density, which is needed below ~150 km at night 
        IF(NNO(J).LT.1.0E8*EXP((100-Z(J))/20)) 
     >    NNO(J)=1.0E8*EXP((100-Z(J))/20)
        IF(NNO(J).GT.1.5E8) NNO(J)=1.5E8      !.. Don't let NO get too big
      ENDIF

      !... O(1D) density
      CALL CO1D(0,J,0,0,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE,PROD(4)
     > ,LOSS(4),O2PLUS,N2D(J),PO1DEL,PO1DSR,O1S(J),NPLUS,O1D(J))
      IF(LOSS(4).GT.0) O1D(J)=PROD(4)/LOSS(4)

      !... Chemical He+ density
      IF(IHEPLS.LT.0) CALL CHEP(J,0,JP,IHEPLS,Z(J),RTS,ON(J),O2N(J),
     >  N2N(J),ZNE,HEPLUS,PRHEP,HE(J),NNO(J),N(1,J),N(2,J),HN(J))
      IF(IHEPLS.LE.0) XIONN(3,J)=HEPLUS

      !... Chemical N+ density
      IF(INPLS.LT.0) CALL CNPLS(J,0,JP,INPLS,Z(J),RTS,ON(J),O2N(J)
     > ,N2N(J),ZNE,DISNP,NPLUS,N(1,J),N2D(J),OP2P,HEPLUS,OTHPR2(3,J)
     > ,O2PLUS,N4S(J),OP2D,N2PLUS,NNO(J))
      IF(INPLS.LE.0) XIONN(4,J)=NPLUS

      !.. EQN2D(J) used in PE2S only
      EQN2D(J)=N2D(J)*RTS(8)*ZNE

      !... calc o+ prod. from minor ions      !$$$
      FD(5)=OP2P*(RTS(26)*ON(J)+RTS(14)*ZNE+0.047)
      FD(6)=OP2D*(RTS(12)*ZNE+RTS(28)*ON(J))
      FD(7)=N2PLUS*RTS(99)*ON(J)
      FD(8)=NPLUS*ON(J)*RTS(31)
      FD(9)=FD(5)+FD(6)+FD(7)+FD(8)

      !.. Store minor ion densities
      XIONN(5,J)=NOPLUS
      XIONN(6,J)=O2PLUS
      XIONN(7,J)=N2PLUS
      XIONN(8,J)=OP2D
      XIONN(9,J)=OP2P

      IF(I.LT.2) RETURN

      !.. this next section for printing densities and production rates
      DISN2S=SUMION(3,6,J)

      IF(I.EQ.2) CALL CION(J,ID,JP,Z(J),N(1,J),O2PLUS,NOPLUS,N2PLUS
     > ,NPLUS,OP2D,OP2P,ZNE,TI(3,J),TI(1,J),HEPLUS
     > ,N(2,J),NNO(J))

      IF(I.EQ.3) CALL CNEUT(J,ID,JP,Z(J),O1S(J),N2A(J),N2D(J),O1D(J),
     > NNO(J),N4S(J),N2P(J))

      IF(I.EQ.5) CALL CO2EM(J,ID,JP,Z(J),O2B1,O2DEL,O2SS,O2ADEL,O2ASIG)

      IF(I.EQ.6) CALL CN2D(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),NOPLUS
     > ,ZNE,PROD(1),LOSS(1),N2PLUS,DISN2D,UVDISN,NPLUS,N2P(J)
     > ,N2D(J),N(1,J),NNO(J),N2A(J))

      !.. --- O(1D) production and loss rates
      IF(I.EQ.7) CALL CO1D(0,J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J)
     > ,ZNE,PROD(4),LOSS(4),O2PLUS,N2D(J),PO1DEL,PO1DSR,O1S(J)
     > ,NPLUS,O1D(J))

      !.. --- O(1D) emission rates
      IF(I.EQ.36) CALL CO1D(1,J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J)
     > ,ZNE,PROD(4),LOSS(4),O2PLUS,N2D(J),PO1DEL,PO1DSR,O1S(J)
     > ,NPLUS,O1D(J))

      IF(I.EQ.8) CALL CNO(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,PROD(3),LOSS(3),N2D(J),N4S(J),N2P(J),NNO(J),O2PLUS,N(1,J)
     >  ,PSRNO,PLYNOP,N2A(J),NPLUS)

      IF(I.EQ.9) CALL CN4S(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,PROD(2),LOSS(2),N4S(J),DISN4S,N2D(J),N2P(J),N(1,J),N2PLUS
     > ,UVDISN,NOPLUS,NPLUS,NNO(J),O2PLUS,PSRNO)

      IF(I.EQ.10) CALL CO1S(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,O1S(J),O2PLUS,N4S(J),PO1SEL,NPLUS,N2A(J),PO1DSR,O2ADEL,O2DEL)

      IF(I.EQ.11) CALL CN2PLS(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J)
     > ,ZNE,N2PLUS,EUVION(3,1,J),EUVION(3,2,J),EUVION(3,3,J)
     > ,PEPION(3,1,J)+PAUION(3,1,J),PEPION(3,2,J)+PAUION(3,2,J)
     > ,PEPION(3,3,J)+PAUION(3,3,J),OP2D,OP2P,HEPLUS
     > ,NPLUS,NNO(J),N4S(J))

      IF(I.EQ.12) CALL CNOP(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,PNOP,NOPLUS,N(1,J),N2PLUS,O2PLUS,N4S(J),NNO(J),NPLUS
     > ,N2P(J),PLYNOP,N2D(J),OP2D,mp,lp)

      IF(I.EQ.13) CALL CO2P(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,PO2P,O2PLUS,TOTO2I,N(1,J),OP2D,N2PLUS,NPLUS,N4S(J)
     > ,NNO(J),OP2P,mp,lp)

      IF(I.EQ.14) THEN
        EUVP4S=EUVION(1,7,J)
        PEOP4S=PEPION(1,7,J)+PAUION(1,7,J)
        CALL COP4S(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,N(1,J),EUVP4S,OP2D,OP2P,PEOP4S,PDISOP,N2PLUS
     > ,N2D(J),NNO(J),HEPLUS,NPLUS)
      ENDIF

      IF(I.EQ.15) CALL COP2D(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     >,OP2D,EUVION(1,8,J),PEPION(1,8,J),OP2P,HEPLUS,N4S(J),NNO(J))

      IF(I.EQ.16) THEN
        PEOP2P=PEPION(1,9,J)+PAUION(1,9,J)
        CALL COP2P(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     >   ,OP2P,TOP2PI,PEOP2P,HEPLUS,N4S(J),NNO(J))
      ENDIF

      IF(I.EQ.17) CALL CNPLS(J,ID,JP,INPLS,Z(J),RTS,ON(J),O2N(J),N2N(J)
     > ,ZNE,DISNP,NPLUS,N(1,J),N2D(J),OP2P,HEPLUS,OTHPR2(3,J)
     > ,O2PLUS,N4S(J),OP2D,N2PLUS,NNO(J))

      IF(I.EQ.18) CALL CN2A(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,N2A(J),SUMEXC(3,1,J),SUMEXC(3,2,J),SUMEXC(3,3,J),SUMEXC(3,4,J))

      IF(I.EQ.19) CALL CN2P(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,PROD(5),LOSS(5),N2P(J),PN2PEL,UVDISN,O2PLUS,NNO(J),N2PLUS)

      IF(I.EQ.20) CALL CHEP(J,ID,JP,IHEPLS,Z(J),RTS,ON(J),O2N(J),N2N(J)
     > ,ZNE,HEPLUS,PRHEP,HE(J),NNO(J),N(1,J),N(2,J),HN(J))

      IF(I.EQ.25) CALL CO2(J,ID,JP,Z(J),RTS,ON(J),O2N(J),N2N(J),ZNE
     > ,BOD3,O1D(J),O2SS,TN(J),O2B1,O2DEL,O2ADEL,O2ASIG)

      !.. IF(I.EQ.26) N2+(v) no longer available

      IF(I.EQ.27) CALL PRDPRB(J,ID,JP,Z(J))


      RETURN
      END
