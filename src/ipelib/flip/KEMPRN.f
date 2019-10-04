C.................... KEMPRN.FOR;22        26-APR-1994 14:06:10.11 
      SUBROUTINE CION(J,I,JP,Z,OPLS,O2P,NOP,N2PLS,NPLUS,OP2D
     > ,OP2P,NE,TE,TI,HEPLUS,HPLUS,NNO)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      IF(JP.EQ.1) WRITE(I,90)
 90   FORMAT(/'   ALT    [O+]    [O2+]    [NO+]    [N2+]    [N+]'
     >, '    [O+2D]   [O+2P]   [He+]     [H+]    [e]     Te   Ti'
     >,'   [NO]    TBRACE  R_TE_TB')
      CALL BRACE(Z,NE,TE,TB,TEOTB)
      WRITE(I,5) Z,OPLS,O2P,NOP,N2PLS,NPLUS,OP2D,OP2P,HEPLUS,HPLUS
     >  ,NE,NINT(TE),NINT(TI),NNO,TB,TEOTB
      RETURN
 5    FORMAT(F7.2,1P,10E9.2,2I5,0P,F5.2,1P,9E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE BRACE(H,NI,TE,TB,RAT)
C...... Empirical Te model from Brace and Theis GRL p275, 1978
      DOUBLE PRECISION H,NI,TE,TB,RAT
      IF(H.LT.120.OR.H.GT.500) RETURN
      TB=1051+(17.07*H-2746)
     > *EXP(-5.122E-4*H+6.094E-6*NI-3.353E-8*H*NI)
      RAT=TE/TB
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CNEUT(J,I,JP,Z,O1S,N2A,N2D,O1D,NNO,N4S,N2P)
C.......neutral densities
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      IF(JP.EQ.1) WRITE(I,91)
 91   FORMAT(/4X,'ALT',3X,'[O1S]',4X,'[N2A]',4X,'[N2D]',4X,'[O1D]',4X
     > ,'[NO]',4X,'[N4S]',4X,'[N2P]')
      WRITE(I,7) Z,O1S,N2A,N2D,O1D,NNO,N4S,N2P
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CEMM(J,I,JP,Z,RTS,O1D,N2D,OP2D,OP2P,O1S,DISNP,NE,ON
     > ,N2N,N4S,P3X12,P3X11,N2A,P3X2,P3X3,DIS2S)
C.......... production rates
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99)
      IF(JP.EQ.1) WRITE(I,92)
 92   FORMAT(/1X,'Some important volume emission rates
     > (photons/cm3/sec)',/4X,'ALT',2X,'E7320'
     > ,4X,'E5577',4X,'E5755',4X,'E6584',4X,'E1493',4X
     > ,'VEGCAP',4X,'E1200',4X,'E3371',4X,'E6300',4X,'E5200',4X,'E2972'
     > ,4X,'E2470',4X,'E2143    E6364  E3726_29 ')
      E6300=O1D*RTS(54) * 0.76
      E6364=O1D*RTS(54) * 0.24
      E5200=N2D*RTS(61)
      E7320=OP2P*0.17
      E5577=O1S*RTS(55)
      E2470=OP2P*0.0477
      E2972=O1S*RTS(56)
      E2143=RTS(70)*DIS2S
      E5755=RTS(72)*DISNP
      E6584=RTS(71)*DISNP*4.0E-3/(4.0E-3+1.0E-7*NE)
      E1493=0.5*P3X12+0.1*P3X11*N4S+0.025*DISNP
      VEGCAP=N2A*0.57
      E1200=P3X12+P3X11*N4S+0.1*DISNP
      PN2B=P3X2
C........ The factor of 2 compensates for dissociation
      E3371=2.0*  0.25*P3X3
	E3726=7.7E-5*OP2D  !.. O+(2D) radiation
      WRITE(I,7) Z,E7320,E5577,E5755,E6584,E1493,VEGCAP
     > ,E1200,E3371,E6300,E5200,E2972,E2470,E2143,E6364,E3726
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CO2EM(J,I,JP,Z,O2B1,O2DEL,O2SS,O2ADEL,O2ASIG)
C........ Emissions from O2
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      IF(JP.EQ.1) WRITE(I,93)
 93   FORMAT(/5X,'Emissions from O+O+M -> O2'/4X,'ALT     O2atm   O2del
     >   HzII     Cham     HzI')
      WRITE(I,7) Z,(.083*O2B1),(2.83E-3*O2DEL),(5E-3*O2SS)
     > ,(0.02*O2ADEL),(6.25*O2ASIG)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN2D(J,I,JP,Z,RTS,ON,O2N,N2N,NOP,NE,P1,L1
     > ,N2PLS,DISN2D,UVDISN,NPLUS,N2P,N2D,OPLS,NNO,N2A)
C.......n(2d)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=NOP*NE*RTS(50)
      PR(2)=N2PLS*NE*RTS(32)*RTS(11)
      PR(3)=N2PLS*ON*RTS(10)
      PR(4)=DISN2D
      PR(5)=RTS(63)*UVDISN
      PR(6)=RTS(65)*NPLUS*O2N
      PR(7)=N2P*RTS(57)
      PR(8)=RTS(27)*N2A*ON
      LR(1)=ON*RTS(15)
      LR(2)=O2N*RTS(16)
      LR(3)=NE*RTS(8)
      LR(4)=OPLS*RTS(29)
      LR(5)=RTS(61)
      LR(6)=RTS(41)*NNO
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
      L1=LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)
C....... EF is used to convert production rates to volume emission rates
      EF=1.0
C....... This line in for volume emission rates
C...      EF=RTS(61)*0.76/L1
      IF(JP.EQ.1.AND.I.GT.0.AND.INT(EF+0.1).NE.1) WRITE(I,189)
      IF(JP.EQ.1.AND.I.GT.0.AND.INT(EF+0.1).EQ.1) WRITE(I,191)
 189  FORMAT(/2X,'N(2D)',25X,'EMISSION',28X,':',20X,'Loss rate')
 191  FORMAT(/2X,'N(2D)',25X,'Production',36X,':',20X,'Loss rate')
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,96)
 96   FORMAT(2X,'ALT   [N2D]   [N2D]c   NO++e   N2++e   N2++O    e+N2'
     >  ,3X,'hv+N2   N++O2   N(2P)   N2A+O    +O     +O2      +e'
     >  ,5X,'+O+     RAD     +NO')
      IF(I.GT.0) WRITE(I,7) Z,N2D,P1/L1,(PR(K)*EF,K=1,8)
     > ,(LR(K)*N2D,K=1,6)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CO1D(IEMS,J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1,L1
     > ,O2P,N2D,P1X1,PO1DSR,O1S,NPLUS,O1D)
C....... O(1D)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=O2P*NE*RTS(51)
      PR(2)=N2D*O2N*RTS(16)*RTS(74)
      PR(3)=P1X1
      PR(4)=PO1DSR
      PR(5)=RTS(18)*NE*ON
      PR(6)=O1S*RTS(55)
      PR(7)=NPLUS*O2N*RTS(66)
      PR(8)=O1S*O2N*0.31*RTS(48)    !... Fox 2001
      LR(1)=RTS(54)
      LR(2)=N2N*RTS(33)
      LR(3)=O2N*RTS(34)
      LR(4)=ON*RTS(69)
      LR(5)=RTS(98)*NE              !..Fox
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
      L1=LR(1)+LR(2)+LR(3)+LR(4)+LR(5)
C....... EF is used to convert production rates to volume emission rates
C....... EF=0 gives production rates else volume emission rates
      IF(IEMS.EQ.0) THEN
         EF=1.0
      ELSE
         EF=RTS(54)*0.76/L1
      ENDIF
      IF(JP.EQ.1.AND.I.GT.0.AND.INT(EF+0.1).NE.1) WRITE(I,189)
      IF(JP.EQ.1.AND.I.GT.0.AND.INT(EF+0.1).EQ.1) WRITE(I,191)
 189   FORMAT(/2X,'O(1D)',30X,'EMISSION',40X,':',7X,'LOSS RATES')
 191   FORMAT(/2X,'O(1D)',30X,'PRODUCTION',40X,':',7X,'LOSS RATES')
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,193)
 193   FORMAT(4X,'ALT',2X,'[O1D]',3X,'[O1D]c',2X,'O2++e',2X,'N2D+O2',3X
     > ,'e* +O    S-R    e+O    O1S->  N++O2  O1S+O2   Emiss    +N2'
     >  ,4X,'  +O2     +O     +e')
      IF(I.GT.0) WRITE(I,7) Z,O1D,P1/L1,(PR(K)*EF,K=1,8)
     > ,(LR(K)*O1D,K=1,5)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CNO(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1,L1
     > ,N2D,N4S,N2P,NNO,O2P,OPLS,PDNOSR,PLYNOP,N2A,NPLUS)
C........no
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=RTS(16)*O2N*N2D
      PR(2)=RTS(7)*O2N*N4S
      PR(3)=RTS(38)*N2P*O2N
	PR(4)=RTS(27)*N2A*ON
      PR(5)=RTS(22)*NPLUS*O2N          !.. Fox
      LR(1)=RTS(9)*N4S
      LR(2)=RTS(23)*O2P
      LR(3)=RTS(24)*OPLS
      LR(4)=RTS(41)*N2D
      LR(5)=PDNOSR
      LR(6)=PLYNOP
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)
      L1=LR(1)+LR(2)+LR(3) + LR(4) + (LR(5)+LR(6))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,192)
 192   FORMAT(/2X,'NO',17X,'PRODUCTION',20X,':',10X,'LOSS RATES'/
     > ,4X,'ALT',3X,'[NO]',5X,'[NO]c',3X,'O2+N2D',
     > 3X,'O2+N4S   N2P+O2   N2A+O    N++O2    N4S+NO   O2P+NO   O++NO'
     > ,3X,'N2D+NO   hv<1910   Lyman-a')
      IF(I.GT.0) WRITE(I,7) Z,NNO,P1/L1,(PR(K),K=1,5)
     > ,(LR(K)*NNO,K=1,6)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN4S(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1,L1,N4S,DISN4S,N2D
     >   ,N2P,OPLS,N2PLS,UVDISN,NOP,NPLUS,NNO,O2P,PDNOSR)
C........N(4S)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=DISN4S
      PR(2)=RTS(15)*ON*N2D
      PR(3)=RTS(8)*NE*N2D
      PR(4)=RTS(3)*OPLS*N2N
      PR(5)=RTS(53)*RTS(11)*N2PLS*NE
      PR(6)=RTS(62)*UVDISN
      PR(7)=NOP*NE*RTS(49)
      PR(8)=N2D*RTS(61)
      PR(9)=N2P*RTS(58)
      PR(10)=RTS(25)*NPLUS*O2N
      PR(11)=PDNOSR*NNO
      PR(12)=NPLUS*NNO*RTS(81)      !..Fox
      LR(1)=RTS(7)*O2N
      LR(2)=RTS(9)*NNO
      LR(3)=RTS(21)*O2P
      LR(4)=RTS(79)*N2PLS     !..Fox
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9)+PR(10)
     >     +PR(11)+PR(12)
      L1=LR(1)+LR(2)+LR(3)+LR(4)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,193)
 193   FORMAT(/2X,'N(4S)',38X,'PRODUCTION',46X,':',7X,'LOSS RATES'/
     > ,4X,'ALT',2X,'[N4S]',2X,'hv->N+'
     > ,3X,'O+N2D',2X,'e+N2D',3X,'O++N2',3X,'N2++e',4X,'hv->2N'
     > ,2X,'NO++e',2X,'N(2D)',4X,'N(2P)   N+&X    hv+NO    +O2  '
     > ,2X,' +NO  ',2X,'+O2+ & N2+')
      PR(10)=PR(10)+PR(12)  !... for printing fit
      LR(3)=LR(3)+LR(4)     !... for printing fit
      IF(I.GT.0) WRITE(I,7) Z,N4S,(PR(K),K=1,11),(LR(K)*N4S,K=1,3)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CO1S(J,I,JP,Z,RTS,ON,O2N,N2N,NE
     > ,O1S,O2P,N4S,P1X2,NPLUS,N2A,PO1DSR,O2ADEL,O2DEL)
C.......o(1s)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=O2P*N4S*RTS(35)
      PR(2)=O2P*NE*RTS(52)
      PR(3)=P1X2
      PR(4)=RTS(59)*NPLUS*O2N
      PR(5)=RTS(60)*RTS(36)*ON*N2A
      PR(6)=RTS(67)*PO1DSR
      PR(7)=RTS(68)*3E-11*ON*O2ADEL
      LR(1)=RTS(55)
      LR(2)=RTS(48)*O2N
      LR(3)=RTS(56)
      LR(4)=1.7E-10*O2DEL
      LR(5)=(RTS(46)+RTS(47))*NE
      O1S=(PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7))/
     > (LR(1)+LR(2)+LR(3)+LR(4))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,194)
 194   FORMAT(/2X,'O(1S)',27X,'PRODUCTION',34X,':',14X,'LOSS RATES'/
     > ,4X,'ALT',3X,'[O1S]',4X,'N+O2+',4X,'e+O2+',4X,'e+O',4X
     > ,'N++O2',4X,'N2A+O',4X,'S-R',6X,'O2ADEL',5X,'5577'
     > ,4X,'+ O2',4X,'2972',5X,'O2DEL')
      IF(I.GT.0) WRITE(I,7) Z,O1S,(PR(K),K=1,7),(LR(K)*O1S,K=1,4)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::: CN2PLS :::::::::::::::::::::::::::::::
C..... Simplified chemistry of N2+.  PUN2P* = production of N2+ by euv 
C..... in the (X,A,B states). PEN2P* same for p.e.s (X,A,B states)
      SUBROUTINE CN2PLS(J,I,JP,Z,RTS,ON,O2N,N2N,NE,N2PLS,PUN2PX,PUN2PA
     >  ,PUN2PB,PEN2PX,PEN2PA,PEN2PB,OP2D,OP2P,HEPLUS,NPLUS,NNO,N4S)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      REAL PUN2PX,PUN2PA,PUN2PB,PEN2PX,PEN2PA,PEN2PB
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=PUN2PX
      PR(2)=PUN2PA
      PR(3)=PUN2PB
      PR(4)=PEN2PX
      PR(5)=PEN2PA
      PR(6)=PEN2PB
      PR(7)=RTS(19)*OP2D*N2N
      PR(8)=RTS(20)*OP2P*N2N
      PR(9)=RTS(44)*HEPLUS*N2N
      PR(10)=RTS(82)*NPLUS*NNO   !..Fox
      LR(1)=RTS(10)*ON
      LR(2)=RTS(11)*NE
      LR(3)=RTS(17)*O2N
      LR(4)=RTS(99)*ON      !.. RTS(99) from RVN2PB file for N2+(v)
      LR(5)=RTS(79)*N4S     !.. Fox
      LR(6)=RTS(80)*NNO     !.. Fox
      IF(I.LE.0)  N2PLS=
     >   (PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9)+PR(10))/
     >     (LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,95)
 95   FORMAT(/2X,'N2+',29X,'PRODUCTION',45X,':',12X,'LOSS RATES'/
     > ,4X,'ALT  [N2+]  EUV-X   EUV-A    EUV-B   PE-X'
     > ,5X,'PE-A    PE-B  O+2D+N2  O+2P+N2  He++N2  O+N2+'
     > ,2X,'e+N2+  O2+N2+  N2++O  Other')
      PR(9)=PR(9)+PR(10)     !.. for printing fit
      LR(5)=LR(5)+LR(6)      !.. for printing fit
      IF(I.GT.0) WRITE(I,7) Z,N2PLS,(PR(K),K=1,9)
     > ,(LR(K)*N2PLS,K=1,5)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C:::::::::::::::::::::::::::::: CNOPO ::::::::::::::::::::::::::::::::::
      SUBROUTINE CNOPO(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1,NOP,OPLS
     >  ,N2PLS,O2P,N4S,NNO,NPLUS,N2P,PLYNOP,N2D,OP2D,mp,lp)
C........no+
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=RTS(3)*N2N*OPLS
      PR(2)=N2PLS*ON*RTS(10)
      PR(3)=O2P*N4S*RTS(21)
      PR(4)=O2P*NNO*RTS(23)
      !.. N+ + O2 -> O2+ + N(2D,4S) or NO+ + O(1S)
      PR(5)=(RTS(30)+RTS(66)+RTS(59))*NPLUS*O2N
      PR(6)=RTS(37)*N2P*ON
      PR(7)=RTS(24)*OPLS*NNO
      PR(8)=PLYNOP*NNO
      PR(9)=O2P*N2D*RTS(77)         !..Fox
      PR(10)=N2PLS*NNO*RTS(80)      !..Fox
      PR(11)=NPLUS*NNO*RTS(81)      !..Fox
      PR(12)=RTS(83)*NNO*OP2D       !..Fox
      LR(1)=NE*RTS(5)
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
     >    +PR(9)+PR(10)+PR(11)+PR(12)
      IF(I.LE.0) NOP=P1/LR(1)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,96)
 96   FORMAT(/2X,'NO+',31X,'PRODUCTION',48X,':',2X,'LOSS RATES'/
     > ,4X,'ALT',3X,'[NO+]',4X,'O++N2',3X,'N2++O',3X,'O2++N4S'
     > ,3X,'O2++NO',3X,'N++O2',4X,'N2P+O',3X,'O++NO   hv+NO'
     > ,6X,'Other   NO++e')
      PR(9)=PR(9)+PR(10)+PR(11)+PR(12)
      IF(I.GT.0) WRITE(I,7) Z,NOP,(PR(K),K=1,9),LR(1)*NOP
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CO2P(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1
     > ,O2P,TPROD5,OPLS,OP2D,N2PLS,NPLUS,N4S,NNO,OP2P,
     > mp,ilp)
C........o2+
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
C......... TPROD5=euv @ p.e. production
      PR(1)=TPROD5
      PR(2)=RTS(4)*O2N*OPLS
      PR(3)=RTS(43)*OP2D*O2N
      PR(4)=RTS(17)*O2N*N2PLS
      PR(5)=RTS(25)*NPLUS*O2N
      PR(6)=RTS(86)*OP2P*O2N           !.. Fox
      PR(7)=RTS(65)*NPLUS*O2N
      LR(1)=RTS(6)*NE
      LR(2)=RTS(21)*N4S
      LR(3)=RTS(23)*NNO
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)
      IF(I.LE.0) O2P=P1/(LR(1)+LR(2)+LR(3))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,97)
  97  FORMAT(/2X,'O2+',22X,'PRODUCTION',48X,':',9X,'LOSS RATES'
     > /,4X,'ALT   [O2+]   hv+O2   O++O2   O+(2D)+O2   N2++O2   N++O2'
     >  ,3X,'O+(2P)+O2  N++O2    O2++e    O2++N    O2++NO')
      IF(I.GT.0) WRITE(I,7) Z,O2P,(PR(K),K=1,7),(LR(K)*O2P,K=1,3)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COP4S(J,I,JP,Z,RTS,ON,O2N,N2N,NE,OPLS,TPROD1,OP2D
     >  ,OP2P,PEPION,PDISOP,N2PLS,N2D,NNO,HEPLUS,NPLUS)
C...........o+(4s)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
C.........pr(1)= euv production of o+(4s)
      PR(1)=TPROD1
      PR(2)=OP2D*NE*RTS(12)
      PR(3)=OP2P*ON*RTS(26)
      PR(4)=PEPION
      PR(5)=PDISOP
      PR(6)=RTS(99)*N2PLS*ON      !.. RTS(99) from RVN2PB file for N2+(v)
      PR(7)=OP2P*NE*RTS(14)
      PR(8)=OP2P*0.047
      PR(9)=RTS(28)*ON*OP2D
      PR(10)=RTS(85)*OP2P*O2N              !.. Fox
      PR(11)=HEPLUS*O2N*(RTS(91)+RTS(93))  !..Fox
      PR(12)=RTS(95)*NNO*HEPLUS            !..Fox
      PR(13)=RTS(22)*NPLUS*O2N             !..Fox
      LR(1)=N2N*RTS(3)
      LR(2)=O2N*RTS(4)
      LR(3)=NNO*RTS(24)
      LR(4)=N2D*RTS(29)    !.... small loss?? ..Fox
      !..LR(4)=(LR(1)+LR(2)+LR(3))  !.. total loss for printing
      PR(10)=PR(10)+PR(11)+PR(12)+PR(13)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,98)
 98   FORMAT(/2X,'O+',41X,'PRODUCTION',39X,':',10X,'LOSS RATES'/
     >  ,' ALT    [O+]     [O+]c   hv+O  O+(2D)+e O+(2P)+O   e*+O   '
     >  ,'O2-diss   N2++O  O+(2P)+e O+(2P) O+O+(2D)   Other   +N2'
     >  ,5X,'+O2    +NO    +N2D')
      OPLUSC=(PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+
     >    PR(9)+PR(10))/(LR(1)+LR(2)+LR(3)+LR(4)+1.0E-10)
      IF(I.GT.0) WRITE(I,7) Z,OPLS,OPLUSC,(PR(K),K=1,10)
     > ,(LR(K)*OPLS,K=1,4)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COP2D(J,I,JP,Z,RTS,ON,O2N,N2N,NE,OP2D,PPROD,
     > EPROD,OP2P,HEPLUS,N4S,NNO)
C.......o+(2d)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
	REAL PPROD,EPROD
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=PPROD                  ! EUV prod
      PR(2)=OP2P*NE*RTS(13)
      PR(3)=OP2P*0.171
      PR(4)=HEPLUS*O2N*RTS(76)     !..Fox
	PR(5)=EPROD                  !.. Photoelectron prod
      LR(1)=RTS(19)*N2N
      LR(2)=7.7E-5                 !.. radiation at 
      LR(3)=NE*RTS(12)
      LR(4)=ON*RTS(28)
      LR(5)=RTS(43)*O2N
      LR(6)=RTS(83)*NNO    !..Fox
      LR(7)=RTS(84)*N4S    !..Fox
      OP2D=(PR(1)+PR(2)+PR(3)+PR(4)+PR(5))/
     >   (LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)+LR(7))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,99)
 99   FORMAT(/2X,'O+(2D)',12X,'PRODUCTION',22X,':',18X,'LOSS RATES'/
     > ,4X,'ALT',3X,'[O+2D]',3X,'hv+O',4X,'Phote',4X,'O+2P+e',3X,
     > 'O+2P>hv',2X,'He++O2     +N2   E3726_29    +e       +O      +O2'
     > '      +NO     +N')
      IF(I.GT.0) WRITE(I,7) Z,OP2D,PR(1),PR(5),PR(2),PR(3),PR(4)
     > ,(LR(K)*OP2D,K=1,7)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COP2P(J,I,JP,Z,RTS,ON,O2N,N2N,NE,OP2P,TPROD3,PSEC
     > ,HEPLUS,N4S,NNO)
C.......o+(2p)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=0.0
      IF(TPROD3.GE.PSEC) PR(1)=TPROD3-PSEC
      PR(2)=PSEC
      PR(3)=HEPLUS*O2N*RTS(92)        !..Fox
      LR(1)=RTS(26)*ON
      LR(2)=RTS(20)*N2N
      LR(3)=RTS(13)*NE
      LR(4)=0.218
      LR(5)=RTS(14)*NE
      LR(6)=(RTS(85)+RTS(86))*O2N    !..Fox
      LR(7)=RTS(87)*N4S              !..Fox
      LR(8)=RTS(88)*NNO              !..Fox
      OP2P=(TPROD3+PR(3))
     >  /(LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)+LR(7)+LR(8))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,100)
 100   FORMAT(/2X,' O+(2P)',6X,'PRODUCTION',10X,':',12X,'LOSS RATES'/
     >,4X,'ALT   [O+2P]    hv+O     e*+O  He++O2      +O',7X,'+N2'
     > ,6x,'+e       RAD      +e      +O2      +N4S     +NO')
      IF(I.GT.0) WRITE(I,7) Z,OP2P,PR(1),PR(2),PR(3),(LR(K)*OP2P,K=1,8)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C......N+
      SUBROUTINE CNPLS(J,I,JP,INPLS,Z,RTS,ON,O2N,N2N,NE,DISNP
     > ,NPLUS,OPLS,N2D,OP2P,HEPLUS,PHOTN,O2P,N4S,OP2D,N2PLS,NNO)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=DISNP
      PR(2)=RTS(29)*OPLS*N2D
      PR(3)=0
      PR(4)=RTS(45)*HEPLUS*N2N
      PR(5)=PHOTN
      PR(6)=O2P*N2D*RTS(78)          !..Fox
      PR(7)=N2PLS*N4S*RTS(79)        !..Fox
      PR(8)=OP2D*N4S*RTS(84)         !..Fox
      PR(9)=RTS(94)*NNO*HEPLUS       !..Fox
      LR(1)=RTS(30)*O2N              !..Fox
      LR(2)=RTS(25)*O2N              !..Fox
      LR(3)=RTS(22)*O2N              !..Fox
      LR(4)=RTS(65)*O2N              !..Fox
      LR(5)=RTS(66)*O2N              !..Fox
      LR(6)=RTS(31)*ON               !..Fox

      CNPLUS=(PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)+PR(9))
     >   /(LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6))
      IF(INPLS.LE.0) NPLUS=CNPLUS
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,101)
 101   FORMAT(/2X,'N+',20X,'PRODUCTION',73X,':',8X,'LOSS RATES'/
     > ,4X,'ALT   [N+]   [N+]c     hv+N2   O++N2D  O+2P+N2',3X
     > ,'He++N2',3X,' hv+N   O2++N2D  N2++N4S O+(2D)+N4S  He++NO'
     > ,3X,'N++O2    N++O2    N++O2    N++O2    N++O2    N++O')
      IF(I.GT.0) WRITE(I,7) Z,NPLUS,CNPLUS
     > ,(PR(K),K=1,9),(LR(K)*NPLUS,K=1,6)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN2A(J,I,JP,Z,RTS,ON,O2N,N2N,NE
     > ,N2A,P3X1,P3X2,P3X3,P3X4)
C........n2(a3sigma) and total LBH
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      REAL P3X1,P3X2,P3X3,P3X4
      DIMENSION RTS(99),LR(22),PR(22)
C....... pr(1,2,3)= electron impact excitation of n2(a,b,c) states
      PR(1)=P3X1
      PR(2)=P3X2
      PR(3)=P3X3
      LR(1)=RTS(36)*ON
      LR(2)=RTS(27)*ON
      LR(3)=0.57
      N2A=(PR(1)+PR(2)+PR(3))/(LR(1)+LR(2)+LR(3))
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,102)
 102   FORMAT(/2X,'N2(A)',12X,'PRODUCTION',13X,':',5X,'LOSS RATES'
     > ,3X,':  Total LBH'
     > /,4X,'ALT',3X,'N2(A)',3X,'e*->N2A',3X,'e*->N2B',3X,'e*->N2C',2X
     >  ,'N2A>O1S',2X,'N2A>NO',2X,'RAD',5X,'LBH')
      IF(I.GT.0) WRITE(I,7) Z,N2A,(PR(K),K=1,3),(LR(K)*N2A,K=1,3)
     > ,P3X4
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CN2P(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1,L1
     > ,N2P,P3X7,UVDISN,O2P,NNO,N2PLUS)
C....... N(2P). the rates are from Zipf et al JGR, 1980 p687
C.... 21-AUG-1992. Added N2+ recombination source
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=P3X7
      PR(2)=RTS(64)*UVDISN
      PR(3)=RTS(73)*RTS(11)*N2PLUS*NE
      LR(1)=RTS(37)*ON
      LR(2)=RTS(38)*O2N
      LR(3)=RTS(39)*O2P
      LR(4)=RTS(40)*NNO
      LR(5)=RTS(57)
      LR(6)=RTS(58)
      LR(7)=(RTS(96)+RTS(97))*NE    !..Fox
      P1=PR(1)+PR(2)+PR(3)
      L1=LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6)+LR(7)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,103)
 103   FORMAT(/2X,'N(2P)',9X,'PRODUCTION',17X,':',20X,'LOSS RATES'/
     > ,4X,'ALT',3X,'[N2P]',3X,'e+N2',5X,'hv+N2',3X,'e+N2+',6X,'+O  '
     > ,3X,'+O2      +O2+      +NO       +2D     +4S      +e')
      IF(I.GT.0) WRITE(I,7) Z,N2P,(PR(K),K=1,3),(LR(K)*N2P,K=1,7)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C....... He+
      SUBROUTINE CHEP(J,I,JP,IHEPLS,Z,RTS,ON,O2N,N2N,NE,HEPLUS,PRHEP,
     >   HE,NNO ,OPLUS,HPLUS,HN)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
      PR(1)=PRHEP
      PR(2)=HE
      LR(1)=RTS(44)*N2N
      LR(2)=RTS(45)*N2N
      LR(3)=RTS(75)*O2N
      LR(4)=(RTS(76)+RTS(91)+RTS(92)+RTS(93))*O2N    !..Fox
      LR(5)=RTS(94)*NNO                              !..Fox
      LR(6)=RTS(95)*NNO                              !..Fox
      !...... He+ from chemical equilibrium
      CHEPLUS=PR(1)/(LR(1)+LR(2)+LR(3)+LR(4)+LR(5)+LR(6))
      IF(IHEPLS.LE.0) HEPLUS=CHEPLUS

      !... H+ production and loss
      PHPLUS=RTS(2)*OPLUS*HN
      LHPLUS=RTS(1)*ON*HPLUS
      CHPLUS=PHPLUS/(RTS(2)*ON)

      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,104)
 104   FORMAT(/2X,'He+',9X,'PRODUCTION',9X,':',27X,'LOSS RATES'
     > ,17X,':',10X,'H+ PROD and LOSS rates'
     > ,/4X,'ALT    [He+]   [He+]c   hv+He     +N2    +N2->N+'
     > ,3X,'+O2   +O2->O+      +NO      +NO'
     > ,3X,'<<>>',4X,'[H+]   [H+]c   O+ + H   H+ + O')
      IF(I.GT.0) WRITE(I,7) Z,HEPLUS,CHEPLUS,PR(1)
     > ,(LR(K)*HEPLUS,K=1,6),HPLUS,CHPLUS,PHPLUS,LHPLUS
      RETURN
 7    FORMAT(F7.2,1P,9E9.2,5X,1P,9E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CRT1(J,I,JP,Z,RTS)
C............ altitude dependent reaction rates .........
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,105)
 105  FORMAT(/2X,'FIRST 13 NON-CONSTANT REACTION RATES FROM SUBROUTINE
     > RATS'/,4X,'ALT',3X,'O+H+',3X,'O++H',4X,'O++N2',3X,'O++O2',
     > 3X,'NO++e',3X,'O2++e',3X,'O2+N4S',2X,'N2D+e',3X,
     > 'NO+N4S',2X,'N2++O',2X,'N2++e',2X,'O+2D+e',2X,'O+2P+e')
      IF(I.GT.0) WRITE(I,7) Z,(RTS(IS),IS=1,13)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CRT2(J,I,JP,Z,RTS)
C............ altitude dependent reaction rates .........
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,106)
 106  FORMAT(/2X,'MORE NON-CONSTANT RATES FROM SUBROUTINE RATS'/
     >,4X,'ALT',2X,'O+2P+e',2X,'N2++O2',2X,'O1D+e'
     >,2X,'O++NO',2X,'N2P+O',2X,'N2++O',2X,'N2+(V)+O')
      IF(I.GT.0) WRITE(I,7) Z,RTS(14),RTS(17),RTS(18),RTS(24)
     > ,RTS(37),RTS(42)*RTS(10),RTS(99)
      RETURN
 7    FORMAT(F7.2,1P,22E8.1)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE C1200(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P3X12,P3X11,N4S,DISNP)
C....... NI 1200 A and 1493 A. Production from e* + N2 (P3X12), e* + N
C....... (P3X11) and hv + N2 -> N+ + N + e. Refs - Wu etal JGR 1983, 2163
C....... Doering and Gumboel, JGR 1993, 16021, Beyer .. J.Ch.Phys. 1969, 5323
C....... FLIP agrees with 1493 from Meier et al. 1980, JGR p2177
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),PR(22)
      PR(1)=P3X12
      PR(2)=P3X11*N4S
      PR(3)=0.1*DISNP
C...... 1493. 0.5 and 0.1 multipliers in PR(4) and PR(5) from X-sec ratios
      PR(4)=0.5*P3X12
      PR(5)=0.1*P3X11*N4S
      PR(6)=0.025*DISNP
      P1200=PR(1)+PR(2)+PR(3)
      P1493=PR(4)+PR(5)+PR(6)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,108)
 108  FORMAT(/9X,'NI 1200 PRODUCTION RATES',8X,': 1493 PRODUCTION RATES'
     > /,'  ALT   T1200    e*+N2     e*+N    hv+N2    T1493    e*+N2
     >    e*+N     hv+N2')
      IF(I.GT.0) WRITE(I,7) Z,P1200,(PR(K),K=1,3),P1493,(PR(K),K=4,6)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CO2(J,I,JP,Z,RTS,ON,O2N,N2N,NE,BOD3,O1D,O2SS,TN
     > ,O2B1,O2DEL,O2ADEL,O2ASIG)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99),LR(22),PR(22)
C...... Avoid underflow at high altitudes
      IF(Z.GT.700.0) RETURN
C...... densities for O2 atmospheric bands
      BOD3=(300/TN)**2*ON**2*(O2N+N2N)*4.7E-33
      O2SS=(0.8*BOD3)/(5E-13*O2N+3E-11*ON+5E-3)
      O2B1=(5E-13*O2SS*O2N+RTS(34)*O1D*O2N)/
     > (0.0*ON+2.2E-15*N2N+0.083)
      O2DEL=(0.25*BOD3)/
     > (1.3E-16*ON+2.6E-20*TN**0.78*O2N+1E-20*N2N+2.83E-3)
      O2ADEL=(0.05*BOD3)/(0.02+2E-13*N2N)
      O2ASIG=0.05*BOD3/(3E-13*N2N+3E-13*ON+6.25)
C
      PR(1)=BOD3
      PR(2)=RTS(34)*O2N*O1D
      PR(3)=5E-13*O2SS*O2N
      PR(4)=0.25*BOD3
      PR(5)=0.05*BOD3
      PR(6)=0.05*BOD3
      LR(1)=(5E-13*O2N+3E-11*ON+5E-3)
      LR(2)=(0.0*ON+2.2E-15*N2N+0.083)
      LR(3)=1.3E-16*ON+2.6E-20*TN**0.78*O2N+1E-20*N2N+2.83E-3
      LR(4)=(0.02+2E-13*N2N)
      LR(5)=(3E-13*N2N+3E-13*ON+6.25)
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,109)
 109  FORMAT(/15X,'Production rates of O2',21X,':',12X,'Total loss rates
     > of O2 excited states'/,4X,'ALT',3X
     > ,'O+O+M   O1D+O2   O2**+O2    O2DEL   O2ADEL   O2ASIG  :   O2**
     >   O2b1   O2del   O2ADEL   O2ASIG')
      IF(I.GT.0) WRITE(I,7) Z,(PR(K),K=1,6),(LR(K),K=1,5)
      RETURN
 7    FORMAT(F7.2,1P,22E9.2)
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PRDPRA(J,I,JP,ALT)
      DOUBLE PRECISION ALT
      IF (JP.EQ.1.AND.I.GT.0) WRITE(I,110)
 110  FORMAT(/2X,'*** O+, N2+, O2+ Secondary to primary production'
     > ,1X,'ratios not implemented on this file, see Field line grid'
     > ,1X,'files instead ****')
      RETURN
      END
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PRDPRB(J,I,JP,ALT)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      DOUBLE PRECISION ALT
      DIMENSION UVP(3,12),RAT(3,12)
      IF (JP.EQ.1.AND.I.GT.0) WRITE(I,110)
 110  FORMAT(/22X,'Secondary to primary production ratios'
     > ,/4X,'ALT  O+(4S)   O+(2D)   O+(2P)   N2+(X)   N2+(A)   N2+(B) 
     > EOX/UOX  EO2/UO2  EN2/UN2  U2S/UN2+  UN+/UN2  UO+/UO2  UVOX
     >  UVO2 UVN2')
      DO 20 IS=1,3
      DO 20 IK=1,12
      RAT(IS,IK)=0.0
      UVP(IS,IK)=EUVION(IS,IK,J)
      IF(UVP(IS,IK).LE.0.0) UVP(IS,IK)=0.0
      IF(UVP(IS,IK).LE.0) GO TO 20
      RAT(IS,IK)=PEPION(IS,IK,J)/UVP(IS,IK)
 20   CONTINUE
C........ Take into account difference between Kirby and Jackman
      IF(UVP(2,2).GT.0) RAT(2,2)=(PEPION(2,2,J)+PEPION(2,3,J))/UVP(2,2)
      IF(UVP(2,3).GT.0) RAT(2,3)=PEPION(2,4,J)/UVP(2,3)
C...... Calculate total production rates ....
      UVOX=UVP(1,1)+UVP(1,2)+UVP(1,3)+UVP(1,4)+UVP(1,5)
      PEOX=PEPION(1,1,J)+PEPION(1,2,J)+PEPION(1,3,J)+PEPION(1,4,J)
     > +PEPION(1,5,J)
      UVO2=UVP(2,1)+UVP(2,2)+UVP(2,3)+UVP(2,4)
      PEO2=PEPION(2,1,J)+PEPION(2,2,J)+PEPION(2,3,J)+PEPION(2,4,J)
      UVN2=UVP(3,1)+UVP(3,2)+UVP(3,3)
      PEN2=PEPION(3,1,J)+PEPION(3,2,J)+PEPION(3,3,J)
      UVNP=UVP(3,4)+UVP(3,5)+UVP(3,6)
      PENP=PEPION(3,4,J)+PEPION(3,5,J)+PEPION(3,6,J)
C....... calculate production ratios
      RATOX=0.0
      RATO2=0.0
      RATN2=0.0
      RAT2S=0.0
      RNP=0.0
      IF(UVOX.GT.0) RATOX=PEOX/UVOX
      IF(UVO2.GT.0) RATO2=PEO2/UVO2
      IF(UVN2.GT.0) RATN2=PEN2/UVN2
      IF(UVO2.GT.0) RNOP=(UVP(2,4)+UVP(2,5)+UVP(2,6))/UVO2
      IF(UVN2.GT.0) RAT2S=(UVP(3,6))/UVN2
      IF(UVN2.GT.0) RNP=UVNP/UVN2
      TPOX=PEOX+UVOX*1.0
      TPO2=PEO2+UVO2*1.0
      TPN2=PEN2+UVN2*1.0
      WRITE(I,90) ALT,RAT(1,1),RAT(1,2),RAT(1,3),RAT(3,1)
     > ,RAT(3,2),RAT(3,3),RATOX,RATO2,
     > RATN2,RAT2S,RNP,RNOP,UVOX,UVO2
C.....     > TPOX,TPO2,TPN2,UVOX,UVO2,UVN2
 90   FORMAT(F7.2,1P,22E9.2)
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CEMOX(J,I,JP,Z,RTS,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,NE,
     >   OPLUS)
C.......... Electron impact excitation rates for atomic oxygen
C...... P1...P6 are the excitation rates for 1D, 1S, 1304, 989,
C...... 1027, 1356. P7-> not used
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION RTS(99)
      IF(JP.EQ.1) WRITE(I,92)
 92   FORMAT(/2X,'Photoelectron excitation rates for atomic oxygen
     > (photons/cm3/sec): 8446 = (0.5*1304 + 1027), 7774 = 0.5*1356'
     > ,/4X,'ALT',2X,'O(1D)'
     > ,4X,'O(1S)',4X,'T1304',4X,'T989',4X,' T1027',4X,'T1356'
     > ,4X,'E8446',4X,'E7774  e+O+ >7774 He10830')
C...... The 8446 excitation rate is estimated using the results
C...... of Zipf and Erdman JGR, page 11089, 1985. in which they
C...... contend that the direct excitation of 1304 and 8446 are
C...... approximately equal. Since 8446 cascades to 1304, the 8446
C...... excitation rate is 0.5 of 1304. According to Meier Space Sci.
C...... Rev. 1991, the 1027 is trapped and cascades through 8446
      E8446=0.5*P3 + P5
C..... The 7774 emission is estimated in a similar manner to 8446
C..... The measurements of Gulcicek and Doering, JGR, p5885, 1988 
C..... indicate that 20% of the observed 1356 emission is due to direct
C..... excitation of 7774 followed by cascade. Most of the other quintet
C..... states cascade through 7774 to 1356 so a value of 0.5 has been 
C..... adopted for the proportion of the 1356 cascading via 7774
      E7774=0.5*P6
      WRITE(I,7) Z,P1,P2,P3,P4,P5,P6,E8446,E7774,7.3E-13*NE*OPLUS,P10
      RETURN
 7    FORMAT(F6.1,1P,22E9.2)
      END

C:::::::::::::::::::::::::::::::::: HOTOXY ::::::::::::::::::::::::::::
        SUBROUTINE HOTOXY (IFILE, J, I, JP, Z, RTS, NOP, O2P, ON, N4S,
     >                     N2N, O2N, ZNE, N2P, N2D, O1D, NNO, NPLUS,
     >                     OP2D, OP2P, NR, VN2O, HN, TI,HEPLUS)
C...Calculate production rates for non-thermal (hot) atomic oxygen as
C...well as exothermicity production rates. Calculate normalized energy
C...spectrum for every altitude.
        IMPLICIT NONE
        INTEGER  IFILE, J, I, JP, JEXO, IDUM, NBIN(29), MAXBIN,VD
	  PARAMETER (VD=181)    !.. VD = dimension of vertical grid
        REAL*8   PHOTO(28), POEXO(28), EXODAT(28), Z, RTS(99), NOP(VD), 
     >           O2P(VD), ON(VD), N4S(VD), N2N(VD), O2N(VD), ZNE(VD), 
     >           N2P(VD), N2D(VD), O1D(VD), NNO(VD), NPLUS(VD), OP2D, 
     >           OP2P, NR(2,VD),VN2O(6),HN(VD),TI(VD), PARTITION(28),
     >           XENERGY(60), SPECAMP(60), SPECRES, TOTHOTO,
     >           TOTAMP,HEPLUS
        DATA  EXODAT / 0.38, 2.75,  6.97, 5.02, 3.05, 2.38, 1.96, 5.00,
     >                 3.31, 1.46, 1.55, 1.33,  1.96, 3.76, 3.58, 3.25,
     >                 1.385, 6.67, 4.865, 3.02, 1.96, 4.20, 5.63, 0.30,
     >                 0.60, 0.90, 1.200, 1.00 /
C  The following PARTITION array accounts for conservation of momentum
C  and provides fraction of exothermic energy available to hot O.
	DATA  PARTITION / 0.467, 0.467, 0.500, 0.500, 0.500, 0.467,
     >                    0.500, 0.500, 0.500, 0.467, 0.667, 0.636,
     >                    0.636, 0.652, 0.467, 0.636, 0.652, 0.652,
     >                    0.667, 0.636, 0.667, 0.652, 0.636,
     >                    4*0.636, 0.056 /

        IF (JP.EQ.1) THEN
          WRITE(IFILE,1)
          IF(I.EQ.40) THEN
             WRITE(IFILE,*) ' Non-Thermal O Production Rates (large?)'
          ELSE IF(I.EQ.41) THEN 
             WRITE(IFILE,*) ' Non-Thermal O Production Rates (mod. ?)'
          ELSE IF(I.EQ.42) THEN 
             WRITE(IFILE,*) ' Kinetic Energy Production Rates (large?)'
          ELSE IF(I.EQ.43) THEN
             WRITE(IFILE,*) ' Kinetic Energy Production Rates (mod. ?)'
          ENDIF
          IF(I.EQ.40.OR.I.EQ.42) WRITE(IFILE,2) 
          IF(I.EQ.41.OR.I.EQ.43) WRITE(IFILE,3) 
          IF(I.EQ.44) WRITE(IFILE,31) 
        END IF
 1      FORMAT (/)
 2      FORMAT(' ALT   E+NO+    E+NO+    E+O2+    E+O2+    E+O2+',
     >   '    O+N(2D)  O+O(1D) O+O+(2P) O+O+(2D) N(2D)+O+  O2+O+',
     >   '  N2+O+(2D) O(1D)+N2 O2+N(2D)')
 3      FORMAT (' ALT   O+N(2P)  N+NO     N+O2     O2+N+'
     >   ,'   O2+O+(2D) N2+O+(2P) O2+O(1D) N+O2+    NO+N(2D)'
     >   ,' N2*(v=1,4) + O',22X,'O2+He+')
 31      FORMAT(4X,' Energy spectrum of hot O'
     >        ,/2X,'ENERGY   ALTIT.    AMP.')
 4      FORMAT (F6.1,1P,22E9.2)
 5      FORMAT (F6.1,1P,22E9.2)
 6	FORMAT (1P,E9.2,2X,0P,F6.1,2X,1P,E9.2)

        IF (I.EQ.40.OR.I.EQ.42.OR.I.EQ.44)  THEN
C...GROUP I CALCULATIONS:
C...The first 11 reactions can potentially produce large hot O populations
C...NO+ + e-  --> N(2D) + O + 0.38eV   (f=0.78)
          PHOTO(1) = RTS(5)*NOP(J)*ZNE(J)*0.78
C...NO+ + e-  --> N + O + 2.75eV  (f=0.22)
          PHOTO(2) = RTS(5)*NOP(J)*ZNE(J)*0.22

C...We neglect O(1S) production due to O2+ DR (=4%).
C...2 O atoms are produced in each reaction, hence the factor of 2
C...O2+ + e-  --> O + O + 6.97eV, f=33%
          PHOTO(3) = RTS(6)*O2P(J)*ZNE(J)*0.33*2.0
C...O2+ + e-  --> O + O(1D) + 5.02eV, f=21%
          PHOTO(4) = RTS(6)*O2P(J)*ZNE(J)*0.21*2.0
C...O2+ + e-  --> O(1D) + O(1D) + 6.97eV, f=42%
          PHOTO(5) = RTS(6)*O2P(J)*ZNE(J)*0.42*2.0

C...N(2D) + O  --> N + O + 2.38eV
          PHOTO(6) = RTS(15)*N2D(J)*ON(J)
C...O(1D) + O  --> O + O + 1.96eV    Note: 2 O atoms produced
          PHOTO(7) = RTS(69)*O1D(J)*ON(J)*2.0
C...O+(2P) + O  --> O + O + 5eV
          PHOTO(8) = RTS(26)*OP2P*ON(J)

C...O+(2D) + O  --> O+ + O + 3.31eV
          PHOTO(9) = RTS(28)*OP2D*ON(J)
C...O+ + N(2D)  --> O + N+ + 1.46eV
          PHOTO(10) = RTS(29)*NR(1,J)*N2D(J)

C...O+ + O2  --> O2+ + O + 1.55eV
          PHOTO(11) = RTS(4)*NR(1,J)*O2N(J)
C...O+(2D) + N2  --> N2+ + O + 1.33eV
          PHOTO(12) = RTS(19)*OP2D*N2N(J)

C...O(1D) + N2  --> O + N2 + 1.96eV
          PHOTO(13) = RTS(33)*O1D(J)*N2N(J)

C...The following hot O is not created chemically. It is now deleted
C...from the calculations, but left commented here for reference only.
C... O+ + H  --> O + H + kTi eV. Only if Ti > 3000; partitioning=0.059
C...energy added for O+ + H reaction 8.63E-5 is Boltzman k in eV units.
C...          PHOTO(14)=0.0
C...          IF(TI(J).GT.3000.) PHOTO(14) = RTS(2)*NR(1,J)*HN(J)
C...          POEXO(14)=PHOTO(14)*8.63E-5*TI(J)

C...N(2D) + O2  --> NO + O + 3.76eV
          PHOTO(14) = RTS(16)*N2D(J)*O2N(J)

          IF(I.EQ.40) WRITE(IFILE,4) Z, (PHOTO(IDUM),IDUM=1,14)
        ENDIF
C
C...GROUP II CALCULATIONS:
C...The next several (15 to 20) reactions may produce moderate
C...populations of hot O:
C
        IF (I.EQ.41.OR.I.EQ.43.OR.I.EQ.44)  THEN
C...N(2P) + O  --> N + O + 3.58eV
          PHOTO(15) = RTS(37)*N2P(J)*ON(J)
C...NO + N  --> N2 + O + 3.25eV
          PHOTO(16) = RTS(9)*N4S(J)*NNO(J)

C...N + O2  --> NO + O + 1.385eV
          PHOTO(17) = RTS(7)*N4S(J)*O2N(J)
C...N+ + O2  --> NO+ + O + 6.67eV
          PHOTO(18) = RTS(30)*NPLUS(J)*O2N(J)

C...O+(2D) + O2  --> O2+ + O + 4.865eV
          PHOTO(19) = RTS(43)*OP2D*O2N(J)
C...O+(2P) + N2  --> N2+ + O + 3.02eV
          PHOTO(20) = RTS(20)*OP2P*N2N(J)

C The next 3 reactions will only produce small populations of hot O:
C...O(1D) + O2  --> O + O2 + 1.96eV
          PHOTO(21) = RTS(34)*O1D(J)*O2N(J)
C...O2+ + N  --> NO+ + O + 4.2eV ?
          PHOTO(22) = RTS(21)*O2P(J)*N4S(J)

C...NO + N(2D)  --> N2 + O + 5.63eV
          PHOTO(23) = RTS(41)*NNO(J)*N2D(J)

C... VN2O(I) = N2*(v=I) + O --> N2(V=0) + O + (I-1)*0.3 eV (I>1)
	  DO  JEXO = 24, 27
	    PHOTO(JEXO) = VN2O(JEXO-22)
	  END DO

C...He+ + O2 --> He + O + O+ ASSUME + 1eV (although He+ gets most of it)
	  PHOTO(28) = RTS(76)*O2N(J)*HEPLUS

          IF(I.EQ.41) WRITE(IFILE,5) Z, (PHOTO(IDUM), IDUM=15,28)
        END IF

C...Calculate production rate of exothermic energy available to O 
C...(POEXO)
C... first block of reactions
        IF (I.EQ.42)  THEN
          DO  JEXO = 1, 14
            POEXO(JEXO) = PARTITION(JEXO)*EXODAT(JEXO)*PHOTO(JEXO)
          END DO
          WRITE(IFILE,4) Z, (POEXO(IDUM), IDUM=1,14)

C... second block of reactions
C... VN2O(I) = N2*(v=I) + O --> N2(V=0) + O + (I-1)*0.3 eV (I>1)
C... This corresponds to array elements 24 to 27 in following arrays.
C... Note complete quenching is assumed - upper limit
        ELSE IF (I.EQ.43)  THEN
          DO  JEXO = 15, 28
            POEXO(JEXO) = PARTITION(JEXO)*EXODAT(JEXO)*PHOTO(JEXO)
          END DO
          WRITE(IFILE,5) Z, (POEXO(IDUM), IDUM=15,28)
        END IF

	IF (I.LT.44)  RETURN
C---- This section computes spectral energy distribution. 

	  SPECRES = 0.2
	  MAXBIN  = 0
C... Determine energy bin to which each reaction contributes
	  DO  JEXO = 1, 28
	    NBIN(JEXO) = 1+INT(PARTITION(JEXO)*EXODAT(JEXO)/SPECRES)
	  END DO
C... Find maximum number of energy bins
	  DO  JEXO = 1, 28
	    MAXBIN = MAX(MAXBIN,NBIN(JEXO))
	  END DO
C... Fill XENERGY array (plots) = energy at midpoint of each energy bin
	  DO  JEXO = 1, MAXBIN
	    XENERGY(JEXO) = SPECRES*(REAL(JEXO-1)+0.5)
	  END DO
C... Determine total number density of hot O
	  TOTHOTO = 0.0
	  DO  JEXO = 1, 28
	    TOTHOTO = TOTHOTO+PHOTO(JEXO)
	  END DO
C... Determine normalized spectral amplitude, SPECAMP
	  DO  JEXO = 1, 28
	    SPECAMP(NBIN(JEXO)) = 0.0
	  END DO
	  DO  JEXO = 1, 28
	    SPECAMP(NBIN(JEXO)) = SPECAMP(NBIN(JEXO))+PHOTO(JEXO)/TOTHOTO
	  END DO
C... Now write every energy, altitude and spectral amplitude to file
C... Note: Write a line for each bin, as required by AXUM plot
	  TOTAMP = 0.0
	  DO  IDUM = 1, MAXBIN
	    WRITE(IFILE,6) XENERGY(IDUM), Z, SPECAMP(IDUM)
	    TOTAMP = TOTAMP+SPECAMP(IDUM)
	  END DO

C... Output several variables here for checking output
C...	IF (Z.LT.203.) THEN
C...	  WRITE(99,*)' Eqn.#    PHOTO      Bin#     SPECAMP'
C...	  DO  JEXO = 1, 28
C...	    WRITE(99,99) JEXO, PHOTO(JEXO), NBIN(JEXO),
C...     >                   SPECAMP(NBIN(JEXO))
C...	  END DO
C...	  WRITE(99,*)' TOTHOTO = ',TOTHOTO,'  TOTAMP = ',TOTAMP
C...	END IF
99	FORMAT(3X,I2,3X,1P,E9.2,5X,0P,I2,1P,6X,E9.2)

        RETURN
        END
C:::::::::::::::::::::::::::::: CNOP ::::::::::::::::::::::::::::::::::
C...... This routine calculates the vibrational distribution of no+
C...... Uses AFRL report by Winick et al. AFGL-TR-87-0334, Environmental
C...... Research papers, NO. 991, "An infrared spectral radiance code 
C...... for the auroral thermosphere (AARC)
C...... Written by P. Richards in February 2004
      SUBROUTINE CNOP(J,I,JP,Z,RTS,ON,O2N,N2N,NE,P1,NOP,OPLS
     >  ,N2PLS,O2P,N4S,NNO,NPLUS,N2P,PLYNOP,N2D,OP2D,mp,lp)
      IMPLICIT NONE
      INTEGER mp,lp
      INTEGER J,I,JP,INV,IJ,IV,K
      PARAMETER (INV=20)            !.. NO+(v) array dimensions
      DOUBLE PRECISION Z,RTS(99),   !.. Altitude, rate coefficients
     >  ON,O2N,N2N,NE,              !.. O, N2, O2, electron densities
     >  P1,                         !.. total source output for finding [e] 
     >  NOP,OPLS,N2PLS,O2P,         !.. NO+, )+,N2+,O2+ densities
     >  N4S,NNO,NPLUS,N2P,          !.. N(4S), NO, N+, N(2P) densities
     >  PLYNOP,                     !.. Lyman-a source
     >  N2D,OP2D,                   !.. N(2D), O+(2D) densities                
     >  NOPV(INV),NOPTOT,           !.. NO+(v) densities and total NO+ density
     >  LR(22),PR(22),              !.. storage for NO+ sources and sinks
     >  EINSCO1(INV),EINSCO2(INV),  !.. Einstein coeffs for delv=1,2
     >  LRV(INV),                   !.. NO+(v) + e rate factors
     >  PNOPV,LNOPV,                !.. Sources and sinks of NO+(v)
     >  PCASC,LRAD,                 !.. Temp total cascade source, sink
     >  K_N2_Q,P_N2_Q,L_N2_Q        !.. N2 queching rate :- coeff, source, sink

	!.. Fractions of each source going to each vib. level. Assume
	!.. N+ + O2 fractions for each source. NEED UPDATE
      REAL PRV1(INV),PRV2(INV),PRV3(INV),PRV4(INV),PRV5(INV),PRV6(INV),
     >  PRV7(INV),PRV8(INV),PRV9(INV),PRV10(INV),PRV11(INV),PRV12(INV)
      DATA PRV1/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV2/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV3/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV4/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV5/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV6/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV7/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV8/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV9/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV10/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV11/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
      DATA PRV12/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/

      DATA K_N2_Q/7.0E-12/  !.. Quenching rate coeff. by N2
      DATA EINSCO1/0.0,10.9,20.2,28.4,35.5,41.5,46.6,50.9,54.3,57.0,
     >     58.9,60.2,60.8,60.7,59.9,5*60.0/  !.. Einstein coeff delv=1
      DATA EINSCO2/0.0,0.0,.697,1.93,3.61,5.74,8.24,11.1,14.2,17.7,
     >     21.3,25.1,29.0,33.2,37.4,5*40.0/ !.. Einstein coeff delv=1
      !.. rate factors for NO+(v)+e -> N + O. Sheehan and St-Maurice 2004
      DATA LRV/1.0,19*0.3333/    
      
      !... Evaluate total production and loss rates
      PR(1)=RTS(3)*N2N*OPLS        !.. N2 + O+
      PR(2)=N2PLS*ON*RTS(10)            !.. N2+ + O
      PR(3)=O2P*N4S*RTS(21)             !.. O2+ + N(4S)
      PR(4)=O2P*NNO*RTS(23)             !.. O2+ + NO
      !.. N+ + O2 -> O2+ + N(2D,4S) or NO+ + O(1S)
      PR(5)=(RTS(30)+RTS(66)+RTS(59))*NPLUS*O2N
      PR(6)=RTS(37)*N2P*ON              !.. N2+ + O
      PR(7)=RTS(24)*OPLS*NNO            !.. O+ + NO
      PR(8)=PLYNOP*NNO                  !.. Lyman-a + NO
      PR(9)=O2P*N2D*RTS(77)             !.. Fox: O2+ + N(2D)
      PR(10)=N2PLS*NNO*RTS(80)          !.. Fox: N2+ + NO
      PR(11)=NPLUS*NNO*RTS(81)          !.. Fox: N+ + NO
      PR(12)=RTS(83)*NNO*OP2D           !.. Fox: O+(2D) + NO
      LR(1)=NE*RTS(5)                   !.. NO+ + e

      !..Total source term used in main program to calculate [e]
      P1=PR(1)+PR(2)+PR(3)+PR(4)+PR(5)+PR(6)+PR(7)+PR(8)
     >    +PR(9)+PR(10)+PR(11)+PR(12)
      IF(I.LE.0) NOP=P1/LR(1)         !.. NO+ density

      DO IJ=1,INV
        NOPV(IJ)=0.0
      ENDDO
      NOPTOT=0.0

      !.. loop down evaluating the vibrational populations. Must start 
      !.. less than INV-1 because need cascade from v+2
      DO IV=INV-4,1,-1
        !... chemical production for v=IV = total source * fraction 
        PNOPV=PR(1)*PRV1(IV)+PR(2)*PRV2(IV)+PR(3)*PRV3(IV)+
     >    PR(4)*PRV4(IV)+PR(5)*PRV5(IV)+PR(6)*PRV6(IV)+
     >    PR(7)*PRV7(IV)+PR(8)*PRV8(IV)+PR(9)*PRV9(IV)+
     >    PR(10)*PRV10(IV)+PR(11)*PRV11(IV)+PR(12)*PRV12(IV)

        !.. cascade production from v+1, v+2
        PCASC=NOPV(IV+1)*EINSCO1(IV+1)+NOPV(IV+2)*EINSCO2(IV+2)
	  LRAD=EINSCO1(IV)+EINSCO2(IV)   !.. total radiative loss

	  L_N2_Q=K_N2_Q*N2N              !.. sink of quanta by N2 quenching
	  IF(IV.EQ.1) L_N2_Q=0.0

	  P_N2_Q=K_N2_Q*N2N*NOPV(IV+1)   !.. source of quanta by N2 quenching
	  LNOPV=LR(1)*LRV(IV)            !.. recombination rate for level IV

        !.. evaluate NO+(iV) population
        NOPV(IV)=(PNOPV+PCASC+P_N2_Q)/(LNOPV+LRAD+L_N2_Q)
        NOPTOT=NOPTOT+NOPV(IV)  !.. total NO+ concentration

      ENDDO

      PR(9)=PR(9)+PR(10)+PR(11)+PR(12)         !.. for printing only
      IF(JP.EQ.1.AND.I.GT.0) WRITE(I,96)
 96   FORMAT(/5X,'NO+',34X,'PRODUCTION',93X,':LOSS RATE'/
     > ,4X,'ALT NO+(v=0) NO+(v=1) NO+(v=2)  NO+(v=3) O++N2    N2++O',
     >  4X,'O2++N4S   O2++NO   N++O2    N2P+O   O++NO   hv+NO',
     >  4X,'O2++N2D   N2++NO   N++NO   OP2P+NO  NO++e')
      IF(I.GT.0) WRITE(I,'(F6.1,1P,22E9.2)')  Z,
     >   (NOPV(K),K=1,4),(PR(K),K=1,12),LR(1)*NOP

      RETURN
      END

