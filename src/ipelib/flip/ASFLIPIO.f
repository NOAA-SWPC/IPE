C.................... FLIPIO.FOR;2 ..........12-JUL-1993 16:58:36.77 
C..... This file contains all the subroutines for reading and writing 
C..... the data files for the FLIP model. As from February 1993, there 
C..... are no longer separate IO routines for the VAX and the CRAY.
C..... They were written by Teri Statum of Boeing at MSFC toward
C..... the end of 1987. Modified by Phil Richards 1993.
C
C:::::::::::::::::::::::::::::: WF8DAT ::::::::::::::::::::::::::::::::
C..... This routine is called to write the geophysical data in the top
C..... of the main FLIP data file F8.DAT. The file is written at the
C..... MONRUN phase.
      SUBROUTINE WF8DAT(IVLPS,JQ,PCO,Z0,SCAL,GRD,GLONN,GLONF,GLATN,
     >  GLATF,ZLBDY,FPAS,UV,HFAC,TIMEON,TIMEOF,AURFAC,PLAT,PLON,UVFAC)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL RATFAC(99)
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,EUV,UVFAC(59),GLONN,GLONF
     >  ,GLATN,GLATF,GRD,TIMEON,TIMEOF,AURFAC,PLAT,PLON
      COMMON/PALT/ZLO,ZHI,JWR,JTI
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/ALT/Z(401),DT,DH,THF,EPS,TF,ITF,JMAX,JMAX1,ITER,ION
      COMMON/LPS/SZA(401),EPSN,DC,ISPC,I1,I2,IPV,JTIMAX,IPP,ISKP
      COMMON/CONFAC/IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP
      DATA DUMMY/0.0/ !.. not used - was PCOFM for FAIM
      EUV=0.0  !.. unused variable
C...... Note second I1 is dummy - not used
      WRITE(8,90)  JMAX,IDAY,IVLPS,I1,I2,JQ,I1,ISKP,IWIND,IVIBN2
     >  ,IODDN,INPLUS,IHEP,ITEMP,IHPOP,TIMEON,TIMEOF,AURFAC,DUMMY
      WRITE(8,91) PCO,Z0,SCAL,GRD,F107,F107A,AP(1),BLON,EUV,GLONN
     > ,GLONF,GLATN,GLATF,ZLBDY,FPAS,UV,HFAC,(Z(J),J=1,JMAX,ISKP),
     >  (UVFAC(J),J=1,59),(RATFAC(J),J=1,99),PLAT,PLON
 90   FORMAT(I4,I8,13I5,4F7.2)
 91   FORMAT(11F11.4)
      RETURN
      END
C:::::::::::::::::::::::::::::: WF8D2 ::::::::::::::::::::::::::::::::
C..... This routine is called to write the actual densities and temperatures
C..... in the F8.DAT file.
      SUBROUTINE WF8D2(JK,JC,UNN,UNC,IVLPS,JQ)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,EXFACN,EXFACS
      REAL NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS,TMSN,TMSS
      INTEGER CONDEN,CONVEL,TPRIN,TMODON
      COMMON/PALT/ZLO,ZHI,JWR,JTI
      COMMON/VN/U(2,401),BG(401),BM(401),GR(2,401),R(2,401),SL(401)
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/ALT/Z(401),DT,DH,THF,EPS,TF,ITF,JMAX,JMAX1,ITER,ION
      COMMON/STAW/STR(12,401),RCON(401)
      COMMON/LPS/SZA(401),EPSN,DC,ISPC,J1,J2,IPV,JTIMAX,IPP,ISKP
      COMMON/CONFAC/IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP
      COMMON/VAGRID/JBV,JVB,JNB
      COMMON/NDFAC/TPRIN,NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS
     >  ,TMSN,TMSS,EXFACN,EXFACS,TMODON
      TPRIN=0
C
      WRITE(8,90) JTI,JK,JC,(AP(I),I=1,7),F107,TINFN,TINFS,NDFNOR,NDFSOU
     >   ,EXFACN,EXFACS
      WRITE(8,91) SEC,N(1,JK),N(1,JC),N(1,JQ),N(2,JQ),TI(2,JQ),TI(3,JQ)
     >  ,UNN,UNC,XIONN(3,JQ),XIONN(4,JQ)
C----- The densities and velocities are converted to integers for storing
      WRITE(8,92) ((CONDEN(N(I,J)),J=1,JMAX,ISKP),I=1,3)
     >  ,((CONDEN(TI(I,J)),J=1,JMAX,ISKP),I=2,3)
      WRITE(8,93) ((CONVEL(XIONV(I,J)),J=1,JMAX,ISKP),I=1,2)
      IF(IODDN.EQ.1) WRITE(8,92) 
     >  ((CONDEN(STR(I,J)),J=1,JBV,ISKP),I=8,10)
      IF(IODDN.EQ.1) WRITE(8,92) 
     >  ((CONDEN(STR(I,JMAX+1-J)),J=1,JBV,ISKP),I=8,10)
      IF(IVIBN2.EQ.1) WRITE(8,92) 
     >  ((CONDEN(STR(I,J)),J=1,JBV,ISKP),I=1,6)
      IF(IVIBN2.EQ.1) WRITE(8,92) 
     >  ((CONDEN(STR(I,JMAX+1-J)),J=1,JBV,ISKP),I=1,6)
      IF(IVIBN2.GE.1) WRITE(8,92) (CONDEN(RCON(J)),J=1,JMAX,ISKP)
      IF(IABS(IVLPS).GE.9) THEN
        WRITE(8,92) (CONDEN(XIONN(3,J)),J=1,JMAX,ISKP)
        WRITE(8,93) (CONVEL(XIONV(3,J)),J=1,JMAX,ISKP)
      ENDIF
      IF(IABS(IVLPS).GE.11) THEN
        WRITE(8,92) (CONDEN(XIONN(4,J)),J=1,JMAX,ISKP)
        WRITE(8,93) (CONVEL(XIONV(4,J)),J=1,JMAX,ISKP)
      ENDIF
      RETURN
 90   FORMAT(3I5,22F7.2)
 91   FORMAT(F11.0,1P,11E12.4)
 93   FORMAT(15I6)
 92   FORMAT(15I6)
      END
C:::::::::::::::::::::::::::: RF8DAT :::::::::::::::::::::::::::::::::
C....... This routine is used for reading the geophysical parameters off
C....... the top of the main data file (F8.DAT) that is created by the FLIP
C....... model. This file is used in the MONPRN phase when printing selected
C....... data.
      SUBROUTINE RF8DAT(IVLPS,JQ,IW,PCO,Z0,SCAL,GRD,GLONN,GLONF,
     >  GLATN,GLATF,ZLBDY,FPAS,UV,HFAC,JMAX,IDAY,I1,I2,F107,
     >  F107A,AP,BLON,EUV,Z,UVFAC,RATFAC,TIMEON,TIMEOF,AURFAC
     >  ,PLAT,PLON)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL F107,F107A,AP(7),BLON,EUV,UVFAC,GLONN,GLONF,
     > GLATN,GLATF,GRD,RATFAC,TIMEON,TIMEOF,AURFAC,PLAT,PLON
      DIMENSION RATFAC(99),Z(401),UVFAC(59)
      COMMON/LPS/SZA(401),EPSN,DC,ISPC,J1,J2,IPV,JTIMAX,IPP,ISKP
      COMMON/CONFAC/IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP
      DATA DUMMY/0.0/ !.. not used - was PCOFM for FAIM
C
C...... Note IW is not used anywhere
      READ(1,90,ERR=253,END=254) JMAX,IDAY,IVLPS,I1,I2,JQ,IW,ISKP
     >  ,IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP,TIMEON,TIMEOF
     >  ,AURFAC,DUMMY
      READ(1,91,ERR=253) PCO,Z0,SCAL,GRD,F107,F107A,AP(1)
     > ,BLON,EUV,GLONN,GLONF,GLATN,GLATF,ZLBDY,FPAS,UV,HFAC
     > ,(Z(J),J=1,JMAX,ISKP),(UVFAC(J),J=1,59),(RATFAC(J),J=1,99)
     > ,PLAT,PLON
 254  RETURN
 253  CONTINUE
      WRITE(6,*) ' ** ERROR: NO valid data on the FLIP print file ***'
      CALL RUN_ERROR    !.. print ERROR warning in output files
      STOP
 90   FORMAT(I4,I8,13I5,4F7.2)
 91   FORMAT(11F11.4)
      END
C:::::::::::::::::::::::::::::: RF8D2 ::::::::::::::::::::::::::::::::
C....... This routine is used for reading the actual data off  the main
C....... data file (F8.DAT) that is created by the FLIP model. This file
C....... is used in the MONPRN phase when printing selected data.
      SUBROUTINE RF8D2(JTX,JK,JC,NMF2N,NMF2F,NOPEQ,NHPEQ,TIEQ,TEEQ,UNN,
     >  UNC,IVLPS,IFLAG,JQ,SEC,N,JMAX,TI,U,STR,RCON,AP,NHEPEQ,F107
     >  ,NNPEQ)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,AP(7),F107,DENCON,VELCON,EXFACN,EXFACS
      REAL NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS,TMSN,TMSS
      INTEGER RN(4,401),RTI(3,401),RU(2,401),RSTR(12,401),RRCON(401)
     >  ,RUE(2,401),RXIONN(2,401),RXIONV(2,401),TPRIN,TMODON
      DIMENSION N(4,401),TI(3,401),U(2,401),STR(12,401),
     >  RCON(401)
      COMMON/LPS/SZA(401),EPSN,DC,ISPC,J1,J2,IPV,JTIMAX,IPP,ISKP
      COMMON/CONFAC/IWIND,IVIBN2,IODDN,INPLUS,IHEP,ITEMP,IHPOP
      COMMON/VAGRID/JBV,JVB,JNB
      COMMON/NDFAC/TPRIN,NDFNOR,NDFSOU,TINFN,TINFS,DELTIN,DELTIS
     >   ,TMSN,TMSS,EXFACN,EXFACS,TMODON
      TPRIN=1
C
C------ Make sure those variables not read every time are zero
       DO 10 J=1,JMAX
          RCON(J)=1.0
          XIONN(4,J)=0.0
          XIONV(4,J)=0.0
          XIONN(3,J)=0.0
          XIONV(3,J)=0.0
 10    CONTINUE
       DO 15 J=1,JMAX
       DO 15 I=1,11
          STR(I,J)=0.0
 15    CONTINUE
C
C------- Read the data into real variables
      READ(1,90,END=233,ERR=233) JTX,JK,JC,(AP(I),I=1,7),F107
     > ,TINFN,TINFS,NDFNOR,NDFSOU,EXFACN,EXFACS
      READ(1,91,END=233,ERR=233) SEC,NMF2N,NMF2F,NOPEQ,NHPEQ,TIEQ
     > ,TEEQ,UNN,UNC,NHEPEQ,NNPEQ
C----- densities and velocities are read as integers and then converted
      READ(1,92,ERR=233) ((RN(I,J),J=1,JMAX,ISKP),I=1,3)
     >  , ((RTI(I,J),J=1,JMAX,ISKP),I=2,3)
      READ(1,93,ERR=233)  ((RU(I,J),J=1,JMAX,ISKP),I=1,2)
      IF(IODDN.EQ.1) READ(1,92,ERR=233) 
     >    ((RSTR(I,J),J=1,JBV,ISKP),I=8,10)
      IF(IODDN.EQ.1) READ(1,92,ERR=233) 
     >    ((RSTR(I,JMAX+1-J),J=1,JBV,ISKP),I=8,10)
      IF(IVIBN2.EQ.1) READ(1,92,ERR=233) 
     >    ((RSTR(I,J),J=1,JBV,ISKP),I=1,6)
      IF(IVIBN2.EQ.1) READ(1,92,ERR=233) 
     >    ((RSTR(I,JMAX+1-J),J=1,JBV,ISKP),I=1,6)
      IF(IVIBN2.GE.1) READ(1,92,ERR=233) (RRCON(J),J=1,JMAX,ISKP)
      IF(IABS(IVLPS).GE.9) READ(1,92,ERR=233) (RN(4,J),J=1,JMAX,ISKP)
      IF(IABS(IVLPS).GE.9) READ(1,93,ERR=233) (RUE(1,J),J=1,JMAX,ISKP)
      IF(IABS(IVLPS).GE.11) READ(1,92,ERR=233) 
     >   (RXIONN(1,J),J=1,JMAX,ISKP)
      IF(IABS(IVLPS).GE.11) READ(1,93,ERR=233)
     >   (RXIONV(1,J),J=1,JMAX,ISKP)
C
C------- Transfer densities and velocities into double precision variables
       DO 20 J=1,JMAX
          N(1,J)=DENCON(RN(1,J))
          N(2,J)=DENCON(RN(2,J))
          N(3,J)=DENCON(RN(3,J))
          IF(IABS(IVLPS).GE.9) XIONN(3,J)=DENCON(RN(4,J))
          TI(2,J)=DENCON(RTI(2,J))
          TI(3,J)=DENCON(RTI(3,J))
          XIONV(1,J)=VELCON(RU(1,J))
          XIONV(2,J)=VELCON(RU(2,J))
          IF(IVIBN2.GE.1) RCON(J)=DENCON(RRCON(J))
          IF(IABS(IVLPS).GE.9) XIONV(3,J)=VELCON(RUE(1,J))
          IF(IABS(IVLPS).GE.9) XIONN(4,J)=DENCON(RXIONN(1,J))
          IF(IABS(IVLPS).GE.9) XIONV(4,J)=VELCON(RXIONV(1,J))
 20    CONTINUE
       IF(IVIBN2.EQ.1) THEN
          DO 30 J=1,JBV
          DO 40 I=1,6
             STR(I,J)=DENCON(RSTR(I,J))
             STR(I,JMAX+1-J)=DENCON(RSTR(I,JMAX+1-J))
 40       CONTINUE
 30       CONTINUE
       ENDIF
       IF(IODDN.EQ.1) THEN
          DO 55 J=1,JBV
          DO 50 I=8,10
             STR(I,J)=DENCON(RSTR(I,J))
             STR(I,JMAX+1-J)=DENCON(RSTR(I,JMAX+1-J))
 50       CONTINUE
 55       CONTINUE
       ENDIF
C
 221  CONTINUE
      RETURN
 233  CONTINUE
C-----  Return error flag to main program
      IFLAG=1
      RETURN
 90   FORMAT(3I5,22F7.2)
 91   FORMAT(F11.0,1P,11E12.4)
 93   FORMAT(15I6)
 92   FORMAT(15I6)
      END
C:::::::::::::::::::::::::::::: RRSTRT ::::::::::::::::::::::::::::::::
C..... This routine is responsible for reading the file containing the
C..... initial densities and temperatures created in a previous run as an
C..... F2.DAT file.
      SUBROUTINE RRSTRT(IU,PCO,Z0,SCAL,GRD,GLONN,GLONF,GLATN,GLATF,ISOL
     >   ,UNN,UNC,HMF2N,HMF2S,NMF2N,NMF2S,PLAT,PLON,UVFAC)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,EUV,UVFAC(59),GLONN,GLONF,
     >  GLATN,GLATF,GRD,PLAT,PLON
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/VN/U(2,401),BG(401),BM(401),GR(2,401),R(2,401),SL(401)
      COMMON/ALT/Z(401),DT,DH,THF,EPS,TF,ITF,JMAX,JMAX1,ITER,ION
      COMMON/STAW/STR(12,401),RCON(401)

      READ(IU,END=391,ERR=391) JMAX,IDAY,PCO,Z0,SCAL,SEC,GRD,F107
     >   ,F107A,AP(1),BLON,EUV,GLONN,GLONF,GLATN,GLATF,UNN,UNC
     >   ,HMF2N,HMF2S,NMF2N,NMF2S,PLAT,PLON
	
      READ(IU) ((N(I,J),J=1,JMAX),I=1,4),((TI(I,J),J=1,JMAX),I=1,3)
     >  ,((XIONV(I,J),J=1,JMAX),I=1,4),((STR(I,J),J=1,JMAX),I=1,12)
     >  ,((U(I,J),J=1,JMAX),I=1,2),(RCON(J),J=1,JMAX)
     >  , ((XIONN(I,J),J=1,JMAX),I=1,4)
      READ(IU) (UVFAC(J),J=1,59)
C----- Return ISOL=-1 if file is OK else =0 for new field line
      ISOL=-1
      RETURN
 391  CONTINUE
      ISOL=0
      RETURN
      END
C:::::::::::::::::::::::::::::: WRSTRT ::::::::::::::::::::::::::::::::
C..... This routine is called to write the F2.DAT file after each successful
C..... time step. It contains the necessary information to continue the
C..... run, including the densities and temperatures.
      SUBROUTINE WRSTRT(PCO,Z0,SCAL,GRD,GLONN,GLONF,GLATN,GLATF
     >   ,UNN,UNC,HMF2N,HMF2S,NMF2N,NMF2S,PLAT,PLON,UVFAC)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL SEC,DEC,ETRAN,F107,F107A,AP,BLON,EUV,UVFAC(59),GLONN,
     >  GLONF,GLATN,GLATF,GRD,PLAT,PLON
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/VN/U(2,401),BG(401),BM(401),GR(2,401),R(2,401),SL(401)
      COMMON/ALT/Z(401),DT,DH,THF,EPS,TF,ITF,JMAX,JMAX1,ITER,ION
      COMMON/STAW/STR(12,401),RCON(401)
      EUV=0.0  !.. unused variable
      WRITE(2) JMAX,IDAY,PCO,Z0,SCAL,SEC,GRD,F107,F107A,AP(1),BLON
     > ,EUV,GLONN,GLONF,GLATN,GLATF,UNN,UNC,HMF2N,HMF2S,NMF2N,NMF2S
     > ,PLAT,PLON
      WRITE(2) ((N(I,J),J=1,JMAX),I=1,4),((TI(I,J),J=1,JMAX),I=1,3)
     >  ,((XIONV(I,J),J=1,JMAX),I=1,4),((STR(I,J),J=1,JMAX),I=1,12)
     >  ,((U(I,J),J=1,JMAX),I=1,2),(RCON(J),J=1,JMAX)
     >  , ((XIONN(I,J),J=1,JMAX),I=1,4)
      WRITE(2) (UVFAC(J),J=1,59)
      REWIND 2
      RETURN
      END
C:::::::::::::::::::::::::::::::: TFILE :::::::::::::::::::::::::
C------- This file checks to see if a file exists for the FORTRAN unit
C------- JUNIT and returns 0 for no and 1 for yes
       SUBROUTINE TFILE(JUNIT,IOPEN)
       CHARACTER ACHAR
       IOPEN=1
       READ(JUNIT,90,END=70,ERR=70) ACHAR 
 90    FORMAT(A)
       REWIND JUNIT
       RETURN
 70    CONTINUE
C------- File did not exist
       IOPEN=0 
       RETURN
       END
C:::::::::::::::::::::::::::::: VELCON ::::::::::::::::::::::::::::::::
C--- This function is the twin of function CONVEL. It converts an integer
C--- back into a velocity. See CONVEL for details
C--- P. Richards February 1996
      REAL FUNCTION VELCON(INTNUM)
      INTEGER INTNUM
      IF(INTNUM.LT.0.0) THEN
         VELCON=-10**(REAL(-INTNUM-22000)/1000.0)
      ELSE
         VELCON=10**(REAL(INTNUM-22000)/1000.0)
      ENDIF
      RETURN
      END
C:::::::::::::::::::::::::::::: CONVEL ::::::::::::::::::::::::::::::::
C--- This function and VELCON are used to convert a velocity into an
C--- small integer. First, the number is converted to a log then
C--- multiplied by 1000 to get the right precision. Negative numbers
C--- have to be treated differently so 22 is added onto the log
C--- to ensure that part is positive, then the sign of the velocity is 
C--- attached. The 1.0E-22 takes care of very small numbers which are
C--- effectively 0.00
C--- P. Richards February 1996 
      INTEGER FUNCTION CONVEL(DPNUM)
      DOUBLE  PRECISION DPNUM
      IF(DPNUM.LT.0.0) THEN
         CONVEL=-NINT((22+DLOG10(DABS(DPNUM-1.0E-22)))*1000.0)
      ELSE
         CONVEL=NINT((22+DLOG10(DABS(DPNUM+1.0E-22)))*1000.0)
      ENDIF
      RETURN
      END
C:::::::::::::::::::::::::::::: DENCON ::::::::::::::::::::::::::::::::
C--- This function with CONDEN is used to convert a density to an integer
C--- for efficient storage. See CONDEN for details
C--- P. Richards, February 1996
      REAL FUNCTION DENCON(INTNUM)
      INTEGER INTNUM
         DENCON=10**(REAL(INTNUM)/1000.0)
      RETURN
      END
C:::::::::::::::::::::::::::::: CONDEN ::::::::::::::::::::::::::::::::
C--- The density is converted to a log, then to a sensible sized integer
C--- to preserve sufficient precision by multiplying by 1000. Function
C--- DENCON is used to retrieve the actual densities.
C--- P. Richards, February 1996
      INTEGER FUNCTION CONDEN(DPNUM)
      DOUBLE  PRECISION DPNUM
         CONDEN=NINT(DLOG10(DPNUM+1.0E-30)*1000.0)
      RETURN
      END
C
C::::::::::::::::: WRITE_Q ::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE WRITE_Q(JMIN, JMAX, SCALK, Q)
      DOUBLE PRECISION Q(401), SCALK

      WRITE(8, *) JMIN, JMAX, SCALK
      WRITE(8, '(1p, 9D14.6)') (Q(J), J=JMIN, JMAX)

      RETURN
      END
C
C::::::::::::::::: READ_Q ::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE READ_Q(SCALK, Q)
      DOUBLE PRECISION Q(401), SCALK

      READ(1, *) JMIN, JMAX, SCALK
      READ(1, *) (Q(J), J=JMIN, JMAX)

      RETURN
      END
C
C::::::::::::::::: WRITE_BQ ::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE WRITE_BQ(JMIN, JMAX, SCALK, Q)
      DOUBLE PRECISION Q(401), SCALK

      WRITE(2) JMIN, JMAX, SCALK
      WRITE(2) (Q(J), J=JMIN, JMAX)

      RETURN
      END
C
C::::::::::::::::: READ_BQ ::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE READ_BQ(SCALK, Q)
      DOUBLE PRECISION Q(401), SCALK
      
      READ(99,END=100,ERR=100) JMIN, JMAX, SCALK
      READ(99) (Q(J), J=JMIN, JMAX)

 100  RETURN
      END
