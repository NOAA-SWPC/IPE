C.................... RSTEMC_EXB.FOR .............................
C.... this subroutine sets up the error functions of the time-dependent
C.... ion and electron energy equations for subsequent solution by the 
C.... Newton iterative procedure. It calls a function for Ti and one for
C.... Te. This version was created by Silvy John to implement the
C.... flux preserving method (1996). The flux preserving approach
C.... is explained in the paper Torr, et al., J. Geophys. Res., 95, 
C.... 21,147-21,168, 1990.
C.... The routines on this file were documented and cleaned up for CTIPe
C....  by P. Richards, July 2011
      SUBROUTINE TFIJ(J,    !.. point on the field line
     >              ILJ,    !.. used to print diagnostics if JSJ NE 0 or 1
     >               DT,    !.. Time step in seconds
     >                N,    !.. O+, H+, He+, minor ion densities array
     >               TI,    !.. Ion and electron temperatures
     >                F,    !.. OUTPUT: The integrated continuity equation
     >            TISAV,mp,lp)    !.. temperatures at previous time step
      USE FIELD_LINE_GRID !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      IMPLICIT NONE
      INTEGER mp,lp
      INTEGER J,ILJ,ID
      DOUBLE PRECISION DT,N(4,FLDIM),TI(3,FLDIM),F(20),TISAV(3,FLDIM)
      ID=3                         !..unit for diagnostics in North
      IF(J.GT.(JMAX-JMIN)/2) ID=9  !..unit for diagnostics in South

      !.... set up the ion and electron energy equations
      CALL GET_ION_TEMP(J,ILJ,ID,DT,N,TI,TISAV,F)
      CALL GET_ELECTRON_TEMP(J,ILJ,ID,DT,N,TI,TISAV,F,mp,lp)
      RETURN
      END
C:::::::::::::::::::::::::: GET_ION_TEMP ::::::::::::::::::::::::::
C.. This subroutine is similar to GET_ELECTRON_TEMP. Look at the Te function
C.. for more explanation. It determines the interpolated production and
C.. loss processes for ions. It calls TERLIN to do the interpolation
      SUBROUTINE GET_ION_TEMP(
     >                  J,    !.. point on the field line
     >                ILJ,    !.. used to print diagnostics if JSJ NE 0 or 1
     >                 ID,    !.. unit to print diagnostics
     >                 DT,    !.. Time step in seconds
     >                  N,    !.. O+, H+, He+, minor ion densities array
     >                 TI,    !.. Ion and electron temperatures
     >              TISAV,    !.. temperatures at previous time step
     >                  F)    !.. OUTPUT: The integrated continuity equation
      USE FIELD_LINE_GRID  !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE     !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc.
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER J,ILJ,ID,I,K
      DOUBLE PRECISION DT,N(4,FLDIM),TI(3,FLDIM),TISAV(3,FLDIM),F(20)
      DOUBLE PRECISION BOLTZ,NE           !.. Boltzmann constant, e density
      DOUBLE PRECISION DS_UP,DS_LO        !.. distances between grid points
      DOUBLE PRECISION CE_ION(3)          !.. e-i heat transfer 
      DOUBLE PRECISION RCheat(3)          !.. Ion ring current heating
      DOUBLE PRECISION KI_UP,KI_LO,KI(3)  !.. Thermal conductivity coeffs
      DOUBLE PRECISION B_UP,B_LO          !.. Magnetic field at half points
      DOUBLE PRECISION HFLUX_UP,HFLUX_LO  !.. Heat flux components
      DOUBLE PRECISION DTI(3),DERIVT      !.. dT/dt and integrated dT/dt
      DOUBLE PRECISION ElossI             !.. Integrated e - i energy exchange
      DOUBLE PRECISION IlossN             !.. Integrated i - n cooling rate
      DOUBLE PRECISION RCrate             !.. Integrated additional ion heating rate
      DOUBLE PRECISION NI                 !.. Ion density
      DOUBLE PRECISION HLOSI(3),ICOOL(22) !.. Ion- neutral cooling rates
      DOUBLE PRECISION DIVV_UP,DIVV_LO    !.. Div of velocity cmpts
      DOUBLE PRECISION IVEL(2,3),VLIMIT   !.. Velocity cmpts and limiting velocity
      DATA BOLTZ/1.3807E-16 /             !.. Boltzmann constant
      DATA VLIMIT/1.0E+6 /                !.. Velocity limit in cm/sec          
  
      !..evaluate terms at J-1, J, J+1 for integration, J is the
      !..actual grid point. All terms divided by magnetic field (BM)
      DO I=1,3
        K=I+J-2
        !.. Electron density = total ion density
        NE=XIONN(1,K)+XIONN(2,K)+XIONN(3,K)+XIONN(4,K)+XIONN(5,K)+
     >       XIONN(6,K)+XIONN(7,K)+XIONN(8,K)+XIONN(9,K)
        DTI(I)=1.5*BOLTZ*NE*(TI(2,K)-TISAV(2,K))/BM(K)/DT
        !.. direct heating rate for ions (e.g. ring current)
        RCheat(I)=EHT(1,K)*1.6E-12/BM(K)
        !.. Electron - ion energy exchange
        CALL GET_EcoolI(K,FLDIM,TI(3,K),TI(2,K),BM(K),CE_ION(I))
        !.. Cooling rate coefficient to neutrals
        CALL GET_IcoolN(N(1,K),N(2,K),ON(K),N2N(K),O2N(K)
     >     ,HN(K),HE(K),TN(K),TI(2,K),ICOOL)
        HLOSI(I)=ICOOL(1)*(TI(2,K)-TN(K))/BM(K)
        !.. determine the thermal conductivity KI 
        KI(I)=1.84E-8*(N(1,K)+4*N(2,K))*(TI(2,K)**2.5)/NE
        IVEL(1,I)=XIONV(1,K)
        IVEL(2,I)=XIONV(2,K)
        !.. Make sure velocities cannot get too high
        IF(IVEL(1,I).GT.VLIMIT) IVEL(1,I)=VLIMIT
        IF(IVEL(1,I).LT.-VLIMIT) IVEL(1,I)=-VLIMIT
        IF(IVEL(2,I).GT.VLIMIT) IVEL(2,I)=VLIMIT
        IF(IVEL(2,I).LT.-VLIMIT) IVEL(2,I)=-VLIMIT
      ENDDO

      !.. distances between grid points
      DS_UP=SL(J+1)-SL(J)
      DS_LO=SL(J)-SL(J-1)

      !.. Now perform the integration to get components of the energy equation
      !.. DERIVT, ElossI, IlossN , RCrate
      CALL TERLIN(DTI(1),DTI(2),DTI(3),DS_LO,DS_UP,DERIVT)
      CALL TERLIN(CE_ION(1),CE_ION(2),CE_ION(3),DS_LO,DS_UP,ElossI)
      CALL TERLIN(HLOSI(1),HLOSI(2),HLOSI(3),DS_LO,DS_UP,IlossN)
      CALL TERLIN(RCheat(1),RCheat(2),RCheat(3),DS_LO,DS_UP,RCrate)

      B_UP=(BM(J)+BM(J+1))*0.5    !.. Upper magnetic field
      B_LO=(BM(J)+BM(J-1))*0.5    !.. Lower magnetic field

      !.. Upper and lower ion thermal conductivity coefficient
      KI_UP=(KI(3)+KI(2))/2.0      !.. Upper conductivity coefficient
      KI_LO=(KI(2)+KI(1))/2.0      !.. Lower conductivity coefficient

      !.. Upper and lower ion thermal conduction heat flux
      HFLUX_UP=KI_UP*(TI(2,J+1)-TI(2,J))/DS_UP/B_UP
      HFLUX_LO=KI_LO*(TI(2,J)-TI(2,J-1))/DS_LO/B_LO

      !.. Calculate the divergence of the velocity terms (integrated). 
      !.. DIVV_UP and DIVV_LO are the values at the upper and lower .
      !.. bounds of the integral.
      !.. 2nd last term on RHS in Bailey, PSS 1983 equation 15
      !.. The O+ and H+ contributions are treated separately in each term.
      ! DIVV_UP=BOLTZ*0.5*(XIONN(1,J+1)*TI(2,J+1)+XIONN(1,J)*TI(2,J))*
      !>     (IVEL(1,3)-IVEL(1,2))/B_UP
      !>   +BOLTZ*0.5*(XIONN(2,J+1)*TI(2,J+1)+XIONN(2,J)*TI(2,J))*
      !>     (IVEL(2,3)-IVEL(2,2))/B_UP
      ! DIVV_LO=BOLTZ*0.5*(XIONN(1,J)*TI(2,J)+XIONN(1,J-1)*TI(2,J-1))*
      !>     (IVEL(1,2)-IVEL(1,1))/B_UP
      !>   +BOLTZ*0.5*(XIONN(2,J)*TI(2,J)+XIONN(2,J-1)*TI(2,J-1))*
      !>     (IVEL(2,2)-IVEL(2,1))/B_UP

      DIVV_UP=0.00
      DIVV_LO=0.00

      !.. Integrated energy equation for ions
      F(1)=DERIVT-ElossI+IlossN-HFLUX_UP+HFLUX_LO-RCrate
     >   +DIVV_UP-DIVV_LO  !.. div.V term
      
      !.. diagnostic print
      IF(ILJ.EQ.4) THEN
        WRITE(ID,501)
        WRITE(ID,5555) Z(J),TI(3,J),TI(2,J),TN(J),F(1),DERIVT,ElossI,
     >    IlossN,HFLUX_UP,HFLUX_LO,DIVV_UP,DIVV_LO,
     >    HFLUX_UP-HFLUX_LO,DIVV_UP-DIVV_LO,
     >    (HFLUX_UP-HFLUX_LO)/(1.0E-33+DIVV_UP-DIVV_LO)
      ENDIF
 501  FORMAT(' Ti: ALT      Te       Ti       Tn      F        dT/dt'
     > ,6X,'C_ie      C_in    HFLUX_UP  HFLUX_LO   DIVV_UP  DIVV_LO'
     > ,2X,'  DIVCOND  DIVCONV')
 5555   FORMAT(4F9.1,1P,22E10.2)
      RETURN
      END
C:::::::::::::::: GET_IcoolN ::::::::::::::::::::::
C.... Ion loss rate coefficients to neutrals (KIN) Rees & Roble(1975) P220
      SUBROUTINE GET_IcoolN(OPLUS,   !.. O+ density
     >                      HPLUS,   !.. H+ density
     >                         ON,   !.. O density
     >                        N2N,   !.. N2 density
     >                        O2N,   !.. O2 density
     >                         HN,   !.. H density
     >                        HEN,   !.. He density
     >                         TN,   !.. Ion temperature
     >                         TI,   !.. Neutral temperature
     >                      ICOOL)   !.. OUTPUT: Ion cooling to neutrals
      IMPLICIT NONE
      DOUBLE PRECISION OPLUS,HPLUS,ON,N2N,O2N,HN,HEN,TN,TI,ICOOL(22)
        !.. Resonant charge exchange
        ICOOL(2)=1.6E-26*(0.21*OPLUS*ON)*SQRT(TI+TN) !.. O+ - O
        ICOOL(3)=1.6E-26*(1.40*HPLUS*HN)*SQRT(TI+TN) !.. H+ - H
        !..Polarization interactions
        ICOOL(4)=1.6E-26*OPLUS*(6.6*N2N)             !.. O+ - N2
        ICOOL(5)=1.6E-26*OPLUS*(2.8*HEN)             !.. O+ - He
        ICOOL(6)=1.6E-26*OPLUS*(5.8*O2N)             !.. O+ - O2
        ICOOL(7)=1.6E-26*OPLUS*(5.6*ON)              !.. O+ - O
        ICOOL(8)=1.6E-26*HPLUS*(5.5*HEN)             !.. H+ - He
        ICOOL(9)=1.6E-26*HPLUS*(1.9*HN)              !.. H+ - H
        ICOOL(10)=1.6E-26*(0.36*OPLUS*HN)*SQRT(TN)   !.. O+ - H
        ICOOL(11)=1.6E-26*(0.40*HPLUS*ON)*SQRT(TN)   !.. H+ - O
        !..total cooling rate coefficient
        ICOOL(1)=ICOOL(2)+ICOOL(3)+ICOOL(4)+ICOOL(5)+ICOOL(6)+ICOOL(7)+
     >    ICOOL(8)+ICOOL(9)+ICOOL(10)+ICOOL(11)
      RETURN
      END
C:::::::::::::::::::::::::::::::: GET_EcoolI ::::::::::::::::::::::
C.. Electron-ion energy exchange rate. Note that electrons usually 
C.. heat ions but ions can heat electrons at night
      SUBROUTINE GET_EcoolI(K,   !.. point on field line
     >                  FLDIM,   !.. Field line dimension
     >                     TE,   !.. Electron temperature
     >                     TI,   !.. Ion temperature
     >                     BM,   !.. Ion temperature
     >                 CE_ION)   !.. OUTPUT: Ion cooling to neutrals
      USE ION_DEN_VEL        !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER K,FLDIM                   !.. See I/O parameter comments
      DOUBLE PRECISION TE,TI,BM,CE_ION  !.. See I/O parameter comments
      DOUBLE PRECISION NE,NI            !.. Electron and ion density
      DOUBLE PRECISION DEBYE,COULOG     !.. Factors for i-e cooling rate
	
        !.. Electron density = total ion density
        NE=XIONN(1,K)+XIONN(2,K)+XIONN(3,K)+XIONN(4,K)+XIONN(5,K)+
     >       XIONN(6,K)+XIONN(7,K)+XIONN(8,K)+XIONN(9,K)
        !.. Mass weighted ion density
        NI=XIONN(1,K)/16+XIONN(2,K)+XIONN(3,K)/4+XIONN(4,K)/14+
     >       XIONN(5,K)/30+XIONN(6,K)/32+XIONN(7,K)/28+XIONN(8,K)/16+
     >       XIONN(9,K)/16

        !.. Electron- ion cooling rate KEI. The Coulomb logarithm from Itikawa
        !.. JATP 1975, page 1001. Previously assumed to be 15.0
        DEBYE=1.602E-12*SQRT(4.0*3.1415926*NE/(1.381E-16*TE))
        COULOG=0.43086+DLOG(TE)-DLOG(DEBYE)
        CE_ION=1.232e-17*NE*NI*(TE-TI)/BM/(TE**1.5)
        CE_ION= CE_ION*COULOG/15.0    !.. Itikawa correction
      RETURN
      END
C::::::::::::::::::::::::::: GET_ELECTRON_TEMP :::::::::::::::::::::
C.. This subroutine sets up the integrated electron energy equation F(2) for
C.. the Newton solver. It evaluates the electron cooling rate,  dT/dt, 
C.. heating rate, heat_flux, and cooling rate for three adjacent grid points 
C.. j,j-1 and j+1 and uses them to determine the values at the upper and lower
C.. midpoints by interpolation. The upper midpoint corresponds to the
C.. point midway between j and j+1 whereas lower midpoint corresponds to
C.. that between j-1 and j. The other terms are HFLUX_UP and HFLUX_LO which
C.. are the integrated values of the heat flow between j+1 and j and j and
C.. j-1 respectively.The basic equation is 
C.. F(2)= DT/DT(DERIVT)+DIV(HEAT_FLUX)+ELEC_ION_COOLING_RATE(ElossI)
C..      +ELEC_NEUT_COOLING(ElossN)-HEATING_RATE(EHrate)=0
C.. TERLIN interpolates and the integration is done between the midpoints.
      SUBROUTINE GET_ELECTRON_TEMP(
     >                 J,    !.. point on the field line
     >               ILJ,    !.. used to print diagnostics if JSJ NE 0 or 1
     >                ID,    !.. unit to print diagnostics
     >                DT,    !.. Time step in seconds
     >                 N,    !.. O+, H+, He+, minor ion densities array
     >                TI,    !.. Ion and electron temperatures
     >             TISAV,    !.. temperatures at previous time step
     >                 F,mp,lp)    !.. OUTPUT: The integrated continuity equation
      USE THERMOSPHERE         !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE FIELD_LINE_GRID      !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL          !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER mp,lp
      INTEGER J,ILJ,ID,I,K
      !.. I/O parameters, see comments above
      DOUBLE PRECISION DT,N(4,FLDIM),TI(3,FLDIM),TISAV(3,FLDIM),F(20)
      DOUBLE PRECISION BOLTZ,NE             !.. Boltzmann constant, e density
      DOUBLE PRECISION DS_UP,DS_LO          !.. distances between grid points
      DOUBLE PRECISION CE_ION(3),CE_NEUT(3) !.. e-i heat transfer, and neut cool
      DOUBLE PRECISION EHEAT(3)             !.. Electron heating
      DOUBLE PRECISION KE_UP,KE_LO,KE(3)    !.. Thermal conductivity coeffs
      DOUBLE PRECISION B_UP,B_LO            !.. Magnetic field at half points
      DOUBLE PRECISION HFLUX_UP,HFLUX_LO    !.. Heat flux components
      DOUBLE PRECISION DTE(3),DERIVT        !.. dT/dt and Integrated dT/dt
      DOUBLE PRECISION ElossI       !.. Integrated e-i energy exchange rate
      DOUBLE PRECISION ElossN       !.. Integrated e-n cooling rate
      DOUBLE PRECISION EHrate       !.. Integrated photoelectron heating rate
      DOUBLE PRECISION kNT(3),DIVV_UP,DIVV_LO  !.. Div of velocity cmpts
      DOUBLE PRECISION EVEL(3),VLIMIT     !.. Velocity cmpts, limiting velocity
      DATA BOLTZ/1.3807E-16 /             !.. Boltzmann constant
      DATA VLIMIT/1.0E+6 /                !.. Velocity limit in cm/sec          

      !.. Loop to evaluate terms at J-1, J, J+1 for integration, J is 
      !.. the actual grid point. All terms divided by magnetic field (BM)
      DO I=1,3
        K=I+J-2
        !.. Electron density = total ion density
        NE=XIONN(1,K)+XIONN(2,K)+XIONN(3,K)+XIONN(4,K)+XIONN(5,K)+
     >       XIONN(6,K)+XIONN(7,K)+XIONN(8,K)+XIONN(9,K)
        !..temperature derivative dT/dt 
        DTE(I)=1.5*BOLTZ*NE*(TI(3,K)-TISAV(3,K))/BM(K)/DT
        !..electron heating rate in ergs
        EHEAT(I)=EHT(3,K)*1.6E-12/BM(K)
        !.. determine the thermal conductivity KE 
        CALL GET_CONDUCTIVITY(TI(3,K),NE,ON(K),O2N(K),N2N(K),HE(K)
     >      ,HN(K),KE(I))
        !.. ion and neutral cooling rates
        CALL GET_EcoolN(K,TI(3,K),NE,CE_NEUT(I),mp,lp)
        !.. Electron - ion energy exchange
        CALL GET_EcoolI(K,FLDIM,TI(3,K),TI(2,K),BM(K),CE_ION(I))
        !.. kNT is used in convection term div(vel)/ds coefficient
        kNT(I)=BOLTZ*NE*TI(3,K)
        !.. Average electron velocity from ion flux = flux(O+) + flux(H+)
        EVEL(I)=(XIONV(1,K)*XIONN(1,K)+XIONV(2,K)*XIONN(2,K))/
     >    (XIONN(1,K)+XIONN(2,K))
        IF(EVEL(I).GT.VLIMIT) EVEL(I)=VLIMIT
        IF(EVEL(I).LT.-VLIMIT) EVEL(I)=-VLIMIT
      ENDDO

      !.. distances along field line between grid points for the integration 
      DS_UP=SL(J+1)-SL(J)
      DS_LO=SL(J)-SL(J-1)

      !..Now perform the integration to get components of the energy equation
      !..DERIVT, ElossI, ElossN, EHrate
      CALL TERLIN(DTE(1),DTE(2),DTE(3),DS_LO,DS_UP,DERIVT)
      CALL TERLIN(CE_ION(1),CE_ION(2),CE_ION(3),DS_LO,DS_UP,ElossI)
      CALL TERLIN(CE_NEUT(1),CE_NEUT(2),CE_NEUT(3),DS_LO,DS_UP,ElossN)
      CALL TERLIN(EHEAT(1),EHEAT(2),EHEAT(3),DS_LO,DS_UP,EHrate)

      B_UP=(BM(J)+BM(J+1))*0.5    !.. Upper magnetic field
      B_LO=(BM(J)+BM(J-1))*0.5    !.. Lower magnetic field

      !.. Upper and lower ion thermal conductivity coefficient
      KE_UP=(KE(3)+KE(2))/2.0      !.. Upper conductivity coefficient
      KE_LO=(KE(2)+KE(1))/2.0      !.. lower conductivity coefficient

      !.. Upper and lower ion thermal conduction heat flux
      HFLUX_UP=KE_UP*(TI(3,J+1)-TI(3,J))/DS_UP/B_UP
      HFLUX_LO=KE_LO*(TI(3,J)-TI(3,J-1))/DS_LO/B_LO

      !.. Calculate the divergence of the velocity term (integrated). 
      !.. DIVV_UP and DIVV_LO are the values at the upper and lower .
      !.. bounds of the integral.
      !.. 2nd last term on RHS in Bailey, PSS 1983 equation 15
      !DIVV_UP=0.5*(kNT(3)+kNT(2))*(EVEL(3)-EVEL(2))/B_UP
      !DIVV_LO=0.5*(kNT(2)+kNT(1))*(EVEL(2)-EVEL(1))/B_LO

      DIVV_UP=0.00
      DIVV_LO=0.00

      !.. Integrated energy equation for ions
      F(2)=DERIVT+ElossI+ElossN-EHrate- HFLUX_UP+HFLUX_LO
     >   +DIVV_UP-DIVV_LO  !.. div.V term (small)
 
      !..This section for diagnostic print
      IF(ILJ.EQ.3) THEN
        WRITE(ID,501)
        WRITE(ID,5555) Z(J),TI(3,J),TI(2,J),TN(J),F(2),DERIVT,ElossI,
     >    ElossN,EHrate,HFLUX_UP,HFLUX_LO,DIVV_UP,DIVV_LO,
     >    HFLUX_UP-HFLUX_LO,DIVV_UP-DIVV_LO
      ENDIF
 501  FORMAT(' Te: ALT      Te       Ti       Tn        F      dT/dt'
     > ,6X,'C_ei      C_en     EHrate   HFLUX_UP  HFLUX_LO'
     > ,3X,'DIVV_UP   DIVV_LO  DIVCOND  DIVCONV')
 5555 FORMAT(4F9.1,1P,22E10.2)
      RETURN
      END
C:::::::::::::::::::::::::::::: GET_EcoolN ::::::::::::::::::::::::::::::::::::::::::::
C... Calculates the ion and neutral cooling rates for electrons
C... Rees and Roble, Rev. Geophys.(1975) P220 
C... Schunk and Nagy Rev Geophys (1978) p366
      SUBROUTINE GET_EcoolN(K,   !.. actual point on the field line
     >                     TE,   !.. Ion and electron temperatures
     >                     NE,   !.. Electron density
     >                CE_NEUT,mp,lp)   !.. Cooling rate to neutrals
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE       !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE ION_DEN_VEL        !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER mp,lp
      INTEGER K                        !.. actual point on the field line
      DOUBLE PRECISION TE,NE,CE_NEUT   !.. I/O parameters, see above
      DOUBLE PRECISION SQTE,TDIF,TDIFRT  !.. Computer efficiency parameters
      !.. Parameters for elec-neut elas coll loss rates
      DOUBLE PRECISION LEN2,LEO2,LEO,LEH,LRN2,LRO2,LVN2,LVO2
      !.. Parameters for O fine structure cooling (LFO)
      DOUBLE PRECISION D1,D2,E1,E2,E3,LF1,LF2,LF3,ZFO,LFO
      DOUBLE PRECISION EF,GE,HS      !.. factors for N2(v) and O2(v)
      !.. Parameters for O'(1D) cooling
      DOUBLE PRECISION DE,EXH,EXH2,LF1D,KEN
	
        !.. O2 density is used to ignore loss to neutrals at high altitudes
        IF(O2N(K).LT.1.0E-2) THEN
           CE_NEUT=0.0
           RETURN
        ENDIF

        !.. These terms are computed here for efficiency because
        !.. they are used in many places
        SQTE=SQRT(TE)
        TDIF=TE-TN(K)
        TDIFRT=TDIF/TE/TN(K)

        !..  Electron-neutral elastic cooling rates S&N p366
        LEN2=1.77E-19*NE*N2N(K)*(1.0-1.2E-4*TE)*TE*TDIF
        LEO2=1.21E-18*NE*O2N(K)*(1.0+3.6E-2*SQTE)*SQTE*TDIF
        LEO=7.9E-19*NE*ON(K)*(1.0+5.7E-4*TE)*SQTE*TDIF
        LEH=9.63E-16*NE*HN(K)*(1.0-1.35E-4*TE)*SQTE*TDIF

        !.. rotational loss rates B&K 268 ####
        LRN2=2.0E-14*NE*N2N(K)*TDIF/SQTE
        LRO2=7.0E-14*NE*O2N(K)*TDIF/SQTE

        !.. N2 vib loss rates B&K p268 and s&n p364  +++++
        EF=1.06E+4+7.51E+3*TANH(1.1E-3*(TE-1800.0))
        GE=3300.0+1.233*(TE-1000.0)-2.056E-4*(TE-1000.0)*
     >    (TE-4000.0)
        LVN2=-2.99E-12*NE*N2N(K)*EXP(5.0E-4/TE*EF*(TE-2000.0))
     >    *(EXP(-GE*TDIFRT)-1.0)

        !..O2 vib loss rate S&N p364  $$$$
        HS=3300.0-839.0*SIN(1.91E-4*(TE-2700.0))
        LVO2=-5.196E-13*NE*O2N(K)*EXP(1.4286E-3/TE*HS*(TE-700.0))
     >        *(EXP(-2770.0*TDIFRT)-1.0)

        !..fine structure cooling by atomic oxygen
        D1=EXP(-228.0/TN(K))
        D2=EXP(-326.0/TN(K))
        E1=EXP(-228.0/TE)
        E2=EXP(-326.0/TE)
        E3=EXP(-98.0/TE)*D1
        LF1=8.49E-6*TE**0.519*(0.02*(D1-E1)-5.91E-9*
     >     TDIF*(2.019*D1+(228.0/TE+2.019)*E1))
        LF2=7.7E-6*TE**0.3998*(0.028*(D2-E2)-5.91E-9*
     >     TDIF*(1.8998*D2+(326.0/TE+1.8998)*E2))
        LF3=2.22E-7*TE**0.768*(0.008*(D2-E3)-5.91E-9*
     >     TDIF*(2.268*D2+(98.0/TE+2.268)*E3))
        ZFO=5.0+3.0*D1+D2
        LFO=-8.629E-6*NE*ON(K)*(LF1+LF2+LF3)/ZFO
 
        !.. Cooling to O(1D)  S&N P365. 
        !.. EXH and EXH2 are to avoid underflow problems on the 
        !.. MSFC IBM computer. Probably no longer needed
        DE=2.4E4+(TE-1500.0)*(0.3-1.947E-5*(TE-4000.0))
        EXH=3.3333E-4*DE*(TE-3000.0)/TE
        IF(EXH.GT.70.0) EXH=70.0
        EXH2=22713.0*TDIFRT
        IF(EXH2.GT.70.0) EXH2=70.0
! GHGM
!       if((mp.eq.65).and.(lp.eq.30)) then
!       write(6,88) k,ne,on(k),exh,exh2,te,tn(k)
!88     format('GHGM ECOOL N ', i6,6e12.4)
!       endif
! GHGM
        LF1D=-1.57E-12*NE*ON(K)*EXP(EXH)*(EXP(-EXH2)-1.0)

 14     KEN=LEN2+LEO2+LEO+LEH+LRN2+LRO2
        CE_NEUT=(KEN+LVN2+LVO2+LFO+LF1D)*1.6E-12/BM(K)
       RETURN
       END

C:::::::::::::::::::::::::: GET_CONDUCTIVITY ::::::::::::::::::::::::::::::::::::::::::::::::
C... calculates the thermal conductivity co-efficient (KE)
C... Schunk and Nagy Rev Geophys (1978) p366
      SUBROUTINE GET_CONDUCTIVITY(TE, !.. Electron temperature
     >                            NE, !.. Electron density
     >                            ON, !.. Atomic oxygen density
     >                           O2N, !.. Molecular oxygen density
     >                           N2N, !.. Nitrogen density
     >                            HE, !.. Helium density
     >                            HN, !.. Hydrogen density
     >                            KE) !.. OUTPUT: Conductivity  
      DOUBLE PRECISION TE,NE,ON,O2N,N2N,HE,HN,KE !.. see I/O comments
      !.. Factors for evaluating KN
      DOUBLE PRECISION SQTE,KN2,KO2,KO,KHE,KH,KN

        !.. Evaluate factor for neutral effects on conductivity (KN).
        !.. Banks Ann. Geophys., 22, 577, [1966]. Also Schunk & Nagy eqn 59
        SQTE=SQRT(TE)   !.. saves evaluating this 3 times
        KN2=(2.82E-17*SQTE-3.41E-21*SQTE*TE)*N2N
        KO2=(2.2E-16+7.92E-18*SQTE)*O2N
        KO=3.4E-16*ON
        KHE=5.6E-16*HE
        KH=(5.47E-15-7.45E-19*TE)*HN
        !.. 1.0E+4 added to NE to avoid numerical problems for low NE
        KN=3.22E+4*TE**2*(KN2+KO2+KO+KHE+KH)/(NE+1.0E+4)
        KE=1.232E-6*(TE**2.5)/(1.0+KN)   !-- Thermal conductivity 
      RETURN
      END
