C.................... ELECXS.FOR;2 ..........30-JUL-1993 15:41:37.35 
C----- This file contains the electron impact cross sections and related
C----- subroutines for both the auroral model and the FLIP model. It 
C----- replaces AURXS.FOR for the auroral model and R2SSUB.FOR for the FLIP
C----- model. Some of the routines not connected with cross sections on 
C----- R2SSUB.FOR were appended to RSPE2B.FOR. These changes were made
C----- by Phil Richards on 23 FEB 1993.
C
C----- Glynn, 
C----- the most significant change in this program is the change of
C----- the name of PRN2XS to EPN2XS. You will need to change the CALL
C----- statement. I noticed that there is a variable by that name and did 
C----- not want the confusion. 
C----- I also updated the O(1D) and the 1356 cross sections.
C----- In EPRNXS, SIG(1,12) was erroneously going to PEXCIT(3,3,J) instead of
C----- PEXCIT(3,2,J).
C----- I'm not sure what you are doing with the two NI subroutines at the
C----- bottom of this file. They are not being called at present. Mt NI
C----- routines have recently been updated.
C
C:::::::::::::::::::::::::: ELASTC ::::::::::::::::::::::::::::::::::::::::
C..... This subroutine returns the elastic electron impact cross sections
C..... (SIGEL) and backscatter ratios (PEBSC) for electron impact on O,
C..... O2, and N2 from SOLOMON ET AL JGR p9867, 1988
C..... 2008-01-07, Cross sections updated to better agree with Solomon et al
      SUBROUTINE ELASTC(ELECEN,PEBSC,SIGEL)
      DIMENSION PEBSC(3),SIGEL(3),SGELN2(10)
C
C....... Elastic cross sections from Solomon et al. JGR, 9867, 1988
      IF(ELECEN.LT.10.0) THEN
         SIGEL(1)=5.3E-16*ELECEN**(0.146)
         SIGEL(2)=5.0E-16*ELECEN**(0.301)
	   SIGEL(3)=1.0E-15+1.5E-15*EXP(-((ELECEN-2.2)/0.8)**2)
      ELSEIF(ELECEN.LT.1000.) THEN
         SIGEL(1)=6.0E-15*(1-0.1/ELECEN)**20.0*ELECEN**(-0.68)
         SIGEL(2)=1.1E-14*(1-0.1/ELECEN)**70.0*ELECEN**(-0.68)
         SIGEL(3)=SIGEL(2)
         IF(SIGEL(1).GT.7.5E-16) SIGEL(1)=7.5E-16
         IF(SIGEL(3).GT.1.0E-15) SIGEL(3)=1.0E-15
      ELSE
         SIGEL(1)=3.1E-14*(1-0.1/ELECEN)**70.0*ELECEN**(-0.93)
         SIGEL(2)=6.2E-14*(1-0.1/ELECEN)**70.0*ELECEN**(-0.93)
         SIGEL(3)=SIGEL(2)
         IF(SIGEL(3).GT.1.0E-15) SIGEL(3)=1.0E-15

      ENDIF
      !.. backscatter coefficient
      IF(ELECEN.LT.1000) THEN
        PEBSC(1)=1.36*(1-0.1/ELECEN)**4.0*ELECEN**(-0.50)
        PEBSC(2)=1.65*(1-0.1/ELECEN)**4.0*ELECEN**(-0.48)
      ELSE
        PEBSC(1)=21.5*(1-0.1/ELECEN)**4.0*ELECEN**(-0.90)
        PEBSC(2)=21.3*(1-0.1/ELECEN)**4.0*ELECEN**(-0.85)
      ENDIF
      PEBSC(3)=PEBSC(2)  !.. N2 same as O2

      !..WRITE(6,291) ELECEN,SIGEL(1),SIGEL(3),PEBSC(1),PEBSC(3)

 291  FORMAT(2X,' Elastic XS ',F10.2,1P,9E10.2)
C......... To avoid dividing by zero  in 2-stream program
      DO 55 IS=1,3
         IF(PEBSC(IS).GT.0.5) PEBSC(IS)=0.5
         IF(SIGEL(IS).LT.1.0E-25) SIGEL(IS)=1.0E-25
 55   CONTINUE

      RETURN
      END
C:::::::::::::::::::::::::SIGEXS ::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SIGEXS(E,SIGEX,SIGION,SIGO1D)
C..... program for evaluating the total inelastic cross sections
      DIMENSION SIGEX(3),SIGION(3),PARSIG(22)
C......... factors for adjusting O, N2 excitation cross sections
      DATA XSFACO,XSFACN,XSFO2/1.0,1.0,1.0/
C
C.......... O cross sections
      CALL OXSIGS(E,PARSIG,TSIG)
      SIGO1D=PARSIG(1)
C.......... total cross section for O
      SIGEX(1)=PARSIG(1)+
     > (PARSIG(2)+PARSIG(3)+PARSIG(4)+PARSIG(5)+PARSIG(6))
     > * XSFACO
C
C...... total excitation cross section for o2:
C...... Pitchford and Phelps (Stamnes and Rees, JGR p6301, 1983)
C...      IF(E.LE.10.0) SIGEX(2)=2.2E-20*E**3.74
C...      IF(E.GT.10.0) SIGEX(2)=3.79E-16/SQRT(E)
C....... Solomon et al. 1988 JGR p9867
      IF(E.LT.10.0) SIGEX(2)=2.0514E-19*E**2.29
      IF(E.GE.10.0) SIGEX(2)=3.0476E-15*(1-5.0/E)**3.75 *E**(-0.75)
      SIGEX(2)=XSFO2*SIGEX(2)
C
C...... total excitation cross section for n2...... from ae fluxes
      IF(E.LE.15) SIGEX(3)=(1.94*E-14.14)*1.0E-17
      IF(E.LT.5.0) SIGEX(3)=5.0E-9*(1.-1.4/E)**9*(1.4/E)**16
      IF(E.GE.4.AND.E.LT.7) SIGEX(3)=1.0E-19
      IF(E.GT.15.0) SIGEX(3)=8.4711E-16*(1-7.0/E)**0.6/SQRT(E)
      IF(SIGEX(3).LE.0.0) SIGEX(3)=1.0E-19
      SIGEX(3)=XSFACN*SIGEX(3)
C........... Ionization cross sections
      CALL TXSION(E,SIGION)
      RETURN
      END
C:::::::::::::::::::::: OXSIGS :::::::::::::::::::::::::::::::::::::
      SUBROUTINE OXSIGS(E,SIGEX,SIGEXT)
C....... Inelastic cross sections for electron impact on atomic oxygen
C....... E=electron energy, SIGEX(22)=array of partial cross sections,
C....... SIGEXT=total excitation cross section, and S
      DIMENSION SIGEX(22),SO1D(7)
      DO 9 I=1,22
 9    SIGEX(I)=0.0
C------ CROSS SECTION FOR O(1D) - New Doering cross section from JGR
C------ p19531, 1992. Increases production by a factor of 1.13
      DATA SO1D/0.0,0.0,15.0,30.0,44.0,54.0,38.0/
C..      IF(E.GT.6.5) SIGEX(1)=3.77E-15*(1-2.0/E)**4.25/E**1.7
C..      IF(E.LE.7.3) SIGEX(1)=SO1D(NINT(E))*1.0E-18
C
C..... Old cross section of Henry
      IF(E.GT.1.96) SIGEX(1)=4E-16*(1-1.96/E)**2/E
C........ O(1S) cross section: may be double Shyn et al. JGR 1986, 13751
      IF(E.GT.4.17) SIGEX(2)=6.54E-17*(1-SQRT(4.17/E))/E
C....... 1304, 1027 A, Zipf and Erdman JGR 1985, 11088 include cascade.
C....... Direct excitation is half for  1304 (Vaughan and Doering,
C........ JGR 1986, 13755 and 1987 in press)
      IF(E.GE.10) SIGEX(3)=67.6E-17*(E-10)/E**2
C....... 989 cross section from Doering 1987 (1/2 of Zipf)
      IF(E.GE.14) SIGEX(4)=7.0E-17*(1-14/E)/SQRT(E)
      SIGEX(5)=0.38*SIGEX(4)
C....... O(5S) 1356 A Stone And Zipf Corrected By Zipf And Erdman 1985
C------ reparameterized 1 May 92 using FITXS.FOR (PGR)
      IF(E.GT.10.0) SIGEX(6)=4.867E-12*(1.0-9.0/E)**2.67/ E**4.0
      SIGEXT=SIGEX(1)+(SIGEX(2)+SIGEX(3)+SIGEX(4)+SIGEX(5)+SIGEX(6))
      RETURN
      END
C::::::::::::::::::::::: TXSION ::::::::::::::::::::::::::::::::::
C..... total ionization cross sections for O, O2, and N2
C..... ionization cross sections keiffer and dunn ........
C..... The N2+ and O2+ cross sections were modified in April 99 to agree
C..... with the Schram et al. cross sections at high energies
      SUBROUTINE TXSION(E,SIGIT)
      DIMENSION SIGIT(3)

      !... SIGTMP is used for N2+ and O2+ at the high energies
      SIGTMP=1.0E-13*EXP(-2.303*ALOG10(E))

      !... N2+ cross section
      SIGIT(3)=0.0
      IF(E.GT.15.0) SIGIT(3)=1.42E-14*(1-9.0/E)**7.1*E**(-0.7) 
      IF(SIGTMP.LT.SIGIT(3)) SIGIT(3)=SIGTMP
      !... This correction to convert units to cm**2. Keiffer and Dunn page 10
      SIGIT(3)=0.87972*SIGIT(3)

      !... O2+ cross section
      SIGIT(2)=0.0
      IF(E.GT.12.0) SIGIT(2)=1.08E-14*(1-7.0/E)**8.6*E**(-0.65)
      IF(SIGTMP.LT.SIGIT(2)) SIGIT(2)=SIGTMP
      !... This correction to convert units to cm**2. Keiffer and Dunn page 10
      SIGIT(2)=0.87972*SIGIT(2)

      !... O+ cross section from Brook et al. J. Phys. B. Vol 11 p 3115, 1978
      SIGIT(1)=0.0
      IF(E.GT.12.0) SIGIT(1)=7.33E-15*(1-2.0/E)**34.3*E**(-0.7)
      RETURN
      END
C:::::::::::::::::::::::::::::::: EPN2XS ::::::::::::::::::::::::::::::::::::
C........ these programs evaluate the photoelectron excitation cross
C........ sections using the parameterisation of Porter et al j. chem. phys.
C........ 1976, p154 (pjg). and jackman et al jgr 1977, p5081 (jgg)
C........ they also evaluate the heat absorption cross sections(sig(2,j))
      SUBROUTINE EPN2XS(J,FLDIM,XN,PHISUM,ALT,EP,SIGEX,SIGION,PEXCIT)
      INTEGER J,IKS,IS,FLDIM
      DIMENSION SIGO(22),SIG(2,20),XN(3,FLDIM),SIGEX(3),SIGION(3)
      DIMENSION PEXCIT(3,12,FLDIM)
      SAVE IKS,ESAVE,SIGO,SIG
      !-- common to hold the EUV and photoelectron production rates 
      DATA ESAVE/0/

      !.. ESAVE is used avoid recalculating cross sections when electron
      !.. energy is the same.
      IF(ABS(EP-ESAVE).GT.1.0E-5) THEN
         CALL N2FOB(EP,SIG,IKS)
         CALL N2OTH(EP,SIGO)
         !.. normalize partial cross sections to total
         SIGFAC=1.0
	   IF(SIG(1,IKS+1).GT.1.0E-18) SIGFAC=SIGEX(3)/SIG(1,IKS+1)
         SIGFAC=1.0    !.. Don't normalize cross sections
	   DO IS=1,IKS
	   SIG(1,IS)=SIGFAC*SIG(1,IS)
	ENDDO
      ENDIF
      ESAVE=EP

      PHIBYN=PHISUM*XN(3,J)
C
C---- excitation rates for the A,B,C,a1, N(4S),N(2D) and N(2P) states of N2
C---- the fractions represent percentage dissociation from Porter et al
C--- Note that there is some uncertainty over the percentage of the B' state
C--- SIG(1,12) that goes to the B state. Some goes to the ground state
      PEXCIT(3,1,J)=PEXCIT(3,1,J)+SIG(1,3)*PHIBYN
      PEXCIT(3,2,J)=PEXCIT(3,2,J)+(SIG(1,4)+SIG(1,11)+SIG(1,12))*PHIBYN
      PEXCIT(3,3,J)=PEXCIT(3,3,J)+(0.5*SIG(1,5))*PHIBYN
      PEXCIT(3,4,J)=PEXCIT(3,4,J)+1.0*SIG(1,7)*PHIBYN
C....... apportion production rates for states that dissociate. note that
C....... fractions are doubled because 2 atoms produced
      PEXCIT(3,5,J)=PEXCIT(3,5,J)
     >   +(SIG(1,5)+0.4*SIG(1,7)+SIG(1,9)+SIG(1,10))*PHIBYN
      PEXCIT(3,6,J)=PEXCIT(3,6,J)+(SIG(1,9)+SIG(1,10))*PHIBYN
C...... store production of n2 vib use 
      IF(EP.GE.1.0.AND.EP.LE.5) 
     >    PEXCIT(3,8,J)=PEXCIT(3,8,J)+SIGEX(3)*PHIBYN
C....... add in dissociation production rates from Rydbergs for n2. 
C....... The Rydberg cross section is assumed to be 0.1 of ionization 
C....... cross section
      SIGRYD=0.2*SIGION(3)
      PEXCIT(3,5,J)=PEXCIT(3,5,J)+0.8*SIGRYD*PHIBYN
      PEXCIT(3,6,J)=PEXCIT(3,6,J)+0.8*SIGRYD*PHIBYN
      PEXCIT(3,7,J)=PEXCIT(3,7,J)+0.4*SIGRYD*PHIBYN
C
C......... 1200a production rates
      PEXCIT(3,12,J)=PEXCIT(3,12,J)+SIGO(1)*PHIBYN
      PEXCIT(3,11,J)=PEXCIT(3,11,J)+SIGO(2)*PHISUM
C......... 1493a production rates
      PEXCIT(3,10,J)=PEXCIT(3,10,J)+SIGO(3)*PHIBYN
C
      SIGIT=SIG(1,IKS+1)+SIGRYD
      IF(EP.LT.5) SIGIT=SIGIT+SIGEX(3)
      RETURN
      END
C::::::::::::::::::::::::::::::: N2FOB ::::::::::::::::::::::::::::::::::::::
      SUBROUTINE N2FOB(EP,SIG,IKS)
C...... this routine combines forbidden and allowed transitions for n2
C...... from porter et al 1976 with some from cartwright et al 1978
C...... note that pot for n2 vib is set to 0.99 to conserve energy
C...... may not be a bad assumption because e may excite more than 1 quanta
C...... the species are in the following order:- v1,v2,A3,B3pi,C3
C...... ,E3,a1pi,a1sigg,b1pi,b1sig,W3,B3sig,a1sigu,w1
C...... note that the e3 and a1sigg states are omitted 
C
      DIMENSION W(15),ALP(15),BET(15),OMEG(15),F(15),SIG(2,20),C(15)
     >  ,POT(15),ISW(15)
      DATA W/1.85,2.15,8.0,8.5,11.05,11.9,8.6,12.25,12.5,13.3,4*8.,0/
     >     POT/0.99,0.99,6.14,7.3,11.03,11 .9,8.5,12.3,12.3,12.7,4*8.,0/
     >     ALP/3*1.0,3*3.0,1.0,2.3,1.24,1.29,5*0/
     >     BET/8*1.0,3.66,3.72,5*0/
     >     OMEG/7.0,9.0,4*3.0,2*1.0,2*0.0,5*0/
     >     F/.273,.241,.226,.178,.28,.048,.136,.027,.666,.321,5*0/
     >     C/8*0.0,.056,.061,5*0/
     >     ISW/0,0,1,2,-2,0,-3,0,-1,-1,3,4,5,7,0/
      IKS=14
      IKSP1=IKS+1
      SIG(1,IKSP1)=0.0

      DO 25 IK=1,IKS
        SIG(1,IK)=0.0
        IF(IK.LE.2) GO TO 20
        IF(IK.EQ.6) GO TO 20
        IF(IK.EQ.8) GO TO 20
        IF(EP.LT.POT(IK)) GO TO 20
        IF(IK.LE.2.AND.EP.LE.W(IK)) GO TO 12
        IF(IK.GT.2.AND.EP.LT.W(IK)) GO TO 20
        IF(ISW(IK).EQ.-1)
     >  CALL PJGXS(EP,W(IK),F(IK),ALP(IK),BET(IK),SIG(1,IK),C(IK))
        IF(ISW(IK).EQ.0)
     >  CALL PJGXS2(EP,W(IK),F(IK),ALP(IK),BET(IK),SIG(1,IK),OMEG(IK))
        IF(ISW(IK).GT.0) CALL CARTXS(EP,ISW(IK),SIG(1,IK))
        !.. C3PI cross section from Imami and Borst 1971
        IF(ISW(IK).EQ.-2) CALL C3PIXS(EP,SIG(1,IK))
        IF(ISW(IK).EQ.-3) CALL A1PIXS(EP,SIG(1,IK))
        GO TO 20
        !.. the n2 vib and a and b states are interpolated to threshold
        !.. because the jackman et al formula fails
 12     IF(IK.LE.2) SIG(1,IK)=2.0E-17
 20     CONTINUE
        IF(SIG(1,IK).LT.0) SIG(1,IK)=0.0
        SIG(1,IKSP1)=SIG(1,IKSP1)+SIG(1,IK)
        SIG(2,IK)=SIG(1,IK)*POT(IK)
 25   CONTINUE

      RETURN
      END
C::::::::::::::::::::::::::::: PJGXS :::::::::::::::::::::::::::::::::::
      SUBROUTINE PJGXS(EP,W,F,ALP,BET,SIG,C)
      S=6.513E-14*F*(1-(W/EP)**ALP)**BET
      SIG=(S/EP/W)*ALOG(4*EP*C/W+2.7183)
      RETURN
      END
C::::::::::::::::::::::::::::::: PJGXS2 ::::::::::::::::::::::::::::::
      SUBROUTINE PJGXS2(EP,W,F,ALP,BET,SIG,OMEG)
      S=6.513E-14*F*(1-(W/EP)**ALP)**BET
      SIG=(S/W/W)*(W/EP)**OMEG
      RETURN
      END
C:::::::::::::::::::::::::::::: CARTXS :::::::::::::::::::::::::::
      SUBROUTINE CARTXS(EP,IK,SIG)
      DIMENSION A3SIG(5),B3PI(5),W3DEL(5),B3SIG(5),A1SIGU(5)
     > ,A1PI(5),W1DEL(5),C3PI(5),S(5,8),C3PISG(20)
      EQUIVALENCE (A3SIG,S(1,1)),(B3PI,S(1,2)),(W3DEL,S(1,3))
     > ,(B3SIG,S(1,4))
     > ,(A1SIGU,S(1,5)),(A1PI,S(1,6)),(W1DEL,S(1,7)),(C3PI,S(1,8))
      DATA A3SIG/2.374E-05,-1.532E-03,3.345E-02,-2.804E-01,8.697E-01/
      DATA B3PI/-5.737E-05,4.282E-03,-1.160E-01,1.336E+00,-5.228E+00/
      DATA W3DEL/8.831E-05,-5.868E-03,1.357E-01,-1.271E+00,4.221E+00/
      DATA B3SIG/8.985E-06,-4.608E-04,6.218E-03,2.410E-03,-2.324E-01/
      DATA A1SIGU/-6.431E-06,5.761E-04,-1.853E-02,2.491E-01,-1.109E+00/
      DATA A1PI/4.904E-05,-3.291E-03,7.662E-02,-7.073E-01,2.256E+00/
      DATA W1DEL/-3.221E-05,2.438E-03,-6.709E-02,7.872E-01,-3.214E+00/
      DATA C3PI/-4.431E-05,3.549E-03,-1.040E-01,1.305E+00,-5.609E+00/
C
      SIG=0
      IF(EP.LT.6.4) RETURN
      EP2=EP*EP
      EP3=EP2*EP
      EP4=EP3*EP
C......... polynomial fit to cartwright cross sections
      IF(EP.LE.26)
     > SIG=S(1,IK)*EP4+S(2,IK)*EP3+S(3,IK)*EP2+S(4,IK)*EP+S(5,IK)
      IF(SIG.LT.0.0) SIG=0.0
C......... blend in 1/e**2 dependenCe above 17 ev
      IF(EP.GT.17) ABC=EXP(0.9*(17-EP))
      IF(SIG.LE.0) ABC=0.0
      IF(EP.GT.17)
     > S17=S(1,IK)*17**4+S(2,IK)*17**3+S(3,IK)*17**2+S(4,IK)*17+S(5,IK)
      IF(EP.GT.17) SIG=(S17*(17/EP)**2+ABC*SIG)/(1+ABC)
C...... Cartwright C3 cross sections at 1 eV Intervals below 20 eV
      DATA C3PISG/11*0,.146,.289,.443,.389,.284,.234,.202,.181,.165/
      IF(IK.EQ.8.AND.EP.LE.20) SIG=C3PISG(NINT(EP))
C.......... note the 6% reduction in Cartwright's Cross sections
C.... Trajmar et al. 1983 Phys. Reports
      SIG=0.94*SIG*1.0E-16
      RETURN
      END
C:::::::::::::::::::::::::::::: C3PIXS ::::::::::::::::::::::::::::::
      SUBROUTINE C3PIXS(EP,SIG)
C....... The second positive cross section from Imami and Borst
C....... J. Chem. Phys. p1115, 1974
      DIMENSION C3PISG(20)
C
C...... Imami and Borst C3 cross sections at 1 eV Intervals below 20 eV
      DATA C3PISG/11*0,.01,.17,0.4,0.44,0.38,0.32,0.27,0.22,0.20/
      IF(EP.LE.20) THEN
          SIG=1.0E-16*C3PISG(NINT(EP))
      ELSE
          SIG=1.4E-14*(EP)**(-2.2)
      ENDIF
      RETURN
      END
C::::::::::::::::::::::::::::: A1PIXS :::::::::::::::::::::::::::::::::
      SUBROUTINE A1PIXS(EP,SIG)
C....... The LBH cross section from Ajello and Shemansky, JGR, 1985
C....... page 9845. See table 5a page 9855
      SIG=0.0
      IF(EP.GT.11)  SIG=9.4E-15*(1.0-11.0/EP)**1.2 *EP**(-1.6)
      RETURN
      END
C:::::::::::::::::::::::: N2OTH :::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE N2OTH(EP,SIGO)
C........ Cross sections for N (1200, 1394) and minor N2 emissions
      DIMENSION SIGO(22)
      DO 10 I=1,22
 10   SIGO(I)=0.0
C
C......  N 1200A Strickland and Meier value replaced by Ajello value
C......  in February 1993. PGR
C...      IF(EP.GT.25.0) SIGO(1)=1.0E-18*(0.115*EP-2.05)
C...      IF(EP.GT.60.0) SIGO(1)=1.42E-15*(1.0-25.0/EP)**2.6/ EP
      CALL NI1200(EP,SIGO(1))
C
C...... e* + N -> 1200A Zipf flouresence x-s (taken from Doering an Gumbel
C...... JGR (1991) page 16021). Cross section from FITXS.FOR
      IF(EP.GT.11.0) SIGO(2)=9.4E-16*(1.0-11.0/EP)**0.64/ EP**0.5
C
C......... CROSS SECTION FOR e*-N2 --> 1493A  
      IF(EP.GT.22.0) SIGO(3)=5.82E-16*(1.0-10.0/EP)**7.5/ EP
      RETURN
      END
C:::::::::::::::::::::::::::::: XSO2N2 :::::::::::::::::::::::::::::::::
      SUBROUTINE XSO2N2(EJ,IS,IK,SIG)
C...... THIS PROGRAM USES THE METHOD OF GREEN AND SAWADA JATP 1972,
C...... P1719, AND JACKMAN ET AL JGR 1977,P5081 TO CALCULATE PHOTOELECTRON
C...... IMPACT IONIZATION CROSS SECTIONS FOR N2 AND O2 AND O
      DIMENSION POT(3,7),RK(3,7),RKB(3,7),RJ(3,7),RJB(3,7),RJC(3,7),
     > GS(3,7),GB(3,7),TS(3,7),TA(3,7),TB(3,7),SIG(3,10),INI(3)
      DATA INI/3,7,6/
C......... POT=0 FOR N2 C, D AND 40EV STATES (THEY DISSOCIATE) PORTER ET AL
      DATA POT/13.6,12.1,15.58,16.9,16.1,16.73,18.5,16.9,18.75,0,18.2,
     >        22.,0,20.,23.6,0.,23.,40.,0.,37.,25./
     >     RK/1.03,.475,2.42,1.31,1.129,1.06,0.78,1.129,0.551,0,1.01,
     >      .371,0,0.653,0.371,0,0.95,0.53,0,0.594,0/
     >     RKB,RJB,RJC/63*0/
     >     RJ/1.81,3.76,1.74,1.79,3.76,1.74,1.78,3.76,1.74,0,3.76,1.74
     >        ,0,3.76,1.74,0,3.76,1.74,0,3.76,0/
     >     GS/13.,18.5,13.8,13.,18.5,13.8,13.,18.5,13.8,0,18.5,13.8,0,
     >        18.5,13.8,0,18.5,13.8,0,18.5,0/
     >     GB/-.815,12.1,15.58,-.815,16.1,16.73,-.815,16.9,18.75,0,
     >       18.2,22.0,0,20.3,23.6,0,23.0,40.0,0,37.0,0/
     >     TS/6.41,1.86,4.71,6.41,1.86,4.71,6.41,1.86,4.71,0,1.86,4.71,
     >        0,1.86,4.71,0,1.86,4.71,0,1.86,0/
     >     TA/3450.,1000.,1000.,3450.,1000.,1000.,3450.,1000.,1000.
     >      ,12*1000.0/
     >     TB/162.,24.2,31.16,162.,32.2,33.46,162.,33.8,37.5,0,36.4,44.
     >       ,0,40.6,47.2,0,46.0,80.0,0,74.0,0/
C
      SIG(IS,IK)=0.0
      IF(IK.GT.INI(IS)) RETURN
      IF(POT(IS,IK).LE.0.0) RETURN
      IF(EJ.LT.POT(IS,IK)) RETURN
      TM=(EJ-POT(IS,IK))/2
      AC=(RK(IS,IK)/EJ+RKB(IS,IK))
      A=AC*ALOG(EJ/RJ(IS,IK)+RJB(IS,IK)+RJC(IS,IK)/EJ)
      GE=GS(IS,IK)*EJ/(EJ+GB(IS,IK))
      T0=TS(IS,IK)-TA(IS,IK)/(EJ+TB(IS,IK))
      SIG(IS,IK)=1E-16*A*GE*(ATAN((TM-T0)/GE)+ATAN(T0/GE))
      RETURN
      END
C:::::::::::::::::::::::::::::::::::: OXRAT ::::::::::::::::::::::::::::::::
      SUBROUTINE OXRAT(E,R4S,R2D,R2P)
C....... This subroutine returns the electron impact branching ratios
C....... for atomic oxygen from Burnett and Rountree Phys. Rev. A. 20
C....... 1979 page 1468
      R4S=1.0
      R2D=0.0
      R2P=0.0
      EV=E
      IF(E.GT.100.0) EV=100.0
      IF(EV.GT.17) R4S=-1.6E-3*EV+0.56
      IF(EV.GT.17) R2D=1.067E-3*EV+0.2933
      R2P=1-R4S-R2D
         IF(EV.LT.22) THEN
         R2P=0.0
         RTOT=R4S+R2D
         R4S=R4S/RTOT
         R2D=R2D/RTOT
         ENDIF
      RETURN
      END
C-GG-03DEC91-:::::::::::::: O21356 :::::::::::::::::::::::::::::::
C
C
C.... Routine to supply the e+O2->1356 cross section
C
C     References:
C
C       Original optical estimate...
C               J.M.Ajello, J.Chem.Phys., v55, p3156, 1971.
C
C       Revised to reflect non-thermal velocity distributions...
C               W.C.Wells, W.L.Borst, and E.C.Zipf,
C               Chem.Phys. Lett., v12, p288, 1971.
C
C       New radiative lifetime determinations which require
C         further revisions...
C               C.E.Johnson, Phys.Rev., vA5, p2699, 1972.
C               W.C.Wells & E.C.Zipf, EOS Trans.AGU, v53, p459, 1972.
C               W.C.Wells & E.C.Zipf, Phys.Rev., vA9, p568, 1974.
C
C       The 7774 cascade contribution was measured and found
C         to be larger than the 1356 cross section.  Work
C         still needs to be done.
C               P.W.Erdman & E.C.Zipf, J.Chem.Phys., v87, p4540, 1987.
C
C......................................................... 
C
      SUBROUTINE O21356(E,SIG)
      DATA W/14.67/
C    
      SIG=0
      IF(E.LE.W) RETURN
      IF(E.LE.70) SIG=PARAM(4.47E-16,W,E,0.251,1.796,-0.515)
      IF(E.GT.70) SIG=PARAM(5.85E-16,W,E,0.584,2.728,-0.73 )

C-GG-15FEB93- Revised cross section to reflect change in
C             1356 radiative lifetime determination since
C             original measurement.

C             (Changed from 565 microsecond to 170 microsecond.)

C             see- Johnson, 1972
C                  Wells & Zipf, 1972
C                  Wells & Zipf, 1974

      SIG = SIG * 170 / 565

      RETURN
      END
C
C-GG-03DEC91-:::::::::::::::: PARAM :::::::::::::::::::::::::::::
C
      FUNCTION PARAM(A0,W,E,ALPHA,NU,OMEGA)
      REAL NU,OMEGA
C
C	this parameterization is similar to eq.(10)
C	Banks, Chappel, Nagy JGR vol 79, 1459, 1974
C
C	This is the form used by PGR throughout most of the program.
C					-glynn germany (1/16/90)
C
      PARAM=A0*((1-(W/E)**ALPHA)**NU)*E**OMEGA
      RETURN
      END
C
C-GG-15FEB93-:::::::::::::::::::::::::: NIFOB ::::::::::::::::::::::::::::::::::::::
C
      SUBROUTINE NIFOB(EP,SIG)
C
C.... Routine to supply cross sections for minor N2 emissions
C
C.............................................................
C
      DIMENSION SIG(22)
      DO 10 I=1,22
 10   SIG(I)=0.0
C
C....  NI 1200A  from N2 dissociation
C
      CALL NI1200(EP,SIG(1))
      RETURN
      END
C
C-GG-15FEB93 ::::::::::::::::::: NI1200 :::::::::::::::::::::::::
C
      SUBROUTINE NI1200(EP,SIG)
C
      DATA W,ALPHA/20.1,1./
      DATA E1,E2 /26.5,66.5/ OMG1,OMG2 /-.485,-.737/
      DATA A01,RNU1 /6.139E-17,1.6866/ A02,RNU2 /2.381E-16,2.5630/
C
      SIG=0.0
      Y1=PARAM(A01,W,E1,ALPHA,RNU1,OMG1)
      Y2=PARAM(A02,W,E2,ALPHA,RNU2,OMG2)
      SLOPE=(Y2-Y1)/(E2-E1)
      B=Y1-SLOPE*E1
      EMID=(E1+E2)/2.5
C
      IF(EP.LE.W) RETURN
      IF(EP.LE.E1) SIG=PARAM(A01,W,EP,ALPHA,RNU1,OMG1)
      IF((EP.GT.E1).AND.(EP.LT.E2)) THEN
       SIG=SLOPE*EP+B
       IF(EP.LE.EMID) THEN
      	C=1.+0.265*(E1-EP)/(EMID-E1)
       ELSE
        C=1.-0.265*(E2-EP)/(E2-EMID)
       ENDIF
       SIG=SIG*C
      ENDIF
      IF(EP.GE.E2) SIG=PARAM(A02,W,EP,ALPHA,RNU2,OMG2)
      RETURN
      END
C:::::::::::::::::::::::::::: BACKSCAT ::::::::::::::::::::::::
      !.. This subroutine calculates the inelastic backscatter
      !.. coefficient for N2 (Solomon et al. JGR, 1988, page 9878
      !.. Note that there is little lab data on backscatter so
      !.. and all angles are folded into one for the 2 stream model
      !.. these are rough. You can use the same backscatter coeff
      !.. for O and O2. Programmed by P. Richards, July 1998
      SUBROUTINE BACKSCAT(E,        !... energy of electron
     >                    BSCAT)    !... backscatter coeff
      IMPLICIT NONE
      REAL E,BSCAT
      IF(E.GT.500.0) THEN
         BSCAT=10**(-0.76862*ALOG10(E)+0.55161)
      ELSE IF(E.GT.100.0) THEN
         BSCAT=10**(-0.31739*ALOG10(E)-0.666674)
      ELSE IF(E.GT.20.0) THEN
         BSCAT=10**(-1.43068*ALOG10(E)+1.56032)
      ELSE IF(E.GE.0.0) THEN
         BSCAT=0.5
      ENDIF
      RETURN
      END


C:::::::::::::::::::::::::::::::: XS_TEST ::::::::::::::::::::::::::::::::::::
C........ For testing electron excitation cross
      SUBROUTINE XS_TEST(EP)
      DIMENSION SIGO(22),SIG(2,20),SIGEX(3),SIGION(3)
      DATA JSAV/0/

      !.. ESAVE is used avoid recalculating cross sections when electron
      !.. energy is the same.
      !write(6,*) '  ELECXS line 555'
      CALL N2FOB(EP,SIG,IKS)
      CALL N2OTH(EP,SIGO)

	SIGFOB=0.0
	SIGOTH=0.0
      CALL SIGEXS(EP,SIGEX,SIGION,SIGO1D)

	SIGFAC=1.0  !.. normalize partial cross sections to total
	IF(SIG(1,IKS+1).GT.1.0E-18) SIGFAC=SIGEX(3)/SIG(1,IKS+1)

	DO IS=1,IKS
	  SIG(1,IS)=SIGFAC*SIG(1,IS)
	  SIGFOB=SIGFOB+SIG(1,IS)
	ENDDO

	SIGOTH=SIGO(1)+SIGO(2)
	SIGTOT=SIGFOB+SIGOTH

      IF(MOD(NINT(EP),15).EQ.0.OR.JSAV.EQ.0) WRITE(6,89) 
	JSAV=1
 89   FORMAT(9X,' E      tot1     fob      oth       A3'
     > ,5X,'B3pi      C3      E3       a1pi     a1sigg    b1pi'
     > ,5X,'b1sig      W3   B3sig   a1sigu    w1')
	WRITE(6,90) EP,SIGEX(3),SIGFOB,SIGFAC,(SIG(1,IS),IS=3,14)
	!..(SIGO(IS),IS=1,3)
 90   FORMAT(2X,F10.1,1P22E9.1)
      RETURN
	END
