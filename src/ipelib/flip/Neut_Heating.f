!nm20111118: added heating rate output
C............................... NEUT_HEATING.FOR ...................................
C..... This program calculates the neutral gas heating rate for the CTIPe model
C..... Written by P. Richards in September 2010 
C..... Modified in April 2011 to output the O2 dissociation frequency and the
C..... heating rates in the format of Roble et al JGR 1987, p8745
C..... O3, SRB, particle, and Joule heating are not included
C..... The chemistry includes additional reactions (e.g. O+ metastables)
C..... The input parameters JMIN, EHT were also added
      SUBROUTINE NEUT_HEATING(JPR,  !.. Input: Turns printing on and off
     >                       JMIN,JMAX,  !.. Input: Minimum index on the field line
     >                         IJ,  !.. Input: Altitude index
     >                        ALT,  !.. Input: Altitude (km)
     >                        RTS,  !.. Input: Reaction rates
     >                   TE,TI,TN,  !.. Input: Electron, ion, and neutral temperatures
     >            OXN,O2N,N2N,HEN,  !.. Input: O, O2, N2, and He densities (cm-3)
     >                        N4S,  !.. Input: N4S should be 0.5*MSIS N density (cm-3)
     >                         NE,  !.. Input: electron density (cm-3)
     >       OXPLUS,O2PLUS,NOPLUS,  !.. Input: O+, O2+, NO+ densities (cm-3)
     >         HPLUS,N2PLUS,NPLUS,  !.. Input: N2+ and N+ densities (cm-3)
     >                NNO,N2D,N2P,  !.. Input: NO, N(2D), and N(2P) density (cm-3)
     >              N2A,OP2D,OP2P,  !.. Input: N2(A3), O+(2D), O+(2D)
     >                    O1D,O1S,  !.. Input: O(1D), O(1S)
     >                        EHT,  !.. Input: Electron heating rate
     >                      NHEAT,  !.. OUTPUT: Total neutral heating rate
     >                     O2DISF,  !.. OUTPUT: O2 dissociation frequency
     > hrate ) !.. OUTPUT: !nm20121020
!     >                   hrate(2),  !.. 
!     >                   hrate(3),  !.. 
!     >                   hrate(4),  !.. 
!     >                   hrate(5),  !.. 
!     >                   hrate(6),  !.. 
!     >                   hrate(7))  !.. 
      USE PRODUCTION    !.. EUV, photoelectron, and auroral production, PHION
! save each component of heating rate for output
!JFM  USE module_IPE_dimension,ONLY: NLP
!JFM  USE module_FIELD_LINE_GRID_MKS,ONLY: mp_save,lp_save,JMIN_IN
!JFM >,                         JMAX_IS,MaxFluxTube,hrate_cgs_save
!nm20121020      USE module_input_parameters,ONLY: sw_neutral_heating_flip
!nm20121020     >,start_time,ip_freq_output,parallelBuild
!JFM  USE module_heating_rate,ONLY: get_neutral_heating_rate
!nm20121020      USE module_precision
!JFM  USE module_PLASMA,ONLY: utime_save
      IMPLICIT NONE
      INTEGER IJ,K          !.. loop control variables
      INTEGER JPR,JMIN,JMAX      !.. Turns on printing of production and loss
      DOUBLE PRECISION ALT              !.. Altitude (km)
      DOUBLE PRECISION RTS(99)          !.. Reaction rates
      DOUBLE PRECISION TE,TN,TI         !.. Electron and ion temperatures
	!..  H+, He+, O+, N2+, NO+, O2+, N+, 
      DOUBLE PRECISION HPLUS,OXPLUS,N2PLUS,NOPLUS,O2PLUS,NPLUS
      !.. O2,O,N2,NO,He,N4S
      DOUBLE PRECISION O2N,OXN,N2N,NNO,HEN,N4S
      !.. Ne, N(2P),N(2D),O+(2P),O+(2D), total minor ion densities
      DOUBLE PRECISION NE,N2P,N2D,OP2D,OP2P,NMINOR
      DOUBLE PRECISION DISNP            !.. Photodissociation of N2 to give N+
      DOUBLE PRECISION UVDISN           !.. Photodissociation of N2 to give N
	DOUBLE PRECISION N2A,O1D,O1S      !.. N2(A), O(1D), O(1S) densities    
      DOUBLE PRECISION HR(99)           !.. Individual heating rates
      DOUBLE PRECISION HRATE(22)        !.. Total heating rates
      DOUBLE PRECISION PO1DSR           !.. Schumann-Runge production of O(1D)
      DOUBLE PRECISION PEO1D,PO1SEL     !.. Photoelectron production of O(1D), O(1S)
      DOUBLE PRECISION PO1D,LO1D        !.. O(1D) production and loss rates
      DOUBLE PRECISION COOL(22)         !.. Electron cooling rates
      DOUBLE PRECISION FOL              !.. Electron OX fine structure cooling rate
      DOUBLE PRECISION TLSS             !.. Electron cooling rate to neutrals
      DOUBLE PRECISION EHT              !.. Electron heating rate
      DOUBLE PRECISION NHEAT            !.. Neutral heating rate
      DOUBLE PRECISION O2DISF           !.. O2 dissociation frequency
      DOUBLE PRECISION CF               !.. Converts heating rates to per unit mass
      DOUBLE PRECISION HN2D             !.. Electron heating from N2D
      
!nm20110404: save each component of the heating rate for output
!nm20121020      INTEGER :: j2d,jth,lun,lp  !J converted to the 2Dsystem in a meridional plain.

!JFM  REAL(KIND=real_prec), DIMENSION(7,MaxFluxTube,NLP) :: hrate_mks !.. each component of the Neutral heating rate (eV/kg/s) 
!nm20121020      REAL(KIND=real_prec) :: min_hrate,max_hrate

      PO1DSR=OTHPR1(3,IJ)      !.. Schumann-Runge production of O(1D)
      PEO1D=PEXCIT(1,1,IJ)     !.. Photoelectron production of O(1D)
      PO1SEL=PEXCIT(1,2,IJ)    !.. Photoelectron production of O(1S)
      UVDISN=OTHPR1(1,IJ)      !.. Dissociation rate of N2

!dbg20120206:
!dbg      print *,'sub-Neut_Heating: JPR',JPR
      !.. O2 dissociation rate frequency (s-1) See Roble et al. JGR 1987, p8745
!     IF(IJ.EQ.JMIN.AND.JPR.GT.0) WRITE(121,150)
!150  FORMAT(3X,' O2 dissociation rate frequency (s-1)'
!    > ,/3X,'ALT   S-R    N2D+O2  N4S+O2  N++O2   E+O2+   N+O2+  TOTAL')
      HR(1)=PO1DSR/O2N
      HR(2)=N2D*RTS(16)
      HR(3)=RTS(7)*N4S
      HR(4)=(RTS(22)+RTS(25)+RTS(30)+RTS(65)+RTS(66))*NPLUS
      HR(5)=RTS(6)*NE*O2PLUS/O2N
      HR(6)=RTS(21)*O2PLUS*N4S/O2N
      O2DISF=HR(1)+HR(2)+HR(3)+HR(4)+HR(5)+HR(6)
!     IF(JPR.GT.0) WRITE(121,88) ALT,(HR(K),K=1,6),O2DISF

      !.. neutral-neutral kinetic heating
!     IF(IJ.EQ.JMIN.AND.JPR.GT.0) WRITE(122,154)
!154  FORMAT(2X,' Neutral-Neutral kinetic heating'
!    > ,/2X,'ALT   N2D+O   N2D+O2   O+N2A   O+N2P   N+NO'
!    > ,4X,'N+O2   O1D+O   O1D+N2  O1D+O2   total    O1D'
!    > ,5X,'O1S     N2D     N2P     N2A     NNO     N4S     O2N')
      !.. Note that N(2D)+e is not included in total n-n heating because it 
      !.. is included in the photoelectron flux and appears as thermal
      !.. electron energy - Richards, Planet. Space Sci., 34, 689, 1996.
      HR(1)=2.38*N2D*OXN*RTS(15)
      HR(2)=1.84*N2D*O2N*RTS(16)
      HR(3)=6.14*RTS(36)*OXN*N2A   !.. PGR
      HR(4)=3.0*RTS(37)*N2P*OXN    !.. PGR
      HR(5)=1.85*RTS(9)*N4S*NNO
      HR(6)=1.42*RTS(7)*O2N*N4S

      !.. Total neutral-neutral kinetic heating
      HRATE(1)=HR(1)+HR(2)+HR(3)+HR(4)+HR(5)+HR(6)

      !.. Heating from O(1D)
      HR(7)=1.96*O1D*OXN*RTS(69)   !.. PGR
      HR(8)=1.96*O1D*N2N*RTS(33)
      HR(9)=1.96*O1D*O2N*RTS(34)
      HRATE(2)=HR(7)+HR(8)+HR(9)   !.. Total heating from O(1D)

!     IF(JPR.GT.0) WRITE(122,88) ALT,(HR(K),K=1,9),HRATE(1)+HRATE(2)
!    >  ,O1D,O1S,N2D,N2P,N2A,NNO,N4S,O2N

      !.. ion-neutral kinetic heating
!     IF(IJ.EQ.JMIN.AND.JPR.GT.0) WRITE(123,155)
!155   FORMAT(2X,'ion-neutral kinetic heating'
!    > ,/2X,'ALT   O+N2+   E+N2+   O2+N2+  E+NO+   E+O2+   N+O2+'
!    > ,3X,'NO+O2+  N2+O+   O2+O+    O2+N+   O2+N+  O+N+   O+2D+N2'
!    > ,1X,'O+2D+e  O+2P+N2 O+2P+e   N2D+O+   total')
      HR(1)=0.70*RTS(10)*OXN*N2PLUS
      HR(2)=3.44*RTS(11)*NE*N2PLUS
      HR(3)=3.53*RTS(17)*O2N*N2PLUS
      HR(4)=0.90*NOPLUS*NE*RTS(5)
      HR(5)=4.42*RTS(6)*NE*O2PLUS
      HR(6)=4.19*RTS(21)*O2PLUS*N4S
      HR(7)=2.81*RTS(23)*O2PLUS*NNO
      HR(8)=1.10*OXPLUS*N2N*RTS(3)
      HR(9)=1.55*OXPLUS*O2N*RTS(4)
      HR(10)=6.67*(RTS(22)+RTS(30))*O2N*NPLUS
      HR(11)=0.11*(RTS(25)+RTS(65)+RTS(66))*O2N*NPLUS
      HR(12)=0.93*RTS(31)*OXN*NPLUS
      HR(13)=1.33*OP2D*RTS(19)*N2N        !.. PGR
      HR(14)=3.31*OP2D*NE*RTS(12)         !.. PGR
      HR(15)=3.02*RTS(20)*N2N*OP2P        !.. PGR
      HR(16)=1.69*RTS(13)*NE*OP2P         !.. PGR
      HR(17)=1.57*RTS(29)*OXPLUS*N2D      !.. PGR
      !.. Total ion-neutral kinetic heating
      HRATE(3)=HR(1)+HR(2)+HR(3)+HR(4)+HR(5)+HR(6)+HR(7)+HR(8)+HR(9)
     > +HR(10)+HR(11)+HR(12)+HR(13)+HR(14)+HR(15)+HR(16)+HR(17)
!     IF(JPR.GT.0) WRITE(123,88) ALT,(HR(K),K=1,17),HRATE(3)

      !.. Electron cooling - heats the neutrals
      NMINOR=NOPLUS+O2PLUS+N2PLUS
      HN2D=2.38*N2D*NE*RTS(8)
      CALL EHCRATS(125,IJ,JMIN,JPR,ALT,TE,TI,TN,OXPLUS,HPLUS
     > ,NMINOR,OXN,N2N,O2N,HPLUS,EHT,HN2D,TLSS,FOL,COOL)
      !.. Electron cooling - heats the neutrals
      HRATE(4)=COOL(5)
      IF(HRATE(4).LT.1E-22) HRATE(4)=1E-22

      !.. O2 S-R dissociative kinetic heating = 
      !.. total S-R energy absorbed - energy going to
      !.. dissociation (5.08 eV) and O(1D) energy (1.96 eV)
      HRATE(5)=OTHPR1(5,IJ)-7.04*PO1DSR
	IF(HRATE(5).LT.0.0) HRATE(5)=0.0

      !.. EUV and photoelectron N2 dissociation rate to form N and N+
      DISNP= EUVION(3,4,IJ)+EUVION(3,5,IJ)+EUVION(3,6,IJ)
     >   +0.1*(PEPION(3,1,IJ)+PEPION(3,2,IJ)+PEPION(3,3,IJ))  !.. Rydberg diss       
      !.. kinetic heating from dissociation of N2 to N + N and N + N+
      !.. by EUV and photoelectrons. EUV dissociation assumed to release
	!.. 1 eV of KE and photoelectron dissociation assumed to release 2 eV
      HRATE(6)=1.0*UVDISN+2.0*DISNP

      !.. Heating from 3 body collisions> See Roble et al. JGR 1987, p8745
      HRATE(7)=9.59E-34*EXP(480/TN)*OXN*OXN*(O2N+N2N)
     >   +2.15E-34*EXP(345/TN)*OXN*O2N*(OXN+O2N)
     >   +8.82E-35*EXP(575/TN)*OXN*O2N*N2N

      !.. Make sure rates are greater than zero for logs
!JFM  DO K=1,22
      DO K=1,7
        IF(HRATE(K).LT.1.0E-22) HRATE(K)=1.0E-22
      ENDDO

      !.. Total neutral heating rate
      NHEAT=HRATE(1)+HRATE(2)+HRATE(3)+HRATE(4)+HRATE(5)+HRATE(6)+
     >   HRATE(7)

!nm20110404: save each component of heating rate for output
!JFM  IF ( sw_neutral_heating_flip==1 .AND.
!JFM &  MOD( (utime_save-start_time),ip_freq_output)==0) THEN
!JFM    if(parallelBuild) then
!JFM      print*,'sw_neutral_heating_flip=1 does not work in parallel'
!JFM      print*,'Stopping in Neut_heating'
!JFM      stop
!JFM    endif

!JFM    j2d=ij+JMIN_IN(lp_save)-1
!JFM    DO jth=1,7
!JFM       hrate_cgs_save(jth,j2d,lp_save)=hrate(jth) !!(1) PGR neu_neu
!      hrate_cgs_save(2,j2d,mp_save)=hrate(2) !!PGR O1D
!      hrate_cgs_save(3,j2d,mp_save)=hrate(3) !!PGR ion_neu
!      hrate_cgs_save(4,j2d,mp_save)=hrate(4) !!PGR elec_ion
!      hrate_cgs_save(5,j2d,mp_save)=hrate(5) !!PGR SRO2dis
!      hrate_cgs_save(6,j2d,mp_save)=hrate(6) !!PGR UVN2dis
!      hrate_cgs_save(7,j2d,mp_save)=hrate(7) !!PGR 3bod CHECK dimension
!JFM    END DO


!JFM    IF ( IJ==JMAX .AND. lp_save==NLP ) THEN
!20111118 output
!JFM      IF ( mp_save==1 ) write(5000,*)utime_save
!JFM      CALL get_neutral_heating_rate ( hrate_mks )
!JFM      DO jth=1,7
!JFM        lun=5000+jth  
!JFM        min_hrate =  huge(min_hrate)
!JFM        max_hrate = -huge(max_hrate)
!JFM        do lp=1,NLP
!JFM          min_hrate=min(min_hrate,MINVAL(
!JFM >                      hrate_mks(jth,JMIN_IN(lp):JMAX_IS(lp),lp)))
!JFM          max_hrate=max(max_hrate,MAXVAL(
!JFM >                      hrate_mks(jth,JMIN_IN(lp):JMAX_IS(lp),lp)))
!JFM        enddo
!JFM  print *,'hrate',lun,jth,mp_save,lp_save,IJ,min_hrate,max_hrate
!JFM        write(lun,*)mp_save
!JFM        do lp=1,NLP
!JFM          write(lun,*)hrate_mks(jth,JMIN_IN(lp):JMAX_IS(lp),lp)
!JFM        enddo
!JFM      END DO !jth=1,7
!JFM    END IF !( lp_save==NLP ) THEN

!JFM  END IF !( sw_neutral_heating_flip==1 ) THEN

      !.. Conversion factor for converting heating rates from eV/cm3/s to ergs/gm/s
      CF=(16*OXN+28*N2N+32*O2N)*1.6726E-24/1.6022E-12

      !.. Total heating rates. The heating per unit mass is for comparing with
      !.. Roble et al JGR 1987, p8745. Note that their SRC seems to include O(1D)
      !.. that is produced by SR dissociation.
!     IF(IJ.EQ.JMIN.AND.JPR.GT.0) WRITE(126,159)
!159   FORMAT(9X,'Total neutral heating rates and O2 dissociation freq'
!    > ,/9X,'<......................... heating rates (eV/cm3/s .......'
!    > ,'...................>|<,,,,,,,,,,,,, log heating rates per mass'
!    > ,' (ergs/gm/s) ,,,,,,,,,,,,,,,,,,,,,,,>'
!    > ,/3X,'ALT   neu_neu     O1D     ion_neu   elec_ion  SRO2dis'
!    > ,3X,'UVN2dis     3bod     Total    neu_neu     O1D     ion_neu'
!    > ,3X,'elec_ion  SRO2dis   UVN2dis    3bod     Total     O2_diss')
!     IF(JPR.GT.0) WRITE(126,'(F7.2,1P,22E10.2)') ALT,
!    >  (HRATE(K),K=1,7),NHEAT
!    > ,(DLOG10(HRATE(K)/CF),K=1,7),DLOG10(NHEAT/CF),O2DISF

!88   FORMAT(F6.1,1P,22E8.1)

      RETURN
       END
C::::::::::::::::::::::: EHCRATS :::::::::::::::::::::::::::::::::::
C.... Subroutine for printing thermal electron cooling rates
C.... This function needs to be eliminated eventually and replaced
C.... by GET_loss. It is used to print cooling rates in the Hyyyyddd
C.... file of the FLIP run and also in the FLIPRINT phase 
      SUBROUTINE EHCRATS(ID,IJ,JMIN,JPR,ALT,TE,TI,TN,OPLS,HPLUS,NMIN,
     > OX,N2,O2,HN,EHT,HN2D,TLSS,LFO,COOL)
      IMPLICIT DOUBLE PRECISION(A-H,K,L,N,O-Z)
      DIMENSION FD(9),COOL(22) !.. heating and cooling rates
C
      
!     IF(IJ.EQ.JMIN.AND.JPR.GT.0) WRITE(ID,757)
!757  FORMAT(/,8X,'Thermal electron cooling(C_) and heating(H_) rates'
!    > /,6X,'ALT  Te   Ti   Tn     C_ei      C_rN2     C_rO2    C_vN2'
!    > ,6X,'C_vO2     C_fsO     C_O1D     C_tot    H_tot     H_N2D')
        TLSS=0.0
        LFO=0.0

        !.. Don't calculate cooling rates if neutral densities too low
        IF(O2.LT.1.0E-10) RETURN     
        SQTE=SQRT(TE)
        NE=OPLS+HPLUS+NMIN
        NI=OPLS/16.0+HPLUS+NMIN/30.0

        !.. Electron- ion cooling rate KEI. The Coulomb logarithm from Itikawa
        !..  JATP 1975, page 1001. Previously assumed to be 15.0
        DEBYE = 1.602E-12 * SQRT(4.0*3.1415926*NE/(1.381E-16*TE))
        COULOG = 0.43086 + LOG(TE) - LOG(DEBYE)
        KEI=1.232E-17*NE*NI*(TE-TI)/(TE**1.5)/1.6E-12    !.. e - i cooling
        KEI= KEI*COULOG/15.0                     !.. Itikawa correction

        RTE2=1.0/TE
        RTNJ=1.0/TN
 13     TDIF=TE-TN
        TDIFRT=TDIF*RTE2*RTNJ
        LEN2=1.77E-19*NE*N2*(1.0-1.2E-4*TE)*TE*TDIF
        LEO2=1.21E-18*NE*O2*(1.0+3.6E-2*SQTE)*SQTE*TDIF
        LEO=7.9E-19*NE*OX*(1.0+5.7E-4*TE)*SQTE*TDIF
        LEH=9.63E-16*NE*HN*(1.0-1.35E-4*TE)*SQTE*TDIF
C
C###  ROTATIONAL loss RATES B&K 268 ####
        LRN2=2.0E-14*NE*N2*TDIF/SQTE
        LRO2=7.0E-14*NE*O2*TDIF/SQTE
C
C+++  N2 VIB loss RATES B&K P268 AND S&N P364  +++++
        EF=1.06E+4+7.51E+3*TANH(1.1E-3*(TE-1800.0))
        GE=3300.0+1.233*(TE-1000.0)-2.056E-4*(TE-1000.0)*(TE-4000.0)
        EXH=5.0E-4*RTE2*EF*(TE-2000.0)
        IF(EXH.GT.70.0)EXH=70.0
        EXH2=GE*TDIFRT
        IF(EXH2.GT.70.0)EXH2=70.0
        LVN2=-2.99E-12*NE*N2*EXP(EXH)
     > *(EXP(-EXH2)-1.0)

C$$$$ O2 VIB loss RATE S&N P364  $$$$
        HS=3300.0-839.0*SIN(1.91E-4*(TE-2700.0))
        EXH=1.4286E-3*RTE2*HS*(TE-700.0)
        IF(EXH.GT.70.0)EXH=70.0
        EXH2=2770.0*TDIFRT
        IF(EXH2.GT.70.0)EXH2=70.0
        LVO2=-5.196E-13*NE*O2*EXP(EXH2)
     >        *(EXP(-EXH2)-1.0)
C
C;;;;;  FINE STRUCTURE EXCITATIONS S&N P365 ;;;;;;
C;;;;; NOTE THAT D1,D2,E1 ETC. MAY NEED TO BE CHANGED
C;;;; NOTE THAT A TERM IS ADDED TO NE AT LOW ALTS FOR O2+ , NO+
        HD1=228.0*RTNJ
        IF(HD1.GT.70.0)HD1=70.0
        D1=EXP(-HD1)
        HD2=326.0*RTNJ
        IF(HD2.GT.70.0)HD2=70.0
        D2=EXP(-HD2)
        HE1=228.0*RTE2
        IF(HE1.GT.70.0)HE1=70.0
        E1=EXP(-HE1)
        HE2=326.0*RTE2
        IF(HE2.GT.70.0)HE2=70.0
        E2=EXP(-HE2)
        HE3=98.0*RTE2
        IF(HE3.GT.70.0)HE3=70.0
        E3=EXP(-HE3)*D1
        LF1=8.49E-6*TE**0.519*(0.02*(D1-E1)-5.91E-9*
     >     TDIF*(2.019*D1+(228.0*RTE2+2.019)*E1))
        LF2=7.7E-6*TE**0.3998*(0.028*(D2-E2)-5.91E-9*
     >     TDIF*(1.8998*D2+(326.0*RTE2+1.8998)*E2))
        LF3=2.22E-7*TE**0.768*(0.008*(D2-E3)-5.91E-9*
     >     TDIF*(2.268*D2+(98.0*RTE2+2.268)*E3))
        ZFO=5.0+3.0*D1+D2
        LFO=-8.629E-6*NE*OX*(LF1+LF2+LF3)/ZFO
C....... LFO=0.333*LFO
C++++  EXCITATION OF O(1D) BY PE'S S&N P365. - SMALL ONLY CALC FOR PRINT
        DE=2.4E4+(TE-1500.0)*(0.3-1.947E-5*(TE-4000.0))
        EXH=3.3333E-4*DE*(TE-3000.0)*RTE2
        IF(EXH.GT.70.0)EXH=70.0
        EXH2=22713.0*TDIFRT
        IF(EXH2.GT.70.0)EXH2=70.0
        LF1D=-1.57E-12*NE*OX*EXP(EXH)*(EXP(-EXH2)-1.0)

C......... add other losses to electron-ion
      TLSS=LEN2+LEO2+LEO+LEH+LRO2+LRN2+LVO2+LVN2+LF1D
	IF(TLSS.LT.0) TLSS=0.0

C......... electron heating from N(2D) .....
      DO III=1,8
        FD(III)=0.0
      ENDDO
!     IF(JPR.GT.0) WRITE(ID,54) ALT,INT(TE),INT(TI),INT(TN),KEI
!    > ,LRN2,LRO2,LVN2,LVO2,LFO,LF1D,TLSS+LFO+KEI,EHT,HN2D
!54    FORMAT(F9.1,3I5,1P,22E10.2)

      COOL(1)=EHT
      COOL(2)=KEI
      COOL(3)=LVN2
      COOL(4)=LFO+LF1D
      COOL(5)=TLSS+LFO+KEI
      RETURN
      END
