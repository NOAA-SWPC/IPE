C........................ INIT-PROFILES.FOR.................................
C.... This file contains routines for setting up initial profiles and 
C.... Adjusting the H+ He+ depleted flux tube profiles
C::::::::::::::::::::::::::::::: PROFIN ::::::::::::::::::::::::::::::::::
C....... set up rough initial O+, H+, and temperature profiles
      SUBROUTINE PROFIN(IHEPLS,INPLS,PCO,F107,N,TI,HPEQ,HEPRAT)
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P) ISPEC
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      IMPLICIT NONE
      INTEGER J,JEQ
      INTEGER IHEPLS,INPLS  !.. switches He+ and N+ diffusive solutions on if > 0
      REAL F107
      DOUBLE PRECISION ALT,N(4,IDIM),TI(3,IDIM),HPEQ,XHOLD,HEPRAT,
     >  PCO,HP_FULL,HP_D_EQ

      !..... Temperatures approximated by a log function
      DO J=JMIN,JMAX
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        TI(3,J)=904.0*DLOG(ALT)-3329.0
        IF(TI(3,J).GT.3000) TI(3,J)=3000
        TI(1,J)= TI(3,J)
        TI(1,J)=0.5*(TN(J)+TI(3,J))        !$$$
        TI(2,J)= TI(1,J)
      ENDDO
      !.. Northern hemisphere O+
      JEQ=(JMAX+1)/2
      DO J=JMIN,JEQ
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        !...... Insert initial O+ profile: first low altitudes
        XHOLD=-5.5E-4*(ALT-250.0)**2
        IF(ALT.LE.500) N(1,J)=(F107/74.0)*5E5*EXP(XHOLD)
        !....... now high altitudes
        IF(ALT.GE.275) THEN 
          XHOLD=(SL(J)-SL(J-1))*1.9E-7*GR(J)/(TI(3,J)+TI(1,J))
          N(1,J)=N(1,J-1)*EXP(XHOLD)
        ENDIF 
      ENDDO

      !.. Southern hemisphere O+ 
      DO J=JMAX,JEQ+1,-1
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        !...... Insert initial O+ profile: first low altitudes
        XHOLD=-5.5E-4*(ALT-250.0)**2
        IF(ALT.LE.500) N(1,J)=(F107/74.0)*5E5*EXP(XHOLD)
        !....... now high altitudes
        IF(ALT.GE.275) THEN
          XHOLD=-(SL(J+1)-SL(J))*1.9E-7*GR(J)/(TI(3,J)+TI(1,J))
          N(1,J)=N(1,J+1)*EXP(XHOLD)
        ENDIF
      ENDDO

      !... set default equatorial H+ density for a full flux tube
      HP_FULL=5.0E4/PCO**3
      HP_D_EQ=DABS(HPEQ)*HP_FULL !.. initial maximum H+ density

      !.. Now calculate H+, He+, and N+
      DO J=JMIN,JMAX
        ALT=Z(J)
        IF(ALT.LT.50) ALT=50
        !.. Chemical equilibrium at low altitudes 
        IF(ON(J).GT.0) N(2,J)=N(1,J)*HN(J)/ON(J)
        IF(N(2,J).GT.HP_D_EQ.OR.Z(J).GE.2000) N(2,J)=HP_D_EQ
        IF(N(2,J).LT.5.0) N(2,J)=5.0

        !-- He+ density is fraction of H+
        XIONN(3,J) = 0.0
        IF(IHEPLS.GT.0) XIONN(3,J)=HEPRAT*N(2,J)
        !-- N+ density is 10% of O+
        XIONN(4,J) = 0.0
        IF(INPLS.GT.0) XIONN(4,J) = 0.1 * N(1,J)
        XIONN(1,J)=N(1,J)   !.. store O+ 
        XIONN(2,J)=N(2,J)   !.. store H+
      ENDDO
      HPEQ=0.0   !.. switches off initial profiles
      RETURN
      END
C:::::::::::::::::: NEW_HP ::::::::::::::::::::::::
C... This routine adjusts the H+ and He+ density for storm depletions.
C... O+ and N+ are unchanged.
C... Written by P. Richards September 2010
      SUBROUTINE NEW_HP(JMIN,     !... First index on field line
     >                  JMAX,     !... Last index on field line
     >                   PCO,     !... L-shell
     >                  HPEQ,     !... New H+ value
     >                     N,     !... ion densities
     >                 EFLAG)     !.. Error flag array
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER J,JMIN,JMAX
      INTEGER EFLAG(11,11)         !.. error flags
      DOUBLE PRECISION ALPHA,HP_MIN,HPEQ,PCO,N(4,IDIM)

      !.. Equatorial density for a depleted flux tube (~20% full).
      HP_MIN=300/PCO**2
      IF(HP_MIN.LT.5.0) HP_MIN=5.0

      ALPHA=HP_MIN/N(2,(JMIN+JMAX)/2)   !.. reduction factor at equator

      !.. Don't delete an already depleted flux tube
      IF(ALPHA.GT.0.8) HPEQ=0.0
      IF(ALPHA.GT.0.8) RETURN

      !.. Reduce H+ and He+ densities. 
      DO J=JMIN,JMAX
        N(2,J)=N(2,J)*ALPHA
        XIONN(3,J)=XIONN(3,J)*ALPHA
      ENDDO
      !.. Debug write if EFLAG(11,11)=1
      IF(EFLAG(11,11).EQ.1) WRITE(6,661) N(2,(JMAX+1)/2),ALPHA,HPEQ
 661  FORMAT(' H+ and He+ reduced; new equatorial [H+]=',F10.2
     >   ,1X,' for large Kp, ALPHA=',F8.1,' HPEQ=',F8.1)
      HPEQ=0.0   !.. switches off reduced densities
      RETURN
      END 
C:::::::::::::::::::::::::::::::::::::: NEW_HPEQ ::::::::::::::::
C... This routine determines if flux tube H+ density needs to be depleted based on Kp.
C... Different formulas are provided. There are 2 additional features
C...   1. You can specify by how much Kp has to change for depletion to occur.
C...      Probably should be between 0.0 and 1.0
C...   2. KPSAVE is used to avoid unnecessary depletions when the flux tube density 
C...      is already very low. 
C... Binsack JGR 1967, page 5321 gives Lpp=6.0-0.6*Kp
C... Moldwin et al. JGR 2002 Lpp=4.7-0.31*Kp (day);=5.5-0.37*Kp (night)
C... Carpenter and Anderson JGR 1991 Lpp=5.6-0.46Kp
C... Richards KPPP is based on FLIP comparisons with Chappell 1972 profiles
C... Horwitz et al. JGR 1990, page 7939 gives similar values at different MLTs
C      KPPP=(5.5-PCO)/0.37    ! Kp for plasmapause: Moldwin Night
C      KPPP=(4.7-PCO)/0.31    ! Kp for plasmapause: Moldwin day
C      KPPP=(5.6-PCO)/0.46    ! Kp for plasmapause: Carpenter
C      KPPP=(6.0-PCO)/0.60    ! Kp for plasmapause: Binsack
C      KPPP=(5.6-PCO)/0.41    ! Kp for plasmapause: Horwitz 0-6 MLT
C      KPPP=(6.2-PCO)/0.45    ! Kp for plasmapause: Horwitz 0-12 MLT
C      KPPP=(10.3-PCO)/1.0    ! Kp for plasmapause: Horwitz 12-18 MLT
C      KPPP=(7.4-PCO)/0.65    ! Kp for plasmapause: Horwitz 18-24 MLT
C      KPPP=(6.0-PCO)/0.50    ! Kp for plasmapause: Richards
      SUBROUTINE NEW_HPEQ(KP,  !..  IN: Current Kp value
     >                   PCO,  !..  IN: L-shell
     >             DEN_HP_EQ,  !..  IN: Equatorial H+ density
     >                KPSAVE,  !..  OUT: Saved Kp value
     >                  HPEQ)  !..  OUT: =0/1 Signal if H+ depletion
      REAL DEN_HP_EQ,KP,KPSAVE
      DOUBLE PRECISION PCO,HPEQ
      REAL KPPP    !.. Kp of plasmapause location
      REAL KPINC   !.. How much Kp has to increase before depletion occurs
      DATA KPINC/0.7/

      !.. Make sure flux tube is not over full. Max tube content ~ 5E14 ions for L>3.0
      !.. Flux tubes < 3 will never reach this limit. Tube volume goes as ~1.0/L**4 
      IF(DEN_HP_EQ.GT.6E5/PCO**4) HPEQ=-1   !.. Not needed if Kp depletions

      !.. Plasmapause location. If desired, replace with different KPPP above
      KPPP=(6.0-PCO)/0.5        !.. Kp for plasmapause Richards

      !.. Determine if flux tube needs to be deleted
      IF(KP-KPSAVE.GT.KPINC.AND.KP.GE.KPPP) HPEQ=-1    ! PGR 

      KPSAVE=KP  !.. Save Kp for changing Lpp only when Kp increases
      RETURN
      END
