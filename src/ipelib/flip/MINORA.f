C.................... MINORA.FOR .......... .................
C.... This function is the main sequencing control for the He+ and N+
C.... routines. It is similar to subroutine LOOPS on RSLPSD.FOR. Look 
C.... there for comments.
C.... Cleaned up and commented by P. Richards in April 2000
      SUBROUTINE XION(TI,    !.. O+,H+ & Ti,Te
     >              DTIN,    !.. Time step in
     >             DTMIN,    !.. Minimum time step
     >            IHEPNP,    !.. He+ - N+ switch
     >             EFLAG,mp,lp,i_which_call)    !.. OUTPUT: Error flag array
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE SOLVARR        !... DELTA RHS WORK S, Variables for solver
      USE ION_DEN_VEL    !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      IMPLICIT NONE
      integer mp,lp,i_which_call
      INTEGER NION,IHEPNP,J,ITER,IRHS,IBW,I
      INTEGER EFLAG(11,11),NFLAG                  !.. solution procedure error flags
      INTEGER JBNN,JBNS,JEQ,JC                    !.. boundary indices
      INTEGER IDIV,KR,ION,IEQ,MIT                 !.. solution variables
      DOUBLE PRECISION N(4,FLDIM),NMSAVE(2,FLDIM) !.. Solution and saved densities 
      DOUBLE PRECISION TI(3,FLDIM)                !.. Ti,Te
      !.. solution variables
      DOUBLE PRECISION DCR,DCRP,DCRQ(2,FLDIM),ADCR,DINC,F(20),DUM(9)
      DOUBLE PRECISION DT,DTIN,DTMIN,DTINC !.. Time step variables
      DOUBLE PRECISION ZLBHE,ZLBNP         !.. He+ & N+ lower boundary
      DOUBLE PRECISION XMAS,XMASS(9)       !.. Ion mass
      DOUBLE PRECISION NMORIG(2,FLDIM)     !.. Original density at previous time step
      integer ret
      DATA XMASS/23.4164D-24,26.7616D-24,6.6904D-24,6*0.0/
      DATA NION/1/      ! # species: do not change

      IF(IABS(IHEPNP).EQ.9)  XMAS=XMASS(3)  !.. FOR He+
      IF(IABS(IHEPNP).EQ.11) XMAS=XMASS(1)  !.. FOR N+

      ZLBHE=200         !.. set He+ lower boundary
      ZLBNP=115         !.. set N+ lower boundary
      DT=DTIN           !.. Set time step for dN/dt
      DTINC=0.0         !.. Used for reduced timestep
      JBNN=JMIN         !.. lower boundary point for O+ and H+
      JBNS=JMAX         !.. lower boundary point for O+ and H+
      JEQ=(JMIN+JMAX)/2 !.. Equatorial point

      !..transferring minor ion densities from storage in XION to the 
      !.. solution variable N
      DO J=JMIN,JMAX
         !... for N+
         IF(IABS(IHEPNP).EQ.11) THEN
           N(1,J)=XIONN(4,J)
           N(3,J)=XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J)
         ENDIF
         !... for He+
         IF(IABS(IHEPNP).EQ.9) THEN
           N(1,J)=XIONN(3,J)
           N(3,J)=XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J)
         ENDIF
         !..store minor ion density to calculate dn/dt
         NMSAVE(1,J)=N(1,J)    !.. density at previous intermediate time step
         NMORIG(1,J)=N(1,J)    !.. Original density at previous time step
         DCRQ(1,J)=1.0
         DCRQ(2,J)=1.0
      ENDDO

!dbg20120301: this part was commented out on 20110815, but un-comment again to solve for N+ problem: "IN BDSLV &&&&&&& BANDWIDTH IS TOO LARGE " i don't remember why we decided to comment out here and i cannot find a program to setup the local chem equil anywhere else???
!dbg20110815      !.. Use local equilibrium for densities if flux tube apex height < 200 km
!      IF(Z(JEQ).LT.200.0) THEN 
!        DO J=JMIN,JMAX
!          IF(IABS(IHEPNP).EQ.9) CALL MCHEMQ(J,DUM,N,TI,FLDIM)
!          IF(IABS(IHEPNP).EQ.11) CALL NCHEMQ(J,DUM,N,TI,FLDIM)
!          N(1,J)=DUM(1)
!        ENDDO
!        RETURN
!      ENDIF
!dbg20110815
      !.. Use local equilibrium for He+ densities if apex height < ZLBHE
      IF(IHEPNP.EQ.9.AND.Z(JEQ).LE.ZLBHE) THEN
        DO J=JMIN,JMAX
          CALL MCHEMQ(J,DUM,N,TI,FLDIM)
          XIONN(3,J)=DUM(1)      !.. He+
        ENDDO
        RETURN
      ENDIF

      !.. Use local equilibrium for N+ densities if apex height < ZLBNP
      IF(IHEPNP.EQ.11.AND.Z(JEQ).LE.ZLBNP) THEN
        DO J=JMIN,JMAX
          CALL NCHEMQ(J,DUM,N,TI,FLDIM)
          XIONN(4,J)=DUM(1)      !.. N+
        ENDDO
        RETURN
      ENDIF 
!dbg20110815


      !... calculate average values for quantities that don't
      !... change in Newton iteration
      CALL DENAVE(JMAX-1,TI)
      CALL AVDEN2(JMAX-1,TI,IHEPNP)

C- OUTER LOOP Return here on Non-Convergence with reduced time step
  10  CONTINUE

        !.. Update indices for the lower boundaries North. 
        DO J=JMIN,JEQ-1
          IF(IABS(IHEPNP).EQ.9.AND.Z(J).LE.ZLBHE) JBNN=J
          IF(IABS(IHEPNP).EQ.11.AND.Z(J).LE.ZLBNP) JBNN=J
        ENDDO
        !.. Update indices for the lower boundaries South. 
        DO J=JMAX,JEQ,-1
          IF(IABS(IHEPNP).EQ.9.AND.Z(J).LE.ZLBHE) JBNS=J
          IF(IABS(IHEPNP).EQ.11.AND.Z(J).LE.ZLBNP) JBNS=J
        ENDDO

        !.. Set up boundary indices for the solution procedure
        MIT=JBNS-JBNN+1       !.. Number of points on field line
        IEQ=MIT-2             !.. Number of equations to set up      

        !************* Main Newton Solver Iteration Loop begins *****************
        DO 220 ITER=1,20
          !.. set boundary conditions on density
          DO J=JMIN,JMAX
            IF(J.LE.JBNN.OR.J.GE.JBNS) THEN
              IF(IABS(IHEPNP).EQ.9) CALL 
     >          MCHEMQ(J,DUM,N,TI,FLDIM)
              IF(IABS(IHEPNP).EQ.11) 
     >         CALL NCHEMQ(J,DUM,N,TI,FLDIM)
              N(1,J)=DUM(1)
            ENDIF
          ENDDO

          !.. call MDFIJ to get unperturbed value to calculate dFij/dn
          DO J=2,MIT-1
            KR=NION*(J-2)
            JC=J+JBNN-1
            CALL MDFIJ(JC,0,NION,DT,N,TI,F,JBNN,JBNS,IHEPNP,XMAS,NMSAVE)
            DO IRHS=1,NION
              RHS(KR+IRHS)=F(IRHS)
            ENDDO
          ENDDO
          !.. Now set up the Jacobian Matrix dFij/dn
!nm20130111: distinguished he+ v.s. n+ according to phil's suggestion
          IF(IABS(IHEPNP).EQ.9) !he+  
     >     CALL HMATRX(FLDIM,S,RHS,IEQ,DT,N,TI,JBNN,JBNS,MIT,IHEPNP,
     >      XMAS,NION,NMSAVE)
          IF(IABS(IHEPNP).EQ.11)  !n+ 
     >     CALL NMATRX(FLDIM,S,RHS,IEQ,DT,N,TI,JBNN,JBNS,MIT,IHEPNP,
     >      XMAS,NION,NMSAVE)


          !.. invert the jacobian matriX *s* in the inversion routine *bdslv*.
          !.. the increments are stored in array delta in this order
          !.. X(1...n,j),X(1...n,j+1),X(1...n,j+2),....X(1...n,jmaX-1)
          IBW=2*NION-1
          CALL BDSLV(IEQ,IBW,S,0,RHS,DELTA,WORK,NFLAG,
     >                mp,lp,i_which_call)

          !.. Check for problems in band solver
          EFLAG(3,2)=0     
          EFLAG(4,2)=0     
          IF(NFLAG.NE.0) THEN
            IF(EFLAG(11,11).EQ.1) WRITE(6,*) ' '
            IF(EFLAG(11,11).EQ.1) WRITE(6,*) 
     >        ' *** Problem in band solver ****'
            IF(IABS(IHEPNP).EQ.9) EFLAG(3,2)=-1    !.. Report problem to calling routine
            IF(IABS(IHEPNP).EQ.11) EFLAG(4,2)=-1   !.. Report problem to calling routine
            RETURN
          ENDIF

          IDIV=0             !... convergence indicator
          DCR=1.0            !... convergence indicator

          !.. test density increments 'dinc' to ensure density > 0 
          !.. (mod steepest descent)
          DO  142 I=1,NION
            ION=2*NION-I
            DO 142 J=2,MIT-1
              DCRP=1.0
              JC=JBNN+J-1
              DCRQ(I,JC)=1.0
              DINC=DELTA(NION*J-ION)
              IF(DINC.LE.0) GO TO 142
              IF(ABS(DINC/N(I,JC)).GT.0.9999) DCRP=0.5*ABS(N(I,JC)/DINC)
              IF((DCRP.LT.DCR).AND.(Z(JC).GT.0.0)) DCR=DCRP
              IF(ITER.GT.0) DCRQ(I,JC)=DCRP
 142      CONTINUE
 
          ADCR=DCR
          !.. add iterative increment to the density. 
          !.. convergence when IDIV=0. 
          DO 42 I=1,NION
            ION=2*NION-I
            DO 42 J=2,MIT-1
              JC=JBNN+J-1
              DINC=DELTA(NION*J-ION)
              N(I,JC)=N(I,JC)-DINC*ADCR
              IF(ABS(DINC/N(I,JC)).GT.1E-3)  IDIV=IDIV+1
 42       CONTINUE

          !.. test to see if convergence has occured.
          IF(IDIV.EQ.0)  GO TO 230
 220    CONTINUE

 230    CONTINUE    !*************** end main loop ***************

        !.. Testing for convergence to see if need to reduce time step
        IF(IDIV.EQ.0) THEN
          !============== Convergence success ========
          EFLAG(3,1)=0     
          EFLAG(4,1)=0     
          DTINC=DTINC+DT   !.. Used for reduced timestep
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)')  
     >       ' He+ N+ ',ITER,DTINC,DTIN,DT
          IF(DTINC.GE.DTIN-1) THEN
            DO J=JMIN,JMAX
              IF(IABS(IHEPNP).EQ. 9) XIONN(3,J)=N(1,J)  ! He+
              IF(IABS(IHEPNP).EQ.11) XIONN(4,J)=N(1,J)  ! N+
            ENDDO
            RETURN
          ENDIF
          !.. increase time step if convergence is easy
          IF(ITER.LT.5.AND.DTINC+2*DT.LE.DTIN) DT=2*DT
          IF(ITER.LT.5.AND.DTINC+2*DT.GT.DTIN) DT=DTIN-DTINC

          !-- Save current densities for dN/dt. 
          DO J=JMIN,JMAX
            NMSAVE(1,J)=N(1,J)
          ENDDO
        ELSE
        !============== Convergence failure ========
        !-- Non-Convergence: Reduce time step and restore densities. 
          DT=DT/2  !.. reduce time step for non-convergence
          !.. Raise lower boundary to maximum of 300 km
          IF(IHEPNP.EQ.9.AND.DT.LT.DTIN/3.0)  ZLBHE=(ZLBHE+350)/2   
          IF(IHEPNP.EQ.11.AND.DT.LT.DTIN/3.0)  ZLBNP=(ZLBNP+350)/2   
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,2I5,9F14.2)')  
     >       ' He+ N+ 2nd ',-IHEPNP,ITER,DTINC,DTIN,DT,ZLBHE,ZLBNP
          DO J=JMIN,JMAX
            N(1,J)=NMSAVE(1,J)
          ENDDO

          !.. Check that DT is not too small
          IF(DT.LT.DTMIN) THEN
            IF(IABS(IHEPNP).EQ.9) EFLAG(3,1)=-1    !.. Report problem to calling routine
            IF(IABS(IHEPNP).EQ.11) EFLAG(4,1)=-1   !.. Report problem to calling routine
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,9I5)') 
     >        '  ERR FLAGS MINA',IHEPNP
            !.. Restore density to original input value
            DO J=JMIN,JMAX
              IF(IABS(IHEPNP).EQ.9) XIONN(3,J)=NMORIG(1,J)
              IF(IABS(IHEPNP).EQ.11) XIONN(4,J)=NMORIG(1,J)
            ENDDO
            RETURN
          ENDIF 
        ENDIF

      GOTO 10
C- END OF OUTER LOOP ----------------------------
      RETURN
      END
C::::::::::::::::::::::::::::: MDFIJ ::::::::::::::::::::::::::::::::::::
C.... This routine sets up the He+ and H+ continuity equations to be 
C.... minimized by the Newton solver. See program DFIJ on RSDENA.FOR 
C.... for more details

      SUBROUTINE MDFIJ(J,JSJ,NION,DT,N,TI,F,JBNN,JBNS,IHEPNP,XMAS,
     >  NMSAVE)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL    !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER I,J,JSJ,NION,IHEPNP,ID
      INTEGER JBNN,JBNS,JEQ                !.. boundary indices
      DOUBLE PRECISION VL(2),VU(2),FLU(2),FLL(2),ANM(2),PM(2),Q(2)
      DOUBLE PRECISION TINCR(2),N(4,FLDIM),TI(3,FLDIM),F(20)
      DOUBLE PRECISION NMSAVE(2,FLDIM),L(2),XMAS
      DOUBLE PRECISION BU,BL,BD,DT,FGR
      DATA FLU/0.0,0.0/

      !... ID is the unit number for writing diagnostics in
      !... both hemispheres
      ID=3  
      IF(J.GT.JMAX/2) ID=9
      JEQ=(JMAX+1)/2
      IF(JSJ.NE.9) ID=0

      !.. MINVEL to get the velocities(fluxes) at the lower 1/2 pt
      FLL(1)=FLU(1)
      FLL(2)=FLU(2)
      CALL MINVEL(J-1,VL,FLL,N,JSJ,0,XMAS,IHEPNP)

      !..  CALL CHEMO to evaluate the source (Q) and sink(L) terms 
      IF(IABS(IHEPNP).EQ.9) 
     >   CALL CHEMO9(NION,J,Q,N,NMSAVE,TI,L)
      IF(IABS(IHEPNP).EQ.11) CALL CHEM11(NION,J,Q,N,NMSAVE,TI,1,L)

      !..  Call DAVE for dn/dt. PM(AM)=future(present) cpt of dn/dt
      CALL DAVE(NION,J,ANM,PM,N,NMSAVE)

      !.. MINVEL to get the velocities(fluxes) at the upper 1/2 pt
360   CALL MINVEL(J,VU,FLU,N,JSJ,ID,XMAS,IHEPNP)

      !.. magnetic field strength at midpoints
      BU=(BM(J)+BM(J+1))*0.5
      BL=(BM(J)+BM(J-1))*0.5

      !.. Set up F=prod-loss-dn/dt-div(flux)
      DO 380 I=1,NION
      TINCR(I)=(PM(I)-ANM(I))/DT     !.. dn/dt
      FGR=(FLU(I)/BU-FLL(I)/BL)      !.. div(flux)

      F(I)=Q(I)-L(I)-TINCR(I)-FGR

      !.... average velocity for He+ and N+ ........
      IF(JSJ.EQ.0.AND.IABS(IHEPNP).EQ.9)  XIONV(3,J)=.5*(VL(I)+VU(I))
      IF(JSJ.EQ.0.AND.IABS(IHEPNP).EQ.11) XIONV(4,J)=.5*(VL(I)+VU(I))
      VL(I)=VU(I)

      !.. printing of factors in the continuity equations ...........
      !.. since all quantities have been divided by bm and integrated,
      !.. multiply by bd to get actual magnitudes
      IF(ID.EQ.0) GO TO 380
        BD=2.0*BM(J)/(SL(J+1)-SL(J-1))
        WRITE(ID,666) J,BD,FLU(I),FLL(I),FGR,Q(I),L(I),F(I),XIONV(3,J)
     >   ,TINCR(I),XIONV(4,J)
  666    FORMAT(2X,'FIJ2',I4,1P,22E11.3)
         IF(I.EQ.2) WRITE(ID,*) '   '
  380 CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::: MCHEMQ :::::::::::::::::::::::::::::::::::::
C... Gets the chemical equilibrium density for He+
      SUBROUTINE MCHEMQ(J,DUM,N,TI,FLDIM)
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      IMPLICIT NONE
      INTEGER J,FLDIM     !.. Field line index & dimension
      DOUBLE PRECISION N(4,FLDIM),TI(3,FLDIM)    !.. He+,N+ & Ti,Te
      DOUBLE PRECISION RTS(99),DUM(9),HELOSS
      CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
        HELOSS=(RTS(44)+RTS(45))*N2N(J)+(RTS(75)+RTS(76))*O2N(J)
       !.. small He+ production added to avoid convergence problems
        DUM(1)=(OTHPR1(2,J)+1.0E-07)/HELOSS 
        DUM(2)=0.0
      RETURN
      END
C::::::::::::::::::::::::::::: NCHEMQ :::::::::::::::::::::::::::::::::::::
C... Gets the chemical equilibrium density for N+
      SUBROUTINE NCHEMQ(J,DUM,N,TI,FLDIM)
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE MINORNEUT !.. N4S N2D NNO N2P N2A O1D O1S
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      IMPLICIT NONE
      INTEGER J,FLDIM     !.. Field line index & dimension
      DOUBLE PRECISION N(4,FLDIM),TI(3,FLDIM)    !.. He+,N+ & Ti,Te
      DOUBLE PRECISION RTS(99),DUM(9),QNP,PNP
        CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
        !... Evaluate loss rates
        QNP=(RTS(30)+RTS(25)+RTS(59)+RTS(22)+RTS(65)+RTS(66))*O2N(J)
     >     +RTS(31)*ON(J)
        !... evaluate production rates
        PNP=RTS(45)*XIONN(3,J)*N2N(J)+NPLSPRD(J)+RTS(29)*N2D(J)*
     >    XIONN(1,J)
        DUM(1)=PNP/QNP        !.. chemical equilibrium density for N+
        DUM(2)=0.0
      RETURN
      END
C::::::::::::::::::::::::::::: MINVEL ::::::::::::::::::::::::::::::::::
C... This subroutine returns the ion fluxes for the continuity equation.
      SUBROUTINE MINVEL(J,V,FLUX,N,JSJ,ID,XMAS,IHEPNP)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc. 
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER J,IHEPNP,JSJ,ID
      DOUBLE PRECISION PP(4),QSIGN(2),V(2),FLUX(2),GAMMA(2),B(3)
      DOUBLE PRECISION N(4,FLDIM),HX(2),GRA
      DOUBLE PRECISION NE,NEU,NEL,D,SUMN,XMAS,GRADNE,DLNI,QION,Q2,Q3,Q4
      DATA QSIGN/-1.,1./

      !... evaluate total ion density
      NEU=XIONN(1,J+1)+XIONN(2,J+1)
      NEL=XIONN(1,J)+XIONN(2,J)
      NE=DSQRT(XIONN(1,J+1)*XIONN(1,J))+DSQRT(XIONN(2,J)*XIONN(2,J+1))
      IF(IABS(IHEPNP).GT.9) THEN
         NEU=NEU+XIONN(3,J+1)+N(1,J+1)+N(3,J+1)
         NEL=NEL+XIONN(3,J)+N(1,J)+N(3,J)
         NE=NE+DSQRT(XIONN(3,J)*XIONN(3,J+1))+DSQRT(N(1,J+1)*N(1,J))
         NE=NE+DSQRT(N(3,J)*N(3,J+1))
      ENDIF
      !.. interpolate minor ion to midpoint
      PP(1)=DSQRT(N(1,J+1)*N(1,J))
      PP(2)=0.0
      !.. CALCULATE ELECTRON DENSITY AT UPPER, LOWER AND MID POINT
      IF (IABS(IHEPNP).LE.9) THEN
        NEU=NEU+N(1,J+1)+N(3,J+1)
        NEL=NEL+N(1,J)+N(3,J)
        NE=NE+PP(1)+DSQRT(N(3,J)*N(3,J+1))
      ENDIF

      !.. CALL DCOEFF to obtain d(i),hX,beta1,beta2,betaX, @ b(i) used
      !.. to calc v(i)  all from st. maurice and schunk
      CALL DCOEFF(J,PP,D,B,HX,GAMMA,JSJ,ID,SUMN,XMAS)    ! f90 added

      !-----calculate altitude derivatives----------
      !.. note that TEJ(J)=0.5(TI(3,J+1)+TI(3,J))/DS(J)/TIJ(J) see AVDEN
      GRADNE=TEJ(J)*(NEU-NEL)/NE

      DLNI=(N(1,J+1)-N(1,J))/DS(J)/PP(1)              !.. dn/ds
      GRA=XMAS*GRAV(J)                                !.. gravity
      QION=QSIGN(2)*GRADTE(J)

      !.... neutral wind term SUMN=sum of collision terms for ions on neutrals
      Q2=UNJ(J)*SUMN
      !.. Q3 is the sum of FRICTION terms.  HX(1)=X - H+, HX(2)=X - O+
      Q3=0.5*(XIONV(1,J+1)+XIONV(1,J))*HX(2)
      Q3=Q2+Q3+0.5*(XIONV(2,J+1)+XIONV(2,J))*HX(1)
      Q4=GRADTI(J)*(1.0D00-(B(1)+B(2)+B(3)))

      !.. Velocity from the momentum equation.
       V(1)=Q3-D*(DLNI-GRA+Q4+GRADNE+QION)

      !.......  print routine for momtm eqn terms ..........
      IF(ID.NE.0) WRITE(ID,605) J,Z(J),DLNI,GRADTI(J),GRADNE
     > ,GRADTE(J),GRA,QION,Q4,Q2,Q3
 605   FORMAT(1X,'MINVEL',I4,F9.0,1P,22E11.3)

       !.. ion flux
       FLUX(1)=V(1)*PP(1)
       !... print velocities ...........
      IF(JSJ.EQ.-4) WRITE(3,606) J,Z(J),V(1),V(2)
 606  FORMAT(1X,'J=',I4,3X,'ALT=',F9.0,2X,'V1=',1P,E13.6,2X,'V2=',E13.6)
      RETURN
      END
C::::::::::::::::::::::::::::: DCOEFF ::::::::::::::::::::::::::::::::::::
C... This routine calculates the diffusion velocities from the minor ion
C... equation 34 as given by St. Maurice and Schunk Planet and Spac. 
C... Sci. 1976.  Units are cgs
      SUBROUTINE DCOEFF(J,PP,D,B,HX,GAMMA,JSJ,ID,SUMN,XMAS)
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc.
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT NONE
      INTEGER I,J,K,JSJ,MM,ID
      DOUBLE PRECISION NU(3,3),NUPI(3,2),N(3),NUP(3),MU(3,3),KB
      DOUBLE PRECISION D1(3,3),D4(3,3),A(3),AM(3),PP(4),AMASS,RE,XMAS
      DOUBLE PRECISION B(3),DELT(2),HX(2),R(2,2),Y(3,2),GAMMA(2)
      DOUBLE PRECISION AI2,AK2,AIK2,SUMN,AST,EPS,X1,DTOP,NUSUM,D,BTOPC
      DATA A/1.0D00,16.0D00,0.0D00/,KB/1.38062D-16/
      DATA AMASS/1.67D-24/,RE/6.371D08/

      !.. set up constants: A= amu of mass A1=hydrogen,A2=oxygen,
      !.. C A3=helium or N+. amass=1.67e-24 grams
      AM(1)=A(1)*AMASS
      AM(2)=A(2)*AMASS
      AM(3)=XMAS       !.. either He+ or N+
      A(3)=XMAS/AMASS  !.. either He+ or N+
      !.. load densities H+,O+,and minor ions.
      N(1)=DSQRT(XIONN(2,J+1)*XIONN(2,J))
      N(2)=DSQRT(XIONN(1,J+1)*XIONN(1,J))
      N(3)=PP(1)
      !..  Mass dependent constants
      DO I=1,3
        AI2=A(I)*A(I)
          DO K=1,3
            IF (K.EQ.3.OR.K.EQ.I) GO TO 11
            AK2=A(K)*A(K)
            AIK2=(A(I)+A(K))*(A(I)+A(K))
            D1(I,K)=(3.0*AI2+.1*A(I)*A(K)-.2*AK2)/AIK2
            D4(I,K)=(1.2*AK2-1.5*A(I)*A(K))/AIK2
11          MU(I,K)=A(I)*A(K)*AMASS/(A(I)+A(K))
          END DO
      END DO

      !.. Evaluate collision terms. Neutral terms from AVDEN2
      SUMN=NUX(1,J)   !.. sum of neutral collision terms

      DO I=1,3
        DO K=1,2
          AST=A(I)*A(K)/(A(I)+A(K))
          NU(I,K)=1.27D00*DSQRT(AST)*N(K)/(A(I)*TIJ(J)**1.5)
        END DO
      END DO

      NUP(1)=NU(1,1)+1.25*NU(1,2)*(D1(1,2)+1.5*MU(1,2)/AM(1))
      NUP(2)=NU(2,2)+1.25*NU(2,1)*(D1(2,1)+1.5*MU(2,1)/AM(2))
      NUP(3)=1.25*NU(3,1)*(D1(3,1)+1.5*MU(3,1)/AM(3))
      NUP(3)=NUP(3)+1.25*NU(3,2)*(D1(3,2)+1.5*MU(3,2)/AM(3))
      NUPI(1,2)=1.25*NU(1,2)*(D4(1,2)+1.5*MU(1,2)/AM(1))
      NUPI(2,1)=1.25*NU(2,1)*(D4(2,1)+1.5*MU(2,1)/AM(2))
      NUPI(3,1)=1.25*NU(3,1)*(D4(3,1)+1.5*MU(3,1)/AM(3))
      NUPI(3,2)=1.25*NU(3,2)*(D4(3,2)+1.5*MU(3,2)/AM(3))
      EPS=1.0/(1.0-(NUPI(1,2)/NUP(1))*(NUPI(2,1)/NUP(2)))

      DO I=1,2
        DO K=1,2
          IF (K.EQ.I) THEN
          MM=3-K
          Y(I,K)=1.0-NUPI(3,I)/NUP(3)-NUPI(MM,I)/NUP(MM)*
     >          NUPI(3,MM)/NUP(3)
          ELSE
          Y(I,K)=NUPI(K,I)/NUP(K)-NUPI(3,I)/NUP(3)-
     >          NUPI(K,I)/NUP(K)*NUPI(3,K)/NUP(3)
          ENDIF
        END DO
      END DO

      Y(3,1)=-1.0/EPS
      Y(3,2)=Y(3,1)
      X1=MU(3,1)/AM(3)*NU(3,1)/NUP(3)+MU(3,2)/AM(3)*NU(3,2)/NUP(3)
      DELT(1)=9.0*(MU(3,1)/AM(3))*X1/8.0
      DELT(2)=9.0*(MU(3,2)/AM(3))*X1/8.0
      DTOP=KB*TIJ(J)/AM(3)
      NUSUM=NU(3,1)*(1.0-DELT(1))+NU(3,2)*(1.0-DELT(2))+SUMN

      !.. divide diffusion and neutral collision terms by total of
      !.. of collision terms nusum. redefine neutral collision term sumn
      D=DTOP/NUSUM
      SUMN=SUMN/NUSUM
      BTOPC=15.0*AM(3)*EPS/8.0
      B(1)=BTOPC/AM(1)*(MU(3,1)/AM(1)*NU(3,1)/NUP(1)*Y(1,1)
     >       +MU(3,2)/AM(1)*NU(3,2)/NUP(1)*Y(1,2))
      B(2)=BTOPC/AM(2)*(MU(3,1)/AM(2)*NU(3,1)/NUP(2)*Y(2,1)
     >       +MU(3,2)/AM(2)*NU(3,2)/NUP(2)*Y(2,2))
      B(3)=BTOPC/AM(3)*(MU(3,1)/AM(3)*NU(3,1)/NUP(3)*Y(3,1)
     >       +MU(3,2)/AM(3)*NU(3,2)/NUP(3)*Y(3,2))

      !...  UNEUT=SUMN*VNEUT*1.0D02/NUSUM
      DO I=1,2
        DO K=1,2
          IF (K.NE.I) THEN
            R(I,K)=MU(I,K)/AM(I)*NU(I,K)/AM(I)*EPS/NUP(I)
            X1=MU(3,I)*NU(3,I)*Y(I,I)+MU(3,K)*NU(3,K)*Y(I,K)
            R(I,K)=9.0*R(I,K)*X1/8.0
          ENDIF
        ENDDO
      END DO

      !..HX=velocity coefficients (friction)
        HX(1)=(NU(3,1)*(1.0-DELT(1))-R(1,2)+R(2,1))/NUSUM
        HX(2)=(NU(3,2)*(1.0-DELT(2))-R(2,1)+R(1,2))/NUSUM
      RETURN
      END
C:::::::::::::::::::::::::::::::: CHEMO9 :::::::::::::::::::::::::::::::
C.... this program determines the interpolated production and loss
C.... processes for He+. It CALLs rates to get the rate constants and TERLIN
C.... to do the interpolation. 
      SUBROUTINE CHEMO9(NION,JI,SOURCE,N,NMSAVE,TI,SINK)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production, PHION
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc. 
      IMPLICIT NONE
      INTEGER I,J,K,JI,NION
      DOUBLE PRECISION N(4,FLDIM),TI(3,FLDIM),NMSAVE(2,FLDIM),Q(2,3)
     > ,L(2,3),SOURCE(2),SINK(2),RTS(99)
      DOUBLE PRECISION HELOSS,QM,LM
      !-- Q and L are chemical source and sink
      DO 300 K=1,3
        J=K+JI-2
        CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
        HELOSS=(RTS(44)+RTS(45))*N2N(J)+(RTS(75)+RTS(76))*O2N(J)
        !.. small He+ production added to avoid convergence problems
        Q(1,K)=(OTHPR1(2,J)+1.0E-07)/BM(J)
        L(1,K)=HELOSS*N(1,J)/BM(J)
 300  CONTINUE

      !.... terd is called to interpolate
      DO 400 I=1,NION
        CALL TERLIN(Q(I,1),Q(I,2),Q(I,3),DS(JI-1),DS(JI),QM)
        CALL TERLIN(L(I,1),L(I,2),L(I,3),DS(JI-1),DS(JI),LM)
        SOURCE(I)=QM
        SINK(I)=LM
 400  CONTINUE
      RETURN
      END
C:::::::::::::::::::::::::: CHEM11 ::::::::::::::::::::::::::::::::::::::
C.... this program determines the interpolated production and loss
C.... processes for N+. It CALLs rates to get the rate constants and TERLIN
C.... to do the interpolation.
      SUBROUTINE CHEM11(NION,JI,SOURCE,N,NMSAVE,TI,ISW,SINK)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE MINORNEUT !.. N4S N2D NNO N2P N2A O1D O1S
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc. 
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      INTEGER J,K,JI,NION
      DIMENSION N(4,FLDIM),TI(3,FLDIM),NMSAVE(2,FLDIM),Q(2,3),L(2,3)
     > ,SOURCE(2),SINK(2),RTS(99)

      DO 300 K=1,3
        J=K+JI-2
        CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
        PNP1=RTS(45)*XIONN(3,J)*N2N(J)
        PNP2=NPLSPRD(J)
        PNP3=RTS(29)*N2D(J)*XIONN(1,J)
        LNP1=RTS(30)*O2N(J)
        LNP2=RTS(25)*O2N(J)
        LNP3=RTS(59)*O2N(J)
        LNP4=RTS(31)*ON(J)
        LNP5=RTS(22)*O2N(J)
        LNP6=RTS(65)*O2N(J)
        LNP7=RTS(66)*O2N(J)

        Q(1,K)=(PNP1+PNP2+PNP3)/BM(J)
        Q(2,K)=0.0
        L(1,K)=(LNP1+LNP2+LNP3+LNP4+LNP5+LNP6+LNP7)*N(1,J)/BM(J)
        L(2,K)=0.0
 300  CONTINUE
      !... terd is called to interpolate
      DO 400 I=1,NION
        CALL TERLIN(Q(I,1),Q(I,2),Q(I,3),DS(JI-1),DS(JI),QM)
        CALL TERLIN(L(I,1),L(I,2),L(I,3),DS(JI-1),DS(JI),LM)
        SOURCE(I)=QM
        SINK(I)=LM
 400  CONTINUE
      RETURN
      END
C::::::::::::::::::::::::::::: AVDEN2 ::::::::::::::::::::::::::::::::::::
C.....................................................................
C   this program evaluates the interpolated densities at the midpoints
C   it also evaluates the ion-neutral collision frequencies nuX(i,j)
C   and o+ and h+ production and loss rates
C........................................................................
      SUBROUTINE AVDEN2(JMAX1,TI,IHEPNP)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc. 
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      REAL OA,HA,N2A,O2A,TNJ,TR,HEA
      DIMENSION TI(3,FLDIM)

      !.. ion-neutral coll freqs taken from schunk and nagy rev. geophys. v18
      !.. p813, 1980. first,  resonant charge eXchange from table 5 p823
      DO 10 J=JMIN,JMAX-1
        OA=SQRT(ON(J)*ON(J+1))
        HA=SQRT(HN(J)*HN(J+1))
        N2A=SQRT(N2N(J)*N2N(J+1))
        O2A=SQRT(O2N(J)*O2N(J+1))
        HEA=SQRT(HE(J)*HE(J+1))
        TNJ=0.5*(TN(J)+TN(J+1))
        TR=(TIJ(J)+TNJ)*0.5
        IF(IABS(IHEPNP).EQ.9) CALL ADEN9(J,OA,HA,N2A,O2A,HEA,TIJJ
     >  ,TR,TNJ,NUX)
        IF(IABS(IHEPNP).EQ.11) CALL ADEN11(J,OA,HA,N2A,O2A,HEA
     >  ,TIJJ,TR,TNJ,NUX)
10    CONTINUE
      RETURN
      END
C:::::::::::::::::: ADEN9::::::::::::::::::::::::::::::
C.........For HE+
      SUBROUTINE ADEN9(J,OA,HA,N2A,O2A,HEA,TIJJ,TR,TNJ,NUX)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      REAL OA,HA,N2A,O2A,TNJ,TR,SQT,HEA
      DIMENSION NUX(2,FLDIM)
      SQT=SQRT(TR)
      RHEPHE=8.73E-11*HEA*SQT*(1-0.093*ALOG10(TR))**2

      !.. non-resonant ion neutral interactions table 6
      CHEPN=(4.71*HA+10.1*OA+16*N2A+15.3*O2A)*1.0E-10
      NUX(1,J)=RHEPHE+CHEPN
      NUX(2,J)=0.0

      !... print diagnostics
      !  ID=3
      !  IF(J.GT.JMAX/2) ID=9
      !  WRITE(ID,99) J,Z(J),OA,HA,N2A,O2A,HEA,TIJJ,TNJ,TR,RHEPHE
      !  WRITE(ID,99) J,Z(J),CHEPN,NUX(1,J),NUX(2,J)
 99   FORMAT(1X,'AVDEN9',I4,F7.0,1P,22E10.3)
      RETURN
      END
C:::::::::::::::::: ADEN11::::::::::::::::::::::::::::::
C........For N+
      SUBROUTINE ADEN11(J,OA,HA,N2A,O2A,HEA,TIJJ,TR,TNJ,NUX)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      REAL OA,HA,N2A,O2A,TNJ,TR,SQT,HEA
      DIMENSION NUX(2,FLDIM)
      SQT=SQRT(TR)

      !..  resonant ion neutral collision. No neutral N at present(8/17/89)
      NNA=0.0
      CRIN=3.83D-11*NNA*SQT*(1.0-0.063*ALOG10(TR))**2
      CNRIN=(1.45*HA+1.49*HEA+4.42*OA+7.47*N2A+7.25*O2A)
      CNRIN=CNRIN*1.0D-10
      NUX(1,J)=CRIN+CNRIN
      NUX(2,J)=0.0

      !... print diagnostics
      !  ID=3
      !  IF(J.GT.JMAX/2) ID=9
      !  WRITE(ID,99) J,Z(J),OA,HA,N2A,O2A,HEA,TIJJ,TNJ,TR,CRIN
      !  WRITE(ID,99) J,Z(J),CNRIN,NUX(1,J),NUX(2,J)
 99   FORMAT(1X,'AVDEN11',I4,F7.0,1P,22E10.3)
      RETURN
      END
C::::::::::::::::::::::: HMATRX ::::::::::::::::::::::::::::::::::::
C.... HMATRX evaluates df/dn using Stephansons method for He+ and N+.
C.... S(L,M) = DEL(FIJ) / DEL(N)
C.... where N is the density of ion IV at grid point JV. L and M
C.... are computed from JF...IV. the array RHS contains values of
C.... FIJ saved from previous computation
C.... Consult file RSLPSD-Algorithm.doc for detailed explanation
      SUBROUTINE HMATRX(FLDIM,   !.. field line grid dimension
     >                      S,   !.. Array to store the Jacobian
     >                    RHS,   !.. Stored values of function at time t
     >                   INEQ,   !.. # of equations
     >                     DT,   !.. time step for dn/dt
     >                      N,   !.. Ion density array
     >                     TI,   !.. Ion and electron temperatures
     >                   JBNN,   !.. Lower boundary index in north
     >                   JBNS,   !.. Lower boundary index in south
     >                    MIT,   !.. # of points on field line
     >                 IHEPNP,   !.. Switch for He+ or N+
     >                    XMA,   !.. mass of ion
     >                   NSPC,   !.. # of species
     >                 NMSAVE)   !.. saved density N at time t (for dn/dt)
      IMPLICIT NONE
      INTEGER FLDIM,INEQ,JBNN,JBNS,MIT,IHEPNP,NSPC   !.. see I/O comments above
      INTEGER KZS,JZS,JF,J1,J2,IV,JV,L,M,KRV,JVC,JFC,IS !.. Loop control variables
      DOUBLE PRECISION F(20)   !.. Function values at time t + delt
      DOUBLE PRECISION RHS(INEQ),S(INEQ,8),N(4,FLDIM),TI(3,FLDIM),XMA
      DOUBLE PRECISION NMSAVE(2,FLDIM),DT,H
!nm20130111
      integer ret

      DO JZS=1,INEQ        !.. loop over # of equations
      DO KZS=1,4*NSPC-1    !.. loop over band width
        S(JZS,KZS)=0.0
      ENDDO
      ENDDO

      !.. Loop to fill up the band matrix with the numerical derivative
      !.. outer loop is over the grid points
      DO JF=2,MIT-1
        J1=MAX0(2,JF-1)     !.. first column index in band
        J2=MIN0(JF+1,MIT-1) !.. last column index in band
        !.. loop over species
        DO IV=1,NSPC
          !.. loop over densities at grid points j-1, j, and j+1
          DO JV=J1,J2
            L=NSPC*(JF-2)    !.. band matrix row index
            !.. Band matrix column index
            IF(JF.LE.3) M=NSPC*(JV-2)+IV
            IF(JF.GT.3) M=NSPC*(JV-JF)+IV+2*NSPC
            KRV=2*(JV-2)+IV  !.. not used

            JVC=JBNN+JV-1   !.. grid index of the density 
            JFC=JBNN+JF-1   !.. grid index of the function

            H=1.E-4*N(IV,JVC)     !.. delta for the derivative
            N(IV,JVC)=N(IV,JVC)+H !.. increment density

            !.. Obtain the function at the new density
            CALL MDFIJ(JFC,1,NSPC,DT,N,TI,F,JBNN,JBNS,IHEPNP,XMA,NMSAVE)
            N(IV,JVC)=N(IV,JVC)-H   !.. restore density

            !.. Store derivatives in the band matrix
            DO IS=1,NSPC
              IF(JF.LE.3) S(L+IS,M)=(F(IS)-RHS(L+IS))/H
              IF(JF.GT.3) S(L+IS,M-IS)=(F(IS)-RHS(L+IS))/H
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
Cnm20130111: copied from HMATRX to distinguish between he+ v.s. n+
C::::::::::::::::::::::: NMATRX ::::::::::::::::::::::::::::::::::::
C.... NMATRX evaluates df/dn using Stephansons method for He+ and N+.
C.... S(L,M) = DEL(FIJ) / DEL(N)
C.... where N is the density of ion IV at grid point JV. L and M
C.... are computed from JF...IV. the array RHS contains values of
C.... FIJ saved from previous computation
C.... Consult file RSLPSD-Algorithm.doc for detailed explanation
      SUBROUTINE NMATRX(FLDIM,   !.. field line grid dimension
     >                      S,   !.. Array to store the Jacobian
     >                    RHS,   !.. Stored values of function at time t
     >                   INEQ,   !.. # of equations
     >                     DT,   !.. time step for dn/dt
     >                      N,   !.. Ion density array
     >                     TI,   !.. Ion and electron temperatures
     >                   JBNN,   !.. Lower boundary index in north
     >                   JBNS,   !.. Lower boundary index in south
     >                    MIT,   !.. # of points on field line
     >                 IHEPNP,   !.. Switch for He+ or N+
     >                    XMA,   !.. mass of ion
     >                   NSPC,   !.. # of species
     >                 NMSAVE)   !.. saved density N at time t (for dn/dt)
      IMPLICIT NONE
      INTEGER FLDIM,INEQ,JBNN,JBNS,MIT,IHEPNP,NSPC   !.. see I/O comments above
      INTEGER KZS,JZS,JF,J1,J2,IV,JV,L,M,KRV,JVC,JFC,IS !.. Loop control variables
      DOUBLE PRECISION F(20)   !.. Function values at time t + delt
      DOUBLE PRECISION RHS(INEQ),S(INEQ,8),N(4,FLDIM),TI(3,FLDIM),XMA
      DOUBLE PRECISION NMSAVE(2,FLDIM),DT,H

      DO JZS=1,INEQ        !.. loop over # of equations
      DO KZS=1,4*NSPC-1    !.. loop over band width
        S(JZS,KZS)=0.0
      ENDDO
      ENDDO

      !.. Loop to fill up the band matrix with the numerical derivative
      !.. outer loop is over the grid points
      DO JF=2,MIT-1
        J1=MAX0(2,JF-1)     !.. first column index in band
        J2=MIN0(JF+1,MIT-1) !.. last column index in band
        !.. loop over species
        DO IV=1,NSPC
          !.. loop over densities at grid points j-1, j, and j+1
          DO JV=J1,J2
            L=NSPC*(JF-2)    !.. band matrix row index
            !.. Band matrix column index
            IF(JF.LE.3) M=NSPC*(JV-2)+IV
            IF(JF.GT.3) M=NSPC*(JV-JF)+IV+2*NSPC
            KRV=2*(JV-2)+IV  !.. not used

            JVC=JBNN+JV-1   !.. grid index of the density 
            JFC=JBNN+JF-1   !.. grid index of the function

            H=1.E-4*N(IV,JVC)     !.. delta for the derivative
            N(IV,JVC)=N(IV,JVC)+H !.. increment density

            !.. Obtain the function at the new density
            CALL MDFIJ(JFC,1,NSPC,DT,N,TI,F,JBNN,JBNS,IHEPNP,XMA,NMSAVE)
            N(IV,JVC)=N(IV,JVC)-H   !.. restore density

            !.. Store derivatives in the band matrix
            DO IS=1,NSPC
              IF(JF.LE.3) S(L+IS,M)=(F(IS)-RHS(L+IS))/H
              IF(JF.GT.3) S(L+IS,M-IS)=(F(IS)-RHS(L+IS))/H
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
Cnm20130111
C::::::::::::::::::::::::::: DENAVE :::::::::::::::::::::::::::::::::::::
C   this program evaluates the interpolated densities at the midpoints
C   it also evaluates the ion-neutral collision frequencies nux(i,j)
C   and o+ and h+ production and loss rates
C........................................................................
      SUBROUTINE DENAVE(JMAX1,TI)
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !.. TEJ TIJ NUX UNJ DS GRADTE GRADTI GRAV OLOSS HLOSS HPROD
      USE AVE_PARAMS    !.. midpoint values - TEJ TIJ NUX UNJ etc. 
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      DIMENSION TI(3,FLDIM),RTS(99)
      DATA BK/1.3807E-16/
      DO 10 J=JMIN,JMAX-1
      DS(J)=SL(J+1)-SL(J)
      TIJ(J)=0.5*(TI(1,J)+TI(1,J+1))
      TEJ(J)=0.5*(TI(3,J)+TI(3,J+1))/DS(J)/TIJ(J)
      GRAV(J)=0.5*(GR(J)+GR(J+1))/BK/TIJ(J)
      GRADTI(J)=(TI(1,J+1)-TI(1,J))/DS(J)/TIJ(J)
      GRADTE(J)=(TI(3,J+1)-TI(3,J))/DS(J)/TIJ(J)
      UNJ(J)=0.5*(UN(J+1)+UN(J))

      !... ion-neutral coll freqs taken from schunk and nagy rev. geophys. v18
      !... p813, 1980. first,  resonant charge exchange from table 5 p823
      OA=SQRT(ON(J)*ON(J+1))
      HA=SQRT(HN(J)*HN(J+1))
      N2A=SQRT(N2N(J)*N2N(J+1))
      O2A=SQRT(O2N(J)*O2N(J+1))
      TNJ=0.5*(TN(J)+TN(J+1))
      TR=(TIJ(J)+TNJ)*0.5
      SQT=SQRT(TR)
      RHPH=2.65E-10*HA*SQT*(1-0.083*DLOG10(TR))**2

      !.. COLFAC is the scaling factor for O+ - O collision frequency  
      !.. should be 1.7 according to Burnside 1987 (in press)
      ROPO= COLFAC * 3.67E-11*OA*SQT*(1-0.064*DLOG10(TR))**2 
      RHPO=6.61E-11*OA*SQRT(TIJ(J))*(1-0.047*DLOG10(TIJ(J)))**2
      ROPH=4.63E-12*HA*SQRT(TNJ+TIJ(J)/16.0)

      !.. non-resonant ion neutral interactions table 6
      CHPN=(33.6*N2A+32*O2A)*1.0E-10
      COPN=(6.28*N2A+6.64*O2A)*1.0E-10
      NUX(1,J)=ROPO+ROPH+COPN
      NUX(2,J)=RHPH+RHPO+CHPN

      !..... O+ and H+ production and loss rates
      CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
      OLOSS(J)=RTS(2)*HN(J)+RTS(3)*N2N(J)+RTS(4)*O2N(J)
      HLOSS(J)=RTS(1)*ON(J)
      HPROD(J)=RTS(2)*HN(J)

      !... printing results
      !  IUN=34
        !... IF(J.GT.JMAX/2) IUN=9
      !  WRITE(IUN,99) J,Z(J),OA,HA,N2A,O2A,TIJ(J),TNJ,TR,RHPH,ROPO,RHPO
      !  WRITE(IUN,99) J,Z(J),ROPH,CHPN,COPN,NUX(1,J),NUX(2,J),OLOSS(J)
      !>   ,HLOSS(J),HPROD(J)
 10   CONTINUE
 99   FORMAT(1X,'AVDEN',I4,F7.0,1P,22E10.3)
      RETURN
      END
