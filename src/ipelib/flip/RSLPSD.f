C.................... RSLPSD.FOR ........................................... 
C:::::::::::::::::::::::::::::: DLOOPS :::::::::::::::::::::::::::::::::;
C--- subroutine loops is the main sequencing control program. It calls sub-
C--- prog DFIJ to obtain the error functions FIJ for the density D.E.
C--- It sets up the Jacobian matrix in subroutine DMATRIX and solves for
C--- the increments using BDSLV. Increments from the solver are
C--- tested to ensure non -ve densities (modified steepest descent).
C--- Consult file RSLPSD-Algorithm.doc for detailed explanation
      SUBROUTINE DLOOPS(JMIN,   !.. first point on the field line
     >                  JMAX,   !.. last point on the field line
     >                 FLDIM,   !.. Field line grid array dimension
     >                     Z,   !.. Altitude array
     >                     N,   !.. O+, H+, He+, minor ion densities array
     >                    TI,   !.. Ion and electron temperatures
     >                  DTIN,   !.. Time step from calling program
     >                 DTMIN,   !.. Minimum time step
     >                 EFLAG,mp,lp,nflag)   !.. OUTPUT: Error flag array
      USE SOLVARR       !... DELTA RHS WORK S, Variables for solver
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE ION_DEN_VEL   !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
!dbg20120301: N+ problem: "IN BDSLV &&&&&&& BANDWIDTH IS TOO LARGE "
      IMPLICIT NONE
      integer mp,lp,i_which_call
      INTEGER FLDIM                      !.. Field line grid array dimension
      INTEGER NFLAG,EFLAG(11,11)         !.. error flags
      INTEGER I,J,JC,ITER                !.. Loop control variables
      INTEGER IDIV,KR,ION,IEQ,MIT        !.. solution variables
      INTEGER JBNN,JBNS,NION             !.. boundary indices, # of ions
      INTEGER JMIN,JEQ,JMAX,JK,JCON      !.. spatial grid indices
      INTEGER DCLON,DCLOS,DCUPP          !.. tests for convergence region
      DOUBLE PRECISION DT,DTIN,DTMIN,DTINC !.. Time step variables
      DOUBLE PRECISION ZLBDY               !.. lower boundary altitudes for N
      DOUBLE PRECISION DCRQ(2,FLDIM),DCR(2)      !.. solution variables
      DOUBLE PRECISION DCRP,DCRT,ADCR,DINC,F(20) !.. solution variables
      !.. see above for description of Z,N,TI,HE
      DOUBLE PRECISION Z(FLDIM),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION RATIN,RATIS,FACIRI,ALPHA
      !.. N and V at previous time step for convergence failure
      DOUBLE PRECISION NORIG(2,FLDIM),VORIG(2,FLDIM)
      !.. N and T at previous intermediate time step for dN/dt and dT/dt
      DOUBLE PRECISION NSAVE(2,FLDIM),TISAV(3,FLDIM)
      DATA NION/2/  !.. # ions (i.e. O+ and H+)

      DT=DTIN           !.. Set time step for dN/dt
      DTINC=0.0         !.. Used for reduced timestep
      JBNN=JMIN         !.. lower boundary point for O+ and H+
      JBNS=JMAX         !.. lower boundary point for O+ and H+
      JEQ=(JMIN+JMAX)/2 !.. Equatorial point
      EFLAG(2,1)=0      !.. Initialize error flag
      EFLAG(2,1)=0      !.. Initialize error flag
!dbg20121130:debug note 
!i vaguely remembered ZLBDY_flip=115 caused crash (or at least very bad results!)
!i can't remember introducing ZLBDY_flip was the cause of the problem? 
! or the bad value, 115???
!therefore, the use of ZLBDY_flip is pending.
!dbg20121130      ZLBDY=115.!ZLBDY_flip        !.. Lower boundary altitude
      ZLBDY=120.        !.. Lower boundary altitude

!dbg20120301: this part was commented out on 20110815, but un-comment again to solve for N+ problem: "IN BDSLV &&&&&&& BANDWIDTH IS TOO LARGE " i don't remember why we decided to comment out here and i cannot find a program to setup the local chem equil anywhere else???
!dbg20110815      !.. Use local equilibrium for densities if flux tube apex height < 200 km
!	IF(sw_LCE.AND.Z(JEQ).LT.ht_LCE) THEN !ht_LCE=200.
!         DO J=JMIN,JMAX
!           CALL HOEQ(FLDIM,J,N,TI)
!           XIONV(1,J)=0.0
!           XIONV(2,J)=0.0
!         ENDDO
!         RETURN
!	ENDIF 
!dbg20110815:
      !.. Use local equilibrium for densities if flux tube
      !.. apex height < ZLBDY + some increment (1.0 km?)
      IF(Z(JEQ).LE.ZLBDY+1.0) THEN
         DO J=JMIN,JMAX
           CALL HOEQ(FLDIM,J,N,TI)
           XIONN(1,J)=N(1,J)
           XIONN(2,J)=N(2,J)
           XIONV(1,J)=0.0
           XIONV(2,J)=0.0
         ENDDO
         RETURN
      ENDIF 



      !-- Save current densities for dn/dt.
      !-- scaling FLIP to the IRI nmF2.
      DO J=JMIN,JMAX
      DO I=1,2
        NSAVE(I,J)=N(I,J)      !.. N for calculating dn/dt
        NORIG(I,J)=N(I,J)      !.. N to restore to previous value
        VORIG(I,J)=XIONV(I,J)  !.. Velocity at previous time step
        DCRQ(I,J)=1.0
      ENDDO
      ENDDO

      !.. AVDEN is called to obtain average midpoint densities and temps
      CALL AVDEN(TI,0,0.0D0,0.0D0)

C*** OUTER LOOP: Return here on Non-Convergence with reduced time step
  10  CONTINUE
        !.. Update indices for the lower boundaries. 
        DO J=JMIN,JEQ-1
          IF(Z(J).LE.ZLBDY) JBNN=J      !.. North
        ENDDO
        DO J=JMAX,JEQ,-1
          IF(Z(J).LE.ZLBDY) JBNS=J      !.. South
        ENDDO

        MIT=JBNS-JBNN+1       !.. Number of points on field line
        IEQ=2*(MIT-2)         !.. Number of equations to set up      
!dbg20120304:
      if ( IEQ<=2 ) then
        STOP 'STOP! sub-DLOOPS:MIT'
      end if

        !.. Main loop: On each iteration the Jacobian is formed and solved for
        !.. the increments of to add to N. 
        DO ITER=1,20
          !.. boundary conditions on density. 
          DO J=JMIN,JBNN
            CALL HOEQ(FLDIM,J,N,TI)
          ENDDO
          DO J=JBNS,JMAX
            CALL HOEQ(FLDIM,J,N,TI)
          ENDDO

          !.. Compute current FIJ values to use in calculating dF/dn
          DO J=2,MIT-1
            KR=2*(J-2)
            JC=J+JBNN-1
            CALL DFIJ(JC,0,DT,N,TI,F,NSAVE)
            RHS(KR+1)=F(1)
            RHS(KR+2)=F(2)
          ENDDO

          !.. Create the Jacobian matrix
          CALL DMATRX(FLDIM,S,RHS,IEQ,DT,N,TI,JBNN,JBNS,MIT,NION,NSAVE)

          !.. Solve the linear system with the band solver
          !.. invert the jacobian matrix 'S' in the inversion routine BDSLV.
          !.. the increments are stored in array delta in this order
          !.. x(1...n,j),x(1...n,j+1),x(1...n,j+2),....x(1...n,jmax-1)
          i_which_call = 3
          CALL BDSLV(IEQ,3,S,0,RHS,DELTA,WORK,NFLAG,
     >               mp,lp,i_which_call)

          IF(NFLAG.NE.0) THEN
            EFLAG(2,2)=-1   !.. Report problem to calling routine
            RETURN
          ENDIF

          IDIV=0
          DCR(1)=1
          DCR(2)=1
          DCRP=1.0

          !.. Test density increments 'dinc' to ensure N>0 (modified 
          !.. steepest descent)
          DO I=1,NION
          DO J=2,MIT-1
            DCRP=1.0
            JC=JBNN+J-1
            DCRQ(I,JC)=1.0
            ION=3
            IF(I.EQ.NION) ION=2
            DINC=DELTA(2*J-ION)  !.. density increment
            !.. Stop density going negative
            IF(DINC.GT.0) THEN
              IF(ABS(DINC/N(I,JC)).GT.0.9999) DCRP=0.5*ABS(N(I,JC)/DINC)
              IF((DCRP.LT.DCR(I)).AND.(Z(JC).GT.0.0)) DCR(I)=DCRP
              IF(ITER.GT.0) DCRQ(I,JC)=DCRP
            ENDIF
          ENDDO
          ENDDO

          !.. Add density increment to the array 'N' and test for
          !.. convergence when IDIV=0. That is, when every density 
          !.. increment is small
          ADCR=DMIN1(DCR(1),DCR(2))    !.. Take smallest DCR value
          DO I=1,NION
            ION=3
            IF(I.EQ.NION) ION=2
            DO J=2,MIT-1
              JC=JBNN+J-1
              DINC=DELTA(2*J-ION)
              N(I,JC)=N(I,JC)-DINC*ADCR
              IF(ABS(DINC/N(I,JC)).GT.1E-3)  IDIV=IDIV+1
            ENDDO
          ENDDO

          IF(IDIV.EQ.0) GOTO 224    !.. Has convergence has occured?
      ENDDO

        !.. END OF SOLUTION LOOP     
 224    CONTINUE

        DCLON=1  !.. North Flag for non-convergence below 200 km
        DCLOS=1  !.. South Flag for non-convergence below 200 km
        DCUPP=1  !.. Flag for non-convergence above 200 km

        !.. Testing for non-convergence to reduce time step
        IF(IDIV.EQ.0) THEN
          !============== Convergence success ========
          EFLAG(2,1)=0     
          !.. DTINC is the accumulated time since DFIJ was called with
          !.. time step DTIN
          DTINC=DTINC+DT
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)') 
     >    '  O+,H+ ',ITER,DTINC,DTIN,DT,Z(JBNN)
          !.. Final convergence, store values and return
          IF(DTINC.GE.DTIN-0.1) THEN
             DO J=JMIN,JMAX
               XIONN(1,J)=N(1,J)
               XIONN(2,J)=N(2,J)
            ENDDO
            !.. CALL DFIJ(40,4,DT,N,TI,F,NSAVE) !diagnostic
            RETURN
          ENDIF
          !.. increase time step if convergence is easy
          IF(ITER.LT.5.AND.DTINC+2*DT.LE.DTIN) DT=2*DT
          IF(ITER.LT.5.AND.DTINC+2*DT.GT.DTIN) DT=DTIN-DTINC

          !-- Save current densities for dN/dt. 
          DO J=JMIN,JMAX
            NSAVE(1,J)=N(1,J)
            NSAVE(2,J)=N(2,J)
          ENDDO
        ELSE
          !============== Convergence failure ========
          !.. Test to see if difficulty only below 200 km
          DO J=JBNN,JBNS
            IF(Z(J).LT.200.0) THEN
              IF(J.LE.JEQ.AND.DCRQ(1,J).LT.0.99999999) DCLON=0
              IF(J.GT.JEQ.AND.DCRQ(1,J).LT.0.99999999) DCLOS=0
            ELSE
              IF(DCRQ(1,J).LT.0.99999999) DCUPP=0
            ENDIF
          ENDDO
          !-- Non-Convergence: Reduce time step and restore densities. 
          DT=DT/2  !.. reduce time step for non-convergence
          DO J=JMIN,JMAX
            N(1,J)=NSAVE(1,J)
            N(2,J)=NSAVE(2,J)
!20121130 Phil suggested to try this! it might help in finding a solution for convergence error...
! GHGM the below 0.9 stuff was switched on in this version
! but Naomi always has it switched off - so I've commented it out....
! GHGM
!            if ( sw_init_guess_flip ) then
!          !.. Try a different initial guess. Could use random increment
!              IF(DT .LT. 60) THEN
!                 N(1,J)=0.9*NSAVE(1,J)
!                 N(2,J)=0.9*NSAVE(2,J) 
!              END IF
!            end if !( sw_init_guess_flip ) then
!20121130
          ENDDO
          !.. Raise lower boundary if the problem is only below 200 km
          IF(DT.LE.DTIN/4.0.AND.DCUPP.EQ.1.AND.DCLON*DCLOS.EQ.0)
     >       ZLBDY=(ZLBDY+200)/2   

          !.. Check that DT is not too small
          IF(DT.LT.DTMIN) THEN
            EFLAG(2,1)=-1   !.. Report problem to calling routine
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,9I5)') 
     >        '  ERR FLAGS LPSD',DCLON,DCLOS,DCUPP,NINT(ZLBDY)
            !-- Restore densities and velocities to original saved values. 
            DO J=JMIN,JMAX
              N(1,J)=XIONN(1,J)
              N(2,J)=XIONN(2,J)
              XIONV(1,J)=VORIG(1,J)
              XIONV(2,J)=VORIG(2,J)
            ENDDO
            RETURN
          ENDIF 
        ENDIF

      GOTO 10
C*** END OF OUTER LOOP: 

      RETURN
      END
C:::::::::::::::::::::::::: HOEQ :::::::::::::::::::::::::::::::::::::::::::::::::
C.....  finds the new chemical equilibrium densities
C.....  of O+ and H+ at point j
      SUBROUTINE HOEQ(FLDIM,   !.. Field line grid array dimension
     >                    J,   !.. point on the field line
     >                    N,   !.. O+, H+, He+, minor ion densities array
     >                   TI)   !.. Ion and electron temperatures
      USE THERMOSPHERE  !.. ON HN N2N O2N HE TN UN EHT COLFAC
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION !.. EUV, photoelectron, and auroral production, PHION
      IMPLICIT NONE
      INTEGER FLDIM   !.. field line grid dimension
      INTEGER J       !.. Field line grid array dimension
      DOUBLE PRECISION RTS(99),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION LOPLS,LOPLS2,LHPLS   !.. O+, H+ loss rates

      !.. Get reaction rates
      CALL RATS(J,TI(3,J),TI(1,J),TN(J),RTS)
      !....... evaluate loss rates
      LOPLS=(RTS(3)*N2N(J)+RTS(4)*O2N(J))
      LOPLS2=+RTS(2)*HN(J)
      LHPLS=RTS(1)*ON(J)
      N(1,J)=PHION(J)/LOPLS
      N(2,J)=N(1,J)*LOPLS2/LHPLS
       RETURN
       END
C::::::::::::::::::::::: DMATRX ::::::::::::::::::::::::::::::::::::
C.... DMATRX evaluates df/dn using Stephansons method.
C.... S(L,M) = DEL(FIJ) / DEL(N)
C.... where N is the density of ion IV at grid point JV. L and M
C.... are computed from JF...IV. the array RHS contains values of
C.... FIJ saved from previous computation
C.... Consult file RSLPSD-Algorithm.doc for detailed explanation
      SUBROUTINE DMATRX(FLDIM,   !.. field line grid dimension
     >                      S,   !.. Array to store the Jacobian
     >                    RHS,   !.. Stored values of function at time t
     >                   INEQ,   !.. # of equations
     >                     DT,   !.. time step for dn/dt
     >                      N,   !.. O+, H+, He+, minor ion densities array
     >                     TI,   !.. Ion and electron temperatures
     >                   JBNN,   !.. Lower boundary index in north
     >                   JBNS,   !.. Lower boundary index in south
     >                    MIT,   !.. # of points on field line
     >                   NSPC,   !.. # of species
     >                  NSAVE)   !.. saved density N at time t (for dn/dt)
      IMPLICIT NONE
      INTEGER FLDIM,INEQ,JBNN,JBNS,MIT,NSPC   !.. see I/O comments above
      INTEGER KZS,JZS,JF,J1,J2,IV,JV,L,M,KRV,JVC,JFC,IS !.. Loop control variables
      DOUBLE PRECISION F(20)   !.. Function values at time t + delt
      DOUBLE PRECISION RHS(INEQ),S(INEQ,8),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION NSAVE(2,FLDIM),DT,H

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
            CALL DFIJ(JFC,1,DT,N,TI,F,NSAVE)
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
C:::::::::::::::::::::::::::: SMOOTH :::::::::::::::::::::::
C.. Smoothing out O+ profile when boundary above 120 km by 
C.. fitting an exponential from absolute lower boundary to the O+ 
C.. lower boundary. This profile is then blended with the
C.. local equilibrium values from subroutine HOEQ. The second
C.. call to HOEQ is to get the H+ density consistent with O+
      SUBROUTINE SMOOTH(FLDIM,JMAX,JBNN,JBS,Z,N,DT,TI)
      IMPLICIT DOUBLE PRECISION(A-H,L,N-Z)
      INTEGER FLDIM   !.. field line grid dimension
      DIMENSION N(4,FLDIM),TI(3,FLDIM),NSAVE(4,FLDIM),Z(FLDIM)
        JINC=1
        JBS=JMAX-JBNN+1
        JBSS=JMAX-JBS+1
        AL1=DLOG(N(1,JBNN+JINC)/N(1,JBS))/(Z(JBNN+JINC)-Z(JBS))
        AL2=DLOG(N(1,JBS-JINC)/N(1,JBSS))/(Z(JBS-JINC)-Z(JBSS))
      !.......... Northern hemisphere
        DO J=JBS,JBNN+JINC
           N(1,J)=N(1,JBS)*EXP(AL1*(Z(J)-Z(JBS)))
           CALL HOEQ(FLDIM,J,NSAVE,TI)
           N(1,J)=(N(1,J)+NSAVE(1,J))*0.5
           !.. second call to get H+ with new O+ 
           CALL HOEQ(FLDIM,J,NSAVE,TI)
           N(2,J)=NSAVE(2,J)
        ENDDO

      !.......... Southern hemisphere
        DO J=JBS-JINC,JBSS
           N(1,J)=N(1,JBSS)*EXP(AL2*(Z(J)-Z(JBSS)))
           CALL HOEQ(FLDIM,J,NSAVE,TI)
           N(1,J)=(N(1,J)+NSAVE(1,J))*0.5
           !.. second call to get H+ with new O+ 
           CALL HOEQ(FLDIM,J,NSAVE,TI)
           N(2,J)=NSAVE(2,J)
        ENDDO
      RETURN
      END
