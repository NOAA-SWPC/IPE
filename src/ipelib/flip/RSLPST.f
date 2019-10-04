C.................... RSLPST.FOR ........................................... 
C:::::::::::::::::::::::::::::: TLOOPS :::::::::::::::::::::::::::::::::;
C--- subroutine loops is the main sequencing control program for the temperatures.
C--- It calls subroutine TFIJ to obtain the error functions FIJ for the density
C--- and temperature d.e.'s. It sets up the Jacobian matrix in subroutine
C--- TMATRIX and solves for the increments using BDSLV. increments from
C--- the solver are tested to ensure non -ve densities (modified steepest
C--- descent). P. Richards May 2010
C.... Consult file RSLPST-Algorithm.doc for detailed explanation
      SUBROUTINE TLOOPS(JMIN,   !.. first point on the field line
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
      IMPLICIT NONE
      integer mp,lp,i_which_call
      INTEGER FLDIM                !.. Field line grid array dimension
      INTEGER NFLAG,EFLAG(11,11)         !.. error flags
      INTEGER I,J,JC,NTI,ITER            !.. Loop control variables
      INTEGER IDIV,KR,ION,IEQ,MIT        !.. solution variables
      INTEGER JBNN,JBNS,JBTN,JBTS        !.. boundary indices
      INTEGER JMIN,JEQ,JMAX              !.. spatial grid indices
      DOUBLE PRECISION DT,DTIN,DTMIN,DTINC !.. Time step variables
      DOUBLE PRECISION bound_alt_den,bound_alt_temp       !.. lower boundary altitudes for N & T
      DOUBLE PRECISION DCRP,DCRT,DINC,F(20) !.. solution variables
      !.. see above for description of Z,N,TI,HE
      DOUBLE PRECISION Z(FLDIM),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION NSAVE(2,FLDIM),TISAV(3,FLDIM),TIORIG(3,FLDIM)
      DATA bound_alt_den/120.0/,bound_alt_temp/101.0/
      DATA NTI/3/             !.. # of temperatures note T(O+)=T(H+)

      DT=DTIN           !.. Set time step for dT/dt
      DTINC=0.0         !.. Used for reduced timestep
      JBTN=1            !.. lower boundary point for Te and Ti
      JBTS=JMAX         !.. lower boundary point for Te and Ti
      JBNN=1            !.. lower boundary point for O+ and H+
      JBNS=JMAX         !.. lower boundary point for O+ and H+
      JEQ=(JMIN+JMAX)/2 !.. Equatorial point
      EFLAG(1,1)=0      !.. Initialize error flag
      EFLAG(1,2)=0      !.. Initialize error flag

     
      !.. Update indices for the lower boundaries. 
      DO J=JMIN,JEQ-1
        IF(Z(J).LE.bound_alt_den) JBNN=J
        IF(Z(J).LE.bound_alt_temp) JBTN=J
      ENDDO
      DO J=JMAX,JEQ,-1
        IF(Z(J).LE.bound_alt_den) JBNS=J
        IF(Z(J).LE.bound_alt_temp) JBTS=J
      ENDDO

      !.. Use Tn for Te and Ti if flux tube apex height < lower boundary.
	IF(Z(JEQ).LE.bound_alt_temp) THEN 
        DO J=JMIN,JMAX
          TI(1,J)=TN(J)
          TI(2,J)=TN(J)
          TI(3,J)=TN(J)
        ENDDO
        RETURN
      ENDIF

      !-- Save current temperatures for dT/dt. 
      DO J=JMIN,JMAX
      DO I=1,3
        TISAV(I,J)=TI(I,J)    !.. Temp at previous intermediate time step
        TIORIG(I,J)=TI(I,J)   !.. Temp at previous original time step
      ENDDO
      ENDDO

      !.. Set up boundary indices for the solution procedure
      MIT=JBTS-JBTN+1       !.. Number of points on field line
      IEQ=2*(MIT-2)         !.. Number of equations to set up      

C*** OUTER LOOP: Return here on Non-Convergence with reduced time step
  10  CONTINUE
        !.. Main loop: On each iteration the Jacobian is formed and solved for
        !.. the increments of to add to TI. 
        DO 220 ITER=1,20
          !... boundary conditions on temperature in North
          DO J=JMIN,JBTN
             TI(1,J)=TN(J)
             TI(2,J)=TN(J)
             TI(3,J)=TN(J)
          ENDDO
          !... boundary conditions on temperature in South
          DO J=JBTS,JMAX
             TI(1,J)=TN(J)
             TI(2,J)=TN(J)
             TI(3,J)=TN(J)
          ENDDO
  
          !.. boundary conditions on density. Needed here for temperatures
          DO J=JMIN,JBNN
            CALL HOEQ(FLDIM,J,N,TI)
          ENDDO
          DO J=JBNS,JMAX
            CALL HOEQ(FLDIM,J,N,TI)
         ENDDO

          !.. Compute FIJ values to use in calculating dF/dT
          DO J=2,MIT-1
            KR=2*(J-2)
            JC=J+JBTN-1
            CALL TFIJ(JC,0,DT,N,TI,F,TISAV,mp,lp)
            RHS(KR+1)=F(1)
            RHS(KR+2)=F(2)
          ENDDO

          !.. Create the Jacobian matrix. Note NTI-1 because only TE 
          !.. and one TI solved
          CALL TMATRX(FLDIM,S,RHS,IEQ,DT,N,TI,
     >                 JBTN,JBTS,MIT,NTI-1,TISAV,mp,lp)

          !.. Solve the linear system with the band solver
          !.. invert the jacobian matrix 'S' in the inversion routine BDSLV.
          !.. the increments are stored in array delta in this order
          !.. x(1...n,j),x(1...n,j+1),x(1...n,j+2),....x(1...n,jmax-1)
          i_which_call = 4
          CALL BDSLV(IEQ,3,S,0,RHS,DELTA,WORK,NFLAG,mp,lp,i_which_call)
 

          IF(NFLAG.NE.0) THEN
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(/A,I5,A,I5)')
     >        '  ITERATION =',ITER,'  RETURN FROM BDSLV =',NFLAG
            EFLAG(1,2)=-1   !.. Report problem to calling routine
            RETURN
          ENDIF

          IDIV=0
          DCRP=1.0
          DCRT=1.0

          !.. test TI increments 'DINC' to ensure TI>0 (mod steep descent)
          DO 137 I=1,NTI
          DO 137 J=2,MIT-1
            DCRP=1.0
            ION=3
            IF(I.EQ.NTI) ION=2
            DINC=DELTA(2*J-ION)
            !.. if DINC exceeds fraction of TI set the factor DCRT so that TI>0
            JC=JBTN+J-1
            IF(ABS(DINC/TI(I,JC)).GT.0.9999) DCRP=0.5*ABS(TI(I,JC)/DINC)
            IF(DCRP.LT.DCRT) DCRT=DCRP
 137      CONTINUE

          !. add iterative increment to the array 'TI' and test for
          !. convergence when idiv=0. 
          DO I=1,NTI
            ION=3
            IF(I.EQ.NTI) ION=2
            DO J=2,MIT-1
              JC=JBTN+J-1
              DINC=DELTA(2*J-ION)
              TI(I,JC)=TI(I,JC)-DINC*DCRT
              IF(ABS(DINC/TI(I,JC)).GT.1E-3)   IDIV=IDIV+1
            ENDDO
          ENDDO
          IF(IDIV.EQ.0) GO TO 224   !. test for convergence
 220    CONTINUE

        !.. END OF SOLUTION LOOP     
 224    CONTINUE

        !.. Testing for convergence to see if need to reduce time step
        IF(IDIV.EQ.0) THEN
          !============== Convergence success ========
          EFLAG(1,1)=0     
          DTINC=DTINC+DT   !.. Used for reduced timestep
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)') 
     >       '  Te Ti ',ITER,DTINC,DTIN,DT,Z(JBTN)
          IF(DTINC.GE.DTIN-1) RETURN
          !.. increase time step if convergence is easy
          IF(ITER.LT.5.AND.DTINC+2*DT.LE.DTIN) DT=2*DT
          IF(ITER.LT.5.AND.DTINC+2*DT.GT.DTIN) DT=DTIN-DTINC

          !-- Save current temperatures for dt/dt. 
          DO J=JMIN,JMAX
            TISAV(1,J)=TI(1,J)
            TISAV(2,J)=TI(2,J)
            TISAV(3,J)=TI(3,J)
          ENDDO
        ELSE
        !============== Convergence failure ========
          IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,I5,9F14.2)') 
     >     ' Temp ',ITER,DTINC,DTIN,DT
          !-- Reduce time step and restore TI. 
          DT=DT/2  !.. reduce time step for non-convergence
          DO J=JMIN,JMAX
            TI(1,J)=TISAV(1,J)
            TI(2,J)=TISAV(2,J)
            TI(3,J)=TISAV(3,J)
          ENDDO
          !.. Check that DT is not too small
          IF(DT.LT.DTMIN) THEN
            EFLAG(1,1)=-1   !.. Report problem to calling routine
            IF(EFLAG(11,11).EQ.1) WRITE(6,'(A,9I5)') '  ERROR IN LPST'
            !-- Restore temperatures to original stored values. 
            DO J=JMIN,JMAX
              TI(1,J)=TIORIG(1,J)
              TI(2,J)=TIORIG(2,J)
              TI(3,J)=TIORIG(3,J)
            ENDDO
            RETURN
          ENDIF 
        ENDIF

      GOTO 10   !... END OF OUTER LOOP ....................

      RETURN
      END
C::::::::::::::::::::::: TMATRX ::::::::::::::::::::::::::::::::::::
C.... TMATRX evaluates dF/dn using Stephansons method.
C.... S(L,M) = DEL(FIJ) / DEL(N)
C.... where N is the density of ion IV at grid point JV. L and M
C.... are computed from JF...IV. the array RHS contains values of
C.... FIJ saved from previous computation
C.... Consult file RSLPST-Algorithm.doc for detailed explanation
      SUBROUTINE TMATRX(FLDIM,   !.. field line grid dimension
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
     >                  TISAV,mp,lp)   !.. saved density N at time t (for dn/dt)
      IMPLICIT NONE
      INTEGER mp,lp
      INTEGER FLDIM,INEQ,JBNN,JBNS,MIT,NSPC   !.. see I/O comments above
      INTEGER KZS,JZS,JF,J1,J2,IV,JV,L,M,KRV,JVC,JFC,IS !.. Loop control variables
      DOUBLE PRECISION F(20)   !.. Function values at time t + delt
      DOUBLE PRECISION RHS(INEQ),S(INEQ,8),N(4,FLDIM),TI(3,FLDIM)
      DOUBLE PRECISION TISAV(3,FLDIM),DT,H

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
        DO IV=1,2
          !.. loop over temperatures at grid points j-1, j, and j+1
          DO JV=J1,J2
            L=2*(JF-2)    !.. band matrix row index
            !.. Band matrix column index
            IF(JF.LE.3) M=2*(JV-2)+IV
            IF(JF.GT.3) M=2*(JV-JF)+IV+4
            KRV=2*(JV-2)+IV  !.. not used

            JVC=JBNN+JV-1   !.. grid index of the temperature
            JFC=JBNN+JF-1   !.. grid index of the function

            !.. Note that IV+1 because solving for TI(2,J) and TI(3,J)
            H=1.E-4*TI(IV+1,JVC)   !.. delta for the derivative
            TI(IV+1,JVC)=TI(IV+1,JVC)+H !.. increment temperature

            !.. Obtain the function at the new temperature
            CALL TFIJ(JFC,1,DT,N,TI,F,TISAV,mp,lp)
            TI(IV+1,JVC)=TI(IV+1,JVC)-H   !.. restore temperature

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
