C.................... RSSUBC.FOR ................................
      SUBROUTINE RFIELD(IWR,PCO,RE,Z0,GL,SCAL,SCALK,ZLBDY,Q,VPERP,
     &             VTOT,COSDIP,JB,JR,JNQ,JTI,ZHI,ZLO,JVERT,JVERC)
C....
C.... this program is responsible for calling the subroutine that sets up
C.... the dipole magnetic field grid and filling the field parameter arrays,
C.... altitude, gravity etc.
C.... Reference: G.J. Bailey, Planet. Space Sci., Vol. 31, No. 4,
C....            pp. 389-409, 1983.
 
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL AP,DEC,ETRAN,BLON,F107,F107A,SEC,RE,GLATD,GLOND
      DIMENSION FD(9),COSDIP(401), X(401),Q(401),GL(401),SINDIP(401)
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      COMMON/VN/U(2,401),BG(401),BM(401),GR(2,401),R(2,401),SL(401)
      COMMON/ALT/Z(401),DT,DH,THF,EPS,TF,ITF,JMAX,JMAX1,ITER,ION
      common/ExB/eq_vperp,div_vem(401)

      RPTS=(JMAX+1)/2-1         ! just converting to the real
      DH=1.0/RPTS               ! distance between adjacent points in x coord.
      ALTR=370.0+130.0*F107/70. ! used in limiting flux calculation
      eq_vperp=vperp

      DO 10 J=1,JMAX
        CALL FIELD(J,JMAX,PCO,RE,SCAL,SCALK,Q(J), X(J),FD)
        JU=JMAX+1-J                ! symetric point of J about equator plane
        IF(J.EQ.1) SREF=FD(8)*1.E5 ! distance along magnetic field
        Z(J)=FD(1)                 ! altitude in km
        IF(Z(J).LE.ZLBDY)  JB =J   ! lower boundary
        IF(Z(J).LT.ALTR)   JR =J   ! index for limiting flux
        IF(Z(J).LT.ZHI) JNQ=J      ! index for plasmaspheric fluxes
        GL(J)=FD(4)                ! latitude
        BM(J)=FD(6)                ! magnetic field strength
        BG(J)=FD(7)                ! coordinate factor
        COSDIP(J)=FD(3)            ! cos(dip)
        SINDIP(J)=FD(2)            ! cos(dip)
      !..     IF(JTI.GT.1) GO TO 10
        GR(2,J)=-FD(5)             ! gravity
        R(2,J)=(RE+FD(1))*1.E5     ! radial distance to grid point
        IF(J.NE.1) SL(J)=SREF-FD(8)*1.E5   ! distance to grid point along B
      !..  Evaluate the ExB term for the continuity equation
        CALL GET_DIVDRIFT_VELOCITY(RE,PCO,X(J),VPERP, DIV_VEM(J))
        IF(DABS(GL(J)).LT.1.0D-4)     GO TO 20
      !..  Set up field line parameters in Southern hemisphere
        GL(JU)=-GL(J)
        Z(JU)=Z(J)
        BM(JU)=BM(J)
        BG(JU)=BG(J)
        COSDIP(JU)=-COSDIP(J)
        SINDIP(JU)=-SINDIP(J)
        GR(2,JU)=-GR(2,J)
        SL(JU)=SREF+FD(8)*1.E5
        R(2,JU)=R(2,J)
        DIV_VEM(JU) = DIV_VEM(J)
  10  CONTINUE
  20  CONTINUE
C
      JQ=(JMAX+1)/2
      IF(Z(JQ).LT.80) THEN
        WRITE(6,696) NINT(Z(JQ)) 
 696    FORMAT(/'   ** ERROR: Field line apex altitude is',I3,' km.'
     >      ,1X,'It does not reach the ionosphere ***'/)
        CALL RUN_ERROR    !.. print ERROR warning in output files
        STOP 
      ENDIF

      IF(ZLO.LT.Z0.OR.Z0.GE.120.) ZLBDY=Z0

       !--- Find altitude points for setting lat, long for vertical grid
      JVERT=1
      DO 245 IL=1,JMAX/2
        IF(Z(IL).LT.250) JVERT=IL
 245  CONTINUE

      JVERC=JMAX+1-JVERT   !... conjugate index

      !.. Calculate total flux tube volume
      VTOT=0.0
      DO 30 J=JB+1,JMAX-JB
        VTOT=VTOT+BM(JB)*(SL(J)-SL(J-1))/BM(J)
  30  CONTINUE
C
      IF(IWR.NE.7) RETURN
      !... overwrite previous values if flux tube moves
      REWIND (17)
      !-- write the flux tube parameters to a file
      WRITE(17,100) IDAY/1000,MOD(IDAY,1000),F107,F107A,AP(1)
     >   ,SEC/3600.0,PCO,BLON,NINT(SL(JMAX)*1.0E-5),VTOT
C
      TOTVOL=0.0
      DELS=0.0
      DELZ=0.0
      DO J=1,JMAX
        IF(J.GT.1) THEN
          DELS=(SL(J)-SL(J-1))/1.0E5
          DELZ=ABS(Z(J)-Z(J-1))
          IF(Z(J).GT.ZLBDY) 
     >     TOTVOL=TOTVOL+BM(JB)*(SL(J)-SL(J-1))/BM(J)
          ENDIF
      !...    transform magnetic to geographic coordinates---
        CALL BTOG(IEXP,BLON,SNGL(GL(J)*57.296),GLATD,GLOND)
        VELEM=0.0
        IF(J.GT.1) VELEM=BM(JB)*(SL(J)-SL(J-1))/BM(J)
        WRITE(17,110) J, Z(J), (SL(J)*1.0E-5), DELZ ,DELS,
     >   (ACOS(COSDIP(J))*57.3),COSDIP(J),BM(J),TOTVOL,GR(2,J),
     >   (R(2,J)*1.0E-5),DIV_VEM(J),GL(J)*57.296,GLATD,GLOND,VELEM
     >   ,SINDIP(J)
      ENDDO
        !... WRITE(6,'(1P,2E11.2)') div_vem(JR),TOTVOL
 110  FORMAT(I5,5F10.2,2F9.2,1P,E10.2,0P,2F9.1,1P,E10.2,0P,3F7.1,
     >  1P,E15.6,0P,F8.3)

      RETURN
 100  FORMAT(/9X,' * * * * File for field line grid parameters,
     > EUV fluxes cross sections * * * * *'
     > ,//,'  YEAR= ',I4,'    DAY= ',I3,'   F107= ',F5.1,'   F107A= '
     >   ,F5.1,'     Ap= ',F5.1
     > ,//9X,' **** The field line grid parameters at UT=',F8.2,' ****'
     > ,//,'  L-shell =',F7.3,';     Magnetic longitude= ',F5.1
     > ,/,' Length of flux tube = ',I8,' km'
     > ,'     Total flux tube volume = ',1P,E10.2,' cm-3'
     > ,//,'    pt     alt    arc_len   alt_dif   arc_dif   dip_ang'
     > ,2X,'cos_dip  mag_fld  Tube_vol  gravity   Re+alt  Div_Vem'
     > ,3X,'B_lat  G_lat G_long   Vol_elem      SinDip')
      END
C:::::::::::::::::::: GET_DIVDRIFT_VELOCITY ::::::::::::::::::::::::::::::
      SUBROUTINE GET_DIVDRIFT_VELOCITY (RE,PCO,X,VPERP, DIV_VEM)
C
C.... This subroutine calculates the ExB term in the continuty equation.
C.... This expression comes from Bailey, PSS 1983, page 392. Note that Re in
C.... his equation (5) is the equatorial radius to the flux tube
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      REAL RE
 
      NUMERATOR = 6.0000 * ((DSIN(X))**2.0)*(1.0+(DCOS(X))**2.0)
      DENOMINATOR = RE*PCO*1.0E5*((1.0+3.0*(DCOS(X))**2.0))**2.0
      DIV_VEM = NUMERATOR/DENOMINATOR
      RETURN
      END
 
C:::::::::::::::::::::::::: QFIELD ::::::::::::::::::::::::::::::::::::::::
        SUBROUTINE QFIELD(JMAX,PCO,RE,Z0,SCAL, SCALK,Q)
C
C.... This program determines Q values.
C.... JMAX = # of grid points, PCO = lshell, Z0 = lower boundary,
C.... scal = scaling factor to ensure that the x-coord ranged from 0 to 1.
C.... SCALK & Q are returned.
C.... Reference: G.J. Bailey, Planet. Space Sci., Vol. 31, No. 4,
C               pp. 389-409, 1983.
C
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      DIMENSION Q(401)
      REAL RE
 
      !... R0 = equatorial radius to flux tube.
      R0=RE*PCO
      PTS=(JMAX+1)/2
      IPTS=PTS
      RAT=(RE+Z0)/RE
      QMAX=(SQRT(1-RAT/PCO))/(RAT**2)
      SCALK=1/SINH(SCAL*(QMAX))
      DO J = 1, JMAX
        DX=1-(J-1)/(PTS-1)
      !..  q=the dipole coord  determined from -- sinh(kq)=dx
        SCX=DX/SCALK
        Q(J) = DLOG(SCX+DSQRT(SCX**2+1))/SCAL
      ENDDO
 
      RETURN
      END
 
C:::::::::::::::::::::::::: FIELD ::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE FIELD(J,JMAX,PCO,RE,SCAL,SCALK,Q, XR,FD)
C..........................................................
C  this program determines the grid point spacing given just jmax=
C   # of grid points, pco=lshell, z0=lower boundary, scal=scaling
C   factor;; x=colatitude, the field parameters are transferred
C   through fd(i);; r0=equatorial radius to flux tube, initial values
C   for x and r are set. dh=distance between points in the x coordinate
C  Reference: G.J. Bailey, Planet. Space Sci., Vol. 31, No. 4,
C               pp. 389-409, 1983.
C--------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
C     temp_store is used to store previous computed values
      DIMENSION FD(9)
      REAL RE
      X=.5
      R0=RE*PCO
      !... Values for equatorial point
      IF(J .EQ. (JMAX+1)/2) THEN
         X=1.570796327
         SHX=1
         CHX=0
      ELSE
         JITER=0

         !... loop for Newton solver
  10     CONTINUE
         SHX=DSIN(X)
         CHX=DCOS(X)
      !.. the next 6 lines are a newton solver for the equation f(x)=0
      !.. this determines the colatitude x
         FEX=(R0**2)*(SHX**4)-(RE**2)*CHX/Q
         DEX=(R0**2)*4*(SHX**3)*CHX+(RE**2)*SHX/Q
         X=X-FEX/DEX
         JITER=JITER+1
         IF(JITER.GT.99) THEN
           WRITE(6,*) 
     >      ' ** ERROR: Terminate after 100 iterations by Newton solver'
           CALL RUN_ERROR    !.. print ERROR warning in output files
           STOP
         ENDIF
         IF(ABS(FEX/(DEX*X)).GT.1.0E-6) GO TO 10
      ENDIF
      !..end newton solver

      !... XR is a returned value other than X in order to keep previous value of
      !... X, generated during last CALL, from being modified by current CALL with
      !... some random value.
      XR = X
      !..  r=radial distance to grid point. sqth=numerator of mag field
      R=R0*SHX**2
      SQTH=DSQRT(3*(CHX**2)+1)
      !..  fd(1)=altitude. fd(2)=sin(dip). fd(3)=cos(dip). fd(4)=latitude
      FD(1)=R-RE
      FD(2)=2*CHX/SQTH
      FD(3)=SHX/SQTH
      FD(4)=1.570796327-X
      !..  fd(5)=gravity. fd(6)=mag field strength. fd(7)= coord factor.
      FD(5)=FD(2)*3.98E+10/((RE+FD(1))**2)
      FD(6)=8.271E+25*SQTH/(R*1.E+5)**3
      FD(7)=SCAL*DCOSH(SCAL*Q)*SCALK*1.E-5*RE**2*SQTH/R**3
      !..  ds & fd(8)= arc length along field calc by 2 methods
      XX=DLOG(1.732*CHX+SQTH)
      FD(8)=R0*.28868*(XX+DSINH(XX)*DCOSH(XX))
 
      RETURN
      END 
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE BTOG(IEXP,BLOND,BLATD,GLATD,GLOND)
C....... Conversion from centered dipole to geographic coordinates. (PLAT
C....... ,PLON)= coords of north magnetic pole. The geomagnetic longitude
C....... is measured east from the geomagnetic prime meridian at 291 degrees
C....... east geographic.
      COMMON/GPAR/RE,PLAT,PLON,PI,ANGVEL
      !....... degrees to radians
      BLONR=BLOND/57.29578
      BLATR=BLATD/57.29578
      !...... Position of point in geomagnetic coords
      XM=COS(BLATR)*COS(BLONR)
      YM=COS(BLATR)*SIN(BLONR)
      ZM=SIN(BLATR)
      !...... rotate coords in north-south direction
      XG=XM*SIN(PLAT)+ZM*COS(PLAT)
      YG=YM
      ZG=-XM*COS(PLAT)+ZM*SIN(PLAT)
C
      !..... geographic latitude and longitude converted to degrees
      GLATD=ASIN(ZG)*57.29578
      GLOND=(PLON+ATAN2(YG,XG))*57.29578
      IF(GLOND.GE.360) GLOND=GLOND-360
      RETURN
      END
C::::::::::::::::::::::::: EPHEM ::::::::::::::::::::::::::::::::::::::
      SUBROUTINE EPHEM(IDAY,ETRAN)
      DIMENSION EPTRAN(26)
      DATA EPTRAN/4338,4374,4398,4404,4392,
     >  4374,4344,4320,4302,4296,4302
     > ,4320,4338,4356,4356,4344,4320
     > ,4290,4260,4236,4218,4224,4248,4284
     > ,4320,4350/
      ANUM=(MOD(IDAY,1000)+15.)/15.
      INUM=ANUM
      FRAC=ANUM-INUM
      ETRAN=(EPTRAN(INUM)*(1.0-FRAC)+EPTRAN(INUM+1)*FRAC)*10.
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ANFLX(JTIMAX,JR,JQ,Z,N,TI,BM,GR,ON,HN,TN,SL,FYLIM,ALSS,
     >  ZTRAN)
C........ This subroutine calculates the limiting H+ flux using the method
C........ described by Richards and Torr JGR 1985 page 5261
      DOUBLE PRECISION RTS(99),ON(401),HN(401),TN(401),BM(401),SL(401)
     > ,Z(401),N(4,401),TI(3,401),GR(2,401),H1,H2,H0,ACON,ZTRAN,HS,
     > FYLIM,ALSS
      !........ determine h+ flux by simple method
      H1=-51.5*(TI(2,JR)+TI(3,JR))/GR(2,JR)
      H2=-51.5*TN(JR)/GR(2,JR)
      CALL RATS(JR,TI(3,JR),TI(2,JR),TN(JR),RTS)
      H0=1.0E5*H1*H2/(H1-H2)
      !...... may need a factor of 2.15 to account for error in paper in 
      !...... thermal diffusion
      ACON=9E7*TI(2,JR)**2.5/RTS(1)/H0**2
      ZTRAN=Z(JR)-DLOG(ACON/N(1,JR)/ON(JR))*H1*H2/(H1+H2)
      DO 717 J=1,JQ
 717  IF(Z(J).LT.ZTRAN) JS=J
      IF(ZTRAN.LE.Z(1).OR.ZTRAN.GT.3000.) JS=JR
      HS=-51.5E5*(TI(2,JS)+TI(3,JS))/GR(2,JS)
      FYLIM=RTS(2)*HN(JS)*N(1,JS)*HS*BM(1)/BM(JS)
      ALSS=0.0D0
      DO 718 J=JS,JQ
      ALSS=ALSS+RTS(1)*N(2,J)*ON(J)*(SL(J+1)-SL(J))
 718  CONTINUE
C
      IF(JTIMAX.GT.1) RETURN
      !.. WRITE(6,95) Z(JR),TI(2,JR),TI(3,JR),H1,H2,TN(JR),H0
      !.. > ,ON(JR),N(1,JR),ACON,RTS(1)
 95   FORMAT(/'  ZJR',3X,'TI',4X,'TE',5X,'HO+',3X,'HO',3X,'TN'
     > ,6X,'H0',5X,'(O)',5X,'(O+)',4X,'ACON',5X,'R1',/6F6.0,1P,5E8.1)
C
      !..    WRITE(6,97) Z(JS),TI(2,JS),TI(3,JS),GR(2,JS),BM(1),BM(JS)
      !..   > ,RTS(2),HN(JS),N(1,JS),HS
  97   FORMAT('  ZJS   TI    TE    GR    BM1   BMJS',3X,'R2',5X,'(H)'
     > ,5X,'(O+)',5X,'HS',/4F6.0,2F6.2,1P,9E8.1)
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SYM(JMAX,HE,EHT,UN,Z,SZA)
C........... This routine gives symmetric ambient parameters
C....... NOTE WELL ! This program does not produce symmetry at night
C....... because the density at grazing incidence is not symmetric
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
      DIMENSION UN(401),Z(401),HE(401),EHT(3,401),SZA(401)
      COMMON/FJS/N(4,401),TI(3,401),F(20)
      COMMON/ND/ON(401),HN(401),N2N(401),O2N(401),PHION(401),TN(401)
      DO 25 J=1,JMAX/2
      JU=JMAX+1-J
      ON(JU)=ON(J)
      O2N(JU)=O2N(J)
      N2N(JU)=N2N(J)
      HN(JU)=HN(J)
      HE(JU)=HE(J)
      PHION(JU)=PHION(J)
      SZA(JU)=SZA(J)
      TN(JU)=TN(J)
      UN(JU)=-UN(J)
      DO 25 I=1,4
      IF(I.LE.3) TI(I,JU)=TI(I,J)
      IF(I.LE.3) EHT(I,JU)=EHT(I,J)
      IF(I.GT.0) N(I,JU)=N(I,J)
 25   CONTINUE
      UN(JMAX/2+1)=0.0D0
      RETURN
      END
C:::::::::::::::::::::::::PHEAD::::::::::::::::::::::::::::::::::::::::::
C....... Subroutine for printing basic information in each of the
C....... output files
      SUBROUTINE PHEAD(IFILE,PCO,GLONN,GLONF,GLATN,GLATF,ALT,UVFAC)
      DOUBLE PRECISION PCO,GL,ALT
      COMMON/MSIS/AP(7),DEC,ETRAN,BLON,F107,F107A,IDAY,SEC
      REAL UVFAC(59)
      CHARACTER*60 HEMIS
      IF(IFILE.EQ.3) THEN
         HEMIS = 'NORTHERN hemisphere * * * * * * * * * *'
      ELSE IF(IFILE.EQ.4) THEN
         HEMIS = 'EQUATORIAL PLANE * * * * * * * * * *'
      ELSE IF(IFILE.EQ.9) THEN
         HEMIS = 'SOUTHERN hemisphere * * * * * * * * * *'
      ELSE IF(IFILE.EQ.6) THEN
         HEMIS = 'BOTH hemispheres * * * * * * * * * *'
      ELSE
         HEMIS = ' * * * * * * * * * *'
      ENDIF
      WRITE(IFILE,112) HEMIS
 112  FORMAT(5X,'* * * * * Initial Parameters for this FLIP run in '
     > ,A60)
 114  FORMAT(/1X,' YEAR=',I4,',  DAY=',I3,',  L-SHELL=',F8.3,
     >  ',  F10.7=' ,F6.1,',  F10.7A=',F6.1,',  AP=',F6.1)
      WRITE(IFILE,114) IDAY/1000,MOD(IDAY,1000),PCO,F107,F107A,AP(1)
      
      !... print lats and longs as integers or else reals
      IF(IDAY.GT.999) THEN
         WRITE(IFILE,115) NINT(GLATN),NINT(GLONN),-NINT(GLATF)
     >  ,NINT(GLONF),NINT(ALT)
      ELSE
         WRITE(IFILE,915) GLATN,GLONN,-GLATF,GLONF,ALT
      ENDIF
 115  FORMAT(/'  Geographic (lat,long):- North (',I3,'N,',I5,'E)'
     > ,2X,',  South (',I3,'S,',I5,'E)  @',I6,'km altitude')
 915  FORMAT(/'  Geographic (lat,long):- North (',F5.1,' N,',F6.1,' E)'
     > ,',  South (',F5.1,' S,',F6.1,' E)  @',F9.2,' km altitude')

      WRITE(IFILE,116) UVFAC(9),UVFAC(6),UVFAC(18),UVFAC(33),UVFAC(35)
 116  FORMAT(/'   Using HEUVAC EUV. Initial scaling factors: 303.78='
     >  ,F5.1,', 284.15=',F5.1,  ', 584.33=',F5.1, ', 977.02=',F5.1,
     >  ', 1025.72=',F5.1)
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE DATRD1(ISOL,ITF,ITFS,I1,I2,IPP,MAXIT,JTIMAX,DTSET
     > ,DTMAX,JTIDT,TWIDT,ZLBDY,IREAD,HPEQ,FPAS,HFAC,IVLPS,TSTOP,JWR
     > ,ZLO,ZHI,IWIND,IVIBN2,IODDN,INPLUS,IHEP,HEPRAT,ITEMP,IHPOP
     > ,IHPKP)
      DOUBLE PRECISION DTSET,DTMAX,TWIDT,ZLBDY,HPEQ,FPAS,HFAC,TSTOP
     > ,ZLO,ZHI,HEPRAT
C
      
	!.. TSTOP no longer checked to be > 38 hours. December 2008
      READ(5,*) TSTOP     !.. Length of FLIP run in hours. 
      TSTOP=ABS(TSTOP)    !.. kept for safety

C------ Test to see if file exists.
      CALL TFILE(15,IOPEN)
      IF(IOPEN.EQ.0) THEN
         WRITE(6,79) 
 79      FORMAT(/1X,' ** ERROR: FLIPRUN.DDD file does not exist ***'/)
         CALL RUN_ERROR    !.. print ERROR warning in output files
         STOP
      ENDIF
C------ Read from file until a good value for ITEMP - get rid of headers
 7    READ(15,*,ERR=7) ITEMP
         IF(ITEMP.LT.0.OR.ITEMP.GT.1) ITEMP=1
         I1=ITEMP
      READ(15,*) IHPOP
         IF(IHPOP.LT.0.OR.IHPOP.GT.1) IHPOP=1
         I2=I1+IHPOP
         IF(ITEMP.EQ.0.AND.IHPOP.EQ.1) I1=2
         IF(ITEMP.EQ.0.AND.IHPOP.EQ.1) I2=2
      READ(15,*) IHEP
         IF(IHEP.LT.0) IHEP=1
      READ(15,*) INPLUS
         IF(INPLUS.LT.0) INPLUS=1
      READ(15,*) IODDN
         IF(IODDN.LT.0) IODDN=1
      READ(15,*) IVIBN2
         IF(IVIBN2.LT.0) IVIBN2=1
      READ(15,*) IDUNNO
      READ(15,*) IDUNNO

      !....... Now set IVLPS to determine which species to solve
         IVLPS=0
         IF(IODDN.GE.1) IVLPS=1
         IF(IVIBN2.GE.1) IVLPS=2
         IF(IHEP.NE.0.AND.IVLPS.EQ.0) IVLPS=-9
         IF(IHEP.NE.0.AND.IVLPS.GT.0) IVLPS=+9
         IF(INPLUS.NE.0.AND.IVLPS.LE.0) IVLPS=-11
         IF(INPLUS.NE.0.AND.IVLPS.GT.0) IVLPS=+11

      READ(15,*) DTSET
         IF(DTSET.LT.1.OR.DTSET.GT.1800) DTSET=1800
      READ(15,*) DTMAX
         IF(DTMAX.LT.1.OR.DTMAX.GT.7200) DTMAX=1800
         IF(DTSET.GT.DTMAX) DTSET=DTMAX
      READ(15,*) TWIDT
         IF(TWIDT.LT.10.OR.TWIDT.GT.3600) TWIDT=600
         IF(TWIDT.LT.300) WRITE(6,109) 
 109     FORMAT(/' **** WARNING: TWIDT < 300, Pyyddd file may be'
     >   ,' very large'/)
      READ(15,*) IWIND
      READ(15,*) HPEQ
        IHPKP=0
        IF(NINT(HPEQ).EQ.-1) IHPKP=1
        IF(HPEQ.GT.1E4) HPEQ=-HPEQ
        IF(ISOL.NE.0.AND.HPEQ.GE.0.0) HPEQ=-HPEQ
      READ(15,*) HEPRAT
        IF(HEPRAT.LT.0.001) HEPRAT=0.001
        IF(HEPRAT.GT.2.0) HEPRAT=2.0
      READ(15,*) FPAS
        IF(FPAS.LT.0.0) WRITE(6,90)
 90     FORMAT(//2X,'FPAS <0, plasmaspheric heating prop. to TEC')
        IF(FPAS.GT.1.0) WRITE(6,91)
 91     FORMAT(//2X,'FPAS >0, no conjugate e* & no plasmaspheric heat')
      READ(15,*) HFAC
        IF(HFAC.LT.0) WRITE(6,92) ABS(HFAC)
 92     FORMAT(//'  HFAC<0, multiplying ionospheric heating rate by +'
     >   ,F5.1)
      !...
      READ(15,*) JWR
         IF(JWR.LT.0) JWR=0
      READ(15,*) IPP
         IF(IPP.LT.0) IPP=0
         IF(JWR.NE.0.AND.IPP.EQ.0) IPP=1
      READ(15,*) ZLO
         IF(ZLO.LT.70.0) ZLO=300.0
      READ(15,*) ZHI
         IF(ZHI.LT.70.0) ZHI=300.0
      READ(15,*) ITF
         ITF=1         !.. ITF is always 1 since 2010-05-24
      READ(15,*) JTIMAX
         IF(JTIMAX.LT.0) JTIMAX=9999999
         !...IF(ISOL.EQ.0.AND.JTIMAX.LT.333) JTIMAX=333
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C..... This subroutine reads in the geophysical parameters for the run
      SUBROUTINE DATRD2(ISOL,IWIND,AP,BLON,F107,F107A,IDAY,SEC,UVFAC)
      REAL AP(7),UVFAC(59),ZFLUX(37)
	DATA RUNHRS/0.0/

      READ(5,*) AP(1)
         IF(ABS(AP(1)).GT.400) THEN
            WRITE(6,*) ' ** ERROR: AP is too big = ',AP(1)
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ENDIF

      READ(5,*) F107
         IF(F107.LT.60.OR.F107.GT.500) THEN
             WRITE(6,*) ' ** ERROR: F107 is out of range =',F107
             CALL RUN_ERROR    !.. print ERROR warning in output files
             STOP
         ENDIF
      READ(5,*) F107A
         IF(F107A.LT.60.OR.F107A.GT.500) THEN
             WRITE(6,*) ' ** ERROR: F107A is out of range =',F107A
             CALL RUN_ERROR    !.. print ERROR warning in output files
             STOP
         ENDIF
      READ(5,*) JDAY
         IF(ISOL.EQ.0.OR.JDAY.GT.0) IDAY=IABS(JDAY)
         IF(IDAY.GT.0.AND.IDAY.LT.367) IDAY=IDAY+2000000
         IF(IDAY.LT.1900000) IDAY=IDAY+1900000
         IF(IDAY.LT.1000000.OR.IDAY.GT.9999999) THEN
           WRITE(6,*) ' ** ERROR: Day number out of range = ',IDAY
           CALL RUN_ERROR    !.. print ERROR warning in output files
           STOP
         ENDIF
      !...... The universal time (SEC) may be changed in DATRD3
      READ(5,*) STRTHRS
         IF(ISOL.EQ.0)  SEC=STRTHRS*3600

      !....... Read in the UV factors
      READ(15,*) (UVFAC(I),I=1,37)
      READ(15,*) (UVFAC(I),I=38,50)
      READ(15,*) (UVFAC(I),I=51,59)

      CALL FACEUV(F107,F107A,UVFAC,ZFLUX)

      CALL FACSR(UVFAC,F107,F107A)

      !... check the magnitude of the euv fluxes ....
      DO 67 I=1,50
 67   IF(UVFAC(I).GT.99.OR.UVFAC(I).LT.0) IUV=1
      IF(IUV.EQ.1) THEN
        WRITE(6,*)
     >  '  **** EUV or S-Runge factors may be wrong - check FLIPRUN.DDD'
        WRITE(6,751) (UVFAC(I),I=1,50)
 751    FORMAT(10F6.1)
      ENDIF

      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE DATRD3(ISOL,BLON,SEC,DTSET,HPEQ,F107,JMAX,PCO,Z0
     > ,SCAL,GRD,L_STAG,L_SHORT,L_PHASE,IDAY)
      INTEGER IDAY
      DOUBLE PRECISION PCO,Z0,SCAL,DTSET,HPEQ,XHOLD
      REAL L_STAG,L_SHORT,L_PHASE,GLOND,GLATD,GLONP,DLAT,DLON
      COMMON/GPAR/RE,PLAT,PLON,PI,ANGVEL
      !...... If not a new field line then read dummy variables so they
      !...... are not altered
      IF(ISOL.NE.0) THEN
         READ(15,*) DUMVAR
         READ(5,*)  L_SHORT
         READ(5,*)  DUMVAR
         READ(5,*)  L_STAG
         READ(5,*)  L_PHASE
         READ(15,*)  DUMVAR
         READ(15,*)  DUMVAR
         READ(15,*)  DUMVAR
         RETURN
      ENDIF
C 
      !....... else setting up a new field line, input parameters
      !...... JMAX is the number of grid points on field line
      READ(15,*) JMAX
      !..... PCO is the L-SHELL. short for P coordinate
      READ(5,*) PCO
      !..... BLON is magnetic longitude
      READ(5,*) BLON
C
      !.. If PCO large it must be geographic latitude need to convert
      !.. to magnetic coordinates using GTOB
      IF(PCO.LT.1.02.OR.PCO.GT.9) THEN
         GLATD=PCO
         GLOND=BLON
         PLAT=-1
         READ(37,*,END=199,ERR=199) PLAT,PLON  !... for testing only
         WRITE(6,'(2F6.1)') PLAT,PLON
         PLAT=PLAT/57.296
         PLON=PLON/57.296
 199     IF(PLAT.LT.0) THEN
            !.. locate magnetic pole
            CALL MPOLE_LOC(GLOND,GLATD,IDAY,PLAT,PLON) 
         ENDIF
         CALL GTOB(IEXP,BLON,BLAT,GLATD,GLOND)
         PCO=((RE+250.)/RE)/(COS(BLAT/57.296))**2
      ELSE
         GLONP=BLON-74.5
         IF(GLONP.LT.0) GLONP=GLONP+360
         !.. locate magnetic pole
         CALL MPOLE_LOC(GLONP,50.0,IDAY,PLAT,PLON)
      ENDIF

      IF((BLON-55.0)**2.0+(GLATD+30.0)**2.0.LT.25.0*25.0)
     >  WRITE(6,'(4X,A,/4X,A,/4X,A)') 
     >  'WARNING: This flux Tube is in the South Atlantic Anomaly.',
     >  'The displaced dipole may not be a good representation.',
     >  'Not a serious problem for the local ionosphere.'
     
      !.. stagnation point and phase for tear drop model
      READ(5,*) L_STAG
      !.. L_SHORT is where the trajectory cuts the short axis
      IF(L_STAG.GT.1.1) L_SHORT=PCO
      READ(5,*) L_PHASE

      !-- correct JMAX if illegal value and make sure it is an odd number
      IF(JMAX.LT.150) JMAX=55*PCO+125
      JMAX= 2 * (JMAX/2) + 1
      IF(JMAX.LT.299) JMAX=299
      IF(JMAX.GT.399) JMAX=399
      IF(NINT(L_STAG).NE.0) JMAX=399   ! max pts for ExB drift

      !.. Z0 is the absolute lower boundary altitude
      READ(15,*) Z0
         !.. reset Z0 for ExB drift
         IF(INT(L_STAG).EQ.0) THEN
            IF(Z0.LT.80) WRITE(6,*) ' *** Z0 reset to 80 km'
            IF(Z0.LT.80) Z0=80
         ELSE
            Z0=90
            WRITE(6,97) NINT(Z0)
 97         FORMAT('  *** Z0 reset to',I4,' km for ExB drift')
         ENDIF
         IF(Z0.GT.190.0) THEN
            Z0=190.0
            WRITE(6,*) ' Z0 too high - reset to 190 km'
         ENDIF
C
      !....... SCAL distributes the grid points along B
      READ(15,*) SCAL
         IF(SCAL.LT.0.5.OR.SCAL.GT.11) SCAL=(8.0/PCO+2.0)
         IF(NINT(L_STAG).NE.0) SCAL=3.5   !.. SCAL for ExB drift
      !........ switch for variable (-1) or constant altitude steps in vertical
      READ(15,*) GRD
      !..   IF(GRD.LT.0) GRD=-1
C
      !...... Set up default starting time at noon by finding the approx. UT
      !...... IF  ABS(SEC) < 24 assume it is LT in hours to start  
       IF(SEC.LT.0) CALL SET_NOON(SEC,PCO,BLON)
C
 95    FORMAT(I6,9F10.3)
      !...... Default initial time step for new field line
      XHOLD=40*(250/F107)**1.5
      IF(DTSET.GT.XHOLD) DTSET=INT(40*(250/F107)**1.5)

      RETURN
      END

C::::::::::::::::::::::::::::: DATRD4 ::::::::::::::::::::::::::::::::
C........ Used for reading additional run control parameters
      SUBROUTINE DATRD4(ISOL,TIMEON,TIMEOF,AURFAC)
      IMPLICIT DOUBLE PRECISION(A-H,N,O-Z)
	INTEGER VD
	PARAMETER (VD=181)    !.. VD = dimension of vertical grid
      REAL TIMEON,TIMEOF,AURFAC,AMP_RING,RC_TIME_CON,TMAX,SIGT
     >   ,LMAX,SIGL,PROP
      !....... COLFAC= O+ - O collision frequency factor. CFAC, and ISCH
      !....... are for future expansion
      COMMON/CFACT/DHDU_FAC,COLFAC,CFAC(9),ISCH(9)
      COMMON/RC/AMP_RING,TMAX,SIGT,LMAX,SIGL,RC_TIME_CON,PROP
	!... common DD transfers EDDY to CALODDN and CALVN2
      COMMON/DD/D(7,VD),E(7,VD),B1(7,VD),B2(7,VD),RDZ(VD),EDDY

      !...... This is the O+ - O collision factor
      READ(15,*) COLFAC
         IF(COLFAC.LT.0) COLFAC=1.0
         IF(COLFAC.LT.0.5.OR.COLFAC.GT.4.0) THEN
            WRITE(6,*) ' ** ERROR: COLFAC out of order = ',COLFAC
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ENDIF

      !..READ(15,*) DHDU_FAC   NOW set in routine HMF2W
      READ(15,*) EDDY
         IF(EDDY.GT.1.0E7) THEN
            WRITE(6,*) ' ** ERROR: EDDY diff. coefficient too big',EDDY
            CALL RUN_ERROR    !.. print ERROR warning in output files
            STOP
         ENDIF

      !.. Amplitude of the ring current heating
      READ(15,*) AMP_RING
      AMP_RING=AMP_RING*5.0E-05
      IF(AMP_RING.LE.0.0) AMP_RING=0.0

      !.. time of maximum and sigmas for ring current heating
      READ(15,*) TMAX,SIGT
      READ(15,*) LMAX,SIGL

      !... proportion of ring current heating going to electrons
      READ(15,*) PROP
      IF(PROP.GT.1.0) PROP=1.0
      IF(PROP.LT.0.0) PROP=0.0

      READ(15,*) RC_TIME_CON     !.. decay time for ring current (hrs)

      !.. input on time and off time for aurora and energy factor
      READ(15,*) TIMEON,TIMEOF
      READ(15,*) AURFAC
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CENFOR(IEXP,BLOND,BLATD,ALT,CF)
C....... Calculation of centrifugal force. Note that it is asymmetrical with
C....... respect to the geomagnetic equator
      DOUBLE PRECISION ALT,BLATD,CF,RP,CL,DIP,CHI
      COMMON/GPAR/RE,PLAT,PLON,PI,ANGVEL
      !....... Find geographic latitude and longitude
      CALL BTOG(IEXP,BLOND,SNGL(BLATD),GLATD,GLOND)
      RP=(RE+ALT)*1.0E5
      !..... rotate the centrifugal force vector in longitude
      XG=COS(GLATD/57.296)*COS(GLOND/57.296-PLON)
      YG=COS(GLATD/57.296)*SIN(GLOND/57.296-PLON)
      ZG=0.0
      !...... rotate in latitude
      XM=XG*SIN(PLAT)-ZG*COS(PLAT)
      YM=YG
      ZM=XG*COS(PLAT)+ZG*SIN(PLAT)
C
      CL=(90.0-BLATD)/57.296
      DIP=ASIN(2*DCOS(CL)/DSQRT(3.0*DCOS(CL)**2+1.0))
      CHI=CL-DIP
      !....... Centrifugal force =CF
      CF=RP*ANGVEL*(DCOS(CHI)*(COS(BLOND/57.296)*XM+SIN(BLOND/
     > 57.296)*YM)-DSIN(CHI)*ZM)
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ZATOLT(CHI,G,D,E,TIME)
      !...... This routine converts solar zenith angle to local time
      H=((COS(CHI)-SIN(G)*SIN(D))/(COS(G)*COS(D)))
      IF(ABS(H).LE.1.0) THEN
         TIME=((ACOS(H)-E)*57.296/15.0)+12.0
      ELSE IF(H.LT.1.0) THEN
         TIME=24.0
      ELSE
         TIME=12.0
      ENDIF
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LTTOZA(TIME,GLAT,DECL,CHI)
      !...... This routine converts solar apparent time to solar zenith angle
      HH=(TIME-12.)*15*3.14159/180.
      CHI=ACOS(COS(GLAT)*COS(DECL)*COS(HH)+SIN(GLAT)*SIN(DECL))
      RETURN
      END
C:::::::::::::::::::::::::::: F2PEAK :::::::::::::::::::::::::::::::
C-- This subroutine determines the F2 peak electron density (NMF2)
C-- and height (HMF2) by fitting a parabola to the 3 points near
C-- the peak D1, D2, and D3 are the densities at H1, H2, and H3
 
      SUBROUTINE F2PEAK(D1,D2,D3,H1,H2,H3,NMF2,HMF2)
      DOUBLE PRECISION D1,D2,D3,H1,H2,H3,NMF2,HMF2
C--    Determine the coefficients A, B, C of the equation
C--    DEN = A * H**2 + B * H + C
C
C--  The first part takes care of low L values where the highest density
C--  may be at the equatorial point
      IF(NINT(H1).EQ.NINT(H3)) THEN
          NMF2=D2
          HMF2=H2
          RETURN
      ENDIF
C-- determine NMF2 and HMF2 for good data
      ANUMER = (D1-D2)*(H2-H3) - (D2-D3)*(H1-H2)
      ADENOM = (H1*H1-H2*H2)*(H2-H3) - (H2*H2-H3*H3)*(H1-H2)
      A= ANUMER/ADENOM
      B= (D1-D2 - A*(H1*H1-H2*H2))/(H1-H2)
      C= D1 - A*H1*H1 - B*H1
C--     Now the peak height and density
      HMF2 = -0.5*B/A
      NMF2 = A*HMF2*HMF2 + B*HMF2 + C
      RETURN
      END       
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GTOB(IEXP,BLOND,BLATD,GLATD,GLOND)
C........ Conversion from geographic to centered dipole coordinates
      COMMON/GPAR/RE,PLAT,PLON,PI,ANGVEL
      !....... convert to radians
      GLAT=GLATD/57.29578
      GLON=GLOND/57.29578
      !..... rotate in longitude
      XG=COS(GLAT)*COS(GLON-PLON)
      YG=COS(GLAT)*SIN(GLON-PLON)
      ZG=SIN(GLAT)
      !...... rotate in latitude
      XM=XG*SIN(PLAT)-ZG*COS(PLAT)
      YM=YG
      ZM=XG*COS(PLAT)+ZG*SIN(PLAT)
      !...... geomagnetic latitude and longitude in degrees
      BLATD=ASIN(ZM)*57.29578
      BLOND=(ATAN2(YM,XM))*57.29578
      IF(BLOND.LT.-0.01) BLOND=BLOND+360
      RETURN
      END
C::::::::::::::::::::::::: AMBBLK ::::::::::::::::::
      BLOCK DATA AMBBLK
      COMMON/AMBCOM/SECSAV,NY,NDY,LY,SX,ID,DECL
      DATA SECSAV,NY,NDY,LY,  SX,ID,DECL
     >     /0.0  , 0, 0 ,0 , 0.0, 0, 0.0  /
      END
C::::::::::::::::::::::::::: SOLDEC :::::::::::::::::::::::::
C...... This routine calculates the solar declination (DELTA radians)
C...... and ephemeris transit time (ETRAN in secs). ETRAN is the time
C...... the sun is overhead in Greenwich
C...... REFERENCE - page 484 "Explanatory Supplement to the Astronomical
C...... Almanac" Kenneth Seidelmann, University Science Books, 20 Edgehill 
C...... Road, Mill Valley, CA 94941, 1992
C...... IDAY is yyyyddd (eg 1988015 = feb 15), UT is the universal time
C...... in hours (0-24) 
       SUBROUTINE SOLDEC(IDAY,UT,DELTA,ETRAN)
      !.... Find day of year (0-366), and year (1996), then Julian day
       INTEGER JD,DAYNUM
       REAL T,YEAR,UT,L,G,LAMBDA,EPSIL,E,GHA,DELTA,SD,ETRAN,DTOR
       DOUBLE PRECISION DJD,DUT
       DATA DTOR/57.29578/   ! degrees to radians
      !..... Recover year and date and make sure UT is less than 24
       DAYNUM=MOD(IDAY,1000)
       YEAR=(IDAY-DAYNUM)/1000
       JD=INT(365.25*(YEAR-1900)+DAYNUM)+2415020      ! Julian day
       DJD=JD                                ! extra precision needed         
       DUT=UT
       T=(DJD+DUT/24.0-2451545.0)/36525.0         ! # of centuries
      !..     WRITE(6,*) DAYNUM,YEAR,JD,UT
       L=AMOD(280.460+36000.770*T,360.0)              ! aberration
       G=AMOD(357.528+35999.050*T,360.0)              ! mean anomaly
      !...... LAMBDA= ecliptic longitude. DELTA=obliquity of the ecliptic.
       LAMBDA=AMOD(L+1.915*SIN(G/DTOR)+0.020*SIN(2.0*G/DTOR),360.0)
       EPSIL=23.4393-0.01300*T
      !...... Equation of time. Time difference between noon and overhead sun
       E=-1.915*SIN(G/DTOR)-0.020*SIN(2.0*G/DTOR)
     >   +2.466*SIN(2*LAMBDA/DTOR)-0.053*SIN(4*LAMBDA/DTOR)
       GHA=15*UT-180+E                    ! Greenwich hour angle
       DELTA=ASIN(SIN(EPSIL/DTOR)*SIN(LAMBDA/DTOR))  ! solar declination
       SD=0.267/(1-0.017*COS(G/DTOR))    ! 
       ETRAN=(12-E/15)*3600 ! Ephemeris transit time in secs for FLIP
       RETURN
       END

C::::::::::::::::::::::::::::::::: SET_NOON ::::::::::::::::
C...... This subroutine resets the UT so that a run begins at noon
          SUBROUTINE SET_NOON(SEC     !.. Universal time in secs
     >                        ,PCO    !.. L shell
     >                        ,BLON)  !.. magnetic longitude
          IMPLICIT NONE
          DOUBLE PRECISION PCO
          REAL SEC,BLON,STIME,BLAT,GLATD,GLON,HRS
          INTEGER IEXP
          STIME=12.0
          IF(ABS(SEC).LT.24) STIME=ABS(SEC)
      !.. geographic longitude first using BTOG
          BLAT=57.3*ACOS(1.0/SQRT(PCO))
          CALL BTOG(IEXP,BLON,BLAT,GLATD,GLON)
          HRS=(STIME-(GLON)/15)
          IF(HRS.GT.24) HRS=HRS-24
          IF(HRS.LT.0) HRS=HRS+24
          IF(HRS.GE.23.0) HRS=0.1    ! ensures a longer settling period   
          SEC=HRS*3600
          !...WRITE(6,'(9F7.2)') SEC/3600.,BLON,BLAT,GLATD,GLON
          RETURN
          END