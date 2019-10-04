!nm20110817: commented out the driver part to use just the subroutine FLIP_GRID
!C.................... RSMNSD.FOR ..................................... 
!C.... this is the main module for testing the FLIP field line grid
!      IMPLICIT NONE
!      INTEGER J,JMAX,IDAY,FLDIM
!      PARAMETER(FLDIM=9001)     !.. Array dimensions
!      DOUBLE PRECISION PCO,BLON,Z0,SCAL,RE,VTOT
!      DOUBLE PRECISION GLATD(FLDIM),GLOND(FLDIM)
!      DOUBLE PRECISION COSDIP(FLDIM),Q(FLDIM),GL(FLDIM)
!      DOUBLE PRECISION Z(FLDIM),SL(FLDIM),BM(FLDIM),GR(FLDIM),R(FLDIM)
!      DOUBLE PRECISION SINDIP(FLDIM)
!
!      CALL OPEN_FILE() !.. Used to assign files on the PC only
!
!      !.. # of grid points on field line (set negative for default)
!      READ(5,*) JMAX
!      READ(5,*) IDAY      !.. Day of year in the form YYYYddd	 
!      READ(5,*) PCO       !.. L-shell or geographic latitude
!      READ(5,*) BLON      !.. Magnetic or geographic longitude
!      READ(5,*) Z0        !.. Lower boundary altitude
!      !.. # of grid points on field line (set negative for default)
!      READ(5,*) SCAL      !.. SCAL distributes the grid points along B
!
!      !..  Call the Field line grid routine
!      CALL FLIP_GRID(FLDIM,JMAX,IDAY,PCO,BLON,Z0,SCAL,RE,GL,Q,COSDIP,
!     >   SINDIP,Z,SL,GLATD,GLOND,BM,GR,R,VTOT)
!
!      !.. Write the field line grid parameters
!      CALL DATAWRITE(JMAX,IDAY,BLON,PCO,RE,GL,Q,COSDIP,SINDIP,
!     >   Z,SL,GLATD,GLOND,BM,GR,R,VTOT)
!
!      STOP
!      END
C:::::::::::::::::::::::: FLIP_GRID :::::::::::::::::::::::::::
C.... The main sequencing module for the FLIP field line grid
      SUBROUTINE FLIP_GRID(FLDIM,JMAX,IDAY,PCO,BLON,Z0,SCAL,RE,GL,Q,
     >   COSDIP,SINDIP,Z,SL,GLATD,GLOND,BM,GR,R,VTOT)
      IMPLICIT NONE
      INTEGER J,JMAX,IDAY,FLDIM
      DOUBLE PRECISION BLON,RE,Z0,VTOT
      DOUBLE PRECISION GLATD(FLDIM),GLOND(FLDIM),PLAT,PLON
      DOUBLE PRECISION COSDIP(FLDIM),Q(FLDIM),GL(FLDIM),SCAL,SCALK,PCO
      DOUBLE PRECISION Z(FLDIM),SL(FLDIM),BM(FLDIM),GR(FLDIM),R(FLDIM)
      DOUBLE PRECISION SINDIP(FLDIM)

      !...... Validate data and set defaults for setting up a new field line
      CALL VALID_DATA(FLDIM,BLON,JMAX,PCO,Z0,SCAL,IDAY,RE,PLAT,PLON)

      !.. Set up the grid points in the orthogonal q coordinates.
      CALL QFIELD(JMAX,PCO,RE,Z0,SCAL,SCALK,Q)

      !... set up the magnetic field line spatial grid .....
      CALL RFIELD(JMAX,BLON,PLAT,PLON,PCO,RE,Z0,GL,SCAL,SCALK,
     >  Q,COSDIP,SINDIP,Z,SL,GLATD,GLOND,BM,GR,R,VTOT)

      RETURN
      END
C::::::::::::::::::::::::::::::: VALID_DATA :::::::::::::::::::::::::::::::::::
C..... This validates the data for setting up a new field line and establishes  
C..... default values. It also sets up the magnetic pole location
      SUBROUTINE VALID_DATA(FLDIM,BLON,JMAX,PCO,Z0,SCAL,IDAY,RE,
     >   PLAT,PLON)
      IMPLICIT NONE
      INTEGER IDAY,JMAX,FLDIM
      DOUBLE PRECISION PCO,Z0,SCAL
      DOUBLE PRECISION GLOND,GLATD,GLONP,DLAT,DLON,BLON,BLAT
      DOUBLE PRECISION RE,PLAT,PLON

      RE=6370.0  !.. Earth radius

      !.. Check day of year is good
      IF(IDAY.LT.1900000.OR.IDAY.GT.2100000) THEN
        WRITE(6,*) ' ** ERROR: Day # YYYYddd out of range = ',IDAY
        STOP
      ENDIF

      !.. If PCO large or too small for L-shell assume geographic latitude
      !.. and longitude and convert to magnetic coordinates using GTOB
!dbg20110817
 !     IF(PCO.LT.1.01.OR.PCO.GT.9) THEN
 !        GLATD=PCO
 !        GLOND=BLON
 !        !.. locate magnetic pole
 !        CALL MPOLE_LOC(GLOND,GLATD,IDAY,PLAT,PLON) 
 !        CALL GTOB(BLON,PLAT,PLON,BLAT,GLATD,GLOND)
 !        PCO=((RE+250.)/RE)/(COS(BLAT/57.296))**2
 !     ELSE
         GLONP=BLON-74.5
         IF(GLONP.LT.0) GLONP=GLONP+360
         CALL MPOLE_LOC(GLONP,50.0D0,IDAY,PLAT,PLON)  !.. locate magnetic pole
      print *,'GLONP=',GLONP,PLAT,PLON

 !     ENDIF

      !-- Apply the default value if JMAX has an inappropriate value
      IF(JMAX.LT.55*PCO+125) JMAX=55*PCO+125
      IF(JMAX.LT.301) JMAX=301
      IF(JMAX.GE.FLDIM) JMAX=FLDIM    !.. Array boundary check
      JMAX= 2*(JMAX/2)-1   !.. Make sure JMAX is an odd number
!dbg20110817
      JMAX=401 !dbg

      !.. Default scaling if unrealistic value
      IF(SCAL.LT.0.5.OR.SCAL.GT.11) SCAL=(8.0/PCO+2.0)

      RETURN
      END
C::::::::::::::::::::::::::::::: DATAWRITE ::::::::::::::::::::::::::
C.... this routine writes the field line data
      SUBROUTINE DATAWRITE(JMAX,IDAY,BLON,PCO,RE,GL,Q,COSDIP,SINDIP,
     >   Z,SL,GLATD,GLOND,BM,GR,R,VTOT)
 
      IMPLICIT NONE
      INTEGER J,JMAX,IDAY
      DOUBLE PRECISION BLON,RE,GLATD(JMAX),GLOND(JMAX)
      DOUBLE PRECISION COSDIP(JMAX),Q(JMAX),GL(JMAX)
      DOUBLE PRECISION BM(JMAX),GR(JMAX),R(JMAX),SL(JMAX)
      DOUBLE PRECISION Z(JMAX),SINDIP(JMAX),PCO,VTOT
      DOUBLE PRECISION TOTVOL,DELS,DELZ,VELEM

      WRITE(6,100) IDAY/1000,MOD(IDAY,1000),PCO,BLON,
     >   NINT(SL(JMAX)*1.0E-5),VTOT

 100  FORMAT(/4X,'  YEAR= ',I4,'    DAY= ',I3
     > ,//9X,' **** The field line grid parameters ****'
     > ,//,'  L-shell =',F7.3,';     Magnetic longitude= ',F5.1
     > ,//,' Length of flux tube = ',I8,' km'
     > ,'     Total flux tube volume = ',1P,E10.2,' cm-3'
     > ,//,'    pt     alt    arc_len   alt_dif   arc_dif   dip_ang'
     > ,2X,'cos_dip  mag_fld  Tube_vol  gravity   Re+alt'
     > ,2X,'B_lat  G_lat G_lon     Vol_elem    SinDip  Q-value')

      !.... Write out the field line parameters
      TOTVOL=0.0  !.. total flux tube volume to grid point
      DELS=0.0    !.. distance between 
      DELZ=0.0    !.. altitude difference of grid points
      DO J=1,JMAX
        IF(J.GT.1) THEN
          DELS=(SL(J)-SL(J-1))/1.0E5
          DELZ=ABS(Z(J)-Z(J-1))
          TOTVOL=TOTVOL+BM(1)*(SL(J)-SL(J-1))/BM(J)
        ENDIF
        VELEM=0.0
        IF(J.GT.1) VELEM=BM(1)*(SL(J)-SL(J-1))/BM(J)
        WRITE(6,110) J, Z(J), (SL(J)*1.0E-5), DELZ ,DELS,
     >   (ACOS(COSDIP(J))*57.3),COSDIP(J),BM(J),TOTVOL,GR(J),
     >   (R(J)*1.0E-5),GL(J)*57.296,GLATD(J),GLOND(J),VELEM
     >   ,SINDIP(J),Q(J)
      ENDDO

 110  FORMAT(I5,5F10.2,2F9.2,1P,E10.2,0P,2F9.1,1P,0P,3F7.1,
     >  1P,E15.6,0P,2F8.3)
      RETURN
      END
C::::::::::::::::::::::::::::::: RFIELD :::::::::::::::::::::::::::::::::
C.... this program is responsible for calling the subroutine that sets up
C.... the dipole magnetic field grid and filling the field parameter arrays,
C.... altitude, gravity etc.
C.... Reference: G.J. Bailey, Planet. Space Sci., Vol. 31, No. 4,
C....            pp. 389-409, 1983.
      SUBROUTINE RFIELD(JMAX,BLON,PLAT,PLON,PCO,RE,Z0,GL,SCAL,
     > SCALK,Q,COSDIP,SINDIP,Z,SL,GLATD,GLOND,BM,GR,R,VTOT) 
      IMPLICIT NONE
      INTEGER J,JU,JMAX,JQ
      DOUBLE PRECISION BLON,PLAT,PLON,RE,GLATD(JMAX),GLOND(JMAX)
      DOUBLE PRECISION FD(9),COSDIP(JMAX), X(JMAX),Q(JMAX),GL(JMAX)
      DOUBLE PRECISION BM(JMAX),GR(JMAX),R(JMAX),SL(JMAX)
      DOUBLE PRECISION Z(JMAX),SINDIP(JMAX),Z0,PCO,SCAL,SCALK
      DOUBLE PRECISION SREF,VTOT

      JQ=(JMAX+1)/2    !.. Equatorial grid point

      DO J=1,JQ
        CALL FIELD(J,JMAX,PCO,RE,SCAL,SCALK,Q(J), X(J),FD)
        IF(J.EQ.1) SREF=FD(8)*1.E5         ! distance along magnetic field
        Z(J)=FD(1)                         ! altitude in km
        GL(J)=FD(4)                        ! latitude
        BM(J)=FD(6)                        ! magnetic field strength
        COSDIP(J)=FD(3)                    ! cos(dip)
        SINDIP(J)=FD(2)                    ! sin(dip)
        GR(J)=-FD(5)                       ! gravity
        R(J)=(RE+FD(1))*1.E5               ! radial distance to grid point
        IF(J.NE.1) SL(J)=SREF-FD(8)*1.E5   ! distance to grid point along B
        !-----  Set up field line parameters in Southern hemisphere
        JU=JMAX+1-J              !.. symetric point of J about equator plane
        GL(JU)=-GL(J)
        Z(JU)=Z(J)
        BM(JU)=BM(J)
        COSDIP(JU)=-COSDIP(J)
        SINDIP(JU)=-SINDIP(J)
        GR(JU)=-GR(J)
        SL(JU)=SREF+FD(8)*1.E5
        R(JU)=R(J)
      ENDDO

      IF(Z(JQ).LT.80) THEN
        WRITE(6,696) NINT(Z(JQ)) 
 696    FORMAT(/'   ** ERROR: Field line apex altitude is',I3,' km.'
     >      ,1X,'It does not reach the ionosphere ***'/)
        STOP 
      ENDIF

      !.. Calculate total flux tube volume and geographic coordinates
      VTOT=0.0
      DO J=1,JMAX
        IF(J.GT.1) VTOT=VTOT+BM(1)*(SL(J)-SL(J-1))/BM(J)
        CALL BTOG(BLON,PLAT,PLON,GL(J)*57.296,GLATD(J),GLOND(J))
      ENDDO

      RETURN
      END
C:::::::::::::::::::::::::: QFIELD ::::::::::::::::::::::::::::::::::::::::
C.... This program determines Q values.
C.... JMAX = # of grid points, PCO = lshell, Z0 = lower boundary,
C.... scal = scaling factor to ensure that the x-coord ranged from 0 to 1.
C.... SCALK & Q are returned.
C.... Reference: G.J. Bailey, Planet. Space Sci., Vol. 31, No. 4,
C               pp. 389-409, 1983.
C
        SUBROUTINE QFIELD(JMAX,PCO,RE,Z0,SCAL,SCALK,Q)
      IMPLICIT NONE
      INTEGER JMAX,J
      DOUBLE PRECISION Q(JMAX),PCO,Z0,SCAL,SCALK,R0,RAT,QMAX,DX,SCX,NPTS
      DOUBLE PRECISION RE
 
      R0=RE*PCO        !.. equatorial radius to flux tube.
      NPTS=(JMAX+1)/2  !.. # points on half field line
      RAT=(RE+Z0)/RE
      QMAX=(SQRT(1-RAT/PCO))/(RAT**2)
      SCALK=1/SINH(SCAL*(QMAX))
      DO J=1,JMAX
        DX=1-(J-1)/(NPTS-1)
        !..  Q=the dipole coord  determined from -- sinh(kq)=dx
        SCX=DX/SCALK
        Q(J)=DLOG(SCX+DSQRT(SCX**2+1))/SCAL
      ENDDO
 
      RETURN
      END
 
C:::::::::::::::::::::::::: FIELD ::::::::::::::::::::::::::::::::::::::::
C  this program determines the grid point spacing given just JMAX=
C   # of grid points, PCO=lshell, Z0=lower boundary, SCAL=scaling
C   factor;; , the field parameters are transferred through FD(i)
C  Reference: G.J. Bailey, Planet. Space Sci., Vol. 31, No. 4,
C               pp. 389-409, 1983.
C--------------------------------------------------
      SUBROUTINE FIELD(J,JMAX,PCO,RE,SCAL,SCALK,Q, XR,FD)
      IMPLICIT NONE
      INTEGER J,JMAX,JITER
      DOUBLE PRECISION RE
      DOUBLE PRECISION Q,PCO,SCAL,SCALK
      DOUBLE PRECISION X,XX,XR,R,R0,SHX,CHX,FD(9),SQTH,FEX,DEX

      X=.5        !.. colatitude
      R0=RE*PCO   !.. equatorial radius to flux tube

      !... Values for equatorial point
      IF(J.EQ.(JMAX+1)/2) THEN
         X=1.570796327
         SHX=1
         CHX=0
      ELSE          !... use Newton solver
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
C::::::::::::::::::::::: MAGDEC :::::::::::::::::::::::::::
C...... This routine uses a table to calculate the magnetic declination
C...... (DECPT) for use with a neutral wind model. The table was constructed 
C...... using " The Earth's magnetic Field" by Robert Merrill and
C...... Michael McElhinny, Academic Press p 18. (UAH library # QC816.M47) 
C...... The declination at the specified location (degrees) is obtained
C...... by bilinear interpolation of the table values
C...... Written by Phil Richards, August 1995
C...... This code is not accurate near the magnetic poles
C...... GLAT, GLONG, and DECPT are in degrees
       SUBROUTINE MAGDEC(GLAT,GLONG,DECPT)
       IMPLICIT NONE
      !..--- NLONG=# longs, NLAT= # lats. I,J,K are array indices
       INTEGER NLONG,NLAT,I,J,K
       PARAMETER (NLONG=13,NLAT=11) 
       DOUBLE PRECISION GLAT,GLONG,DECPT
       DOUBLE PRECISION ALONG(NLONG),ALAT(NLAT),DECLS(13,11),RLAT,RLONG
      !... Fill up the longitude, latitude, and magnetic declination arrays
       DATA ALONG/0,30,60,90,120,150,180,210,240,270,300,330,360/
       DATA ALAT/-90,-80,-60,-40,-20,0,20,40,60,80,90/
      !..-- The declinations are loaded in latitude slices. The first NLONG
      !..-- values are for the first latitude in the ALAT array
       DATA DECLS / -20,-45,-67,-95,-125,170,110,80,60,40,20,0,-20
     > ,-20,-40,-63,-82,-120,140,82,65,52,35,16,-2,-20, -20,-35,-55,-70
     > ,-60,35,45,40,38,30,10,-7.5,-20, -27,-28,-42,-35,-5,13,20,21,21
     > ,20,0,-21,-27, -21,-15,-22,-10,2,8,13,13,12,11,-7,-25,-21, -10
     > ,-2,-4,-3.5,2,6,12,9,8,6,-10,-20,-10, -5,1,0,-2,-2,1,10,12,11,5
     > ,-14,-16,-5, -5,3,5,2,-7,-5,7,17,17,2,-20,-16,-5, -7,7,16,10,-12
     > ,-11,6,25,30,-7.5,-36,-15,-7, -10,10,25,20,-11,-11,7,30,42,-40
     > ,-52,-33,-10, -12,10,27,25,-7,-7,7,32,50,-70,-60,-35,-12/
      !.. Find the index in the longitude array
       DO 20 I=1,NLONG
         IF(GLONG.LT.ALONG(I)) GO TO 30
  20   CONTINUE
  30   CONTINUE
       J=I-1
      !.. Find the index in the latitude array
       DO 40 I=1,NLAT
       print *,'MAGDEC:i',i,'nlat'  ,nlat,ALAT(I),GLAT,GLONG  
       IF(GLAT.LT.ALAT(I)) GO TO 50
 40    CONTINUE
 50    CONTINUE
       K=I-1
      !..-- Now do the bilinear fit (Press et al. Num. Rec. page 96)
       RLAT=(GLAT-ALAT(K))/(ALAT(K+1)-ALAT(K))
       RLONG=(GLONG-ALONG(J))/(ALONG(J+1)-ALONG(J))
       DECPT=(1-RLAT)*(1-RLONG)*DECLS(J,K)+RLONG*(1-RLAT)*DECLS(J+1,K)
     >   +RLAT*RLONG*DECLS(J+1,K+1)+(1-RLONG)*RLAT*DECLS(J,K+1)
       END
C:::::::::::::::::::::::: MPOLE_LOC ::::::::::::::::::::::::::::::::
C... This function calls the routines that determine the latitude and 
C... longitude of the pseudo north magnetic pole. These are different
C... for the two hemispheres
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_LOC(GLONG,     !... Geographic long. in degrees
     >                     GLAT,      !... Geographic lat. in degrees
     >                     IDAY,      !... year
     >                     PLAT,PLON) !... pole position in radians
      IMPLICIT NONE
      INTEGER IDAY
      DOUBLE PRECISION GLAT,GLONG,PLAT,PLON
      IF(GLAT.GE.0) CALL MPOLE_NTH(GLONG,IDAY,PLAT,PLON)
      IF(GLAT.LT.0) CALL MPOLE_STH(GLONG,IDAY,PLAT,PLON)
      RETURN
      END
C:::::::::::::::::::::::: MPOLE_NTH ::::::::::::::::::::::::::::::::
C... This function returns the latitude and longitude of the
C... pseudo north magnetic pole for a tilted, centred, dipole. The 
C... pole position varies with geographic location on the Earth so 
C... as to provide fair agreement with L values from the IGRF model.
C... This routine is for the northern latitudes. There is a separate
C... routine for southern latitudes. It is most accurate for L=3 and 
C... not very accurate at the Atlantic longitudes.
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_NTH(GLONG,      !... Geographic long. in degrees
     >                     IDAY,      !... year
     >                     PLAT,PLON)  !... pole position in radians
      IMPLICIT NONE
      INTEGER I,DIM,IDAY
      PARAMETER (DIM=19)
      DOUBLE PRECISION GLONG,PLAT,PLON,DLAT,DLON,GEOLON(DIM),MAGLON(DIM)
      DOUBLE PRECISION MAGLAT(DIM),GRAD,XCUT

      DATA GEOLON/0,20,40,60,80,100,120,140,160,180,200,220,240
     > ,260,280,300,320,340,360/

      DATA MAGLAT/78.0,78.0,80.0,81.0,83.0,83.5,83.5,82.0,81.0
     > ,80.0,80.0,79.0,78.0,77.0,77.5,79.0,80.0,79.0,78.0/

      DATA MAGLON/270.0,280.0,284.0,295.0,300.0,295.0,295.0,295.0,295.0
     > ,295.0,295.0,293.0,290.0,290.0,290.0,280.0,267.0,265.0,270.0/


      IF(GLONG.LT.0) GLONG=GLONG+360
      IF(GLONG.GT.360) GLONG=GLONG-360
      PLAT=MAGLAT(1)
      PLON=MAGLON(1)

      !..Loop backwards to locate the GLONG within the GLONG array
      I=DIM+1
  10  I=I-1
      IF(GLONG.LT.GEOLON(I)) GO TO 10

      !... Now find the longitude by linear interpolation
      GRAD=(MAGLON(I+1)-MAGLON(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLON(I)-GEOLON(I)*GRAD
      PLON=GLONG*GRAD+XCUT

      !... now find the latitude by linear interpolation
      GRAD=(MAGLAT(I+1)-MAGLAT(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLAT(I)-GEOLON(I)*GRAD
      PLAT=GLONG*GRAD+XCUT

      CALL MPOLE_INC(IDAY,DLAT,DLON)   !..find change since 1980

      PLON=(DLON+PLON)/57.2958   !.. convert to radians
      PLAT=(DLAT+PLAT)/57.2958

      RETURN
      END
C:::::::::::::::::::::::: MPOLE_STH ::::::::::::::::::::::::::::::::
C... This function returns the latitude and longitude of the
C... pseudo north magnetic pole for a tilted, centred, dipole. The 
C... pole position varies with geographic location on the Earth so 
C... as to provide fair agreement with L values from the IGRF model.
C... This routine is for the southern latitudes. There is a separate
C... routine for northern latitudes. It is most accurate for L=3 and 
C... not very accurate at the Atlantic longitudes.
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_STH(GLONG,      !... Geographic long. in degrees
     >                     IDAY,      !... year
     >                     PLAT,PLON)  !... pole position in radians
      IMPLICIT NONE
      INTEGER I,DIM,IDAY
      PARAMETER (DIM=19)
      DOUBLE PRECISION GLONG,PLAT,PLON,DLAT,DLON,GEOLON(DIM),MAGLON(DIM)
      DOUBLE PRECISION MAGLAT(DIM),GRAD,XCUT

      DATA GEOLON/0,20,40,60,80,100,120,140,160,180,200,220,240
     > ,260,280,300,320,340,360/

      DATA MAGLON/285.0,280.0,280.0,280.0,280.0,290.0,300.0,310.0,310.0
     > ,305.0,300.0,297.0,290.0,290.0,290.0,285.0,285.0,285.0,285.0/

      DATA MAGLAT/78.0,78.0,79.0,79.0,78.50,77.0,77.0,78.0,79.0,80.0
     > ,80.0,79.0,78.0,76.0,75.0,75.0,78.0,78.0,78.0/


      IF(GLONG.LT.0) GLONG=GLONG+360
      IF(GLONG.GT.360) GLONG=GLONG-360
      PLAT=MAGLAT(1)
      PLON=MAGLON(1)

      !.. Loop backwards to locate the GLONG within the GLONG array
      I=DIM+1
  10  I=I-1
      IF(GLONG.LT.GEOLON(I)) GO TO 10

      !... Now find the longitude by linear interpolation
      GRAD=(MAGLON(I+1)-MAGLON(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLON(I)-GEOLON(I)*GRAD
      PLON=GLONG*GRAD+XCUT

      !... now find the latitude by linear interpolation
      GRAD=(MAGLAT(I+1)-MAGLAT(I))/(GEOLON(I+1)-GEOLON(I))
      XCUT=MAGLAT(I)-GEOLON(I)*GRAD
      PLAT=GLONG*GRAD+XCUT

      !... Year correction not needed in South
      PLON=PLON/57.2958   !.. convert to radians
      PLAT=PLAT/57.2958

      RETURN
      END
C:::::::::::::::::::: MPOLE_INC ::::::::::::::::::::::::::::::
C... This function finds returns the change in latitude and longitude 
C... of the north magnetic pole referred to its location in 1980.
C... This is from the NSSDC web site
C... Written by P. Richards in June 2001
      SUBROUTINE MPOLE_INC(IDAY,DLAT,DLON)
      IMPLICIT NONE
      INTEGER I,DIM,IDAY
      PARAMETER (DIM=12)
      DOUBLE PRECISION YEAR,PLAT,PLON,MAGYR(DIM),MAGLON(DIM),MAGLAT(DIM)
      DOUBLE PRECISION DLAT,DLON,GRAD,XCUT
      !... Years that data are location specified
      DATA MAGYR/1950,1955,1960,1965,1970,1975,1980,1985,1990,1995,
     >           2000,2005/
      !... Magnetic longitude of magnetic pole
      DATA MAGLON/281.05,280.87,280.66,280.41,280.14,279.84,279.44
     >            ,279.11,278.76,278.35,278.06,277.89/
      !... Magnetic latitude of magnetic pole
      DATA MAGLAT/79.85,79.95,80.05,80.13,80.22,80.34,80.55,80.80
     >             ,81.04,81.28,81.64,81.99/

      YEAR=IDAY/1000
      !... make sure the year is within the legal range
      IF(YEAR.LT.MAGYR(1)) THEN
        WRITE(6,90) NINT(YEAR),NINT(MAGYR(1)),NINT(MAGYR(1))
        PLAT=MAGLAT(1)
        PLON=MAGLON(1)
        RETURN
      ELSE IF(YEAR.GE.MAGYR(DIM)) THEN
        WRITE(6,91) NINT(YEAR),NINT(MAGYR(DIM)),NINT(MAGYR(DIM))
        PLAT=MAGLAT(DIM)
        PLON=MAGLON(DIM)
        RETURN
      ENDIF

 90   FORMAT(1X,' *** Requested year',I5,' is less than',I5,
     >  '; the magnetic pole location has been set to the',I5,
     >  ' value')
 91   FORMAT(1X,' *** Requested year',I5,' is greater than',I5,
     >  '; the magnetic pole location has been set to the',I5,
     >  ' value')

      !.. Loop backwards to locate the year within the year array
      I=DIM+1
  10  I=I-1
      IF(YEAR.LT.MAGYR(I)) GO TO 10

      !... Now find the longitude by linear interpolation
      GRAD=(MAGLON(I+1)-MAGLON(I))/(MAGYR(I+1)-MAGYR(I))
      XCUT=MAGLON(I)-MAGYR(I)*GRAD
      PLON=YEAR*GRAD+XCUT

      !... now find the latitude by linear interpolation
      GRAD=(MAGLAT(I+1)-MAGLAT(I))/(MAGYR(I+1)-MAGYR(I))
      XCUT=MAGLAT(I)-MAGYR(I)*GRAD
      PLAT=YEAR*GRAD+XCUT

      !.. determine the change in pole position since 1980 in radians
      DLON=(PLON-MAGLON(7))
      DLAT=(PLAT-MAGLAT(7))
      RETURN
      END
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE BTOG(BLOND,PLAT,PLON,BLATD,GLATD,GLOND)
C....... Conversion from centered dipole to geographic coordinates. (PLAT
C....... ,PLON)= coords of north magnetic pole. The geomagnetic longitude
C....... is measured east from the geomagnetic prime meridian at 291 degrees
C....... east geographic.
      IMPLICIT NONE
      DOUBLE PRECISION BLOND,BLATD,GLATD,GLOND
      DOUBLE PRECISION BLONR,BLATR,XM,YM,ZM,XG,YG,ZG
      DOUBLE PRECISION RE,PLAT,PLON,PI,ANGVEL

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
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GTOB(BLOND,PLAT,PLON,BLATD,GLATD,GLOND)
C........ Conversion from geographic to centered dipole coordinates
      IMPLICIT NONE
      DOUBLE PRECISION BLOND,BLATD,GLATD,GLOND
      DOUBLE PRECISION XM,YM,ZM,XG,YG,ZG,GLAT,GLON
      DOUBLE PRECISION RE,PLAT,PLON,PI,ANGVEL

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
