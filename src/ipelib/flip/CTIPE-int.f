C............................ CTIPE-int.for .........................
C.. This file contains the basic interface routines for the CTIP interface
C::::::::::::::::::::::::::: Module FIELD_LINE_GRID :::::::::::::::::
C.... Contains parameters associated with the field line grid
C.... Written by P. Richards June 2010.
      MODULE FIELD_LINE_GRID
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER JMIN,JMAX   !.. first and last indices on field line grid
        INTEGER, PARAMETER :: FLDIM = 1115     !.. Field line grid dimension
        DOUBLE PRECISION Z(FLDIM),SZA(FLDIM)  !.. altitude, Solar zenith angle
        DOUBLE PRECISION BM(FLDIM)            !.. magnetic field strength
        DOUBLE PRECISION SL(FLDIM)            !.. distance along the field line
        DOUBLE PRECISION GR(FLDIM),GL(FLDIM)  !.. Gravity, magnetic latitude
      END MODULE FIELD_LINE_GRID
C::::::::::::::::::::::::::: ION_DEN_VEL :::::::::::::::::
C.... Ion densities and velocities
C.... Written by P. Richards June 2010.
C.... O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      MODULE ION_DEN_VEL
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: IDIM = 1115  !.. Field line grid dimension
        INTEGER, PARAMETER :: ISPEC = 9   !.. Species dimension
      !.. Ion densities & Velocities
      DOUBLE PRECISION XIONN(ISPEC,IDIM),XIONV(ISPEC,IDIM)
      END MODULE ION_DEN_VEL
C:::::::::::::::::::::::::::Module SOLVAR :::::::::::::::::
C.... This module is used to conserve space used by the Newton solver.
C.... The arrays are workspace for the solver
C.... Written by P. Richards June 2010.
      MODULE SOLVARR
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: SDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION DELTA(3*SDIM),RHS(3*SDIM)
        DOUBLE PRECISION WORK(15*SDIM),S(15*SDIM)
      END MODULE SOLVARR
C::::::::::::::::::::::::::: Module THERMOSPHERE :::::::::::::::::
C.... Contains the neutral densities and temperatures, O+ production, electron 
C.... heating
C.... Written by P. Richards June 2010.
      MODULE THERMOSPHERE
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: TDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION ON(TDIM),HN(TDIM),N2N(TDIM),O2N(TDIM),HE(TDIM)
        DOUBLE PRECISION TN(TDIM),UN(TDIM),TINF(TDIM)
        DOUBLE PRECISION EHT(3,TDIM)
        DOUBLE PRECISION COLFAC  !.. O+ - O collision frequency Burnside factor 1-1.7
      END MODULE THERMOSPHERE
C::::::::::::::::::::::::::: AVE_PARAMS :::::::::::::::::
C.... Computes midpoint values associated with the field line grid
C.... Written by P. Richards June 2010.
      MODULE AVE_PARAMS
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: AVDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION TEJ(AVDIM),TIJ(AVDIM),NUX(2,AVDIM), RBM(AVDIM)
        DOUBLE PRECISION UNJ(AVDIM),DS(AVDIM)
        DOUBLE PRECISION GRADTE(AVDIM),GRADTI(AVDIM),GRAV(AVDIM)
        DOUBLE PRECISION OLOSS(AVDIM),HLOSS(AVDIM),HPROD(AVDIM)
      END MODULE AVE_PARAMS
C::::::::::::::::::::::::::: PRODUCTION :::::::::::::::::
C.... Holds the EUV, photoelectron, and auroral production
C.... Written by P. Richards August 2010.
      !.. Total ionization and excitation rates (EUV + PE + Aurora)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      MODULE PRODUCTION
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: PDIM = 1115  !.. Field line grid dimension
      !.. EUV and photoelectron production rates 
      REAL EUVION(3,12,PDIM),PEXCIT(3,12,PDIM),PEPION(3,12,PDIM)
      DOUBLE PRECISION OTHPR1(6,PDIM),OTHPR2(6,PDIM)
      !.. Total ionization and excitation rates (EUV + PE + Aurora)
      REAL SUMION(3,12,PDIM),SUMEXC(3,12,PDIM)
      !.. Auroral ionization and excitation rates
      REAL PAUION(3,12,PDIM),PAUEXC(3,12,PDIM)
      DOUBLE PRECISION NPLSPRD(PDIM)  !.. N+ production
      DOUBLE PRECISION PHION(PDIM)    !.. Total O+ production
      END MODULE PRODUCTION
C::::::::::::::::::::::::::: Module MINORNEUT :::::::::::::::::
C.... Contains the minor neutral densities 
C.... Written by P. Richards September 2010.
      !.. USE MINORNEUT N4S N2D NNO N2P N2A O1D O1S
      MODULE MINORNEUT
        IMPLICIT NONE
        !... Check dimensions are the same in all FLIP modules
        INTEGER, PARAMETER :: NDIM = 1115  !.. Field line grid dimension
        DOUBLE PRECISION N4S(NDIM),N2D(NDIM),NNO(NDIM),N2P(NDIM),
     >    N2A(NDIM),O1D(NDIM),O1S(NDIM),EQN2D(NDIM)
      END MODULE MINORNEUT
C::::::::::::::::::::::::::::::::: PHOTOIONIZATION_DATA :::::::::::::::
C..... This module contains data for the PRIMPR photoionization routine
C..... Written by P. Richards in August 2010
      MODULE PHOTOIONIZATION_DATA  !.. NPOT LMAX PROB ZLAM SIGION SIGABS TPOT 
      IMPLICIT NONE
      INTEGER, PARAMETER :: LMAX =37  !.. # of EUV wavelengths
      INTEGER NPOT(3)                 !.. # of ionization potentials
      REAL PROB(3,6,LMAX)             !.. Ionization probability,
      REAL ZLAM(LMAX)                 !.. EUV wavelengths,
      REAL SIGION(3,LMAX)             !.. Ionization cross sections
      REAL SIGABS(3,LMAX)             !.. Absorption cross sections
      REAL TPOT(3,10)                 !.. Ionization potentials

      DATA NPOT/5,5,6/   !.. # of O, O2, N2 ionization potentials 
      !.. Ionization potentials
      DATA TPOT/13.6,12.1,15.6,16.9,16.5,16.7,18.6,18.2,18.8,
     >   28.5,20.0,25.0,40.0,25.0,29.0,0.0,0.0,37.0, 12*0.0/
      !.... wavelength data. average is taken for bands
      DATA ZLAM/1025.,1031.91,1025.72,975.,977.02,925.,875.,825.,775.,
     > 789.36,770.41,765.15,725.,703.36,675.,625.,629.73,609.76,575.,
     > 584.33,554.31,525.,475.,465.22,425.,375.,368.07,325.,303.78,
     > 303.31,275.,284.15,256.3,225.,175.,125.,75./

      DATA SIGION/0.0,2.59000E-19,0.00000E+00,0.00000E+00,0.0,
     > 0.00000E+00,0.00000E+00,1.05000E-18,0.00000E+00,0.00000E+00,
     > 1.39400E-17,0.00000E+00,0.00000E+00,1.55400E-17,0.00000E+00,
     > 1.31500E-18,9.37400E-18,0.00000E+00,4.55400E-18,5.49400E-18,
     > 0.00000E+00,3.49800E-18,6.41300E-18,0.00000E+00,5.09100E-18,
     > 1.05970E-17,1.42740E-17,3.74900E-18,1.01910E-17,8.86000E-18,
     > 3.89000E-18,8.47000E-18,8.50000E-18,4.00000E-18,1.17200E-17,
     > 6.58000E-17,1.07360E-17,2.38050E-17,1.50600E-17,1.14600E-17,
     > 2.37500E-17,2.54800E-17,1.72450E-17,2.13060E-17,2.92350E-17,
     > 1.33650E-17,2.49370E-17,2.33390E-17,1.34000E-17,3.11000E-17,
     > 2.33700E-17,1.34000E-17,2.63900E-17,2.27900E-17,1.30240E-17,
     > 2.66100E-17,2.27870E-17,1.30900E-17,2.27900E-17,2.24000E-17,
     > 1.25900E-17,2.60400E-17,2.41300E-17,1.20590E-17,2.46060E-17,
     > 2.45010E-17,1.21270E-17,2.31010E-17,2.34710E-17,1.19300E-17,
     > 2.19100E-17,2.31600E-17,1.14960E-17,2.03100E-17,2.16750E-17,
     > 9.68700E-18,1.81180E-17,1.63950E-17,9.84000E-18,1.83200E-17,
     > 1.69100E-17,8.69300E-18,1.74380E-17,1.38570E-17,7.70000E-18,
     > 1.68100E-17,1.17000E-17,7.68000E-18,1.68000E-17,1.16700E-17,
     > 6.46100E-18,1.43870E-17,1.04930E-17,7.08000E-18,1.57900E-17,
     > 1.09000E-17,6.05000E-18,1.33700E-17,1.02100E-17,5.20200E-18,
     > 1.09000E-17,8.39200E-18,3.73200E-18,7.50900E-18,4.95800E-18,
     > 1.83900E-18,3.80600E-18,2.26100E-18,7.30000E-19,1.31600E-18,
     > 7.20000E-19/
      DATA SIGABS/0.0000,1.34600E-18,0.000,0.00000E+00,1.00000E-18,
     > 0.00000E+00,0.00000E+00,1.63000E-18,0.00000E+00,0.00000E+00,
     > 2.11080E-17,5.09880E-17,0.00000E+00,1.87300E-17,2.24000E-18,
     > 1.31500E-18,1.28170E-17,9.68000E-18,4.55400E-18,8.56200E-18,
     > 2.02490E-17,3.49800E-18,1.66310E-17,1.69920E-17,5.09100E-18,
     > 2.21450E-17,3.35780E-17,3.74900E-18,2.66680E-17,1.64870E-17,
     > 3.89000E-18,1.89100E-17,1.41800E-17,4.00000E-18,2.08000E-17,
     > 1.20490E-16,1.07360E-17,2.85350E-17,2.46620E-17,1.14600E-17,
     > 2.74400E-17,2.65400E-17,1.72450E-17,2.19190E-17,3.17550E-17,
     > 1.33650E-17,2.60170E-17,2.33390E-17,1.34000E-17,3.20600E-17,
     > 2.33700E-17,1.34000E-17,2.80700E-17,2.27900E-17,1.30240E-17,
     > 2.66100E-17,2.27870E-17,1.30900E-17,2.27900E-17,2.24000E-17,
     > 1.25900E-17,2.60400E-17,2.41300E-17,1.20590E-17,2.46060E-17,
     > 2.45010E-17,1.21270E-17,2.31010E-17,2.34710E-17,1.19300E-17,
     > 2.19100E-17,2.31600E-17,1.14960E-17,2.03100E-17,2.16750E-17,
     > 9.68700E-18,1.81180E-17,1.63950E-17,9.84000E-18,1.83200E-17,
     > 1.69100E-17,8.69300E-18,1.74380E-17,1.38570E-17,7.70000E-18,
     > 1.68100E-17,1.17000E-17,7.68000E-18,1.68000E-17,1.16700E-17,
     > 6.46100E-18,1.43870E-17,1.04930E-17,7.08000E-18,1.57900E-17,
     > 1.09000E-17,6.05000E-18,1.33700E-17,1.02100E-17,5.20200E-18,
     > 1.09000E-17,8.39200E-18,3.73200E-18,7.50900E-18,4.95800E-18,
     > 1.83900E-18,3.80600E-18,2.26100E-18,7.30000E-19,1.31600E-18,
     > 7.20000E-19/
       DATA PROB/0.0,1.0,1.0,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,1.000000,1.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,0.957794,1.000000,0.000000,0.042206,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 1.000000,0.895954,1.000000,0.000000,0.104046,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.627701,0.565000,0.527407,0.372299,0.435000,0.472593,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.426702,0.668824,0.388933,0.573298,0.331176,0.611067,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.297362,0.395273,0.335986,0.664135,0.463591,0.664014,
     > 0.038504,0.141136,0.000000,0.000000,0.000000,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.317995,0.240422,0.308000,0.460831,0.365353,0.589000,
     > 0.221175,0.340922,0.103000,0.000000,0.053304,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.289552,0.242130,0.308000,0.469403,0.359162,0.589000,
     > 0.241045,0.352241,0.103000,0.000000,0.046467,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.269403,0.234886,0.308000,0.460448,0.383177,0.589000,
     > 0.270149,0.306381,0.103000,0.000000,0.075556,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.277948,0.354208,0.332000,0.453778,0.280858,0.568250,
     > 0.268274,0.214870,0.099750,0.000000,0.150064,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.299465,0.303890,0.323043,0.459893,0.328831,0.575994,
     > 0.240642,0.213666,0.100963,0.000000,0.153613,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.289913,0.355053,0.351862,0.450357,0.235358,0.551077,
     > 0.259730,0.216569,0.097060,0.000000,0.193020,0.000000,
     > 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
     > 0.283689,0.359411,0.380000,0.450037,0.282907,0.526750,
     > 0.266274,0.128627,0.093250,0.000000,0.219594,0.000000,
     > 0.000000,0.009461,0.000000,0.000000,0.000000,0.000000,
     > 0.279871,0.265374,0.405764,0.450153,0.160558,0.469567,
     > 0.269976,0.077351,0.082915,0.000000,0.084977,0.039289,
     > 0.000000,0.411741,0.002465,0.000000,0.000000,0.000000,
     > 0.279966,0.242785,0.414980,0.450126,0.139718,0.465742,
     > 0.269908,0.068087,0.081994,0.000000,0.064217,0.034222,
     > 0.000000,0.485192,0.003062,0.000000,0.000000,0.000000,
     > 0.267310,0.419160,0.439757,0.425887,0.235695,0.446402,
     > 0.259743,0.120885,0.078528,0.047060,0.122618,0.028560,
     > 0.000000,0.101641,0.006687,0.000000,0.000000,0.000066,
     > 0.260452,0.397023,0.361839,0.400846,0.223099,0.480085,
     > 0.250026,0.122412,0.099384,0.088676,0.144374,0.026681,
     > 0.000000,0.113092,0.030068,0.000000,0.000000,0.001944,
     > 0.260061,0.393955,0.346750,0.400000,0.221354,0.478991,
     > 0.259959,0.122624,0.101074,0.079980,0.147389,0.029676,
     > 0.000000,0.114679,0.040767,0.000000,0.000000,0.002741,
     > 0.259864,0.374885,0.267385,0.396411,0.210504,0.462571,
     > 0.249971,0.123939,0.107338,0.093754,0.166130,0.022626,
     > 0.000000,0.124542,0.114657,0.000000,0.000000,0.025423,
     > 0.250000,0.365000,0.256053,0.370000,0.205000,0.440509,
     > 0.250000,0.125000,0.103843,0.090000,0.055006,0.021138,
     > 0.040000,0.249993,0.089033,0.000000,0.000000,0.089423,
     > 0.250000,0.365000,0.255808,0.370052,0.205000,0.440035,
     > 0.250000,0.125000,0.103766,0.089974,0.055006,0.021076,
     > 0.039974,0.249993,0.088240,0.000000,0.000000,0.091075,
     > 0.251973,0.365000,0.235293,0.359851,0.205000,0.404410,
     > 0.230305,0.125000,0.095588,0.098592,0.055006,0.011581,
     > 0.059279,0.249993,0.069486,0.000000,0.000000,0.183642,
     > 0.250000,0.365000,0.242093,0.370056,0.205000,0.416098,
     > 0.230226,0.125000,0.098350,0.100282,0.055006,0.014699,
     > 0.049435,0.249993,0.078089,0.000000,0.000000,0.150670,
     > 0.254380,0.365000,0.221395,0.353388,0.205000,0.380522,
     > 0.225289,0.125000,0.089942,0.098843,0.055006,0.006278,
     > 0.068099,0.249993,0.037670,0.000000,0.000000,0.264193,
     > 0.257513,0.365000,0.207040,0.346101,0.205000,0.355850,
     > 0.224710,0.125000,0.084110,0.099784,0.055006,0.000000,
     > 0.071892,0.249993,0.000000,0.000000,0.000000,0.353000,
     > 0.270685,0.365000,0.204800,0.332954,0.205000,0.352000,
     > 0.218368,0.125000,0.083200,0.098948,0.055006,0.000000,
     > 0.079045,0.249993,0.000000,0.000000,0.000000,0.360000,
     > 0.294011,0.365000,0.204800,0.320024,0.205000,0.352000,
     > 0.208711,0.125000,0.083200,0.098609,0.055006,0.000000,
     > 0.078645,0.249993,0.000000,0.000000,0.000000,0.360000,
     > 0.296412,0.365000,0.204800,0.321373,0.205000,0.352000,
     > 0.209048,0.125000,0.083200,0.096724,0.055006,0.000000,
     > 0.076443,0.249993,0.000000,0.000000,0.000000,0.360000/ 
      END MODULE PHOTOIONIZATION_DATA
      

C:::::::::::::::::::::::::::::::::: CTIPINT :::::::::::::::::::::::::::::
C.... Interface for CTIP model. This routine uploads the CTIPe variables to 
C.... modules and calls the different FLIP routines.
C.... The parameters ending in X are dummy variables to be transferred to 
C.... the FLIP module variables. 
C.... Dummy variables are also used because it is assumed that the values 
C.... of FLIP variables are not preserved from call to call
C---- Additional mods Change TI to (temp_ti_te(2,DIM)
C---- Change SCOLUMN for grazing incidence
C.... Written by P. Richards June-September 2010.
      SUBROUTINE CTIPINT(
     >             JMINX,  !.. index of the first point on the field line
     >             JMAXX,  !.. index of the last point on the field line
     >           CTIPDIM,  !.. CTIPe array dimension, must equal to FLDIM
     >                ZX,  !.. array, altitude (km)
     >               PCO,  !.. p coordinate (L-shell)
     >               SLX,  !.. array, distance of point from northern hemisphere
     >               GLX,  !.. array, magnetic latitude (radians)
     >               BMX,  !.. array, magnetic field strength, (nT)
     >               GRX,  !.. array, gravity, cm2 s-1
     >                OX,  !.. array, O density (cm-3)
     >                HX,  !.. array, H density (cm-3)
     >               N2X,  !.. array, N2 density (cm-3)
     >               O2X,  !.. array, O2 density (cm-3)
     >               HEX,  !.. array, He density (cm-3)
     >              N4SX,  !.. array, N(4S) density (cm-3)
     >              INNO,  !.. switch
     >              NNOX,  !.. array, NO density (cm-3)
     >               TNX,  !.. array, Neutral temperature (K)
     >             TINFX,  !.. array, Exospheric Neutral temperature (K)
     >               UNX,  !.. array, Neutral wind (cm/s)
     >                DT,  !.. CTIPe time step (secs)
     >             DTMIN,  !.. Minimum time step allowed (>=10 secs?)
     >              F107,  !.. Daily F10.7
     >             F107A,  !.. 81 day average F10.7
     >              SZAX,  !.. Solar Zenith angle (radians)
     >              FPAS,  !.. Pitch angle scattering fraction
     >              HPEQ,  !.. Sets initial equatorial H+ density. See declaration below
     >            HEPRAT,  !.. Intial He+/H+ ratio (.01 to 1.0)
     >           COLFACX,  !.. O+ - O collision frequency Burnside factor (1.0 to 1.7)
     >            IHEPLS,  !..
     >             INPLS,  !..
     >             UTHRS,  !.. Universal time in hours !$$$
     >              EHTX,  !.. IN/OUT 2D array, Electron & ion heating rate (eV cm-3 s-1)
     >          AUR_PROD,  !.. IN 2D, array, ionization rates 
     >            TE_TIX,  !.. IN/OUT: 2D array, Electron and ion temperatures (K) (see below)
     >     XIONNX,XIONVX,  !.. IN/OUT: 2D array, Storage for ion densities and velocities 
     >             NHEAT,  !.. OUT: array, Neutral heating rate (eV/cm^3/s) 
     >             EFLAG,mp,lp,nflag_t,nflag_d)  !.. OUT: 2D array, Error Flags
      USE THERMOSPHERE       !.. ON HN N2N O2N HE TN UN EHT COLFAC
      USE MINORNEUT          !.. N4S N2D NNO N2P N2A O1D O1S
      USE FIELD_LINE_GRID    !.. FLDIM JMIN JMAX FLDIM Z BM GR SL GL SZA
      USE ION_DEN_VEL        !.. O+ H+ He+ N+ NO+ O2+ N2+ O+(2D) O+(2P)
      !..EUVION PEXCIT PEPION OTHPR1 OTHPR2 SUMION SUMEXC PAUION PAUEXC NPLSPRD
      USE PRODUCTION         !.. EUV, photoelectron, and auroral production
      USE SOLVARR
      USE AVE_PARAMS
      USE PRODUCTION

      IMPLICIT NONE
      integer mp,lp,i_which_call,nflag_t,nflag_d,test_te
      INTEGER IHEPLS, INPLS, INNO
      INTEGER CTIPDIM         !.. CTIPe array dimension, must equal to FLDIM
      INTEGER JTI             !.. Dummy variable to count the number of calls to this routine
      INTEGER I,J,JMINX,JMAXX !.. lcv + spatial grid indices
      INTEGER EFLAG(11,11)    !.. error flags, check =0 on return from FLIP
      INTEGER DEBUG           !.. switch to turn on debug writes 0=off, 1=on
      !.. sw_HEplus,sw_Nplus turn on diffusive solutions if > 0. no solution if 0, 
      !.. chemical equilibrium if < 0
      REAL F107,F107A         !.. daily and 81 day average F10.7      !.. 
      REAL EUVFLUX(37),UVFAC(59)  !.. Solar EUV fluxes and mult. factors
      DOUBLE PRECISION M_to_CM,M3_to_CM3,nT_TO_Gauss !.. Unit conversion factors
      !.. Dummy variables ending in X for transferring CTIP to FLIP modules, 
      !.. see above for documentation
      DOUBLE PRECISION N4SX(CTIPDIM),NNOX(CTIPDIM),ZX(CTIPDIM)
      DOUBLE PRECISION UNX(CTIPDIM),COLFACX,BMX(CTIPDIM)
      DOUBLE PRECISION GRX(CTIPDIM),SLX(CTIPDIM),GLX(CTIPDIM)
      DOUBLE PRECISION OX(CTIPDIM),HX(CTIPDIM),N2X(CTIPDIM),O2X(CTIPDIM)
      DOUBLE PRECISION HEX(CTIPDIM),TNX(CTIPDIM),SZAX(CTIPDIM)
      DOUBLE PRECISION XIONNX(9,CTIPDIM),XIONVX(9,CTIPDIM)
      !.. TE_ti(3,J) = Te, TE_TIX(2,J) = Ti = TE_TIX(2,J)
      DOUBLE PRECISION TE_TIX(3,CTIPDIM)
      !.. EHTX(3,J) = e heating rate, EHTX(1,J) = ion heating rate, EHTX(2,J) unused
      DOUBLE PRECISION EHTX(3,CTIPDIM)
      DOUBLE PRECISION AUR_PROD(3,CTIPDIM)
      !.. TINFX has to be an array for grazing incidence column densities
      DOUBLE PRECISION TINFX(CTIPDIM)  !.. Exospheric temperature
      !.. End dummy variable declarations 
      !.. HPEQ is equatorial H+ density = HPEQ * density of a full flux tube.
      !.. If 0.0 the densities from the previous time step are used.
      !.. If positive, initial densities and temperatures are set. 
      !.. If negative, H+ and He+ densities are reset for flux tube depletion
      DOUBLE PRECISION HPEQ              !.. HPEQ is equatorial H+ density
      DOUBLE PRECISION DT,DTMIN,FD(9),BCKPRD,FPAS,HEPRAT,PCO,UTHRS !$$$
      DOUBLE PRECISION COLUM(3,CTIPDIM)  !.. Neutral column densities for PEPRIM
      DOUBLE PRECISION N(4,CTIPDIM)      !.. FLIP variable for O+ H+ & total ions
      DOUBLE PRECISION temp_ti_te(3,CTIPDIM)     !.. FLIP variable for Te and Ti
      DOUBLE PRECISION N_SAVE(4,CTIPDIM)      !.. FLIP variable for O+ H+ & total ions
      DOUBLE PRECISION TI_SAVE(3,CTIPDIM)     !.. FLIP variable for Te and Ti
      DOUBLE PRECISION NHEAT(CTIPDIM)    !.. Neutral heating rate
      DOUBLE PRECISION RTS(99)           !.. Reaction rates
      DOUBLE PRECISION electron_density(CTIPDIM)     !.. electron density
      DOUBLE PRECISION O2DISF(CTIPDIM)   !.. O2 dissociation frequency
      DOUBLE PRECISION sza_north, sza_south !solar zenith angles at the north and south end of a tube                         
      INTEGER IWR         !$$$  for writing in files
      integer istop

      DATA M_TO_CM,M3_TO_CM3,nT_TO_Gauss/1.0E+2,1.0E-6,1.0E+4/    !.. Unit conversion factors
      DATA DEBUG/0/  !.. turn on debug writes if DEBUG=1

      DATA JTI/0/

      JTI=JTI+1

        nflag_t = 0
        nflag_d = 0

        Z(:)=0.0D0               !(FLDIM)
        SZA(:)=0.0D0
        BM(:)=0.0D0
        SL(:)=0.0D0
        GR(:)=0.0D0
        GL(:)=0.0D0

!ION_DEN_VEL
        XIONN(:,:)=0.0D0!(ISPEC,IDIM)
        XIONV(:,:)=0.0D0

!(2) SOLVARR
        DELTA(:)=0.0D0 !(10*SDIM)
        RHS(:)=0.0D0   !(10*SDIM)
        WORK(:)=0.0D0   !(50*SDIM)
        S(:)=0.0D0   !(50*SDIM)

!THERMOSPHERE
        ON(:)=0.0D0              !(TDIM)
        HN(:)=0.0D0
        N2N(:)=0.0D0
        O2N(:)=0.0D0
        HE(:)=0.0D0
        TN(:)=0.0D0
        UN(:)=0.0D0
        TINF(:)=0.0D0
        EHT(:,:)=0.0D0           !(3,TDIM)
        COLFAC=0.0D0 

!(3)AVE_PARAMS
      TEJ(:)=0.0D0 !(AVDIM)
      TIJ(:)=0.0D0 !(AVDIM)
      NUX(:,:)=0.0D0 !(2,AVDIM)
      RBM(:)=0.0D0 !(AVDIM)
      UNJ(:)=0.0D0 !(AVDIM)
      DS(:)=0.0D0 !(AVDIM)
      GRADTE(:)=0.0D0 !(AVDIM)
      GRADTI(:)=0.0D0 !(AVDIM)
      GRAV(:)=0.0D0 !(AVDIM)
      OLOSS(:)=0.0D0 !(AVDIM)
      HLOSS(:)=0.0D0 !(AVDIM)
      HPROD(:)=0.0D0 !(AVDIM)

!(1) production
      EUVION(:,:,:)=0.0D0
      PEXCIT(:,:,:)=0.0D0
      PEPION(:,:,:)=0.0D0
      PAUION(:,:,:)=0.0D0
      OTHPR1(:,:)=0.0D0
      OTHPR2(:,:)=0.0D0
      SUMION(:,:,:)=0.0D0 !(3,12,PDIM)
      SUMEXC(:,:,:)=0.0D0 !(3,12,PDIM)
      PAUEXC(:,:,:)=0.0D0 !(3,12,PDIM)
      NPLSPRD(:)=0.0D0 !(PDIM)
      PHION(:)=0.0D0 !(PDIM)

!(4)MINORNEUT
      N2D(:)=0.0D0 !(NDIM)
      N2P(:)=0.0D0 !(NDIM),
      N2A(:)=0.0D0 !(NDIM)
      O1D(:)=0.0D0 !(NDIM)
      O1S(:)=0.0D0 !(NDIM)
      EQN2D(:)=0.0D0 !(NDIM)
      N4S(:)=0.0D0 !(NDIM)
      NNO(:)=0.0D0 !(NDIM)
      !.. Load debug flag into EFLAG
      EFLAG(11,11) =  0 !DEBUG 

      !.. Check that dimensions agree
      IF(CTIPDIM.GT.FLDIM) THEN
        EFLAG(11,1) =-1
        RETURN
      ENDIF

      !.. Upload field line grid parameters to FIELD_LINE_GRID module
      JMIN=JMINX
      JMAX=JMAXX
      DO J=JMIN,JMAX
        Z(J)=ZX(J)
        SZA(J)=SZAX(J)
        BM(J)=BMX(J)*nT_TO_Gauss
        SL(J)=SLX(J)*M_to_CM
        GL(J)=GLX(J)
        GR(J)=GRX(J)*(M_to_CM)
      ENDDO

      !.. Upload densities and velocities to ION_DEN_VEL module
      DO J=JMIN,JMAX
      DO I=1,ISPEC
        XIONN(I,J)=XIONNX(I,J)*M3_to_CM3
        XIONV(I,J)=XIONVX(I,J)*M_to_CM
      ENDDO
!       XIONN(3,J)=0.0
!       XIONV(3,J)=0.0
!       XIONN(4,J)=0.0
!       XIONV(4,J)=0.0
! GHGM
      ENDDO

      !.. Transfer Te and Ti to FLIP variable TI
      DO J=JMIN,JMAX
      DO I=1,3
        temp_ti_te(I,J)=TE_TIX(I,J)
        EHT(I,J)=EHTX(I,J)
      ENDDO
      ENDDO

      !.. Upload thermosphere parameters to THERMOSPHERE module
      COLFAC=COLFACX  !.. O+ - O collision frequency Burnside factor 1-1.7 
      DO J=JMIN,JMAX
! GHGM - set N+ and He+ to zero
!       if(uthrs.lt.2.0E-2) then
!       XIONN(3,J)=0.0
!       XIONV(3,J)=0.0
!       XIONN(4,J)=0.0
!       XIONV(4,J)=0.0
!       endif
! GHGM
        ON(J)=OX(J)*M3_to_CM3
        HN(J)=HX(J)*M3_to_CM3
! GHGM
        if(hn(j).lt.1.e3) hn(j) = 1.e3
        N2N(J)=N2X(J)*M3_to_CM3
        O2N(J)=O2X(J)*M3_to_CM3
        HE(J)=HEX(J)*M3_to_CM3
        TN(J)=TNX(J)
        TINF(J)=TINFX(J)
        N4S(J)=N4SX(J)*M3_to_CM3
        NNO(J)=NNOX(J)*M3_to_CM3
        UN(J)=UNX(J)*M_to_CM
        !.. transfer densities from storage to FLIP solution variable N
        N(1,J)=XIONN(1,J)
        N(2,J)=XIONN(2,J)
        N(3,J)=(XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J))
        N(4,J)=XIONN(3,J)
        NHEAT(J)=0.0
        O2DISF(J)=0.0
! GHGM - intialize ti and te ......
!       if(uthrs.lt.2.0E-2) then
!       temp_ti_te(3,J)=904.0*DLOG(z(j))-3329.0
!       IF(temp_ti_te(3,J).GT.3000.) temp_ti_te(3,J)=3000.
!       temp_ti_te(1,J)=0.5*(TN(J)+temp_ti_te(3,J))
!       endif
! GHGM
      ENDDO

      !.. Set up initial temperature and density profiles.
      !.. 0.1 < HPEQ < 1.0.
      IF(HPEQ.GT.0.001) THEN
        CALL PROFIN(IHEPLS,INPLS,PCO,F107,N,temp_ti_te,HPEQ,HEPRAT)
      ENDIF

      !.. This routine adjusts the H+ and He+ densities for depleted flux tubes      
      !..  if HPEQ is negative. 0.1 < -HPEQ < 1.0
      IF(HPEQ.LT.-0.001) CALL NEW_HP(JMIN,JMAX,PCO,HPEQ,N,EFLAG)

      !.. Update solar EUV flux factors
      CALL FACEUV(F107,F107A,UVFAC,EUVFLUX)
!      CALL FACSEE(F107,F107A,UVFAC,EUVFLUX)
      
      !.. Update Schumann-Runge UV flux factors
      CALL FACSR(UVFAC,F107,F107A)

      !.... evaluate primary EUV production
      CALL PRIMPR(F107,F107A,UVFAC,COLUM,EUVFLUX)

      !.. electron density for photoelectron routine
      DO J=JMIN,JMAX
        electron_density(J)=XIONN(1,J)+XIONN(2,J)+XIONN(3,J)+XIONN(4,J)+
     >                      XIONN(5,J)+XIONN(6,J)
      ENDDO

      !.. 2-stream photoelectron routine to get electron heating 
      !.. rate and secondary ion production

      sza_north = sza(5) * 180.0 / 3.14159
      sza_south = sza(jmax - 4) * 180.0 / 3.14159

      IF(sza_north.LT.103.0.OR.sza_south.LT.103.0) then

      CALL PE2S(F107,F107A,N,temp_ti_te,FPAS,electron_density,UVFAC,
     >          COLUM,IHEPLS,INPLS,INNO,mp,lp)

      ENDIF

      !-- Sum the EUV, photoelectron, and auroral production rate
      CALL SUMPRD(JMIN,JMAX,AUR_PROD)

      !.. Loop to calculate O+(4S) total ionization rate
      !.. PHION=total O+(4S) prod, including EUV, e*, dissoc of O2 and
      !.. minor ions. BCKPRD = small background production to eliminate 
      !,, instability below 200 km
!     if((mp.eq.2).and.(lp.eq.27)) then
!       write(6,*) 'GHGM UT ',UTHRS
!       write(6,*) '**************'
!     endif
      DO J=JMIN,JMAX 
         PHION(J)=1.0E-22
         N(3,J)=1.0E-22
         DO I=1,9
            FD(I)=0.0
         ENDDO
         N2D(J)=0.0
         IF(Z(J).GE.80.AND.Z(J).LE.700) THEN
           !.. CALL cminor to get NO+, O+(2D), O2+, N2+ & O+(2P) densities
           CALL CMINOR(0,J,0,IHEPLS,INPLS,INNO,FD,7,N,temp_ti_te,
     >                 Z,EFLAG,mp,lp)
           BCKPRD=2.0E-10*N2N(J)*EXP(-5.0E-14*N2N(J)*TN(J))
           PHION(J)=SUMION(1,7,J)+SUMION(2,4,J)+SUMION(2,5,J)+FD(9)
           PHION(J)=PHION(J)+BCKPRD
         ELSE
           !.. Make sure minor ions and neutrals are 0.0 above upper boundary
           DO I=5,ISPEC
             XIONN(I,J)=0
           ENDDO
           N2D(J)=0.0
           N2A(J)=0.0
           O1D(J)=0.0
           O1S(J)=0.0
           N2P(J)=0.0
           NNO(J)=0.0
         ENDIF
         !.. Sum minor ions N+, NO+, O2+, N2+ for electron density at low altitudes
         N(3,J)=XIONN(4,J)+XIONN(5,J)+XIONN(6,J)+XIONN(7,J)+XIONN(8,J)
      ENDDO

      !..  electron and ion temperature solution
      n_save = n
      ti_save = temp_ti_te
      CALL TLOOPS(JMIN,JMAX,CTIPDIM,Z,N,temp_ti_te,DT,DTMIN,EFLAG,
     >            mp,lp,nflag_t) 
! GHGM - simple attempt to stop Te climbing above 5,000K
!     test_te = 0
!     DO J=JMIN,JMAX
!       if(temp_ti_te(3,j).ge.5000.) then
!         test_te = 1
!         goto 2345
!       endif
!     ENDDO
!2345 continue
!     if(test_te.eq.1) then
!!     write(6,*) 'GHGM TE greater than 5000K, resetting ',mp,lp
!     DO J=JMIN,JMAX
!       temp_ti_te(3,J)=904.0*DLOG(z(j))-3329.0
!       IF(temp_ti_te(3,J).GT.3000.) temp_ti_te(3,J)=3000.
!       temp_ti_te(1,J)=0.5*(TN(J)+temp_ti_te(3,J))
!     ENDDO
!     endif
! GHGM - simple attempt to stop Te climbing above 5,000K
!     if((mp.eq.2).and.(lp.eq.27)) then
!       do j = jmin,jmax
!         write(6,11) j,z(j),n(1,j),n(2,j),n(3,j),n(4,j),
!    >    temp_ti_te(1,j),temp_ti_te(2,j),temp_ti_te(3,j)           
!         write(6,11) j,z(j),XIONN(4,J),XIONN(5,J),XIONN(6,J), 
!    >                XIONN(7,J),XIONN(8,J)    
!         write(6,*) '*******************************'
!11       format(i6,f10.0,7e12.4)
!       enddo
!     endif
!     if(nflag_t.ne.0) then
!       write(6,*) 'GHGM flagged Temperature ', mp,lp,nflag_t
!       n = n_save
!       ti = ti_save
!       return
!     endif

      !.. O+, H+ solution
      n_save = n
      ti_save = temp_ti_te
      CALL DLOOPS(JMIN,JMAX,CTIPDIM,Z,N,temp_ti_te,DT,DTMIN,EFLAG,
     >            mp,lp,nflag_d)
!     if(nflag_d.ne.0) then
!       write(6,*) 'GHGM flagged Density ', mp,lp,nflag_d
!       n = n_save
!       ti = ti_save
!     endif

      !.. Recalculate the minor ion densities with the new O+ density
      !.. Added by PGR 2012-11-29
      DO J=JMIN,JMAX 
       IF(Z(J).GE.80.AND.Z(J).LE.700) THEN
          CALL CMINOR(0,J,0,IHEPLS,INPLS,INNO,FD,7,N,temp_ti_te,
     >                Z,EFLAG,mp,lp)
       ENDIF
      ENDDO

      !.. He+ solution
      IF(EFLAG(2,1).EQ.0.AND.IHEPLS.GT.0) THEN
        i_which_call = 1
        CALL XION(temp_ti_te,DT,DTMIN,9,EFLAG,mp,lp,i_which_call)
      ENDIF
      !.. N+ solution
      IF(EFLAG(2,1).EQ.0.AND.INPLS.GT.0) THEN
        i_which_call = 2
        CALL XION(temp_ti_te,DT,DTMIN,11,EFLAG,mp,lp,i_which_call)
      ENDIF

        !.. transfer densities from FLIP to CTIP variable
      DO J=JMIN,JMAX
        NNOX(J)=NNO(J)/M3_to_CM3
        DO I=1,ISPEC
          XIONNX(I,J)=XIONN(I,J)/M3_to_CM3
          XIONVX(I,J)=XIONV(I,J)/M_to_CM
        ENDDO
      ENDDO

      !.. Transfer Te and Ti and e heating to CTIPe variable
      DO J=JMIN,JMAX
        DO I=1,3
          TE_TIX(I,J)=temp_ti_te(I,J)
          EHTX(I,J)=EHT(I,J)
        ENDDO
      ENDDO

!        !.. Get neutral gas heating rate NHEAT
!        DO J=JMIN,JMAX
!          IF(Z(J).GE.80.AND.Z(J).LE.700) THEN
!            !.. electron density for photoelectron routine
!            electron_density(J)=XIONN(1,J)+XIONN(2,J)+XIONN(3,J)+XIONN(4,J)+
!     >      XIONN(5,J)+XIONN(6,J)
!            CALL RATS(J,temp_ti_te(3,J),temp_ti_te(1,J),TN(J),RTS)  !.. Reaction rates
!            !.. Neutral heating rate
!            CALL NEUT_HEATING(0,JMIN,J,Z(J),RTS,temp_ti_te(3,J),temp_ti_te(2,J),TN(J),
!     >        ON(J),O2N(J),N2N(J),HE(J),N4S(J),electron_density(J),N(1,J),XIONN(6,J)
!     >       ,XIONN(5,J),N(2,J),XIONN(7,J),XIONN(4,J),NNO(J),N2D(J)
!     >       ,N2P(J),N2A(J),XIONN(8,J),XIONN(9,J),O1D(J),O1S(J)
!     >       ,EHT(3,J),NHEAT(J),O2DISF(J))
!          ENDIF
!      ENDDO

      RETURN
      END
C:::::::::::::::::: WRITE_EFLAG ::::::::::::::::::::::::
C... This routine prints the information about error flags
C... Written by P. Richards September 2010
      SUBROUTINE WRITE_EFLAG(PRUNIT,   !.. Unit number to print results
     >                        EFLAG)   !.. Error flag array
      IMPLICIT NONE
      INTEGER PRUNIT,EFLAG(11,11)         !.. error flags
      IF(EFLAG(1,1).NE.0) WRITE(PRUNIT,11)
 11   FORMAT(/'  Convergence failure in Temperature solution (TLOOPS).'
     >  ,2X,'Time step less than minimum.')
      IF(EFLAG(1,2).NE.0) WRITE(PRUNIT,12)
 12   FORMAT(/'  Convergence failure in Temperature solution (TLOOPS).'
     >  ,2X,'Incorrect input to the band solver BDSLV.')

      IF(EFLAG(2,1).NE.0) WRITE(PRUNIT,21)
 21   FORMAT(/'  Convergence failure in O+ - H+ solution (DLOOPS).'
     >  ,2X,'Time step less than minimum.')
      IF(EFLAG(2,2).NE.0) WRITE(PRUNIT,22)
 22   FORMAT(/'  Convergence failure in O+ - H+ solution (DLOOPS).'
     >  ,2X,'Incorrect input to the band solver BDSLV.')

      IF(EFLAG(3,1).NE.0) WRITE(PRUNIT,31)
 31   FORMAT(/'  Convergence failure in He+ solution (XION).'
     >  ,2X,'Time step less than minimum.')
      IF(EFLAG(3,2).NE.0) WRITE(PRUNIT,32)
 32   FORMAT(/'  Convergence failure in He+ solution (XION).'
     >  ,2X,'Incorrect input to the band solver BDSLV.')

      IF(EFLAG(4,1).NE.0) WRITE(PRUNIT,41)
 41   FORMAT(/'  Convergence failure in He+ solution (XION).'
     >  ,2X,'Time step less than minimum.')
      IF(EFLAG(4,2).NE.0) WRITE(PRUNIT,42)
 42   FORMAT(/'  Convergence failure in N+ solution (XION).'
     >  ,2X,'Incorrect input to the band solver BDSLV.')

      IF(EFLAG(5,1).NE.0) WRITE(PRUNIT,51)
 51   FORMAT(/'  Convergence failure in CMINOR.')

      IF(EFLAG(11,1).NE.0) WRITE(PRUNIT,111)
 111  FORMAT(/3X,'** CTIP dimension not equal to FLIP dimensions'
     >  /3X,'** Check dimensions in all FLIP modules')
      RETURN
      END 
