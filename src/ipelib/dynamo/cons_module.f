!t #include "defs.h"
!
      module cons_module
      use params_module,only: kmlat,kmlon,kmlonp1
       
      implicit none
      
      real ::                                                           &
     &  pi,                                                             &! set with 4*atan(1)    C(110)
     &  dtr,                                                            &! degrees-to-radians (pi/180.)
     &  rtd                              ! radians-to-degrees (180./pi)
     
      real,parameter ::                                                 &
     &  pi_dyn=3.14159265358979312 ! pi used in dynamo calculations
      real,parameter ::                                                 &
     &  re   = 6.37122e8,                                               &! earth radius (cm)                
     &  h0 =9.0e6, r0 =re+h0    ! use mean earth radius
      real ::                                                           &
     &  dlatm, dlonm,                                                   &! grid spacing
     &  xlatm(kmlat),                                                   &! magnetic latitudes (radians)
     &  xlonm(kmlonp1),                                                 &! magnetic longitudes (radians)
     &  xlatm_deg(kmlat),                                               &! magnetic latitudes (degree)
     &  xlonm_deg(kmlonp1),                                             &! magnetic longitudes (degree)
     &  rcos0s(kmlat),                                                  &! cos(theta0)/cos(thetas)
     &  dt1dts(kmlat)  ! dt0dts/abs(sinim) (non-zero at equator)
!
! Critical colatitude limits for use of Heelis potential in dynamo:
      real,parameter ::                                                 &
     &  crit(2) = (/0.261799387, 0.523598775/) ! original values
!    &  crit(2) = (/0.523598775, 0.61086524/)  ! plasmasphere has zero  
!nm031407:     &  crit(2) = (/0.523598775, 0.525298775/)  !nm041106: test
                                          ! conductances aboce &lam&>60 deg therefore I set the pcb to 60deg
                                          ! and the auroral zone equator boundary to 55 deg
! now above &lam&>60 set the conductances to some value which shouldn't
! matter and below it's an interpolation between high and low latitude
! potential
!nm040607:     &  crit(2) = (/0.5314827561378479, 0.61086524/)  !crit(1) corresponds to max lat gip(83)=59.54828119329674  !nm040607:

      integer, public :: idyn_save(kmlat) !correspondance between lp_plas & lp_dyn

!-----------------------------------------------------------------------
      end module cons_module
