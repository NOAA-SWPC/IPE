!
      module heelis_module
!
! Module used to calculcate the Heelis model potential in both hemispheres
! Byimf, Ctpoten and Power at a minimum using paramaters from aurora_cons
!
!      use params_module,only: kmlat,kmlonp1,kmlon,kmlonp1
!      use dynamo_module,only: kmlat0,pfrac,phihm,potential_model
!      use cons_module,only: dtr  ! degrees to radians (pi/180)
!
      implicit none
      integer,parameter :: isouth=1, inorth=2

! am 10/04 defined in original code in input_module       
!     real, parameter ::  ctpoten = 45. ! cross-cap potential (kV)      (e.g., 45.)
      real ::  ctpoten ! cross-cap potential (kV)      (e.g., 45.)
!
! am 10/04 defined in aurora.F
! Additional auroral parameters (see sub aurora_cons):
! (dimension 2 is for south, north hemispheres)
      real ::                                                           &
     &  theta0(2),                                                      &! convection reversal boundary in radians
     &  offa(2),                                                        &! offset of oval towards 0 MLT relative to magnetic pole (rad)
     &  dskofa(2),                                                      &! offset of oval in radians towards 18 MLT (f(By))
     &  phid(2),                                                        &! dayside convection entrance in MLT converted to radians (f(By))
     &  phin(2)      ! night convection entrance in MLT converted to radians (f(By))
!
! am 10/04 defined in aurora.F
! The following parameters are used only by heelis module for dynamo:
      real ::                                                           &
     &  offc(2),                                                        &!
     &  dskofc(2),                                                      &!
     &  psim(2),                                                        &! 
     &  psie(2),                                                        &!
     &  pcen(2),                                                        &!
     &  phidp0(2),                                                      &!
     &  phidm0(2),                                                      &!
     &  phinp0(2),                                                      &!
     &  phinm0(2),                                                      &!
     &  rr1(2)
!
!-----------------------------------------------------------------------
      end module heelis_module
