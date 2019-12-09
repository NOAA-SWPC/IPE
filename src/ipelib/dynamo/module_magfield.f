! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
!20120829: copied from magfield.f and splitted  
!
      module module_magfield
!
! Sub sunloc is called once per timestep from advance to determine
!   sun's longitudes for current ut model time.
! sunlons: sun's longitude in dipole coordinates (see sub sunloc)
! (this was dlons in earlier versions)
!
      real :: sunlons
      real(8) :: sangle,bt,stilt,swvel,swden
!-----------------------------------------------------------------------
      end module module_magfield
