!20111107: copied originally from tiegcm1.8_dynamo_lres
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
      module module_threed
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: threed
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: threed
! !INTERFACE:
      subroutine threed
      use params_module,only: kmlat,kmlon,kmlonp1,kmlath
      use dynamo_module,only: ed1dy,ed2dy,phim
!     
! !USES:
      use cons_module,only: dlatm,dlonm,rcos0s,r0,pi_dyn,dt1dts
!t      use module_sub_ncplot,ONLY:ncplot
      IMPLICIT NONE
!     
! !DESCRIPTION: 
! calculates electric field 
!      E_d1 = -1/(R cos lam_m) d Phi/ d phi_m
!      E_d2 = 1/ (R |sinI_m|)  d Phi/ d lam_m
! calculate 3D electric potential and electric field assuming
! that the fieldlines are dipolar and with constant electric potential
! along the field line
!
! !REVISION HISTORY:
! 07.03.05  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: i,j,ip1f,ip2f,ip3f
      real :: pi
      real :: csth0
!nm101306:
      character :: fname*10,labl*56,units*12
!
! Externals:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
!
! Calculate E_d1, E_d2 components of electric field:
! E_d1 = -1/(R cos lam_0)* (cos lam_0/cos lam_m) d Phi/ d phi_m
      do j=2,kmlat-1
        csth0 = cos(-pi/2.+(j-1)*dlatm)
        do i=2,kmlon
          ed1dy(i,j) = -(phim(i+1,j)-phim(i-1,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
        enddo ! i=2,kmlon
        ed1dy(1,j) = -(phim(2,j)-phim(kmlon,j))/(2.*dlonm*csth0)*
     |      rcos0s(j)/(r0*1.e-2)
        ed1dy(kmlonp1,j) = ed1dy(1,j)
      enddo ! j=2,kmlat-1
!
! E_d2 = +/- 1/(R |sinI_m|) d Phi/d lam_0* (d lam_0/ (d lam_m*|sin I_m|))
! Southern hemisphere:
      do j=2,kmlath-1
        do i=1,kmlonp1
          ed2dy(i,j) = -(phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
        enddo ! i=1,kmlonp1
      enddo ! j=2,kmlath-1
!
! Northern hemisphere:
      do j=kmlath+1,kmlat-1
        do i=1,kmlonp1
          ed2dy(i,j) = (phim(i,j+1)-phim(i,j-1))/(2.*dlatm)*dt1dts(j)/
     |      (r0*1.e-2)
        enddo ! i=1,kmlonp1
      enddo ! j=kmlath+1,kmlat-1
!
! Poles:
      do i = 1,kmlonp1
        ip1f = i + kmlon/4
        if (ip1f > kmlonp1) ip1f = ip1f - kmlon
        ip2f = i + kmlon/2
        if (ip2f > kmlonp1) ip2f = ip2f - kmlon
        ip3f = i + 3*kmlon/4
        if (ip3f > kmlonp1) ip3f = ip3f - kmlon
        ed1dy(i,1) = .25*(ed1dy(i,2) - ed1dy(ip2f,2) +
     |                  ed2dy(ip1f,2) - ed2dy(ip3f,2))
        ed1dy(i,kmlat) = .25*(ed1dy(i,kmlat-1) - ed1dy(ip2f,kmlat-1) +
     |                      ed2dy(ip1f,kmlat-1) - ed2dy(ip3f,kmlat-1))
        ed2dy(i,1) = .25*(ed2dy(i,2) - ed2dy(ip2f,2) -
     |                  ed1dy(ip1f,2) + ed1dy(ip3f,2))
        ed2dy(i,kmlat) = .25*(ed2dy(i,kmlat-1) - ed2dy(ip2f,kmlat-1) -
     |                      ed1dy(ip1f,kmlat-1) + ed1dy(ip3f,kmlat-1))
!
! Equator
       ed2dy(i,kmlath) = (4.*phim(i,kmlath+1)-phim(i,kmlath+2)
     |    -3.*phim(i,kmlath))/(2.*dlatm)/(R0*1.e-2)
       ed2dy(i,kmlath+1) = (4.*phim(i,kmlath+2)-phim(i,kmlath+3)
     |    -3.*phim(i,kmlath+1))/(2.*dlatm)/(R0*1.e-2)
       ed2dy(i,kmlath-1) = (4.*phim(i,kmlath-2)-phim(i,kmlath-3)
     |    -3.*phim(i,kmlath-1))/(2.*dlatm)/(R0*1.e-2)
      enddo ! i = 1,kmlonp1
!
!
      write(unit=4022,FMT='(20E12.4)') ed1dy
      write(unit=4023,FMT='(20E12.4)') ed2dy

      end subroutine threed
!-----------------------------------------------------------------------
      end module module_threed
!-----------------------------------------------------------------------
