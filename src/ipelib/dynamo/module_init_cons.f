!t #include "defs.h"
!
      module module_init_cons
      PRIVATE
      PUBLIC :: init_cons     
      contains
!-----------------------------------------------------------------------
      subroutine init_cons
      use cons_module      
      IMPLICIT NONE
! Local:
      integer :: j,i
      real,parameter :: r1=1.06e7, alfa=1.668
      real ::                                                           &
     &  tanth0(kmlat),                                                  &
     &  tanths(kmlat),                                                  &
     &  theta0(kmlat),                                                  &
     &  hamh0(kmlat)
      real :: tanths2
!
      pi = 4.*atan(1.)                ! C(110)
      dtr = pi/180.                   ! degrees to radians
      rtd = 180./pi                   ! radians to degrees
      dlatm = pi_dyn/float(kmlat-1) ! note use of pi_dyn
      dlonm = 2.*pi_dyn/float(kmlon) 

! Set magnetic latitudes xlatm and magnetic longitudes xlonm:
!
! xlatm is equally spaced in theta0, but holds corresponding value
!   of thetas.
      do j=1,kmlat
        theta0(j) = -pi_dyn/2.+float(j-1)*dlatm ! note use of pi_dyn
      enddo ! j=1,kmlat
      do j=2,kmlat-1
        tanth0(j) = abs(tan(theta0(j)))
        hamh0(j) = r1*tanth0(j)+r0*tanth0(j)**(2.+2.*alfa)/             &
     &    (1.+tanth0(j)**2)**alfa
        tanths(j) = sqrt(hamh0(j)/r0)
        xlatm(j) = sign(atan(tanths(j)),theta0(j))
        xlatm_deg(j) = xlatm(j)*rtd
        rcos0s(j) = sqrt((1.+tanths(j)**2)/(1.+tanth0(j)**2))
        tanths2  = tanths(j)**2
        dt1dts(j) =                                                     &
     &    (r0*sqrt(1.+4.*tanths2)*(1.+tanths2))/                        &
     &    (r1*(1.+tanth0(j)**2)+2.*r0*tanth0(j)**(2.*alfa+1.)*          &
     &    (1.+alfa+tanth0(j)**2)/(1.+tanth0(j)**2)**alfa)
        
      enddo ! j=2,kmlat-1
      
! Magnetic poles:
      xlatm(1)     = theta0(1)
      xlatm(kmlat) = theta0(kmlat)
      xlatm_deg(1)     = xlatm(1)*rtd    
      xlatm_deg(kmlat) = xlatm(kmlat)*rtd 
      rcos0s(1)    = 1.
      rcos0s(kmlat)= 1.
      
! Magnetic longitudes:
      do i=1,kmlonp1
        xlonm(i) = -pi_dyn+float(i-1)*dlonm
        xlonm_deg(i) = xlonm(i)*rtd
      enddo ! i=1,kmlonp1

! output grid
! GHGM get rid of these write statements for now....
!     print *,'(20) output dyn grid'
!      write(unit=4020,FMT='(I12)')utime
!     write(unit=4020,FMT='(20E12.4)')xlatm_deg
!     write(unit=4020,FMT='(20E12.4)')xlonm_deg
!     write(unit=4020,FMT='(20E12.4)')crit(1),crit(2)
!     print *,'sub-init_con:crit1/2=',(90.-crit(1:2)*rtd) 
!
      
      end subroutine init_cons
!-----------------------------------------------------------------------
      end module module_init_cons
