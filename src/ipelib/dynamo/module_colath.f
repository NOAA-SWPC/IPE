!
      module module_colath
!
      PRIVATE
      PUBLIC :: colath
      contains 
!-----------------------------------------------------------------------
      subroutine colath(offset1_deg,offset2_deg)
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).  (crit is in cons module)
!
      use heelis_module !,ONLY:what?
      use module_magfield,only: sunlons 
      use cons_module,only: rtd,                                        &
     &  crit,pi_dyn,                                                    &! critical colatitudes crit(2)
     &  xlonm,xlatm  ! magnetic grid lons, lats
      use dynamo_module,only: kmlat0,pfrac
      use params_module,only: kmlonp1
      implicit none
!
      real,intent(in) :: offset1_deg,offset2_deg

      real,dimension(kmlonp1,kmlat0) :: colatc
!
! Local:
      integer :: i,j
      real :: sinlat,coslat,aslonc,ofdc,cosofc,sinofc
!
! offc(2), dskofc(2) are for northern hemisphere aurora (see aurora.F)
      ofdc = sqrt(offc(2)**2+dskofc(2)**2)
      cosofc = cos(ofdc)
      sinofc = sin(ofdc)
      aslonc = asin(dskofc(2)/ofdc)
!
! Define colatc with northern convection circle coordinates
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!
      do j=1,kmlat0
        sinlat = sin(abs(xlatm(j+kmlat0-1)))
        coslat = cos(    xlatm(j+kmlat0-1))
        do i=1,kmlonp1
          colatc(i,j) = cos(xlonm(i)-sunlons+aslonc)
          colatc(i,j) = acos(cosofc*sinlat-sinofc*coslat*colatc(i,j))
        enddo ! i=1,kmlonp1

!
! Calculate fractional presence of dynamo equation at each northern
! hemisphere geomagnetic grid point. Output in pfrac(kmlonp1,kmlat0)
!
         crit(1)=theta0(1)+(offset1_deg*pi_dyn/180.)
         crit(2)=theta0(1)+(offset2_deg*pi_dyn/180.)
	do i=1,kmlonp1
          pfrac(i,j) = (colatc(i,j)-crit(1))/(crit(2)-crit(1))
          if (pfrac(i,j) < 0.) pfrac(i,j) = 0.
          if (pfrac(i,j) >= 1.) pfrac(i,j) = 1.
	enddo ! i=1,kmlonp1

!       write(6,"('colath: j=',i3,' pfrac(:,j)=',/,(6e12.4))")
!     &    j,pfrac(:,j)

      enddo ! j=1,kmlat0

      end subroutine colath
!-----------------------------------------------------------------------
      end module module_colath
