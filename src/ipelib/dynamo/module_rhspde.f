!
      module module_rhspde
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: rhspde
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: rhspde
! !INTERFACE:
      subroutine rhspde
!!USES:
      use dynamo_module     
      use cons_module,only: pi_dyn,dlatm,dlonm,r0
      IMPLICIT NONE
!     !DESCRIPTION:
!  differentiate RHS 
!  R_0 * d lam_m/d lam_0 * 1/cos lam_0 [ d K_(m phi)^DT(0) / d phi  +
!  (d [ K_(m lam)^DT(0) * cos(lam_0)]/ d lam_0  ] 
!
! !REVISION HISTORY:
! 05.02.16  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: i,j,jj
      real :: pi
      real,dimension(kmlat) :: tint1,tint2,tint3
      real,dimension(-1:kmlonp1+1) :: tint33
!
! External:
      real,external :: sddot ! in util.F
!
      pi = pi_dyn
!
! cos lam_0:
      do j=1,kmlat
        tint1(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo ! j=1,kmlat
!
! Perform differentiation (1/ cos lam0 ) * d K_(m phi)^DT(0)/ d phi
! central differencing
! rhs(i,j) = 1/ cos lam0(j) *  [K_(m phi)^DT(0)(i+1,j)-K_(m phi)^DT(0)(i-1,j)]/ (2 d lon)
! rim(1) =  K_(m phi)^DT(0)
! copy rim(1) into tint33 and set up wrap around points
!
      do j=2,kmlath-1   !  2,48                            
        jj = j+kmlath-1 ! 49,96
        do i=1,kmlon
          tint33(i) = rim(i,jj,1)         ! tint33(1:kmlon)
        enddo ! i=1,kmlon
        do i=1,2                          ! wrap around points
          tint33(i-2) = tint33(i-2+kmlon) ! -1:0 <= kmlon-1:kmlon
          tint33(i+kmlon) = tint33(i)     ! kmlon+1:kmlon+2 <= 1,2
        enddo ! i=1,2
        do i=1,kmlon                                         
          rhs(i,j) = 1./(dlonm*tint1(kmlath+j-1))*
     |      .5*(tint33(i+1)-tint33(i-1))  ! tint33 2:kmlon+1 and 0:kmlon-1
        enddo ! i=1,kmlon
      enddo ! j=2,kmlath-1
!
! Perform differentiation (1/ cos lam0 ) * d [K_(m lam)^DT(0) cos lam_0]/ d lam_0:
! central differencing
! rhs(i,j) = rhs(i,j) + 1/ cos lam0(j) *  
!  [K_(m lam)^DT(0)(i,j+1)cos lam0(j+1)-K_(m lam)^DT(0)(i,j-1)cos lam0(j-1)]/ (2 d lat)
!
      do j=kmlath+1,kmlat-1 ! 50,96                          
        jj = j-kmlath+1     !  2,48
        do i=1,kmlon
          rhs(i,jj) = rhs(i,jj)+1./(dlatm*tint1(j))*
     |      .5*(rim(i,j+1,2)*tint1(j+1)-rim(i,j-1,2)*tint1(j-1))
        enddo ! i=1,kmlon
      enddo ! j=kmlath+1,kmlat-1
!
! polar value:
      rhs(1,kmlath) = -2./float(kmlon)*
     |  sddot(kmlon,unitvm,rim(1,kmlat-1,2))/tint1(kmlat-1)
!
! equatorial value:
! [K_(m phi)^DT(0)(i+1,j)-K_(m phi)^DT(0)(i-1,j)]/ (2 d lon) +
! [K_(m lam)^DT(0)(i,j+1/2)cos lam0(j+1/2)-0]/ ( d lat/2)
!   with   -  K_(m lam)^DT(0)(i,j_eq) = 0
!          -  cos lam0(j_eq) =1
!
      i = 1
      rhs(i,1) = 0.5/dlonm*(rim(i+1,kmlath,1)-rim(kmlon,kmlath,1))
      rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(kmlath)*rim(i,kmlath,2)+
     |                            tint1(kmlath+1)*rim(i,kmlath+1,2))
      do i = 2,kmlon-1
        rhs(i,1) = 0.5/dlonm*(rim(i+1,kmlath,1)-rim(i-1,kmlath,1))
        rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(kmlath)*rim(i,kmlath,2)+
     |                              tint1(kmlath+1)*rim(i,kmlath+1,2))
      enddo ! i = 2,kmlon-1
      i = kmlon
      rhs(i,1) = 0.5/dlonm*(rim(1,kmlath,1)-rim(i-1,kmlath,1))
      rhs(i,1) = rhs(i,1)+1./dlatm*(tint1(kmlath)*rim(i,kmlath,2)+
     |                            tint1(kmlath+1)*rim(i,kmlath+1,2))
!
! Extend over longitude:
      do i=2,kmlon
        rhs(i,kmlath) = rhs(1,kmlath)
      enddo ! i=2,kmlon
!
! Periodic points:
      do j=1,kmlath
        rhs(kmlonp1,j) = rhs(1,j)
      enddo ! j=1,kmlath
!
! multiply by earth radius + h_0 in meter  = R0*1.E-2 
!  R_0 * d lam_m/d lam_0 * 1/cos lam_0 [ d K_(m phi)^DT(0) / d phi  +
!  (d [ K_(m lam)^DT(0) * cos(lam_0)]/ d lam_0  ] 
!
      do j=1,kmlath
        do i=1,kmlonp1
          rhs(i,j) = rhs(i,j)*r0*1.e-2
        enddo ! i=1,kmlonp1
      enddo ! j=1,kmlath
!
! Save rhs to secondary histories (redundant in vertical and
!  zero in north hem):
!     real,dimension(kmlonp1,kmlath)  :: rhs ! right-hand side 
!     real,dimension(kmlonp1,-2:nlev) :: rhs_plt ! diag
!
!     rhs_plt = 0.
!     do j=1,kmlat
!       if (j <= kmlath) then
!         do i=1,kmlonp1
!           rhs_plt(i,:) = rhs(i,j)
!         enddo
!       endif
!       call addfsech_ik('RHS',' ',' ',rhs_plt,
!    |    1,kmlonp1,nmlev,nmlev-1,j)
!     enddo ! j=1,kmlat
!
      end subroutine rhspde
!-----------------------------------------------------------------------
      end module module_rhspde
!-----------------------------------------------------------------------
