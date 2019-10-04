!
      module module_stenmod
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: stenmod
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: stenmod
! !INTERFACE:
      subroutine stenmod(inlon,inlat,c,phihm,pfrac)
      use dynamo_module,ONLY: kmlon0,kmlat0
!     
! !USES:
      use cons_module,only: dlatm,dtr
      IMPLICIT NONE
!
! !DESCRIPTION:
! Modify stencil c and cofum to set potential to heelis value within auroral circle.
!    c     is the stencil w   upwinding
! pfrac set up in heelis module: fractional presence of dynamo equation
!        pfrac = 1  for |colam_m| < crit(2)
!        pfrac = 0  for |colam_m| > crit(1)
!        pfrac = (colat_m-crit(1))/(crit(2)-crit(1))  
!                                for crit(2) < |lam_m| < crit(1)
! to ensure diagonal dominance only the diogonal element c(9) is used
! factor of (dlatm/(10.*dtr))**2 to eliminate dependence on grid spacing
! for c(9)
!
! !PARAMETERS: 
      integer,intent(in) :: inlon,inlat
      real,dimension(kmlon0,kmlat0),intent(in) :: 
     |  phihm,  ! heelis potential (from subs potm, flwv32)
     |  pfrac   ! fractional presence of dynamo (from sub colath)
! !RETURN VALUE:
      real,intent(inout) :: c(inlon,inlat,*)
! !REVISION HISTORY:
! 05.03.10  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: nint,i0,j0,i,j,n,jj
!
! Compute separation of grid points for this resolution:
      nint = (kmlon0-1)/(inlon-1)
      i0 = 1-nint
      j0 = 1-nint
!
! If nint==1, then we are at the highest resolution.
! Correct RHS, which is in c(10)
!
      if (nint==1) then
        do j=1,inlat
          do i=1,inlon
            c(i,j,10) = pfrac(i,j)*c(i,j,10)+(1.-pfrac(i,j))*c(i,j,9)*
     |        (dlatm/(10.*dtr))**2*phihm(i,j)
          enddo ! i=1,inlon
        enddo ! j=1,inlat
      endif
!
! Modify stencil, c(i,j,n),n=1,9:
!
      do j=1,inlat
        jj = j0+j*nint
        do n = 1,8
          do i = 1,inlon
            c(i,j,n) = c(i,j,n)*pfrac(i0+i*nint,jj)
          enddo ! i = 1,inlon
        enddo ! n = 1,8
        do i = 1,inlon
          c(i,j,9) = c(i,j,9)*pfrac(i0+i*nint,jj)+
     |      (1.-pfrac(i0+i*nint,jj))*c(i,j,9)*
     |      (dlatm*float(nint)/(10.*dtr))**2
        enddo ! i = 1,inlon
      enddo ! j=1,inlat
!
      end subroutine stenmod
!-----------------------------------------------------------------------
      end module module_stenmod
!-----------------------------------------------------------------------
