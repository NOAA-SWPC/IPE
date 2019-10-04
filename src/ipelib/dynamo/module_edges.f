!
      module module_edges
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: edges
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: edges
! !INTERFACE:
      subroutine edges(c,nlon,nlat)
      use dynamo_module
      IMPLICIT NONE
!     
! !USES:
! !DESCRIPTION: 
! Insert polar boundary conditions in stencil c(nlon,nlat,9)
!    c(:,nlat,1:8) = 0 & diagonal term c(:,nlat,9) = 1.
!
! !ARGUMENTS:
      integer,intent(in) :: nlon,nlat    ! dimension of coeff. array
! !RETURN VALUE:
      real,intent(out) :: c(nlon,nlat,*) ! coefficient array
!        
! !REVISION HISTORY:
! 21.02.05  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: n,i
!
      do n=1,8
        do i=1,nlon
          c(i,nlat,n) = 0.
        enddo
      enddo
      do i=1,nlon
        c(i,nlat,9) = 1.
      enddo
      end subroutine edges
!-----------------------------------------------------------------------
      end module module_edges
!-----------------------------------------------------------------------
