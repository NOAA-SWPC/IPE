!
      module module_divide
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: divide
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: divide
! !INTERFACE:
      subroutine divide(c,nlon,nlat,nlon0,nlat0,cs,igrid)
      use dynamo_module
      IMPLICIT NONE
!     
! !USES:
! !DESCRIPTION: 
! Divide stencil C by cos(theta(i,j))
!
! !ARGUMENTS:
      integer,intent(in) :: nlon,nlat,   ! dim. of finest grid
     |    nlon0,nlat0,                   ! dim. of actual grid
     |    igrid                          ! flag for divide by cos lam_0
      real,intent(in) :: cs(*)           ! cos lam_0
! !RETURN VALUE:
      real,intent(out) :: c(nlon,nlat,*) ! coefficient array for grid level
!        
! !REVISION HISTORY:
! 21.02.05  <Astrid Maute> <include header> 
!                          <remove 1/cos of RHS at the equator> 
! 
! EOP
!
! Local:
      integer :: nint,j0,n,j,i
!
      nint = (nlon0-1)/(nlon-1)
      j0 = 1-nint
      do n = 1,9
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,n) = c(i,j,n)/(cs(j0+j*nint)*nint**2)
          enddo ! i = 1,nlon
        enddo ! j = 1,nlat-1
      enddo ! n = 1,9
!
      end subroutine divide
!-----------------------------------------------------------------------
      end module module_divide
!-----------------------------------------------------------------------
