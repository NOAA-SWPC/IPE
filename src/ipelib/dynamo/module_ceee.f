!
      module module_ceee
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: ceee
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: ceee
! !INTERFACE:
      subroutine ceee(cee,nx,ny,cf)
!t      use dynamo_module
      IMPLICIT NONE
!     !USES:
!
! !DESCRIPTION:
! Called from mudpack solvers to transfer coefficients.
!
      save
! !PARAMETERS: 
      integer,intent(in) :: nx,ny
      real,intent(in) :: cee(nx,ny,9)
! !RETURN VALUE:
      real,intent(out) :: cf(nx,ny,9)
!
! !REVISION HISTORY:
! 05.02.16  <Astrid Maute> <include header> 
! 
! EOP
!
      integer :: i,j,n
      do n = 1,9
        do j = 1,ny
          do i = 1,nx
            cf(i,j,n) = cee(i,j,n)
          enddo
        enddo
      enddo
      end subroutine ceee
!-----------------------------------------------------------------------
      end module module_ceee
!-----------------------------------------------------------------------
