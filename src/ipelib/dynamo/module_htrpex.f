!
      module module_htrpex
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: htrpex
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: htrpex
! !INTERFACE:
      subroutine htrpex(coeff,kmlon0,kmlat0)
      use dynamo_module,ONLY: wkarray
      IMPLICIT NONE
!     
! !DESCRIPTION: 
! copy coefficients into array and
! extend array over 16 grid points for the wrap around points
! on the 5 different grid levels. 
! Result returned in array(1:kmlath) - equator to pole.
!
! !ARGUMENTS:
      integer,intent(in) :: kmlon0,kmlat0        ! kmlon, kmlath
      real,intent(in) :: coeff(kmlon0,kmlat0)    ! Sigma(SH) pole to equator
!        
! !REVISION HISTORY:
! 18.02.05  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: i,j,jj
!
! Copy coeff into positions in array:
! MODIFY: NO SYM. USE NH VALUES
!
      do j=1,kmlat0      ! 1:kmlath
        do i=1,kmlon0
          wkarray(i,j) = coeff(i,j)
        enddo ! i=1,kmlon0
      enddo ! j=1,kmlat0
!
! Extend over 32 grid spaces to allow for wrap around points for
! a total of 5 grid levels:
!
      do i=1,16
        do j=1,kmlat0
          wkarray(1-i,j)      = wkarray(kmlon0-i,j) 
          wkarray(kmlon0+i,j) = wkarray(1+i,j)
        enddo ! j=1,kmlat0
      enddo ! i=1,16
!
      end subroutine htrpex
!-----------------------------------------------------------------------
      end module module_htrpex
!-----------------------------------------------------------------------
