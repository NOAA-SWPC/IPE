!
      module module_clearcee
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: clearcee
      contains
!-----------------------------------------------------------------------
 ! BOP
! !IROUTINE: clearce
! !INTERFACE:
      subroutine clearcee(cee,nlon0,nlat0)
      IMPLICIT NONE
!     
! !DESCRIPTION:
! Zero C arrays for stencil coefficients.
! Cee will contain:
!   c0(kmlon0,kmlat0,10), c1(kmlon1,kmlat1,9), c2(kmlon2,kmlat2,9),
!   c3(kmlon3,kmlat3,9),  c4(kmlon4,kmlat4,9)
!
! !PARAMETERS: 
      integer,intent(in) :: nlon0,nlat0
! !RETURN VALUE:
      real,intent(out) :: cee(*)
!
! !REVISION HISTORY:
! 05.03.8  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: nlon,nlat,n,m,i
!
! Compute total size of cee
      nlon = nlon0
      nlat = nlat0
      n = 0
      do m=1,5 ! 5 resolution levels
        n = n+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! m=1,5 ! 5 resolution levels
      n = 9*n+nlon0*nlat0
!
! Clear cee:
      do i=1,n
        cee(i) = 0.
      enddo
      end subroutine clearcee
!-----------------------------------------------------------------------
      end module module_clearcee
!-----------------------------------------------------------------------
