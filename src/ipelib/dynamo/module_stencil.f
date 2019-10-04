!
      module module_stencil
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: stencil
      contains
!-----------------------------------------------------------------------
 ! BOP
! !IROUTINE: stencil
! !INTERFACE:
      subroutine stencil(zigm,nlon0,nlat0,cee,ncoef)
!t      use dynamo_module
      use module_htrpex,ONLY: htrpex
      use module_cnm,ONLY: cnm
      IMPLICIT NONE
!     
! !DESCRIPTION:
! Calculate contribution fo 3 by 3 stencil from coefficient zigm
! at each grid point and level.
!
! !PARAMETERS: 
      integer,intent(in) :: 
     |  nlon0, ! longitude dimension of finest grid level
     |  nlat0, ! latitude dimension of finest grid level
     |  ncoef   ! integer identifier of coefficient
      real,intent(in) :: 
     |  zigm(nlon0,nlat0)  ! coefficients (nlon0+1/2,(nlat0+1)/2) 
!     
! !RETURN VALUE:
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!
! Local:
      integer :: nc,nlon,nlat,n
!
! Perform half-way interpolation and extend zigm in array:
!
      call htrpex(zigm,nlon0,nlat0)
!
! Calculate contribution to stencil for each grid point and level:
      nc = 1
      nlon = nlon0
      nlat = nlat0
      do n=1,5 ! 5 levels of resolution
        call cnm(nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencil
!-----------------------------------------------------------------------
      end module module_stencil
!-----------------------------------------------------------------------
