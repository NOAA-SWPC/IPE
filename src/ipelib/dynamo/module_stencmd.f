!
      module module_stencmd
!----------------------------------------------------------------------- 
      PRIVATE
      PUBLIC :: stencmd
      CONTAINS
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: stencmd
! !INTERFACE:
      subroutine stencmd(zigm,nlon0,nlat0,cee,ncoef)
!t      use dynamo_module
      use module_htrpex,ONLY: htrpex
      use module_cnmmod,ONLY: cnmmod
      use module_cnm,ONLY: cnm
      IMPLICIT NONE
!     
! !DESCRIPTION: 
! Calculate contribution fo 3 by 3 stencil from coefficient zigm
! at each grid point and level.
!
! !RETURN VALUE:
      integer,intent(in) :: 
     |  nlon0, ! longitude dimension of finest grid level
     |  nlat0, ! latitude dimension of finest grid level
     |  ncoef  ! integer identifier of coefficient
      real,intent(in) :: 
     |  zigm(nlon0,nlat0)  ! coefficients (nlon0+1/2,(nlat0+1)/2) 
!     
! !RETURN VALUE:
      real,intent(inout) :: 
     |  cee(*)  ! output stencil array consisting of c0,c1,c2,c3,c4
!       
! !REVISION HISTORY:
! 18.02.05  <Astrid Maute> <include header> 
! 
! EOP
! Local:
      integer :: nc,nlon,nlat,n
!
! Perform half-way interpolation and extend zigm in array:
!
      call htrpex(zigm,nlon0,nlat0)
!
! Calculate contribution to stencil for each grid point and level:
!
      nc = 1
      nlon = nlon0
      nlat = nlat0
!
! Calculate modified and unmodified stencil on finest grid
!
      call cnmmod(nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
!
! Stencils on other grid levels remain the same.
      nc = nc+10*nlon*nlat 
      nlon = (nlon+1)/2
      nlat = (nlat+1)/2
!
      do n=2,5
        call cnm(nlon0,nlat0,nlon,nlat,cee(nc),ncoef)
        nc = nc+9*nlon*nlat
        if (n==1) nc = nc+nlon*nlat
        nlon = (nlon+1)/2
        nlat = (nlat+1)/2
      enddo ! n=1,5
      end subroutine stencmd
!-----------------------------------------------------------------------
      end module module_stencmd
!-----------------------------------------------------------------------
