!
      module module_cnmmod
!----------------------------------------------------------------------- 
!
      PRIVATE
      PUBLIC :: cnmmod
      contains
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: cnmmod
! !INTERFACE:
      subroutine cnmmod(nlon0,nlat0,nlon,nlat,c,ncoef)
      use dynamo_module
      IMPLICIT NONE
!
!     !DESCRIPTION: 
! Compute contribution to stencil from zigm(ncoef) on finest grid nlon0 by nlat0.
! calculates the stencil with (c) and without (cofum)
! the use of upwinding method if diagonal dominance is not preserved
! output: c     ( upwinding)
!         cofum (no upwinding)
!
! ncoef: integer id of Sigma coefficient:
! ncoef = 1 for zigm11 Sigma_(phi phi)
! ncoef = 2 for zigmc  Sigma_(phi lam)
! ncoef = 3 for zigm2  Sigma_(lam phi)
! ncoef = 4 for zigm22 Sigma_(lam lam)
!
! !ARGUMENTS:
      integer,intent(in) :: 
     |  nlon0,nlat0, ! finest grid dimensions
     |  nlon,nlat    ! output grid dimensions
      integer,intent(in) :: ncoef
! !RETURN VALUE:
      real,intent(inout) :: 
     |  c(nlon,nlat,*)    ! output array for grid point stencils at
                          ! resolution nlon x nlat
!        
! !REVISION HISTORY:
! 18.02.05  <Astrid Maute> <include header> 
! 
! EOP
!
! Local:
      integer :: i,j,nint,i0,j0
      real,parameter :: pi=3.141592654
      real :: wk(nlon0,3)
!
! Compute separation of grid points of resolution nlon x nlat within
! grid of resolution nlon0,nlat0. Evaluate dlon and dlat, grid spacing 
! of nlon x nlat.
!
      nint = (nlon0-1)/(nlon-1)
!
! Scan wkarray nlon x nlat calculating and adding contributions to stencil
! from zigm(ncoef)
      i0 = 1-nint
      j0 = 1-nint
!
! zigm11: Sigma_(phi phi)
      if (ncoef==1) then
        do j = 1,nlat-1
          do i = 1,nlon
            c(i,j,1) = c(i,j,1)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            c(i,j,5) = c(i,j,5)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+(i+1)*nint,j0+j*nint)+
     |        2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
            cofum(i,j,1)=cofum(i,j,1)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            cofum(i,j,5)=cofum(i,j,5)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            cofum(i,j,9) = cofum(i,j,9)-
     |        .5*(wkarray(i0+(i+1)*nint,j0+j*nint)+
     |         2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigmc Sigma_(phi lam)
      elseif (ncoef==2) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,6) = c(i,j,6)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
            c(i,j,8) = c(i,j,8)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+(i+1)*nint,j0+j*nint))
            wk(i,1) = .5*(wkarray(i0+(i+1)*nint,j0+j*nint)-
     |        wkarray(i0+(i-1)*nint,j0+j*nint))
!
! Unmodified:
            cofum(i,j,2) = c(i,j,2)
            cofum(i,j,4) = c(i,j,4)
            cofum(i,j,6) = c(i,j,6)
            cofum(i,j,8) = c(i,j,8)
            cofum(i,j,3) = cofum(i,j,3)+wk(i,1)
            cofum(i,j,7) = cofum(i,j,7)-wk(i,1)
! test for upwinding
            wk(i,2) = (c(i,j,3)+wk(i,1))*(c(i,j,7)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,3)+c(i,j,7))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,3) = c(i,j,3)+wk(i,1)+wk(i,3)
            c(i,j,7) = c(i,j,7)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! zigm2 Sigma_(lam phi)
      elseif (ncoef==3) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,6) = c(i,j,6)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            c(i,j,8) = c(i,j,8)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
            wk(i,1) = .5*(wkarray(i0+i*nint,j0+(j+1)*nint)-
     |        wkarray(i0+i*nint,j0+(j-1)*nint))
!
! Unmodified:
            cofum(i,j,2) = c(i,j,2)
            cofum(i,j,4) = c(i,j,4)
            cofum(i,j,6) = c(i,j,6)
            cofum(i,j,8) = c(i,j,8)
            cofum(i,j,1) = cofum(i,j,1)+wk(i,1)
            cofum(i,j,5) = cofum(i,j,5)-wk(i,1)
! test for upwinding
            wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
            wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
            if (wk(i,2) >= 0.) wk(i,3) = 0.
            c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
            c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
            c(i,j,9) = c(i,j,9)-2.*wk(i,3)
          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,2) = c(i,j,2)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,4) = c(i,j,4)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))
          wk(i,1) = .5*(wkarray(i0+i*nint,j0+j*nint)+
     |      wkarray(i0+i*nint,j0+(j+1)*nint))

          cofum(i,j,2) = c(i,j,2)
          cofum(i,j,4) = c(i,j,4)
          cofum(i,j,1) = cofum(i,j,1)+wk(i,1)
          cofum(i,j,5) = cofum(i,j,5)-wk(i,1)

          wk(i,2) = (c(i,j,1)+wk(i,1))*(c(i,j,5)-wk(i,1))
          wk(i,3) = sign(wk(i,1),c(i,j,1)+c(i,j,5))
          if (wk(i,2) >= 0.) wk(i,3) = 0.
          c(i,j,1) = c(i,j,1)+wk(i,1)+wk(i,3)
          c(i,j,5) = c(i,j,5)-wk(i,1)+wk(i,3)
          c(i,j,9) = c(i,j,9)-2.*wk(i,3)
        enddo ! i=1,nlon
!
! zigm22: Sigma_(lam lam)
      elseif (ncoef==4) then
        do j = 2,nlat-1
          do i = 1,nlon
            c(i,j,3) = c(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
            c(i,j,7) = c(i,j,7)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j-1)*nint))
            c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+i*nint,j0+(j-1)*nint)
     |        +2.*wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
!
! Unmodified:
            cofum(i,j,3) = cofum(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j+1)*nint))
            cofum(i,j,7) = cofum(i,j,7)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |        +wkarray(i0+i*nint,j0+(j-1)*nint))
            cofum(i,j,9) = cofum(i,j,9)-.5*(wkarray(i0+i*nint,j0+(j-1)*
     |        nint)+2.*wkarray(i0+i*nint,j0+j*nint)+
     |        wkarray(i0+i*nint,j0+(j+1)*nint))

          enddo ! i = 1,nlon
        enddo ! j = 2,nlat-1
!
! Low latitude boundary condition:
        j = 1
        do i=1,nlon
          c(i,j,3) = c(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)
     |      +wkarray(i0+i*nint,j0+(j+1)*nint))
          c(i,j,9) = c(i,j,9)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                            wkarray(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,3) = cofum(i,j,3)+.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                                 wkarray(i0+i*nint,j0+(j+1)*nint))
          cofum(i,j,9) = cofum(i,j,9)-.5*(wkarray(i0+i*nint,j0+j*nint)+
     |                                 wkarray(i0+i*nint,j0+(j+1)*nint))

        enddo ! i=1,nlon
      endif ! ncoef
!
      end subroutine cnmmod
!-----------------------------------------------------------------------
      end module module_cnmmod
!-----------------------------------------------------------------------
