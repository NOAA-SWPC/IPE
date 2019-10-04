!nm022007: copied originally from ./tiegcm1.8_dynamo_lres/util.F
!
! Utility subprograms for tgcm:
!
!-------------------------------------------------------------------
      real function sddot(n,x,y)
      implicit none
!
! Call sdot (single precision) if on Cray, or ddot (double precision) 
!   if on SGI. (ddot must be called even if -r8 on sgi compiler command 
!   line). Ddot is from -lblas on the sgi.
! On IBM AIX use dot_product()
!
! 2/10/00: removing incx,incy args (i.e., incx=incy=1 on unicos
!   and irix -- IBM dot_product does not use increment args --
!   this function must be called with stride-1 vectors 
!   (see bndry.f, bndry2.f, bndrya.f, threed.f, transf.f)
!
      integer,intent(in) :: n
      real,intent(in) :: x(n),y(n)
!
!nm022007: #ifdef UNICOS
!nm022007:       real,external :: sdot
!nm022007:       sddot = sdot(n,x,1,y,1)
!nm022007: #elif SUN
!nm022007:       real,external :: sdot
!nm022007:       sddot = dot_product(x,y)
!nm022007: #elif IRIX
!nm022007:       real,external :: ddot
!nm022007:       sddot = ddot(n,x,1,y,1)
!nm022007: #elif AIX
!nm022007:       sddot = dot_product(x,y)
!nm022007: #elif OSF1
!nm022007:       sddot = dot_product(x,y)
!nm022007: #elif LINUX
      sddot = dot_product(x,y)
!nm022007: #else
!nm022007:       write(6,"('>>> WARNING sddot: unresolved OS pre-processor',
!nm022007:      |  ' directive.')")
!nm022007: #endif
      end function sddot
