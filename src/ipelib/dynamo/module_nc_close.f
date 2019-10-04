!
      module module_nc_close
!      
      PRIVATE
      PUBLIC :: nc_close
      contains
!--------------------------------------------------------------------------
      subroutine nc_close
      use nc_module
      use module_check_err,ONLY:check_err
      implicit none
!     
! Local 
      integer :: istat 
!#ifdef SUN
!#include "/opt/local/include/netcdf.inc"
!#elif AIX
!#include "netcdf.inc"
!#elif LINUX
!#include "netcdf.inc"
      include "netcdf.inc"
!#endif
!     
! close the new created file
      istat = nf_close(noid)
      if (istat /= NF_NOERR) call check_err(istat,'closing new file ')
!      
      end subroutine nc_close
!--------------------------------------------------------------------------
      end module module_nc_close
