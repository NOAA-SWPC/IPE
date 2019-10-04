!test022208: v2: how much degree can the polar cap conductanc affect the potential solution??? 
!
      module module_check_err
      contains
!-----------------------------------------------------------------------
      subroutine check_err(status,name)
      implicit none
!
!nm20140326
!#ifdef SUN
!#include "/opt/local/include/netcdf.inc"
!#else
!#include "netcdf.inc"
      include "netcdf.inc"
!#endif
!
      integer,intent(in) :: status
      character(len=*),intent(in) :: name
!
      write(6,*) 'error in ',name
      write(6,*) NF_STRERROR(status)
      write(6,*) ' '
!
      stop
! 
      end subroutine check_err
!-----------------------------------------------------------------------
      end module module_check_err

