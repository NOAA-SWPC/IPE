!
      module module_sub_heelis
!
      PRIVATE
      PUBLIC :: heelis
      contains 
!-----------------------------------------------------------------------
      subroutine heelis(offset1_deg,offset2_deg)
      use heelis_module !,ONLY: what?
      use module_colath,ONLY:colath
      use module_potm,ONLY:potm
!     
! Heelis driver, called from sub dynamo (dynamo module, dynamo.F).
! These routines return pfrac and phihm to the dynamo.
!   (see argument descriptions below). 
!
      implicit none
      real,intent(in) :: offset1_deg,offset2_deg
!
! Args:
! pfrac:  Fractional presence of dynamo equation given critical 
!           convection colatitudes crit(2).
! phihm:  Heelis potential in magnetic coordinates (single level).
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).  (crit is in cons module)
!
      call colath(offset1_deg,offset2_deg)
!
! Calculate  the heelis potential phihm in geomagnetic coordinates:
! (potm calls sub flwv32)
!
      call potm
!
      end subroutine heelis
!-----------------------------------------------------------------------
      end module module_sub_heelis
