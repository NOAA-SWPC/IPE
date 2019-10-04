      module module_highlat
!
      PRIVATE
      PUBLIC :: highlat
      
      contains 
!-----------------------------------------------------------------------
      subroutine highlat
      use params_module,ONLY: kmlonp1
      use dynamo_module,ONLY: kmlat0,phihm,potential_model
      use module_sub_heelis,ONLY: heelis
      use module_colath,ONLY:colath
      implicit none
      integer :: i,j
      
!   
! am 10/04 remove weimer part: first test without any potential model
! Dynamo calls Heelis (heelis.F), Weimer (wei01gcm.F), or neither
!   for high latitude electric potential, according to user-provided
!   "model_potential".
! Get high latitude (Heelis or other) colatitudes, NH pfrac, and poten phihm.
!  If Weimer is used, then theta0,phid etc is changed before use in aurora
!   in dynamics.
!
      if (potential_model == 'HEELIS') then
        call heelis
      else  !  potential_model='NONE'
        do j=1,kmlat0
          do i=1,kmlonp1
	    phihm(i,j) = 0.
          enddo ! i=1,kmlonp1
        enddo ! j=1,kmlat0
        call colath
      endif
      
      end subroutine highlat
!-----------------------------------------------------------------------
      end module module_highlat
