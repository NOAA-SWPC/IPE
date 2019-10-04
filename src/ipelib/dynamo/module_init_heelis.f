!
      module module_init_heelis
!
      PRIVATE
      PUBLIC :: init_heelis
      contains 
!-----------------------------------------------------------------------
      subroutine init_heelis
      use heelis_module !,ONLY:offc, dskofc,phid, phin, psim,psie, pcen &
!     &,phidp0, phidm0, phinp0, phinm0,rr1,theta0,offa,dskofa
      use cons_module,only: dtr  ! degrees to radians (pi/180) 
      implicit none
!      
! The following parameters (offc through rr1) are used only in heelis 
!   potential calculation for the dynamo (see heelis.F)
!
      offc(isouth) = 1.*dtr
      offc(inorth) = 1.*dtr
      dskofc(isouth) = 0.
      dskofc(inorth) = 0.
      phid(isouth) = 0.
      phid(inorth) = 0.
      phin(isouth) = 180.*dtr
      phin(inorth) = 180.*dtr
      psim(:) =  0.50 * ctpoten * 1000.
      psie(:) = -0.50 * ctpoten * 1000.
      pcen(isouth) = 0.
      pcen(inorth) = 0.
      phidp0(:) = 90.*dtr
      phidm0(:) = 90.*dtr
      phinp0(:) = 90.*dtr
      phinm0(:) = 90.*dtr
      rr1(:) = -2.6
!
      theta0(isouth) = (-3.80+8.48*(ctpoten**0.1875))*dtr
      theta0(inorth) = (-3.80+8.48*(ctpoten**0.1875))*dtr
      offa(isouth) = 1.0*dtr
      offa(inorth) = 1.0*dtr
      dskofa(isouth) = 0.
      dskofa(inorth) = 0.
!
      end subroutine init_heelis
!-----------------------------------------------------------------------

      end module module_init_heelis
