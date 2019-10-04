!
      module module_potm
!
      PRIVATE
      PUBLIC :: potm
!     
      contains 
!-----------------------------------------------------------------------
      subroutine potm
      use heelis_module !,ONLY:
      use module_flwv32,ONLY: flwv32
      use module_magfield,only: sunlons 
      use cons_module,only: 
     |  xlonm,xlatm, ! magnetic grid lons, lats
     |  pi_dyn       ! pi used in dynamo calculations
      use params_module,ONLY: kmlat,kmlon,kmlonp1
      use module_flwv32,ONLY: flwv32
      use dynamo_module,ONLY: phihm
      implicit none
!
! Calculate heelis potential in geomagnetic coordinates.
!
! Local:
      integer :: i,j
      real,dimension(kmlon) :: dlat,dlon,ratio
      integer,dimension(kmlon) :: iflag
!
      ratio(:) = 1.
      do j=1,kmlat
        iflag(:) = 1 ! must be updated at each j
        dlat(:) = xlatm(j)
        dlon(:) = xlonm(1:kmlon)-sunlons
!
! flwv32 returns single-level Heelis potential in geomag coords:
!
        
        if (abs(xlatm(j)) > pi_dyn/6.) then
          call flwv32(dlat,dlon,ratio,iflag,kmlon,phihm(:,j),j)

!          write(6,"('potm: j=',i3,' phihm(:,j)=',/,(6e12.4))")
!     |       j,phihm(:,j)
        else
          phihm(1:kmlon,j) = 0.
        endif
      enddo ! j=1,kmlat
!
! Periodic points:
      do j=1,kmlat
        phihm(kmlonp1,j) = phihm(1,j)
      enddo ! j=1,kmlat
      end subroutine potm
!-----------------------------------------------------------------------
      end module module_potm
