      module module_highlat
      
      PRIVATE
      PUBLIC :: highlat  
      
      contains 
!-----------------------------------------------------------------------
      subroutine highlat(offset1_deg_r8,offset2_deg_r8,potential_model,
     &           rc)   
!    
! am_2023.03 modify to include Weimer potential of both hemispheres & use 
! convection reversal boundary from Weimer; should use the average NH+SH
! potential in the dynamo equation & later overwrite the solved for potential
! in the NH & SH hemisphere 
!       
      use ipe_error_module
      use params_module,ONLY: kmlonp1,kmlat,kmlon,kmlath
      use dynamo_module,ONLY: kmlat0,phihm
      use module_sub_heelis,ONLY: heelis
      use module_colath,ONLY:colath
      use module_magfield
      use module_weimer2005Ipe,ONLY:setmodel2005Ipe,epotval2005Ipe,
     &  fileLocation,weictpoten
      use cons_module,ONLY:xlatm_deg,xlonm_deg,pi,dtr,xlatm
      use module_init_heelis
      use heelis_module,ONLY:ctpoten,theta0
      use cons_module,only: rtd,dtr

      implicit none
      real*8,intent(in) :: offset1_deg_r8,offset2_deg_r8
      integer,intent(in) :: potential_model
      real :: offset1_deg,offset2_deg
      integer,optional,intent(out) :: rc

      integer :: i,j,lrc,ihem
      real*8,parameter::fill=0. !fill in value outside the boundary
      real*8 ::mlat,mlt,epot,hstilt,hsangle
      real ::offset1_loc,offset2_loc               ! am2023.03 real*8 not specified in heelis & colat

      offset1_deg = real(offset1_deg_r8)
      offset2_deg = real(offset2_deg_r8)
!   
! am 10/04 remove weimer part: first test without any potential model
! Dynamo calls Heelis (heelis.F), Weimer (wei01gcm.F), or neither
!   for high latitude electric potential, according to user-provided
!   "model_potential".
! Get high latitude (Heelis or other) colatitudes, NH pfrac, and poten phihm.
!  If Weimer is used, then theta0,phid etc is changed before use in aurora
!   in dynamics.
!
      if (present(rc)) rc = IPE_SUCCESS

      if ( potential_model == 1 ) then ! 1 is HEELIS
        call init_heelis
        call heelis(offset1_deg,offset2_deg)

      else if ( potential_model == 2 ) then ! 2 is WEIMER2005

      ! am_2023.03 include southern hemisphere ihem = -1 
        ihem    = -1
        hstilt  = ihem*stilt      ! flip sign for southern hemisphere
        hsangle = ihem*sangle 
        ! write(6,"('highlat_setSH: ',(2e12.4))")  hstilt, hsangle      
       
        call setmodel2005Ipe(hsangle,bt,hstilt,swvel,swden,
     &       fileLocation,'epot',rc=lrc)
        if (ipe_error_check(lrc,msg="call to setmodel2005Ipe SH failed",
     &    rc=rc)) return
!       
        mlatLoopSH:do j=1,kmlath  ! loop over SH CHECK THAT THIS IS THE SH

          mlat=xlatm_deg(j)      !mlat[deg]

          mlonLoopSH: do i=1,kmlon
              !mltrad = mlonRad - sunlons + pi !mlt(rad)
              !sunlons is time dependent!
              mlt = (xlonm_deg(i)*dtr -sunlons + pi)*12./pi !mlt(hr)
              if (mlt>=24.) mlt=mod(mlt,24.0d0)
              if (mlt<0.) mlt=mlt+24.
              call epotval2005Ipe(abs(mlat),mlt,fill,epot)
              phihm(i,j) = epot*1.0E+3                      ! epot[kV]-->[Volt] 
              !write(6,"('highlatSH: ',(3e12.4))")xlonm_deg(i),mlat,epot             
          end do mlonLoopSH   !i=1,kmlon

          phihm(kmlonp1,j) = phihm(1,j) ! set periodic point

        end do mlatLoopSH      !j=1,kmlath SH
        
        ! determine cross polar cap potential & convection reversal boundary in SH 1
        call wei05loc(1,hsangle)

      ! am_2023.03 northern hemisphere ihem = 1 
        ihem    = 1
        hstilt  = ihem*stilt      ! no sign flip for northern hemisphere
        hsangle = ihem*sangle
       
        call setmodel2005Ipe(hsangle,bt,hstilt,swvel,swden,
     &       fileLocation,'epot',rc=lrc)
        if (ipe_error_check(lrc,msg="call to setmodel2005Ipe NH failed",
     &    rc=rc)) return
!       
        mlatLoopNH:do j=kmlath+1,kmlat  ! loop over NH -CHECK THAT THIS IS THE NH

          mlat=xlatm_deg(j)      !mlat[deg]

          mlonLoopNH: do i=1,kmlon
              !mltrad = mlonRad - sunlons + pi !mlt(rad)
              !sunlons is time dependent!
              mlt = (xlonm_deg(i)*dtr -sunlons + pi)*12./pi !mlt(hr)
              if (mlt>=24.) mlt=mod(mlt,24.0d0)
              if (mlt<0.) mlt=mlt+24.
              call epotval2005Ipe(abs(mlat),mlt,fill,epot)
              phihm(i,j) = epot*1.0E+3 ! epot[kV]-->[Volt]                       
          end do mlonLoopNH   !i=1,kmlon

          phihm(kmlonp1,j) = phihm(1,j) ! set periodic point

        end do mlatLoopNH      !j=1,kmlath SH

        ! determine cross polar cap potential & convection reversal boundary in SH & NH
        call wei05loc(2,hsangle)
        
        ctpoten     = 0.5*(weictpoten(1)+weictpoten(2))  ! avg SH& NH potential drop
        offset1_loc = -0.5*(theta0(1)+theta0(2))*rtd     ! avg convection reversal Bnd deg?
        offset2_loc = 15.
        call colath(offset1_loc,offset2_loc)
!        write(6,"('highlat_set: ',(7(x,e12.4)))") hstilt,hsangle,bt,
!     |       weictpoten(1),weictpoten(2),theta0(1),theta0(2)
        
      else  !  0 - potential_model='NONE'
        do j=1,kmlat0
          do i=1,kmlonp1
            phihm(i,j) = 0.
          enddo ! i=1,kmlonp1
        enddo ! j=1,kmlat0
        call colath(offset1_deg,offset2_deg)
      endif
      
      end subroutine highlat
!-----------------------------------------------------------------------
      subroutine wei05loc (ih,sangle)
! ih=1,2 for SH,NH called from subroutine highlat
! calculates convection reversal boundar theta(1/2)
!       offc(1/2)
!       ctpoten(1/2)
!
      use params_module,only: kmlat,kmlonp1,kmlath
      use heelis_module,only: dskofc,dskofa,offc,offa,psim,psie,pcen,
     |    phidp0,phidm0,phinp0,phinm0,rr1,phid,phin,theta0,ctpoten
      use module_weimer2005Ipe,only: weictpoten,bndyfitr
      use module_magfield,only: bt
      use cons_module,only: rtd,dtr
      use dynamo_module,only: phihm 
      
      implicit none
! Args:
      integer,intent(in) :: ih
      real*8,intent(in)  :: sangle  ! IMF angle considering hemisphere
      
      integer :: i,j,j1,j2
      integer :: jnx(2,2),inx(2,2)
      real*8 :: vnx(2,2),hem,mltd,mltn,byloc,crad,
     | offcdegp(2),dskofcp(2)
      
! set lat. loop boundaries ih=1 is SH, ih=2 is NH
	if (ih .eq. 1) then
	  j1 = 1
	  j2 = kmlath
	  hem = -1.
	else
	  j1 = kmlath + 1
	  j2 = kmlat
	  hem = 1.
	endif
! find cross polar cap potential drop in each hemisphere
!  Find min/max
      vnx(ih,1) = 0.
      vnx(ih,2) = 0.
      do j=j1,j2
        do i=1,kmlonp1-1
          if (phihm(i,j) .gt. vnx(ih,2)) then
            vnx(ih,2) = phihm(i,j)
            jnx(ih,2) = j
            inx(ih,2) = i
          endif
          if (phihm(i,j) .lt. vnx(ih,1)) then
            vnx(ih,1) = phihm(i,j)
            jnx(ih,1) = j
            inx(ih,1) = i
          endif
        enddo  !  i=1,nmlonp1-1
      enddo  !  j=j1,j2
      
      weictpoten(ih) = 0.001*(vnx(ih,2) - vnx(ih,1))
!      
      crad       = bndyfitr/2.  ! bndyfitr in setboundary2005Ipe calculate
      theta0(ih) = crad*dtr  
      !write(6,"('weiloc_tehta: ',(i4,x,e12.4))") ih,theta0(ih)
                    
      byloc = tan(sangle*dtr)/sqrt(1.+sangle**2*dtr**2)*bt     ! sangle [deg]; asymmetric wrt hemispheres - am2023.03 make it symmetric
      !write(6,"('weiloc_by: ',(i4,3(x,e12.4)))") ih,sangle,bt,byloc
      if (byloc .gt.  11.) byloc =  11.     ! this was byloc=7 before
      if (byloc .lt. -11.) byloc = -11.
!  day-night-time convection entrance in MLT converted to radians (f(By))
!  Use parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                                and phin (phin(MLT)=23.50 +/- 0.15By - 12)
      mltd = 9.39  - hem*0.21*byloc
      mltn = 23.50 - hem*0.15*byloc
      phid(ih) = (mltd-12.) * 15.*dtr 
      phin(ih) = (mltn-12.) * 15.*dtr 
      
!  Use default constant value of offcdegp from setboundary in Weimer 2005
      offcdegp(ih) = 4.2
      offc(ih) = offcdegp(ih)*dtr 
      offa(ih) = offcdegp(ih)*dtr 
     
      dskofcp(ih)= 0.
      dskofc(ih) = dskofcp(ih)*dtr
      dskofa(ih) = (dskofcp(ih)-2.5)*dtr  !  USED? oval offset is 2.5 deg towards dawn (more neg dskof)

! the following is just set since init_heelis will not be called
! check if any of these are used for Weimer
      psim(:) =  0.50 * ctpoten * 1000.
      psie(:) = -0.50 * ctpoten * 1000.
      pcen(:) = 0.
      phidp0(:) = 90.*dtr
      phidm0(:) = 90.*dtr
      phinp0(:) = 90.*dtr
      phinm0(:) = 90.*dtr
      rr1(:) = -2.6

      end subroutine wei05loc
!-----------------------------------------------------------------------
      end module module_highlat
