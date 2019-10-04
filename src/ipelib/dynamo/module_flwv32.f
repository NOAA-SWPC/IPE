!
      module module_flwv32
!
      PRIVATE
      PUBLIC :: flwv32
      contains 
!-----------------------------------------------------------------------
      subroutine flwv32(dlat,dlon,ratio,iflag,kmlon,poten,mlat)
!
! Calculate heelis potential at current magnetic latitude mlat.
!
      use heelis_module !,only: what?
      use cons_module,only: pi_dyn
      implicit none
!
! Args:
      integer,intent(in) :: mlat,kmlon
      integer,intent(inout) :: iflag(kmlon)
      real,dimension(kmlon),intent(in)  :: dlat,dlon,ratio
      real,dimension(kmlon+1),intent(out) :: poten
!
! Local:
      integer :: i,n,ihem
      real,parameter :: eps=1.e-6
      real ::                                                           &
     &  pi,pi2,pih,sinthr1,psi(8),phirc,sinth0,                         &
     &  ofda,cosofa(2),sinofa(2),aslona(2),                             &
     &  ofdc,cosofc(2),sinofc(2),aslonc(2),                             &
     &  phdpmx(2),phnpmx(2),phnmmx(2),phdmmx(2)
      real,dimension(kmlon) :: sinlat,coslat,sinlon,coslon,alon,colat,  &
     &  wk1,wk2,wk3,phifun,phifn2
      integer :: ifn(kmlon)
      real :: phi(kmlon,8)
!
      pi = pi_dyn
      pi2 = 2.*pi
      pih = .5*pi
      do n=1,2
        ofda = sqrt(offa(n)**2+dskofa(n)**2)
        cosofa(n) = cos(ofda)
        sinofa(n) = sin(ofda)
        aslona(n) = asin(dskofa(n)/ofda)
!
        ofdc = sqrt(offc(n)**2+dskofc(n)**2)
        cosofc(n) = cos(ofdc)
        sinofc(n) = sin(ofdc)
        aslonc(n) = asin(dskofc(n)/ofdc)
!
        if (phin(n) < phid(n)) phin(n) = phin(n)+pi2  ! modifies aurora phin
        phdpmx(n) = .5*min(pi,(phin(n)-phid(n)))
        phnpmx(n) = .5*min(pi,(phid(n)-phin(n)+pi2))
        phnmmx(n) = phdpmx(n)
        phdmmx(n) = phnpmx(n)
      enddo ! n=1,2

!     write(6,"('flwv32: mlat=',i3,' cosofa=',2e12.4)") mlat,cosofa
!     write(6,"('flwv32: mlat=',i3,' sinofa=',2e12.4)") mlat,sinofa
!     write(6,"('flwv32: mlat=',i3,' aslona=',2e12.4)") mlat,aslona
!     write(6,"('flwv32: mlat=',i3,' cosofc=',2e12.4)") mlat,cosofc
!     write(6,"('flwv32: mlat=',i3,' sinofc=',2e12.4)") mlat,sinofc
!     write(6,"('flwv32: mlat=',i3,' aslonc=',2e12.4)") mlat,aslonc
!     write(6,"('flwv32: mlat=',i3,' phdpmx=',2e12.4,' phnpmx=',
!    |  2e12.4)") mlat,phdpmx,phnpmx
!     write(6,"('flwv32: mlat=',i3,' phnmmx=',2e12.4,' phdmmx=',
!    |  2e12.4)") mlat,phnmmx,phdmmx

!
! Set ihem=1,2 for South,North hemisphere:
!  am 10/04 this doesn't make sense dlat(:) =  xlatm(j)
!      ihem = int(dlat(max0(1,nlon/2))*2./3.1416+2.) 
      ihem = int(dlat(1)*2./3.1416+2.) 
      sinth0 = sin(theta0(ihem))
!
! Average amie results show r1=-2.6 for 11.3 degrees
!   (0.1972 rad) beyond theta0.
!
      sinthr1 = sin(theta0(ihem)+0.1972)
      psi(1) = psie(ihem)
      psi(3) = psim(ihem)
      do n=2,4,2
        psi(n) = psi(n-1)
      enddo ! n=2,4,2
      do n=1,4
        psi(n+4) = psi(n)
      enddo ! n=1,4
!
! Transform to auroral circle coordinates:
!
      do i=1,kmlon
        sinlat(i) = sin(abs(dlat(i)))
        coslat(i) = cos(dlat(i))
        sinlon(i) = sin(dlon(i)+aslonc(ihem))
        coslon(i) = cos(dlon(i)+aslonc(ihem))
        colat(i) = cosofc(ihem)*sinlat(i)-sinofc(ihem)*coslat(i)*       &
     &    coslon(i)
        alon(i) = amod(atan2(+sinlon(i)*coslat(i),sinlat(i)*            &
     &    sinofc(ihem)+cosofc(ihem)*coslat(i)*coslon(i))-               &
     &    aslonc(ihem)+3.*pi,pi2)-pi
        colat(i) = acos(colat(i))*sqrt(ratio(i))
!
! Boundaries for longitudinal function:
!
        wk1(i) = ((colat(i)-theta0(ihem))/theta0(ihem))**2
        phi(i,4)=phid(ihem)+eps-min(phidm0(ihem)+wk1(i)*                &
     &    (pih-phidm0(ihem)),phdmmx(ihem))
        phi(i,5)=phid(ihem)-eps+min(phidp0(ihem)+wk1(I)*                &
     &    (pih-phidp0(ihem)),phdpmx(ihem))
        phi(i,6)=phin(ihem)+eps-min(phinm0(ihem)+wk1(i)*                &
     &    (pih-phinm0(ihem)),phnmmx(ihem))
        phi(i,7)=phin(ihem)-eps+min(phinp0(ihem)+wk1(i)*                &
     &    (pih-phinp0(ihem)),phnpmx(ihem))
        phi(i,1)=phi(i,5)-pi2
        phi(i,2)=phi(i,6)-pi2
        phi(i,3)=phi(i,7)-pi2
        phi(i,8)=phi(i,4)+pi2
        phifun(i)=0.
        phifn2(i) = 0.
        if (colat(i)-theta0(ihem) >= 0.) then
          ifn(i) = 3
        else
          ifn(i) = 2
        endif
        if (iflag(i) == 1) iflag(i) = ifn(i)
!
! Add ring current rotation to potential (phirc)
!
        phirc = 0.
        wk2(i) = amod(alon(i)+phirc+2.*pi2+pi,pi2)-pi
        wk3(i) = amod(alon(i)+phirc+3.*pi2,pi2)-pi
      enddo ! i=1,kmlon
!
! Longitudinal variation:
!
      do n=1,7
        do i=1,kmlon
          phifun(i)=phifun(i)+.25*(psi(n)+psi(n+1)+(psi(n)-             &
     &      psi(n+1))*cos(amod(pi*(wk2(i)-phi(i,n))/(phi(i,n+1)-        &
     &      phi(i,n)),pi2)))*(1.-sign(1.,(wk2(i)-phi(i,n))*(wk2(i)-     &
     &      phi(i,n+1))))
          phifn2(i)=phifn2(i)+.25*(psi(n)+psi(n+1)+(psi(n)-             &
     &      psi(n+1))*cos(amod(pi*(wk3(i)-phi(i,n))/(phi(i,n+1)-        &
     &      phi(i,n)),pi2)))*(1.-sign(1.,(wk3(i)-phi(i,n))*(wk3(i)-     &
     &      phi(i,n+1))))
        enddo
      enddo
!
! Evaluate total potential:
!
      do i=1,kmlon
        if (iflag(i)==2) then
          poten(i) = (2.*(pcen(ihem)-phifun(i))+(phifun(i)-phifn2(i))*  &
     &      0.75)*(colat(i)/theta0(ihem))**3 +                          &
     &      (1.5*(phifun(i)+phifn2(i))-3.*pcen(ihem))*(colat(i)/        &
     &      theta0(ihem))**2 + 0.75*(phifun(i)-phifn2(i))*(colat(i)/    &
     &      theta0(ihem)) + pcen(ihem)
        else
          poten(i) = phifun(i)*(max(sin(colat(i)),                      &
     &      sinth0)/sinth0)**rr1(ihem)*exp(7.*(1.-max(sin(colat(i)),    &
     &      sinthr1)/sinthr1))
        endif
      enddo

!      write(6,"(/'flwv32: j=',i2,' ihem=',i2)") mlat,ihem
!      write(6,"('  theta0(ihem)=',e12.4,' pcen(ihem)=',e12.4,
!     &  ' rr1(ihem)=',e12.4)") theta0(ihem),pcen(ihem),rr1(ihem)
!      write(6,"('  sinth0=',e12.4,' sinthr1=',e12.4)") sinth0,sinthr1
!      write(6,"('  iflag=',/,(20i3))") iflag
!      write(6,"('  colat=',/,(6e12.4))") colat
!      write(6,"('  phifun=',/,(6e12.4))") phifun
!      write(6,"('  phifn2=',/,(6e12.4))") phifn2
!      write(6,"('  poten =',/,(6e12.4))") poten
!      write(6,"(/'flwv32: j=',i2,' ihem=',i2,' poten=',/,(6e12.4))") 
!     &  mlat,ihem,poten   

      end subroutine flwv32
!-----------------------------------------------------------------------
      end module module_flwv32
