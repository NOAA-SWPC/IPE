!-----------------------------------------------------------------------
      module module_weimer2005Ipe
!module w05sc
      implicit none
!
! Data read from W05scEpot.dat or W05scBpot.dat:
      integer,parameter :: csize=28, d1_pot=15, d2_pot=18
      integer :: ab(csize), ls(csize), ms(csize)
      real*8 :: alschfits(d2_pot,csize), schfits(d1_pot,csize)
      real*8 :: ex_pot(2)
      integer :: maxl_pot,maxm_pot
!
! Data read from SCHAtable.dat:
      integer,parameter :: d1_scha=19, d2_scha=7, d3_scha=68
      real*8 :: allnkm(d1_scha,d2_scha,d3_scha)
      integer :: maxk_scha,maxm_scha
      real*8 :: th0s(d3_scha)
!
! Data read from W05scBndy.dat:
      integer,parameter :: na=6, nb=7
      real*8 :: bndya(na),bndyb(nb),ex_bndy(2)
!
! Location of W05scEpot.dat, W05SCHAtable.dat, W05scBndy.dat
!  character*7  :: file_location = 'moddta/'
      character(len=*),parameter  :: fileLocation = './'
!
      real*8 :: rad2deg,deg2rad ! set by setmodel
      real*8 :: bndyfitr        ! calculated by setboundary
      real*8 :: esphc(csize),bsphc(csize) ! calculated by setmodel
      real*8 :: tmat(3,3),ttmat(3,3) ! from setboundary
      integer,parameter :: mxtablesize=200
      real*8 :: plmtable(mxtablesize,csize)
      real*8 :: colattable(mxtablesize)
      real*8 :: nlms(csize)
!
      PRIVATE
      PUBLIC :: setmodel2005Ipe,epotval2005Ipe,fileLocation
      contains
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine setmodel2005Ipe(angle,bt,tilt,swvel,swden,file_path    &
     &,model)
      implicit none
!
! Args:
!
! file_path: directory in which to find data file (must have "/" at end)
! model: must be either 'epot' or 'bpot' for electric or magnetic potential
!
      real*8,intent(in) :: angle,bt,tilt,swvel,swden
      character(len=*),intent(in) :: file_path,model
!
! Local:
      integer :: i,j
      real*8 :: pi,stilt,stilt2,sw,swp,swe,c0
      real*8 :: rang,cosa,sina,cos2a,sin2a
      real*8 :: cfits(d1_pot,csize),a(d1_pot)
!
      if (trim(model) /= 'epot'.and.trim(model) /= 'bpot') then
         write(6,"('>>> model=',a)") trim(model)
         write(6,"('>>> setmodel: model must be either epot or bpot')")
         stop 'setmodel' 
      endif

!   write(6,"('setmodel: angle=',f8.2,' bt=',f8.2,' tilt=',f8.2,' swvel=',&
!     &f8.2,' swden=',f8.2)") angle,bt,tilt,swvel,swden
!
! Read data:
      if (trim(model) == 'epot') then
      call read_potential2005Ipe('global_idea_coeff_W05scEpot.dat')
      else
      call read_potential2005Ipe('global_idea_coeff_W05scBpot.dat')
      endif
      call schatable2005Ipe('global_idea_coeff_W05SCHAtable.dat')
!
      pi = 4.*atan(1.)
      rad2deg = 180./pi
      deg2rad = pi/180.
!     
      call setboundary2005Ipe(angle,bt,tilt,swvel,swden,file_path)
!
      stilt = sin(tilt*deg2rad)
      stilt2 = stilt**2
      sw = bt*swvel/1000.
      swe = (1.-exp(-sw*ex_pot(2)))*sw**ex_pot(1)
      c0 = 1.
      swp = swvel**2 * swden*1.6726d-6
      rang = angle*deg2rad
      cosa = cos(rang)
      sina = sin(rang)
      cos2a = cos(2.*rang)
      sin2a = sin(2.*rang)
      if (bt < 1.) then         ! remove angle dependency for IMF under 1 nT
         cosa = -1.+bt*(cosa+1.)
         cos2a = 1.+bt*(cos2a-1.)
         sina = bt*sina
         sin2a = bt*sin2a
      endif
      cfits = schfits           ! schfits(d1_pot,csize) is in module w05read_data
      a = (/c0      , swe       , stilt      , stilt2     , swp,        &
     &     swe*cosa, stilt*cosa, stilt2*cosa, swp*cosa,                 &
     &     swe*sina, stilt*sina, stilt2*sina, swp*sina,                 &
     &     swe*cos2a,swe*sin2a/)
      if (trim(model) == 'epot') then
         esphc(:) = 0.
         do j=1,csize
            do i=1,int(d1_pot)
               esphc(j) = esphc(j)+cfits(i,j)*a(i)
            enddo
         enddo
!     write(6,"('setmodel: esphc=',/,(6e12.4))") esphc
      else
         bsphc(:) = 0.
         do j=1,csize
            do i=1,int(d1_pot)
               bsphc(j) = bsphc(j)+cfits(i,j)*a(i)
            enddo
         enddo
!     write(6,"('setmodel: bsphc=',/,(6e12.4))") bsphc
      endif
      end subroutine setmodel2005Ipe
!-----------------------------------------------------------------------
      subroutine setboundary2005Ipe(angle,bt,tilt,swvel,swden,file_path &
     & )
      implicit none
!
! Args:
!
! file_path: directory in which to find data file (must have "/" at end)
      character(len=*),intent(in) :: file_path
      real*8,intent(in) :: angle,bt,tilt,swvel,swden
!
! Local:
      integer :: i
      real*8 :: swp,xc,theta,ct,st,tilt2,cosa,btx,x(na),c(na)

!   write(6,"('Enter setboundary: angle=',f8.3,' bt=',f8.3)") angle,bt
!
! Read data:
      call read_bndy2005Ipe('global_idea_coeff_W05scBndy.dat')
!
! Calculate the transformation matrix to the coordinate system
! of the offset pole.
!
      xc = 4.2
      theta = xc*(deg2rad)
      ct = cos(theta)
      st = sin(theta)
!
      tmat(1,:) = (/ ct, 0.d0, st/) 
      tmat(2,:) = (/ 0.d0, 1.d0, 0.d0/) 
      tmat(3,:) = (/-st, 0.d0, ct/)
!     
      ttmat(1,:) = (/ct, 0.d0,-st/)
      ttmat(2,:) = (/ 0.d0,1.d0, 0.d0/)
      ttmat(3,:) = (/st, 0.d0, ct/)
!
      swp = swden*swvel**2*1.6726d-6 ! pressure
      tilt2 = tilt**2
      cosa = cos(angle*deg2rad)
      btx = 1.-exp(-bt*ex_bndy(1))
      if (bt > 1.) then
         btx = btx*bt**ex_bndy(2)
      else
         cosa = 1.+bt*(cosa-1.) ! remove angle dependency for IMF under 1 nT
      endif
      x = (/1.d0, cosa, btx, btx*cosa, swvel, swp/)
      c = bndya
      bndyfitr = 0.
      do i=1,na
         bndyfitr = bndyfitr+x(i)*c(i)
      enddo

!   write(6,"('setboundary: cosa=',f8.3,' btx=',f8.3)") cosa,btx
!   write(6,"('setboundary: x=',/,(6e12.4))") x
!   write(6,"('setboundary: c=',/,(6e12.4))") c
!   write(6,"('setboundary: bndyfitr=',e12.4)") bndyfitr

      end subroutine setboundary2005Ipe
!-----------------------------------------------------------------------
      subroutine epotval2005Ipe(lat,mlt,fill,epot)
      implicit none
!
! Args:
      real*8,intent(in) :: lat,mlt,fill
      real*8,intent(out) :: epot
!     
! Local:
      integer :: inside,j,m,mm,skip
      real*8 :: z,phir,plm,colat,nlm
      real*8 :: phim(2),cospm(2),sinpm(2)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!     
      call checkinputs2005Ipe(lat,mlt,inside,phir,colat)
      if (inside == 0) then
         epot = fill
         return
      endif
!
! IDL code: 
! phim=phir # replicate(1,maxm) * ((indgen(maxm)+1) ## replicate(1,n_elements(phir)))
!   where the '#' operator multiplies columns of first array by rows of second array,
!   and the '##' operator multiplies rows of first array by columns of second array.
! Here, maxm == maxm_pot == 2 (from w05read_data module), and phir is a scalar. The 
!   above IDL statement then becomes: phim = ([phir] # [1,1]) * ([1,2] ## [phir]) where
!   phim will be dimensioned [1,2]
!
      phim(1) = phir
      phim(2) = phir*2.
!   write(6,"('epotval: phir=',1pe12.4,' phim=',2(1pe12.4))") phir,phim
      cospm(:) = cos(phim(:))
      sinpm(:) = sin(phim(:))
!
      z = 0.
      do j=1,csize
         if (skip == 1) then
            skip = 0
            cycle
         endif
         m = ms(j)
         if (ab(j)==1) then

            plm = scplm2005Ipe(j,colat,nlm) ! scplm function is in this module
!       write(6,"('epotval: j=',i3,' plm=',1pe12.4,' m=',i3)") &
!        j,plm,m

            skip = 0
            if (m == 0) then
               z = z+plm*esphc(j)
            else
               z = z+plm*(esphc(j)*cospm(m)+esphc(j+1)*sinpm(m))
               skip = 1
            endif

         endif                  ! ab(j)
      enddo
      epot = z 
      end subroutine epotval2005Ipe
!-----------------------------------------------------------------------
      subroutine mpfac2005Ipe(lat,mlt,fill,fac)
      implicit none
!
! Args:
      real*8,intent(in) :: lat,mlt,fill
      real*8,intent(out) :: fac
!     
! Local:
      integer :: j,m,inside,skip
      real*8 :: phim(2),cospm(2),sinpm(2),cfactor
      real*8 :: re,z,phir,plm,colat,nlm,pi
!     
      re = 6371.2 + 110.        ! km radius (allow default ht=110)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
      call checkinputs2005Ipe(lat,mlt,inside,phir,colat)
      if (inside == 0) then
         fac = fill
         return
      endif
!
      phim(1) = phir
      phim(2) = phir*2.
      cospm(:) = cos(phim(:))
      sinpm(:) = sin(phim(:))
!     
      z = 0.
      jloop: do j=1,csize
      if (skip == 1) then
         skip = 0
         cycle
      endif
      if (ls(j) >= 11) exit jloop
      m = ms(j)
      if (ab(j) == 1) then
         plm = scplm2005Ipe(j,colat,nlm) ! colat and nlm are returned (both reals)
         plm = plm*(nlm*(nlm+1.))
!
! bsphc was calculated in setmodel (when setmodel called with 'bpot')
         if (m==0) then
            z = z-plm*bsphc(j)
         else
            z = z-(plm*(bsphc(j)*cospm(m)+bsphc(j+1)*sinpm(m)))
            skip = 1
         endif
      endif
      enddo jloop               ! j=1,csize
      pi = 4.*atan(1.)
      cfactor = -1.d5/(4.*pi*re**2) ! convert to uA/m2
      z = z*cfactor
      fac = z
!   write(6,"('mpfac: lat=',f8.3,' mlt=',f8.3,' fac=',1pe12.4)") lat,mlt,fac
      end subroutine mpfac2005Ipe
!-----------------------------------------------------------------------
      real function scplm2005Ipe(indx,colat,nlm)
      implicit none
!
! Args:
      integer,intent(in) :: indx
      real*8,intent(in) :: colat
      real*8,intent(out) :: nlm
!
! Local:
      integer,save :: tablesize
      integer :: istat,i,j,l,m,skip
      real*8 :: th0,output(1),colata(1),plm1
      real*8 :: cth(mxtablesize)
      real*8,save :: prevth0=1.d36

      scplm2005Ipe = 0.
      th0 = bndyfitr
      if (prevth0 /= th0) then
         tablesize = 3*nint(th0)
         if (tablesize > mxtablesize) then 
            write(6,"('>>> tablesize > mxtablesize: tablesize=',i5,     &
     &' mxtablesize=',i5,' tn0=',e12.4)") tablesize,mxtablesize,th0
            stop 'tablesize'
         endif
!     write(6,"('scplm: indx=',i3,' colat=',f8.3,' th0=',e12.4,&
!       &' tablesize=',i3)") indx,colat,th0,tablesize
 
         do i=1,tablesize
            colattable(i) = float(i-1)*(th0/float(tablesize-1))
            cth(i) = cos(colattable(i)*deg2rad)
         enddo

!     write(6,"('scplm: tablesize=',i4,' colattable=',/,(6f8.3))") &
!       tablesize,colattable(1:tablesize)
!     write(6,"('scplm: tablesize=',i4,' cth=',/,(6(1pe12.4)))") &
!       tablesize,cth(1:tablesize)

         prevth0 = th0
         nlms = 0.              ! whole array init 
         do j=1,csize
            if (skip == 1) then
               skip = 0
               cycle
            endif
            l = ls(j)
            m = ms(j)
            
            nlms(j) = nkmlookup2005Ipe(l,m,th0) ! nkmlookup in this module

!       write(6,"('scplm after nkmlookup: j=',i3,' l=',i3,' m=',i3,&
!         &' nlms(j)=',f8.4)") j,l,m,nlms(j)

! real*8 :: plmtable(mxtablesize,csize)
            call pm_n2005Ipe(m,nlms(j),cth,plmtable(1:tablesize,j)      &
     &           ,tablesize)
!     write(6,"('scplm: j=',i3,' indx=',i3,' plmtable(:,j)=',/,(6e12.4))") &
!         j,indx,plmtable(1:tablesize,j)

            skip = 0
            if (m /= 0 .and. ab(j) > 0) then
               plmtable(1,j+1) = plmtable(1,j)
               nlms(j+1) = nlms(j)
               skip = 1
            endif

         enddo                  ! j=1,csize

      endif                     ! prevth0
      nlm = nlms(indx)
      colata(1) = colat
      
      call interpol_quad2005Ipe(plmtable(1:tablesize,indx)              &
     &     ,colattable(1:tablesize),colata,output)
      scplm2005Ipe = output(1)

!   write(6,"('scplm: indx=',i3,' scplm=',e12.4,' plmtable=',/,(6e12.4))") &
!     indx,scplm,plmtable(1:tablesize,indx)
!   write(6,"('scplm returning: indx=',i3,' scplm=',1pe12.4)") indx,scplm

      end function scplm2005Ipe
!-----------------------------------------------------------------------
      subroutine pm_n2005Ipe(m,r,cth,plmtable,tablesize)
      implicit none
!
! Args:
      integer,intent(in) :: m,tablesize
      real*8,intent(in) :: r
      real*8,intent(in) :: cth(tablesize)
      real*8,intent(out) :: plmtable(tablesize)
!     
! Local:
      integer :: i,k,ii
      real*8 :: rm,rk,div,ans,xn
      real*8,dimension(tablesize) :: a,x,tmp,table
!
      if (m == 0) then 
         a = 1.                 ! whole array op
      else
         do i=1,tablesize
            a(i) = sqrt(1.-cth(i)**2)**m
         enddo
      endif
      xn = r*(r+1.)
      x(:) = (1.-cth(:))/2.

!   write(6,"('pm_n: a=',/,(6(1pe12.4))") a
!   write(6,"('pm_n: xn=',1pe12.4") xn
!   write(6,"('pm_n: x=',/,(6(1pe12.4))") x

      table = a                 ! whole array init
!
      k = 1
!   write(6,"(/)")
      pmn_loop: do              ! repeat-until loop in idl code
      do i=1,tablesize
         rm = float(m)
         rk = float(k)
         a(i) = a(i)*(x(i)*((rk+rm-1.)*(rk+rm)-xn)/(rk*(rk+rm)))
         table(i) = table(i)+a(i) ! "result" in idl code
      enddo

!     write(6,"('pm_n: k=',i3,' a=',/,(6(1pe12.4)))") k,a
!     write(6,"('pm_n: k=',i3,' table=',/,(6(1pe12.4)))") k,table

      k = k+1
      do i=1,tablesize
         div = abs(table(i))
         if (div <= 1.d-6) div = 1.d-6
         tmp(i) = abs(a(i)) / div
      enddo

!     write(6,"('pm_n: k=',i3,' abs(a)=',/,(6(1pe12.4)))") k,abs(a)
!     write(6,"('pm_n: k=',i3,' abs(table)=',/,(6(1pe12.4)))") k,abs(table)
!     write(6,"('pm_n: k=',i3,' tmp=',/,(6(1pe12.4)))") k,tmp
!     write(6,"('pm_n: k=',i3,' max(tmp)=',1pe12.4)") k,maxval(tmp)

!     write(6,"('pm_n: k=',i5,' min,max table=',2(1pe12.4),' min,max tmp=',2(1pe12.4))")&
!       k,minval(table),maxval(table),minval(tmp),maxval(tmp)

      if (maxval(tmp) < 1.d-6) exit pmn_loop
      enddo pmn_loop
      ans = km_n2005Ipe(m,r)

!   write(6,"('pm_n: ans=',1pe12.4,' table=',/,(6(1pe12.4)))") ans,table(1:tablesize)

      plmtable(:) = table(:)*ans

!   write(6,"('pm_n returning: tablesize=',i4,' plmtable=',/,6(1pe12.4))") &
!     tablesize,plmtable

      end subroutine pm_n2005Ipe
!-----------------------------------------------------------------------
      real function km_n2005Ipe(m,rn)
      implicit none
!
! Args:
      integer,intent(in) :: m
      real*8,intent(in) :: rn
!
! Local:
      integer :: i,n
      real*8 :: rm
!
      if (m == 0) then 
         km_n2005Ipe = 1.
         return
      endif
    
      rm = float(m)
      km_n2005Ipe = sqrt(2.*exp(lngamma2005Ipe(rn+rm+1.)                &
     & -lngamma2005Ipe(rn-rm+1.))) / (2.**m*factorial2005Ipe(m))
!   write(6,"('km_n: m=',i3,' rn=',f8.4,' km_n=',e12.4)") m,rn,km_n

      end function km_n2005Ipe
!-----------------------------------------------------------------------
      real function nkmlookup2005Ipe(k,m,th0)
      implicit none
!
!     Args:
      integer,intent(in) :: k,m
      real*8,intent(in) :: th0
!
! Local:
      integer :: kk,mm
      real*8 :: th0a(1),out(1)

      if (th0 == 90.) then
         nkmlookup2005Ipe = float(k)
         return
      endif
      th0a(1) = th0
      kk = k+1
      mm = m+1
      if (kk > maxk_scha) then
         write(6,"('>>> nkmlookup: kk > maxk: kk=',i4,' maxk=',i4)") kk &
     &    ,maxk_scha
         call interpol_quad2005Ipe(allnkm(maxk_scha,mm,:),th0s,th0a,out &
     & )
      endif
      if (mm > maxm_scha) then
         write(6,"('>>> nkmlookup: mm > maxm: kk=',i4,' maxm=',i4)") kk &
     &    ,maxm_scha
         call interpol_quad2005Ipe(allnkm(kk,maxm_scha,:),th0s,th0a,out &
     & )
      endif
      if (th0 < th0s(1)) then
         write(6,"('>>> nkmlookup: th0 < th0s(1): th0=',e12.4,          &
     & ' th0s(1)=',e12.4)") th0,th0s(1)
      endif

!   write(6,"('nkmlookup call interpol: kk=',i3,' mm=',i3,' th0=',e12.4,&
!     &' allnkm=',/,(6(1pe12.4)))") kk,mm,th0a,allnkm(kk,mm,:)

      call interpol_quad2005Ipe(allnkm(kk,mm,:),th0s,th0a,out)

      nkmlookup2005Ipe = out(1)

      end function nkmlookup2005Ipe
!-----------------------------------------------------------------------
      subroutine checkinputs2005Ipe(lat,mlt,inside,phir,colat)
      implicit none
!
! Args:
      real*8,intent(in) :: lat,mlt
      integer,intent(out) :: inside
      real*8,intent(out) :: phir,colat
!
! Local:
      real*8 :: lon,tlat,tlon,radii
      
      lon = mlt*15.
      call dorotation2005Ipe(lat,lon,tlat,tlon)
      radii = 90.-tlat
      inside = 0
      if (radii <= bndyfitr) inside = 1 ! bndyfitr from setboundary
      phir = tlon*deg2rad
      colat = radii

      end subroutine checkinputs2005Ipe
!-----------------------------------------------------------------------
      subroutine dorotation2005Ipe(latin,lonin,latout,lonout)
      implicit none
!
! Args:
      real*8,intent(in) :: latin,lonin
      real*8,intent(out) :: latout,lonout
!
! Local:
      real*8 :: latr,lonr,stc,ctc,sf,cf,a,b,pos(3)
      integer :: i

      latr = latin*deg2rad
      lonr = lonin*deg2rad
      stc = sin(latr)
      ctc = cos(latr)
      sf = sin(lonr)
      cf = cos(lonr)
      a = ctc*cf
      b = ctc*sf
!
! IDL code: Pos= TM ## [[A],[B],[STC]]
! The ## operator multiplies rows of first array by columns of second array.
! Currently, TM(3,3) = Tmat (or TTmat if "reversed" was set)
! If called w/ single lat,lon, then a,b,stc are dimensioned (1), and
!   Pos is then (1,3)
!
      do i=1,3
         pos(i) = tmat(1,i)*a + tmat(2,i)*b + tmat(3,i)*stc
      enddo
  
      latout = asin(pos(3))*rad2deg
      lonout = atan2(pos(2),pos(1))*rad2deg

!   write(6,"('dorotation: latin,lonin=',2f9.4,' latout,lonout=',2f9.4)") &
!     latin,lonin,latout,lonout
  
      end subroutine dorotation2005Ipe
!-----------------------------------------------------------------------
      subroutine interpol_quad2005Ipe(v,x,u,p)
!
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
      implicit none
!
! Args:
      real*8,intent(in) :: v(:),x(:),u(:)
      real*8,intent(out) :: p(:)
!
! Local:
      integer :: nv,nx,nu,i,ix
      real*8 :: x0,x1,x2
!
      nv = size(v)
      nx = size(x)
      nu = size(u)
      if (nx /= nv) then
         write(6,"('>>> interpol_quad: nx /= nv: nx=',i4,' nv=',i4)")   &
     &    nx,nv
         p(:) = 0.
         return
      endif
      do i=1,nu
         ix = value_locate2005Ipe(x,u(i))
         if (ix <= 1.or.ix >= nx) then
!       write(6,"('>>> interpol_quad: ix out of range: nu=',i4,' ix=',i4)") nu,ix
            p(i) = 0.
            continue
         endif
         x1 = x(ix)
         x0 = x(ix-1)
         x2 = x(ix+1)
         if(x0.eq.0..and.x1.eq.0.) then
            p(i) = 0.0
         else
            p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2)/((x0-x1) * (x0-x2))+ &
     &             v(ix)   * (u(i)-x0) * (u(i)-x2)/((x1-x0) * (x1-x2))+ &
     &             v(ix+1) * (u(i)-x0) * (u(i)-x1)/((x2-x0) * (x2-x1))
         endif
      enddo
!   write(6,"('interpol_quad: nu=',i4,' p=',/,(1pe12.4)") nu,p

      end subroutine interpol_quad2005Ipe
!-----------------------------------------------------------------------
      integer function value_locate2005Ipe(vec,val)
!
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val >= vec(i+1)
! Input vec must be monotonically increasing
!
      implicit none
!
! Args:
      real*8,intent(in) :: vec(:),val
!
! Local:
      integer :: n,i
!
      value_locate2005Ipe = 0
      n = size(vec)
      if (val < vec(1)) return
      if (val > vec(n)) then
         value_locate2005Ipe = n
         return
      endif
      do i=1,n-1
         if (val >= vec(i) .and. val <= vec(i+1)) then
            value_locate2005Ipe = i
            return
         endif
      enddo

      end function value_locate2005Ipe
!-----------------------------------------------------------------------
      real function lngamma2005Ipe(xx)
!
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
!
      implicit none
      real*8,intent(in) :: xx
      real*8 :: x,y,tmp,ser
      real*8 :: cof(6) = (/76.18009172947146, -86.50532032941677        &
     &,24.01409824083091                                                &
     &,-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5/)
      integer :: j
!
      y = xx
      x = xx
      tmp = x+5.5
      tmp = tmp-(x+0.5)*log(tmp)
      ser = 1.000000000190015
      do j=1,5
         y = y+1
         ser = ser+cof(j)/y
      enddo
      lngamma2005Ipe = -tmp+log(2.5066282746310005*ser/x)
      end function lngamma2005Ipe
!-----------------------------------------------------------------------
      real function factorial2005Ipe(n)
      implicit none
      integer,intent(in) :: n
      integer :: m
      if (n <= 0) then
         write(6,"('>>> factorial: n must be positive: n=',i4)") n
         factorial2005Ipe = 0.
         return
      endif
      if (n == 1) then
         factorial2005Ipe = 1.
         return
      endif
      factorial2005Ipe = float(n)
      do m = n-1,1,-1
         factorial2005Ipe = factorial2005Ipe * float(m)
      enddo
      end function factorial2005Ipe
!-----------------------------------------------------------------------
      subroutine read_potential2005Ipe(infile)
!     
! Read ascii data file W05scEpot.dat or W05scBpot.dat, written by 
!   pro write_potential (write_data.pro)
!
      implicit none
!
! Args:
      character(len=*),intent(in) :: infile
!
! Local:
!
      character(len=16) :: fname
      integer :: i,lu=20
      integer :: csize_rd,d1_rd,d2_rd
!
      open(lu,file=infile,status='old', ACCESS ='SEQUENTIAL')
      read(lu,"(a)") fname
      read(lu,"(28i3)") ab
      read(lu,"(3i3)") csize_rd,d1_rd,d2_rd
      if (csize_rd /= csize) then
         write(6,"('>>> read_potential: file ',a                        &
     &,': incompatable csize: ','csize_rd=',i4,' csize=',i4)")          &
     &        fname,csize_rd,csize
         stop 'csize'
      endif
      if (d1_rd /= d1_pot) then
         write(6,"('>>> read_potential: file ',a,': incompatable d1: '  &
     &,'d1_rd=',i4,' d1_pot=',i4)") fname,d1_rd,d1_pot
         stop 'd1'
      endif
      if (d2_rd /= d2_pot) then
         write(6,"('>>> read_potential: file ',a,': incompatable d2: '  &
     &,'d2_rd=',i4,' d2_pot=',i4)") fname,d2_rd,d2_pot
         stop 'd2'
      endif
      do i=1,csize
         read(lu,"(6e20.9)") alschfits(:,i)
      enddo
      read(lu,"(2f10.3)") ex_pot
      read(lu,"(28i3)") ls
      read(lu,"(2i3)") maxl_pot,maxm_pot
      read(lu,"(28i3)") ms
      do i=1,csize
         read(lu,"(6e20.9)") schfits(:,i)
      enddo
      close(lu)
      end subroutine read_potential2005Ipe
!-----------------------------------------------------------------------
      subroutine schatable2005Ipe(infile)
!
! Read ascii data file SCHAtable.dat, written by pro write_scha
!   (write_data.pro)
!
      implicit none
!
! Args:
      character(len=*),intent(in) :: infile
!
! Local:
!
      character(len=16) :: fname
      integer :: i,j,lu=20
!
      open(lu,file=infile,status='old', ACCESS = 'SEQUENTIAL')
      read(lu,"(a)") fname
      read(lu,"(2i3)") maxk_scha,maxm_scha
      do i=1,d3_scha
         do j=1,d2_scha
            read(lu,"(6e20.9)") allnkm(:,j,i)
         enddo
      enddo
      read(lu,"(8f10.4)") th0s
      close(lu)
      end subroutine schatable2005Ipe
!-----------------------------------------------------------------------
      subroutine read_bndy2005Ipe(infile)
!
! Read ascii data file W05scBndy.dat, written by pro write_bndy
!   (write_data.pro)
!
      implicit none
!
! Args:
      character(len=*),intent(in) :: infile
!
! Local:
!
      character(len=16) :: fname
      integer :: rd_na,rd_nb,lu=20
!
      open(lu,file=infile,status='old', ACCESS = 'SEQUENTIAL')
      read(lu,"(a)") fname
      read(lu,"(2i3)") rd_na,rd_nb
      if (rd_na /= na) then
         write(6,"('>>> read_potential: file ',a,': incompatable na: ', &
     &'rd_na=',i4,' na=',i4)") fname,rd_na,na
         stop 'na'
      endif
      if (rd_nb /= nb) then
         write(6,"('>>> read_potential: file ',a,': incompatable nb: ', &
     &'rd_nb=',i4,' nb=',i4)") fname,rd_nb,nb
         stop 'nb'
      endif
      read(lu,"(8e20.9)") bndya
      read(lu,"(8e20.9)") bndyb
      read(lu,"(8e20.9)") ex_bndy
      close(lu)
      end subroutine read_bndy2005Ipe
!-----------------------------------------------------------------------
      end module module_weimer2005Ipe
!end module w05sc
