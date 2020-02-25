!20111107: copied originally from tiegcm1.8_dynamo_lres
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!-------------------------------------------- 
!
      module dynamo_module
!----------------------------------------------------------------------- 
! BOP
! !MODULE: dynamo_module
! !DESCRIPTION: 
!
! Module for electrodynamic dynamo (the "e" in tiegcm). 
!
! 5/02 btf: This is the "new" dynamo code, adapted by Astrid Maute
! and Art Richmond from the original ("old") dynamo in tgcm15,
! written by E. Cicely Ridley completed in about 1992. 
!
! Sub prep_dynamo is called from advance. Prep_dynamo prepares dynamo 
! input fields and gathers them to the root task (mp_dynamo_gather).
! Next, advance calls sub dynamo (the dynamo driver), which controls
! the dynamo subroutine calls. Prep_dynamo and the dynamo driver are
! called once per timestep, after sub dynamics.
!
! Since the dynamo is serial (as of 4/02), the master task must receive
! subdomain data from other tasks for input fields to the dynamo 
! (on the geographic grid). These fields are exchanged via sub 
! mp_dynamo_gather which is called from prep_dynamo below.
!
! !USES:
!
      use params_module,only:                                           &
     &  kmlon,                                                          &! number of geomagnetic grid longitudes
     &  kmlonp1,                                                        &! kmlon+1
     &  kmlat,                                                          &! number of geomagnetic grid latitudes
     &  kmlatp1,                                                        &! kmlat+1
     &  kmlath  ! (kmlat+1)/2 (index to magnetic equator)
!
! Dimensions of the 5 grid resolutions for the multi-grid PDE:
      integer,parameter ::                                              &
     &  kmlon0=kmlon+1,                                                 &
     &  kmlat0=(kmlat+1)/2,                                             &
     &  kmlon1=(kmlon0+1)/2,                                            &
     &  kmlat1=(kmlat0+1)/2,                                            &
     &  kmlon2=(kmlon1+1)/2,                                            &
     &  kmlat2=(kmlat1+1)/2,                                            &
     &  kmlon3=(kmlon2+1)/2,                                            &
     &  kmlat3=(kmlat2+1)/2,                                            &
     &  kmlon4=(kmlon3+1)/2,                                            &
     &  kmlat4=(kmlat3+1)/2
!
! Space needed for descretized coefficients of dynamo pde at all 
!   5 levels of resolution:
!
      integer,parameter ::                                              &
     &  ncee=10*kmlon0*kmlat0+9*(kmlon1*kmlat1+kmlon2*kmlat2+kmlon3*    &
     &    kmlat3+kmlon4*kmlat4)
!
! Coefficients are stored in 1-d array cee(ncee)
! cee transmits descretized dynamo PDE coefficients to the multi-grid 
!   mudpack solver. (cee was formerly in ceee.h)
! The common block /cee_com/ is retained from earlier versions because
!   of the equivalencing below of coefficient arrays c0, c1, etc.
!
      real :: cee(ncee) 
      common/cee_com/ cee
!
! The following parameters nc0,nc1,... are pointers to the beginning of 
!   the coefficients for each level of resolution.
!
      integer,parameter ::                                              &
     & nc0=1,                                                           &
     & nc1=nc0+10*kmlon0*kmlat0,                                        &
     & nc2=nc1+9 *kmlon1*kmlat1,                                        &
     & nc3=nc2+9 *kmlon2*kmlat2,                                        &
     & nc4=nc3+9 *kmlon3*kmlat3
!
! nc(1:6) are pointers to beginning of coefficient blocks at each of 
!   5 levels of resolution: 
! nc(1) = nc0, pointer to coefficients for highest resolution.
! nc(2) = nc1, pointer to coefficients at half the resolution of nc0, 
!   and so on for nc(3), nc(4), nc(5), etc.
! nc(6) = ncee, the dimension of the entire cee array, containing
!   coefficients for all 5 levels of resolution.
! 
      integer :: nc(6)
!
      real ::                                                           &
     &  c0(kmlon0,kmlat0,10),                                           &
     &  c1(kmlon1,kmlat1,9),                                            &
     &  c2(kmlon2,kmlat2,9),                                            &
     &  c3(kmlon3,kmlat3,9),                                            &
     &  c4(kmlon4,kmlat4,9)
      equivalence                                                       &
     &  (cee,c0),                                                       &
     &  (cee(nc1),c1),                                                  &
     &  (cee(nc2),c2),                                                  &                         
     &  (cee(nc3),c3),                                                  &
     &  (cee(nc4),c4)
!
! Unmodified coefficients for using modified mudpack:
      real,dimension(kmlon0,kmlat0,9) :: cofum
!
! phim: Single level dynamo potential in geomagnetic coordinates, as
!   output by PDE solver mud (formerly in dynphi.h):
! The dynamo is symmetric about the magnetic equator, but the high latitude
!  is anti-symmetric in both hemispheres.  However, since Mudpack uses the
!  NH potential pattern, then the SH potential pattern must be added
!  back into the 2-D phim before the call threed, and before it is
!  transformed to geographic coordinates.
! jn index for NH part of potential (kmlat down to ~kmlat0)
! jp index for NH pfrac (kmlat0 down to 1), the fraction of the dynamo
!
      real,dimension(kmlonp1,kmlat) :: phim
      integer :: jn,jp
!      
! electric field
      real,dimension(kmlonp1,kmlat) :: ed1dy,ed2dy
!
! Coefficients and RHS terms for PDE on geomagnetic grid:
! (formerly in coefm.h)
!
      real,dimension(kmlonp1,kmlat) ::                                  &
     &  zigm11,                                                         &! sigma11*cos(theta0)
     &  zigmc,                                                          &! sigmac
     &  zigm2,                                                          &! sigma2
     &  zigm22   ! sigma22/cos(theta0)
!
! rim(1)=id(1), rim(2)=id(2)/cos(theta0)
      real,dimension(kmlonp1,kmlat,2) :: rim 
      real,dimension(kmlonp1,kmlath)  :: rhs ! right-hand side 
!
! isolve = 0 -> original mud version 5.
! isolve = 1 -> muh hybrid solver (only as direct solver -- slow)
! isolve = 2 -> modified mudpack solver (mudmod)
! Default is isolve=2 for DYNAMO=2 (new dynamo), however under SGI/IRIX,
!   the hybrid solver is used because SGI crashes in mudmod if isolve==2. 
!
!nm20140804 #ifdef IRIX
!nm20140804      integer,parameter :: isolve = 1
!nm20140804 #else
!      integer,parameter :: isolve = 2 ! default is 2 for new dynamo
      integer,parameter :: isolve = 1 
!      integer,parameter :: isolve = 0
!nm20140804 #endif
!
! For dot products:
      real,parameter :: unitvm(kmlon)=1.
!
! am 10/04 define potential model : NONE or HEELIS - first test NONE
!       character(len=*),parameter :: potential_model='NONE'
!      character(len=*),parameter :: potential_model='HEELIS'
       character(len=*),parameter :: potential_model='weimer2005'
!
! Electric potential from heelis or weimer:
      real,dimension(kmlonp1,kmlat0) :: pfrac  ! NH fraction of potential
      real,dimension(kmlonp1,kmlat)  :: phihm  ! potential in magnetic
!      
! Work array (private to this module)
      real :: wkarray(-15:kmlon0+16,kmlat0)
!
! !REVISION HISTORY:
! 24.06.05  Astrid Maute: merge TIEGCM version1.8 into coupling with CTIP
!           for George Millward
! 05.02.04  Astrid Maute: include headers
! 
! comments from the original coupled code with CTIP
! am 10/04: adapte the electrodynamo from tiegcm1.6 for the
! coupling with a plasmasphere model
! -electro-dynamo is serial
! -high latitude potential is still in the code but shouldn't be 
! specified (set crit(1) and crit(2) to zero or small value
! crit(1), crit(2) colatitude for boundary of specified high lat. 
! potential (lat_m < crit(1) ) and transitions zone (crit(1)< lat_m < crit(2)
! - output eletric potential in magnetic coordinates
!          electric field in magnetic coordinates
!
! Module for electrodynamic dynamo (the "e" in tiegcm). 
!
! 5/02 btf: This is the "new" dynamo code, adapted by Astrid Maute
! and Art Richmond from the original ("old") dynamo in tgcm15,
! written by E. Cicely Ridley completed in about 1992. 
! 
! EOP
!
      end module dynamo_module
!-----------------------------------------------------------------------
