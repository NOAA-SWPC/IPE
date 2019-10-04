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
!Aug2011: the original code was provided from Fei Wu from WAM version
!nm20110906: modified to implement to IPE
!copied from /.../testall.f
!ylonm(1:nmlon=180)
!ylatm(1:nmlat=90)
!      program ts_efield
      MODULE module_eldyn
      USE module_precision
!----------------------
!c idea
!      subroutine idea_geteb(im,ix,dayno,utsec,f107,kp,maglat,maglon,
!     &essa,ee1,ee2)
      USE efield
!c     use date_def
!c     use physcons, pi => con_pi
      IMPLICIT NONE
!      integer, intent(in) :: im  ! number of data points in efield 
!      integer, intent(in) :: ix  ! max data points in efield
!      integer :: dayno=254  ! calender day
!      real :: utsec=0.0  ! second
!      real :: f107=70.  ! 
!      real, intent(in) :: maglat(im)  ! magnetic latitude (rad)
!      real, intent(in) :: maglon(im)  ! magnetic longitude (rad)
!      real, intent(in) :: essa(im)  ! degree
!      real, intent(out)   :: ee1(im)    ! electric field x direction mV/m
!      real, intent(out)   :: ee2(im)    ! electric field y direction mV/m
!c     character*(*), intent(in) ::   dir    ! directory located coef files
!c local
!      integer i,k,iref,jref
!       real :: utsec_last
!      real utsec_last,dx,dy,aa,bb,maglond,maglatd,
!     &ed11(0:nmlon,0:nmlat),ed22(0:nmlon,0:nmlat),ylatm1(0:nmlat),
!     &ylonm1(0:nmlon),pi
!      real, parameter :: pi=3.141592653
!      logical first
!c
!      data first/.true./
!      data utsec_last/-1./
!      save first,utsec_last,ed11,ed22,ylatm1,ylonm1

      REAL   (KIND=real_prec),DIMENSION(0:nmlat  ),PUBLIC :: theta90_rad

!SMS$DISTRIBUTE(dh,2) BEGIN
      INTEGER(KIND=int_prec ),allocatable,public :: j0      (:,:) !1:NH; 2:SH
      INTEGER(KIND=int_prec ),allocatable,public :: j1      (:,:) !1:NH; 2:SH
      REAL   (KIND=real_prec),allocatable,public :: coslam_m(:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,2,3) BEGIN
      REAL   (KIND=real_prec),allocatable,public :: Ed1_90  (:,:,:)
      REAL   (KIND=real_prec),allocatable,public :: Ed2_90  (:,:,:)
!nm20131025:...ideally, this should be located for sw_eldyn=0, should not be mixed with waccm efield!!!
      REAL   (KIND=real_prec),allocatable,public :: plas_fli (:,:,:,:) !(2,NLP,NMP,6)
!t      REAL   (KIND=real_prec),allocatable,public ::sigma_ped(:,:,:)
!t      REAL   (KIND=real_prec),allocatable,public ::sigma_hall(:,:,:)
!SMS$DISTRIBUTE END

!nm20121003:subroutine init_eldyn, eldyn are separated into module_sub_eldyn.f

      END MODULE module_eldyn
