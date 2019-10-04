!
      module module_nc_create      
!
      PRIVATE
      PUBLIC :: nc_create
      contains
!--------------------------------------------------------------------------
      subroutine nc_create
!
      use nc_module
!     nc-output: define plot dimensions & attributes
      use cons_module,only:                                             &
     &  xlatm_deg,                                                      &! magnetic latitudes (deg)
     &  xlonm_deg   ! magnetic longitudes (deg)
      use params_module,only:                                           &
     &  kmlonp1,                                                        &! kmlon+1
     &  kmlat   ! number of geomagnetic grid latitudes
      use module_check_err,ONLY: check_err
!
      implicit none
      
!#ifdef SUN
!#include "/opt/local/include/netcdf.inc"
!#elif AIX
!#include "netcdf.inc"
!#elif LINUX
!#include "netcdf.inc"
      include "netcdf.inc"
!#endif
!
      integer :: istat,idtime,idlon,idlat,itime,ilon,ilat
!      
      write(6,*) 'output nc-file created:',nc_path

! create output file
      istat = nf_create(nc_path,NF_CLOBBER,noid) 
      if (istat /= NF_NOERR) call check_err(istat,'create nc-file ')

! define dimensions
      istat = nf_def_dim(noid,'time',NF_UNLIMITED,idtime)
      if (istat /= NF_NOERR) call check_err(istat,'defdim time ')

      istat = nf_def_dim(noid,'mlon',kmlonp1,idlon)
      if (istat /= NF_NOERR) call check_err(istat,'defdim mlon ')

      istat = nf_def_dim(noid,'mlat',kmlat,idlat)
      if (istat /= NF_NOERR) call check_err(istat,'defdim mlat ')

! ids of 3dim array on irregular grid
      dim3(1) = idlon
      dim3(2) = idlat
      dim3(3) = idtime

! define variables      
      istat = nf_def_var(noid,'time',NF_REAL,1,idtime,itime)
      if (istat /= NF_NOERR) call check_err(istat,'defvar time ')
      
      istat = nf_def_var(noid,'mlon',NF_REAL,1,idlon,ilon)
      if (istat /= NF_NOERR) call check_err(istat,'defvar mlon ')

      istat = nf_def_var(noid,'mlat',NF_REAL,1,idlat,ilat)
      if (istat /= NF_NOERR) call check_err(istat,'defvar mlat ')

! end definition part      
      istat = nf_enddef(noid)
      if (istat /= NF_NOERR) call check_err(istat,'enddef ') 

! put latitude/longitude in new file
      istat = nf_put_var_double(noid,ilat,xlatm_deg)
      if (istat /= NF_NOERR) call check_err(istat,'putvar mlat ')
      
      istat = nf_put_var_double(noid,ilon,xlonm_deg)
      if (istat /= NF_NOERR) call check_err(istat,'putvar mlon ')
!
      end subroutine nc_create
!--------------------------------------------------------------------------
      end module module_nc_create

