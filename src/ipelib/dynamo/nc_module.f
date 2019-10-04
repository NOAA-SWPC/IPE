!
      module nc_module
      
      use params_module,only:                                           &
     &  kmlonp1,                                                        &! nmlon+1
     &  kmlat   ! number of geomagnetic grid latitudes
!      
      integer :: noid,                                                  &! id of output netcdf-file
     &   dim3(3),                                                       &! id of  defining 3dim-arrays on irregular grid
     &   dim1 	           ! id of  defining 1dim-arrays on irregular grid

!      
      integer,parameter ::                                              &
     &     count3(3)= (/kmlonp1,kmlat,1/),                              &! only for read in 3D fields 
     &     count1(1)= (/              1/),                              &! only for read in 1D fields 
     &     start3_out(3)= (/1,1,1/),                                    &! only for put out 3D fields 
     &     start1_out(1)= (/    1/),                                    &! only for put out 3D fields 
     &     start3(3)= (/1,1,2/)             ! only for read in 3D fields  
         
      character(len=*),parameter ::                                     &
     &   nc_path =                                                      &
     & 'potent_lres.nc'  !nm021207:
!nm021207:     & '/home/naomi/astrid/tiegcm1.8_dynamo_lres/potent_lres.nc'  !nm082306
!nm082306:   & '/iapetus/i/naomi/tmp/RCM_bnd/potent_rcmbnd.nc' !nm041106
!nm041106:   & '/iapetus/i/naomi/tmp/tiegcm1.8_dynamo_lres/potent.nc' !nm040706
!nm040706:   & '/suncat/e/maute/plasmasphere/potent.nc'

      integer :: n_time=0      
!--------------------------------------------------------------------------
      end module nc_module
