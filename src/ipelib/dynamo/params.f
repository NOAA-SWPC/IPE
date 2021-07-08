!
!t #include "defs.h"
      
      module params_module
!
! Magnetic grid:
      integer,parameter ::                                              &
!    &  kmlat =97,                                                      &! KMLAT,         ! number of magnetic latitudes
     &  kmlat =193,                                                     &! KMLAT,         ! number of magnetic latitudes
!    &  kmlon =80,                                                      &! KMLON,         ! number of magnetic longitudes
!    &  kmlon =160,                                                     &! KMLON,         ! number of magnetic longitudes
     &  kmlon =320,                                                     &! KMLON,         ! number of magnetic longitudes
     &  kmlonp1=kmlon+1,                                                &
     &  kmlatp1=kmlat+1,                                                &
     &  kmlath=(kmlat+1)/2     ! index to magnetic equator
!
      end module params_module
