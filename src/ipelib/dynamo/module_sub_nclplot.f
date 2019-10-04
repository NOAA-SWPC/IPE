!
!t      module module_sub_ncplot
!t      
!t      PRIVATE
!t      PUBLIC :: ncplot
!t      contains
!--------------------------------------------------------------------------
      subroutine ncplot(noid,fname,labl,start3,count3,dim3,             &
     &                  psi,nw,units,nw_u)

!      use nc_module
      use module_check_err,ONLY:check_err
!     nc-output: define plot variable & attributes

      implicit none
!     
!#ifdef SUN
!#include "/opt/local/include/netcdf.inc"
!#else
!#include "netcdf.inc"
      include "netcdf.inc"
!#endif
      
      integer, intent(in) :: noid,start3(3),count3(3),nw
      integer, intent(in) :: dim3(3),nw_u
      real, intent(in) :: psi
      character, intent(in) :: fname*10,labl*56,units*12

      integer :: istat,id
      character :: msgerr*19

      
! check if field is already defined
      istat = nf_inq_varid(noid,fname,id)
      if (istat /= NF_NOERR) then
	
        istat = nf_redef(noid)
        if (istat /= NF_NOERR) call check_err(istat,'redef ')

        istat = nf_def_var(noid,fname,NF_REAL,3,dim3,id)
        msgerr = 'defvar '//fname
        if (istat /= NF_NOERR) call check_err(istat,msgerr)

        istat = nf_put_att_text(noid,id,'long_name',nw,labl)
        msgerr = 'put att_'//fname//' name'
        if (istat /= NF_NOERR) call check_err(istat,msgerr)
	
        istat = nf_put_att_text(noid,id,'units',nw_u,units)
        msgerr = 'put att_'//fname//' units'
        if (istat /= NF_NOERR) call check_err(istat,msgerr)

        istat = nf_enddef(noid)
        if (istat /= NF_NOERR) call check_err(istat,'enddef ')
	 
      endif
      
! put values in new file
      msgerr = 'putvar '//fname
      istat = nf_put_vara_double(noid,id,start3,count3,psi)
      if (istat /= NF_NOERR) call check_err(istat,msgerr) 

      end subroutine ncplot
!--------------------------------------------------------------------------
!t      end module module_sub_ncplot
