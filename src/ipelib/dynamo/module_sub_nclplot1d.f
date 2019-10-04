!
!t      module module_sub_ncplot1d
!t      
!t      PRIVATE
!t      PUBLIC :: ncplot
!t      contains
!--------------------------------------------------------------------------
!t      note: 1D version of sub_ncplot
      subroutine ncplot1d(noid,fname,labl,start1_in,count1_in,dim1,     &
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
      
      integer, intent(in) :: noid,start1_in(:),count1_in(:),nw
      integer, intent(in) :: dim1,nw_u
      real, intent(in) :: psi(:)
      character (LEN=*), intent(in) :: fname
      character (LEN=*), intent(in) :: labl
      character (LEN=*), intent(in) :: units

      integer :: istat,id
      character :: msgerr*19

      
! check if field is already defined
      istat = nf_inq_varid(noid,fname,id)
      print *,'id=', id,' istat=', istat
      if (istat /= NF_NOERR) then
	
        istat = nf_redef(noid)
        if (istat /= NF_NOERR) call check_err(istat,'redef ')

        istat = nf_def_var(noid,fname,NF_REAL,1,dim1,id)
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
      

!dbg20150615
      print *,'fname=', fname,' start1=',start1_in,' count1=',count1_in


! put values in new file
      msgerr = 'putvar '//fname
      istat = nf_put_vara_double(noid,id,start1_in,count1_in,psi)
      if (istat /= NF_NOERR) call check_err(istat,msgerr) 

!dbg20150615
      print *, 'sub-ncplot1D: noid',noid


      end subroutine ncplot1d
!--------------------------------------------------------------------------
!t      end module module_sub_ncplot
