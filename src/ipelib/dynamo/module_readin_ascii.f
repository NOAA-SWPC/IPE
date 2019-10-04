!test022208: v2: how much degree can the polar cap conductanc affect the potential solution??? 
!
      module module_readin_ascii      
!      
      PRIVATE
      PUBLIC :: readin_ascii
      contains
!-----------------------------------------------------------------------
      subroutine readin_ascii
      use params_module,only:                                           &
     &  kmlonp1,                                                        &! kmlon+1
     &  kmlat   ! number of geomagnetic grid latitudes

      use read_module,only:path_ascii

! read in data file with integrals from George Millward
!
      use dynamo_module,only:                                           &
     &  zigm11,                                                         &! zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
     &  zigm22,                                                         &! zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
     &  zigmc,                                                          &! zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
     &  zigm2,                                                          &! zigm2  is int[sigma_h] ds,	      i.e. Sigma_h
     &  rim      ! K_(m phi)^D/abs(sin Im),  K_(m lam)^D
      
!t      use module_sub_ncplot,ONLY:ncplot
!t      use cons_module,only: secs
      implicit none
!
! Local 
      integer:: lp,mp,i,is,ie,itotal,j
      integer, parameter :: iunit= 179,ounit=4001
      character :: fname*10,labl*56,units*12
      real,dimension(kmlonp1,kmlat) ::                                  &
     &  part1,part2,part3  ! for comparison with k_mlam terms
      integer,parameter:: sw_part1=0
      integer :: jth,utime
      real,dimension(kmlonp1,kmlat)::dum
!      
      open(iunit,file=path_ascii,status='OLD',err=99)
      write(6,*) 'open ascii file ',path_ascii
      read(unit=iunit,FMT="(I12)",err=199)utime
      print *,'utime=',utime
!t      secs = REAL(utime)
!t      print *,'secs=',secs
!
      jth_loop: do jth=1,6
      write(6,*) 'read in jth=',jth
        read(unit=iunit,FMT="(20E12.4)")dum





        if (jth==1) then 
           zigm11=dum
        else if (jth==2) then 
           zigm22=dum
        else if (jth==3) then 
           zigm2 =dum
        else if (jth==4) then 
           zigmc =dum
        else if (jth==5) then 
           rim(:,:,1)=dum
        else if (jth==6) then 
           rim(:,:,2)=dum
        end if
      end do jth_loop


      close(iunit,err=100)  
	
! am 10/04 so far no value at the equator therefore set it
! but this shouldn't be in the code
      j = kmlat/2+1
      do i = 1,kmlonp1
	zigm11(i,j)= .125*(zigm11(i,j-1)+ zigm11(i,j+1))
	zigm22(i,j)= .125*(zigm22(i,j-1)+ zigm22(i,j+1))
	zigmc(i,j) = .125*(zigmc(i,j-1) + zigmc(i,j+1))
	zigm2(i,j) = .06 *(zigm2(i,j-1) + zigm2(i,j+1))
	rim(i,j,1) = .06 *(rim(i,j-1,1) + rim(i,j+1,1))
	rim(i,j,2) = .06 *(rim(i,j-1,2) + rim(i,j+1,2))
      enddo ! i = 1,kmlon
      
! am 10/04 change sign of K_(m lam)^D in the SH- that's what TIEGCM dynamo
! expects
      do lp = 1,(kmlat+1)/2
!        rim(:,lp,2) = rim(:,lp,2)*1.e8
        rim(:,lp,2) = -rim(:,lp,2)
      enddo
!      
! am 10/04 set high latitude values since these are no calculated
! in the plasmasphere model region of Heelis pattern
! value itself should not matter since potential is prescribed  
!t      do lp = 1,14
!t        zigm11(:,lp) = 0.01
!t        zigm22(:,lp) = 0.01
!t        zigm11(:,kmlat+1-lp) = 0.01
!t        zigm22(:,kmlat+1-lp) = 0.01
!t      enddo
      
!dbg20140801(2)
      print *,'output dyn fli at utime=',utime
      write(unit=4002,FMT='(I12)')utime
      write(unit=4002,FMT='(20E12.4)')zigm11

      return
!      
 2435 format(10e12.4)
 2436 format(e12.4)
   99 write(6,*) 'Error opening input file',path_ascii
      STOP
  199 write(6,*) 'Error reading input file'
      STOP
  100 write(6,*) 'Error closing input file'
!      
      end subroutine readin_ascii 
!--------------------------------------------------------------------------  
      end module module_readin_ascii
