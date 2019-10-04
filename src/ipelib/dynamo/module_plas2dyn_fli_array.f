      module module_plas2dyn_fli_array
        IMPLICIT NONE
        PRIVATE
        PUBLIC :: plas2dyn_fli_array
      contains
        subroutine plas2dyn_fli_array ( utime )
          USE module_precision
          USE module_eldyn,ONLY: plas_fli !t,Je_3d
!t      USE module_IPE_dimension,ONLY: NMP-->kmlon
          USE params_module,ONLY:kmlat,kmlon,kmlonp1
          USE cons_module,ONLY:idyn_save 
          use dynamo_module,only:zigm11,zigm22,zigmc,zigm2,rim
          use module_input_parameters,ONLY: stop_time
          IMPLICIT NONE
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      REAL (KIND=real_prec)    ::  dum(kmlon+1,kmlat) !dynamo FLI
      INTEGER (KIND=int_prec) :: jth,lp_dyn,lp_plas,ihem
      INTEGER (KIND=int_prec) :: ilat_dyn,ilon_dyn,mp,lp
      INTEGER (KIND=int_prec),parameter ::  lp_dyn_eq=47 !the lowest latitude index for FLI
!
      print *,'convert plas2dyn fli at utime=',utime,stop_time
      if (utime==stop_time) then
      open(4030,file='output_fli',form='formatted',status='new')
      write(4030,FMT='(I12)')utime
      endif


!convert from plas_fli to dynamo fli array: dum(kmlon+1,kmlat)
!(1) NH; (2) SH
!SMS$SERIAL(<plas_fli,IN>:default=ignore) BEGIN
      jth_loop: do jth=1,6
         print *, '!dbg20140407: jth=',jth,' mp',mp

         lp_dyn_loop: do lp_dyn=2,lp_dyn_eq !from SH toward eq
            lp_plas = idyn_save(lp_dyn)
            if (jth==1) print *, '!dbg20140407:lp_dyn=',lp_dyn          &
     &,' lp_plas',lp_plas

            ihem_loop: do ihem=1,2               
               if (ihem==1) then !NH
                  ilat_dyn = kmlat - lp_dyn +1
               else if ( ihem==2 ) then !SH
                  ilat_dyn = lp_dyn
               end if
               if (jth==1) print *, '!dbg20140407:ilat_dyn=',ilat_dyn   &
     &,ihem,lp_plas               

!mp4dyn is 180deg shifted from IPE
               mp_loop: do mp=1,kmlon

                  ilon_dyn = mp + (kmlon / 2)
                  if ( ilon_dyn > kmlon+1 ) ilon_dyn = ilon_dyn - kmlon
                  if (jth==1) print *, '!dbg20140407: mp',mp            &
     &,' ilon_dyn',ilon_dyn

                  dum(ilon_dyn,ilat_dyn) = plas_fli(ihem,lp_plas,mp,jth)
               end do mp_loop!: do mp=1,kmlon
            end do ihem_loop
         end do lp_dyn_loop
         
!perdiodic solution: kmlon+1<--mp=1
!dbg20140804         dum(kmlon+1,:)=dum(1,:)
         dum(1,:)=dum(kmlon+1,:)!dbg20140804

! magnetic poles: is this correct? all longitudes should be collapsed to one same value?
         do mp=1,kmlon+1
            dum(mp,     1)=dum(mp,       2) !South
            dum(mp,kmlat)=dum(mp,kmlat-1) !North
         end do

!no value at equator lp_dyn=49
!dbg20140407:lp_dyn=          47  lp_plas         170
 !dbg20140407:ilat_dyn=          51           1         170
 !dbg20140407:ilat_dyn=          47           2         170
 !dbg20140407:lp_dyn=          48  lp_plas           0
 !dbg20140407:ilat_dyn=          50           1           0
 !dbg20140407:ilat_dyn=          48           2           0
 !dbg20140407:lp_dyn=          49  lp_plas           0
 !dbg20140407:ilat_dyn=          49           1           0
 !dbg20140407:ilat_dyn=          49           2           0
         do mp=1,kmlon+1
            dum(mp,49) = ( dum(mp,47)+dum(mp,51) )*0.5 !magnetic eq
            dum(mp,48) = ( dum(mp,47)+dum(mp,49) )*0.5
            dum(mp,50) = ( dum(mp,49)+dum(mp,51) )*0.5
         end do

!dbg20140804: NOTE! this quick fix is made only for the original old grid!
! i must comment out these lines when using the new grid201207 etc!!!
!dbg20140407: bug in the original old grid at lp_dyn=39
 !dbg20140407:lp_dyn=          39  lp_plas           0
 !dbg20140407:ilat_dyn=          59           1           0
 !dbg20140407:ilat_dyn=          39           2           0
!         do mp=1,kmlon+1
!!SH
!            dum(mp,39)=( dum(mp,38)+dum(mp,40) )*0.5
!!NH
!            dum(mp,59)=( dum(mp,58)+dum(mp,60) )*0.5
!          end do

        if (utime==stop_time) then
        write(*,*) 'within pla2dyn',utime,stop_time
        write(4030,FMT='(20E12.4)')dum
        endif

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

!UNDERCONSTRUCTION!!!
!dbg20150608: copied from module_readin_ascii.f
! am 10/04 so far no value at the equator therefore set it
! but this shouldn't be in the code
      lp = kmlat/2+1
      do mp = 1,kmlonp1
	zigm11(mp,lp)= .125*(zigm11(mp,lp-1)+ zigm11(mp,lp+1))
	zigm22(mp,lp)= .125*(zigm22(mp,lp-1)+ zigm22(mp,lp+1))
	zigmc(mp,lp) = .125*(zigmc(mp,lp-1) + zigmc(mp,lp+1))
	zigm2(mp,lp) = .06 *(zigm2(mp,lp-1) + zigm2(mp,lp+1))
	rim(mp,lp,1) = .06 *(rim(mp,lp-1,1) + rim(mp,lp+1,1))
	rim(mp,lp,2) = .06 *(rim(mp,lp-1,2) + rim(mp,lp+1,2))
      enddo ! i = 1,kmlon
!
! am 10/04 change sign of K_(m lam)^D in the SH- that's what TIEGCM dynamo
! expects
      do lp = 1,(kmlat+1)/2
!        rim(:,lp,2) = rim(:,lp,2)*1.e8
        rim(:,lp,2) = -rim(:,lp,2)
      enddo
!dbg20150608: copy end

!SMS$SERIAL END
      !t           IF ( sw_3DJ==1 )  Je_3d(IN:IS,mp,1:2) = zero
      end subroutine plas2dyn_fli_array
      end module module_plas2dyn_fli_array
