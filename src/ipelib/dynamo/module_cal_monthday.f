      MODULE module_cal_monthday
      implicit none
      PRIVATE
      PUBLIC :: cal_monthday
      CONTAINS
!---------------------
      SUBROUTINE cal_monthday (year,sum,month,day)
      USE module_input_parameters,ONLY:sw_debug
        implicit none

        integer,intent(in) :: year
        integer,intent(in) :: sum !=NDAY
        integer,intent(out) :: month
        integer,intent(out) :: day

!most of years
        integer, parameter :: nmonth=12
        integer, parameter, dimension(nmonth) :: day_of_month0=(/       &
     & 31,28,31,30,31,30                                                &
     &,31,31,30,31,30,31 /)

!leap years (29days in Feb): 2000,2004,2008,2012
        integer, parameter, dimension(nmonth) :: day_of_month1=(/       &
     & 31,29,31,30,31,30                                                &             
     &,31,31,30,31,30,31 /)

        integer :: m,nday4m
!-----------------


      if ( sw_debug ) then
         print *, "input YEAR=",year
         print *, "day of year=",sum
      end if


        day=sum
        month_loop1: do m=1,nmonth

!leap years
        if( MOD(year, 4)==0 ) then
           nday4m = day_of_month1(m)
        else 
           nday4m = day_of_month0(m)
        endif

        if ( day <= nday4m ) then
           month = m
           EXIT month_loop1
        else                    ! if ( day > nday4m ) then
           day = day - nday4m
        end if
     
      enddo month_loop1
!     print *, "output DATE: year=",year," month=",month," day=",day   

      if ( sw_debug )       print *,'subroutine cal_monthday finished!'
      END SUBROUTINE cal_monthday
      END MODULE module_cal_monthday
