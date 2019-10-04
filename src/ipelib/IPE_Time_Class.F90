MODULE IPE_Time_Class

USE IPE_Precision

IMPLICIT NONE

  TYPE IPE_Time
    REAL(prec) :: utime
    REAL(prec) :: elapsed_sec
    INTEGER    :: year, month, day, day_of_year, hour, minute
    INTEGER(8) :: day_number

    CONTAINS

    PROCEDURE :: Build => Build_IPE_Time

    PROCEDURE :: Set_Date  => Set_Date_IPE_Time
    PROCEDURE :: Get_Date  => Get_Date_IPE_Time
    PROCEDURE :: Set_UTime => Set_UTime_IPE_Time
    PROCEDURE :: Get_UTime => Get_UTime_IPE_Time
    PROCEDURE :: Set_HourMinute => Set_HourMinute_IPE_Time

    PROCEDURE :: Update => Update_IPE_Time
    PROCEDURE :: Increment => Increment_IPE_Time
    PROCEDURE :: Calculate_Hour_and_Minute => Calculate_Hour_and_Minute_IPE_Time
    PROCEDURE :: Calculate_UTime           => Calculate_UTime_IPE_Time

    PROCEDURE :: DateStamp => DateStamp_IPE_Time
    PROCEDURE :: Time_From_DateStamp => Time_From_DateStamp_IPE_Time

    PROCEDURE :: Calculate_DayNumber    => Calculate_DayNumber_IPE_Time
    PROCEDURE :: Calculate_YearMonthDay => Calculate_YearMonthDay_IPE_Time
    PROCEDURE :: Calculate_Date_Difference
    PROCEDURE :: Calculate_Day_Difference

  END TYPE IPE_Time



CONTAINS

  SUBROUTINE Build_IPE_Time( time_tracker, time_stamp, rank_id )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    CHARACTER(12), INTENT(in)        :: time_stamp
    INTEGER, INTENT(in)              :: rank_id

      CALL time_tracker % Time_From_DateStamp( time_stamp )

      CALL time_tracker % Calculate_DayNumber( )
      CALL time_tracker % Calculate_UTime( )
      ! The day of the year is calculated as the difference between the current
      ! date and January 1 of the same year.
      time_tracker % day_of_year = time_tracker % Calculate_Day_Difference( time_tracker % year, 1, 1 ) + 1

      time_tracker % elapsed_sec = 0.0_prec

      if(rank_id.eq.0) then
      write(6,*) ''
      write(6,*) '***********************************************'
      write(6,*) '*                                             *'
      write(6,*) '*    IPE : Initial DateTime : ',time_stamp,'    *'
      write(6,*) '*                                             *'
      write(6,*) '***********************************************'
      write(6,*) ''
      endif

  END SUBROUTINE Build_IPE_Time

  SUBROUTINE Set_Date_IPE_Time( time_tracker, year, month, day )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    INTEGER, INTENT(in)              :: year, month, day

      time_tracker % year  = year
      time_tracker % month = month
      time_tracker % day   = day

      CALL time_tracker % Calculate_DayNumber( )

  END SUBROUTINE Set_Date_IPE_Time

  SUBROUTINE Get_Date_IPE_Time( time_tracker, year, month, day )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(in) :: time_tracker
    INTEGER, INTENT(out)          :: year, month, day

      year  = time_tracker % year
      month = time_tracker % month
      day   = time_tracker % day

  END SUBROUTINE Get_Date_IPE_Time

  SUBROUTINE Set_UTime_IPE_Time( time_tracker, utime )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    REAL(prec), INTENT(in)           :: utime

      time_tracker % utime = utime

  END SUBROUTINE Set_UTime_IPE_Time

  SUBROUTINE Set_HourMinute_IPE_Time( time_tracker, hour, minute )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    INTEGER, INTENT(in)              :: hour, minute

      time_tracker % hour = hour
      time_tracker % minute = minute

      CALL time_tracker % Calculate_UTime( )

  END SUBROUTINE Set_HourMinute_IPE_Time

  SUBROUTINE Get_UTime_IPE_Time( time_tracker, utime )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(in) :: time_tracker
    REAL(prec), INTENT(out)       :: utime

      utime = time_tracker % utime

  END SUBROUTINE Get_UTime_IPE_Time

  SUBROUTINE Update_IPE_Time( time_tracker, utime )
   ! Uses the utime to advance the minute, hour, day, month, year
   IMPLICIT NONE
   CLASS( IPE_Time ), INTENT(inout) :: time_tracker
   REAL(prec), INTENT(in), OPTIONAL :: utime
   ! Local
   INTEGER :: additional_hours, additional_days

      IF( PRESENT( utime ) )THEN
        time_tracker % utime = utime
      ENDIF

      CALL time_tracker % Calculate_Hour_and_Minute( )

      IF( time_tracker % minute >= 60 )THEN

        additional_hours = MOD( time_tracker % minute, 59 )
        time_tracker % minute = time_tracker % minute - additional_hours*60
        time_tracker % hour   = time_tracker % hour + additional_hours

      ENDIF


      IF( time_tracker % hour >= 24 )THEN

        additional_days = MOD( time_tracker % hour, 23 )
        time_tracker % hour = time_tracker % hour - additional_days*24
        time_tracker % day_number  = time_tracker % day_number + additional_days
        time_tracker % utime = time_tracker % utime - 86400.0_prec

      ENDIF

      CALL time_tracker % Calculate_UTime( )
      CALL time_tracker % Calculate_YearMonthDay( )

  END SUBROUTINE Update_IPE_Time

  SUBROUTINE Increment_IPE_Time( time_tracker, delta_utime )
   ! Uses the utime to advance the minute, hour, day, month, year
   IMPLICIT NONE
   CLASS( IPE_Time ), INTENT(inout) :: time_tracker
   REAL(prec), INTENT(in)           :: delta_utime
   ! Local
   INTEGER :: additional_hours, additional_days
   INTEGER :: delta_day, delta_hour, delta_minute

      time_tracker % elapsed_sec = time_tracker % elapsed_sec + delta_utime
      time_tracker % utime = time_tracker % utime + delta_utime

      delta_day    = INT(delta_utime/86400.0_prec)
      delta_hour   = INT((delta_utime-REAL(delta_day,prec)*24.0_prec)/3600.0_prec)
      delta_minute = INT((delta_utime-REAL(delta_hour,prec)*3600.0_prec-REAL(delta_day,prec)*86400.0_prec)/60.0_prec)

      time_tracker % day_number = time_tracker % day_number + delta_day
      time_tracker % hour       = time_tracker % hour + delta_hour
      time_tracker % minute     = time_tracker % minute + delta_minute

      IF( time_tracker % minute >= 60 )THEN

        additional_hours      = time_tracker % minute/60
        time_tracker % minute = time_tracker % minute - additional_hours*60
        time_tracker % hour   = time_tracker % hour + additional_hours

      ENDIF


      IF( time_tracker % hour >= 24 )THEN

        additional_days = time_tracker % hour/24
        time_tracker % hour = time_tracker % hour - additional_days*24
        time_tracker % day_number  = time_tracker % day_number + additional_days
        time_tracker % utime = time_tracker % utime - 86400.0_prec

      ENDIF

      CALL time_tracker % Calculate_YearMonthDay( )
      time_tracker % day_of_year = time_tracker % Calculate_Day_Difference( time_tracker % year, 1, 1 ) + 1

  END SUBROUTINE Increment_IPE_Time

  SUBROUTINE Calculate_Hour_and_Minute_IPE_Time( time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    ! Local
    INTEGER :: additional_hours

      time_tracker % hour   = INT( time_tracker % utime/3600.0_prec )
      time_tracker % minute = INT( (time_tracker % utime - REAL(time_tracker % hour,prec)*3600.0_prec)/60.0_prec )

  END SUBROUTINE Calculate_Hour_and_Minute_IPE_Time

  SUBROUTINE Calculate_UTime_IPE_Time( time_tracker )
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker

      time_tracker % utime = REAL( time_tracker % hour*3600 + time_tracker % minute*60 )


  END SUBROUTINE Calculate_UTime_IPE_Time

  FUNCTION DateStamp_IPE_Time( time_tracker ) RESULT( timestamp )
    CLASS( IPE_Time ) :: time_tracker
    CHARACTER(12)     :: timestamp

      WRITE( timestamp, '(I4,I2.2,I2.2,I2.2,I2.2)' ) time_tracker % year, &
                                                     time_tracker % month, &
                                                     time_tracker % day, &
                                                     time_tracker % hour, &
                                                     time_tracker % minute



  END FUNCTION DateStamp_IPE_Time

  SUBROUTINE Time_From_DateStamp_IPE_Time( time_tracker, timestamp )
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    CHARACTER(12), INTENT(in)        :: timestamp

      READ( timestamp, '(I4,I2.2,I2.2,I2.2,I2.2)' ) time_tracker % year, &
                                                    time_tracker % month, &
                                                    time_tracker % day, &
                                                    time_tracker % hour, &
                                                    time_tracker % minute



  END SUBROUTINE Time_From_DateStamp_IPE_Time

  SUBROUTINE Calculate_DayNumber_IPE_Time( time_tracker )
  ! Algorithm adapted from
  ! https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    ! Local
    INTEGER :: month, year

      ! Adjust the calendar to start with march 1 at the "beginning of the year"
      month = MOD( time_tracker % month + 9, 12 )
      year  = time_tracker % year - month/10

      time_tracker % day_number = 365*year + &
                                  year/4 - &
                                  year/100 + &
                                  year/400 + &
                                 (month*306 + 5)/10 + &
                                 (time_tracker % day - 1 )

  END SUBROUTINE Calculate_DayNumber_IPE_Time

  SUBROUTINE Calculate_YearMonthDay_IPE_Time( time_tracker )
  ! Algorithm adapted from
  ! https://alcor.concordia.ca/~gpkatch/gdate-algorithm.html
    IMPLICIT NONE
    CLASS( IPE_Time ), INTENT(inout) :: time_tracker
    ! Local
    INTEGER(8) :: month, year, day, adj_day, mi

      year = (10000*time_tracker % day_number)/3652425
      adj_day = time_tracker % day_number - (365*year + year/4 - year/100 + year/400)

      IF( adj_day < 0 )THEN
       year = year - 1
       adj_day = time_tracker % day_number - (365*year + year/4 - year/100 + year/400)
      ENDIF

      mi    = (100*adj_day + 52)/3060

      month = MOD( (mi + 2), 12 ) + 1
      year  = year + (mi + 2)/12
      day   = adj_day - (mi*306 + 5)/10 + 1

      time_tracker % year  = year
      time_tracker % month = month
      time_tracker % day   = day

  END SUBROUTINE Calculate_YearMonthDay_IPE_Time

  FUNCTION Calculate_Date_Difference( time_tracker, year, month, day, hour, minute ) RESULT( diff_minutes )
  ! Calculates the difference between the date in time_tracker and the given
  ! date. The result is reported in units of minutes
    IMPLICIT NONE
    CLASS( IPE_Time ) :: time_tracker
    INTEGER           :: year, month, day, hour, minute
    REAL(prec)        :: diff_minutes
    ! Local
    INTEGER :: month_adj, year_adj, day_number, day_diff
    REAL(prec) :: ut

      ! Adjust the calendar to start with march 1 at the "beginning of the year"
      month_adj = MOD( month + 9, 12 )
      year_adj  = year - month_adj/10

      day_number = 365*year_adj + &
                   year_adj/4 - &
                   year_adj/100 + &
                   year_adj/400 + &
                   (month_adj*306 + 5)/10 + &
                   (day - 1 )

      CALL time_tracker % Calculate_DayNumber( )

      day_diff = time_tracker % day_number - day_number

      ut = REAL(hour*60 + minute,prec) ! [ minutes ]


      diff_minutes = REAL( day_diff*86400, prec ) + time_tracker % utime/60.0_prec - ut

  END FUNCTION Calculate_Date_Difference

  FUNCTION Calculate_Day_Difference( time_tracker, year, month, day ) RESULT( day_diff )
  ! Calculates the difference between the date in time_tracker and the given
  ! date. The result is reported in units of minutes
    IMPLICIT NONE
    CLASS( IPE_Time ) :: time_tracker
    INTEGER           :: year, month, day
    REAL(prec)        :: day_diff
    ! Local
    INTEGER :: month_adj, year_adj, day_number
    REAL(prec) :: ut

      ! Adjust the calendar to start with march 1 at the "beginning of the year"
      month_adj = MOD( month + 9, 12 )
      year_adj  = year - month_adj/10

      day_number = 365*year_adj + &
                   year_adj/4 - &
                   year_adj/100 + &
                   year_adj/400 + &
                   (month_adj*306 + 5)/10 + &
                   (day - 1 )

      CALL time_tracker % Calculate_DayNumber( )

      day_diff = time_tracker % day_number - day_number



  END FUNCTION Calculate_Day_Difference

END MODULE IPE_Time_Class
