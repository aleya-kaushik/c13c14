subroutine time_manager( )

!----------------------------------------------------------------------
! updates the time variables
!

use kinds
use module_pparams, only: &
    eqnx, secs_per_day, secs_per_hr, &
    days_per_year, days_per_month,  &
    doy1_month
use module_sibconst, only: &
    subcount, &
    leapyr_switch, &
    spinup, spinup_done, &
    spinup_numyrs, spinup_maxiter, &
    spinup_lnum, spinup_ynum, &
    tm5mr_switch
use module_io, only: &
    driver_seccur, driver_readf, &
    driver_updatef, &
    fire_seccur, fire_readf, &
    fire_updatef,   &
    requib_writef, &
    tm5_seccur, tm5_readf, &
    tm5_updatef
use module_time


implicit none


! local variables
integer(i4) :: x, dayperyr, yrnow
real(r8) :: dtdayfrac, hrtmp

!---------------------------------------------------------------------------
! TIME VARIABLES
!---------------------------------------------------------------------------

    !...update time variables
    new_year = .false.
    sec_tot = sec_tot + dtsib
    sec_year = sec_year + dtsib
    sec_month = sec_month + dtsib
    sec_day = sec_day+dtsib

    if (mod(sec_year,secs_per_day)==0) then
           doy = doy+1
           day = day+1
           sec_day = 0
    endif
    hrtmp=sec_day/secs_per_hr
    hour = floor(hrtmp)

    if ((sec_day == 0) .and. &
        (day > curday_per_mon)) then
        day = 1
        month = month + 1
        sec_month = 0

        if (month .gt. 12) then
           doy = 1
           month = 1
           year = year + 1
           sec_year = 0
           new_year = .true.
        endif


        curday_per_mon = days_per_month(month)
        if ((leapyr_switch) .and. &
             (month .eq. 2) .and. &
             (((mod(year,4) .eq. 0) .and. &
             (mod(year,100) .ne. 0)) .or. &
             (mod(year,400) .eq. 0))) then
            curday_per_mon = curday_per_mon + 1
        endif
    endif

    nmonth = month + 1
    if ( nmonth > 12 ) then
       nmonth = 1
       nyear = year + 1
    endif

    if (sec_month == 0) then
       new_month = .true.
    else
       new_month = .false.
    endif
    if (sec_month == secs_per_day) then
       new_month_day_delay = .true.
    else
       new_month_day_delay = .false.
    endif
 
    nday = day + 1
    if ( nday .eq. curday_per_mon) then
       curday_per_mon_day_ahead = days_per_month(nmonth)
       if ((leapyr_switch) .and. &
             (nmonth .eq. 2) .and. &
             (((mod(year,4) .eq. 0) .and. &
              (mod(year,100) .ne. 0)) .or. &
              (mod(year,400) .eq. 0))) then
            curday_per_mon_day_ahead = curday_per_mon_day_ahead + 1
        endif
    endif
    if ( nday > curday_per_mon ) nday=1

    if ( sec_day == 0 ) then
        new_day = .true.
    else
        new_day = .false.
    endif

    !----------------------
    ! UPDATE DAYS PER YEAR
    !-----------------------
    IF (leapyr_switch) THEN
       cureqnx=eqnx
       curday_per_year=days_per_year
       if (((mod(year,4) .eq. 0) .and. &
           (mod(year,100) .ne. 0)) .or. &
           (mod(year,400) .eq. 0)) then
           cureqnx=cureqnx+1
           curday_per_year=curday_per_year+1
       endif
   ENDIF

   IF (sec_year .EQ. secs_per_day) THEN
      curday_per_year_day_delay = curday_per_year
   ENDIF

    !----------------------
    ! UPDATE LOCAL STANDARD TIME HOUR FRACTION
    !-----------------------
    dtdayfrac = dble(dtsib) / dble(secs_per_day)
    call update_lst(subcount, dtdayfrac, &
                    dayfrac_lst, new_day_lst)

!---------------------------------------------------------------------------
! DRIVER DATA
!---------------------------------------------------------------------------
if ((driver_updatef) .and. &
    (sec_tot .eq. driver_seccur)) then
   driver_readf = .true.
endif

!---------------------------------------------------------------------------
! TM5 DATA  !ARACHO
!---------------------------------------------------------------------------
if ((tm5_updatef) .and. &
    (sec_tot .eq. tm5_seccur)) then
   tm5_readf = .true.
endif

!---------------------------------------------------------------------------
! FIRE DATA
!---------------------------------------------------------------------------
if ((fire_updatef) .and. &
    (sec_tot .eq. fire_seccur)) then
   fire_readf = .true.
endif

!---------------------------------------------------------------------------
! EQUILIBRIUM/SPINUP
!---------------------------------------------------------------------------
!...update spinup years ran
IF (sec_year .eq. 0) THEN 
   spinup_ynum = spinup_ynum + 1

   IF (year .GT. endyear) year = startyear
ENDIF

IF (sec_tot .ge. end_second) THEN

   IF (new_year) requib_writef = .true.
   IF (.not. spinup) THEN
       spinup_done = .true.  
   ELSE
       spinup_ynum = 1
       spinup_lnum = spinup_lnum + 1
       sec_tot=0
       start_second=0
       end_second=0

       yrnow=year
       DO x=1,spinup_numyrs
          dayperyr = days_per_year
          if ((leapyr_switch) .and. &
              (((mod(yrnow,4) .eq. 0) .and. &
              (mod(year,100) .ne. 0)) .or. &
              (mod(year,400) .eq. 0))) then
              dayperyr = dayperyr+1
          endif
          end_second=end_second + dayperyr*secs_per_day
          yrnow = yrnow + 1
          if (yrnow .gt. end_year) yrnow=start_year
        ENDDO

        IF (spinup_lnum > spinup_maxiter) THEN
            spinup_done = .true.
        ELSE
            call driver_init()
            if (tm5mr_switch) then
                call tm5_init()  !ARACHO
            endif
        ENDIF
   ENDIF  !spinup simulation
ENDIF  !end_second


end subroutine time_manager
