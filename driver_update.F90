!--------------------------------------------------------
subroutine driver_update()
!--------------------------------------------------------
!Updates the time/control variables for the driver data.
!--------------------------------------------------------

use kinds
use module_pparams, only:  &
    secs_per_day, secs_per_hr, &
    hrs_per_day
use module_io
use module_sib, only: sib
use module_sibconst, only: single_pt, subcount
use module_time, only: &
    startyear, endyear, end_second, &
    curday_per_mon_day_ahead

implicit none

! local variables
integer(long) :: local_secahead

!-----------------------------------------------------

if (driver_readf) then
    !.....read in driver data needed
    if (single_pt) then
       call driver_read_single(sib%g(1)%gprogt)
    else
       call driver_read_global()
    endif

    driver_recnum = driver_recnum + 1
    driver_seccur = driver_seccur + driver_step
    driver_secnext = driver_secnext + driver_step
    local_secahead = driver_secnext + driver_step

    driver_hour = driver_hour + floor(driver_step/dble(secs_per_hr))
    if (driver_hour == hrs_per_day) then
        driver_day = driver_day + 1
        driver_hour = 0
    endif

    driver_readf = .false.
    driver_switchf = .false.
    if ( driver_recnum >  driver_permon ) then
        !....switch driver data for new month
        if (local_secahead .lt. end_second) then
            driver_recnum = 1
            driver_month = driver_month + 1
            driver_day = 1
            if (driver_month > 12) then
                driver_month=1
                driver_year = driver_year+1
                if (driver_year .GT. endyear) driver_year=startyear
           endif
           driver_permon = curday_per_mon_day_ahead * driver_perday
           driver_switchf = .true.
        else
           driver_updatef = .false.
        endif
    endif
endif !driver_readf

end subroutine driver_update
