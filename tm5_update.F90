!--------------------------------------------------------
subroutine tm5_update(gprogt)
!--------------------------------------------------------
!Updates the time/control variables for the tm5 data.
!Written by Ara Cho (2020)
!--------------------------------------------------------

use kinds
use module_pparams, only:  &
    secs_per_day, secs_per_hr, &
    hrs_per_day
use module_io
use module_sib, only: gprog_type
use module_sibconst, only: single_pt, subcount
use module_time, only: &
    startyear, endyear, end_second, &
    curday_per_mon_day_ahead

implicit none

! input variables
type(gprog_type), dimension(subcount), intent(inout) :: gprogt

! local variables
integer(long) :: local_secahead

!-----------------------------------------------------

if (tm5_readf) then
    !.....read in tm5 data needed
    call tm5_read(gprogt)

    tm5_recnum = tm5_recnum + 1
    tm5_seccur = tm5_seccur + tm5_step
    tm5_secnext = tm5_secnext + tm5_step
    local_secahead = tm5_secnext + tm5_step

    tm5_hour = tm5_hour + floor(tm5_step/dble(secs_per_hr))
    if (tm5_hour == hrs_per_day) then
        tm5_day = tm5_day + 1
        tm5_hour = 0
    endif

    tm5_readf = .false.
    tm5_switchf = .false.
    if ( tm5_recnum >  tm5_permon ) then
        !....switch tm5 data for new month
        if (local_secahead .lt. end_second) then
            tm5_recnum = 1
            tm5_month = tm5_month + 1
            tm5_day = 1
            if (tm5_month > 12) then
                tm5_month=1
                tm5_year = tm5_year+1
                if (tm5_year .GT. endyear) tm5_year=startyear
           endif
           tm5_permon = curday_per_mon_day_ahead * tm5_perday
           tm5_switchf = .true.
        else
           tm5_updatef = .false.
        endif
    endif
endif !tm5_readf

end subroutine tm5_update
