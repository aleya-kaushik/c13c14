subroutine time_init( )

!--------------
!Initializes the time variables
!-------------

use module_oparams, only: &
     clim_len, seas_len, seas_len_precip
use module_pparams, only: &
     eqnx, secs_per_day, &
     days_per_year, days_per_month, &
     doy1_month
use module_sibconst
use module_io
use module_time

implicit none

! local variables
integer(i4) :: x, yrstorun
integer(i4) :: dayperyr, doy1now, yrnow
integer(i4) :: starttimet, endtimet
real(r8) :: dayfrac
integer(i4) :: ppos
integer, parameter :: nlen = 256
character(len=nlen) :: str

!----------------------------------------------
print*,''
print*,'Setting Time Variables'

    dtisib = done/dtsib
    steps_per_day = secs_per_day / dtsib

    if (leapyr_switch) then
       if (((mod(endyear,4) .eq. 0) .and. &
          (mod(endyear,100) .ne. 0)) .or. &
          (mod(endyear,400) .eq. 0)) then
          if (endtime .eq. 365) then
              endtime = endtime + 1
          endif
       else
          if (endtime .eq. 366) then
              endtime = endtime - 1
          endif
       endif
    endif !leapyear

    if (spinup) then
       endtime = days_per_year
       if ((leapyr_switch) .and. &
            (((mod(endyear,4) .eq. 0) .and. &
            (mod(endyear,100) .ne. 0)) .or. &
            (mod(endyear,400) .eq. 0))) then
          endtime = endtime + 1
       endif
    endif
    starttimet = (starttime-1)*secs_per_day
    endtimet = endtime*secs_per_day

    !...make sure endtime doesn't occur before starttime
    if ( endyear == startyear ) then
      if ((starttimet >= endtimet) )  &
          stop 'simulation ends before it starts, check starttime,'  &
              //' endtime, startyear, and endyear'
    endif
    
    !...make sure endyear doesn't occur before startyear
    if ( endyear < startyear )  &
        stop 'simulation ends before it starts, check startyear and endyear'
    
    !...make sure number of seconds in simulation is evenly divisible by dtsib
    if ( mod( endtimet-starttimet, dtsib) /= 0 )  &
        stop 'dtsib does not divide evenly into the total simulation'

    !...make sure dtsib divides evenly into a day
    if (  mod( secs_per_day, dtsib ) /= 0 )  &
        stop 'dtsib does not divide evenly into a day'
        
    !...make sure restart_dtsib is evenly divisible by dtsib
    if ( restart_dtsib > 0 .and. mod( restart_dtsib, dtsib ) /= 0 )  &
        stop 'Restart Error: restart_dtsib not evenly divisible by dtsib'
    
    !...make sure qp_dtsib is evenly divisible by dtsib
    if ( qp_dtsib > 0 .and. mod( qp_dtsib, dtsib ) /= 0 )  &
        stop 'QP Output Error: qp_dtsib not evenly divisible by dtsib'
    
    !...make sure pbp_dtsib is evenly divisible by dtsib
    if ( pbp_dtsib > 0  .and. mod( pbp_dtsib, dtsib) /= 0) &
         stop 'PBP Output Error: pbp_dtsib not evenly divisible by dtsib'
    
    !...make sure hr_dtsib is evenly divisible by dtsib
    if ( hr_dtsib > 0  .and. mod( hr_dtsib, dtsib) /= 0) &
         stop 'HR Output Error: hr_dtsib not evenly divisible by dtsib'

    !...set initial values
    init_year = startyear
    init_doy = starttime
    init_second = starttimet

    new_day = .true.
    do x = 1, 12
       doy1now = doy1_month(x)
       if (((MOD(startyear,4) .EQ. 0) .AND. &
           (MOD(startyear,100) .ne. 0)) .or. &
           (MOD(startyear,400) .eq. 0)) then
          IF (x .gt. 2) THEN
             doy1now = doy1now + 1
           ENDIF
        ENDIF
        if ( init_doy >= doy1now ) then
            init_month = x
            init_day = init_doy - doy1now + 1
        endif
    enddo
    
    if (init_doy == doy1_month(init_month)) then
       new_month = .true.
    else
       new_month = .false.
    endif
    if (init_day == 2) then
       new_month_day_delay = .true.
    else 
       new_month_day_delay = .false.
    endif

    start_year = init_year
    start_month = init_month
    start_day = init_day
    start_doy = init_doy
    start_second = init_second
    
    if (new_day .and. new_month) then
       if (start_month .eq. 1) then
           new_year = .true.
       endif
    endif

    !...set ending values
    end_year = endyear
    end_doy  = endtime
    
    dayperyr = days_per_year
    if ((leapyr_switch) .and. &
        (((mod(end_year,4) .eq. 0) .and. &
        (mod(end_year,100) .ne. 0)) .or. &
        (mod(end_year,400) .eq. 0))) then
         dayperyr = dayperyr + 1
    endif

    dayfrac = dtsib/dble(secs_per_day)
    end_dtime = endyear + (end_doy - dayfrac)/dble(dayperyr)

    end_second = 0
    if (spinup) then
         yrnow=start_year
         yrstorun = spinup_numyrs
         do x=1,yrstorun
            dayperyr = days_per_year
             if ((leapyr_switch) .and. &
                 (((mod(yrnow,4) .eq. 0) .and. &
                 (mod(yrnow,100) .ne. 0)) .or. &
                 (mod(yrnow,400) .eq. 0))) then
                dayperyr = dayperyr + 1
             endif
             end_second = end_second + dayperyr * secs_per_day
             yrnow = yrnow + 1
             if (yrnow .gt. end_year) yrnow=start_year
         enddo
     else 
         yrstorun = (end_year - start_year)
         do x = start_year, start_year+yrstorun-1
             dayperyr = days_per_year
             if ((leapyr_switch) .and.  &
                 (((mod(x,4) .eq. 0) .and. &
                 (mod(x,100) .ne. 0)) .or. &
                 (mod(x,400) .eq. 0))) then
                  dayperyr = dayperyr + 1
              endif 
              end_second = end_second + dayperyr * secs_per_day
         enddo
         end_second = end_second + endtimet
    endif

    !...set current times to get ready for simulation
    year = start_year
    month = start_month
    doy = start_doy
    day = start_day
    hour = 0
    sec_day = 0
    sec_month = 0 
    sec_year = start_second
    sec_tot = start_second

    nmonth = month
    nyear = year   
    nday = day

    cureqnx = eqnx
    curday_per_year  = days_per_year
    curday_per_mon=days_per_month(month)
    curday_per_mon_day_ahead = curday_per_mon

    if ((leapyr_switch) .and. &
        (((mod(year,4) .eq. 0) .and. &
        (mod(year,100) .ne. 0)) .or. &
        (mod(year,400) .eq. 0))) then
       cureqnx = cureqnx + 1
       curday_per_year = curday_per_year + 1
 
       if (month .eq. 2) then
         curday_per_mon = curday_per_mon + 1
       endif
    endif
    curday_per_year_day_delay = curday_per_year

    !...initialize local standard time information
    allocate(dayfrac_lst(subcount), new_day_lst(subcount))
    dayfrac = dble(sec_day) / dble(secs_per_day)
    call calc_lst(dayfrac, subcount, sublonsib, dayfrac_lst)
    !new_day_lst(:) = .true.

    !...initialize sun position based on doy
    call init_solar_dec(doy, cureqnx, curday_per_year, lonearth)
    call solar_dec(doy, cureqnx, lonearth, sin_dec, cos_dec, tan_dec)
    
    !...set maximum daylength
    allocate(daylenmax(subcount))
    call max_day_length(subcount, sublatsib, daylenmax)
    
    !...initialize spinup variables
    if (spinup .and. (.not. spinup_continue)) then
        spinup_ynum = 1
        spinup_lnum = 1
        if (spinup_maxiter .eq. -1) then
           spinup_maxiter = 1
           spinup_numyrs = spinup_numyrs*((endyear-startyear)+1)
        endif
    elseif ((.not. spinup) .and. (.not. spinup_continue)) then
        spinup_ynum = 1
        spinup_lnum = spinup_maxiter
        spinup_numyrs = (endyear-startyear)+1
    endif

    !...initialize spinup variables when continuing
    if (spinup .and. (spinup_continue)) then
        spinup_ynum = 1
        !spinup_lnum = 1
        str=trim(ic_file)
        ppos=scan(str,".",BACK=.true.)
        !print*,"ppos: ", ppos
        read(str(ppos-1:ppos),*) spinup_prevnum
        spinup_lnum = 1+spinup_prevnum
        if (spinup_maxiter .eq. -1) then
           spinup_maxiter = 1
           spinup_numyrs = spinup_numyrs*((endyear-startyear)+1)
        endif
    elseif ((.not. spinup) .and. (spinup_continue)) then
        spinup_ynum = 1
        spinup_lnum = spinup_maxiter
        spinup_numyrs = (endyear-startyear)+1
    endif

    !...set weights
    wt_clim = 1./MAX(done,real(clim_len)*real(steps_per_day))
    wt_daily = 1./real(steps_per_day)
    wt_seas = 1./MAX(done,real(seas_len)*real(steps_per_day))
    wt_seas_precip = 1./MAX(done,real(seas_len_precip)*real(steps_per_day))
    
    wtd_clim = 1./MAX(done,real(clim_len))
    wtd_seas = 1./MAX(done,real(seas_len))

end subroutine time_init

