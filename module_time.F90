module module_time

!------------------------------------------------------------
! 
!   SiB4 Time Variable Module
!
!------------------------------------------------------------

use kinds
implicit none

    ! specified times from namel_sibdrv
    integer(i4) :: &
        starttime,   & ! day of year to start simulation
        startyear,   & ! year to start simulation
        endtime,     & ! day of year to end simulation
        endyear        ! year to end simulation

    integer(i4) :: dtsib    !  # seconds in simulation time step
    real(r8)    :: dtisib   ! inverse model timestep (1/s)
    integer(i4) :: steps_per_day   ! number of time steps per day

    ! time constants
    integer(i4) :: init_year     ! initial year of simulation
    integer(i4) :: init_month    ! initial month of simulation
    integer(i4) :: init_day      ! initial day of month of simulation
    integer(i4) :: init_doy      ! initial day of year of simulation
    integer(i4) :: init_second   ! initial second of simulation

    integer(i4) :: start_year    ! year of restart file
    integer(i4) :: start_month   ! month of restart file
    integer(i4) :: start_day     ! day of month of restart file
    integer(i4) :: start_doy     ! day of year of restart file
    integer(i4) :: start_second  ! second of restart file
    
    integer(i4)   :: end_year    ! last year of simulation
    integer(i4)   :: end_doy     ! last day of year of simulation
    integer(long) :: end_second  ! last second of simulation
    real(r8)      :: end_dtime   ! last timestep of simulation,
                                 !    in fractional year
    
    integer(i4) :: total_years   ! total number of years in simulation
    integer(i4) :: total_months  ! total number of months in simulation
    integer(i4) :: total_days    ! total number of days in simulation

    ! time variables
    integer(i4) :: year     ! current year in the simulation
    integer(i4) :: month    ! current month in the simulation
    integer(i4) :: hour     ! current hour of day in the simulation
    integer(i4) :: day           ! current day of the current month 
    integer(i4) :: doy           ! current day of current year
    integer(i4) :: sec_day       ! current second in the current day
    integer(i4) :: sec_month     ! current second in the current month
    integer(i4) :: sec_year      ! current second in the current year
    integer(long) :: sec_tot       ! current second in the whole simulation
    
    integer(i4) :: nyear         ! year of next month of simulation
    integer(i4) :: nmonth        ! next month of simulation
    integer(i4) :: nday          ! next day of simulation

    integer(i4) :: curday_per_year  ! number of days in current year
    integer(i4) :: curday_per_year_day_delay  ! number of days in year,
                                              !    switched on January 2nd
    integer(i4) :: curday_per_mon   ! number of days in current month
    integer(i4) :: curday_per_mon_day_ahead    ! number of days in month, 
                                               !   switched on last day
                                               !   of previous month
    integer(i4) :: cureqnx   ! day of vernal equinox in current year

    ! Flags
    logical :: new_day               ! new day ?
    logical :: new_month             ! new month ?
    logical :: new_month_day_delay   ! second day of month?
    logical :: new_year              ! new year ?

    ! Local Time Information (subcount)
    real(r8), dimension(:), allocatable :: &
          dayfrac_lst ! fraction of day in local time
    logical, dimension(:), allocatable ::  &
          new_day_lst  ! new day in local time?
 
    ! Day Length Information
    real(r4), dimension(:), allocatable :: &  !(subcount)
           daylenmax

    ! Weights
    real(r8) :: wt_clim  !time-step weight for climatological running-mean
    real(r8) :: wt_daily !time-step weight for daily running-mean
    real(r8) :: wt_seas  !time-step weight for seasonal running-mean
    real(r8) :: wt_seas_precip !time-step weight for precip seas mean

    real(r8) :: wtd_clim !daily weight for climatological running-mean
    real(r8) :: wtd_seas !daily weight for seasonal running-mean

end module module_time
