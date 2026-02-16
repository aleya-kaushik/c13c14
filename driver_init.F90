!---------------------------------------------------------------------
subroutine driver_init( )
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Sets initial drive data information
! Reads initialization information
!
! Expects standard driver output:
!   Year, DOY, Hour
!   Temperature (at 2m, in K)
!   Specific Humidity (kg/kg)
!   Wind Speed (m/s)
!   Surface Pressure (mb)
!   Surface Thermal Ratiaion Downwards (W/m2)
!   Surface Solar Radiation Downwards (W/m2)
!   Large Scale Precipitation (mm/s)
!   Convective Precipitation (mm/s)
!
! For time specification in the file:
!   Global: Mid-point of valid time range 
!   Site: Beginning of valid time range
!
! For format:
!   Global: Netcdf File
!   Site: Text File
!
!---------------------------------------------------------------------

use module_pparams, only: &
    secs_per_day, secs_per_hr
use module_sibconst, only: &
    single_pt, spinup, spinup_lnum
use module_io
use module_sib, only: sib
use module_time, only:  &
    day, hour, month, year, &
    curday_per_mon, sec_tot

use netcdf
implicit none

!driver data step setting variables
integer(i4) :: status, varid
integer(i4) :: yr, doy
real(r4), dimension(2) :: hr
character(len=256) :: record
logical :: exist

!misc variables
integer :: drnum

!-------------------------------------------

!...print message
if ((.not. spinup) .or. (spinup_lnum .eq. 1)) then
    print *, 'Setting Driver Data'
endif

!...set driver time step
if (single_pt) then
    write(unit=driver_filename,fmt='(a,i4.4,i2.2,a4)') &
          trim(driver_path), year, month, '.dat'
    inquire(file=trim(driver_filename),exist=exist)
    IF (.not. exist) THEN
        print*,'Missing Driver Data!'
        print*,'Please check dr_path in namel_sibdrv,'
        print*,'   as well as starting time.'
        STOP
    ENDIF

    driverid = 87
    open(unit=driverid, file=trim(driver_filename), &
         status='old', form='formatted', iostat=status)
    DO  ! read until not a comment
        read(driverid,'(a)',iostat=status) record
        IF (status .ne. 0) THEN
           print*,'Stopping due to driver data error!'
           stop
        ENDIF
        IF (record(1:1) .ne. '#') THEN
           exit
        ENDIF
    ENDDO

    read(unit=record,fmt=*) yr,doy,hr(1)
    read(driverid,'(a)',iostat=status) record
    read(unit=record,fmt=*) yr,doy,hr(2)
    close(driverid)
else
    write(driver_filename,fmt='(a,i4.4,i2.2,a3)') &
          trim(driver_path), year, month, '.nc'
    inquire(file=trim(driver_filename), exist=exist)
    IF (.not. exist) THEN
        print*,'Missing Global Driver Data!'
        print*,'Filename: ',trim(driver_filename)
        print*,'Please check dr_path in namel_sibdrv,'
        print*,'   as well as starting time.'
        STOP
    ENDIF

    status = nf90_open(trim(driver_filename), nf90_nowrite, &
                        driverid)
    status = nf90_inq_varid(driverid,'hour',varid)
    status = nf90_get_var(driverid,varid,hr) 
    IF (status .ne. nf90_noerr) THEN
        print*,'Error Reading Hour in Global Driver Data!'
        STOP
    ENDIF
    status = nf90_close(driverid)
endif

!.....check driver step
driver_step = int((hr(2)-hr(1))*3600)
if (driver_step == 0) then
    print*,'Time initialization for driver type not set.'
    print*,'Please check file and method in driver_init.F90'
    stop
elseif (mod (secs_per_day,driver_step) /= 0) then
    print*,'Driver Time Step Does Not Divide Evenly Into Day.'
    stop
endif

!...set time driver data variables
driver_perday = secs_per_day/driver_step
driver_permon = curday_per_mon * driver_perday
driver_switchf = .true.
driver_readf = .true.
driver_updatef = .true.

driver_year = year
driver_month = month
driver_day = day
driver_hour = hour
driver_recnum = (day-1)*driver_perday + 1

!...read in initial driver data
if ( single_pt ) then
    driver_seccur = sec_tot
    driver_secnext = sec_tot + driver_step

    do drnum=1, driver_recnum
       call driver_read_single(sib%g(1)%gprogt)
       driver_switchf = .false.
    enddo
    driver_recnum = driver_recnum + 1
else
    driver_seccur = sec_tot - (driver_step/2)
    driver_secnext = sec_tot + (driver_step/2)
    do drnum=1, driver_recnum
       call driver_read_global()
       driver_switchf = .false.
    enddo
endif

!...set variables
sib%g(:)%gprogt%cupr1 = sib%g(:)%gprogt%cupr2
sib%g(:)%gprogt%cupr  = sib%g(:)%gprogt%cupr2
sib%g(:)%gprogt%cuprt = sib%g(:)%gprogt%cupr * 0.001 !(m/s)
   
sib%g%gprogt%lspr1 = sib%g%gprogt%lspr2
sib%g%gprogt%lspr  = sib%g%gprogt%lspr2
sib%g%gprogt%lsprt = sib%g%gprogt%lspr * 0.001 !(m/s)
 
sib%g%gprogt%dlwbot1 = sib%g%gprogt%dlwbot2
sib%g%gprogt%dlwbot  = sib%g%gprogt%dlwbot2

sib%g%gprogt%ps1 = sib%g%gprogt%ps2
sib%g%gprogt%ps  = sib%g%gprogt%ps2

sib%g%gprogt%sh1 = sib%g%gprogt%sh2
sib%g%gprogt%sh  = sib%g%gprogt%sh2

sib%g%gprogt%spdm1 = sib%g%gprogt%spdm2
sib%g%gprogt%spdm  = sib%g%gprogt%spdm2

sib%g%gprogt%sw_dwn1 = sib%g%gprogt%sw_dwn2
sib%g%gprogt%sw_dwn  = sib%g%gprogt%sw_dwn2

sib%g%gprogt%tm1 = sib%g%gprogt%tm2
sib%g%gprogt%tm  = sib%g%gprogt%tm2

end subroutine driver_init
