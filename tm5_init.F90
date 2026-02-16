!---------------------------------------------------------------------
subroutine tm5_init( )
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Sets initial mixing ratio data information  
! Reads initialization information
!
! Expects standard tm5 output:
!   
!   Pressure (Pa)
!   COS mixing_ratio (molCOS/molAir)
!
! For format:
!   Netcdf File
!
! Written by Ara Cho (2020)
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

!tm5 data step setting variables
integer(i4) :: status, varid
integer(i4) :: yr !, doy
real(r4), dimension(2) :: hr
character(len=256) :: record
logical :: exist

!misc variables
integer :: drnum

!-------------------------------------------

!...print message
if ((.not. spinup) .or. (spinup_lnum .eq. 1)) then
    print *, 'Setting TM5 Mixing Ratio Data'
endif

!...set tm5 time step
write(tm5_filename,fmt='(a,a,i4.4,i2.2,a3)') &
      trim(tm5_path),'mix_TM5_', year, month, '.nc'
inquire(file=trim(tm5_filename), exist=exist)
IF (.not. exist) THEN
    print*,'Missing TM5 Data!'
    print*,'Filename: ',trim(tm5_filename)
    print*,'Please check tm5_path in namel_sibdrv,'
    print*,'   as well as starting time.'
    STOP
ENDIF

status = nf90_open(trim(tm5_filename), nf90_nowrite, &
                        tm5id)
status = nf90_inq_varid(tm5id,'hour',varid)
status = nf90_get_var(tm5id,varid,hr) 
IF (status .ne. nf90_noerr) THEN
    print*,'Error Reading Hour in TM5 Data!'
    STOP
ENDIF
status = nf90_close(tm5id)


!.....check tm5 step
tm5_step = int((hr(2)-hr(1))*3600)
if (tm5_step == 0) then
    print*,'Time initialization for tm5 type not set.'
    print*,'Please check file and method in tm5_init.F90'
    stop
elseif (mod (secs_per_day,tm5_step) /= 0) then
    print*,'TM5 Time Step Does Not Divide Evenly Into Day.'
    stop
endif

!...set time tm5 data variables
tm5_perday = secs_per_day/tm5_step
tm5_permon = curday_per_mon * tm5_perday
tm5_switchf = .true.
tm5_readf = .true.
tm5_updatef = .true.

tm5_year = year
tm5_month = month
tm5_day = day
tm5_hour = hour
tm5_recnum = (day-1)*tm5_perday + 1

!...read in initial tm5 data
tm5_seccur = sec_tot - (tm5_step/2)
tm5_secnext = sec_tot + (tm5_step/2)
do drnum=1, tm5_recnum
   call tm5_read(sib%g%gprogt)
   tm5_switchf = .false.
enddo


!...set variables
sib%g%gprogt%cosm_tm51 = sib%g%gprogt%cosm_tm52
sib%g%gprogt%cosm_tm5  = sib%g%gprogt%cosm_tm52

end subroutine tm5_init
