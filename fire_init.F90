!---------------------------------------------------------------------
subroutine fire_init( )
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Sets initial fire information
!
! Expects standard driver output:
!   Year, DOY, Hour
!   Emis_C (micromoles C/m2/s)
!   Emis_CO2 (micromoles C/m2/s)
!
! For time specification in the file:
!   Mid-point of valid time range 
!
! For format:
!   Global: Netcdf File
!   Site: Text File
!
!---------------------------------------------------------------------

use kinds
use module_io
use module_pparams, only: &
    days_per_month, secs_per_day, mol_to_umol
use module_sib, only: sib
use module_sibconst
use module_time

use netcdf
implicit none

!parameters
logical, parameter :: fire_required = .true.

!fire data step setting variables
integer(i4) :: status, varid
integer(i4) :: fyr, fdoy, fmon, fmonref, dpermon
integer(i4) :: fire_recnumstart
integer(long) :: fsec_tot
real(r4), dimension(2) :: fhr
character(len=256) :: record
logical :: exist

!local variables
integer :: frnum
character(len=iofnlen) :: fire_tempname

!----------------------------------------------------------------------

if ((spinup) .or. (.not. fire_switch)) then
   !Only use fires for production runs
   fire_step = izero
   RETURN
else
   print*,'Setting Fire Data'
endif

!...find first fire file
fire_year = izero
fire_month = izero
do fyr=init_year, end_year
   if (fyr .eq. init_year) then
      fmonref = init_month
      if (fmonref .gt. 1) then
          dpermon = days_per_month(fmonref-1)
          if ((mod(fyr,4) .eq. 0) .and. (fmonref .eq. 2)) then
              dpermon = dpermon + 1
           endif
      else
          dpermon = izero
      endif
      fsec_tot = dpermon*secs_per_day
   else
      fmonref = 1
   endif
         
   do fmon=fmonref,12
      if (single_pt) then
          write(fire_tempname,fmt='(a,i4.4,i2.2,a4)') &
               trim(fire_path), fyr, fmon, '.dat'
      else
         write(fire_tempname,fmt='(a,i4.4,i2.2,a3)') &
              trim(fire_path), fyr, fmon, '.nc'
      endif
      inquire(file=trim(fire_tempname), exist=exist)

      dpermon = days_per_month(fmon)
      if ((mod(fyr,4) .eq. 0) .and. (fmon .eq. 2)) then
         dpermon = dpermon + 1
      endif
      
      if ((exist) .AND. (fire_year .eq. izero)) then
         fire_year = fyr
         fire_month = fmon
         fire_seccur = fsec_tot
         fire_filename = fire_tempname
      else
         fsec_tot = fsec_tot + dpermon*secs_per_day
      endif
   enddo
enddo


IF (fire_year .eq. izero) THEN
    if (fire_required) THEN
        print*,''
        print*,'Missing/Non-Existant Fire Emissions File!'
        print*,'Directory: ',trim(fire_path)
        print*,'Please check fr_path in namel_sibdrv.'
        print*,''
        STOP
     else
        print*,'  Missing Fire Data, Using Zero Emissions.'
        fire_step = 0
        RETURN
     endif
ENDIF


!...set fire timestep
if (single_pt) then !site simulation
    fireid = 43
    open(unit=fireid, file=trim(fire_filename), &
         status='old', form='formatted', iostat=status)
    DO   ! read until not a comment
         read(fireid,'(a)',iostat=status) record
         IF (status .ne. 0) THEN
            print*,'Stopping due to fire data error!'
            stop
         ENDIF
         IF (record(1:1) .ne. '#') THEN
            exit
         ENDIF
     ENDDO

     read(unit=record,fmt=*) fyr, fdoy, fhr(1)
     read(fireid,'(a)',iostat=status) record
     read(unit=record,fmt=*) fyr, fdoy, fhr(2)
     close(fireid)

else !global simulation
     status = nf90_open(trim(fire_filename), nf90_nowrite, &
                  fireid)
     status = nf90_inq_varid(fireid,'hour',varid)
     status = nf90_get_var(fireid,varid,fhr) 
     IF (status .ne. nf90_noerr) THEN
         print*,'Error Reading Hour in Global Fire Data!'
         STOP
     ENDIF
     status = nf90_close(fireid)
endif


!...set and check fire step
fire_step = int((fhr(2)-fhr(1))*3600)
if (fire_step .le. 0) then
   print*,'Missing Fire Time Step'
   STOP
elseif (mod(secs_per_day,fire_step) /= 0) then
   print*,'Fire Time Step Does Not Divide Evenly Into Day.'
   STOP
endif


!...set fire variables
fire_perday = secs_per_day/fire_step
fire_permon = curday_per_mon * fire_perday
fire_hour = fhr(1)

if ((fire_year .eq. init_year) .and. &
    (fire_month .eq. init_month)) then
     fire_recnumstart = (day-1)*fire_perday + 1
     fire_day = day
else
     fire_recnumstart = 1
     fire_day = 1
endif
fire_secnext = fire_seccur + fire_step
fire_updatef = .true.
fire_readf = .false.


!...read in initial fire data
fire_switchf = .true.
fire_recnum = 1
do frnum=1, fire_recnumstart
    if (single_pt) then
       call fire_read_single(sib%g(1)%gprogt)
    else
       call fire_read_global()
    endif
    fire_hour = fire_hour + real(fire_step/3600.)
    fire_recnum = fire_recnum + 1
    fire_seccur = fire_seccur + fire_step
enddo

sib%g(:)%gprogt%firec1 = sib%g(:)%gprogt%firec2
sib%g(:)%gprogt%firec  = sib%g(:)%gprogt%firec2
sib%g(:)%gprogt%fireco21 = sib%g(:)%gprogt%fireco22
sib%g(:)%gprogt%fireco2  = sib%g(:)%gprogt%fireco22
fire_switchf = .false.

end subroutine fire_init


