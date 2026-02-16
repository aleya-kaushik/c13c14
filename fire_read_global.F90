#include "nc_util.h"

subroutine fire_read_global()

!****--------------------------------------------------------------------
!    This subroutines reads the fire emissions data for the current
!    time step for global (group/multi-point) simulations.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.
!****--------------------------------------------------------------------

use kinds
use netcdf

use module_io
use module_pparams, only: mol_to_umol
use module_sib, only: sib
use module_sibconst, only: &
   nsib, subset, subcount, &
   print_fire, print_stop
use module_time, only: &
   month, year, day, hour

implicit none

!...data variables
real(r8), dimension(nsib) :: firecin
real(r8), dimension(nsib) :: fireco2in

!...netcdf variables
integer :: status
integer :: dimid, dimlen
integer :: ncyid, ncmid, ncdid
integer :: fcid, fco2id
character(len=20) :: dim_name
integer, dimension(2) :: mstart, mcount

!...local variables
integer :: i
integer :: xyear, xmonth, xday
logical :: fire_existf

!****--------------------------------------------------------------------
! storing previous time steps data
do i=1,subcount
    sib%g(i)%gprogt%firec1 = sib%g(i)%gprogt%firec2
    sib%g(i)%gprogt%fireco21 = sib%g(i)%gprogt%fireco22
enddo

! switch files if needed
fire_existf = .false.
if (fire_switchf) then
   status = nf90_close(fireid)

   write(fire_filename,fmt='(a,i4.4,i2.2,a3)') trim(fire_path), &
         fire_year, fire_month, '.nc'
   inquire(file=trim(fire_filename), exist=fire_existf)
   if (fire_existf) then
      STATUS = nf90_open(trim(fire_filename), nf90_nowrite, fireid)
   else
      if (firefile_stop) then
          print*,''
          print*,'Missing/Non-Existant Fire Emissions File!'
          print*,'File: ',trim(fire_filename)
          print*,'Please check fr_path in namel_sibdrv.'
          print*,''
          STOP
        else
          fire_step = 0
          fireid = 0
          sib%g(:)%gprogt%firec2 = dzero
          sib%g(:)%gprogt%fireco22 = dzero

          RETURN
        endif
    endif
endif !fire_switchf

! read data 
!...check nsib
CHECK(nf90_inq_dimid(fireid, trim(fnsibname), dimid))
CHECK(nf90_inquire_dimension(fireid, dimid, dim_name, dimlen))
if (dimlen /= nsib) then
   print*,'Mismatching Fire File!'
   print*,' File nsib: ',dimlen,' Sim nsib: ',nsib
   STOP
endif

!...check time values
ENSURE_VAR(fireid,'year',ncyid)
ENSURE_VAR(fireid,'month',ncmid)
ENSURE_VAR(fireid,'day',ncdid)
mstart(1) = fire_recnum

CHECK(nf90_get_var(fireid, ncyid, xyear, mstart(1:1)))
if (xyear .ne. fire_year) then
   print*,'Fire year in file does not match simulation.  Stopping.'
   print*,'   File: ',xyear,' Sim: ',year
   stop
endif

CHECK(nf90_get_var(fireid, ncmid, xmonth, mstart(1:1)))
if (xmonth .ne. fire_month) then
   print*,'Fire month in file does not match simulation.  Stopping.'
   print*,'   File: ',xmonth,' Sim: ',month
   stop
endif

CHECK(nf90_get_var(fireid, ncdid, xday, mstart(1:1)))
if (xday .ne. fire_day) then
   print*,'Fire day in file does not match simulation. Stopping.'
   print*,'  File: ',xday,' Sim: ', day
   stop
endif

!...get variable id's
ENSURE_VAR( fireid, trim(firecname), fcid)   !firec
ENSURE_VAR( fireid, trim(fireco2name), fco2id) !fireco2

!...get data
mstart=(/1,fire_recnum/); mcount=(/nsib,1/)
STATUS = nf90_get_var(fireid, fcid, firecin, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting Fire C Emissions!'
   stop
ENDIF

status = nf90_get_var(fireid, fco2id, fireco2in, mstart, mcount)
IF (status .ne. nf90_noerr) THEN
   print*,'Error Getting Fire CO2 Emissions!'
   stop
ENDIF

!...pull out points in subdomain
do i=1, subcount
   sib%g(i)%gprogt%firec2 = firecin(subset(i))/mol_to_umol
   sib%g(i)%gprogt%fireco22 = fireco2in(subset(i))/mol_to_umol
enddo

!...print out the new data if requested
if (print_fire) then
    print*,'     ------------------------------------'
    print*,'     !!!READING FIRE EMISSIONS!!!'

    if (fire_switchf) then
       print*,'     Opening fire file: '
       print*,'       ',trim(fire_filename)
    endif
    print*,'   REAL TIMES:' 
    print*,'   Year   Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,i4)'), '   ', &
         year,'  ',month,'  ',day,'  ',hour
    print*,''
    print*,'   FIRE TIMES:'
    print*,'   Year    Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,f6.2)'), '   ', &
          fire_year,'  ',fire_month,'  ',fire_day, &
          '  ',fire_hour
    print*,'   FIRE RECNUM/TOTNUM: ', fire_recnum, fire_permon
    print('(a,2e16.6)'),'     Min/Max Fire C Loss (umol C/m2/s)    : ', &
          minval(firecin), maxval(firecin)
    print('(a,2e16.6)'),'     Min/Max Fire CO2 Resp (umol CO2/m2/s): ', &
          minval(fireco2in), maxval(fireco2in)
    print*,'     ------------------------------------'
 
    if (print_stop) stop
endif

end subroutine fire_read_global
