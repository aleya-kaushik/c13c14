subroutine fire_read_single(gprogt)

!****--------------------------------------------------------------------
!    This subroutines reads the fire emissions data 
!    for the current time step for single site (point) simulations.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.
!****--------------------------------------------------------------------

use kinds
use module_io
use module_pparams, only: mol_to_umol
use module_sib, only: gprog_type
use module_sibconst, only: &
    print_fire, print_stop
use module_time, only: &
     year, month, day, hour

implicit none

!...input variables
type(gprog_type), intent(inout) :: gprogt

!...local variables
integer(i4) :: yr, doy
real(r4) :: hr
real(r8) :: firecin, fireco2in

!...misc variables
integer(i4) :: status
logical :: opened, exist
character(len=256) :: record

!...storing previous data
gprogt%firec1 = gprogt%firec2
gprogt%fireco21 = gprogt%fireco22

!...switch files if needed
if (fire_switchf) then !switch and read
    inquire(unit=fireid, exist=exist,opened=opened)
    if (opened) close(fireid, iostat=status)

    write(unit=fire_filename,fmt='(a,i4.4,i2.2,a4)') &
        trim(fire_path), fire_year, fire_month, '.dat'
    inquire(file=trim(fire_filename), exist=exist)

    if (exist) then
       open(unit=fireid,file=trim(fire_filename), status='old', &
            form='formatted', iostat=status)
       if (status > 0) then
          print*,'!!!Error opening fire file!!!'
          stop
       endif
     else !fire file does not exist
         if (firefile_stop) then
             print*,''
             print*,'Stopping due to non-existent fire file: '
             print*,' ',trim(fire_filename)
             stop
         else
             fire_step = 0
             fireid = 0
             gprogt%firec2 = dzero
             gprogt%fireco22 = dzero

             RETURN
         endif
           
     endif 
endif !switch files

if (fire_recnum .gt. 0) then
   DO  !Read until not a comment
       read(fireid,'(a)',iostat=status) record
       IF (status /= 0) THEN
           print*,'Stopping due to error in single fire data.'
           stop
       ENDIF
       IF (record(1:1) .ne. '#') THEN
          exit
       ENDIF
   ENDDO

   read(unit=record,fmt=*) yr, doy, hr, &
        firecin, fireco2in
else
   firecin=dzero
   fireco2in=dzero
endif

!...set the new variables
gprogt%firec2 = firecin/mol_to_umol
gprogt%fireco22 = fireco2in/mol_to_umol

!...print out the new data if requested
if (print_fire) then
    print*,'     ------------------------------------'
    print*,'     !!!READING FIRE EMISSIONS!!!'

    if (fire_switchf) then
        print*,'   Opening fire file: '
        print*,'       ',trim(fire_filename)
    endif

    print*,'   REAL TIMES:'
    print*,'   Year    Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,i4)'), '   ', &
         year,'  ',month,'  ',day,'  ',hour
    print*,''
    print*,'   FIRE TIMES:'
    print*,'   Year    Mon   Day   Hour'
    print('(a,i5,a,i4,a,i4,a,f6.2)'), '   ', &
          fire_year,'  ',fire_month,'  ',fire_day, &
          '  ',fire_hour
    print*,'   FIRE RECNUM/TOTNUM: ', fire_recnum, fire_permon
    print('(a,e16.6)'),'     Fire C Loss (umol C/m2/s)    : ', &
              firecin
    print('(a,e16.6)'),'     Fire CO2 Resp (umol CO2/m2/s): ', &
              fireco2in
    print*,'     ------------------------------------'
 
    if (print_stop) stop
endif


end subroutine fire_read_single
