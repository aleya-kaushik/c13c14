#include "nc_util.h"

subroutine tm5_read(gprogt)

!****--------------------------------------------------------------------
!    This subroutines reads the forcing surface meteorological data
!    for the current time step for non-point simulations (regional/global).
!    If required, it closes the current month's data file and opens the
!    next month's data file.
!
!    Written by Ara Cho (2020)
!****--------------------------------------------------------------------

use kinds
use netcdf

use module_io
use module_sib, only: gprog_type
use module_sibconst, only: &
       nsib, subset, subcount !, &
       !print_tm5, print_stop
use module_time, only: &
       month, year, day

!...input variables
type(gprog_type), dimension(subcount), intent(inout) :: gprogt

!...local variables
real(r4), dimension(nsib) :: cosm_tm5   ! COSM
!real(r4), dimension(nsib) :: pressure_tm5   ! Surface Pressure (mb)

!...netcdf variables
character(len=20) dim_name
integer :: dimid, dimlen, status
integer :: ncyid,ncmid,nctdid
integer, dimension(2) :: mstart,mcount
integer :: psid_tm5, cosmid_tm5
logical :: tm5_existf

!...misc variables
integer :: i
integer :: xyear
integer(byte) :: xmonth,xday

!***--------------------------------------------------------------------
   !...storing previous time steps data

   do i=1,subcount
       gprogt(i)%cosm_tm51    = gprogt(i)%cosm_tm52
   !    gprogt(i)%pressure_tm51    = gprogt(i)%pressure_tm52  !PRES
   enddo
   !CHANGE ----------------_***
   ! switch files if needed
   tm5_existf = .false.
   if ( tm5_switchf ) then
       status = nf90_close( tm5id )

       write(tm5_filename,fmt='(a,a,i4.4,i2.2,a3)') trim(tm5_path), &
             'mix_TM5_', tm5_year, tm5_month,'.nc'
       inquire(file=trim(tm5_filename), exist=tm5_existf)
       if (tm5_existf) then
           STATUS = nf90_open( trim(tm5_filename), nf90_nowrite, tm5id )
           IF (status .ne. nf90_noerr) THEN
               print*,'Error Opening TM5 File: '
               print*,'  ',trim(tm5_filename)
               STOP
           ENDIF
           print*,'   Reading tm5 file: ',trim(tm5_filename)
        else
            tm5id = -1
        endif
    endif

    ! read data if tm5 file is open
    if (tm5id > 0) then

        !...check nsib
        CHECK( nf90_inq_dimid(tm5id, trim(nsibname), dimid ) )
        CHECK( nf90_inquire_dimension(tm5id, dimid, dim_name,dimlen ) )
        if ( dimlen /= nsib ) then
             print*,'File and specified nsib are different.  Stopping.'
             print*,'  File nsib: ',dimlen,'Sim nsib: ', nsib
             stop
        endif

        !...check time values in tm5 data file
        ENSURE_VAR( tm5id,'year', ncyid )
        ENSURE_VAR( tm5id,'month',ncmid )
        ENSURE_VAR( tm5id,'day',  nctdid )
        mstart(1) = tm5_recnum

        CHECK( nf90_get_var( tm5id, ncyid, xyear, mstart(1:1) ) )
        if (xyear .ne. tm5_year) then
           print*,'Year in file does not match simulation. Stopping.'
           print*,'   File: ',xyear,' Sim: ',year
           stop !ARACHO; need to block?
        endif

        CHECK( nf90_get_var( tm5id, ncmid, xmonth, mstart(1:1) ) )
        if (xmonth .ne. tm5_month) then
           print*,'Month in file does not match simulation. Stopping.'
           print*,'   File: ',xmonth,' Sim: ',tm5_month
           print*, ' Driver:',driver_month
           stop
        endif

        CHECK( nf90_get_var( tm5id, nctdid, xday, mstart(1:1) ) )
        if (xday .ne. tm5_day) then
           print*,'Day in file does not match simulation. Stopping.'
           print*,'   File: ',xday,' Sim: ', tm5_day
           stop
        endif

        !...get veriable id's
        !PRES-ARA
        !ENSURE_VAR( tm5id, trim(psname_tm5), psid_tm5 ) ! Surface Pressure
        ENSURE_VAR( tm5id, trim(cosmname_tm5), cosmid_tm5 )  ! COS mixing ratio from TM5

        !...get data
        mstart=(/1,tm5_recnum/); mcount=(/nsib,1/)
        !PRES-ARA
        !STATUS = nf90_get_var(tm5id, psid_tm5, pressure_tm5, mstart, mcount) !Pressure
        !IF (status .ne. nf90_noerr) THEN
        !    print*,'Error Getting Pressure'
        !    STOP
        !ENDIF

        STATUS = nf90_get_var(tm5id, cosmid_tm5,cosm_tm5 , mstart, mcount) !Temperature
        IF (status .ne. nf90_noerr) THEN
            print*,'Error Getting cosm'
            STOP
        ENDIF

        !...pull out points in subdomain
         do i=1, subcount
         ! PRES-ARA
         !  gprogt(i)%pressure_tm52 = pressure_tm5(subset(i))
           gprogt(i)%cosm_tm52 = cosm_tm5(subset(i))
         enddo

        !...print out the new data if requested
        !if (print_tm5) then
        !    print*,'New tm5 data read (yr/mon/day/hr/sec/recnum): '
        !    print*,'   ',tm5_year,tm5_month,tm5_day, &
        !                 tm5_hour,tm5_secnext,tm5_recnum
        !    print*,'------------------------------------------------------'
        !    print*,'Extrema of new input data'
        ! ! PRES-ARA
        ! !   print*, minval(gprogt%pressure_tm52      ),maxval(gprogt%pressure_tm52 ),' Pressure (Pa)'
        !    print*, minval(gprogt%cosm_tm52      ),maxval(gprogt%cosm_tm52 ),' cos mixing ratio'
        !    print*,'-----------------------------------------------------'

        !    if (print_stop) stop
        !endif

   else  !tm5 data file not found
         print*,'!!!TM5 Data File Not Found.'
   endif

end subroutine tm5_read
