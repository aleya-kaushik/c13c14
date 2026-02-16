#include "nc_util.h"

subroutine driver_read_global()

!****--------------------------------------------------------------------
!    This subroutines reads the forcing surface meteorological data 
!    for the current time step for non-point simulations (regional/global).
!    If required, it closes the current month's data file and opens the 
!    next month's data file.
!****--------------------------------------------------------------------

use kinds
use netcdf

use module_io
use module_sib, only: sib
use module_sibconst, only: &
       nsib, subset, subcount, &
       print_driver, print_stop
use module_time, only: &
       month, year

!...local variables
real(r4), dimension(nsib) :: ps   ! Surface Pressure (mb)
real(r4), dimension(nsib) :: tm   ! Temperature (K)
real(r4), dimension(nsib) :: sh   ! Specific Humidity (kg/kg)
real(r4), dimension(nsib) :: spdm ! Wind speed (m/s)
real(r4), dimension(nsib) :: lspr ! Large Scale Precipitation (mm/s)
real(r4), dimension(nsib) :: cupr ! Convective Precipitation (mm/s)
real(r4), dimension(nsib) :: swd  ! Surface solar radiation downwards (W/m2)           
real(r4), dimension(nsib) :: lwd  ! Surface thermal radiation downwards (W/m2)

!...netcdf variables
character(len=20) dim_name
integer :: dimid, dimlen, status
integer :: ncyid,ncmid,nctdid
integer, dimension(2) :: mstart,mcount
integer :: spdmid, tmid, swdid, &
            lwdid, shid, psid, lsprid, cuprid
logical :: driver_existf

!...misc variables
integer :: i
integer :: xyear
integer(byte) :: xmonth,xday

!***--------------------------------------------------------------------
   !...storing previous time steps data
   do i=1,subcount
       sib%g(i)%gprogt%ps1    = sib%g(i)%gprogt%ps2
       sib%g(i)%gprogt%tm1    = sib%g(i)%gprogt%tm2
       sib%g(i)%gprogt%sh1    = sib%g(i)%gprogt%sh2
       sib%g(i)%gprogt%spdm1  = sib%g(i)%gprogt%spdm2
       sib%g(i)%gprogt%lspr1  = sib%g(i)%gprogt%lspr2
       sib%g(i)%gprogt%cupr1  = sib%g(i)%gprogt%cupr2
       sib%g(i)%gprogt%dlwbot1 = sib%g(i)%gprogt%dlwbot2
       sib%g(i)%gprogt%sw_dwn1 = sib%g(i)%gprogt%sw_dwn2
   enddo

   ! switch files if needed
   driver_existf = .false.
   if ( driver_switchf ) then
       status = nf90_close( driverid )

       write(driver_filename,fmt='(a,i4.4,i2.2,a3)') trim(driver_path), &
             driver_year, driver_month,'.nc'
       inquire(file=trim(driver_filename), exist=driver_existf)
       if (driver_existf) then
           STATUS = nf90_open( trim(driver_filename), nf90_nowrite, driverid )
           IF (status .ne. nf90_noerr) THEN
               print*,'Error Opening Driver File: '
               print*,'  ',trim(driver_filename)
               STOP
           ENDIF
           !print*,'   Reading driver file: ',trim(driver_filename)
        else
            driverid = -1
        endif
    endif

    ! read data if driver file is open
    if (driverid > 0) then

        !...check nsib 
        CHECK( nf90_inq_dimid(driverid, trim(nsibname), dimid ) )
        CHECK( nf90_inquire_dimension(driverid, dimid, dim_name,dimlen ) )
        if ( dimlen /= nsib ) then
             print*,'File and specified nsib are different.  Stopping.'
             print*,'  File nsib: ',dimlen,'Sim nsib: ', nsib
             stop
        endif

        !...check time values in driver data file
        ENSURE_VAR( driverid,'year', ncyid )
        ENSURE_VAR( driverid,'month',ncmid )
        ENSURE_VAR( driverid,'day',  nctdid )
        mstart(1) = driver_recnum

        !CHECK( nf90_get_var( driverid, ncyid, xyear, mstart(1:1) ) )
        !if (xyear .ne. driver_year) then
        !   print*,'Year in file does not match simulation. Stopping.'
        !   print*,'   File: ',xyear,' Sim: ',year
        !   stop
        !endif

        CHECK( nf90_get_var( driverid, ncmid, xmonth, mstart(1:1) ) )
        if (xmonth .ne. driver_month) then
           print*,'Month in file does not match simulation. Stopping.'
           print*,'   File: ',xmonth,' Sim: ',month
           stop
        endif

        CHECK( nf90_get_var( driverid, nctdid, xday, mstart(1:1) ) )
        if (xday .ne. driver_day) then
           print*,'Day in file does not match simulation. Stopping.'
           print*,'   File: ',xday,' Sim: ', driver_day
           stop
        endif

        !...get veriable id's
        ENSURE_VAR( driverid, trim(psname), psid ) ! Surface Pressure
        ENSURE_VAR( driverid, trim(tmname), tmid ) ! Temperature at 2 m
        ENSURE_VAR( driverid, trim(shname), shid ) ! Humidity at 2 m
        ENSURE_VAR( driverid, trim(spdmname), spdmid ) ! Wind at 10 m
        ENSURE_VAR( driverid, trim(swdname), swdid ) ! Surface solar rad downwards
        ENSURE_VAR( driverid, trim(lwdname), lwdid ) ! Surface thermal rad down

        status = nf90_inq_varid(driverid, trim(lsprname), lsprid) ! Large Scale Precipitation
        if (status /= nf90_noerr) then
           status = nf90_inq_varid(driverid, trim(lsprnameopt), lsprid)
           if (status /= nf90_noerr) then 
              print*,'Missing Large Scale Precip Data.  Stopping.'
              stop
           endif
        endif

        status = nf90_inq_varid( driverid, trim(cuprname), cuprid ) ! Convective Precipitation
        if (status /= nf90_noerr) then
            status = nf90_inq_varid(driverid, trim(cuprnameopt), cuprid)
            if (status /= nf90_noerr) then
               print*,'Missing Convective Precip Data.  Stopping.'
               stop
            endif
         endif

        !...get data
        mstart=(/1,driver_recnum/); mcount=(/nsib,1/)
        STATUS = nf90_get_var(driverid, psid, ps, mstart, mcount) !Pressure
        IF (status .ne. nf90_noerr) THEN
            print*,'Error Getting Pressure'
            STOP
        ENDIF

        STATUS = nf90_get_var(driverid, tmid, tm, mstart, mcount) !Temperature
        IF (status .ne. nf90_noerr) THEN
            print*,'Error Getting Temperature'
            STOP
        ENDIF

        STATUS = nf90_get_var(driverid, shid, sh, mstart, mcount) !Humidity
        IF (status .ne. nf90_noerr) THEN
            print*,'Error Getting Humidity'
            STOP
        ENDIF

        STATUS = nf90_get_var(driverid, spdmid, spdm, mstart, mcount) !Wind
        IF (status .ne. nf90_noerr) THEN
           print*,'Error Getting Wind'
           STOP
        ENDIF

        STATUS = nf90_get_var(driverid, lsprid, lspr, mstart, mcount) !Large Scale Precip
        IF (status .ne. nf90_noerr) THEN
            print*,'Error Getting Large Scale Precip!'
            STOP
        ENDIF
        where(lspr .lt. 1.e-20) lspr=rzero

        STATUS = nf90_get_var(driverid, cuprid, cupr, mstart, mcount) !Convective Precip
        IF (status .ne. nf90_noerr) THEN
           print*,'Error Getting Convective Precip!'
           STOP
        ENDIF
        where(cupr .lt. 1.e-20) cupr=rzero

        STATUS = nf90_get_var(driverid, lwdid, lwd, mstart, mcount) !Longwave Rad
        IF (status .ne. nf90_noerr) THEN
           print*,'Error Getting Longwave Radiation'
           STOP
        ENDIF

        STATUS = nf90_get_var(driverid, swdid, swd, mstart, mcount) !Shortwave Rad
        IF (status .ne. nf90_noerr) THEN
           print*,'Error Getting Shortwave Radiation'
           STOP
        ENDIF

        !...pull out points in subdomain
        do i=1, subcount
           sib%g(i)%gprogt%ps2 = ps(subset(i))
           sib%g(i)%gprogt%tm2 = tm(subset(i))
           sib%g(i)%gprogt%sh2 = sh(subset(i))
           sib%g(i)%gprogt%spdm2 = spdm(subset(i))
           sib%g(i)%gprogt%lspr2 = lspr(subset(i))
           sib%g(i)%gprogt%cupr2 = cupr(subset(i))
           sib%g(i)%gprogt%sw_dwn2 = swd(subset(i))
           sib%g(i)%gprogt%dlwbot2 = lwd(subset(i))
        enddo

        !...print out the new data if requested
        if (print_driver) then
            print*,'New driver data read (yr/mon/day/hr/sec/recnum): '
            print('(a,i4.4,i5,i6,i5,i12,i6)'),'   ',driver_year,driver_month,driver_day, &
                         driver_hour,driver_secnext,driver_recnum
            print*,'------------------------------------------------------'
            print*,'Extrema of new input data'
            print*, minval(ps),maxval(ps),' Pressure (mb)'
            print*, minval(tm),maxval(tm),' Temperature (K)'
            print*, minval(sh),maxval(sh),' Water Vapor Mixing Ratio (kg/kg)'
            print*, minval(spdm),maxval(spdm),' Wind Speed (m/s)'
            print*, minval(lspr),maxval(lspr),' Large Scale Precip (mm/s)'
            print*, minval(cupr),maxval(cupr),' Convective Precip (mm/s)'
            print*, minval(swd),maxval(swd),' Short wave down (W/m2)'
            print*, minval(lwd),maxval(lwd),' Long wave down (W/m2)'
            print*,'-----------------------------------------------------'

            if (print_stop) stop
         endif

   else  !driver data file not found
         print*,'!!!Driver Data File Not Found.'
         print*,'   File: ',trim(driver_filename)
         print*,'   ID: ',driverid
         !stop
   endif

end subroutine driver_read_global

