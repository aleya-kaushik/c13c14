subroutine driver_read_single(gprogt)

!****--------------------------------------------------------------------
!    This subroutines reads the forcing surface meteorological data 
!    for the current time step for single site (point) simulations.
!    If required, it closes the current month's data file and opens the 
!    next month's data file.
!****--------------------------------------------------------------------

    use kinds
    use module_io
    use module_sib, only: gprog_type
    use module_sibconst, only: &
       print_driver, print_stop

    implicit none

    !...input variables
    type(gprog_type) :: gprogt

    !...local variables
    real(r8) :: psin,tmin,shin,spdmin,lsprin,cuprin,dlwbotin,swdwnin

    !...misc variables
    integer(i4) :: status
    logical :: opened, existed
    real(r8) :: yr, doy, hr
    character(len=256) :: record

    !...storing previous time steps data
    gprogt%ps1       = gprogt%ps2
    gprogt%tm1       = gprogt%tm2
    gprogt%sh1       = gprogt%sh2
    gprogt%spdm1     = gprogt%spdm2
    gprogt%lspr1     = gprogt%lspr2
    gprogt%cupr1     = gprogt%cupr2
    gprogt%dlwbot1   = gprogt%dlwbot2
    gprogt%sw_dwn1   = gprogt%sw_dwn2

    !...switch files if needed
    if ( driver_switchf ) then  !switch and read
        inquire( unit=driverid,exist=existed,opened=opened)
        if (opened) close( driverid, iostat = status )

        write(unit=driver_filename,fmt='(a,i4.4,i2.2,a4)')  trim(driver_path), &
              driver_year, driver_month, '.dat'
        inquire (file=trim(driver_filename), exist=existed)

       if (existed) then
            if (print_driver) then
                 print*,'Opening driver file: '
                 print*,'  ',trim(driver_filename)
            endif
            open( unit=driverid, file=trim(driver_filename), status='old', &
                  form='formatted', iostat=status)
            if ( status > 0 ) then
                 print *, '!!!Error opening driver file!!!'
                 stop
            endif
        else  !driver_existf == .false.
            print*,''
            print*,'Stopping due to non-existant driver file: '
            print*,'  ',trim(driver_filename)
            stop
        endif 
    endif !switch files

    if (driver_recnum .gt. 0) then
        DO  ! Read until not a comment.
            read( driverid,'(a)', iostat=status ) record
            IF ( status /= 0 ) THEN
                print*,'Stopping due to error in single driver data.'
                stop
            ENDIF
            IF (record(1:1) .ne. '#') THEN
               exit
            ENDIF
         ENDDO

         read(unit=record,fmt=*) yr,doy,hr,tmin,shin, &
                spdmin,psin,dlwbotin,swdwnin,lsprin,cuprin
     else
         tmin=nodata
         shin=nodata
         spdmin=nodata
         psin=nodata
         dlwbotin=nodata
         swdwnin=nodata
         lsprin=nodata
         cuprin=nodata
     endif

     !...set the new variables
     !...check for nodata values
     if (tmin==nodata) tmin=gprogt%tm1
     if (shin==nodata) shin=gprogt%sh1
     if (spdmin==nodata) spdmin=gprogt%spdm1
     if (psin==nodata) psin=gprogt%ps1
     if (dlwbotin==nodata) dlwbotin=gprogt%dlwbot1
     if (swdwnin==nodata) swdwnin=gprogt%sw_dwn1
     if (lsprin==nodata) lsprin=gprogt%lspr1
     if (cuprin==nodata) cuprin=gprogt%cupr1

     !...copy data into data structure
     gprogt%tm2 = tmin
     gprogt%sh2 = shin
     gprogt%spdm2 = spdmin
     gprogt%ps2 = psin
     gprogt%dlwbot2 = dlwbotin
     gprogt%sw_dwn2 = swdwnin
     gprogt%lspr2 = lsprin 
     gprogt%cupr2 = cuprin

     !...print out the new data if requested
     if (print_driver) then
        print*,'New driver data read (yr/mon/day/hr): '
        print*,'   ',driver_year,driver_month,driver_day,driver_hour
        print*,'------------------------------------------------------'
        print*,'Extrema of new input data'
        print*, gprogt%ps2,' Pressure (mb)'
        print*, gprogt%tm2,' Temperature (K)'
        print*, gprogt%sh2,' Water Vapor Mixing Ratio (kg/kg)'
        print*, gprogt%spdm2,' Wind Speed (m/s)'
        print*, gprogt%lspr2,' Large Scale Precip (mm/s)'
        print*, gprogt%cupr2,' Convective Precip (mm/s)'
        print*, gprogt%sw_dwn2,' Short wave down (W/m2)'
        print*, gprogt%dlwbot2,' Long wave down (W/m2)'
        print*,'-----------------------------------------------------'

        if (print_stop) stop
     endif

end subroutine driver_read_single
