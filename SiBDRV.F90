program sibdrive
  
    use kinds
    use module_io, only: rank, nchunks
    use module_pparams, only: pi180, month_names
    use module_sib, only: sib
    use module_sibconst
    use module_time
    use module_param, only: physcon

    implicit none

    ! local variables
    integer(i4) :: i,n,pnum
    integer(long) :: t
    integer(i4), external :: iargc
    character(len=200) :: namel_name
    character(len=16) :: buf

    ! variables for timing
    real etime          ! Declare the type of etime()
    real elapsed(2)     ! For receiving user and system time
    real total          ! For receiving total time

    ! parameters
    logical, parameter :: show_year_info = .true.
    logical, parameter :: show_mon_info = .false.

    !-----------------------------------------------------
    !...read in parallelization values from command line
    n = command_argument_count()
    if (n == 0) then
       namel_name = 'namel_sibdrv'
       rank = 1
       nchunks = 1
    elseif (n == 1) then
       call parse_name( namel_name )
       rank = 1
       nchunks = 1
    elseif (n == 2) then
       namel_name = 'namel_sibdrv'
       call parse_args(rank, nchunks)
    elseif (n == 3) then
       call get_command_argument(1, namel_name)
       call get_command_argument(2, buf)
       read (buf, *) rank
       call get_command_argument(3, buf)
       read(buf, *) nchunks
    else
       print *,'Wrong number of arguments to SiB4.'
       print *,'Stopping.'
       stop
    endif

    !...read in namel_sibdrv
    call read_namel(namel_name)

    !...read in PFT and pool information
    call read_pftinfo()
    call read_poolinfo()

    !...read in vegetation structure 
    call read_sibvs()

    !...initialize the grid
    call grid_init()

    !...read in output specifications
    call read_outopts()

    !...read in restart specifications
    call read_routopts()

    !...initialize time
    call time_init()

    !...read in parameters
    call read_ciso()
    call read_aero()
    call read_pgdd()
    call read_pstg()
    call read_phys()
    call read_pool()

    !...initialize sib structure
    call sibtype_init()

    !...initialize driver data
    call driver_init()

    if (tm5mr_switch) then
         !...initialize tm5 data !ARACHO
         call tm5_init()
    endif

    !...initialize fire emissions
    call fire_init()
    
    !...read in restart variables
    call restart_read()

    !...setup sib diagnostic variables
    call sibtype_setup()

    !...initialize output
    call output_init()

    !...reset pools for continuing spinup
    if (spinup_continue) then
      call set_continue_spinup()
    endif

    print*,''
    if (spinup .and. (.not. spinup_continue)) print*,'Starting Spinup Simulation'
    if (spinup .and. spinup_continue) print*,'Continuing Spinup Simulation'

    print*,'Starting timestep loop'
    do while ((spinup_lnum .le. spinup_maxiter) .and. (.not. spinup_done))

        !...timestep loop
        do t = start_second, end_second - dtsib, dtsib
           if (new_year) then
             print('(a,i4)'),' Processing: ', year
           endif

           do i=1,subcount
             !itb...reset accumulated potential ET
             sib%g(i)%l%fluxt%et0a = 0.0
             sib%g(i)%l%fluxt%et1a = 0.0
             !itb...reset growing season diagnostics
             do n=1,7
               sib%g(i)%l%co2t%gs_stress(n) = 0.0
             enddo
           enddo

          if (new_month .and. show_mon_info) then
              !...Print output information once a month
              call print_moninfo()
          endif

          if (new_day) then
               !...Calculate solar declination 
               call solar_dec(doy, cureqnx, lonearth, &
                    sin_dec, cos_dec, tan_dec)

               !...Calculate day length
               call day_length(subcount, pi180, tan_dec, sublatsib, &
                    sib%g%gdiagt%daylen, sib%g%gdiagt%daylendt)
           endif !new_day

           !...Prepare driver meteorology
           call driver_update()

           if (tm5mr_switch) then
               !...Prepare tm5 data !ARACHO
               call tm5_update(sib%g%gprogt)
           endif

           !...Prepare fire emissions
           call fire_update()
           
           !...loop over points    
           do i = 1, subcount
              if (.not. sib%g(i)%gdiagt%gridcell_spunup) then
     
                  !...Interpolate driver data each time step
                  call driver_interp(i, sublonsib(i), sublatsib(i), &
                       sib%g(i)%gdiagt, sib%g(i)%gprogt)

                  if (tm5mr_switch) then
                       !...Interpolate tm5 data each time step !ARACHO
                       call tm5_interp(i, sublonsib(i), sublatsib(i), &
                            sib%g(i)%gdiagt, sib%g(i)%gprogt)
                  endif

                  !...Interpolate fire emissions each time step
                  call fire_interp(i, sublonsib(i), sublatsib(i), &
                       sib%g(i))

                  !call fire_interp_c13(i, sublonsib(i), sublatsib(i), &
                  !     sib%g(i), physcon(pnum))

                  !...Set CO2 values
                  !if (varco2_switch) then
                  if (new_day) then
                    call set_co2( year, nmonth, sib%g(i)%gprogt )
                  endif
                  !endif

                  call set_cos(sib%g(i)%gprogt)
                  
                  !...Call SiB control
                  call sib_control(i, subset(i),   &
                       sublonsib(i), sublatsib(i), &
                       new_day_lst(i), daylenmax(i), &
                       doy, sib%g(i))

                  !...Save gridcell diagnostics
                  call diagnostic_save_gall( &
                       i, sib%g(i), sib%g(i)%gprogt )

                endif !.not. gridcell_spunup
           enddo  !i=1,subcount

           !...Set CO2 values
           !if (varco2_switch) then
           !  if (new_day) then
           !    call set_co2( year, nmonth )
           !  endif
           !endif

           !...Set COS values
           !if (new_day) then
           !call set_cos(sib%g%gprogt)
           !endif

           !...Increment the time
           call time_manager()

           !...Write any requested output
           call output_control()

           !...Calculate and write the equilibrium pools if necessary
           call equipools_control()

        enddo ! timestep loop
    enddo ! spinup loop

    !...Ensure all files have been closed
    call output_closer()

    !...Print final message
    total = etime(elapsed)
    print*, ''
    print*, 'Times: total=', total, ' user=', elapsed(1),' system=', elapsed(2)
    print*, 'End Simulation'
    print*,''

end program sibdrive

!====================================
subroutine parse_args(rank, nchunks)
    integer, intent(out) :: rank, nchunks

    character(len=16) :: buf

    call getarg(1, buf)
    if (buf == '' .or. buf == '>') then
        rank = 1
        nchunks = 1
    else
        read(buf, *) rank
        call getarg(2, buf)
        if (buf == '' .or. buf == '>') then
            stop 'Command line arguments incorrect:  SiBD4 rank nchunks'
        else
            read(buf, *) nchunks
            if (rank > nchunks) stop 'rank greater than nchunks'
            if (rank < 1 .or. nchunks < 1) stop 'rank or nchunks < 1'
        endif
    endif
end subroutine parse_args

!===================================
subroutine parse_name( namel_name )
    character(len=120), intent(out) :: namel_name

    character(len=120) :: buf

    call getarg(1, buf)
    if (buf == '' .or. buf == '>') then
        namel_name = 'namel_sibdrv'
    else
        namel_name = buf
    endif
end subroutine parse_name

!===========================
subroutine print_moninfo()

use module_io
use module_pparams, only: &
    month_names
use module_time

implicit none

integer(i4) :: indx1, indx2, xmoni, xmon
character(len=10) :: lstring

!-----------------------------------------------------------
print('(a,a,i6)'),'     ',month_names(month), year

!.....driver data information
print('(a)'),'     Driver Data: '
indx1 = INDEX(driver_filename,'/',back=.true.)
if (indx1 .gt. 1) then
   indx1 = INDEX(driver_filename(1:indx1-1),'/',back=.true.)
endif
if (indx1 .gt. 1) then
   indx1 = INDEX(driver_filename(1:indx1-1),'/',back=.true.)
endif
if (indx1 .gt. 1) then
   lstring = '       ...'
else
   lstring = '          '
   indx1 = 1
endif
indx2 = len(driver_filename)
print('(a,a)'),lstring,trim(driver_filename(indx1:indx2))

!.....restart information
IF (restart_savef) THEN
    IF (restart_dtsib .GT. 0) THEN
        print('(a)'),' Writing Restart Output'
    ELSE
        xmoni = index(trim(restart_filename),'sib_r')
        read(restart_filename(xmoni+9:xmoni+10),'(i2)') xmon
        if (xmon .EQ. month) then
            print('(a)'),'     Restart Output: '
            indx1 = INDEX(restart_filename,'/',back=.true.)
            if (indx1 .gt. 1) then
               indx1 = INDEX(restart_filename(1:indx1-1),'/',back=.true.)
            endif
            if (indx1 .gt. 1) then
               indx1 = INDEX(restart_filename(1:indx1-1),'/',back=.true.)
            endif
            if (indx1 .gt. 1) then
               lstring = '       ...'
            else
               lstring = '          '
              indx1 = 1
            endif
            indx2 = len(restart_filename)
            print('(a,a)'),lstring,trim(restart_filename(indx1:indx2))
        endif
     ENDIF
ENDIF

!.....diagnostic output information
IF ((hr_savegf) .OR. (hr_saveluf)) THEN
    xmoni = index(trim(hr_filenameg),'hsib_')
    read(hr_filenameg(xmoni+9:xmoni+10),'(i2)') xmon
    if (xmon .EQ. month) then
        print('(a)'),'     HR Output: '
        IF (hr_savegf)  print('(a,a)'),'      ',trim(hr_filenameg)
        IF (hr_saveluf) print('(a,a)'),'      ',trim(hr_filenamelu)
    endif
ENDIF
IF ((pbp_savegf) .OR. (pbp_saveluf)) THEN
    xmoni = index(trim(pbp_filenameg),'psib_')
    read(pbp_filenameg(xmoni+9:xmoni+10),'(i2)') xmon
    if (xmon .EQ. month) then
       print('(a)'),'     PBP Output: '
       IF (pbp_savegf)  print('(a,a)'),'      ',trim(pbp_filenameg)
       IF (pbp_saveluf) print('(a,a)'),'      ',trim(pbp_filenamelu)
    endif
ENDIF
IF ((qp_savegf) .OR. (qp_saveluf)) THEN
    xmoni = index(trim(qp_filenameg),'qsib_')
    read(qp_filenameg(xmoni+9:xmoni+10),'(i2)') xmon
    if (xmon .EQ. month) then
        print('(a)'),'     QP Output: '
        IF (qp_savegf)  print('(a,a)'),'      ',trim(qp_filenameg)
        IF (qp_saveluf) print('(a,a)'),'      ',trim(qp_filenamelu)
     endif
ENDIF
print('(a)'),''

end subroutine print_moninfo
