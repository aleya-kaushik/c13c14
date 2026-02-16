!---------------------------------------------------------------------
subroutine output_init( )
!---------------------------------------------------------------------

!Initializes output information

use module_pparams, only:  &
    secs_per_day, days_per_month
use module_sibconst, only: &
     subcount, nlu
use module_io
use module_time

use kinds
implicit none

!...local variables
integer(i4) :: xmon
integer(i4) :: daysref, daystot, monref, yearref

    !-----------------------------------------
    !...Initialize output variables
    IF (restart_savef .OR. hr_savegf .OR. hr_saveluf .OR. &
         pbp_savegf .OR. pbp_saveluf .OR. &
         qp_savegf .OR. qp_saveluf) THEN
        print *, 'Setting Output'
    ENDIF

    !-----------------------------------------
    !...Restart Files
    !-----------------------------------------
    if (restart_savef) then
       if (restart_dtsib > 0) then
          if (mod(restart_dtsib, 86400) .ne. 0.) then
             !...Hourly (not allowed since SiB4 starts at 0 GMT)
             print('(a)'),'Restart Output Must Be Days or Months.  Stopping.'
             STOP
           else
              !..Daily
              daystot = restart_dtsib/86400.
              print('(a,i6)'), '    SiB Restart Files Written (days) = ', daystot
              restart_ntpersave = daystot * steps_per_day
              restart_countt = 1

              !...update doy for output file name
              daysref = doy + 1
              yearref = year
              if (doy .gt. curday_per_year) then
                 daysref = 1
                 yearref = yearref + 1
              endif

              if (printrank) then
                 write( restart_filename, '(a,i4.4,i3.3,a,i5.5,a)') trim(out_path)//'sib_r', &
                      yearref,daysref,'p',rank,'.nc'
              else
                 write( restart_filename, '(a,i4.4,i3.3,a)') trim(out_path)//'sib_r', &
                      yearref,daysref,'.nc'
              endif
           endif

       else !monthly
           print('(a,i4)'),'    SiB restart files written (month) = ', -restart_dtsib

           daystot = 0
           monref = month-1
           yearref = year
           do xmon=1,abs(restart_dtsib)
              monref = monref + 1
              if (monref .GT. 12) then
                 monref = 1
                 yearref = yearref + 1
              endif

              daysref = days_per_month(monref)
              if ((monref .eq. 2) .and. &
                  (((mod(endyear,4) .eq. 0) .and. &
                   (mod(endyear,100) .ne. 0)) .or. &
                   (mod(endyear,400) .eq. 0))) then
                  daysref = daysref + 1
              endif
              daystot = daystot + daysref 
           enddo

           restart_ntpersave = daystot * steps_per_day
           restart_countt = 1

           !...update mon/year ref for output file name
           monref = monref+1
           if (monref .gt. 12) then
               monref=1
               yearref=yearref+1
           endif

           if (printrank) then
              write( restart_filename, '(a,i4.4,i2.2,a,i5.5,a)') trim(out_path)//'sib_r', &
                     yearref,monref,'p',rank,'.nc'
           else
              write( restart_filename, '(a,i4.4,i2.2,a)') trim(out_path)//'sib_r', &
                     yearref,monref,'.nc'
           endif
       endif !restart_dtsib
    endif !restart_savef

    !-----------------------------------------
    !...HR Diagnostic Output (Typically Hourly)
    !-----------------------------------------
    if ((hr_savegf) .OR. (hr_saveluf)) then
       if (hr_dtsib > secs_per_day) then
           !...Multi-Day HR Output
           print('(a)'),'Not Expecting Multi-Day HR Output.'
           print('(a)'),'Stopping.'
           stop
       elseif (hr_dtsib > izero) then
           !...Daily and Hourly HR Output
           print('(a,i6)'), '    SiB HR files written (s) = ', hr_dtsib
         
           hr_start = doy-1
           hr_step = hr_dtsib / dble(secs_per_day)
           hr_ntpersave = hr_dtsib / dtsib
           hr_nsaveperout = dtoutwrite / hr_dtsib
           hr_nsaveperfile = curday_per_mon*hr_nsaveperout
           if ((year .eq. endyear) .and. &
               ((doy+curday_per_mon) .gt. endtime)) then
               hr_nsaveperfile = (endtime - doy + 1)*hr_nsaveperout
           endif
           monref = month
           yearref = year

       elseif (hr_dtsib > -13) then
           !...Monthly HR Output
           print('(a,i4)'), '    SiB HR files written (mon) = ', abs(hr_dtsib)

            daystot = izero
            monref = month-1
            yearref = year
            do xmon=1,abs(hr_dtsib)
               monref = monref + 1
               if (monref .GT. 12) then
                  monref = 1
                  yearref = yearref + 1
               endif

               daystot = daystot + curday_per_mon
           enddo

           hr_start = (doy-1) + daystot/dble(2.)
           if (hr_start .gt. curday_per_year) &
               hr_start = hr_start - curday_per_year

           hr_step = 0
           hr_ntpersave = daystot*steps_per_day
           hr_nsaveperout = 1
           hr_nsaveperfile = 1

       else
           !...Yearly HR Output
           print('(a)'),'Not Expecting Yearly HR Output.'
           print('(a)'),'Stopping.'
           stop
       endif

       !...General HR counters/holders
       hr_countt = 1
       hr_countsave = 1
       hr_wtpersave = 1. / real(hr_ntpersave)

       call diagnostic_filename('hr', yearref, monref, &
            hr_filenameg, hr_filenamelu)

        IF (hr_savegf) THEN
           allocate( hr_g(subcount,hr_nvarg,hr_nsaveperout))
           allocate( hr_idvarg(hr_nvarg))
           hr_g = dzero
           hr_idvarg = izero
        ENDIF
        IF (hr_saveluf) THEN
           allocate( hr_lu(subcount,nlu,hr_nvarlu,hr_nsaveperout))
           allocate( hr_idvarlu(hr_nvarlu))
           hr_lu = dzero
           hr_idvarlu = izero
        ENDIF
        
        !...General HR file info
        hr_openf = .false.
        hr_idg = izero
        hr_idlu = izero

        hr_writef = .false.
        hr_countwrite = 1

        !...SIF HR file info
        allocate(hr_sif_satcount(subcount,hr_nsaveperout,2))
        hr_sif_satcount = izero
    
    endif  !hrintg or hrintlu

    !-----------------------------------------
    !...QP Diagnostic Output (Typically Monthly)
    !-----------------------------------------
    if ((qp_savegf) .OR. (qp_saveluf)) then
       if (qp_dtsib > secs_per_day) then
           !...Multi-Day QP Output
           print('(a)'),'Not Expecting Multi-Day QP Output.'
           print('(a)'),'Stopping.'
           stop
       elseif (qp_dtsib > izero) then
           !...Daily and Hourly QP Output
           print('(a,i6)'), '    SiB QP files written (s) = ', qp_dtsib
         
           qp_start = doy-1
           qp_step = qp_dtsib / dble(secs_per_day)
           qp_ntpersave = qp_dtsib / dtsib
           qp_nsaveperout = dtoutwrite / qp_dtsib
           qp_nsaveperfile = curday_per_mon*qp_nsaveperout
           if ((year .eq. endyear) .and. &
               ((doy+curday_per_mon) .gt. endtime)) then
               qp_nsaveperfile = (endtime - doy + 1)*qp_nsaveperout
           endif
           monref = month
           yearref = year

       elseif (qp_dtsib > -13) then
           !...Monthly QP Output
           print('(a,i4)'), '    SiB QP files written (mon) = ', abs(qp_dtsib)

            daystot = izero
            monref = month-1
            yearref = year
            do xmon=1,abs(qp_dtsib)
               monref = monref + 1
               if (monref .GT. 12) then
                  monref = 1
                  yearref = yearref + 1
               endif

               daystot = daystot + curday_per_mon
           enddo

           qp_start = (doy-1) + daystot/dble(2.)
           if (qp_start .gt. curday_per_year) &
               qp_start = qp_start - curday_per_year

           qp_step = 0
           qp_ntpersave = daystot*steps_per_day
           qp_nsaveperout = 1
           qp_nsaveperfile = 1

       else
           !...Yearly QP Output
           print('(a)'),'Not Expecting Yearly QP Output.'
           print('(a)'),'Stopping.'
           stop
       endif

       !...General QP counters/holders
       qp_countt = 1
       qp_countsave = 1
       qp_wtpersave = 1. / real(qp_ntpersave)

       call diagnostic_filename('qp', yearref, monref, &
            qp_filenameg, qp_filenamelu)

        IF (qp_savegf) THEN
           allocate( qp_g(subcount,qp_nvarg,qp_nsaveperout))
           allocate( qp_idvarg(qp_nvarg))
           qp_g = dzero
           qp_idvarg = izero
        ENDIF
        IF (qp_saveluf) THEN
           allocate( qp_lu(subcount,nlu,qp_nvarlu,qp_nsaveperout))
           allocate( qp_idvarlu(qp_nvarlu))
           qp_lu = dzero
           qp_idvarlu = izero
        ENDIF

        !...General QP file info
        qp_openf = .false.
        qp_idg = izero
        qp_idlu = izero

        qp_writef = .false.
        qp_countwrite = 1

        !...SIF QP file info
        allocate(qp_sif_satcount(subcount,qp_nsaveperout,2))
        qp_sif_satcount = izero

    endif  !qp_savegf or qp_saveluf

    !-----------------------------------------
    !...PBP Diagnostic Output (Typically Daily)
    !-----------------------------------------
    if ((pbp_savegf) .or. (pbp_saveluf)) then

       if (pbp_dtsib > secs_per_day) then
           !...Multi-Day PBP Output
           print('(a)'),'Not Expecting Multi-Day PBP Output.'
           print('(a)'),'Stopping.'
           stop
       elseif (pbp_dtsib > izero) then
           !...Daily and Hourly PBP Output
           print('(a,i6)'), '    SiB PBP files written (s) = ', pbp_dtsib
         
           pbp_start = doy-1
           pbp_step = pbp_dtsib / dble(secs_per_day)
           pbp_ntpersave = pbp_dtsib / dtsib
           pbp_nsaveperout = dtoutwrite / pbp_dtsib

           daystot = curday_per_mon
           if ((year .eq. endyear) .and. &
               ((doy+curday_per_mon) .gt. endtime)) then
               daystot = (endtime - doy + 1)
           endif
           pbp_nsaveperfile = daystot*pbp_nsaveperout
           monref = month
           yearref = year

       elseif (pbp_dtsib > -13) then
           !...Monthly PBP Output
           print('(a,i4)'), '    SiB PBP files written (mon) = ', abs(pbp_dtsib)

            daystot = izero
            monref = month-1
            yearref = year
            do xmon=1,abs(pbp_dtsib)
               monref = monref + 1
               if (monref .GT. 12) then
                  monref = 1
                  yearref = yearref + 1
               endif

               daystot = daystot + curday_per_mon
           enddo

           pbp_start = (doy-1) + daystot/dble(2.)
           if (pbp_start .gt. curday_per_year) &
               pbp_start = pbp_start - curday_per_year

           pbp_step = 0
           pbp_ntpersave = daystot*steps_per_day
           pbp_nsaveperout = 1
           pbp_nsaveperfile = 1

       else
           !...Yearly PBP Output
           print('(a)'),'Not Expecting Yearly PBP Output.'
           print('(a)'),'Stopping.'
           stop
       endif

       !...General PBP counters/holders
       pbp_countt = 1
       pbp_countsave = 1
       pbp_wtpersave = 1. / real(pbp_ntpersave)

       call diagnostic_filename('pbp', yearref, monref, &
            pbp_filenameg, pbp_filenamelu)

        IF (pbp_savegf) THEN
           allocate( pbp_g(npbp,pbp_nvarg,pbp_nsaveperout))
           allocate(pbp_idvarg(pbp_nvarg))
           pbp_g = dzero
           pbp_idvarg = izero
        ENDIF
        IF (pbp_saveluf) THEN
           allocate( pbp_lu(npbp,nlu,pbp_nvarlu,pbp_nsaveperout))
           allocate(pbp_idvarlu(pbp_nvarlu))
            pbp_lu = dzero
            pbp_idvarlu = izero
         ENDIF
         
        !...General PBP file info
        pbp_openf = .false.
        pbp_idg = izero
        pbp_idlu = izero

        pbp_writef = .false.
        pbp_countwrite = 1

        !...SIF PBP file info
        allocate(pbp_sif_satcount(subcount,pbp_nsaveperout,2))
        pbp_sif_satcount = izero
        
    endif  !pbpintg or pbpintlu

end subroutine output_init

!---------------------------------------------------------------------
subroutine diagnostic_filename( &
           cout_type, yr, mon, filenameg, filenamelu) 
!---------------------------------------------------------------------

use kinds
use module_io

implicit none

!...input variables
character(len=*), intent(in) :: cout_type
integer(i4), intent(in) :: yr, mon
character(len=iofnlen), intent(inout) :: &
    filenameg, filenamelu

!...PBP Files
IF (cout_type .eq. 'pbp') THEN

   if (printrank) then
       write( filenameg, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//pbp_prefix, &
              yr,mon,'p',rank,trim(pbp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//pbp_prefix, &
              yr,mon,'p',rank,trim(pbp_suffixlu)
    else
       write( filenameg, '(a,i4.4,i2.2,a)') &
              trim(out_path)//pbp_prefix, &
              yr,mon,trim(pbp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a)') &
               trim(out_path)//pbp_prefix, &
               yr,mon,trim(pbp_suffixlu)
     endif

ELSEIF (cout_type .eq. 'qp') THEN

   if (printrank) then
       write( filenameg, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//qp_prefix, &
              yr,mon,'p',rank,trim(qp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//qp_prefix, &
              yr,mon,'p',rank,trim(qp_suffixlu)
    else
       write( filenameg, '(a,i4.4,i2.2,a)') &
              trim(out_path)//qp_prefix, &
              yr,mon,trim(qp_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a)') &
               trim(out_path)//qp_prefix, &
               yr,mon,trim(qp_suffixlu)
     endif

ELSEIF (cout_type .eq. 'hr') THEN

   if (printrank) then
       write( filenameg, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//hr_prefix, &
              yr,mon,'p',rank,trim(hr_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a,i5.5,a)') &
              trim(out_path)//hr_prefix, &
              yr,mon,'p',rank,trim(hr_suffixlu)
    else
       write( filenameg, '(a,i4.4,i2.2,a)') &
              trim(out_path)//hr_prefix, &
              yr,mon,trim(hr_suffixg)
       write( filenamelu, '(a,i4.4,i2.2,a)') &
               trim(out_path)//hr_prefix, &
               yr,mon,trim(hr_suffixlu)
     endif

ELSE
    print*,'Unknown Output File Type:',cout_type
    print*,'Stopping.'
    stop
ENDIF


end subroutine diagnostic_filename
