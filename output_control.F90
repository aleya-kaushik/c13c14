!--------------------------------
subroutine output_control( )
!--------------------------------

!Controls all output (restart and diagnostic)
!   - Updates output time information
!   - Writes output to a file
!   - Updates the file names
!   - Switches output files if necessary

use kinds
use module_pparams, only: &
    secs_per_day, days_per_month
use module_io
use module_sibconst
use module_time

implicit none

! local variables
integer(i4) :: secref, daysref, monref, yearref
integer(i4) :: daystot, xmon

!-----------------------------------
!RESTART FILES
!-----------------------------------
IF (restart_savef) THEN
   !...update counters
   restart_countt = restart_countt + 1

   IF (restart_countt .GT. restart_ntpersave) THEN
        !...write the restart file
        call restart_write(.false.)

        !...reset counter and filename
        restart_countt = 1

        IF (restart_dtsib > 0) THEN
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
        ELSE !monthly

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
                   (((mod(yearref,4) .eq. 0) .and. &
                   (mod(yearref,100) .ne. 0)) .or. &
                   (mod(yearref,400) .eq. 0))) then
                    daysref = daysref + 1
               endif
               daystot = daystot + daysref
            enddo
            restart_ntpersave = daystot * steps_per_day

            monref = monref + 1
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

      ENDIF !restart_dtsib
   ENDIF !restart_countt
ENDIF !restart_savef


!--------------------------------
!HR FILES
!--------------------------------
IF ((hr_savegf) .OR. (hr_saveluf)) THEN
    !...update counters 
    hr_countt = hr_countt + 1
    IF (hr_countt .GT. hr_ntpersave) THEN
        hr_countt = 1
        hr_countsave = hr_countsave+1
        IF (hr_countsave .GT. hr_nsaveperout) THEN
            hr_writef = .true.
            if (hr_countwrite .eq. 1) hr_openf = .true.
            hr_countsave = 1
         ENDIf
    ENDIF 
 
    !...create output files
    IF (hr_openf) THEN
        !...gridcell variables
        if (hr_savegf) THEN
            call diagnostic_createg( &
                 year, curday_per_year, &
                 hr_filenameg, hr_nvarg, subcount, &
                 sublonsib, sublatsib, subsitename, &
                 hr_start, hr_step, &
                 hr_nsaveperfile,  &
                 hr_nameoutg, hr_listoutg,   &
                 hr_idg, hr_idvarg, hr_idtimeg)
         endif

         if (hr_saveluf) THEN
             call diagnostic_createlu( &
                   year, curday_per_year, &
                   hr_filenamelu, hr_nvarlu, subcount, &
                   sublonsib, sublatsib, subsitename, &
                   sublarea, subpref,     &
                   hr_start, hr_step,     &
                   hr_nsaveperfile,       &  
                   hr_nameoutlu, hr_listoutlu,    &
                   hr_idlu, hr_idvarlu, hr_idtimelu)
         ENDIF

         hr_openf = .false.
    ENDIF

    !...write output
    IF (hr_writef) THEN
        IF (hr_savegf) THEN
            hr_g = hr_g * hr_wtpersave
            call diagnostic_writeg( &
                 subcount, hr_ntpersave, hr_nsaveperout, &
                 hr_nvarg, hr_idg, hr_countwrite, &
                 hr_vrefg, hr_idvarg, hr_sif_satcount, &
                 hr_g)
            hr_g = dzero
            hr_sif_satcount = izero
        ENDIF

        IF (hr_saveluf) THEN
            hr_lu = hr_lu * hr_wtpersave
            call diagnostic_writelu( &
                 subcount, hr_nsaveperout, hr_idlu, &
                 hr_nvarlu, nlu, hr_idvarlu, &
                 hr_countwrite, hr_lu )
            hr_lu = dzero
        ENDIF

        hr_countwrite = hr_countwrite + hr_nsaveperout
        hr_writef = .false.
    ENDIF

    !...close output files
    IF (hr_countwrite .GT. hr_nsaveperfile) THEN
       IF (hr_savegf)  call diagnostic_close(hr_idg)
       IF (hr_saveluf) call diagnostic_close(hr_idlu)

       !.....prepare for next files
       if (hr_dtsib > izero) then
          !.....daily and hourly 
          monref = month
          yearref = year

          hr_start = doy-1
          hr_nsaveperfile = curday_per_mon*hr_nsaveperout
          if ((year .eq. endyear) .and. &
              ((doy+curday_per_mon) .gt. endtime)) then
              hr_nsaveperfile = (endtime - doy + 1)*hr_nsaveperout
          endif
        else
          !.....monthly
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
          hr_ntpersave = daystot*steps_per_day
          hr_wtpersave = 1./real(hr_ntpersave) 

       endif  !daily/hourly vs. monthly

       hr_countwrite = 1
       call diagnostic_filename('hr', &
            yearref, monref, &
            hr_filenameg, hr_filenamelu)

    ENDIF  !closing file and resetting

ENDIF !HR Files


!--------------------------------
!QP FILES
!--------------------------------
IF ((qp_savegf) .OR. (qp_saveluf)) THEN
    !...update counters 
    qp_countt = qp_countt + 1
    IF (qp_countt .GT. qp_ntpersave) THEN
        qp_countt = 1
        qp_countsave = qp_countsave+1
        IF (qp_countsave .GT. qp_nsaveperout) THEN
            qp_writef = .true.
            if (qp_countwrite .eq. 1) qp_openf = .true.
            qp_countsave = 1
         ENDIf
    ENDIF 
 
    !...create output files
    IF (qp_openf) THEN
        !...gridcell variables
        if (qp_savegf) THEN
            call diagnostic_createg( &
                 year, curday_per_year, &
                 qp_filenameg, qp_nvarg, subcount, &
                 sublonsib, sublatsib, subsitename, &
                 qp_start, qp_step, &
                 qp_nsaveperfile,  &
                 qp_nameoutg, qp_listoutg,   &
                 qp_idg, qp_idvarg, qp_idtimeg)
         endif

         if (qp_saveluf) THEN
             call diagnostic_createlu( &
                   year, curday_per_year, &
                   qp_filenamelu, qp_nvarlu, subcount, &
                   sublonsib, sublatsib, subsitename, &
                   sublarea, subpref,     &
                   qp_start, qp_step,     &
                   qp_nsaveperfile,       &  
                   qp_nameoutlu, qp_listoutlu,    &
                   qp_idlu, qp_idvarlu, qp_idtimelu)
         ENDIF

         qp_openf = .false.
    ENDIF

    !...write output
    IF (qp_writef) THEN
        IF (qp_savegf) THEN
            qp_g = qp_g * qp_wtpersave
            call diagnostic_writeg( &
                 subcount, qp_ntpersave, qp_nsaveperout, &
                 qp_nvarg, qp_idg, qp_countwrite, &
                 qp_vrefg, qp_idvarg, qp_sif_satcount, &
                 qp_g)
            qp_g = dzero
            qp_sif_satcount = izero
        ENDIF

        IF (qp_saveluf) THEN
            qp_lu = qp_lu * qp_wtpersave
            call diagnostic_writelu( &
                 subcount, qp_nsaveperout, qp_idlu, &
                 qp_nvarlu, nlu, qp_idvarlu, &
                 qp_countwrite, qp_lu )
            qp_lu = dzero
        ENDIF

        qp_countwrite = qp_countwrite + qp_nsaveperout
        qp_writef = .false.
    ENDIF

    !...close output files
    IF (qp_countwrite .GT. qp_nsaveperfile) THEN
       IF (qp_savegf)  call diagnostic_close(qp_idg)
       IF (qp_saveluf) call diagnostic_close(qp_idlu)

       !.....prepare for next files
       if (qp_dtsib > izero) then
          !.....daily and hourly 
          monref = month
          yearref = year

          qp_start = doy-1
          qp_nsaveperfile = curday_per_mon*qp_nsaveperout
          if ((year .eq. endyear) .and. &
              ((doy+curday_per_mon) .gt. endtime)) then
              qp_nsaveperfile = (endtime - doy + 1)*qp_nsaveperout
          endif
        else
          !.....monthly
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
          qp_ntpersave = daystot*steps_per_day
          qp_wtpersave = 1./real(qp_ntpersave) 

       endif  !daily/hourly vs. monthly

       qp_countwrite = 1
       call diagnostic_filename('qp', &
            yearref, monref, &
            qp_filenameg, qp_filenamelu)

    ENDIF  !closing file and resetting
ENDIF !QP Files


!--------------------------------
!PBP FILES
!--------------------------------
IF ((pbp_savegf) .OR. (pbp_saveluf)) THEN
    !...update counters 
    pbp_countt = pbp_countt + 1
    IF (pbp_countt .GT. pbp_ntpersave) THEN
        pbp_countt = 1
        pbp_countsave = pbp_countsave+1
        IF (pbp_countsave .GT. pbp_nsaveperout) THEN
            pbp_writef = .true.
            if (pbp_countwrite .eq. 1) pbp_openf = .true.
            pbp_countsave = 1
         ENDIf
    ENDIF 
 
    !...create output files
    IF (pbp_openf) THEN
        !...gridcell variables
        if (pbp_savegf) THEN
            call diagnostic_createg( &
                 year, curday_per_year, &
                 pbp_filenameg, pbp_nvarg, subcount, &
                 sublonsib, sublatsib, subsitename, &
                 pbp_start, pbp_step, &
                 pbp_nsaveperfile,  &
                 pbp_nameoutg, pbp_listoutg,   &
                 pbp_idg, pbp_idvarg, pbp_idtimeg)
         endif

         if (pbp_saveluf) THEN
             call diagnostic_createlu( &
                   year, curday_per_year, &
                   pbp_filenamelu, pbp_nvarlu, subcount, &
                   sublonsib, sublatsib, subsitename, &
                   sublarea, subpref,     &
                   pbp_start, pbp_step,     &
                   pbp_nsaveperfile,       &  
                   pbp_nameoutlu, pbp_listoutlu,    &
                   pbp_idlu, pbp_idvarlu, pbp_idtimelu)
         ENDIF

         pbp_openf = .false.
    ENDIF

    !...write output
    IF (pbp_writef) THEN
        IF (pbp_savegf) THEN
            pbp_g = pbp_g * pbp_wtpersave
            call diagnostic_writeg( &
                 subcount, pbp_ntpersave, pbp_nsaveperout, &
                 pbp_nvarg, pbp_idg, pbp_countwrite, &
                 pbp_vrefg, pbp_idvarg, pbp_sif_satcount, &
                 pbp_g)
            pbp_g = dzero
            pbp_sif_satcount = izero
        ENDIF

        IF (pbp_saveluf) THEN
            pbp_lu = pbp_lu * pbp_wtpersave
            call diagnostic_writelu( &
                 subcount, pbp_nsaveperout, pbp_idlu, &
                 pbp_nvarlu, nlu, pbp_idvarlu, &
                 pbp_countwrite, pbp_lu )
            pbp_lu = dzero
        ENDIF

        pbp_countwrite = pbp_countwrite + pbp_nsaveperout
        pbp_writef = .false.
    ENDIF

    !...close output files
    IF (pbp_countwrite .GT. pbp_nsaveperfile) THEN
       IF (pbp_savegf)  call diagnostic_close(pbp_idg)
       IF (pbp_saveluf) call diagnostic_close(pbp_idlu)

       !.....prepare for next files
       if (pbp_dtsib > izero) then
          !.....daily and hourly 
          monref = month
          yearref = year

          pbp_start = doy-1
          pbp_nsaveperfile = curday_per_mon*pbp_nsaveperout
          if ((year .eq. endyear) .and. &
              ((doy+curday_per_mon) .gt. endtime)) then
              pbp_nsaveperfile = (endtime - doy + 1)*pbp_nsaveperout
          endif
        else
          !.....monthly
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
          pbp_ntpersave = daystot*steps_per_day
          pbp_wtpersave = 1./real(pbp_ntpersave) 

       endif  !daily/hourly vs. monthly

       pbp_countwrite = 1
       call diagnostic_filename('pbp', &
            yearref, monref, &
            pbp_filenameg, pbp_filenamelu)

    ENDIF  !closing file and resetting
ENDIF !PBP Files


end subroutine output_control
