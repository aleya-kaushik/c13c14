!===============================================================================
subroutine read_outopts()
!===============================================================================

use module_sibconst, only: &
    nlu, nsoil, nsnow, ntot
use module_io

implicit none

integer :: numtem
logical :: doqptem, dopbptem, dohrtem
character (len=1) :: luorgtem
character (len=4) :: levtem
character (len=21) :: nametem
character (len=21) :: luorgnametem
character (len=100) :: listtemf, listtem
character (len=100) :: luorglisttem
integer :: i, count

    !---------------------------------------------------------------------------
    ! read sib_outopts and count number of variables to be output
    !---------------------------------------------------------------------------
 
    hr_nvarg = 0
    hr_nvarlu = 0
    pbp_nvarg = 0
    pbp_nvarlu = 0
    qp_nvarg = 0
    qp_nvarlu = 0

    print*,'Reading Output File: '
    print*,'  ',trim(out_info)
    open(unit=2,file=trim(out_info),form='formatted')

    do 
        read(2,*, end=1200) doqptem,dopbptem,dohrtem,luorgtem, &
                           levtem,nametem,numtem,listtem

        if (luorgtem == 'l') then
            if (doqptem) then
                if (levtem == 'sing') then
                    qp_nvarlu = qp_nvarlu + 1
                elseif (levtem == 'soil') then
                    qp_nvarlu = qp_nvarlu + nsoil
                elseif (levtem == 'stot') then
                    qp_nvarlu = qp_nvarlu + ntot
                else
                    print*,'Bad level descriptor in sib_outopts.'
                    print*,'Stopping.'
                    print*,''
                    stop
                endif
            endif

            if (dopbptem) then
                if (levtem == 'sing') then
                    pbp_nvarlu = pbp_nvarlu + 1
                elseif (levtem == 'soil') then
                    pbp_nvarlu = pbp_nvarlu + nsoil
                elseif (levtem == 'stot') then
                    pbp_nvarlu = pbp_nvarlu + ntot
                else
                    print*,'Bad level descriptor in sib_outopts.'
                    print*,'Stopping.'
                    print*,''
                    stop
                endif
             endif

            if (dohrtem) then
                if (levtem == 'sing') then
                    hr_nvarlu = hr_nvarlu + 1
                elseif (levtem == 'soil') then
                    hr_nvarlu = hr_nvarlu + nsoil
                elseif (levtem == 'stot') then
                    hr_nvarlu = hr_nvarlu + ntot
                else
                    print*,'Bad level descriptor in sib_outopts.'
                    print*,'Stopping.'
                    print*,''
                    stop
                endif
            endif

        elseif (luorgtem == 'g') then
            if (doqptem) then
                if (levtem == 'sing') then
                    qp_nvarg = qp_nvarg + 1
                else
                    print*,'Currently Expecting Only Single Levels'
                    print*,'Stopping'
                    print*,''
                    stop
                endif
            endif
            if (dopbptem) then
                if (levtem == 'sing') then
                    pbp_nvarg = pbp_nvarg + 1
                else
                    print*,'Currently Expecting Only Single Levels'
                    print*,'Stopping.'
                    print*,''
                    stop
                endif
            endif
            if (dohrtem) then
                if (levtem == 'sing') then
                    hr_nvarg = hr_nvarg + 1
                else
                    print*,'Currently Expecting Only Single Levels'
                    print*,'Stopping.'
                    print*,''
                    stop
                endif
            endif
        else
            print*,'Invalid Variable Structure Type in sib_outopts.  Stopping'
            stop
        endif
    enddo

    1200  continue
    rewind 2

    IF ((hr_nvarlu .gt. 0) .and. (hr_saveluf)) THEN
       allocate (hr_nameoutlu(hr_nvarlu), hr_listoutlu(hr_nvarlu), &
                 hr_vreflu(hr_nvarlu))
    ELSE
       hr_saveluf = .false.
    ENDIF
    IF ((hr_nvarg .gt. 0) .and. (hr_savegf)) THEN
        allocate (hr_nameoutg(hr_nvarg), hr_listoutg(hr_nvarg), hr_vrefg(hr_nvarg))
    ELSE
        hr_savegf = .false.
    ENDIF

    IF ((qp_nvarlu .gt. 0) .and. (qp_saveluf)) THEN
        allocate (qp_nameoutlu(qp_nvarlu), qp_listoutlu(qp_nvarlu), qp_vreflu(qp_nvarlu))
    ELSE
       qp_saveluf = .false.
    ENDIF
    IF ((qp_nvarg .gt. 0) .and. (qp_savegf)) THEN
       allocate (qp_nameoutg(qp_nvarg), qp_listoutg(qp_nvarg), qp_vrefg(qp_nvarg))
    ELSE
       qp_savegf = .false.
    ENDIF

    IF ((pbp_nvarlu .gt. 0) .and. (pbp_saveluf)) THEN
       allocate (pbp_nameoutlu(pbp_nvarlu), pbp_listoutlu(pbp_nvarlu), &
                 pbp_vreflu(pbp_nvarlu))
    ELSE
       pbp_saveluf = .false.
    ENDIF
    IF ((pbp_nvarg .gt. 0) .and. (pbp_savegf)) THEN
       allocate (pbp_nameoutg(pbp_nvarg), pbp_listoutg(pbp_nvarg), &
                 pbp_vrefg(pbp_nvarg))
    ELSE
       pbp_savegf = .false.
    ENDIF

    hr_nvarlu = 0
    hr_nvarg = 0
    qp_nvarlu = 0
    qp_nvarg = 0
    pbp_nvarlu = 0
    pbp_nvarg = 0

    do 
        read(2,*, end=942) doqptem,dopbptem,dohrtem,luorgtem, &
                           levtem,nametem,numtem,listtemf

        if (luorgtem == 'l') then
           if (doqptem .and. qp_saveluf) then
               if (levtem == 'sing') then
                   qp_nvarlu = qp_nvarlu + 1
                   qp_nameoutlu(qp_nvarlu) = nametem
                   write(listtem,'(a,a2,a)') trim(nametem), &
                           ': ',trim(listtemf)
                   qp_listoutlu(qp_nvarlu) = listtem
                   qp_vreflu(qp_nvarlu) = numtem
                elseif (levtem == 'stot') then
                    count=nsnow
                    do i=-nsnow+1, nsoil
                       qp_nvarlu = qp_nvarlu + 1
                       if (i < 1) then
                           write(luorgnametem,'(a,a1,i2.2)') trim(nametem),'_',count
                           write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                           write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Snow Level ', count
                           count = count - 1
                       else
                           write(luorgnametem,'(a,i2.2)') trim(nametem),i
                           write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                           write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Soil Level ',i
                       endif
                       qp_nameoutlu(qp_nvarlu) = luorgnametem
                       qp_listoutlu(qp_nvarlu) = luorglisttem
                       qp_vreflu(qp_nvarlu) = numtem
                     enddo
                 elseif (levtem == 'soil') then
                    do i=1, nsoil
                       qp_nvarlu = qp_nvarlu + 1
                       write(luorgnametem,'(a,i2.2)') trim(nametem), i
                       write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                       write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Soil Level ', i
                       qp_nameoutlu(qp_nvarlu) = luorgnametem
                       qp_listoutlu(qp_nvarlu) = luorglisttem
                       qp_vreflu(qp_nvarlu) = numtem
                    enddo
                endif
            endif  !doqptem

            if (dopbptem .and. pbp_saveluf) then
               if (levtem == 'sing') then
                   pbp_nvarlu = pbp_nvarlu + 1
                   pbp_nameoutlu(pbp_nvarlu) = nametem
                   write(listtem,'(a,a2,a)') trim(nametem), &
                             ': ',trim(listtemf)
                   pbp_listoutlu(pbp_nvarlu) = listtem
                   pbp_vreflu(pbp_nvarlu) = numtem
                elseif (levtem == 'stot') then
                    count = nsnow
                    do i=-nsnow+1, nsoil
                       pbp_nvarlu = pbp_nvarlu + 1
                       if (i < 1) then
                           write(luorgnametem,'(a,a1,i2.2)') trim(nametem),'_',count
                           write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                           write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Snow Level ', count
                           count = count - 1
                       else
                           write(luorgnametem,'(a,i2.2)') trim(nametem),i      
                           write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                           write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Soil Level ',i
                       endif
                       pbp_nameoutlu(pbp_nvarlu) = luorgnametem
                       pbp_listoutlu(pbp_nvarlu) = luorglisttem
                       pbp_vreflu(pbp_nvarlu) = numtem
                     enddo
                 elseif (levtem == 'soil') then
                    do i=1, nsoil
                       pbp_nvarlu = pbp_nvarlu + 1
                       write(luorgnametem,'(a,i2.2)') trim(nametem), i
                       write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                       write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Soil Level ', i
                       pbp_nameoutlu(pbp_nvarlu) = luorgnametem
                       pbp_listoutlu(pbp_nvarlu) = luorglisttem
                       pbp_vreflu(pbp_nvarlu) = numtem
                    enddo
                 endif 
            endif !dopbptem

           if (dohrtem .and. hr_saveluf) then
               if (levtem == 'sing') then
                   hr_nvarlu = hr_nvarlu + 1
                   hr_nameoutlu(hr_nvarlu) = nametem
                   write(listtem,'(a,a2,a)') trim(nametem), &
                           ': ',trim(listtemf)
                   hr_listoutlu(hr_nvarlu) = listtem
                   hr_vreflu(hr_nvarlu) = numtem
                elseif (levtem == 'stot') then
                    count=nsnow
                    do i=-nsnow+1, nsoil
                       hr_nvarlu = hr_nvarlu + 1
                       if (i < 1) then
                           write(luorgnametem,'(a,a1,i2.2)') trim(nametem),'_',count
                           write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                           write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Snow Level ', count
                           count = count - 1
                       else
                           write(luorgnametem,'(a,i2.2)') trim(nametem),i
                           write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                           write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Soil Level ',i
                       endif
                       hr_nameoutlu(hr_nvarlu) = luorgnametem
                       hr_listoutlu(hr_nvarlu) = luorglisttem
                       hr_vreflu(hr_nvarlu) = numtem
                     enddo
                 elseif (levtem == 'soil') then
                    do i=1, nsoil
                       hr_nvarlu = hr_nvarlu + 1
                       write(luorgnametem,'(a,i2.2)') trim(nametem), i
                       write(listtem,'(a,a2,a)') trim(luorgnametem),': ',trim(listtemf)
                       write(luorglisttem,'(a,a,i2.2)') trim(listtem), ' Soil Level ', i
                       hr_nameoutlu(hr_nvarlu) = luorgnametem
                       hr_listoutlu(hr_nvarlu) = luorglisttem
                       hr_vreflu(hr_nvarlu) = numtem
                    enddo
                endif
            endif  !dohrtem

        elseif (luorgtem == 'g') then
           write(listtem,'(a,a2,a)') trim(nametem), &
                    ': ',trim(listtemf)
           if (dohrtem .and. hr_savegf) then
               hr_nvarg = hr_nvarg + 1
               hr_nameoutg(hr_nvarg) = nametem
               hr_listoutg(hr_nvarg) = listtem
               hr_vrefg(hr_nvarg) = numtem
           endif  !dohrtem

           if (doqptem .and. qp_savegf) then
               qp_nvarg = qp_nvarg + 1
               qp_nameoutg(qp_nvarg) = nametem
               qp_listoutg(qp_nvarg) = listtem
               qp_vrefg(qp_nvarg) = numtem
            endif !doqptem

            if (dopbptem .and. pbp_savegf) then
               pbp_nvarg = pbp_nvarg + 1
               pbp_nameoutg(pbp_nvarg) = nametem
               pbp_listoutg(pbp_nvarg) = listtem
               pbp_vrefg(pbp_nvarg) = numtem
            endif !pbptem
        endif  !luorgtem choices
    enddo !read file

    942 continue
    close(2)

end subroutine read_outopts



