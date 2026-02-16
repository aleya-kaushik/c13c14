#include "nc_util.h"

!--------------------------------------------
!Includes routines to create new output files
!   for both global and land unit output.
!---------------------------------------------

subroutine diagnostic_createg( &
    year, curday_per_year,     &
    filename, numvars, subcount, &
    lonsib, latsib, sitesib, &
    period_start, period_step, &
    nsaveperfile, &
    namevars, listvars, &
    fileid, varids, timeid)

    use module_io, only: iofnlen
    use module_sibconst, only: slen
    use netcdf
    use kinds

    implicit none

    ! input variables
    integer(i4), intent(in) :: year            !current year of simulation
    integer(i4), intent(in) :: curday_per_year !number of days in current year
    character(len=iofnlen), intent(in) :: filename   ! filename for output
    integer(i4), intent(in) :: numvars           ! max # of variables written to file
    integer(i4), intent(in) :: subcount          ! number of landpoints on subset
    real(r4), dimension(subcount), intent(in) :: lonsib   ! array of longitude coordinates
    real(r4), dimension(subcount), intent(in) :: latsib   ! array of latitude coordinates
    character(len=slen), dimension(subcount), intent(in) :: sitesib  ! site names
    real(r8), intent(in) :: period_start    ! start of averaged period
    real(r8), intent(in) :: period_step     ! timestep for averaged period
    integer(i4), intent(in) :: nsaveperfile ! number of times saved per file
    character(len=21), dimension(numvars), intent(in) :: namevars   ! names of vars
    character(len=100), dimension(numvars), intent(in) :: listvars  ! variable descriptions
    integer(i4), intent(out) :: fileid           ! file id#
    integer(i4), dimension(numvars), intent(out):: varids  ! variable id#s 
    integer(i4), intent(out) :: timeid           ! time variable id#

    ! netcdf variables
    integer(i4) :: nscid   ! landpoints dimension id #
    integer(i4) :: slenid  ! site name string length dimension id #
    integer(i4) :: timedid ! time dimension id #

    integer(i4) :: latsibid ! latitude variable id #
    integer(i4) :: lonsibid ! longitude variable id #
    integer(i4) :: snameid  ! site name variable id #

    ! local variables
    logical :: savesite
    integer(i4) :: t, nv, syeari
    integer(i4) :: status          ! netcdf status
    character(len=4)  :: syear     !year string
    character(len=40) :: units          ! variable units
    character(len=80) :: longname       ! variable description
    character(len=60) :: titlename      ! variable title

    integer(i4) :: step
    real(r8) :: ddoy, dtime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! determine if saving site names
    if (sitesib(1) .ne. '') then
       savesite = .true.
    else
       savesite = .false.
    endif

    ! make sure fileid is not tied to any open file
    status = nf90_close( fileid )
    CHECK( nf90_create( trim(filename), nf90_clobber, fileid) )

    ! define dimensions
    CHECK( nf90_def_dim( fileid, 'landpoints', subcount, nscid ) )
    CHECK( nf90_def_dim( fileid, 'time', nsaveperfile, timedid ) )
    IF (savesite) CHECK( nf90_def_dim( fileid, 'slen', slen, slenid ) )    

    ! define control variables
    CHECK( nf90_def_var( fileid, 'lonsib', nf90_float, (/nscid/), lonsibid ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'long_name', 'SiB Point Longitude' ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'units', 'degrees_east' ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'quantity', 'longitude' ) )

    CHECK( nf90_def_var( fileid, 'latsib', nf90_float, (/nscid/), latsibid ) )
    CHECK( nf90_put_att( fileid, latsibid, 'long_name', 'SiB Point Latitude' ) )
    CHECK( nf90_put_att( fileid, latsibid, 'units', 'degrees_north' ) )
    CHECK( nf90_put_att( fileid, latsibid, 'quantity', 'latitude' ) )

    IF (sitesib(1) .NE. '') THEN
       CHECK( nf90_def_var( fileid, 'site_names', nf90_char, (/slenid,nscid/), snameid ) )
    ENDIF

    CHECK( nf90_def_var( fileid, 'time', nf90_double, (/timedid/), timeid ) )
    CHECK( nf90_put_att( fileid, timeid, 'quantity', 'time' ) )
    syeari = index(trim(filename),'_',back=.true.)
    syeari = syeari + 1
    write(syear,'(a4)') filename(syeari:syeari+3)
    CHECK( nf90_put_att( fileid, timeid, 'units', 'days since '//syear//'-01-01'  ) )
    
    ! define data variables
    do nv = 1, numvars
       status = nf90_def_var(fileid,trim(namevars(nv)),nf90_double,(/nscid,timedid/),varids(nv))
       if (status .ne. nf90_noerr) then
          print*,'Error defining variables for diagnostic_createg.F90. Stopping.'
          stop
       endif
       call get_units( listvars(nv), longname, units )
       call get_title( longname, titlename )

       CHECK( nf90_put_att( fileid, varids(nv), 'long_name', trim(longname) ) )
       CHECK( nf90_put_att( fileid, varids(nv), 'title', trim(titlename) ) )
       CHECK( nf90_put_att( fileid, varids(nv), 'units', trim(units) ) )
    enddo

    ! switch from definition mode to data mode
    CHECK( nf90_enddef( fileid ) )

    ! assign values to variables not variant with time
    CHECK( nf90_put_var( fileid, lonsibid, lonsib ) )
    CHECK( nf90_put_var( fileid, latsibid, latsib ) )
    IF (savesite) CHECK( nf90_put_var( fileid, snameid, sitesib) )

    ! write time variables
    ddoy = dble(period_start) + (period_step * 0.5)
    step=0
    do t=1,nsaveperfile
       dtime = year + ddoy / dble(curday_per_year)
       step = step + 1
       CHECK( nf90_put_var( fileid, timeid, ddoy, (/step/) ) )
       ddoy = ddoy + period_step
   enddo
    
end subroutine diagnostic_createg

!------------------------------------------------------------

subroutine diagnostic_createlu( &
    year, cdaypyr, &
    filename, numvars, subcount, &
    lonsib, latsib, sitesib, &
    lgarea, lgref, &
    period_start, period_step, &
    nsaveperfile,  &
    namevars, listvars, & 
    fileid, varids, timeid)

    use netcdf
    use kinds
    use module_io, only: iofnlen
    use module_pftinfo, only: &
         clen, pft_name, pft_ref
    use module_sibconst, only: &
         nlu, npft, slen

    implicit none

    ! input variables
    integer(i4), intent(in) :: year     !current year of simulation
    integer(i4), intent(in) :: cdaypyr  !number of days in current year
    character(len=iofnlen), intent(in) :: filename
    integer(i4), intent(in) :: numvars     ! # of variables to write out
    integer(i4), intent(in) :: subcount    ! # of landpoints simulated
    real(r4), dimension(subcount), intent(in) :: lonsib   ! array of longitude coordinates
    real(r4), dimension(subcount), intent(in) :: latsib   ! array of latitude coordinates
    character(len=slen), dimension(subcount), intent(in) :: sitesib ! site names of landpoints
    real(r4), dimension(subcount,nlu), intent(in) :: lgarea    ! % areal coverage of lu/gridcell
    integer(i4), dimension(subcount,nlu), intent(in) :: lgref  !reference of lu/gridcell 
    real(r8), intent(in) :: period_start    ! start of averaged period
    real(r8), intent(in) :: period_step     ! timestep for averaged period
    integer(i4), intent(in) :: nsaveperfile ! number of times saved per file
    character(len=21), dimension(numvars), intent(in) ::  namevars  ! variable names
    character(len=100), dimension(numvars), intent(in) :: listvars ! variable descriptions
    integer(i4), intent(out) :: fileid     ! file id #
    integer(i4), dimension(numvars), intent(out) :: varids  ! variable id #s
    integer(i4), intent(out) :: timeid       ! time variable id #

    ! netcdf variables
    integer(i4) :: nscid  ! landpoints dimension id #
    integer(i4) :: nluid  ! nlu dimension id #
    integer(i4) :: npftid ! npft dimension id #
    integer(i4) :: clenid ! clen dimension id #
    integer(i4) :: slenid ! slen dimension id #
    integer(i4) :: tdid   ! time dimension id #
    integer(i4) :: pnameid, prefid    ! pft name and reference id #
    integer(i4) :: lonsibid, latsibid ! landpoints lat/lon variable id #
    integer(i4) :: snameid ! sitename variable id #
    integer(i4) :: areaid  ! percent coverage of gridcell
    integer(i4) :: rupid   ! PFT reference number

    ! local variables
    logical :: savesite   ! flag for writing out site names
    integer(i4) :: n, t, syeari      ! index variable
    integer(i4) :: status            ! status of netcdf 
    character(len=4)  :: syear       ! year string
    character(len=40) :: units       ! variable units
    character(len=80) :: longname    ! variable description
    character(len=60) :: titlename   ! variable title

    integer(i4) :: step
    real(r8) :: ddoy, dtime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! determine if saving site names
    if (sitesib(1) .ne. '') then
       savesite = .true.
    else
       savesite = .false.
    endif

    ! make sure fileid is not tied to any open file
    status = nf90_close( fileid )
    CHECK( nf90_create( trim(filename), nf90_clobber, fileid) )

    ! define dimensions
    CHECK( nf90_def_dim( fileid, 'landpoints', subcount, nscid ) )
    CHECK( nf90_def_dim( fileid, 'nlu', nlu, nluid ) )
    CHECK( nf90_def_dim( fileid, 'npft', npft, npftid ) )    
    CHECK( nf90_def_dim( fileid, 'time', nsaveperfile, tdid ) )
    CHECK( nf90_def_dim( fileid, 'clen', clen, clenid ) )
    IF (savesite) CHECK( nf90_def_dim( fileid, 'slen', slen, slenid ) )  
    
    ! define control variables
    CHECK( nf90_def_var( fileid, 'lonsib', nf90_float, (/nscid/), lonsibid ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'long_name', 'SiB Point Longitude' ) )
    CHECK( nf90_put_att( fileid, lonsibid, 'units', 'degrees_east' ) )

    CHECK( nf90_def_var( fileid, 'latsib', nf90_float, (/nscid/), latsibid ) )
    CHECK( nf90_put_att( fileid, latsibid, 'long_name', 'SiB Point Latitude' ) )
    CHECK( nf90_put_att( fileid, latsibid, 'units', 'degrees_north' ) )

    IF (sitesib(1) .NE. '') THEN
       CHECK( nf90_def_var( fileid, 'site_names', nf90_char, (/slenid,nscid/), snameid ) )
    ENDIF

    CHECK( nf90_def_var( fileid, 'pft_names', nf90_char, (/clenid,npftid/), pnameid ) )
    CHECK( nf90_put_att( fileid, pnameid, 'long_name', 'PFT Name'))

    CHECK( nf90_def_var( fileid, 'pft_refnums', nf90_byte, (/npftid/), prefid ) )
    CHECK( nf90_put_att( fileid, prefid, 'long_name', 'PFT Reference Numbers' ) )
    
    CHECK( nf90_def_var( fileid, 'lu_area', nf90_float, (/nscid, nluid/), areaid ) )
    CHECK( nf90_put_att( fileid, areaid, 'long_name', 'Land Unit Areal Coverage' ) )
    CHECK ( nf90_put_att( fileid, areaid, 'units', 'Fraction (0-1)' ) )

    CHECK( nf90_def_var( fileid, 'lu_pref', nf90_int, (/nscid, nluid/), rupid ) )
    CHECK( nf90_put_att( fileid, rupid, 'long_name', 'Land Unit PFT Reference' ) )
    CHECK( nf90_put_att( fileid, rupid, 'units', 'Number') )
    CHECK( nf90_put_att( fileid, rupid, 'min_val', 1 ) )
    CHECK( nf90_put_att( fileid, rupid, 'max_val', maxval(pft_ref) ) )
    
    CHECK( nf90_def_var( fileid, 'time', nf90_double, (/tdid/), timeid ) )
    CHECK( nf90_put_att( fileid, timeid, 'quantity', 'time' ) )
    syeari = index(trim(filename),'_',back=.true.)
    syeari = syeari + 1
    write(syear,'(a4)') filename(syeari:syeari+3)
    CHECK( nf90_put_att( fileid, timeid, 'units', 'days since '//syear//'-01-01'  ) )

    ! define data variables
    do n = 1, numvars
        status = nf90_def_var(fileid,trim(namevars(n)),nf90_double,(/nscid,nluid,tdid/),varids(n))
        if (status .ne. nf90_noerr) then
           print*,'Error defining variables for diagnostic_createlu.F90. Stopping.'
           stop
        endif
        call get_units( listvars(n), longname, units )
        call get_title( longname, titlename )

        CHECK( nf90_put_att( fileid, varids(n), 'long_name', trim(longname) ) )
        CHECK( nf90_put_att( fileid, varids(n), 'title', trim(titlename) ) )
        CHECK( nf90_put_att( fileid, varids(n), 'units', trim(units) ) )
    enddo

    ! switch from definition mode to data mode
    CHECK( nf90_enddef( fileid ) )

    ! assign values to variables not variant with time
    CHECK( nf90_put_var( fileid, latsibid, latsib ) )
    CHECK( nf90_put_var( fileid, lonsibid, lonsib ) )
    IF (savesite) CHECK( nf90_put_var( fileid, snameid, sitesib) )
    CHECK( nf90_put_var( fileid, pnameid, pft_name) )
    CHECK( nf90_put_var( fileid, prefid, pft_ref) )
    CHECK( nf90_put_var( fileid, areaid, lgarea ) )
    CHECK( nf90_put_var( fileid, rupid, lgref ) )

    ! write time variables
    ddoy = dble(period_start) + (period_step * 0.5)
    step=0
    do t=1,nsaveperfile
       dtime = year + ddoy / dble(cdaypyr)
       step = step + 1
       CHECK( nf90_put_var( fileid, timeid, ddoy, (/step/) ) )
       ddoy = ddoy + period_step
    enddo
    
end subroutine diagnostic_createlu

!-----------------------------------------------------------------------
subroutine diagnostic_close(ncid)

    use netcdf
    use kinds

    implicit none

    ! input variables
    integer(i4), intent(inout) :: ncid

    ! local variables
    integer(i4) :: status

    status = nf90_noerr
    IF (ncid > 0) status = nf90_close(ncid)

    IF (status /= nf90_noerr) print*,'Error Closing Diagnostic Output'
    ncid = 0

end subroutine diagnostic_close

!-----------------------------------------------------------------------
subroutine get_units(description, longname, units)

    !-----------------------------------------------------------------------
    !  Purpose:
    !   Extracts a string enclosed within parentheses. 
    !   Used for units which are contained in a general description string.
    !
    !   Returns the units string (units) and
    !   the description string with the units removed (longname).
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: description
    character(len=*), intent(out) :: units
    character(len=*), intent(out) :: longname

    integer :: n, start_paren, end_paren, paren_count

    paren_count = 0
    start_paren = len_trim(description)
    end_paren = len_trim(description)

    do n = len(description), 1, -1
        if (description(n:n)==")") then
            if (paren_count == 0) then
                end_paren = n
            endif
            paren_count = paren_count + 1
        else if (description(n:n) == "(") then
            paren_count = paren_count - 1
            if (paren_count == 0) then
                start_paren = n
                exit
            endif
        end if
    end do

    !   in case of confusion, clear units and return unaltered description
    !   note: start_paren > end_paren should not be possible, but just in case...
    !         start_paren = end_paren occurs when there are no units
    !         start_paren = end_paren-1 occurs when units are "()"
    !   FIXME: n==1 is too limiting - what if I wanted only units?
    if (n == 1 .or. start_paren >= end_paren) then   ! no units
        units = " "
        longname = trim(description)
    else if (start_paren == (end_paren-1)) then      ! "()" case
        units = " "
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
    else                                             ! normal units
        units = description(start_paren+1:end_paren-1)
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
    end if

end subroutine get_units

!-----------------------------------------------------------------------
subroutine get_title(description, titlename)

    !-----------------------------------------------------------------------
    !  Purpose:
    !   Extracts the string following a colon.
    !
    !   Returns the title string (titlename).
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: description
    character(len=*), intent(out) :: titlename

    integer :: n, colon_pos, end_pos

    colon_pos = 0
    end_pos = len_trim(description)

    do n = len(description), 1, -1
        if (description(n:n)==":") then
           colon_pos = n
        end if
    end do

    if (colon_pos == 0) then
        titlename = trim(description)
    else
        titlename = trim(description(colon_pos+2:end_pos))
    end if

end subroutine get_title
