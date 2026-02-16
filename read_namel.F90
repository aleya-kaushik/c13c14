
! Reads in sibdrv control variables.
subroutine read_namel( namel_name )

use module_sibconst
use module_io
use module_time

implicit none

!...input variables
character(len=200) :: namel_name

!...local variables
character(len=256) :: dr_path, fr_path
!character(len=512) :: icfile

!...NAMELISTS
namelist /CONTROL_LIST_SIBDRV/ & ! USER SELECTED TIMES
    nsib, starttime, startyear, endtime, endyear, &
    dtsib, restart_dtsib, qp_dtsib, pbp_dtsib, hr_dtsib
namelist /IO_LIST_SIBDRV/ & ! USER SELETED I/O OPTIONS
    pft_info, pool_info, &
    aero_file, pgdd_file, pstg_file, &
    phys_file, pool_file, isodata_file, &
    vs_file, ic_file, dr_path, tm5_path, fr_path, &
    out_path, out_info, out_rinfo
namelist /SPINUP_LIST_SIBDRV/ & ! USER SELECTED SPINUP OPTIONS
    spinup, spinup_default, spinup_continue, spinup_numyrs, spinup_maxiter, &
    spinup_threshold, spinup_writediag, spinup_writetxtf
namelist /SUBGRID_SIBDRV/ & ! SUBGRID INFO FOR NON-POINT RUNS
    minlon, maxlon, minlat, maxlat
namelist /PBP_LIST_SIBDRV/ & ! USER DEFINED PBP DIAGNOSTIC LOCATIONS
    npbp
namelist /BALAN_LIST_SIBDRV/ & ! USER SELECTED BALANCE CHECKING OPTIONS
    badtc_print, badtc_stop, &
    canb_print, canb_stop, canb_thresh, &
    carbonb_print, carbonb_stop, carbonb_thresh, &
    carbonb_threshc13, &
    fireb_print, fireb_stop, fireb_thresh, &
    snocmbn_print, snocmbn_stop, snocmbn_thresh, &
    bnum_allow, energyb_print, energyb_stop, energyb_thresh, &
    waterb_print, waterb_stop, waterb_thresh
namelist /PRINT_LIST_SIBDRV/ & ! USER SELECTED PRINTING OPTIONS
    print_avec, print_driver, print_fire, print_harvest, &
    print_pftinfo, print_pooll, print_soil, print_sscol, &
    print_veg, print_stop
namelist /SWITCH_LIST_SIBDRV/ & ! USER SELECTED SWITCHES
    cornsoy_switch, fire_switch, grazing_switch, green_switch, &
    eqclear_switch, leapyr_switch, updatelst_switch, tm5mr_switch, &
    soilogee_switch, varco2_switch, varc13_switch, varc13m_switch

!print*,'icfile:',icfile
    !-----------------------------------------------------------------------
    ! read in namel_sibdrv
    !-----------------------------------------------------------------------
    print*,''
    print*,'Reading in namelist file: ',trim(namel_name)
    open(unit=2,file=trim(namel_name),form='formatted')
    read (2,CONTROL_LIST_SIBDRV)
    read (2,IO_LIST_SIBDRV)
    read (2,SPINUP_LIST_SIBDRV)
    read (2,SUBGRID_SIBDRV)
    read (2,PBP_LIST_SIBDRV)
    IF (npbp .GT. 0) THEN
       allocate (pbp_lonlat(2,npbp))
       pbp_lonlat = 0.0
       read(2,*,err=919)pbp_lonlat
       919  continue
    ELSE
       allocate(pbp_lonlat(2,nsib))
       
    ENDIF

    read (2,BALAN_LIST_SIBDRV)
    read (2,PRINT_LIST_SIBDRV)
    read (2,SWITCH_LIST_SIBDRV)

    close(2)

    print ('(a,2i6)'), '   Starting Time (DOY/Year): ',starttime, startyear
    print ('(a,2i6)'), '   Ending Time (DOY/Year): ',endtime, endyear
    print ('(a,i6)'),  '   SiB time step (s) = ',dtsib

    driver_path = dr_path
    fire_path = fr_path
!    ic_file = icfile       

    if ((spinup) .and. (.not. spinup_writediag)) then
       restart_savef = .false.
       hr_savegf   = .false.
       hr_saveluf  = .false.
       qp_savegf   = .false.
       qp_saveluf  = .false.
       pbp_savegf  = .false.
       pbp_saveluf = .false.
    else
       restart_savef = restart_dtsib /= 0
       hr_savegf  = (hr_dtsib /= 0)
       hr_saveluf = hr_savegf
       qp_savegf  = (qp_dtsib /= 0)
       qp_saveluf = qp_savegf
       pbp_savegf = (pbp_dtsib /= 0)
       pbp_saveluf = pbp_savegf
    endif

    single_pt = .false.
    if (nsib == 1) then
        single_pt=.true. 
        printrank = .false.
    endif

end subroutine read_namel
