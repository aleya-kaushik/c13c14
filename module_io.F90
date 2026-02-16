module module_io

  use module_sibconst, &
        only: slen
    use kinds
    implicit none

    integer(i4), parameter :: iofnlen = 512
    real(r4), parameter :: missing = -9999.

    integer(i4) :: nchunks       ! number of processors
    integer(i4) :: rank          ! process number

    ! file path names read in from namel file
    character (len=iofnlen) ::     &
        pft_info,     & ! path to pft information
        pool_info,    & ! path to pool information
        aero_file,    & ! path to morphological/aerovar parameters
        pgdd_file,    & ! path to phenology gdd parameters
        pstg_file,    & ! path to phenology stage parameters
        phys_file,    & ! path to plant physiological parameters
        pool_file,    & ! path to pool parameters
        isodata_file, & ! path to isotope data
        vs_file,      & ! path to PFT/soil file 
        ic_file,      & ! path to initial conditions file
        driver_path,  & ! path to driver data
        tm5_path,     & ! path to tm5 data, ARACHO
        fire_path,    & ! path to fire data
        out_path,     & ! path for output files
        out_info,     & ! path to diagnostic preferences (sib_outopts)
        out_rinfo       ! path to restart preferences (sib_routopts)

    ! constant file id numbers
    integer(i4), parameter :: &
        pftid  = 23, &  !file id for pft information
        piid   = 24, &  !file id for pool information
        pgddid = 25, &  !file id for stage-phenology parameters
        pstgid = 26, &  !file id for gdd-phenology parameters
        physid = 28, &  !file id for physiological parameters
        poolid = 32, &  !file id for pool parameters
        cisoid = 34     !file id for isotope data

    !!!!!DRIVER DATA!!!!!!!!!
    ! parameters
    real(r4), parameter :: nodata = -9999.
    real(r4), parameter :: &
           ztemp = 100., &  ! temperature measurement height (m)
           zwind = 100.     ! wind measurement height (m)

    ! driver file variable names
    character(len=15), parameter :: &
        nsibname = 'nsib',        & ! number of landpoints
        psname = 'ps',            & ! pressure
        tmname = 'tm',            & ! temperature
        shname = 'sh',            & ! sensible heat
        spdmname = 'spdm',        & ! wind speed
        lsprname = 'lspr_scaled', & ! large-scale precip
        cuprname = 'cupr_scaled', & ! convective precip
        swdname = 'swd',          & ! downwelling shortwave radiation
        lwdname = 'lwd',          & ! downwelling longwave radiation
        lsprnameopt = 'lspr',     & ! optional/secondary large-scale precip
        cuprnameopt = 'cupr'        ! optional/secondary convective precip

    ! file access information
    integer(i4) :: driverid  ! netcdf id for driver data
    character (len=iofnlen) :: driver_filename

    ! flags
    logical :: driver_readf   ! read driver data ?
    logical :: driver_switchf ! switch driver data file ?
    logical :: driver_updatef ! update driver data ?

    ! time information
    integer(i4) :: driver_step     ! # seconds in driver data timestep
    integer(i4) :: driver_perday   ! # driver data timesteps per day
    integer(i4) :: driver_permon   ! # driver data timesteps per month

    integer(i4) :: driver_recnum   ! driver data current record number
    integer(i4) :: driver_year     ! year of driver data to read
    integer(i4) :: driver_month    ! month of driver data to read
    integer(i4) :: driver_day      ! day of month of driver data to read
    integer(i4) :: driver_hour     ! hour of driver data to read
    integer(long) :: driver_seccur   ! current total second of driver data
    integer(long) :: driver_secnext  ! next total second of driver data

    !!!!!TM5 DATA!!!!!!!!!    ARACHO
    ! tm5 file variable names
    character(len=15), parameter :: &
        nsibname_tm5 = 'nsib',         & ! number of landpoints
        !psname_tm5 = 'pressure',       & ! pressure
        cosmname_tm5 = 'mixing_ratio'    ! COS mixing ratio

    ! flags
    logical :: tm5_readf   ! read driver data ?
    logical :: tm5_switchf ! switch driver data file ?
    logical :: tm5_updatef ! update driver data ?

    ! time information
    integer(i4) :: tm5_step     ! # seconds in tm5 data timestep
    integer(i4) :: tm5_perday   ! # tm5 data timesteps per day
    integer(i4) :: tm5_permon   ! # tm5 data timesteps per month

    integer(i4) :: tm5_recnum   ! tm5 data current record number
    integer(i4) :: tm5_year     ! year of tm5 data to read
    integer(i4) :: tm5_month    ! month of tm5 data to read
    integer(i4) :: tm5_day      ! day of month of tm5 data to read
    integer(i4) :: tm5_hour     ! hour of tm5 data to read
    integer(long) :: tm5_seccur   ! current total second of tm5 data
    integer(long) :: tm5_secnext  ! next total second of tm5 data
    integer(i4) :: tm5id  ! netcdf id for tm5 data
    character (len=iofnlen) :: tm5_filename

    !!!!!FIRE DATA!!!!!!!!!
    ! switch to stop if fire fires are not present
    logical, parameter :: firefile_stop = .false.

    ! fire file variable names
    character(len=15), parameter :: &
        fnsibname = 'nsib',       & ! number of landpoints
        ftimename = 'time',       & ! time dimension
        firecname = 'emis_C',     & ! fire C emissions
        fireco2name = 'emis_CO2'    ! fire CO2 emissions

    ! fire file access information
    integer(i4) :: fireid
    character (len=iofnlen) :: fire_filename

    ! flags
    logical :: fire_readf   ! read fire data ?
    logical :: fire_switchf ! switch fire data file ?
    logical :: fire_updatef ! update fire data?

    ! time information
    integer(i4) :: fire_step     ! # seconds in fire data timestep
    integer(i4) :: fire_perday   ! # fire data timesteps per day
    integer(i4) :: fire_permon   ! # fire data timesteps per month

    integer(i4) :: fire_recnum   ! fire data current record number
    integer(i4) :: fire_year     ! year of fire data to read
    integer(i4) :: fire_month    ! month of fire data to read
    integer(i4) :: fire_day      ! day of fire data to read
    real(r4)    :: fire_hour     ! hour of fire data to read
    integer(long) :: fire_seccur  ! current total second of firedata
    integer(long) :: fire_secnext ! next total second of fire data
    
    !!!!!RESTART DATA!!!!!
    ! flags
    logical :: restart_savef   ! write restart files at all ?

    ! time information
    integer(i4) :: restart_dtsib  !restart output interval
                                  ! > 0 units of seconds, < 0 units of months
    integer(i4) :: restart_ntpersave ! total # of timesteps per restart output
    integer(i4) :: restart_countt    ! simulated # of timesteps between restart output

    ! file information
    character(len=iofnlen) :: restart_filename

    ! variable information
    !....number of variables
    integer(i4) :: sibr_nvar  

    !...variable names
    character(len=20), dimension(:), allocatable :: & 
          sibr_vname

    !...variable references
    integer(i4), dimension(:), allocatable :: &
          sibr_vref

    !...variable default values
    real(r4), dimension(:), allocatable :: &
          sibr_vd

    !...output flags (for restart and requib)
    logical, dimension(:), allocatable :: &
          sibr_doref

    !!!!!REQUIB DATA!!!!!
    logical :: requib_writef   ! calculate and write equilibrium pools ?

    ! file information
    character(len=iofnlen) :: requib_filename

    !...output flags
    logical, dimension(:), allocatable :: &
         sibreq_doref



    !!!!!!!!!!!!!!!!!!!!!!!OUTPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  HR=Global Time-Series Output (typically hourly)              !!
    !!  PBP=Selected Point Time-Series Output                        !!
    !!       - Typically hourly for global runs with selected sites  !!
    !!       - Typically daily for site/single point runs            !!
    !!  QP=Global Time-Series Output (typically monthly)             !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! parameters
    integer, parameter :: dtoutwrite = 86400

    ! filename information
    logical :: printrank = .true.
    character(len=5), parameter :: &
        hr_prefix = 'hsib_', &
        pbp_prefix = 'psib_', &
        qp_prefix = 'qsib_'

    character(len=12), parameter :: &
        hr_suffixg  = '.g.hr2.nc', &
        hr_suffixlu = '.lu.hr3.nc', &
        pbp_suffixg  = '.g.pbp1.nc', &
        pbp_suffixlu = '.lu.pbp2.nc', &
        qp_suffixg  = '.g.qp2.nc',  &
        qp_suffixlu = '.lu.qp3.nc'

    ! file names
    character (len=iofnlen) :: &
        hr_filenameg,   & ! filename for hr gridcell files
        hr_filenamelu,  & ! filename for hr land unit files
        pbp_filenameg,  & ! filename for pbp gridcell files
        pbp_filenamelu, & ! filename for pbp land unit files
        qp_filenameg,   & ! filename for qp gridcell files
        qp_filenamelu     ! filename for qp land unit files

    ! flags
    logical  ::      &
        hr_savegf,   &  ! flag to save hr g info
        hr_saveluf,  &  ! flag to save hr lu info
        pbp_savegf,  &  ! flag to save pbp g info
        pbp_saveluf, &  ! flag to save pbp lu info
        qp_savegf,   &  ! flag to save qp g info
        qp_saveluf      ! flag to save qp lu info

    logical  ::     &
        hr_openf,   & ! open hourly files?
        hr_writef,  & ! write to hourly files?
        pbp_openf,  & ! open pbp files?
        pbp_writef, & ! write to pbp files?
        qp_openf,   & ! open qp files?
        qp_writef     ! write to qp files?

    ! time information

    ! > 0 units of seconds; < 0 units of months
    integer(i4) :: & ! output intervals (from namel file)                      
        hr_dtsib,  & ! hr output interval 
        pbp_dtsib, & ! pbp output interval
        qp_dtsib     ! qp output interval

    integer(i4) :: &  ! number of timesteps per file
        hr_nsaveperfile, qp_nsaveperfile, &
        pbp_nsaveperfile

    integer(i4) :: &  ! number of saved steps per output
        hr_nsaveperout, qp_nsaveperout, pbp_nsaveperout
 
    integer(i4) :: &  ! number of timesteps per save
        hr_ntpersave, qp_ntpersave, pbp_ntpersave
    real(r8) ::    &  ! weight of timestep output to saved output
        hr_wtpersave, qp_wtpersave, pbp_wtpersave

    real(r8) ::  &    ! (units of days)
         hr_start,  & ! day of hr starting period
         pbp_start, & ! day of pbp starting period
         qp_start     ! day of qp starting period

    real(r8) ::     &   ! (units of day fractions)
         hr_step,   &   ! day fraction of hourly timestep
         pbp_step,  &   ! day fraction of pbp timestep
         qp_step        ! day fraction of qp timestep

    ! variable information
    integer(i4) ::  &
        hr_nvarg,    & ! number of saved hrsib fields on Grid Cell
        hr_nvarlu,   & ! number of saved hrsib fields on Land Units
        pbp_nvarg,   & ! number of saved pbp fields on Grid Cell
        pbp_nvarlu,  & ! number of saved pbp fields on Land Units
        qp_nvarg,    & ! number of saved qpsib fields on Grid Cell
        qp_nvarlu      ! number of saved qpsib fields on Land Units

    integer(i4), dimension(:), allocatable ::    & !(hr/pbp/qp nvar)
        hr_vrefg,   & ! reference # for saved grid cell hr variables
        hr_vreflu,  & ! reference # for hrsib fields on Land Units to be saved
        pbp_vrefg,  & ! reference # for pbpsib fields on Grid Cell to be saved
        pbp_vreflu, & ! reference # for pbpsib fields on Land Units to be saved
        qp_vrefg,   & ! reference # for qpsib fields on Grid Cell to be saved
        qp_vreflu     ! reference # for qpsib fields on Land Units to be saved

    character(len=100), dimension(:), allocatable ::  & !(hr/pbp/qp nvar)
        hr_listoutg,   & ! descriptions of hr variables on Grid Cell
        hr_listoutlu,  & ! descriptions of hr variables on Land Unitse
        pbp_listoutg,  & ! descriptions of pbp variables on Grid Cell
        pbp_listoutlu, & ! descriptions of pbp variables on Land Units
        qp_listoutg,   & ! descriptions of qp variables on Grid Cell
        qp_listoutlu     ! descriptions of qp variables on Land Units

    character(len=21), dimension(:), allocatable ::    & !(hr/pbp/qp nvar)
        hr_nameoutg,    & ! field names of hr variables on Grid Cell
        hr_nameoutlu,   & ! field names of hr variables on Land Units
        pbp_nameoutg,   & ! field names of pbp variables on Grid Cell
        pbp_nameoutlu,  & ! field names of pbp variables on Land Units
        qp_nameoutg,    & ! field names of qp variables on Grid Cell
        qp_nameoutlu      ! field names of qp variables on Land Units

    ! output information
    integer(i4) :: &  ! file ids
        hr_idg, hr_idlu, pbp_idg, pbp_idlu, &
        qp_idg, qp_idlu

    integer(i4) :: &  ! time ids
        hr_idtimeg, hr_idtimelu, pbp_idtimeg, pbp_idtimelu, &
        qp_idtimeg, qp_idtimelu

    integer(i4), dimension(:), allocatable ::  & !(nvars)
        hr_idvarlu,  hr_idvarg,  & !variable ids for hr variables
        pbp_idvarlu, pbp_idvarg, & !variable ids for pbp variables
        qp_idvarlu,  qp_idvarg     !variable ids for qp variables

    integer(i4) :: &  ! counter information
        hr_countt, &  ! count of timestep being saved per output interval
        hr_countsave,  & ! count of saved value per output interval
        hr_countwrite, & ! count of timestep being written out
        pbp_countt, pbp_countsave, pbp_countwrite, &
        qp_countt, qp_countsave, qp_countwrite

    real(r8), dimension(:,:,:), allocatable ::  & !(subcount,nvar,ntperout)
        hr_g,  & ! time series diagnostics all points (typically hourly)
        qp_g,  & ! time series diagnostics all points (typically monthly)
        pbp_g    ! time series diagnostics at selected points

    real(r8), dimension(:,:,:,:), allocatable ::  & !(subcount,nlu,nvar,ntperout)
        hr_lu,   & ! time-series diagnostics all points (typically hourly)
        qp_lu,   & ! time-series diagnostics all points (typically monthly)
        pbp_lu     ! time-series diagnostics at selected points


    !!!!!PBP-Specific Information
    integer(i4) :: npbp  ! number of gridpoints where pbp fields are saved

    integer(i4), dimension(:), allocatable :: & !(subcount)
        pbp_outref

    integer(i4), dimension(:), allocatable ::  & !(npbp)
        pbp_sibpt,  &  ! original sibpt references
        pbp_gref       ! references for PBP output

    real(r4), dimension(:,:), allocatable :: &  !(npbp,npbp)
       pbp_lonlat  !lon and lat of pbp points from namelist

    real(r4), dimension(:), allocatable ::     & !(npbp)
       pbp_lon, & !  Longitudes for PBP output
       pbp_lat    !  Latitudes for PBP output

    integer(i4), dimension(:,:), allocatable :: &  !(npbp,nlu)
        pbp_pref    ! land unit PFT refs for PBP output
    real(r4), dimension(:,:), allocatable ::    &  !(npbp,nlu)
        pbp_larea   ! land unit areas for PBP output

    character(len=slen), dimension(:), allocatable :: & !(npbp)
        pbp_sitename

    !!!!!Satellite Information - Fluorescence
    !Output times for GOME-2 and OCO-2
    real(r4), parameter :: sif_gome2_tstart = 0.375
    real(r4), parameter :: sif_gome2_tstop = 0.4585
    real(r4), parameter :: sif_oco2_tstart = 0.5104
    real(r4), parameter :: sif_oco2_tstop = 0.5521

    ! Count of output times for satellite SIF comparisons
    ! 1=GOME-2  2=OCO-2
    integer(i4), dimension(:,:,:), allocatable :: & !(subcount,ntperout,2)
         hr_sif_satcount,  &
         pbp_sif_satcount, &
         qp_sif_satcount

end module module_io
