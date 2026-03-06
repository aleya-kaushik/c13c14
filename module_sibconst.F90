module module_sibconst

!----------------------------------------------------------------------
!
!   SiB4 Constant Module
!
!----------------------------------------------------------------------

    use kinds
    implicit none

    !--------------------------------
    !User specific variables for SiB
    real (r4), parameter  ::  version = 4.0 ! code version identifier
    character(len=10), parameter :: sib_source = 'SiB4'
    character(len=10), parameter :: grid_source = 'lat/lon'

    !...Land Unit Information
    integer(i4), parameter :: nlu = 10      ! # of land units
    integer(i4), parameter :: nsoil = 10    ! # of soil layers 
    integer(i4), parameter :: nsoiltop = 3  ! # of top soil layers
    integer(i4), parameter :: nsnow = 5     ! max snow layers
    integer(i4), parameter :: ntot = nsoil+nsnow  ! total soil column layers

    !...Phenology/Pool Information
    integer(i4), parameter :: dtpool = 86400 ! # seconds to update

    !...Site Information
    integer(i4), parameter :: slen=6 !site name length

    !---------------------------------
    !...Isotope data from c_iso_time_series.dat
    integer(i4) :: nisodatayr !number of years of isotope data read in

    !-------------------
    !...PFT Information
    integer(i4) :: npft   ! number of PFTs
    integer(i4) :: ntype  ! number of PFT types
    integer(i4) :: ngroup ! number of PFT groups
    integer(i4) :: npmeth ! number of phenology methods
    integer(i4) :: npstgmax ! maximum number of phenology stages

    !---------------------
    !...Pool Information
    integer(i4) :: npoolpft ! number of live carbon pools with PFT traits
    integer(i4) :: npoollu  ! number of dead/soil carbon pools with LU traits
    integer(i4) :: ntpool   ! total number of pools (npoolpft+npoollu)

    integer(i4) :: npoolcan   ! number of live carbon pools in canopy
    integer(i4) :: npoolsfc   ! number of dead carbon pools at surface
    integer(i4) :: npoolsoil  ! number of dead carbon pools in soil
    !...same as above but for C13 pools
    integer(i4) :: npoolcanc13   ! number of live carbon pools in canopy
    integer(i4) :: npoolsfcc13   ! number of dead carbon pools at surface
    integer(i4) :: npoolsoilc13  ! number of dead carbon pools in soil

    !...same as above but for C14 pools
    integer(i4) :: npoolcanc14   ! number of live carbon pools in canopy
    integer(i4) :: npoolsfcc14   ! number of dead carbon pools at surface
    integer(i4) :: npoolsoilc14  ! number of dead carbon pools in soil

    !integer(i4) :: cantot
 !   integer(i4) :: npoolpft_c13 ! number of live carbon pools with PFT traits
 !   integer(i4) :: npoollu_c13  ! number of dead/soil carbon pools with LU traits
 !   integer(i4) :: ntpool_c13   ! total number of pools (npoolpft+npoollu)

 !   integer(i4) :: npoolcan_c13   ! number of live carbon pools in canopy
 !   integer(i4) :: npoolsfc_c13   ! number of dead carbon pools at surface
 !   integer(i4) :: npoolsoil_c13  ! number of dead carbon pools in soil
    
    !-------------------------
    !...Control Switches
    logical  :: cornsoy_switch   !flag: alternate between corn and soybeans
    logical  :: fire_switch      !flag: fire emissions
    logical  :: grazing_switch   !flag: grazing
    logical  :: green_switch     !flag: use greenness fraction
    logical  :: eqclear_switch   !flag: set equilibrium calculation sums
                                 !      to zero at simulation start
    logical  :: leapyr_switch    !flag: use leap years or constant 365 days per year
    logical  :: updatelst_switch !flag: use local standard time for phenology updates
    logical  :: tm5mr_switch     !flag: use tm5 COS mole fraction
    logical  :: soilogee_switch  !flag: use COS soil from Ogee et al., 2016
    logical  :: varco2_switch    !flag: use variable CO2 mole fraction
    logical  :: varciso_switch    !flag: use variable kiecps/d13cca C13 mole fraction    
    logical  :: varcisom_switch    !flag: use variable atm C13 and C14 mole fraction
    !-----------------------
    !...Lat/Lon Information
    logical :: single_pt
    integer(i4) ::     &
        nsib,           & !  number of SiB points in datasets
        subcount          !  actual number of SiB points in simulation
        
    real(r4), dimension(:), allocatable :: &  !(nsib)
        latsib,       & !  SiB point latitude
        lonsib          !  SiB point longitude

    character(len=slen), dimension(:), allocatable :: &  !(nsib)
        sitenamesib     ! SiB point site name

    ! grid variables from namel file - subdomain limits
    real(r4) :: minlon  !minimum longitude
    real(r4) :: maxlon  !maximum longitude
    real(r4) :: minlat  !minimum latitude
    real(r4) :: maxlat  !maximum latitude

    integer(i4), dimension(:), allocatable :: & !(subcount)
        subset   !  array of landpoint indices for subgrid based on nsib

    real(r4), dimension(:), allocatable :: & !(subcount)
        sublatsib,    & !  latitude of subset
        sublonsib       !  longitude of subset

    real(r4), dimension(:,:), allocatable :: & !(subcount, nlvc)
        sublarea  !land unit fractional areas
    integer(i4), dimension(:,:), allocatable :: & !(subcount, nlvc)
        subpref    !PFT references per lvc unit

    character(len=slen), dimension(:), allocatable :: & !(subcount)
        subsitename   !Site Name 

    !----------------------------------
    !...Solar Zenith Angle Information
    real(r8) sin_dec  ! (-) sin solar declination
    real(r8) cos_dec  ! (-) cosine solar declination
    real(r8) tan_dec  ! (-) tangent solar declination
    real(r8) :: lonearth ! (rad) Earth lon about Sun from vernal equinox

   !---------------------------------------
   !...Spinup/Equilibrium Pool Information
   logical :: spinup     ! spinup simulation if set to .true.
   logical :: spinup_default ! uses default initial conditions if set to .true.
                             !    => all pools zero except minimal froot pool
                             !    => soil moisture starts at saturation
   logical :: spinup_continue !spinup but start from a restart file  
   integer(i4) :: spinup_numyrs  ! # of years in equilibrium pool calculation
   integer(i4) :: spinup_maxiter ! maximum number of spinup iterations 
   real(r4)    :: spinup_threshold  ! threshold to determine pool equilibrium (%)
   logical :: spinup_writediag  ! write diagnostic output files?
   logical :: spinup_writetxtf  ! write equilibrium pools to a text file?

   logical :: spinup_done  ! flag for determining if further spin-up is required
   integer(i4) :: spinup_lnum !  count of spinup iterations
   integer(i4) :: spinup_ynum !  count of years in spinup
   integer(i4) :: spinup_prevnum !  count of spinup iterations before restart

   !---------------------------------------
   !...Balance Information
   logical  :: badtc_print    !print canopy temperatures?
   logical  :: badtc_stop     !stop for bad canopy temperatures?
   integer(i4) :: bnum_allow  !number of allowable balance offenses
   logical  :: canb_print     !print canopy balance values?
   logical  :: canb_stop      !stop if canopy balance fails?
   real(r4) :: canb_thresh    !canopy balance threshold
   logical  :: carbonb_print  !print carbon balance values?
   logical  :: carbonb_stop   !stop if carbon balance fails?
   real(r4) :: carbonb_thresh !carbon balance threshold 
   real(r4) :: carbonb_threshc13 !carbon balance threshold for C13
   logical  :: fireb_print    !print fire balance values?
   logical  :: fireb_stop     !stop if fire balance fails?
   real(r4) :: fireb_thresh   !fire balance threshold
   logical  :: snocmbn_print  !print snow combine information?
   logical  :: snocmbn_stop   !stop if snow combine water balance fails?
   real(r4) :: snocmbn_thresh !snow combine balance threshold
   logical  :: energyb_print  !print energy balance values?
   logical  :: energyb_stop   !stop if energy balance fails?
   real(r4) :: energyb_thresh !energy balance threshold
   logical  :: waterb_print   !print water balance values?
   logical  :: waterb_stop    !stop if water balance fails?
   real(r4) :: waterb_thresh  !water balance threshold 

   !---------------------------------------
   !...Printing Information
   logical :: print_avec    !print avec/bvec values?
   logical :: print_driver  !print driver data?
   logical :: print_fire    !print fire data?
   logical :: print_harvest !print harvest information?
   logical :: print_pftinfo !print PFT information?
   logical :: print_pooll   !print live pool values?
   logical :: print_soil    !print soil properties?
   logical :: print_sscol   !print soil/snow layer info?
   logical :: print_veg     !print vegetation values?
   logical :: print_stop    !stop when printed?

end module module_sibconst
