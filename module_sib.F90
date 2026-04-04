Module module_sib

!----------------------------------------------------------------------
!
!   SiB4 Variable Module
!
!   Using derived type hierarchy.  
!   gridcell -> landpft
!
!   ---
!   every point to simulate consists of a 
!      single gridcell type
!   ---
!   landpft types currently have corresponding
!      land, canopy air space (CAS),
!      canopy, and pft variables
!   ---
!
!   pft types can have values 
!   1=> not vegetated (bare ground)
!   2=> needleleaf evergreen
!   4=> needleleaf deciduous
!   5=> broadleaf evergreen
!   8=> broadleaf deciduous
!  11=> shrub (non-tundra)
!  12=> tundra shrub
!  13=> tundra c3 grass
!  14=> c3 grass (non-tundra)
!  15=> c4 grass
!  17=> generic crop
!  20=> maize
!  22=> soybean
!  24=> winter wheat
!----------------------------------------------------------------------

use kinds
implicit none

!=================================================
!Time Invariant Variables
!---------------------------------------------------------------------
! Soil Properties (Soil)
!---------------------------------------------------------------------
type, public :: soil_type

    !...soil properties
    real(r8) :: sandfrac  ! soil sand fraction
    real(r8) :: clayfrac  ! soil clay fraction
    real(r8) :: soref_vis ! soil shortwave reflectance (-)
    real(r8) :: soref_nir ! soil longwave reflectance (-)

    real(r8) :: poros  ! soil porosity (zero to one)
    real(r8) :: satco  ! hydraulic conductivity at saturation (m/s)

    real(r8) :: &
          csolid, &  ! heat capacity, soil solids (J/m3/K)
          tkdry,  &  ! thermal conductivity, dry soil (W/m/K)      
          tkmg,   &  ! thermal conductivity, soil minerals (W/m/K)
          tksat      ! thermal conductivity, saturated soil (W/m/K)

    !...soil variables for plant water stress
    real(r8) :: bee      ! Clapp & Hornberger 'b' exponent (-)    
    real(r8) :: phsat    ! Soil tension at saturation (m)
    real(r8) :: fieldcap ! Field Capacity (1/m3)
    real(r8) :: vwcmin   ! Wilting Point (1/m3)

    !...soil variables for respiration
    real(r8) :: wopt   ! respiration function optimum soil moisture (-)
    real(r8) :: woptzm ! wopt to the zm exponent
    real(r8) :: wsat   ! respiration function value (-)
    real(r8) :: zm     ! respiration function exponent (-)

     !...soil variables for water availability
     real(r8), dimension(:), allocatable :: &  !(nsoil)
          fc_eff     !Effective Field Capacity (1/m3)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
          wp_eff     !Effective Wilting Point (1/m3)

end type soil_type

!===================================================================
!=================================================
!Time-Step Varying Variables
!-------------------------------------------------
! Canopy (Can) and Canopy Air Space (CAS) Variables 
!-------------------------------------------------
type, public :: cas_type

     !...Canopy prognostic variables
     real(r8) :: tc    !canopy temperature (K)
 
     !...CAS prognostic variables
     real(r8) :: eacas    ! CAS water vapor pressure (hPa or mb)
     real(r8) :: shcas    ! CAS water vapor mixing ratio (kg/kg)
     real(r8) :: tcas     ! CAS temperature (K)
     real(r8) :: tkecas   ! CAS turbulent kinetic energy (J/kg)

     !...Canopy environmental variables
     real(r8) :: hcapc    !canopy heat capacity (J/m2/K)
     real(r8) :: tcmin    !frost 

     !...CAS environmental variables
     real(r8) :: thcas    ! CAS potential temperature (K)
     real(r8) :: hcapcas  ! CAS heat capacity (J/m2/K)
     real(r8) :: vcapcas  ! CAS vapor capacity (J/m2/hPa)
     real(r8) :: vpd      ! CAS vapor pressure deficit (hPa or mb)

     !...Additional pressure diagnostics
     real(r8) :: press_addinc
     real(r8) :: press_flxvrbrd
     real(r8) :: press_flxupdate
     real(r8) :: press_flxupdateps

end type cas_type

!-------------------------------------------------
! CO2/Photosynthesis Variables (CAN; CAS; PFT)
!-------------------------------------------------
type, public :: co2_type

     !...assimilation
     real(r8) :: assim  ! gross assimilation (mol C/m2/s)
     real(r8) :: assimd ! daily (24hr running-mean) assimilation (mol C/m2/s)
     real(r8) :: clim_assim ! climatological assimilation (mol C/m2/s)

     !...canopy scaling
     real(r8) :: assimpot  !potential top leaf photosynthesis (mol C/m2/s)
     real(r8) :: apar    ! absorbed photosynthetically active radiation (mol/m2/s)
     real(r8) :: aparkk  ! factor for scaling of leaf radiation (-)
     real(r8) :: gamma   ! CO2 photocompensation point (Pa)
     real(r8) :: par     ! photosynthetically active radiation (mol/m2/s)
     real(r8) :: nspar   ! non-scattered photosynthetically active radiation (mol/m2/s)

     
!...canopy air space (CAS) fluxes
     real(r8) :: casd    ! CAS depth for CO2 (m)
     real(r8) :: cflux   ! CAS to ref height carbon flux (mol C/m2/s)

     !...co2 concentrations
     real(r8) :: pco2cas ! CAS CO2 partial pressure (Pa)
     real(r8) :: pco2c   ! chloroplast CO2 partial pressure (Pa) 
     real(r8) :: pco2i   ! leaf internal CO2 partial pressure (Pa) 
     real(r8) :: pco2s   ! leaf surface CO2 partial pressure (Pa)
     real(r8) :: pco2m   ! mixed layer CO2 partial pressure (Pa)

     real(r8) :: co2cas ! CAS CO2 concentration (mol C/mol air)
     real(r8) :: co2m    ! reference level CO2 concentration (mol C/mol air)
     real(r8) :: co2s    ! leaf surface CO2 concentration (mol C/mol air)
     real(r8) :: co2i    ! leaf internal CO2 concentration (mol C/mol air)
     real(r8) :: co2c    ! leaf chloroplast CO2 concentration (mol C/mol air)
     real(r8) :: co2gamma    ! CO2 compensation point (mol C/mol air)

     !...resistances
     real(r8) :: rst ! prognostic stomatal resistance (s/m)

     !..soil freeze functions
     real(r8) :: soilfrz   !soil freeze function (-)
     real(r8) :: soilfrztg !soil freeze function for top soil layer (-)
     real(r8) :: soilfrztd !soil freeze function for second soil layer (-)

     !.....total stress factors
     real(r8) :: rstfac(4) ! canopy stress factors (-)
                           !  (1) leaf surface RH stress
                           !  (2) rootzone water stress
                           !  (3) temperature stress
                           !  (4) product of factors 1-3
     real(r8) :: vmaxss    ! stressed rubisco velocity (mol/m2/s)

     !itb...some diagnostics to determine accumulated stress. 
     real(r8) :: gs_stress(7) !
     ! (1) high-temperature stress (-) 
     ! (2) low-temperature stress (-) 
     ! (3) frost stress (-)
     ! (4) rstfac2 (-)
     ! (5) vpd (hPa or mb)
     ! (6) gpp (umol m-2 sec-1)
     ! (7) ect (transpiration: W m-2)

     !... Additional pressure diagnostics
     !real(r8) :: press_phosib
     !real(r8) :: press_phosibps
     integer(i4) :: icconv_phosib 
     real(r8) :: eyy_phosib
     real(r8) :: resp_casn    !canopy resp factoring vmax0, rstfac2, qt
     !integer(i4) :: icconv_phosib     

end type co2_type

!-------------------------------------------------
! Carbon Isotopes Variables
!-------------------------------------------------
type, public :: fract_type

    !...carbon isotopes
    real(r8) :: d13cca   ! isotope effects
    real(r8) :: d13cca_updated ! updated for varco2,varciso runs in c13_iso_calc
    real(r8) :: d13cm     ! 
    real(r8) :: c13ca     ! updated in cfrax
    real(r8) :: c12ca     ! updated in cfrax
    real(r8) :: c13cca     ! updated in c13_iso_calc based on recycling
    real(r8) :: c12cca     ! updated in c13_iso_calc based on recycling
    real(r8) :: c13cm     ! updated in cfrax
    real(r8) :: c12cm     ! updated in cfrax
    real(r8) :: kiecps    ! fractionation
    real(r8) :: kiecps_nog    ! fractionation, no gamma term
    real(r8) :: rcassim    ! isotope ratio value of assimilation 
    real(r8) :: rcassimfac
    real(r8) :: rcassim_nog    ! isotope ratio value of assimilation, no gamma term 
    real(r8) :: d13cassim   ! delta value of assimilation
    real(r8) :: d13cassim_nog   ! delta value of assimilation, no gamma term
    real(r8) :: c13assim    ! recently assimilated carbon-13 (mol C/m2/s)
    real(r8) :: c13assim_nog    ! recently assimilated carbon-13 (mol C/m2/s), no gamma
    real(r8) :: c12assim    ! recently assimilated carbon-12 (mol C/m2/s)
    real(r8) :: c13assimd  ! daily (24hr running-mean) C-13 assimilation (mol C/m2/s)
    real(r8) :: c13resptot    ! recently respired carbon-13
    real(r8) :: c12resptot    ! recently respired carbon-12

    real(r8) :: c13assimn   ! recently net assimilated carbon-13 (mol C/m2/s)
    real(r8) :: c12assimn   ! recently net assimilated carbon-12 (mol C/m2/s)

    real(r8) :: kiecps_k1   ! fractionation term 1
    real(r8) :: kiecps_k2   ! fractionation term 2
    real(r8) :: kiecps_k3   ! fractionation term 3
    real(r8) :: kiecps_k4   ! fractionation term 4
    real(r8) :: kiecps_k5   ! fractionation term 5

    real(r8) :: d14cca   ! isotope effects
    real(r8) :: d14cm     ! 
    real(r8) :: c14ca     !
    real(r8) :: c14cm     !
    real(r8) :: c14cca
    
    real(r8) :: rcassimc14    ! isotope ratio value of assimilation 
    real(r8) :: rcassimfacc14
    real(r8) :: d14cassim   ! delta value of assimilation
    real(r8) :: D_14cassim
    real(r8) :: c14assim    ! recently assimilated carbon-13 (mol C/m2/s)
    real(r8) :: c14assimn   ! recently net assimilated carbon-14 (mol C/m2/s)
    real(r8) :: c14assimd   ! daily (24-hour running mean) assimilated carbon-14 (mol C/m2/s)
    real(r8) :: c14resptot 

    real(r8) :: c14alpha  ! 14alpha_ph = (13alpha_ph)^2, alpha_ab = 1 + wtkiecps/1000.
    real(r8) :: normfac

    !...pool variables
    real(r8) :: d13cpool_leafc13
    real(r8) :: d13cpool_frootc13
    real(r8) :: d13cpool_crootc13
    real(r8) :: d13cpool_stwdc13
    real(r8) :: d13cpool_prodc13
    real(r8) :: d13cpool_cdbc13
    real(r8) :: d13cpool_metlc13
    real(r8) :: d13cpool_strlc13
    real(r8) :: d13cpool_slitc13
    real(r8) :: d13cpool_slowc13
    real(r8) :: d13cpool_armc13

    real(r8) :: d14cpool_leafc14
    real(r8) :: d14cpool_frootc14
    real(r8) :: d14cpool_crootc14
    real(r8) :: d14cpool_stwdc14
    real(r8) :: d14cpool_prodc14
    real(r8) :: d14cpool_cdbc14
    real(r8) :: d14cpool_metlc14
    real(r8) :: d14cpool_strlc14
    real(r8) :: d14cpool_slitc14
    real(r8) :: d14cpool_slowc14
    real(r8) :: d14cpool_armc14

    !...resp variables
    real(r8) :: d13cresp_autoc13
    real(r8) :: d13cresp_growc13
    real(r8) :: d13cresp_leafc13
    real(r8) :: d13cresp_mntnc13
    real(r8) :: d13cresp_rootc13
    real(r8) :: d13cresp_hetc13
    real(r8) :: d13cresp_firec13
    real(r8) :: d13cemisfire_pool
    real(r8) :: d13cresp_grzc13
    real(r8) :: d13cresp_hrvstc13
    real(r8) :: d13cresp_nvegc13
    real(r8) :: d13cresp_soilc13
    real(r8) :: d13cresp_totc13

    real(r8) :: d14cresp_autoc14
    real(r8) :: d14cresp_growc14
    real(r8) :: d14cresp_leafc14
    real(r8) :: d14cresp_mntnc14
    real(r8) :: d14cresp_rootc14
    real(r8) :: d14cresp_hetc14
    real(r8) :: d14cresp_firec14
    real(r8) :: d14cemisfire_pool
    real(r8) :: d14cresp_grzc14
    real(r8) :: d14cresp_hrvstc14
    real(r8) :: d14cresp_nvegc14
    real(r8) :: d14cresp_soilc14
    real(r8) :: d14cresp_totc14

    !...met variables
    real(r8) :: co2m_cfrax
    real(r8) :: press_cfrax
    real(r8) :: press_cfraxps

    !...fire variables
    real(r8) :: rcpoolfirec13 !pool ratio for (lp+wp+cdb+metl+strl) for fire
    real(r8) :: rcpoolfirec14
    real(r8) :: poolemistotC !pool sum for totC for (lp+wp+cdb+metl+strl) 
    real(r8) :: poolemisc13 !pool sum for C13 for (lp+wp+cdb+metl+strl) 
    real(r8) :: poolemisc14

    !real(r8) :: d13cassimxassim
    !real(r8) :: kiecpsxassim
    !real(r8) :: kiecpsk1xassim
    !real(r8) :: kiecpsk2xassim
    !real(r8) :: kiecpsk3xassim
    !real(r8) :: kiecpsk4xassim
    !real(r8) :: kiecpsk5xassim

    !real(r8) :: d13crautoxrauto
    !real(r8) :: d13crgrowxrgrow
    !real(r8) :: d13crleafxrleaf
    !real(r8) :: d13crmntnxrmntn
    !real(r8) :: d13crrootxrroot
    !real(r8) :: d13crfirexrfire
    !real(r8) :: d13crgrzxrgrz
    !real(r8) :: d13crhrvstxrhrvst
    !real(r8) :: d13crnvegxrnveg
    !real(r8) :: d13crhetxrhet
    !real(r8) :: d13crsoilxrsoil
    !real(r8) :: d13crtotxrtot

    !real(r8) :: d13cpleafxpleaf
    !real(r8) :: d13cpfrootxpfroot
    !real(r8) :: d13cpcrootxpcroot
    !real(r8) :: d13cpstwdxpstwd
    !real(r8) :: d13cpprodxpprod
    !real(r8) :: d13cpcdbxpcdb
    !real(r8) :: d13cpmetlxpmetl
    !real(r8) :: d13cpstrlxpstrl
    !real(r8) :: d13cpslitxpslit
    !real(r8) :: d13cpslowxpslow
    !real(r8) :: d13cparmxparm

end type fract_type


!-------------------------------------------------
! Carbonyl Sulfide (COS) Variables (CAN; CAS; PFT)
!-------------------------------------------------
type, public :: cos_type

    !...COS CAS
     real(r8) :: cos_casd ! CAS depth for COS (m)
     real(r8) :: cos_flux ! CAS COS flux (mol/m2/s)

     !...COS concentrations
     real(r8) :: coscas   ! CAS COS (mol COS/mol air)
     real(r8) :: coss     ! Leaf Surface COS (mol COS/mol air)
     real(r8) :: cosi     ! Leaf Internal COS (mol COS/mol air)
     real(r8) :: coscasp  ! CAS COS partial pressure (Pa)

     !...COS fluxes
     real(r8) :: cos_assim  ! COS assimilation (mol/m2/s)
     real(r8) :: cos_lru    ! COS leaf relative uptake (-)
     real(r8) :: cos_lru2   ! COS leaf relative uptake, ci/ca calculation (-)
     real(r8) :: cos_lru3   ! COS leaf relative uptake variant 3
     real(r8) :: cos_lru4   ! COS leaf relative uptake variant 4
     real(r8) :: cos_grnd   ! COS uptake by the soil (mol/m2/s)
     real(r8) :: cos_grnd_Ogee   ! COS uptake by the soil (mol/m2/s)
     real(r8) :: cos_soil ! COS uptake by the soil (mol/m2/s)

     !...COS conductances
     real(r8) :: gsh2onew   ! canopy conductance (mol/m2/s)
     real(r8) :: cosm       ! reference level COS concentration (mol COS/mol air)

     !...COS concentration in the soil
     real(r8), dimension(:), allocatable :: & !(nsoil)
          cos_s             ! COS concentration in the soil (xxxx)

     !...Conductance diagnostics
     real(r8) :: cosgm   ! 'apparent' mesophyll conductance (mol/m2/s)
     real(r8) :: cosgt   ! total conductance (mol/m2/s)
     real(r8) :: badcos  ! Diagnostic for bad COS calculations

     !...Pressure diagnostic
     real(r8) :: press_coscalc
     
end type cos_type


!-------------------------------------------------
! Dead Pool Equilibrium Variables (Soil)
!-------------------------------------------------
type, public :: equibd_type

      !...prognostic land pool variables for equilibrium
      real(r8), dimension(:), allocatable :: &  !(npoollu)
           poollu_totgain, & !sum of pool gains (mol/m2)
           poollu_totloss    !sum of pool losses (mol/m2)

     !...equilibrium variables for individual dead pools
     real(r8), dimension(:), allocatable :: &  !(npoollu)
            poollu_init,  & !initial pools (mol C/m2)
            poollu_end,   & !ending pools (mol C/m2)
            poollu_min,   & !minimum pool value (mol C/m2)
            poollu_max,   & !maximum pool value (mol C/m2)
            poollu_gain,  & !net pool gain (mol C/m2)
            poollu_loss,  & !net pool loss (mol C/m2)
            poollu_ratio, & !ratio of input/output (-)
            poollu_equib    !equilibrium pools (mol C/m2)
     
     logical, dimension(:), allocatable :: &  !(npoollu)
            poollu_notdone !flag for if dead pools are spunup

     !...equilibrium variables for surface pools
     !.....cwd + litmet + litstr
     real(r8) :: deadsfc_init  !initial pool total (mol C/m2)
     real(r8) :: deadsfc_end   !ending pool total (mol C/m2)
     real(r8) :: deadsfc_gain  !net gain (mol C/m2)
     real(r8) :: deadsfc_loss  !net loss (mol C/m2)
     real(r8) :: deadsfc_ratio !gain/loss ratio (-)
     logical  :: deadsfc_notdone  !flag for if pools are spunup

     !...equilibrium variables for surface C-13 pools
     !.....cwd + litmet + litstr
     real(r8) :: deadsfc_c13_init  !initial pool total (mol C/m2)
     real(r8) :: deadsfc_c13_end   !ending pool total (mol C/m2)
     real(r8) :: deadsfc_c13_gain  !net gain (mol C/m2)
     real(r8) :: deadsfc_c13_loss  !net loss (mol C/m2)
     real(r8) :: deadsfc_c13_ratio !gain/loss ratio (-)
     logical  :: deadsfc_c13_notdone  !flag for if pools are spunup

     !...equilibrium variables for surface C-14 pools
     !.....cwd + litmet + litstr
     real(r8) :: deadsfc_c14_init  !initial pool total (mol C/m2)
     real(r8) :: deadsfc_c14_end   !ending pool total (mol C/m2)
     real(r8) :: deadsfc_c14_gain  !net gain (mol C/m2)
     real(r8) :: deadsfc_c14_loss  !net loss (mol C/m2)
     real(r8) :: deadsfc_c14_ratio !gain/loss ratio (-)
     logical  :: deadsfc_c14_notdone  !flag for if pools are spunup

    !...equilibrium variables for soil pools
    !.....slit + slow + arm
     real(r8) :: deadsoil_init  !initial pool total (mol C/m2)
     real(r8) :: deadsoil_end   !ending pool total (mol C/m2)
     real(r8) :: deadsoil_gain  !net gain (mol C/m2)
     real(r8) :: deadsoil_loss  !net loss (mol C/m2)
     real(r8) :: deadsoil_ratio !gain/loss ratio (-)
     logical  :: deadsoil_notdone  !flag for if pools are spunup

    !...equilibrium variables for soil C-13 pools
    !.....slit + slow + arm
     real(r8) :: deadsoil_c13_init  !initial pool total (mol C/m2)
     real(r8) :: deadsoil_c13_end   !ending pool total (mol C/m2)
     real(r8) :: deadsoil_c13_gain  !net gain (mol C/m2)
     real(r8) :: deadsoil_c13_loss  !net loss (mol C/m2)
     real(r8) :: deadsoil_c13_ratio !gain/loss ratio (-)
     logical  :: deadsoil_c13_notdone  !flag for if pools are spunup

    !...equilibrium variables for soil C-14 pools
    !.....slit + slow + arm
     real(r8) :: deadsoil_c14_init  !initial pool total (mol C/m2)
     real(r8) :: deadsoil_c14_end   !ending pool total (mol C/m2)
     real(r8) :: deadsoil_c14_gain  !net gain (mol C/m2)
     real(r8) :: deadsoil_c14_loss  !net loss (mol C/m2)
     real(r8) :: deadsoil_c14_ratio !gain/loss ratio (-)
     logical  :: deadsoil_c14_notdone  !flag for if pools are spunup

    !...spin-up variables
    logical :: lupft_spunup

end type equibd_type


!-------------------------------------------------
! Live Pool Equilibrium Variables (PFT)
!-------------------------------------------------
type, public :: equibl_type

    !...prognostic pool variables for equilibrium calculation
    real(r8), dimension(:), allocatable :: & !(npoolpft)
          poolpft_totgain, & !sum of pool gains (mol/m2)
          poolpft_totloss    !sum of pool losses (mol/m2)

     !...equilibrium variables for individual pools
     real(r8), dimension(:), allocatable :: & !(npoolpft)
          poolpft_init,  & !initial pools (mol C/m2)
          poolpft_end,   & !ending pools (mol C/m2)
          poolpft_min,   & !minimum pool value (mol C/m2)
          poolpft_max,   & !maximum pool value (mol C/m2)
          poolpft_gain,  & !net gain (mol C/m2)
          poolpft_loss,  & !net loss (mol C/m2)
          poolpft_ratio, & !ratio of input/output (-)
          poolpft_equib    !equilibrium pools (mol C/m2)
     
     logical, dimension(:), allocatable :: & !(npoolpft)
          poolpft_notdone  !flag for if live pools are spunup

     !...equilibrium variables for live pool sums 
     !.....leaf + root + wood + prod
     real(r8) :: live_init     !initial live pool total (mol C/m2)
     real(r8) :: live_end      !ending live pool total (mol C/m2)
     real(r8) :: live_gain     !live pool net gain (mol C/m2)
     real(r8) :: live_loss     !live pool net loss (mol C/m2)
     real(r8) :: live_ratio    !gain/loss ratio (-)
     logical  :: live_notdone     !flag for if live pools are spunup


     !...equilibrium variables for live pool sums for C-13 pools
     !.....leaf + root + wood + prod
     real(r8) :: live_c13_init     !initial live pool total (mol C/m2)
     real(r8) :: live_c13_end      !ending live pool total (mol C/m2)
     real(r8) :: live_c13_gain     !live pool net gain (mol C/m2)
     real(r8) :: live_c13_loss     !live pool net loss (mol C/m2)
     real(r8) :: live_c13_ratio    !gain/loss ratio (-)
     logical  :: live_c13_notdone     !flag for if live pools are spunup

     !...equilibrium variables for live pool sums for C-14 pools
     !.....leaf + root + wood + prod
     real(r8) :: live_c14_init     !initial live pool total (mol C/m2)
     real(r8) :: live_c14_end      !ending live pool total (mol C/m2)
     real(r8) :: live_c14_gain     !live pool net gain (mol C/m2)
     real(r8) :: live_c14_loss     !live pool net loss (mol C/m2)
     real(r8) :: live_c14_ratio    !gain/loss ratio (-)
     logical  :: live_c14_notdone     !flag for if live pools are spunup

end type equibl_type


!---------------------------------------------------------------------
! Flux Variables (CAN; CAS; PFT)
!---------------------------------------------------------------------
type, public :: flux_type

    !...land-atmosphere exchange info
    real(r8) :: ct       !thermal transfer coefficient (-)
    real(r8) :: cu       !momentum transfer coefficient (-) 
    real(r8) :: drag     !drag (kg/m2/s)
    real(r8) :: ustar    !friction velocity (m/s)
    real(r8) :: ventmf   !ventilation mass flux (kg/m2/s)

     !...latent heat flux
     real(r8) :: ec      ! canopy latent heat flux (J/m2)
     real(r8) :: eci     ! latent heat flux, canopy interception (puddles) (J/m2) 
     real(r8) :: ect     ! latent heat flux, canopy transpiration (J/m2)
     real(r8) :: eg      ! ground latent heat flux (J/m2)
     real(r8) :: egi     ! latent heat flux, ground interception (puddles) (J/m2) 
     real(r8) :: egs     ! latent heat flux, ground evaporation (J/m2)
     real(r8) :: egsmax  ! maximum ground evapotration per timestep (J/m2)
     real(r8) :: es      ! snow latent heat flux (J/m2)
     real(r8) :: fws     ! CAS-BL latent heat flux (W/m2)

     !...sensible heat flux
     real(r8) :: hc      ! canopy sensible heat flux (J/m2)
     real(r8) :: hg      ! ground sensible heat flux (J/m2)
     real(r8) :: hs      ! snow sensible heat flux (J/m2)
     real(r8) :: fss     ! CAS-BL sensible heat flux (W/m2)

     !...storage heat flux
     real(r8) :: storhc  ! canopy heat storage flux (W/m2)
     real(r8) :: storhg  ! ground heat storage flux (W/m2)

     !...resistances
     real(r8) :: ra      ! canopy air space - mixed layer resistance (s/m)
     real(r8) :: rb      ! canopy to canopy air space resistance (s/m)
     real(r8) :: rbc     ! canopy to canopy air space resistance 
                         !    adjusted for snow (-)
     real(r8) :: rc      ! bulk leaf to canopy resistance (s/m)
     real(r8) :: rd      ! ground to canopy air space resistance (s/m)
     real(r8) :: rdc     ! ground to canopy air space resistance
                         !    adjusted for snow (-)
     real(r8) :: rds     ! ground and soil resistance (s/m)

    !...balance checks
    integer(i4) :: ebalnum  !energy balance 
    integer(i4) :: wbalnum  !water balance

    !..potential ET terms
    real(r8) :: et0      !ASCE Penman-Monteith calculation
    real(r8) :: et1      !ASCE Standardized calculation
    real(r8) :: et0a     !Accumulated ASCE Penman-Monteith calculation
    real(r8) :: et1a     !Accumulated ASCE Standardized calculation

end type flux_type


!---------------------------------------------------------------------
! Soil Hydrology Variables (Soil)
!---------------------------------------------------------------------
type, public :: hydros_type

     !...environmental variables
     real(r8) :: rhsoil  ! soil surface relative humidity (-)
     real(r8) :: rsoil   ! soil surface resistance 
                         !  (due to surface tension, s/m)

     !...evapotranspiration
     real(r8) :: ecmass   ! canopy evapotranspiration (kg/m2 or mm water)
     real(r8) :: egmass   ! ground evaporation (kg/m2 or mm water)

     !...precipitation
     real(r8) :: infil     ! water infiltrated into top soil layer (mm)
     real(r8) :: p0        ! ground surface precip (mm)
     real(r8) :: pcpg_rain ! ground surface rain precip (mm/s)
     real(r8) :: pcpg_snow ! ground surface snow precip (mm/s)

     !...runoff
     real(r8) :: roff      ! total subsurface runoff from soil layers (mm)
     real(r8) :: roffo     ! overland runoff (mm)

     !...snow
     real(r8) :: snow_gdepth  ! depth of snow on ground (m)
     real(r8) :: snow_gmass   ! mass of snow on ground (kg/m2)
     real(r8) :: snow_gvfc    ! snow ground cover fraction (0-1)
                              !   (formerly areas)

     !...soil diagnostics
     real(r8) :: www_tot     !total soil water-all layers, water+ice (kg/m2)
     real(r8) :: www_inflow  !water inflow at ground surface (kg/m2/s)
     real(r8) :: satfrac     !total fraction of water saturation in soil column (-)
     
     !...water interception
     real(r8) :: capacc_liq   ! prognostic canopy surface liquid (kg/m2)
     real(r8) :: capacc_snow  ! prognostic canopy surface snow (kg/m2)
                              !   (formerly snow_veg)
     real(r8) :: capacg       ! prognostic ground surface liquid (kg/m2)
     real(r8) :: satcapc      ! canopy wetness storage limit (kg/m2)
     real(r8) :: satcapg      ! ground wetness storage limit (kg/m2)

     real(r8) :: snow_cvfc  ! snow vertical cover fraction (-)
                            !   (formerly canex=1-snow_cvfc)
     real(r8) :: wetfracc  ! canopy wetness fraction (-)
     real(r8) :: wetfracg  ! ground wetness fraction (-)

end type hydros_type


!-------------------------------------------------
! Vegetation-Specific Hydrology Variables (PFT)
!-------------------------------------------------
type, public :: hydrov_type
     !...rooting zone information
     !......Plant Available Water (PAW; liquid only)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
         paw_lay,    &  !PAW per soil layer (kg/m3)
         pawmax_lay, &  !PAW maximum per soil layer (kg/m3)
         pawfrac_lay    !PAW fraction per soil layer (-)
     real(r8) :: &
         pawfrw,  &     !Root-weighted PAW fraction in soil column (kg/m2)
         pawftop, &     !Mean PAW fraction in top 3 soil layers (-)
         pawfzw         !Soil-layer depth-weighted PAW fraction (kg/m2)

     !......Total Available Water (TAW; liquid + ice)
     real(r8), dimension(:), allocatable :: &
         taw_lay,     & !TAW per soil layer (kg/m3)
         tawfrac_lay    !TAW fraction per soil layer (-)

     real(r8) :: &
          tawfrw,  &  !Root-weighted TAW in soil column (kg/m2)
          tawftop, &  !Mean TAW fraction in top 3 soil layers (-)
          tawfzw      !Soil-layer depth weighted TAW fraction (kg/m2)

     !.....Climatological Water Availability
     real(r8) :: clim_pawfrw !climatological root-weighted PAW fraction (-)
     real(r8) :: clim_tawfrw !climatological root-weighted TAW fraction (-)

end type hydrov_type


!-----------------------------------------------------
! Phenology Variables (PFT)
!-----------------------------------------------------
type, public :: phen_type

     !...growing season start determinants
     real(r8) :: phenave_assim     !Running-Mean Assimilation (mol C/m2/s)
     real(r8) :: phenave_assimsm   !Seasonal Maximum Mean Assimilation (mol C/m2/s)
     real(r8) :: phenave_assimpot  !Mean Assimilation Potential (-)
     logical ::  phenflag_assimlow !Assimilation Flag For Growing Season Reset

     real(r8) :: phenave_pr       !Seasonal Mean Precipitation (mm/day)
     real(r8) :: phenave_prsm     !Seasonal Maximum Mean Precipitation (mm/day)
     real(r8) :: phenave_prsdoy   !Seasonal Day of Maximum Precip (doy)
     real(r8) :: phenave_prcdoy   !Climatological Mean Day of Max Precip (doy)
     real(r8) :: phenave_prpot    !Seasonal Precipitation Potential (-)
     logical  :: phenflag_precip  !Precipitation Flag for Growing Season Start
     
     real(r8) :: phenave_tawftop   !Running-Mean TAW in Top 3 Soil Layers (-)
     logical  :: phenflag_moist    !Moisture Flag for Growing Season Start
     
     real(r8) :: phenave_tm        !Running-Mean Temperature (K)
     logical  :: phenflag_temp     !Temperature Flag for Growing Season Start
     
     logical :: phenflag_daylen    !Daylength Flag for Growing Season Start
     logical :: phenflag_gsspass   !Combined Growing Season Start Flag

     !...growing season information
     integer(i4) :: nd_dormant !# of days dormant
     integer(i4) :: nd_gs      !# of days of growing season
     integer(i4), dimension(:), allocatable :: & !(npstg-1)
             nd_stg   !# of days per stage (npstg-1)

     !...phenology stage
     integer(i4) :: phen_istage !Phenology Stage (1-5)
     real(r8) :: phen_pi       !Phenology Stage Index

     !...dynamic phenology stage variables
     real(r8) :: phens_dayl     !Phenology Stage Daylength Potential

     real(r8) :: phenc_climp    !Climatological Suitability (-)
     real(r8) :: phenc_laimax   !Max Potential LAI (m2/m2)
     real(r8) :: phenc_laimin   !Min Potential LAI (m2/m2)
     real(r8) :: phens_grw       !Phenology Stage Growth Potential

     real(r8) :: phenave_env  !Environmental Conditions Potential
     real(r8) :: phenave_wa   !Water Availability Potential
     real(r8) :: phenave_wac  !Combined Environmental and Water Potential
     real(r8) :: phenave_wacsm !Seasonal Maximum Combined Potential
     real(r8) :: phens_wx       !Phenology Stage Weather Potential

     !...defined phenology stage variables
     integer(i4) :: ipd  !planting date (doy)
     integer(i4) :: dapd    !days after planting date (days)
     integer(i4) :: dapdaf  !days after planting above freezing (days)
     real(r8) :: gdd     !growing degree days (-)
     real(r8) :: seed_pool    !seed pool carbon (mol C/m2)

end type phen_type


!---------------------------------------------------------------------
! Dead Pool Variables (Soil)
!---------------------------------------------------------------------
type, public :: poold_type

     !====Dead Pool Gains (per timestep)====!

     !-------------
     !Grazing Gains (per soil layer)
     real(r8), dimension(:,:), allocatable :: &  !(npoollu)
             gain_grz_lay !gain from grazing (mol C/m2/s)

     !--------------
     !Harvest Gains (per soil layer)
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
            gain_hrvst_lay !gain from harvest (mol C/m2/s)

     !-------------------------
     !Live Pool Transfer Gains (per soil layer)
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
            gain_transl_lay  !gain from live pools (mol C/m2/s)
   
     !-------------------------
     !Dead Pool Transfer Gains (per soil layer)
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
            gain_transd_lay  !gain from dead pools (mol C/m2/s)


     !====Dead Pool Losses (per timestep)====!
     !---------------------------------------

     !Fire Loss
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
           loss_fire_lay   !loss from fire (mol C/m2/s)

     !...loss from radioactive decay, zero for everything except C14
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
           loss_raddecay_lay  !loss from radioactive decay (mol C/m2/s)
 
     !Heterotrophic Respiration/Transfer Loss

     !...surface pools
     real(r8) :: mhrt_sfc_assim  !surface assimilation scalar (-)
     real(r8) :: mhrt_sfc_freeze !surface cold/freezing scalar (-)
     real(r8) :: mhrt_sfc_hot    !surface high temperature scalar (-)
     real(r8) :: mhrt_sfc_precip !surface precip scaling factor (-)
     real(r8) :: mhrt_sfc_scale  !surface respiration scaling coefficient (-)

     !...soil pools per soil layer
     real(r8), dimension(:), allocatable ::  &  !(nsoil)
            mhrt_soil_freeze_lay, & !freeze inhibition scalar (-)
            mhrt_soil_hot_lay,   & !high temperature scalar (-)
            mhrt_soil_moist_lay, & !soil moisture scalar (-)
            mhrt_soil_pawf_lay,  & !PAW fraction scalar (-)
            mhrt_soil_scale_lay    !combined soil scalar (-)

     !...soil scalars root-weighted
     real(r8) :: &
            mhrt_soil_assim,  & !assimilation scalar (-)
            mhrt_soil_freeze, & !freeze inhibition scalar (-)
            mhrt_soil_hot,    & !high temperature scalar (-)
            mhrt_soil_moist,  & !soil moisture scaling factor (-) 
            mhrt_soil_pawfrw, & !soil pawfrw scalar (-)
            mhrt_soil_precip, & !soil precipitation scalar (-)
            mhrt_soil_scale     !combined soil scalar (-)

     !...respiration and transfer information per soil layer
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
           kratert_lay,    & !scaled decay rate (1/s)
           loss_resp_lay,  & !loss from respiration (mol C/m2/s)
           loss_trans_lay    !loss from decay transfers (mol C/m2/s)

     !...combined respiration rates
     real(r8) :: resp_het  !heterotrophic respiration (mol C/m2/s)
     real(r8) :: resp_soil !soil respiration (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
            resp_soil_lay
     real(r8) :: resp_soilnr  !soil respiration w/o roots (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
            resp_soilnr_lay

     !... combined respiration rates for C-13 pools
     real(r8) :: resp_hetc13  !heterotrophic respiration (mol C/m2/s)
     real(r8) :: resp_soilc13 !soil respiration (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
            resp_soilc13_lay
     real(r8) :: resp_soilnrc13  !soil respiration w/o roots (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
            resp_soilnrc13_lay

     !... combined respiration rates for C-14 pools
     real(r8) :: resp_hetc14  !heterotrophic respiration (mol C/m2/s)
     real(r8) :: resp_soilc14 !soil respiration (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
            resp_soilc14_lay
     real(r8) :: resp_soilnrc14  !soil respiration w/o roots (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(nsoil)
            resp_soilnrc14_lay

     !------------------------
     !Daily net pool change (per soil layer)
     real(r8), dimension(:,:), allocatable ::   &  !(npoollu,nsoil)
           poollu_dgain,  & !pool gain (mol C/m2/day)
           poollu_dloss     !pool loss (mol C/m2/day)

     !------------------------
     !Prognostic Carbon Pools
     real(r8), dimension(:), allocatable ::   &  !(npoollu)
           poollu  !vertically integrated pool size (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
           poollu_lay !prognostic dead carbon pools (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoollu,nsoil)
           poollu_flay  !fraction of carbon per soil layer (-)

     !---------------
     !Carbon Balance
     real(r8), dimension(:), allocatable :: & !(npoollu)
         poollup !previous carbon pool (mol C/m2)

     !---------------
     !Isotope pool-related variables
     real(r8), dimension(:,:), allocatable ::  & !(npoollu,nsoil)
          rcpoollu_lay   !C13 or C14 pool ratio based on pool size, per layer
     real(r8), dimension(:), allocatable ::  & !(npoollu)
          rcpoollu       !C13 or C14 pool ratio based on pool size
     real(r8), dimension(:), allocatable :: &
          curpoollu      !current pool size computed in c13_iso_calc

end type poold_type


!-------------------------------------------------
! Live Pool Variables (PFT)
!-------------------------------------------------
type, public :: pooll_type

    !====Live Pool Gains (lpg, per timestep)====!
     !===========================================!

     !------------------
     !Assimilation Gains
     !...allocation fractions for live biomass (vary from 0 to 1)
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          alloc !allocation fractions

     !...dynamic allocation contributions (vary from -1 to 1)
     logical :: aadj_moist !allow moisture adjustments?
     logical :: aadj_temp  !allow temperature adjustments?

     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          alloc_phen,  & !allocation based on phenological state
          alloc_moist, & !allocation adjustments due to moisture stress 
          alloc_temp     !allocation adjustments due to temperature stress

     !...pool gains from photosynthesis divied up by allocation fraction
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          gain_assim  !gain from photosynthesis (mol C/m2/s)

     !---------------------
     !Seed (Transfer) Gains
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          gain_seed   !gain from seed (mol C/m2/s)

     !====Live Pool Losses (lpl, per timestep)====!
     !===========================================!

     !-----------------------
     !Autotrophic Respiration
     !....Growth.....
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          loss_gresp   !loss from growth resp (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          resp_nveg_tmp   !non-vegetation respiration (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          resp_nvegc13_tmp   !non-vegetation respiration (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
          resp_nvegc14_tmp   !non-vegetation respiration (mol C/m2/s)

     !.....Maintenance.....
     !...canopy respiration scaling coefficients
     real(r8) :: mcr_assim     !canopy assimilation scalar
     real(r8) :: mcr_freeze    !canopy freeze inhibition
     real(r8) :: mcr_hot       !canopy high temperature exponential
     real(r8) :: mcr_scale     !combined canopy respiration scalar

     !...root respiration scaling coefficients
     real(r8), dimension(:), allocatable :: &  !(per soil layer)
           mrr_freeze_lay, & !roots freeze inhibition scalar
           mrr_hot_lay,    & !roots high temperature exponential
           mrr_scale_lay     !combined roots scalar

     real(r8) :: &        !root-weighted
           mrr_assim,  &  !assimilation
           mrr_freeze, &  !freeze inhibition
           mrr_hot,    &  !high temperature
           mrr_lai,    &  !leaf support
           mrr_scale      !combined root-weighted scalar for respiration

     !...maintenance respiration information per soil layer
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          krater_lay, &  !scaled maintainance loss rate (1/s)
          loss_mresp_lay !loss from maintainence (mol C/m2/s)

     !.....Combined Respiration Rates.....
     real(r8) :: resp_auto !autotrophic respiration (mol C/m2/s)
     real(r8) :: resp_grow !growth respiration (mol C/m2/s)
     real(r8) :: resp_leaf !leaf respiration (mol C/m2/s)
     real(r8) :: resp_mntn !maintainence respiration (mol C/m2/s)
     real(r8) :: resp_nveg !non-vegetation respiration (mol C/m2/s)
     real(r8) :: resp_root !root respiration (mol C/m2/s)

     !.....Combined Respiration Rates for C-13 pools
     real(r8) :: resp_autoc13 !autotrophic respiration (mol C/m2/s)
     real(r8) :: resp_growc13 !growth respiration (mol C/m2/s)
     real(r8) :: resp_leafc13 !leaf respiration (mol C/m2/s)
     real(r8) :: resp_mntnc13 !maintainence respiration (mol C/m2/s)
     real(r8) :: resp_nvegc13 !non-vegetation respiration (mol C/m2/s)
     real(r8) :: resp_rootc13 !root respiration (mol C/m2/s)

     !.....Combined Respiration Rates for C-13 pools
     real(r8) :: resp_autoc14 !autotrophic respiration (mol C/m2/s)
     real(r8) :: resp_growc14 !growth respiration (mol C/m2/s)
     real(r8) :: resp_leafc14 !leaf respiration (mol C/m2/s)
     real(r8) :: resp_mntnc14 !maintainence respiration (mol C/m2/s)
     real(r8) :: resp_nvegc14 !non-vegetation respiration (mol C/m2/s)
     real(r8) :: resp_rootc14 !root respiration (mol C/m2/s)

     real(r8) :: resp_leafc12 !leaf respiration c12 (mol C/m2/s)
     !real(r8) :: resp_casn    !canopy resp factoring vmax0, rstfac2, qt

     !...loss from radioactive decay, zero for everything except C14
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
              loss_raddecay_lay  !loss from radioactive decay (mol C/m2/s)

     !----------------------
     !Live-To-Dead Transfer
     !.....transfers use fractions of the pool that is transferring
     !.....transfer = fraction * pool_size

     !...leaf transfer (fractions per day)
     real(r8) :: tfl_daylen   !shortening daylength
     real(r8) :: tfl_freeze   !freezing
     real(r8) :: tfl_dry      !water deficiency
     real(r8) :: tfl_pstage   !phenology stage
     real(r8) :: tfl_total    !total leaf transfer
     real(r8) :: tfl_totalc13    !total leaf transfer
     real(r8) :: tfl_totalc14    !total leaf transfer

     !.....turnover transfer
     real(r8), dimension(:), allocatable :: &  !(npoolpft)
           tf_turnover !fractions per day

     !...transfer loss information (per soil layer)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          loss_trans_lay  !loss to dead pool transfers (mol C/m2/s)

     !------------
     !Fire Loss
     real(r8) :: nd_fire   !# of days burned
     real(r8) :: resp_fire !fire respiration (mol C/m2/s)
     real(r8) :: rmmd_fire !fire emitted but not removed due to
                           !  fire dataset and SiB4 mismatch (mol C/m2/s)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
              loss_fire_lay  !loss from fire (mol C/m2/s)
     
     real(r8) :: resp_firec13 !fire respiration for C-13 (mol C/m2/s)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
              loss_firec13_lay  !loss from fire for C-13 (mol C/m2/s)
     real(r8) :: rmmd_firec13 !fire emitted but not removed due to
                           !  fire dataset and SiB4 mismatch (mol C/m2/s)

     real(r8) :: resp_firec14 !fire respiration for C-14 (mol C/m2/s)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
              loss_firec14_lay  !loss from fire for C-14 (mol C/m2/s)
     real(r8) :: rmmd_firec14 !fire emitted but not removed due to
                           !  fire dataset and SiB4 mismatch (mol C/m2/s)

     !------------
     !Grazing Loss
     real(r8) :: nd_grz !# of days grazed 
     real(r8) :: resp_grz !grazing respiration (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(npoolcan)
              loss_grz !loss from grazing (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(npoolcanc13)
              loss_grzc13 !loss from grazing (mol C/m2/s)
     real(r8) :: resp_grzc13 !C13 grazing respiration for C-13 (mol C/m2/s)
     real(r8), dimension(:), allocatable :: &  !(npoolcanc14)
              loss_grzc14 !loss from grazing (mol C/m2/s)
     real(r8) :: resp_grzc14 !C14 grazing respiration for C-13 (mol C/m2/s)

     !-------------
     !Harvest Loss
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
         loss_hrvst_lay     !loss from harvest per soil layer (mol C/m2/s)
     real(r8) :: resp_hrvst !C harvest respiration (mol C/m2/s)
     real(r8) :: resp_hrvstc13 !C13 harvest respiration for C-13 (mol C/m2/s)
     real(r8) :: resp_hrvstc14 !C14 harvest respiration for C-13 (mol C/m2/s)
     real(r8) :: rmvd_hrvst !C harvested and removed (mol C/m2)
     real(r8) :: rmvd_hrvstc13 !C13 harvested and removed (mol C/m2)
     real(r8) :: rmvd_hrvstc14 !C14 harvested and removed (mol C/m2)

     !-------------------------------------
     !Daily net pool change (per soil layer)
     real(r8), dimension(:,:), allocatable ::  & !(npoolpft,nsoil)
           poolpft_dgain, & !pool gain (mol C/m2/day)
           poolpft_dloss    !pool loss (mol C/m2/day)

     !------------------------
     !Prognostic Carbon Pools
     real(r8), dimension(:), allocatable ::  & !(npoolpft)
          poolpft      !vertically integrated pool size (mol C/m2)
     real(r8), dimension(:), allocatable ::  & !(npoolpft)
          poolpftmin   !min poolpft (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          poolpft_lay  ! prognostic carbon pools (mol C/m2)
     real(r8), dimension(:,:), allocatable :: & !(npoolpft,nsoil)
          poolpft_flay !fraction of carbon per soil layer (-)

     !---------------
     !Carbon Balance
     real(r8), dimension(:), allocatable :: & !(npoolpft)
         poolpftp   !previous carbon pool (mol C/m2)

     !--------------
     !.. Additional diagnostics for pool C13
     real(r8), dimension(:), allocatable :: & !(npoolpft)
          mrespr   !for check/calc maintentance resp in pool_auto_resp
     real(r8), dimension(:,:), allocatable :: & !(npoolpft,nsoil)
          pftavail_lay !for logging poolpft_avail_lay from pool_auto_tran
     real(r8) :: tfrac_lp
     real(r8) :: tfrac_lpc13
     real(r8) :: tfrac_lpc14
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          assim_dgain  ! gain in pool_assim (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          assim_dloss  ! loss in pool_assim (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          autoresp_dgain  ! gain in pool_auto_resp (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          autoresp_dloss  ! loss in pool_auto_resp (mol C/m2)
     real(r8), dimension(:,:), allocatable :: &  !(npoolpft,nsoil)
          autotran_dloss  ! loss in pool_auto_tran (mol C/m2)
     real(r8), dimension(:), allocatable ::  & !(npoolpft)
          poolpftmin_updated   !min poolpft updated (mol C/m2)

     real(r8), dimension(:), allocatable ::  & !(npoolpft)
          isofactorp   !min poolpft updated (mol C/m2)

     real(r8), dimension(:,:), allocatable ::  & !(npoolpft,nsoil)
          rcpoolpft_lay   !C13 or C14 pool ratio based on pool size, per layer
     real(r8), dimension(:), allocatable ::  & !(npoolpft)
          rcpoolpft       !C13 or C14 pool ratio based on pool size
     !real(r8) :: rcpoolfire !pool ratio for (lp+wp+cdb+metl+strl) for fire
     !real(r8) :: poolemistotC !pool sum for totC for (lp+wp+cdb+metl+strl) 
     !real(r8) :: poolemisc13 !pool sum for C13 for (lp+wp+cdb+metl+strl) 
     real(r8), dimension(:), allocatable :: &
          curpoolpft      !current pool size computed in c13_iso_calc

end type pooll_type


!---------------------------------------------------------------------
! Radiation Variables (CAN, CAS, Soil)
!---------------------------------------------------------------------
type, public :: rad_type

    !...albedos
    real(r8) :: albedo_visb   ! albedo, visible beam (-)
    real(r8) :: albedo_visd   ! albedo, visible diffuse (-)
    real(r8) :: albedo_nirb   ! albedo, nir beam (-)
    real(r8) :: albedo_nird   ! albedo, nir diffuse (-)

     !...radiation variables
    real(r8) :: radfacc(2,2)  ! canopy radiation absorption factors (-)
            !   (1,1) - visible, beam
            !   (1,2) - visible, diffuse
            !   (2,1) - nir, beam
            !   (2,2) - nir, diffuse
    real(r8) :: radfacg(2,2)  ! ground radiation absorption factors (-)
            !   (1,1) - visible, beam
            !   (1,2) - visible, diffuse
            !   (2,1) - nir, beam
            !   (2,2) - nir, diffuse

    real(r8) :: radc3c    ! absorbed radiation by canopy (W/m2)
    real(r8) :: radc3g    ! absorbed radiation by ground (W/m2)
    real(r8) :: radtc     ! canopy net radiation (W/m2)
    real(r8) :: radtg     ! ground net radiation (W/m2)
    real(r8) :: radts     ! snow net radiation (W/m2)
    real(r8) :: effgc     ! effective ground cover for thermal radiation (-)

    !...temperatures
    real(r8) :: tsfc      ! surface temperature (K)

end type rad_type


!-------------------------------------------------
! Fluorescence (SIF) Variables (CAN; CAS; PFT)
!-------------------------------------------------
type, public :: sif_type

     !...electron transports
     real(r8) :: sif_je    !electron transport
     real(r8) :: sif_jo    !max electron transport
     real(r8) :: sif_jejo  !fractional transport (je/jo)

     !...k coefficients, defined as the probability of
     !.....excitons to follow certain pathways (-) 
     real(r8) :: sif_kd    !Heat dissipation
     real(r8) :: sif_kn    !Non-photochemical quenching (NPQ)
     real(r8) :: sif_kp    !Photosynthesis

     !...x factor (0 when GPP=potential, 1 when GPP=0)
     real(r8) :: sif_x

     !...yields (-)
     real(r8) :: phi_d  !Heat dissipation yield
     real(r8) :: phi_f  !SIF yield
     real(r8) :: phi_n  !NPQ yield
     real(r8) :: phi_p  !Photosynthetic yield

     !...resulting sif values 
     real(r8) :: sif  !fluorescence (W m-2 sr-1 nm-1)

end type sif_type


!---------------------------------------------------------------------
! Soil/Snow Column Variables (Soil)
!---------------------------------------------------------------------
type, public :: sscol_type  

     !...prognostic number of snow layers (negative)
     integer(byte) :: nsl     

     !...soil diagnostics for soil column (nsoil)
     real(r8), dimension(:), allocatable :: &
         rootr,     & !effective rooting frac for soil hydrology (-)
         satfrac_lay  !fraction of water saturation (-)

     !...snow/soil diagnostic column variables 
     real(r8), dimension(:), allocatable :: & !(-nsnow+1:nsoil)
         eff_poros, & ! soil/snow liquid effective porosity (-)
         layer_z,   & ! soil/snow layer interface depth (m)     
         node_z,    & ! soil/snow layer node depth (m)
         shcap,     & ! soil/snow total heat capacity (J/m2/K)
         slamda,    & ! soil/snow heat flux term (W/m2/K)
         tksoil,    & ! soil/snow thermal conductivity (W/m/K)
         vol_liq,   & ! soil/snow liquid water volume (kg/m3)
         vol_ice      ! soil/snow ice volume (kg/m3)

     !...snow/soil prognostic column variables
     real(r8), dimension(:), allocatable :: & !(-nsnow+1:nsoil)
         dz, &      ! soil/snow layer thickness (m)
         td, &      ! soil/snow temperature (K)
         www_liq, & ! soil/snow liquid water (kg/m2)
         www_ice    ! soil/snow ice (kg/m2)

end type sscol_type


!-------------------------------------------------
! Vegetation Description and State Variables (PFT)
!-------------------------------------------------
type, public :: veg_type

     !...land-atmos interactions
     real(r8) :: z0       ! canopy snow-adjusted roughness length (m)
     real(r8) :: z0d      ! canopy roughness length (m)
     real(r8) :: zp_dispd ! zero-plane displacement (m)
     real(r8) :: zpd_adj  ! snow-adjusted zero-plane displacement (m)
     real(r8) :: zztemp   ! temperature height for mass flux (m)
     real(r8) :: zzwind   ! wind height for mass flux (m)

     !...resistances
     real(r8) :: cc1      ! bulk pbl resistance coefficient (sqrt(s/m))
     real(r8) :: cc2      ! ground to canopy air space resistance (-)     

     !...root profile
     real(r8), dimension(:), allocatable :: &  !(nsoil)
          rootf  ! root fraction (-)
 
     !...vegetation state 
     real(r8) :: fpar     ! absorbed fraction of PAR (-)
     real(r8) :: green    ! green fraction of LAI (-)
     real(r8) :: lai      ! leaf area index (-)
     real(r8) :: lait     ! canopy total leaf area index (w/ dead, -)
     real(r8) :: vcover   ! fraction of vegetation cover (-)

     !...vegetation properties
     real(r8) :: gmudmu   ! time-mean leaf projection (-)
     real(r8) :: park     ! solar absorption factor / 
                          !      extinction coefficient for PAR (-)
     real(r8) :: vmax     ! rubisco velocity (mol/m2/s)

     !...climatological vegetation states
     real(r8) :: clim_lai  ! climatological LAI (-)

end type veg_type


!=================================================
! Driver Data Variables
!------------------------------------------------------------------
! Grid Cell Diagnostic Variables
! -----------------------------------------------------------------
type, public :: gdiag_type

     !...daytime/sunlight properties
     real(r8) :: cosz     ! cosine of solar zenith angle (-)
     real(r8) :: daylen   ! length of daylight (hrs)
     real(r8) :: daylendt ! change in length of day (hrs)

     !...misc driver forcings
     real(r8) :: tmdf   ! daily (24-hr running mean) temperature (F)
     real(r8) :: thm    ! mixed layer potential temperature (K)
     real(r8) :: bps(2) ! (ps/1000)**kapa - turns theta into temp
     real(r8) :: em     ! mixed layer water vapor pressure (hPa or mb)
     real(r8) :: ros    ! surface air density (kg/m3)
     real(r8) :: psy    ! psycrometric constant (gamma) (hPa/K)

    !...precipitation properties
    real(r8) :: seas_precippot ! seasonal precipitation potential (-)

    !...radiation properties
    real(r8) :: radvbc    ! visible beam radiation (W/m2)
    real(r8) :: radvdc    ! visible diffuse radiation (W/m2)
    real(r8) :: radnbc    ! nir beam radiation (W/m2)
    real(r8) :: radndc    ! nir diffuse radiation (W/m2)

    real(r8) :: toa_solar   ! Top-of-Atmosphere insolation (W/m2)
    real(r8) :: toa_radvbc  ! TOA visible beam radiation (W/m2)
    real(r8) :: toa_radvdc  ! TOA visible diffuse radiation (W/m2)
    real(r8) :: toa_radnbc  ! TOA near-infrared (NIR) beam radiation (W/m2)
    real(r8) :: toa_radndc  ! TOA NIR diffuse radiation (W/m2)

    real(r8) :: toa_par   ! TOA PAR (mol/m2/s)
    real(r8) :: aod       ! aerosol + cloud optical depth (-)

    !...solar-induced fluorescence (SIF) properties
    real(r8) :: sif_atten ! SIF attenuation factor (-)
    logical, dimension(2) :: sif_flag  !conditions suitable for 
                                       !1=GOME2 2=OCO2?

    !...spinup variables
    logical :: gridcell_spunup

end type gdiag_type

!------------------------------------------------------------------
! Grid Cell Prognostic Variables
! -----------------------------------------------------------------
type, public :: gprog_type

    !...prognostic driver data/forcings
    real(r8) :: cupr      ! cumulus precipitation rate (mm/s)
    real(r8) :: cuprt     ! cumulus precip rate (m/s)
    real(r8) :: cupr1     ! cumulus precipitation rate (mm/s)
    real(r8) :: cupr2     ! cumulus precipitation rate (mm/s)
    real(r8) :: dlwbot    ! surface incident longwave 
    real(r8) :: dlwbot1   ! surface incident longwave 
    real(r8) :: dlwbot2   ! surface incident longwave 
    real(r8) :: lspr      ! stratiform precipitation rate (mm/s)
    real(r8) :: lsprt     ! stratiform precip rate (m/s) 
    real(r8) :: lspr1     ! stratiform precipitation rate (mm/s)  
    real(r8) :: lspr2     ! stratiform precipitation rate (mm/s) 
    real(r8) :: ps        ! surface pressure (hPa or mb)
    real(r8) :: ps1       ! surface pressure (hPa or mb)
    real(r8) :: ps2       ! surface pressure (hPa or mb)
    real(r8) :: sh        ! mixed layer water vapor mixing ratio (kg/kg)
    real(r8) :: sh1       ! mixed layer water vapor mixing ratio (kg/kg)
    real(r8) :: sh2       ! mixed layer water vapor mixing ratio (kg/kg)
    real(r8) :: spdm      ! wind speed (m/s) 
    real(r8) :: spdm1     ! wind speed (m/s) 
    real(r8) :: spdm2     ! wind speed (m/s) 
    real(r8) :: sw_dwn    ! surface incident shortwave radiation (W/m2)
    real(r8) :: sw_dwn1   ! surface incident shortwave radiation (W/m2)
    real(r8) :: sw_dwn2   ! surface incident shortwave radiation (W/m2)
    real(r8) :: tm        ! mixed layer temperature (K)
    real(r8) :: tm1       ! mixed layer temperature (K)
    real(r8) :: tm2       ! mixed layer temperature (K)

    !...fire emissions
    real(r8) :: firec     ! C loss from fire (mol C/m2/s)
    real(r8) :: firec1    ! C loss from fire (mol C/m2/s)
    real(r8) :: firec2    ! C loss from fire (mol C/m2/s)
    real(r8) :: fireco2   ! CO2 respiration from fire (mol C/m2/s)
    real(r8) :: fireco21  ! CO2 respiration from fire (mol C/m2/s)
    real(r8) :: fireco22  ! CO2 respiration from fire (mol C/m2/s)    
    
    !...carbon cycle
    real(r8) :: pco2m     ! Mixed Layer (background) CO2 partial pressure (Pa)
    real(r8) :: pcosm     ! Mixed Layer (background) COS partial pressure (Pa)

    real(r8) :: co2m      ! Mixed Layer (background) CO2 concentration (ppm)
    real(r8) :: cosm      ! Mixed Layer (background) COS concentration (ppt)

    !...daily values
    real(r8) :: tmd    ! daily (24-hr running mean) temperature (K)

    !...seasonal values
    real(r8) :: seas_precip !seasonal precipitation (mm/day)
    real(r8) :: seas_tm     !seasonal temperature (K)

    !...climatological values
    real(r8) :: clim_cupr   ! clim-mean convective precipitation (mm/day)
    real(r8) :: clim_precip ! clim-mean precipitation (mm/day)
    real(r8) :: clim_tm     ! clim-mean temperature (K)

    !...tm5 values ARACHO
    real(r8) :: cosm_tm5      ! COS mixing ratio (molCOS/molAir)
    real(r8) :: cosm_tm51     ! COS mixing ratio (molCOS/molAir)
    real(r8) :: cosm_tm52     ! COS mixing ratio (molCOS/molAir)
    !real(r8) :: pressure_tm5  ! pressure (Pa)
    !real(r8) :: pressure_tm51 ! pressure (Pa)
    !real(r8) :: pressure_tm52 ! pressure (Pa) 
     
end type gprog_type


!******************************************************************
!------------------------------------------------------------------
! Begin definition of spatial scaling hierarchy
!------------------------------------------------------------------

!------------
! define the land unit structure
! includes corresponding:
!    - soil column (Soil)
!    - vegetation (PFT)
!    - canopy (CAN) and canopy air space (CAS)
!------------
type, public :: lu_type

    ! vegetation information
    integer(i4) :: ipft   !pft reference number
    real(r4)    :: larea  !fraction of coverage per gridcell (0-1)

    ! time-invariant variables
    type(soil_type) :: soilt   !soil properties

    ! time-variant variables
    type(cas_type)     :: cast      !canopy air space variables (CAS) 
    type(co2_type)     :: co2t      !photosynthetic/CO2 variables
    type(fract_type)   :: fract     !fractionation variables
    type(cos_type)     :: cost      !carbonyl sulfide variables
    type(equibd_type)  :: equibdt   !dead pool equilibrium variables (Soil)
    type(equibl_type)  :: equiblt   !live pool equilibrium variables (PFT)
    type(flux_type)    :: fluxt     !flux variables
    type(hydros_type)  :: hydrost   !soil hydrological variables (Soil)
    type(hydrov_type)  :: hydrovt   !veg hydrological variables (PFT)
    type(phen_type)    :: phent     !phenology variables (PFT)
    type(poold_type)   :: pooldt    !dead pool variables (Soil)
    type(pooll_type)   :: poollt    !live pool variables (PFT)
    type(rad_type)     :: radt      !radiation variables
    type(sif_type)     :: sift      !fluorescence variables
    type(sscol_type)   :: sscolt    !soil/snow column variables (Soil)
    type(veg_type)     :: vegt      !vegetation information (PFT)

end type lu_type

!------------
! define the gridcell structure
!------------
type, public :: gridcell_type

    ! gridcell information
    real(r8) :: lat         !latitude (degrees)
    real(r8) :: lon         !longitude (degrees)

    ! variables defined at the gridcell level
    type(gdiag_type) :: gdiagt !diagnostic variables
    type(gprog_type) :: gprogt !prognostic variables

    ! gridcell -> landunit hierarchy
    integer(i4) :: g_nlu !number of land units per gridcell
                         !...includes corresponding
                         !....soil, veg (PFT), and canopy 

    type(lu_type), dimension(:), allocatable :: &
          l  !land unit data structure

end type gridcell_type

!-----------
! define the top-level structure
!-----------
type, public :: sib_t
   type(gridcell_type),  dimension(:), allocatable :: &
          g  !gridcell data structure
end type sib_t

!------------------------------------------------------------------
! End definition of spatial scaling hierarchy
!------------------------------------------------------------------
!******************************************************************

!***********************************************************************
!-----------------------------------------------------------------------
! Declare single instance of sibtype
    type(sib_t) :: sib
!-----------------------------------------------------------------------

end module module_sib
