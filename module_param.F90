Module module_param

!----------------------------------------------------------------------
!
!   SiB4 Parameter Module
!
!----------------------------------------------------------------------

use kinds
implicit none


!******************************************************************

!------------------------------------------------------------------
! Aerodynamic Parameters
!------------------------------------------------------------------
type, public :: aero_param

     real(r4) :: zo       !canopy roughness coeff
     real(r4) :: zp_disp  !zero plane displacement
     real(r4) :: RbC      !coefficient for canopy to CAS aerodynamic resistance
     real(r4) :: RdC      !coefficient for ground to CAS aerodynamic resistance

end type aero_param

!-------------------------------------------------
! Phenology Parameters
!-------------------------------------------------
type, public :: phen_param

    !-------------------------------------------------
    ! Defined and Dynamic Common Phenology Parameters
    !-------------------------------------------------

    !...growing season start factors
   integer(i4) ::  precip_len,  &  !length of precip running-mean (days)
                   precip_bef,  &  !time length before maximum precipitation (days)
                   precip_aft,  &  !time length after maximum precipitation (days)
                   tawftop_len, &  !time length water must be available (days)
                   tm_len          !time length temperature must be acceptable (days)

    real(r4) :: daylen_mini,   & !minimum day length if day length is increasing (hr)
                daylen_offd,   & !day length growing ability for decreasing day length (hr)
                         ! > 0: Offset from maximum day length
                         ! < 0: Daylength threshold to restrict day lengths (winter wheat)
                tawftop_min,  & !minimum water availability (-)
                tm_min,       & !minimum temperature (K)
                tm_max          !maximum temperature (K)

    !...growing season reset factors
    integer(i4) :: assim_resetl !running-mean length for assim mean
    real(r4) :: assim_resetv    !assimilation factor threshold 
                                !...to reset growing season variables

    !...phenology stage index parameters
    integer(i4) :: npstg  !number of phenological stages
    real(r4), dimension(:), allocatable ::  &  !(npstgmax-1)
         threshp  !thresholds for phenology stages

    !...phenology stage varying parameters
    real(r4), dimension(:,:), allocatable :: & !(npoolpft,npstgmax)
         allocp     !phenology-based allocation fractions
    logical :: &
          adj_moist, & !adjust allocations due to moisture
          adj_temp     !adjust allocations due to temperature

    real(r4), dimension(:), allocatable :: & !(npstgmax)
        lptransfer !phenology-based leaf pool (LAI) transfer
    real(r4), dimension(:), allocatable :: & !(npstgmax)
        vmax0      !phenology-based rubisco velocity of sun leaf (mol/m2/s)

    !---------------------------------
    !Dynamic Phenology Parameters
    !---------------------------------

    !...phenology factor parameters
    !.....day length potential
    real(r4) :: &
         psdayl_ref,  &  !day length reference
         psdayl_mul,  &  !daily change
         psdayl_min      !minimum value
    
    !.....growth potential
    !.....climatological suitability for growth
    !......using y=AB^clim_wa + C*(clim_wa-D)
    real(r4) :: &
         climp_a,   & !climatological suitability (climp) exponential adjustment
         climp_b,   & !climp exponential adjustment base
         climp_c,   & !climp multiplicative adjustment coefficient
         climp_d,   & !climp multiplicative adjustment offset
         climp_min, & !climp minimum
         climp_max    !climp maximum

    integer(i4) :: cwa_type !climatological suitability water availability type
    !   1=CUPR (convective precipitation)
    !   2=Total Precipitation
    !   3=PAWFRW (Plant Available Water Root-Weighted Fraction)
    !   4=TAWFRW (Total Available Water Root-Weighted Fraction)

    real(r4) :: &
         clai_coef, &  !climatological LAI coefficient (-)
         clai_offl, &  !climatological LAI lower offset (-)
         clai_offg     !climatological LAI upper offset (-)
    
    real(r4) :: psg_min !growth potential minimum (-)

    !.....weather potential
    integer(i4) :: pswx_rml  !weather potential running-mean length
    integer(i4) :: pswx_type !weather potential water availability type
       !    1=PAW_FTop
       !    2=PAW_FAll
       !    3=2*PAW_Fall
       !    4=PAW_FAll > 0. ==> 1.0
       !    5=TAW_FTop
       !    6=TAW_FAll
       !    7=RSTFAC2
       !    8=RSTFAC4    

    !-----------------------------
    !Defined Phenology Parameters
    !-----------------------------

    !...phenology stage determination method
    !.....1=Growing Degree Days (GDD)
    !.....2=Days After Planting Date (DAPD)
    integer(byte) :: gdd_or_pd

    !...growing degree day information
    real(r4) :: gdd_tbase  ! Base temperature (F)
    real(r4) :: gdd_tmax   ! Maximum temperature (F)
    
    !...growing season length (GSL) maximum (days)
    integer(i4) :: gslmax     !growing season length (GSL) max (days)

    !...seed information
    real(r4) :: seed_carbon   !carbon in seed (mol C)
    real(r4) :: seed_release  !daily carbon released from seed (mol C)

    !-------------------------
    !Calculated Parameters
    !------------------------

    !.....weights for growing season start
    real(r8) :: wt_assim,   & !time-step contribution for assimilation
                wt_precip,  & !time-step contribution for precipitation
                wt_tawftop, & !time-step contribution for tawftop
                wt_tm         !time-step contribution for temperature

    !.....weights for dynamic phenology
    real(r8) :: wt_pswx  !time-step contribution for weather potential

end type phen_param


!------------------------------------------------------------------
! Physiological Parameters 
!------------------------------------------------------------------
type, public :: phys_param

    !...vegetation-specific characteristics
    real(r4) :: sla      !specific leaf area (m2 leaf area/g leaf)
    real(r4) :: laimin   !m2/m2
    real(r4) :: laisat   !m2/m2
    real(r4) :: fparsat  !-
    logical  :: pftgraze !switch for grazing

    integer(byte)  :: c4flag  ! 0=c3,  1=c4
    real(r4) :: chil    ! leaf angle distribution factor (-)
    real(r4) :: z1      ! canopy bottom (m)
    real(r4) :: z2      ! canopy top (m)
    real(r4) :: kroot   ! root density extinction coefficient (-)
    real(r4) :: rootd   ! maximum rooting depth (m)

    !...temperature stress parameters
    real(r4) :: slti    ! slope of lo-temp inhibition (1/K)
    real(r4) :: shti    ! slope of hi-temp inhibition (1/K)
    real(r4) :: hlti    ! 1/2 point of lo-temp inhibition (K)
    real(r4) :: hhti    ! 1/2 point of hi-temp inhibition (K)
    real(r4) :: hfti    ! 1/2 point of frost inhibition (K)
    real(r4) :: sfti    ! slope of frost inhibition (1/K)

    !...soil moisture stress parameters
    real(r4), dimension(:), allocatable :: & !(nsoil/2)
                fc_min  ! minimum field capacity (1/m3)
    real(r4), dimension(:), allocatable :: & !(nsoil/2)
                wp_min  ! minimum wilting point (1/m3)
    real(r4) :: wssp    ! water stress shape parameter (0.1-1.0)

    !...photosynthesis parameters
    real(r4) :: effcon  ! quantum efficiency (mol CO2/mol quanta)
    real(r4) :: gmeso   ! mesophyll conductance (mol/m^2/sec) 
    real(r4) :: binter  ! conductance-photosynthesis intercept (mol m^-2 sec^-1)
    real(r4) :: gradm   ! conductance-photosynthesis slope parameter (-)
    real(r4) :: atheta  ! WC WE coupling parameter (-)
    real(r4) :: btheta  ! WC WE WS coupling parameter (-)
    real(r4) :: gmin    ! Minimum stomatal conductance (mol m^-2 sec^-1) - Table 1 median L17

    !...radiation parameters
    real(r4) :: tran(2,2) ! leaf transmittance (-)
    !  (1,1) - shortwave, green plants
    !  (1,2) - shortwave, brown plants
    !  (2,1) - longwave, green plants
    !  (2,2) - longwave, brown plants
    real(r4) :: ref(2,2)  ! leaf reflectance (-)
    !  (1,1) - shortwave, green plants
    !  (1,2) - shortwave, brown plants
    !  (2,1) - longwave, green plants
    !  (2,2) - longwave, brown plants

end type phys_param

!------------------------------------------------------------------
! Pool Parameters 
!------------------------------------------------------------------
type, public :: pool_param

    !----Live Pool Info----
    !...autotrophic respiration parameters
    real(r4), dimension(:), allocatable :: &  !(npoolpft)
         gr_frac  !growth respiration coefficients (0-1)

    real(r4), dimension(:), allocatable :: & !(npoolpft)
         lresp_eff  !respiration efficiency (0-1)
    
    real(r4) :: cr_aml,  & !canopy assimilation respiration low multiplier (-)
                cr_amh,  & !canopy assimilation respiration high multiplier (-)
                cr_amin, & !canopy assimilation respiration min (-)
                cr_amax    !canopy assimilation respiration max (-)

    real(r4) :: cr_fmul, & !canopy freeze respiration inhibition multiplier (-)
                cr_fref, & !canopy freeze respiration inhibition ref temperature (K)
                cr_fmin    !canopy freeze respiration inhibition minimum (-)

    real(r4) :: cr_hq10, & !Q10 base for canopy respiration hot scalar (-)
                cr_href, & !reference temperature for canopy respiration hot scalar (K)
                cr_hmax    !max value for canopy respiration high temp scalar (-)
    
    !...leaf transfer parameters
    real(r4) :: lt_fq10,  & !leaf freeze transfer Q10 base (-)
                lt_fref,  & !leaf freeze transfer reference temperature (K)
                lt_fmax     !leaf freeze transfer maximum (fraction per day)

    real(r4) :: lt_dcoef,  & !leaf daylength transfer coefficient (-)
                lt_dref,   & !leaf daylength transfer ref (diff from max daylenth)
                lt_dmax      !leaf daylength transfer maximum (fraction per day)
    
    real(r4) :: lt_wref,   & !leaf water deficiency transfer pawfrw reference (-)
                lt_wbase,  & !leaf water deficiency transfer exponential base (-)
                lt_wcoef,  & !leaf water deficiency transfer coefficient (-)
                lt_wmax      !leaf water deficiency transfer maximum (fraction per day)

    !...root respiration and transfer parameters
    real(r8) :: rrt_aml,  & !root resp/transfer assimilation respiration low multiplier (-)
                rrt_amh,  & !root resp/transfer assimilation respiration high multiplier (-)
                rrt_amin, & !root resp/transfer assimilation respiration min (-)
                rrt_amax    !root resp/transfer assimilation respiration max (-)

    real(r8) :: rrt_fmul, & !root resp/transfer freeze inhibition multiplier (-)
                rrt_fref, & !root resp/transfer freeze inhibition ref temperature (K)
                rrt_fmin    !root resp/transfer freeze inhibition minimum (-)

    real(r4) :: rrt_hq10, & !root resp/transfer hot Q10 (-)
                rrt_href, & !root resp/transfer hot reference temperature (K)
                rrt_hmax    !root resp/transfer max high temperature scalar (-)

    real(r4) :: rrt_laimin, & !root resp/transfer LAI ratio (LAI/Clim_LAI) min (-)
                rrt_laimax    !root resp/transfer LAI ratio max (-)

    !----Dead Pool Info----
    !...surface pool heterotrophic respiration and transfer parameters
    real(r4) :: hrt_sfc_aml,  & !assimilation respiration low multiplier (-)
                hrt_sfc_amh,  & !assimilation respiration high multiplier (-)
                hrt_sfc_amin, & !assimilation respiration min (-)
                hrt_sfc_amax    !assimilation respiration max (-)
    
    real(r4) :: hrt_sfc_fmul, &  !freeze inhibition multiplier (-)
                hrt_sfc_fref, &  !freeze inhibition ref temperature (K)
                hrt_sfc_fmin     !freeze inhibition minimum (-)

    real(r4) :: hrt_sfc_hq10, &  !high temperature Q10 (-)
                hrt_sfc_href, &  !high temperature reference temperature (K)
                hrt_sfc_hmax     !high temperature maximum scalar (-)

    real(r4) :: hrt_sfc_pml, & !precip inhibition low multiplier (-)
                hrt_sfc_pmin   !precip inhibition minimum (-)

    !...soil pool heterotrophic respiration/transfer parameters
    real(r4) :: hrt_soil_aml,  & !assimilation respiration low multiplier (-)
                hrt_soil_amh,  & !assimilation respiration high multiplier (-)
                hrt_soil_amin, & !assimilation respiration min (-)
                hrt_soil_amax    !assimilation respiration max (-)
    
    real(r4) :: hrt_soil_fmul, &  !freeze inhibition multiplier (-)
                hrt_soil_fref, &  !freeze inhibition ref temperature (K)
                hrt_soil_fmin     !freeze inhibition minimum (-)

    real(r4) :: hrt_soil_hq10, &  !high temperature Q10 (-)
                hrt_soil_href, &  !high temperature reference temperature (K)
                hrt_soil_hmax     !high temperature maximum scalar (-)

    real(r4) :: hrt_soil_mmin     !moisture inhibition minimum (-)
    real(r4) :: hrt_soil_pawmin   !PAW inhibition minimum (-)


    !...dead pool respiration efficiency (0-1)
    real(r4), dimension(:,:), allocatable :: &  !(npoollu,npoollu)
         dresp_eff

    !...grazing and harvest transfer parameters (npoollu+2)
    real(r4), dimension(:), allocatable :: &
         graze_trans, & ! transfer fractions for grazing (0-1)
         harvest_trans  ! transfer fractions for harvest (0-1)

    !----Live and Dead Pool Info----
    !...pool turnover times (yr)
    real(r4), dimension(:), allocatable ::  &  !(ntpool)
         turnover
    
    !...pool transfer fractions (0-1)
    real(r4), dimension(:,:), allocatable :: & !(ntpool,ntpool)
         pool_trans_frac    ! transfer fractions between pools (-)


    !...calculated parameters...
    real(r4), dimension(:), allocatable :: & !(ntpool)
         k_rate  !pool decay rate (1/s)

    real(r4), dimension(:), allocatable :: &  !(npoolpft)
          poolpft_min  !minimum pool values (mol C/m2)
    
end type pool_param

!***********************************************************************
!-----------------------------------------------------------------------
! Aerodynamic Parameters
  integer(i4) :: ngrid  !number of points in interpolation table
  type(aero_param), allocatable :: aerovar(:,:,:) !(npft,ngrid,ngrid)
  real(r4), allocatable :: LAIgrid(:)      !(ngrid)
  real(r4), allocatable :: fVCovergrid(:)  !(ngrid)

! Phenological Parameters (npft)
  type(phen_param),  dimension(:), allocatable :: phencon

! Physiological Parameters (npft)
  type(phys_param), dimension(:), allocatable :: physcon

! Pool Parameters (npft)
  type(pool_param), dimension(:), allocatable :: poolcon
!-----------------------------------------------------------------------

end module module_param
