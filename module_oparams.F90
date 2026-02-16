!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	
module module_oparams

!-----------------------------------------------------------------
!  This modules specifices all non-PFT specific parameters:
!     - Climtology
!     - COS
!     - Fluxes
!     - Fire
!     - Grazing
!     - Hydrology
!     - Phenology Allocation Adjustment 
!        (Friedlinstein et al.)
!     - Photosynthesis  
!     - Radiation
!     - SIF
!-----------------------------------------------------------------

       use kinds
       implicit none

       real(r4), parameter :: &
            near_zero= 1.E-10     ! Near zero value

       ! Climatological Parameters
       integer(i4), parameter :: &
            clim_len = 3650, &  !length of climatological averaging (days)
            seas_len = 10,   &  !length of seasonal averaging (days)
            seas_len_precip = 18  !length of precip seasonal averaging (days)
       
       ! COS Parameters
       integer(i4), parameter :: cos_i3=3
       real(r4), parameter :: &
            k_cos_soil = 1.2E4, & ! COS soil decay rate (1/s)
            cos_casd_min = 10.0   ! CAS depth minimum for COS (m)
       
       ! Constants for surface flux functions
       !  from Holtslag and Boville (1993, J. Climate)
       real(r4), parameter :: & 
            bunstablM = 10.0, &
            bunstablT = 15.0, &
            cunstablM = 75.0, &
            cunstablT = 75.0, &
            bstabl =  8.0,    &
            cstabl = 10.0      

      ! Fire Parameters
       real(r4), parameter :: &
            fire_leaff = 0.3, & ! fire emission fraction removed from leaf pool
            fire_stwdf = 0.1, & ! fire emission fraction removed from stwd pool
            fire_cdbf  = 0.2, & ! fire emission fraction removed from cdb pool
            fire_metlf = 0.2, & ! fire emission fraction removed from metl pool
            fire_strlf = 0.2    ! fire emission fraction removed from strl pool
            !fire_leaffc13 = 0.3, & ! fire emission fraction removed from leaf 13C pool
            !fire_stwdfc13 = 0.1, & ! fire emission fraction removed from stwd 13C pool
            !fire_cdbfc13  = 0.2, & ! fire emission fraction removed from cdb 13C pool
            !fire_metlfc13 = 0.2, & ! fire emission fraction removed from metl 13C pool
            !fire_strlfc13 = 0.2    ! fire emission fraction removed from strl 13C pool
       
       ! Grazing Parameters
       real(r4), parameter :: &
            graze_cfracp = 0.008, & ! Fraction of canopy carbon grazed daily 
                                   !   for productive systems
            graze_cfracd = 0.003, & ! Fraction of canopy carbon grazed daily
                                   !   for non-productive systems
            graze_minlai = 0.70, & ! Minimum LAI required for grazing (m2/m2)
            graze_climlai = 1.0    ! LAI threshold to switch between grazing
                                   !   for productive and desert ecosystems (m2/m2)
       
       ! Hydrological Parameters
       real(r4), parameter ::  &
            canlai_min = 0.01     ! Minimum LAI to use canopy water storage

       real(r4), parameter ::  &
            snow_c2 = 23.0e-3, &  ! Snow compaction parameter (m3/kg)
            snow_c3 = 2.77e-6, &  ! Snow compaction parameter (1/s)
            snow_c4 = 0.04,    &  ! Snow compaction parameter (1/K)
            snow_c5 = 2.0,     &  ! Snow compaction parameter (-)
            snow_dm = 100.0,   &  ! Uupper limit on destructive 
                                  !   metamorphism compaction (kg/m3)
            snow_eta0 = 9.0e5     ! Viscosity coefficient (kg/m3)

       real(r4), parameter :: &
            wsat_default = 0.95 ! Default initial water saturation fraction (-)

       ! Phenology Allocation Adjustment Parameters
       real(r4), parameter :: &
            aadjustmin = 0.02, &   ! Adjustment minimum (-)
            lftit = 273.0,     &   ! Leaf adjustment temperature ref (K)
            lftif = 1.3,       &   ! Leaf adjustment factor (-)
            lgrw_min = 0.6,    &   ! Leaf growth minimum (-)
            moistmult = 2.8,   &   ! Moisture adjustment multiplier
            wftit = 278.0,     &   ! Wood/Stem adjustment temperature ref (K)
            wftif = 1.3            ! Wood/Stem adjustment factor (-)

       ! Respiration Parameters
       real(r4) :: &
            rt_moist_exp = 30.0,  &  ! Moisture resp/transfer scalar exponent (-)
            rt_moist_range = 0.98    ! Moisture resp/transfer scalar range (-)

       ! Radiation Parameters
       real(r4), parameter :: &
            rad_c1 = 580., &
            rad_c2 = 464., &
            rad_c3 = 499., &
            rad_c4 = 963., &
            rad_c5 = 1160.
       real(r4), parameter :: &
            tcbot = 253.15

       ! SIF Parameters
       real(r4), parameter :: &
            mu_1 = -0.6,      &   ! SIF aerosol optical depth parameter   
            perih = 1.7963,   &   ! SIF TOA solar parameter
            sif_a1 = 1.93,    &   ! SIF parameter
            sif_a2 = 10.0,    &   ! SIF parameter
            sif_kf0 = 0.05,   &   ! SIF parameter
            sif_kn0 = 5.01,   &   ! SIF parameter
            sif_kp0 = 4.0        ! SIF parameter

end module module_oparams

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
