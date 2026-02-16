module module_phosib

!----------------------------------------------------------------------
!   SiB4 Photosynthesis Variables    
!----------------------------------------------------------------------

    use kinds

    implicit none

    !...Parameters
    integer(i4), parameter :: numic = 7
    real(r8), parameter :: minrad=1.e-6
    real(r8), parameter :: minassim=1.e-7
    real(r8), parameter :: minassimdiff=10.e-6

    real(r8), parameter :: toatosfc = 0.8
    real(r8), parameter :: ctoco2 = 44./12.
    real(r8), parameter :: convertppm = 1.e6
    real(r8), parameter :: convertphoton = 0.219

    !...Parameters - Base values for Q10 relationship
    real(r8), parameter :: vmax_q10=2.1   !VMax temperature sensitivity
    real(r8), parameter :: vmax_tref=298. !VMax reference temperature
    real(r8), parameter :: oms_q10=1.8    !Export-limited assimilation

    !...Parameters - Damping Variables
    real(r8), parameter :: zln2 = 6.9314718e-1
    real(r8), parameter :: ghalf = 1.0257068e1
    real(r8), parameter :: dttin = 3.6e3 
    real(r8), parameter :: dmin = 6.0e1

    !...Parameters - Limits
    real(r4), parameter ::  co2_casd_min = 4.0 ! CAS depth minimum for CO2 (m)
    real(r4), parameter ::  rst_max = 5.E6 !Max stomatal resistance (s/m)
    real(r4), parameter :: &
            rhfac_astart = 0.6,   & ! Humidity stress curvature start (-)
            rhfac_exp = 2.2,      & ! Humidity stress curvature (-)
            rhfac_exp_crop = 0.7, & ! Humidity stress crop exponent (-)
            rhfac_nforest = 0.7,  & ! Humidity stress min for needle forests (-)
            rhfac_tundra = 0.6      ! Humidity stress min for tundra (-)


    !===============================================================================
    !...Assimilation-Limiting Variables
    real(r8) :: assim_omc ! rubisco-limited assimilation (mol C/m2/s)
    real(r8) :: assim_ome ! light-limited assimilation (mol C/m2/s)
    real(r8) :: assim_oms ! sink-limited assimilation (mol C/m2/s)
    real(r8) :: assimpot_omc !potential rubisco-unlimited assimilation (mol C/m2/s)
    real(r8) :: assimpot_ome !potential light-unlimited assimilation (mol C/m2/s)
    real(r8) :: assimpot_oms !potential sink-unlimited assimilation (mol C/m2/s)
    real(r8) :: assimfac(4)  !assimilation rate stress factors (-)
                             ! (1) rubisco-limited stress
                             ! (2) light-limited stress
                             ! (3) sink-limited stress
                             ! (4) product of factors 1-3
    
    real(r8) :: assimdiff ! diff btwn actual & potential assim (mol/m2/s)
    real(r8) :: rrkk      ! export-limited rate portion 
    real(r8) :: omss      ! export-limited rate portion 

    !...Assimilation/Ball-Berry Variables
    real(r8) :: bintc     ! Minimal assimilation (y-intercept of Ball-Berry; mol/m2/s)
    real(r8) :: range     ! Adjusted/First Guess choloroplast CO2 (Pa)
    real(r8) :: assimn    ! Net leaf assimilation (mol/m2/s)
    real(r8) :: eyy(numic)     ! Iterated difference between pco2 estimates (Pa)
    real(r8) :: pco2y(numic)   ! Iterated chloroplast CO2 (Pa)
    real(r8) :: assimc(numic)  ! Iterated rubisco-limited assimilation (mol/m2/s)
    real(r8) :: assime(numic)  ! Iterated light-limited assimilation (mol/m2/s)
    real(r8) :: assims(numic)  ! Iterated sink-limited assimilation (mol/m2/s)
    real(r8) :: assimy(numic)  ! Iterated assimilation (mol/m2/s)
    real(r8) :: gminc      ! Minimum stomatal conductance (mol/m2/s)

    !...Atmospheric Variables
    real(r8) :: co2cap   !air capacity for CO2 exchange (mol air/m2)
    real(r8) :: pressure !surface pressure (mb)

    !...Conductance Variables
    real(r8) :: gxco2     ! Canopy biophysical-maximum conductance (mol/m2/s)
    real(r8) :: gah2o     ! CAS-Mixed Layer conductance (mol/m2/s)
    real(r8) :: gbh2o     ! Leaf-CAS conductance (mol/m2/s)
    real(r8) :: gsh2o     ! Canopy conductance (mol/m2/s) 
    real(r8) :: gsh2onew  ! Updated canopy conductance (mol/m2/s)
    real(r8) :: drst      ! Delta of stomatal resistance (s/m)

    !...CO2 Concentrations
    real(r8) :: co2cas  ! CAS CO2 concentration (mol C/mol air)
    real(r8) :: co2m    ! reference level CO2 concentration (mol C/mol air)
    real(r8) :: co2s    ! leaf surface CO2 concentration (mol C/mol air)
    real(r8) :: co2i    ! leaf internal CO2 concentration (mol C/mol air)
    real(r8) :: co2c    ! leaf chloroplast CO2 concentration (mol C/mol air)
    real(r8) :: co2gamma    ! CO2 compensation point (mol C/mol air)

    real(r8) :: pco2cas ! CAS CO2 concentration (Pa)
    real(r8) :: pco2s   ! leaf surface CO2 partial pressure (Pa)
    real(r8) :: pco2i   ! leaf internal (stomatal) CO2 partial pressure (Pa)
    real(r8) :: pco2c   ! leaf chloroplast CO2 partial pressure (Pa)
    real(r8) :: pco2m   ! mixed layer CO2 partial pressure (Pa)

    !...Damping Factors
    real(r8) :: pdamp
    real(r8) :: qdamp
    real(r8) :: tprcor    ! Temperature correction (K)
 
    !...Leaf Parameters
    real(r8) :: qt        ! Temperature Q10 factor (-)
    real(r8) :: vmax_unscaled ! LK: Vmax not scaled with T, for update in cos_calc.F90
    real(r8) :: vmaxts    ! VMax temperature-scaled (mol/m2/s)
    real(r8) :: zkc       ! Michaelis-Menten coefficient for CO2 (Pa)
    real(r8) :: zko       ! Inhibition coefficient for O2 (Pa)
    real(r8) :: spfy      ! CO2/O2 specificity (-)

    !...Potential Photosynthesis Variables
    real(r8) :: assimnp   ! Top leaf assimilation (mol C/m2/s)
    real(r8) :: antemp    ! Bottom stopped assimimlation (mol C/m2/s)
    real(r8) :: pco2ipot  ! Potential intercellular CO2 (Pa)
    real(r8) :: omcpot    ! Potential rubisco limitation (mol C/m2/s)
    real(r8) :: omepot    ! Potential light limitation (mol C/m2/s)
    real(r8) :: omspot    ! Potential sink or pep limitation (mol C/m2/s) 

    real(r8) :: sqrtin    ! Intermediate potential (mol C/m2/s)
    real(r8) :: omppot    ! Intermediate top leaf photosynthesis (mol C/m2/s)

    !...Radiation and Canopy Scaling Factors
    real(r8) :: pfd       ! Photon flux (mol/m2/s)
    real(r8) :: scatc     ! scattering coefficient (-)
    real(r8) :: toa_pfd   ! Top-of-Atmosphere (TOA) photon flux (mol/m2/s)
    real(r8) :: toa_apar  ! Absorbed TOA PAR (mol/m2/s)

    !...Respiration Rates
    real(r8) :: resp_cas  ! total respiration, auto + heterotrophic (mol C/m2/s)

    !...Vegetation Info
    real(r8) :: c3, c4, atheta, btheta

end module module_phosib
