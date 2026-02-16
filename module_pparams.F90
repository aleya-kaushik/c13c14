!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	
module module_pparams

!-----------------------------------------------------------------
!	This modules specifices physical parameters using MKS.
!-----------------------------------------------------------------

        use kinds
        implicit none

        real(r8), parameter :: & 
            pi           = 3.14159265358979323846, &  !pi
            pi180        = pi / 180.0_r8,       & ! degrees to radians
            pidaypy      = 0.0172142,           & ! constant for zenith angle
            grav         = 9.81,                & ! earth gravity (m/s^2)
            gas_const_R  = 287.0,               & ! gas constant for dry air (J/kg/K)
            spec_heat_cp = 1005.0,              & ! specific heat at constant pressure (J/kg/K)
            kappa  = gas_const_R/spec_heat_cp,  & ! constant for pot temp conversion
            vkrmn        = 0.35,                & ! Von Karmann's constant (unitless)
            rgfac        = 100.0/gas_const_r,   & ! 
            delta        = 0.608,        & ! molecular_weight_air/molecular_weight_water - 1
            stefan       = 5.67e-08,     & ! stefan boltzmann constant (W/m^2/K^4)
            rv           = 4.61e+02,     & ! gas constant for water vapor
            rstar        = 8.3143,       & ! universal gas constant (m3 Pa/mol/K)
            amagat       = 44.6,         & ! reciprocal of molar volume of a gas at p0/tice (mol/m3)
            amagatwv     = 44.032476,    & ! water-vapor adjusted amagat (mol/m3)
            molc_h2o     = 55.56           ! number of moles per liter of water (mol/l)

       real(r8), parameter :: &
            par_conv     = 4.6E-6,  & ! PAR conversion factor (umol quanta/J)
            drytoc       = 1./0.45    ! dry mass to carbon ratio

       real (r8), parameter  ::     &
            lvap  = 2.25e+06,       & ! latent heat of vaporization (J/kg)
            lfus  = 0.3336e+06,     & ! latent heat of fusion
            lsub  = lvap + lfus,    & ! latent heat of sublimation
            snomel  = 3.705185e+08, & ! latent heat of fusion of ice (J/m3) 
            snofac = lvap/(lvap + snomel * 1.E-3), & ! ratio from Sellers (1986)
            cv     = 1952.0,        & ! specific heat of water vapor at constant pressure (J/deg/kg)
            cpice  = 2117.27,       & ! specific heat of ice (J/deg/kg)
            cpliq  = 4188.0,        & ! specific heat of water (J/deg/kg)
            denh2o = 1000.0,        & ! density of water (kg/m^3)
            denice = 917.0,         & ! density of ice (kg/m^3) 
            tkair  = 0.023,         & ! thermal conductivity of air (W/m/K)
            tkwat  = 0.6,           & ! thermal conductivity of water (W/m/K)
            tkice  = 2.29,          & ! thermal conductivity of ice (W/m/K)
            leafhc = 4.186*1000.0*0.2,    & ! leaf heat capacity  (J/deg/m3)
            h2ohc  = 4.186*1000.0*1000.0    ! water heat capacity (J/deg/m3)

      real(r8), parameter ::  &
            phmin  = -1.e8,   & ! minimum value for soil potential (mm)
            wimp   = 00.05,   & ! water impermeable if porosity below this value
            wpotfc = -15.00,  & ! water potential at field capacity (J/kg)
            wpotwp = -1500.,  & ! water potential at wilting point (J/kg)
            wtfact = 00.30,   & ! fraction of area with high water table (-)
            ssi    = 0.033,   & ! irreducible water fraction of snow (-)
            cwlim  = 0.001,   & ! canopy water storage limit (kg/m2)
            gwlim  = 10.00,   & ! ground water storage limit (kg/m2)
                                !.....equivalent to 0.01 m 
            gwctog = 0.25,    & ! adjustment for water storage from point to grid scale
            zlnd   = 0.01,    & ! roughness length for land (m)
            psb    = 50.0       ! pbl mass depth (hPa or mb)

       real (r8), parameter::   &
            psref  = 1013.246,  & ! reference pressure (hPa)
            p0_sfc = 1.0e+05,   & ! surface pressure (Pa)
            tsref  = 373.15,    & ! reference temperature (K)
            tref   = 273.15,    & ! reference temperature (K)
            tice   = 273.15,    & ! freezing temperature of water (K)
            tbgmin = 253.15,    & ! minimum temperature for snow calculations (K)
            tbgmax = 273.15,    & ! maximum temperature for snow calculations (K) 
            twmin  = 173.16,    & ! lowest allowed temperature boundary for water (K)
            twmax  = 373.16,    & ! highest allowed temperature boundary for water (K)
            timin  = 173.16,    & ! lowest allowed temperature boundary for ice (K)
            timax  = 273.16,    & ! highest allowed temperature boundary for ice (K)
            tffrz  = 32.0         ! freezing temperature (F)

       real(r8), parameter :: &
            mwc    = 12.0,       & ! molecular weight of carbon (g/mol)
            mwc13  = 13.0,       & ! molecular weight of carbon-13 (g/mol)
            po2m   = 20900.0,    & ! mixed layer O2 concentration
            bco2m  = 340.0,      & ! mixed layer CO2 concentration (ppm)
            bcosm  = 500.0         ! mixed layer COS concentration (ppt)

       real(r8), parameter :: &
           pdb       = 0.0112372,  & ! 13C/12C ratio of Pee Dee Belemnite
           kieclfbl  = -2.9,       & ! canopy air space to leaf boundary layer
           kiecstom  = -4.4,       & ! leaf boundary layer to stomatal cavity
           kieclphas = -0.7,       & ! liquid phase fractionation
           kiecdis   = -1.1,       & ! dissolution
           kiecrbsco = -28.2,      & ! C3 C-fixation enzyme rubisco
           kiecphtrsp = -12.0,     & ! photorespiration, Keeling2017: 12 +/- 4
                                     ! appropriate for leaves at 25 deg C
           rpoolinitc3 = 0.0108324381949482, &  ! initial c3 pool allocation based on d13C=-26
                                     ! -26 = (R_spl/R_pdb -1)*1000.
                                     ! -26 = (rpoolinit/pdb-1)*1000.
           rpoolinitc4 = 0.0109808123837865   ! initial c4 pool allocation based on d13C=-12.4                 

       integer(i4), parameter :: &
            eqnx  = 80, &             !  day of vernal equinox
            nhsolstice = 172, &       !  longest day of year in NH
            shsolstice = 355          !  longest day of year in SH

       real(r8), parameter :: &
            solar_const = 1367.0,   & ! solar constant (W/m2)
            perhl       = 102.7,    & ! factor for zenith angle (-)
            eccn        = 0.016715, & ! eccentricity
            decmax      = 23.441,   & ! max declination
            cosz_min    = -0.1045     ! minimum cosine of zenith angle value
                                      !(96 deg. which includes civil twilight)

      real(r8), parameter :: &
            mol_to_dw = 24.0,   &  ! conversion factor: mol C/m2 to g DW/m2
            mol_to_mg = 0.12,   &  ! conversion factor: mol C/m2 to Mg C/ha
            molc13_to_mg = 0.13,   &  ! conversion factor: mol C-13/m2 to Mg C/ha
            mol_to_umol = 1.0D6, &  ! conversion factor: mol to micromol
            mol_to_pmol = 1.0D12    ! conversion factor: mol to picomol

      real(r8), parameter :: &
            mol_to_bu_mze = 4.248, & ! crop conversion: mol C/m2 to bu/ac
            mol_to_bu_soy = 3.965, & ! crop conversion: mol C/m2 to bu/ac
            mol_to_bu_wwt = 3.965    ! crop conversion: mol C/m2 to bu/ac

       integer(i4), parameter ::  &
            days_per_year = 365,  &    ! # days in non-leap year
            hrs_per_day   = 24,   &    ! # hours in a day
            secs_per_hr  = 3600,  &    ! # seconds in an hour
            secs_per_day = 86400, &    ! # seconds in a day
            ave_days_per_month = 30    ! # of days per month on average

       integer(i4), parameter, dimension(12) :: &
            days_per_month = (/31,28,31,30,31,30,31,31,30,31,30,31/),  &  
            doy1_month = (/1,32,60,91,121,152,182,213,244,274,305,335/)
       character(len=10), parameter, dimension(12) :: &
            month_names = (/'January   ','February  ','March     ',  &
                            'April     ','May       ','June      ',  &
                            'July      ','August    ','September ',  &
                            'October   ', 'November  ','December  '/)

end module module_pparams

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
