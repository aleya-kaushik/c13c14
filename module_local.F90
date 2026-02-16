module module_local

!----------------------------------------------------------------------
!   SiB4 Local Variables    
!----------------------------------------------------------------------

    use kinds
    use module_sibconst, only: nsoil, nsnow

    implicit none

    !...previous time-step values
    real(r8) :: capaccl_old      ! previous canopy liquid water storage (kg/m2)
    real(r8) :: capaccs_old      ! previous canopy snow water storeage (kg/m2)
    real(r8) :: capacg_old      ! previous ground liquid water storage (kg/m2)
    real(r8) :: snow_cmass_old  ! previous canopy snow cover (kg/m2)
    real(r8) :: snow_gmass_old  ! previous ground snow cover (kg/m2)

    integer(byte) :: nsl_old    ! prev timestep number of snow layers (-)
    real(r8), dimension(-nsnow+1:nsoil) :: &
          wwwliq_old, &  ! prev timestep soil liquid (kg/m2)
          wwwice_old, &  ! prev timestep soil ice (kg/m2)
          frac_iceold, & ! prev timestep ice fraction of total soil liquid
          td_old         ! prev timestep soil temperature (K)

    !...saturation vapor pressures
    real(r8) :: etc       ! saturation vapor pressure at Tc (hPa)
                          !   ('e-star' of Tc)
    real(r8) :: getc      ! derivative of etc with respect to temp
                          !   (d(etc)/dTc (hPa K^-1)
    real(r8) :: etg       ! 'e-star' of ground surface (hPa)
    real(r8) :: getg      ! d(etg)/dTg (hPa K^-1)
    real(r8) :: ets       ! 'e-star' of snow surface (hPa)
    real(r8) :: gets      ! d(ets)/dTs (hPa K^-1)

    !...canopy air space (CAS) and canopy
    real(r8) :: dtg       ! change in ground surface temperature (K)
    real(r8) :: dts       ! change in snow surface temperature (K)
    real(r8) :: dtc       ! change in canopy temperature (K)
    real(r8) :: dth       ! change in ref level temperature (K)
    real(r8) :: dqm       ! change in ref level moisture (hPa)
    real(r8) :: dta       ! change in CAS temperature (K)
    real(r8) :: dea       ! change in CAS moisture (hPa)

    real(r8) :: gect      ! dry fraction of veg / rc
    real(r8) :: geci      ! wet fraction of veg / 2rb
    real(r8) :: gegs      ! dry fraction of ground / rds
    real(r8) :: gegi      ! wet fraction of ground /rd
    real(r8) :: coc       ! gect + geci
    real(r8) :: cog1      ! gegi + gegs*hrr
    real(r8) :: cog2      ! gegi + gegs
    real(r8), dimension(1) :: &
       ppl, ttl, &  ! temp vars for humidity
       esst, dtesst ! temp vars for easat

    !...radiation
    real(r8) :: dtc4      ! d(canopy thermal em)/dT (W/m2/K)
    real(r8) :: dtg4      ! d(ground thermal em)/dT (W/m2/K)
    real(r8) :: dts4      ! d(snow thermal em)/dT   (W/m2/K)
    real(r8) :: lcdtc     ! d(canopy thermal em)/dtc (W/m2/K)
    real(r8) :: lcdtg     ! d(canopy thermal em)/dtg (W/m2/K)
    real(r8) :: lcdts     ! d(canopy thermal em)/dts (W/m2/K)
    real(r8) :: lgdtc     ! d(ground thermal em)/dtc (W/m2/K)
    real(r8) :: lgdtg     ! d(ground thermal em)/dtg (W/m2/K)
    real(r8) :: lsdts     ! d(snow thermal em)/dts   (W/m2/K)
    real(r8) :: lsdtc     ! d(snow thermal em)/dtc   (W/m2/K)
    real(r8) :: hcdtc     ! d(canopy H)/dtc  (W/m2/K)
    real(r8) :: hcdta     ! d(canopy H)/dta  (W/m2/K)
    real(r8) :: hgdta     ! d(ground H)/dta  (W/m2/K)
    real(r8) :: hgdtg     ! d(ground H)/dtg  (W/m2/K)
    real(r8) :: hsdta     ! d(snow H)/dta  (W/m2/K)
    real(r8) :: hsdts     ! d(snow H)/dtsnow  (W/m2/K)
    real(r8) :: hadta     ! d(CAS H)/dta  (W/m2/K)
    real(r8) :: hadth     ! d(CAS H)/dtheta  (W/m2/K)
    real(r8) :: ecdtc     ! d(canopy LE)/dtc (W/m2/K)
    real(r8) :: ecdea     ! d(canopy LE)/dea (W/m2/K)
    real(r8) :: egdtg     ! d(ground LE)/dtg (W/m2/K)
    real(r8) :: egdea     ! d(ground LE)/dea (W/m2/K)
    real(r8) :: esdts     ! d(snow LE)/dtsnow (W/m2/K)
    real(r8) :: esdea     ! d(snow LE)/dea (W/m2/hPa)
    real(r8) :: eadea     ! d(CAS LE)/dea (W/m2/hPa)
    real(r8) :: eadem     ! d(CAS LE)/dem (W/m2/hPa)
    real(r8) :: closs     ! canopy thermal loss     (W/m2)
    real(r8) :: gloss     ! ground thermal loss     (W/m2)
    real(r8) :: sloss     ! snow thermal loss       (W/m2)
    real(r8) :: fac1      ! effective ground cover for thermal radiation (-)

    !...soil/snow variables to save
    real(r8), dimension(-nsnow+1:nsoil) :: &
          dtd, & ! delta soil temperature (K)
          imelt  ! flag for 1=melting 2=freezing in soil column

    !...variables to save for energy balance
    real(r8) :: cas_e_storage   !CAS energy storage
    real(r8) :: cas_w_storage   !CAS water storage

    real(r8) :: radttc ! copy of PFT net radiation (W/m2)
    real(r8) :: radttg ! copy of ground net radiation (W/m2)
    real(r8) :: radtts ! copy of snow net radiation (W/m2)

end module module_local
