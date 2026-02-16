!=====================SUBROUTINE HYDRO_CANOPY===========================

subroutine hydro_nveg( &
           tm, tcas, &
           cuprt, lsprt,  &
           hydrost, sscolt)

use kinds
use module_local, only: frac_iceold
use module_pparams, only: &
    tice, tkice, denice, denh2o, &
    cpice, cpliq
use module_sib, only: &
    hydros_type, sscol_type
use module_time, only: dtsib, dtisib


implicit none

!=======================================================================
!     RESET ALL CANOPY HYDROLOGICAL VARIABLES.
!     SET ALL PRECIPITATION AS THROUGHFALL
!     UPDATE GROUND PRECIPITATION
!=======================================================================


!...input variables
real(r8), intent(in) :: tm, tcas
real(r8), intent(inout) :: cuprt, lsprt
type(hydros_type), intent(inout) :: hydrost
type(sscol_type), intent(inout) :: sscolt

!...local variables
real(r8) :: totalp    ! total precip (meters)
real(r8) :: p0, pinf
real(r8) :: xsc   ! canopy excess water (kg/m2)
real(r8) :: xss   ! canopy excess snow (kg/m2)
real(r8) :: cliqmo   ! veg liquid store in meters
real(r8) :: csnowmo ! veg snow in meters
real(r8) :: satcapm ! canopy store saturation in meters

!...water and vegatation (J/m2/deg)
real(r8) :: tti       ! direct throughfall in meters
real(r8) :: tex       ! canopy drainage in meters
real(r8) :: thru      ! total throughfall (tti + tex)

!...interception
real(r8) :: fliq        ! fraction of liquid in precip (0-1)
real(r8) :: bifall      ! density of falling snow (kg/m3)
real(r8) :: dz_snowfall ! new snow depth (m)
integer(byte) :: pcptype    ! precip type; 1=rain 2=snow

    !----------------------------------------------------
    !...SET INCOMING PRECIP
    totalp  = (cuprt + lsprt) *  dtsib  !meters
    if (totalp < 1.0E-10) then
       totalp = dzero
       cuprt = dzero
       lsprt = dzero
    endif
    p0 = totalp

    !...REMOVE EXCESS LEAF WATER/SNOW FROM PREVIOUS TIMESTEP
    cliqmo = hydrost%capacc_liq/denh2o    !meters
    csnowmo = hydrost%capacc_snow/denice  !meters
    satcapm = dzero

    xsc = max(dzero, cliqmo - satcapm)
    xss = max(dzero, csnowmo - satcapm)
    thru = xsc + xss

    !...TOTAL THROUGHFALL
    tti = p0
    tex = dzero
    thru = thru + tti + tex
    pinf = p0 - tti - tex

    if (pinf .lt. -1.E-10) then
        print*,'Precipitation Imbalance', pinf
        print*,'Incoming Precip  : ', p0
        print*,'Throughfall (tti): ', tti
        print*,'Leaf Runoff (tex): ', tex
        stop
     endif

    !...ASSIGN VARIABLES
    p0 = thru
    p0 = p0 * dtisib * 1000.  !convervion from m to mm/sec
    hydrost%p0 = p0

    hydrost%capacc_snow = dzero
    hydrost%capacc_liq = dzero
    hydrost%wetfracc = dzero

    !...precipitation onto ground (follows clm_hydro_canopy)
    !.....determine precip type
    if (tm >= tice) then
        pcptype = 1
        fliq = 1.0
        hydrost%pcpg_snow = dzero
        hydrost%pcpg_rain = hydrost%p0
    else
        pcptype = 2
        fliq = dzero

        if (tm > tice) then
           fliq = 0.40
           bifall = 189.
        elseif (tm > (tice - 15.0)) then
            fliq = -54.61 + 0.2*tm
            bifall = 50.0 + 1.7 * (tm - tice + 15.0)**1.5
        else
            fliq = dzero
            bifall = 50.0
        endif
        fliq = max(dzero,fliq)

        hydrost%pcpg_rain = hydrost%p0 * fliq
        hydrost%pcpg_snow = hydrost%p0 * (done - fliq)

        !...snowfall rate will be in m/sec
        dz_snowfall = hydrost%pcpg_snow/bifall 

        if (dz_snowfall > dzero) then
           !...snow_depth will be in meters...
           hydrost%snow_gdepth = hydrost%snow_gdepth + dz_snowfall * dtsib
           hydrost%snow_gmass = hydrost%snow_gmass + hydrost%pcpg_snow * dtsib

           if (sscolt%nsl < 0) then  !Add Snow To Existing Layer
              !...add snow depth to top snow layer
              sscolt%www_ice(sscolt%nsl+1) = &
                    sscolt%www_ice(sscolt%nsl+1) + dtsib * hydrost%pcpg_snow
              sscolt%dz(sscolt%nsl+1) = sscolt%dz(sscolt%nsl+1) + &
                    dz_snowfall * dtsib
           else   !Initialize a New Snow Layer
              sscolt%nsl    = -1
              sscolt%dz(0)  = hydrost%snow_gdepth
              sscolt%node_z(0)  = 0.5*sscolt%dz(0)
              sscolt%layer_z(0) = dzero
              sscolt%layer_z(-1) = sscolt%dz(0)
              sscolt%td(0) = MIN(tice, tcas)
              sscolt%www_ice(0) = hydrost%snow_gmass
              sscolt%www_liq(0) = dzero
              frac_iceold(0) = 1.0
              sscolt%shcap(0) = cpliq*sscolt%www_liq(0) + &
                  cpice*sscolt%www_ice(0)
              sscolt%tksoil(0) = tkice
          endif !new snow layer
        endif !snow layer check
    endif !snow precip

    !...end with p0 in mm (per timestep)
    hydrost%p0 = hydrost%p0 * dtsib

end subroutine hydro_nveg
