!=====================SUBROUTINE HYDRO_CANOPY===========================

subroutine hydro_canopy( &
     gref, glon, glat, pref,  &
     chil, tm, tcas,      &
     lai, vcover, &
     cuprt, lsprt, tc,  &
     hydrost, sscolt)

use kinds
use module_local, only: frac_iceold
use module_pparams, only: &
    tice, tkice, denice, denh2o, snomel, &
    cpice, cpliq, h2ohc, leafhc, month_names
use module_sibconst, only: &
    nsoil, nsnow, &
    canb_print, canb_stop, canb_thresh, &
    badtc_print, badtc_stop
use module_sib, only: &
    hydros_type, sscol_type
use module_time, only: dtsib, dtisib, &
    day, month, year, sec_tot


implicit none

!=======================================================================
!
!     CALCULATION OF INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW,
!     INCORPORATING EFFECTS OF PATCHY SNOW COVER AND TEMPERATURE
!     ADJUSTMENTS.
!
!----------------------------------------------------------------------
!
!     NON-UNIFORM PRECIPITATION DESCRIBED BY RELATIONSHIP:
!        F(X) = A*EXP(-B*X)+C
!
!     THROUGHFALL, INTERCEPTION AND INFILTRATION EXCESS 
!         ARE BASED ON THIS RELATIONSHIP.
!
!     REFERENCES:
!         1) SA-96A, APPENDIX D.
!         2) SA-89B: Sato et al, Implementing the Simple Biosphere Model 
!             (SiB) in a General Circulation Model: Methodologies and Results
!             NASA Contractor Report 185509 (1989)
!
!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       TC             CANOPY TEMPERATURE (K)
!       CAPACC_LIQ     CANOPY LIQUID INTERCEPTION STORE (kg m^-2)
!       CAPACC_SNOW    CANOPY SNOW INTERCEPTION STORE (Kg M^-2) or (mm water)
!       P0             THROUGHFALL PRECIPITATION (MM/HR)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!...input variables
integer(i4), intent(in) :: gref, pref
real(r4), intent(in) :: glon, glat
real(r4), intent(in) :: chil
real(r8), intent(in) :: tm, tcas
real(r8), intent(in) :: lai, vcover
real(r8), intent(inout) :: cuprt, lsprt
real(r8), intent(inout) :: tc
type(hydros_type), intent(inout) :: hydrost
type(sscol_type), intent(inout) :: sscolt

!...local variables
real(r8) :: pcoefs(2,2)
real(r8) :: ap, bp, cp
real(r8) :: totalp    ! total precip (meters)
real(r8) :: p0, pinf, fpi
real(r8) :: xsc   ! canopy excess water (kg/m2)
real(r8) :: xss   ! canopy excess snow (kg/m2)
real(r8) :: xsb   ! canopy excess combined snow/water (kg/m2)
real(r8) :: xc    ! canopy excess water (kg/m2)
real(r8) :: xs    ! canopy excess snow (kg/m2)
real(r8) :: cliqm, cliqmo   ! veg liquid store in meters
real(r8) :: csnowm, csnowmo ! veg snow in meters
real(r8) :: satcapm ! canopy store saturation in meters
real(r8) :: spechc   ! specific heat of canopy (intercepted)

!...water and vegatation (J/m2/deg)
real(r8) :: chiv      ! leaf angle dist factor (unitless)
real(r8) :: aa, bb
real(r8) :: exrain, zload
real(r8) :: tti       ! direct throughfall in meters
real(r8) :: tex       ! canopy drainage in meters
real(r8) :: thru      ! total throughfall (tti + tex)
real(r8) :: freeze
real(r8) :: diff
real(r8) :: ccp, cct
real(r8) :: newtc     ! canopy temperature modified for 
                      ! canopy precip coverage (K)
real(r8) :: tsd       ! temperature for freezing/melting fractions

!...interception
real(r8) :: tta, ttb
real(r8) :: cca, ccb, ccc
real(r8) :: arg
real(r8) :: fliq        ! fraction of liquid in precip (0-1)
real(r8) :: bifall      ! density of falling snow (kg/m3)
real(r8) :: dz_snowfall ! new snow depth (m)
integer(byte) :: pcptype    ! precip type; 1=rain 2=snow

!...balance check
real(r8) :: hbal, cbal, pbal
real(r8) :: hcpgaind, hcstord, hclossd
real(r8) :: hgpgaind, hgtgaind, hggaind

    data pcoefs(1,1)/ 20. /, pcoefs(1,2)/ .206E-8 /, &
         pcoefs(2,1)/ 0.0001 /, pcoefs(2,2)/ 0.9999 /, bp /20. /

    !----------------------------------------------------
    !...REMOVE EXCESS LEAF WATER/SNOW FROM PREVIOUS TIMESTEP
    cliqmo = hydrost%capacc_liq/denh2o    !meters
    csnowmo = hydrost%capacc_snow/denice  !meters
    satcapm = hydrost%satcapc/denh2o     !meters

    xsc = max(dzero, cliqmo - satcapm)
    cliqm = cliqmo - xsc

    xss = max(dzero, csnowmo - satcapm)
    csnowm = csnowmo - xss
    zload = cliqm + csnowm

    xsb = max(dzero, zload - satcapm)
    if (xsb .gt. dzero) then
       if (tm < tice) then
          cliqm = cliqm - xsb
       else
          csnowm = csnowm - xsb
       endif
       if (cliqm .lt. dzero) then
          csnowm = csnowm - cliqm
          cliqm = dzero
       elseif (csnowm .lt. dzero) then
          cliqm = cliqm - csnowm
          csnowm = dzero
       endif
       zload = cliqm + csnowm
    endif

    thru = xsc + xss + xsb

    spechc = &  !units of J/m2/K
        min( 0.05, zload)*h2ohc + lai*leafhc

    !...Set precip coefficients
    !......distributing large-scale and convective
    !......precip according to Sato et al. (1989B; eq. C3)
    dz_snowfall = 0.0
    ap = pcoefs(2,1)
    cp = pcoefs(2,2)
    totalp  = (cuprt + lsprt) *  dtsib  !meters
    if (totalp < 1.0E-10) then
       totalp = dzero
       cuprt = dzero
       lsprt = dzero
    endif
    p0 = totalp

    if (hydrost%capacc_snow > 0.0 .or. tm < tice) cuprt = 0.0
    lsprt = totalp/dtsib - cuprt

    if (totalp > 1.e-8) then
        ap = cuprt * dtsib / totalp * pcoefs(1,1) + &
             lsprt * dtsib / totalp * pcoefs(2,1)
        cp = cuprt * dtsib / totalp * pcoefs(1,2) + &
             lsprt * dtsib / totalp * pcoefs(2,2)
    endif

    !----------------------------------------------------
    !...PROPORTIONAL SATURATED AREA (XS) AND
    !...DIRECT THROUGHFALL (TTI) AND LEAF DRAINAGE (TEX)
    !...TTI ( D-D )     : EQUATION (C.4), SA-89B
    !...XS  ( X-S )     : EQUATION (C.7), SA-89B
    !...TEX ( D-C )     : EQUATION (C.8), SA-89B

    !TTI....direct throughfall in meters
    chiv = chil
    if (abs(chiv) <= 0.01) chiv=0.01

    aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv   ! unitless
    bb = 0.877 * ( 1. - 2. * aa )                  ! unitless
    exrain = aa + bb                               ! unitless
    
    fpi   = ( 1.-exp( - exrain*lai/vcover )) * vcover
    tti   = p0 * ( 1. - fpi )  !rainfall through to ground (m)

    !TEX...leaf drainage in meters
    xs    = 1.
    if ( p0 >= 1.e-9 ) then  
        arg = ( satcapm-zload )/ &
                ( p0*fpi*ap ) - cp/ap 

        if ( arg >= 1.e-9 ) then                                 
            xs = -1./bp * log( arg )                                      
            xs = min( xs, done ) 
            xs = max( xs, dzero ) 
        endif   
    endif 

    tex = p0 * fpi * &
        ( ap/bp*( 1.- exp( -bp*xs )) + cp*xs ) - &
        ( satcapm - zload ) * xs                        
    tex = max( tex, dzero ) 

     !  TOTAL THROUGHFALL AND STORE AUGMENTATION
     thru = thru + tti + tex
     pinf = p0 - tti - tex

     if (pinf .lt. -1.E-10) then
        print*,'Precipitation Imbalance', pinf
        print*,'Incoming Precip  : ', p0
        print*,'Throughfall (tti): ', tti
        print*,'Leaf Runoff (tex): ', tex
        stop
     endif

     !--------------------------------------------------------
     !...TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION
     if ( tm > tice ) then
        cliqm = cliqm + pinf
     else
        csnowm = csnowm + pinf
     endif    
     
     !...check storages
     if (cliqm .lt. dzero) then
        csnowm = csnowm + abs(cliqm)
        cliqm = dzero
     endif
     if (csnowm .lt. dzero) then
         cliqm = cliqm + abs(csnowm)
         csnowm = dzero
     endif

    freeze = dzero
    diff = ( cliqm + csnowm - cliqmo - csnowmo )*h2ohc
    ccp = spechc
    cct = spechc + diff
    newtc = ( tc * ccp + tm * diff ) / cct
    if (( tc >  tice .AND. tm <= tice ) .or. &
        ( tc <= tice .AND. tm >  tice ) ) then
        tta = tc
        ttb = tm
        cca = ccp
        ccb = diff
        if ( newtc <= tice ) then

            !...FREEZING OF WATER ON CANOPY 
            ccc = cliqm * 0.001 * snomel
            if ( tc < tm ) ccc = diff * snomel / h2ohc
            tsd = ( tta * cca + ttb * ccb + ccc ) / cct
            freeze = ( tice * cct - ( tta * cca + ttb * ccb ) )
            freeze = (min ( ccc, freeze )) / snomel
            if (tsd > tice) tsd = tice - 0.01
        else

            !...MELTING OF SNOW ON CANOPY 
            ccc = -csnowm * snomel 
            if ( tc > tm ) ccc = - diff * snomel / h2ohc

            tsd = ( tta * cca + ttb * ccb + ccc ) / cct

            freeze = ( tice * cct - ( tta * cca + ttb * ccb ) )
            freeze = (max( ccc, freeze )) / snomel
            if (tsd <= tice) tsd = tice - 0.01
        endif
    endif

    csnowm = csnowm + freeze
    cliqm  = cliqm - freeze
    if (csnowm .lt. dzero) then
       cliqm = cliqm + csnowm
       csnowm = dzero
    elseif (cliqm .lt. dzero) then
       csnowm = csnowm + cliqm
       cliqm = dzero
    endif

    xc = max( dzero, (cliqm - satcapm))
    cliqm = cliqm - xc

    xs = max(dzero, (csnowm - satcapm))
    csnowm = csnowm - xs
    thru = thru + xs + xc

    !...assign variables
    !...update canopy interception storage
    p0 = thru
    p0 = p0 * dtisib * 1000.  !convervion from m to mm/sec
    hydrost%p0 = p0

    hydrost%capacc_snow = csnowm * denice
    hydrost%capacc_liq = cliqm * denh2o

    !...update canopy wetness fraction
    hydrost%wetfracc = MAX(dzero, MIN( done, &
         (hydrost%capacc_liq + hydrost%capacc_snow) &
         / hydrost%satcapc ))
    
    !...check balance
    hcpgaind = pinf
    hcstord = (cliqm - cliqmo) + (csnowm - csnowmo)
    hclossd = xsc + xss + xsb + xc + xs

    hgpgaind = tti + tex
    hgtgaind = thru - tti - tex
    hggaind = hgpgaind + hgtgaind

    hbal = totalp - hcstord - hggaind
    cbal = pinf - hcstord - hclossd
    pbal = totalp - hcpgaind - hgpgaind
    
    if ((canb_print) .or. &
        (abs(hbal) .gt. canb_thresh) .or. &
        (abs(cbal) .gt. canb_thresh) .or. &
        (abs(pbal) .gt. canb_thresh)) then
        print*,'CBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCB'
        print*,'Canopy Hydrological Imbalance'
        print*,'Imbalance (m): ',hbal
        print('(a,a,i3,a,i4)'), &
              'Date: ', trim(month_names(month)), day, ', ', year
        print'(a,i12)','Second=',sec_tot
        print'(a,2f12.4,1i7,2i4)','Point (Lon/Lat/Ref/PFT): ', &
                 glon, glat, gref, pref
        print'(a,f8.4)',    'Current LAI: ',lai
        print*,''
        print*,'Precip Imbalance : ', pbal
        print*,'   totalp, canopy, ground'
        print*,'   ', totalp, pinf, tti + tex
        print*,''
        print*,'Canopy Imbalance: ', cbal
        print*,'Precip Input: ',pinf
        print*,'Storage Changes: ', hcstord
        print*,'  Liq (change/new/orig):',cliqm-cliqmo,cliqm,cliqmo
        print*,'  Ice (change/new/orig):',csnowm-csnowmo,csnowm,csnowmo
        print*,'Loss to Ground: ',hclossd
        print*,'   xsc,xss,xsb: ',xsc,xss,xsb
        print*,'   xc,xs  : ',xc,xs
        print*,'   tex    : ',tex
        print*,'Phase Change Exchange: ',freeze
        print*,''

        print*,'Net Ground Gain: ',hydrost%p0
        print*,'Precip Gain: ',tti + tex
        print*,'Transfer Gain: ', xsc+xss+xsb+xc+xs
        print*,'CBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCB'

        if (canb_stop) stop

    endif

    !.....update canopy temperature (K)
    if ((newtc > 400) .or. (newtc < 200)) then
       if (badtc_print) then
          print*,''
          print('(a)'),'   BAD CANOPY TEMP UPDATE IN HYDRO_CANOPY:'
          print('(a,a,i3,a,i4)'), &
              '   DATE: ', trim(month_names(month)), day, ', ', year
          print'(a,i12)','   SECOND: ',sec_tot
          print'(a,2f12.4,1i7,2i4)','   POINT (LON/LAT/REF/PFT): ', &
                 glon, glat, gref, pref
          print'(a,f8.4)',    '   CURRENT LAI: ',lai
          print*,'  ORIG TC/TM: ', tc, tm
          print*,'  NEW TC    : ', newtc
          print*,'  ORIG CAPACC/CAPACS: ',cliqmo,csnowmo
          print*,'  NEW CAPACC/CAPACS: ',cliqm,csnowm
          print*,''
       endif
       if (badtc_stop) stop
    else
       tc = newtc
    endif

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

end subroutine hydro_canopy
