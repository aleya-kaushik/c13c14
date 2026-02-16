!=====================SUBROUTINE FLUX_UPDATE============================

subroutine flux_update( psy, ros, &
     radtc, radtg, radts, radc3c, radc3g,  &
     poros, paw_lay, cast, &
     fluxt, hydrost, sscolt, press, spdm, dlwbot)

!=======================================================================
! UPDATING OF FLUXES AND ASSOCIATED HYDROLOGICAL PROGNOSTIC VARIABLES.  
!=======================================================================

!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
! FLUXES: 
!    ECI, ECT
!    EGI, EGS, EGSMAX
!    ES
!    HC, HG, HS
!    STORHC, STORHG
!
! HYDROLOGY:
!    ECMASS, EGMASS
!    SNOW_GMASS, SNOW_GDEPTH
!    ROFF
!
! SOIL COLUMN:
!    CAPACC_LIQ, CAPACC_SNOW, CAPACG
!    WWW_ICE, WWW_LIQ
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


use kinds
use module_local
use module_pparams, only: &
    lvap, spec_heat_cp, snofac, &
    denice, denh2o, tice, lfus
use module_sibconst, only:   &
    nsnow, nsoil
use module_sib, only: &
    cas_type, flux_type, &
    hydros_type, sscol_type
use module_time, only: dtsib, dtisib

implicit none

!----------------------------------------------------------------------
!...input variables
real(r8), intent(in) :: psy, ros
real(r8), intent(in) :: radtc, radtg, radts
real(r8), intent(in) :: poros
real(r8), dimension(nsoil), intent(in) :: paw_lay
real(r8), intent(in) :: press
real(r8), intent(in) :: radc3c, radc3g, spdm, dlwbot

type(cas_type), intent(inout) :: cast
type(flux_type), intent(inout) :: fluxt
type(hydros_type), intent(inout) :: hydrost
type(sscol_type), intent(inout) :: sscolt


!...LOCAL VARIABLES
real(r8) :: flux_imb(-nsnow+1:nsoil)
real(r8) :: water_imb(-nsnow+1:nsoil)
real(r8) :: wice0(-nsnow+1:nsoil)
real(r8) :: wliq0(-nsnow+1:nsoil)
real(r8) :: wmass(-nsnow+1:nsoil)
real(r8) :: sflux(-nsnow+1:nsoil)
real(r8) :: sflx1(-nsnow+1:nsoil)
real(r8) :: htstor(-nsnow+1:nsoil)
real(r8) :: tinc(-nsnow+1:nsoil)

real(r8) :: rsnow  ! fraction of snow that is ice
real(r8) :: ecpot
real(r8) :: egpot
real(r8) :: espot
real(r8) :: cogs1
real(r8) :: cogs2
real(r8) :: ectmax(nsoil)  ! layer max transpiration (J m^-2)
real(r8) :: ectmax_tot !column total max transp (J m^-2)
real(r8) :: facks
real(r8) :: cpdpsy
real(r8) :: lvapi
real(r8) :: ecidif
real(r8) :: egidif
real(r8) :: ecit
real(r8) :: egit
real(r8) :: temp
real(r8) :: propor
real(r8) :: heatr
logical :: icemelt

integer(i4) :: i, j
real(r8) :: rnet, t1, t2, t3
!----------------------------------------------------------------------

!...set local variables
imelt(:) = 0
flux_imb(:) = dzero
water_imb(:) = dzero
wice0(:) = dzero
wliq0(:) = dzero
wmass(:) = dzero
sflux(:) = dzero
sflx1(:) = dzero
htstor(:) = dzero
tinc(:) = dzero

lvapi = done / lvap      ! inverse heat of vap, units kg/J
cpdpsy = spec_heat_cp/psy

if (sscolt%nsl < 0) then
   wice0(1) = dzero
   wliq0(1) = dzero
   do i=sscolt%nsl+1,0,-1
      wice0(1) = wice0(1) + sscolt%www_ice(i)
      wliq0(1) = wliq0(1) + sscolt%www_liq(i)
   enddo

  if (wliq0(1) > dzero) then
     rsnow = wice0(1) / (wliq0(1) + hydrost%capacg)
  else
     rsnow = done
  endif
  wice0(1) = dzero
  wliq0(1) = dzero
else
  rsnow = dzero
endif
facks = 1. + rsnow * (snofac-1.)

!------------------------------------------------------
!...GROUND EVAPORATION AND INTERCEPTION FLUXES (J/m2)
!------------------------------------------------------
if (sscolt%nsl == 0 ) then  !NO SNOW
   !...potential evaporative gradient (ground to CAS)
   egpot = (etg + getg*dtd(1)) - cast%eacas
   espot = dzero

   cogs1 =  gegs * hydrost%rhsoil  
   cogs2 =  gegs

   !...Make sure you don't evap more than you have (units are J/m2).
   !.....Limiting evaporation to 1/2 the water in top soil layer. 
   fluxt%egsmax = 0.5_r8 * (sscolt%www_liq(1) + sscolt%www_ice(1)) * lvap
   fluxt%egs = min(fluxt%egsmax, &
               ((etg + getg * dtd(1)) * cogs1 &
                - cast%eacas * cogs2 ) * ros * cpdpsy * dtsib)

    !...interception flux
    fluxt%egi = egpot * gegi * ros * cpdpsy * dtsib

    !...snow evaporation
    fluxt%es = dzero

else  !SNOW CASE

    !...liquid evaporation
    fluxt%egs = dzero

    !...liquid surface interception
    fluxt%egi = dzero

    !...snow surface evaporation
    !.....potential gradient (snow to CAS)
    espot = (ets + gets * dtd(sscolt%nsl+1)) - cast%eacas
    fluxt%es = espot * ros * cpdpsy / fluxt%rd / snofac * dtsib
    !...can't evaporate more snow than you have...
    fluxt%es = min( fluxt%es, lvap / snofac &
                 * (sscolt%www_ice(sscolt%nsl+1) &
                  + sscolt%www_liq(sscolt%nsl+1)))

    !...now subtract it out
    sscolt%www_ice(sscolt%nsl+1) = sscolt%www_ice(sscolt%nsl+1) &
             - fluxt%es * snofac / lvap

endif

!.....adjust ground interception by surface store
egidif = MAX(dzero, &
            fluxt%egi - (hydrost%snow_gmass + hydrost%capacg) * lvap)
egit = MIN(fluxt%egi, &
            (hydrost%snow_gmass + hydrost%capacg) * lvap + fluxt%es)
fluxt%egi = egit

!.....reduce ground surface store (ponding) by ground intercepted flux
hydrost%capacg = hydrost%capacg - fluxt%egi * lvapi
hydrost%egmass = (fluxt%egi + fluxt%egs) * facks * lvapi

!--------------------------------------------------------
!...CANOPY INTERCEPTION and TRANSPIRATION FLUXES (J/m2)
!-------------------------------------------------------

!...potential gradient in hPa (canopy to CAS)
ecpot = (etc + getc*dtc) - cast%eacas

!...canopy interception
fluxt%eci = ecpot * geci * ros * cpdpsy * dtsib

!.....adjust canopy interception by canopy storage
ecidif = MAX(dzero, &
         fluxt%eci - hydrost%capacc_snow * lvap - hydrost%capacc_liq * lvap)
ecit = MIN(fluxt%eci, hydrost%capacc_snow * lvap + hydrost%capacc_liq * lvap)
fluxt%eci = ecit

!.....update canopy storage
rsnow = hydrost%capacc_snow/(hydrost%capacc_liq + &
        hydrost%capacc_snow + 1.0E-10)
facks = 1. + rsnow * (snofac-1.)
hydrost%capacc_liq = hydrost%capacc_liq - (1.-rsnow)*ecit*facks*lvapi
hydrost%capacc_snow = hydrost%capacc_snow - rsnow*ecit*facks*lvapi
if (hydrost%capacc_liq .lt. dzero) then
   hydrost%capacc_snow = hydrost%capacc_snow + hydrost%capacc_liq
   hydrost%capacc_liq = dzero
endif
if (hydrost%capacc_snow .lt. dzero) then
   hydrost%capacc_liq = hydrost%capacc_liq + hydrost%capacc_snow
   hydrost%capacc_snow = dzero
endif
if (abs(hydrost%capacc_liq) .lt. 1.E-10) hydrost%capacc_liq = dzero
if (abs(hydrost%capacc_snow) .lt. 1.E-10) hydrost%capacc_snow = dzero

if ((hydrost%capacc_liq .lt. dzero) .or. &
    (hydrost%capacc_snow .lt. dzero)) then
    print*,'Error in ECI calculation in hydro_update.F90.'
    print*,'capacc_liq: ',hydrost%capacc_liq
    print*,'capacc_snow: ',hydrost%capacc_snow
    stop
endif

!...canopy transpiration
fluxt%ect = ecpot * gect * ros * cpdpsy * dtsib

!.....Make sure you don't transpire more
!.....water than available
ectmax(:)  = dzero
ectmax_tot = dzero
do i=1,nsoil
    ectmax(i)  = paw_lay(i) * sscolt%dz(i) * lvap
    ectmax_tot = ectmax_tot + ectmax(i)
enddo
fluxt%ect = min(ectmax_tot, fluxt%ect)

!...canopy evapotranspiration
hydrost%ecmass = (fluxt%eci + fluxt%ect) * facks * lvapi


!---------------------------------------------------------------------
!     CALCULATION OF SENSIBLE HEAT FLUXES FOR THE END OF THE TIMESTEP.
!        SEE FIGURE (2) OF SE-86.  NOTE THAT INTERCEPTION LOSS EXCESS
!        ENERGIES (ECIDIF, EGIDIF) ARE ADDED.
!
!      HC          (HC)    : EQUATION (63) , SE-86
!      HG          (HGS)   : EQUATION (65) , SE-86
!----------------------------------------------------------------------

 !...one-sided leaf for now...
fluxt%hc = (cast%tc - cast%tcas) / fluxt%rb &
            * ros * spec_heat_cp * dtsib + ecidif
fluxt%storhc = cast%hcapc * dtisib * dtc

!..keep in mind that currently we are not worrying about 
!...partial snow cover issues. with that in mind, there will
!...only exist ground H when the is NO snow, snow H when
!...there IS snow.

!...snow fluxes
if(sscolt%nsl < 0 ) then  !SNOW case
   !...remember that td was updated in addinc
   fluxt%hs = (sscolt%td(sscolt%nsl+1) - cast%tcas) / &
              fluxt%rd * ros * spec_heat_cp * dtsib + egidif
   fluxt%hg = dzero
else
   fluxt%hg = (sscolt%td(1) - cast%tcas) / &
        fluxt%rd * ros * spec_heat_cp * dtsib + egidif
   fluxt%hs = dzero
endif

!----------------------------------------------------------------------
!     CALCULATION OF STORAGE HEAT FLUXES
!---------------------------------------------------------------------- 

if (sscolt%nsl == 0) then
    fluxt%storhg = dtd(1) * sscolt%shcap(1) * dtisib + &
        sscolt%slamda(1) * (sscolt%td(1) - sscolt%td(2))
else
    fluxt%storhg = dtd(sscolt%nsl+1) * sscolt%shcap(sscolt%nsl+1) * dtisib + &
           sscolt%slamda(sscolt%nsl+1) * &
           (sscolt%td(sscolt%nsl+1) - sscolt%td(sscolt%nsl+2))
endif

!----------------------------------------------------------------------
!   CALCULATION OF PHASE CHANGES WITHIN SNOW/SOIL LAYERS
!   This code based on routine CLM_MELTFREEZE from the Common
!   Land Model (CLM)  (Dai et al, submitted)
!   CLM web info:  http://clm.gsfc.nasa.gov
!----------------------------------------------------------------------

!...initialize some values
do j = sscolt%nsl+1,nsoil
   imelt(j) = 0
   flux_imb(j) = 0.
   water_imb(j) = 0.

   wice0(j) = sscolt%www_ice(j)
   wliq0(j) = sscolt%www_liq(j)
   wmass(j) = sscolt%www_ice(j) + sscolt%www_liq(j)
enddo

!...calculate some fluxes. 
do j=sscolt%nsl+1,nsoil-1
   sflux(j) = sscolt%slamda(j)*(td_old(j+1)-td_old(j))
   sflx1(j) = sscolt%slamda(j)*(sscolt%td(j+1) - sscolt%td(j))
enddo
sflux(nsoil) = dzero
sflx1(nsoil) = dzero

!...melting check
do j = sscolt%nsl+1,nsoil

   !...if air temperature greater than freezing,
   !......do not freeze more
   icemelt = (sscolt%td(j) .ge. tice)

   !...if ice exists in warm conditions, melt it
   if (sscolt%www_ice(j) > 0.0 .and. icemelt) then
      imelt(j) = 1
      tinc(j) = sscolt%td(j) - tice
      sscolt%td(j) = tice
   endif

  !...if water exists in cold conditions, freeze it
   if (sscolt%www_liq(j) > 0.0 .and. (.not. icemelt)) then
      imelt(j)  = 2
      tinc(j) = sscolt%td(j) - tice
      sscolt%td(j) = tice

   endif

enddo

!...check for existence of snow, less depth than 0.01 m
if(sscolt%nsl == 0  .and. hydrost%snow_gmass > 0.0) then
   if (sscolt%td(1) > tice) then
       imelt(1) = 1
       tinc(1) = sscolt%td(1) - tice
       sscolt%td(1)    = tice
    endif
endif

!...change in storage
do j = sscolt%nsl+1, nsoil
    htstor(j) = sscolt%shcap(j)*dtisib*   &
       (sscolt%td(j) - td_old(j))
enddo

!...calculate energy surplus or loss if there was melting or freezing
do j = sscolt%nsl+1,nsoil
    if(imelt(j) > 0 ) then    ! did melting/freezing occur?
       flux_imb(j) = tinc(j) * sscolt%shcap(j) ! (J/m^2)
    endif
enddo

!...the CLM boys say these have been checked carefully-they are the 
!...result of computed error in the tridiagonal matrix solver.
do j=sscolt%nsl+1,nsoil
    if (imelt(j) == 1 .and. flux_imb(j) < 0.0 ) then
        imelt(j)    = 0
        flux_imb(j) = dzero
    endif

    if (imelt(j) == 2 .and. flux_imb(j) > 0.0 ) then
        imelt(j)    = 0
        flux_imb(j) = dzero
     endif
enddo 

do j=sscolt%nsl+1,nsoil

   !...melting or freezing occurring, some error present
   if (imelt(j) > 0 .and. abs(flux_imb(j))>0.0) then

      !...water_imb > 0 ==> melting ice
      !...water_imb < 0 ==> freezing water
       water_imb(j) = flux_imb(j) / lfus  !(kg of water)

       !...snow exists, but less than critical depth. CLM boys say this
       !...needs work.

       !...LEVEL 1 - top soil layer, small amt of snow present, 
       !...no snow layers
       if (j == 1) then  
          if (sscolt%nsl == 0 .and. hydrost%snow_gdepth >  &
              0.01_r8 .and. water_imb(j) > 0.0_r8) then 
               temp = hydrost%snow_gmass
               hydrost%snow_gmass = max(dzero,temp-water_imb(j))
               propor = hydrost%snow_gmass/temp
               hydrost%snow_gdepth = hydrost%snow_gdepth * propor

               heatr = flux_imb(j)*dtisib -       &
                        lfus*(temp-hydrost%snow_gmass)*dtisib

               if (heatr > 0.0) then
                  water_imb(j) = heatr*dtsib/lfus 
                  flux_imb(j) = heatr
               else
                   water_imb(j) = dzero
                   flux_imb(j)  = dzero
                endif

            endif  
        endif ! j==1

        heatr = dzero
        if (water_imb(j) > 0.0) then
            sscolt%www_ice(j) = max(0.0_r8,wice0(j)-water_imb(j))
            heatr = flux_imb(j)*dtisib  &
                    - lfus*(wice0(j)-sscolt%www_ice(j))*dtisib
        elseif(water_imb(j) < 0.0_r8) then
            sscolt%www_ice(j) = min(wmass(j),wice0(j)-water_imb(j))
            heatr = flux_imb(j)*dtisib  &
                  - lfus*(wice0(j)-sscolt%www_ice(j))*dtisib
        endif

        sscolt%www_liq(j) = max(dzero,wmass(j)-sscolt%www_ice(j))

        if (abs(heatr) > 0.0) then
            sscolt%td(j) = sscolt%td(j) + dtsib*heatr/sscolt%shcap(j)

            if (sscolt%www_liq(j) * sscolt%www_ice(j) > 0.0) then
                sscolt%td(j) = tice
            endif
        endif 

   endif   ! imelt \= 0 and flux_imb \= 0
enddo 

!----------------------------------------------------------------------
!      UPDATE ENERGY, HYDROLOGY, AND RADIATION
!----------------------------------------------------------------------

!...calculate change in energy and water of CAS
   cas_e_storage = dta * cast%hcapcas * dtisib
   cas_w_storage = dea * cast%vcapcas * dtisib

!...update total soil water
   hydrost%www_tot = dzero
   do i=1,nsoil
      hydrost%www_tot = hydrost%www_tot + &
         sscolt%www_liq(i) + sscolt%www_ice(i)
   enddo

!...update effective porosity
   do i = sscolt%nsl+1,nsoil
      sscolt%vol_ice(i) = min(poros, &
           sscolt%www_ice(i) / (sscolt%dz(i) * denice))
      sscolt%eff_poros(i) = poros - sscolt%vol_ice(i)
      sscolt%vol_liq(i) = min(sscolt%eff_poros(i), &
           sscolt%www_liq(i) / (sscolt%dz(i) * denh2o))

      !.....if no effective porosity,
      !.....lose liquid to runoff
      if (sscolt%vol_liq(i) == 0.0 .and. sscolt%www_liq(i) > 0.0 ) then
           hydrost%roff = hydrost%roff + sscolt%www_liq(i)
           sscolt%www_liq(i) = dzero
      endif
   enddo

!...update radiation
   if (sscolt%nsl == 0) then
      radttc = radtc - lcdtc*dtc &
           - lcdtg*dtd(1)
      
      radttg = radtg - lgdtc*dtc &
           - lgdtg*dtd(1)
      
      radtts = 0.
   else
      radttc = radtc - lcdtc*dtc &
           - lcdts*dtd(sscolt%nsl+1)
      
      radttg = 0.
      
      radtts = radts - lsdts*dtd(sscolt%nsl+1) &
           - lsdtc*dtc
   endif

!----------------------------------------------------------------------
!      CALCULATE VPD, Potential ET (from Ian)
!----------------------------------------------------------------------

!...VPD
ppl(1) = press *100.0
ttl(1) = cast%tcas
call dtess_eau(1,ppl,ttl,esst,dtesst)
cast%vpd = (esst(1) - cast%eacas * 100.0)/100.0
cast%vpd = MAX(cast%vpd,0.0_r8)

!...PET: ASCE Penman-Monteith calculation
rnet = (radc3c * radc3g - dlwbot - fluxt%storhg) * (600.0/1.0E6)
t1 = (dtesst(1)/100.0) * rnet
t2 = ros * spec_heat_cp * (cast%vpd/fluxt%ra)
t3 = (dtesst(1)/100.0) + psy*(1.0_r8 + (fluxt%rc+fluxt%rb)/fluxt%ra)
fluxt%et0 = (((t1 + t2) / t3)/lvap)*3600.0_r8
fluxt%et0a = fluxt%et0a + fluxt%et0/6.0_r8

!...PET: ASCE Standardized calculation (more empirical)
t1 = 0.408_r8 * (dtesst(1)/100.0) * rnet
t2 = psy * (66.0_r8 / cast%tcas) * spdm * cast%vpd
t3 = (dtesst(1)/100.0) + psy * (1.0_r8 + (0.25 * spdm))
fluxt%et1 = (((t1 + t2)/t3)/lvap)*3600.0_r8
fluxt%et1a = fluxt%et1a + fluxt%et1/6.0_r8

end subroutine flux_update
