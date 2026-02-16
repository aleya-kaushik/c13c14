!
!================SUBROUTINE PHOSTRESS=========================
!          Calculation of vegetation stress to be 
!                  used in photosynthesis calculations.
!============================================================= 
 !Output
!
!       PAW            PLANT AVAILABLE WATER
!       TAW            TOTAL AVAILABLE WATER (ice + liquid)
!
!       RWSTRESS       SOIL MOISTURE STRESS FACTOR 
!       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS
!           RSTFAC(1) ( F(H-S) : humidity)      :
!                     EQUATION (17,18b), Sellers et al. [1992]
!           RSTFAC(2) ( F(SOIL): soil moisture) :
!                     EQUATION (1-3), Baker et al. [2008]
!           RSTFAC(3) ( F(TEMP): temperature )  :
!                     EQUATION (5b)   , Jarvis [1976]
!           RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
!==============================================================

subroutine phostress(pnum, physcont, &
     ps, etc, tc, &
     eacas, rb, ecmass, &
     td1, td2,    &
     pawfzw, tawfrw, &
     tcmin, co2t, vpd, ect)

    use kinds
    use module_oparams, only: tcbot
    use module_phosib,only: &
        rhfac_astart, rhfac_exp, &
        rhfac_exp_crop, &
        rhfac_nforest, rhfac_tundra
    use module_param, only: phys_param
    use module_pparams, only: &
        amagatwv, &
        molc_h2o, p0_sfc, &
        secs_per_day,   &
        tice,     &
        denh2o, denice
    use module_pftinfo
    use module_sibconst, only: nsoil
    use module_sib, only: &
        co2_type
    use module_time, only: &
        dtisib, dtsib
 
    implicit none

    !...Input Variables
    integer(i4), intent(in) :: pnum
    type(phys_param), intent(in) :: physcont
    real(r8), intent(in) :: ps, etc, tc
    real(r8), intent(in) :: eacas, rb, ecmass
    real(r8), intent(in) :: td1, td2
    real(r8), intent(in) :: pawfzw, tawfrw

    real(r8), intent(in) :: vpd, ect

    real(r8), intent(inout) :: tcmin
    type(co2_type), intent(inout) :: co2t

    !...Humitidy Stress Variables
    real(r8) :: ecmol   ! water vapor flux from CAS to leaf (mol/m2/s)
    real(r8) :: h2oa    ! CAS water vapor mixing ratio (hPa/Pa)
    real(r8) :: h2os    ! leaf surface water vapor mixing ratio (hPa/Pa)
    real(r8) :: h2oi    ! leaf (internal) water vapor mixing ratio (hPa/Pa)
    real(r8) :: h2osrh  ! leaf humidity stress (h2os/h2oi; -)
    real(r8) :: tprcor  ! temperature correction (K)

    !...Soil Moisture Stress Factors
    real(r8) :: lawf
    
    !...Temperature Factors
    real(r8) :: templ     ! low temperature stress (inhibition function; -)
    real(r8) :: temph     ! high temperature stress (inhibition function; -)
    real(r8) :: tempf     ! frost stress (inhibition function; -)

    !...Local Variables
    logical :: iseforest, isnforest, isshrub, iscrop
    real(r8) :: gstemp

    !-------------------------------------------
    !...local variables
    iseforest = (pft_type(pnum) .eq. type_evg)
    isnforest = (pft_group(pnum) .eq. group_ndlfor)
    isshrub = (pft_group(pnum) .eq. group_shrub)
    iscrop = (pft_group(pnum) .eq. group_crop)

    !...leaf humidity stress
    h2oi   = etc / ps
    h2oa   = eacas / ps

    tprcor = tice*ps*p0_sfc
    ecmol = molc_h2o * ecmass * dtisib 
    h2os = h2oa + ecmol / &
            (0.5/rb * amagatwv*tprcor/tc)

    h2os  = min( h2os, h2oi )
    h2os  = max( h2os, 1.0e-7_r8)
    h2osrh = h2os / h2oi

    !...soft landing: add curvature at low 
    !...relative humidities to 
    !...lessen positive feedback
    if (iscrop) then
        h2osrh = h2osrh ** rhfac_exp_crop
    else
       if (h2osrh < rhfac_astart) then
            h2osrh = h2osrh + (rhfac_astart - h2osrh) ** rhfac_exp
       endif
    endif

    !...set stress factor
    co2t%rstfac(1) = h2osrh

    !...set minimum value to relative humidity
    !...stress in larch forests due to hypothesis
    !...needles and deep roots
    !...counterbalances extreme stress
   ! if (isnforest) &
   !      co2t%rstfac(1) = MAX(h2osrh, rhfac_nforest)

    !...set minimum value to relative humidity
    !...stress in grass tundra due to hypothesis
    !...extra moisture from melting permafrost 
    !...counterbalances extreme stress
    if (pnum .eq. pft_c3a) &
         co2t%rstfac(1) = MAX(h2osrh, rhfac_tundra)

    !------------------------
    !....PFT-dependent root zone stress
    if ((iseforest) .or. (iscrop))then
        lawf = pawfzw
    else
        lawf = tawfrw
    endif

    co2t%rstfac(2) = max(0.1, min(1.0, &
         ((1+physcont%wssp)*lawf) / &
          (physcont%wssp+lawf)))

    !if (((iseforest) .or. (iscrop)) .and.  &
    !     (lawf .gt. 0.0)) then
    !      co2t%rstfac(2) = MAX(0.7, co2t%rstfac(2))
    ! endif
     
    !---------------------
    !...temperature stress

    !.....low and high temperature inhibition functions
    templ = 0.98 + EXP(physcont%slti * (physcont%hlti - tc))
    temph = 0.98 + EXP(physcont%shti * (tc - physcont%hhti))

    !...frost stress factor
    if (tc < tcmin) then
       tcmin = tc
    endif
    tcmin = MAX(tcmin, tcbot) !bottom-stop tcmin at -20C
    if (tc > tice) &  !frost recovery at 2C/day
        tcmin = tcmin + ((4.0*dtsib)/86400.0)
    tempf = 1. + EXP(physcont%sfti * (physcont%hfti - tcmin))   

    !...overall temperature scaling factor
    co2t%rstfac(3) = MIN(1.0, 1./(temph*templ*tempf))

    !...soil freeze factors
    co2t%soilfrztg = 1.+exp(-1.5 * &
                (max(270.0_r8, td1)-tice))
    co2t%soilfrztd = 1.+exp(-1.5 * &
                (max(270.0_r8, td2)-tice))
    co2t%soilfrz = max(1./co2t%soilfrztg, 1./co2t%soilfrztd)
    co2t%soilfrz = max( co2t%soilfrz, 0.05_r8)

   !------------------------
   !...combined plant stress
    co2t%rstfac(4) = co2t%rstfac(1) * &
                      co2t%rstfac(2) * co2t%rstfac(3)

   !itb...calculate some stress diagnostics
    gstemp = 1.0_r8/temph
    gstemp = MIN(gstemp,1.0_r8)
    co2t%gs_stress(1) = co2t%gs_stress(1) + gstemp/6.0_r8

    gstemp = 1.0_r8/templ
    gstemp = MIN(gstemp,1.0_r8)
    co2t%gs_stress(2) = co2t%gs_stress(2) + gstemp/6.0_r8

    gstemp = 1.0_r8/tempf
    gstemp = MIN(gstemp,1.0_r8)
    co2t%gs_stress(3) = co2t%gs_stress(3) + gstemp/6.0_r8

    co2t%gs_stress(4) = co2t%gs_stress(4) + co2t%rstfac(2)/6.0_r8

    co2t%gs_stress(5) = co2t%gs_stress(5) + vpd/6.0_r8

    co2t%gs_stress(6) = co2t%gs_stress(6) + co2t%assim/6.0_r8

    co2t%gs_stress(7) = co2t%gs_stress(7) + ect/6.0_r8

end subroutine phostress
