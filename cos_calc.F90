!!!Calculation of uptake of §carbonyl sulfide by plants and soil.
!!!Equations presented here follow:
!
!  Berry J.A. et al., 2013: A coupled model of the global cycles of carbonyl sulfide
!        and CO2: A possible new window on the carbon cycle. Journal of
!        Geophysical Research, Biogeosciences, 118, doi:10.1002/jgrg.20068.
!
!  Ogée, J., J. Sauze, J. Kesselmeier, B. Genty, H. Van Diest, T. Launois, and
!        L.Wingate (2016). A new mechanistic framework to predict OCS fluxes from soils. 
!        Biogeosciences, 13(8):2221-2240. doi: 10.5194/bg-13-2221-2016.
!
!  Meredith, L.K., Ogée, J., Boye, K. et al. Soil exchange rates of COS and CO18O
!        differ with the diversity of microbial communities and their carbonic anhydrase
!        enzymes. ISME J 13, 290-300 (2019). doi: 10.1038/s41396-018-0270-2
!
!  Meredith, L.K.; Boye, K.; Youngerman, C.; Whelan, M.; Ogée, J.; Sauze, J.;
!        Wingate, L. Coupled Biological and Abiotic Mechanisms Driving Carbonyl
!        Sulfide Production in Soils. Soil Syst. 2018, 2, 37.

subroutine cos_calc(gref, pftref, lonsib, latsib, &
           i3, pcosm, tcas, &
           aparkk, co2_assim, rstfac2, &
           poros, zm, woptzm, wsat, &
           rootf3, dz3, www_ice3, www_liq3, &
           soilresp_lay, hrt_soil_temp, tcan, &
           co2t, cost, pref, sscolt)

    use kinds
    use module_oparams, only: &
        k_cos_soil, near_zero
    use module_phosib, only:  &
        co2m, c4, &
        gah2o, gbh2o, gsh2onew, &
        pressure, vmaxts, vmax_unscaled !, qt
    use module_pparams, only: &
        denh2o, &
        p0_sfc, rstar, tice
    use module_sib, only: &
        cos_type, co2_type, sscol_type
    use module_sibconst, only: &
        nsoil, soilogee_switch
    use module_time, only: dtsib
    use module_pftinfo
    use cos_soil_production
    use fcos_solver_Ogee

    implicit none

    !Parameters
    logical, parameter :: &
         cosbad_print = .false., &
         cosbad_stop = .false.

    !Input Variables
    integer(i4), intent(in) :: gref, pftref
    real(r4), intent(in) :: lonsib, latsib
    integer(i4), intent(in) :: i3
    real(r8), intent(in) :: pcosm, tcas
    real(r8), intent(in) :: aparkk, co2_assim, rstfac2
    real(r8), intent(in) :: poros, zm, woptzm, wsat
    real(r8), dimension(i3), intent(in) :: rootf3
    real(r8), dimension(i3), intent(in) :: dz3, www_ice3, www_liq3
    real(r8), dimension(nsoil) :: soilresp_lay
    real(r8), intent(in) :: hrt_soil_temp, tcan
    type(co2_type), intent(in) :: co2t
    type(cos_type), intent(inout) :: cost
    integer(i4), intent(in) :: pref ! PFT
    type(sscol_type), intent(in) :: sscolt

    !Local Variables
    real(r8) :: cosaprev ! previous CAS COS concentration (mol COS/mol air)
    real(r8) :: coscap   ! air capacity for COS exchange (mol air/m2)
    !real(r8) :: cosm     ! reference level COS concentration (mol COS/mol air)

    ! Variables needed for COS soil flux calculation after Berry et al., 2013
    real(r8),dimension(i3) :: &
                cos_grnd_lay, & ! ground COS flux in top 3 layers (mol COS/m2/s)
                resp_lay        ! soil respiration in top 3 layers

    real(r8) :: freeze  ! term to restrict COS uptake by frozen ground
    real(r8) :: moist   ! term to restrict COS uptake by saturated ground
    real(r8) :: wfract    ! fraction of saturation in top 3 soil layers
    real(r8) :: wet_exp   ! wetness scaling of respiration term
    real(r8) :: resp_tot3   ! respiration from top 3 soil layers
    real(r8) :: resp_totall ! respiration from all soil layers
    real(r8) :: root_tot3   ! total root amount in top 3 soil layers

    ! Variables needed for COS soil flux calculation after Ogee et al., 2016
    real(r8) :: cos_P
    real(r8) :: pKw
    real(r8) :: pH
    real(r8) :: T_s
    real(r8) :: Tref
    real(r8) :: kuncat
    real(r8) :: fca
    integer:: pnum
    real(r8) :: vegtype
    logical :: isdesert, iseforest, isdforest, isgrass, isagri

    !Misc Variables
    integer :: j

    !!!!!!!!!!!!!!!!!!!!!!

    !...if no assimilation, reset cos
    if ((gsh2onew .lt. near_zero) .or. &
        (aparkk .lt. near_zero)) then

        cost%cos_assim = dzero
        cost%cos_flux = dzero
        cost%cos_grnd = dzero
        cost%cos_lru = dzero
        cost%cos_lru2 = dzero
        cost%cos_lru3 = dzero
        cost%cos_lru4 = dzero

        cost%coss = cost%coscas
        cost%cosi = cost%coscas
!        cost%cosgm = (dble(1.40E3) * vmaxts * (1.0 + 5.33 * c4) &
!             * aparkk * rstfac2) * (pressure / p0_sfc) * (tcan / tice)
        cost%cosgm = (dble(1.40E3) * vmax_unscaled * (1.0 + 5.33 * c4) &
             * aparkk * rstfac2) * (pressure / p0_sfc) * (tcan / tice)
        cost%cosgt = dzero

    !...calculate COS
    else

        !...set local variables
        cosaprev = cost%coscas
        coscap = cost%cos_casd * pressure/rstar/tcas
        ! COS mixing ratio is either 500 ppt or varies based on TM5 inversions
        ! depending on tm5mr_switch in namel_sibdrv
        cost%cosm = pcosm / pressure   

        !...calculations following paper
        !...the terms 1.40E3 and 5.33 come from empirical regressions
!        cost%cosgm = dble(1.40E3) * vmaxts * (1.0 + 5.33 * c4) * &
!             aparkk * rstfac2 !* 2.1**qt
        cost%cosgm = dble(1.40E3) * vmax_unscaled * (1.0 + 5.33 * c4) * &
             aparkk * rstfac2   !* 2.1**qt

        cost%cosgm = cost%cosgm * (pressure / p0_sfc) * (tcan / tice)

        cost%gsh2onew = gsh2onew

        cost%cosgt = dble(1.0) / ((dble(1.94) / cost%gsh2onew) &
                + (dble(1.56) / gbh2o) + (dble(1.0) / cost%cosgm))
        cost%cos_assim = cost%cosgt * cost%coscas

        if (co2_assim > dzero) then
           cost%cos_lru = (cost%cos_assim/co2_assim)*(co2m/cost%cosm)

           if ((gsh2onew > dzero) .and. (cost%coscas > dzero) &
                .and. (co2t%pco2cas > dzero) .and. (co2t%pco2i .ne. co2t%pco2cas)) then
                 cost%cos_lru2 = (cost%cos_assim/cost%coscas) * (1.6/gsh2onew) &
                       * (1.0 / (1.0 - (co2t%pco2i/co2t%pco2cas)))
           else
                 cost%cos_lru2 = dzero
           endif

           if ((cost%cosgm > dzero) .and. (co2t%pco2cas > dzero)) then
                cost%cos_lru3 = 1./(0.8333 * (1. + &
                          (gsh2onew * 1./1.6 * 1./1.2)/cost%cosgm) &
                          * (1.0 - co2t%pco2i/co2t%pco2cas))
           else
                cost%cos_lru3 = dzero
           endif

           if (co2t%pco2cas > dzero) then
               cost%cos_lru4 = 1./(0.8333 * (1.0 - co2t%pco2i/co2t%pco2cas))
           else
               cost%cos_lru4 = dzero
           endif

        else
           cost%cos_lru = dzero
           cost%cos_lru2 = dzero
           cost%cos_lru3 = dzero
           cost%cos_lru4 = dzero
        endif

        !...ground uptake of COS

        !!!!!!!!!!!!!!!
        !!! Berry et al., 2013 soil flux calculation, coupled to soil respiration
        !!!!!!!!!!!!!!!

        !...calculations limited to top 3 soil layers, but
        !...scaling contributions to total soil respiration
        cos_grnd_lay(:) = dzero
        resp_lay(1) = soilresp_lay(1)
        resp_lay(2) = soilresp_lay(2)
        resp_lay(3) = soilresp_lay(3)
        resp_tot3 = sum(resp_lay)
        resp_totall = sum(soilresp_lay)
        root_tot3 = rootf3(1) + rootf3(2) + rootf3(3)
        do j=1,3  !top 3 soil layers
           if (www_liq3(j) < 1.E-8 .and. www_ice3(j) < 1.E-8) then
              freeze = 1.0
           else
              freeze = www_liq3(j) / (www_liq3(j) + www_ice3(j))
           endif
           wfract = &
             ((www_liq3(j)/ (dz3(j) * poros * denh2o)) &
             * (rootf3(j) / root_tot3))
           moist = MIN(1.0,1.42 - (wfract * 1.42))
           !...adjust moist to utilize soil moisture respiration scaling factor
           if (wfract == 0.0 .AND. zm<0.0) then
               wet_exp = 0.01
           else
               wet_exp = ((wfract**zm - woptzm) / (1.0 - woptzm))**2.0
               wet_exp = MIN(wet_exp,dble(10.0))
           endif
           moist = moist * (0.8 * wsat**wet_exp + 0.2)
           if (resp_tot3 > 0.) then
              cos_grnd_lay(j) = &
                 resp_lay(j)/resp_tot3 * resp_totall * cost%coscas * k_cos_soil * &
                 freeze * moist * hrt_soil_temp
           else
              cos_grnd_lay(j) = 0.
           endif
        enddo
        cost%cos_grnd = sum(cos_grnd_lay)


        !!!!!!!!!!!!!!!
        !!! Ogee et al., 2016 soil flux calculation
        !!! Added by Linda Kooijmans, July 2020.
        !!!!!!!!!!!!!!!
        
        pnum = pft_num(pref)
        vegtype = pft_type(pnum)
        isdesert = (vegtype .eq. type_bare)
        iseforest = (vegtype .eq. type_evg)
        isdforest = (vegtype .eq. type_decid)
        isgrass  = (vegtype .eq. type_grass)
        isagri   = (vegtype .eq. type_crop)
        

        ! COS soil production (cos_P in mol m-3 s-1) based on Meredith et al., 2018 SS and 
        ! COS soil uptake (kuncat in s-1 and fca) based on Meredith et al., 2018 ISME

        pKw = 14.0
        pH  = 4.5        
        T_s = sscolt%td(1)
        Tref = 293
        kuncat = 2.15d-5 * exp(-10450.0*(1./T_s - 1./Tref)) + &
               12.7*10**(-pKw+pH)*exp(-6040.0*(1./T_s - 1./Tref))


        ! if (iswetland) then 
        ! .... other production term can be implemented. 
        if (iseforest) then
           cos_P = COS_enf_production(sscolt%td(3)-273.15)
           kuncat = kuncat !1.2d-5
           fca = 32000
        elseif (isdforest) then
           cos_P = COS_dbf_production(sscolt%td(3)-273.15)
           kuncat = kuncat !1.2d-5
           fca = 32000
        elseif (isgrass) then
           cos_P = COS_gra_production(sscolt%td(3)-273.15)
           kuncat = kuncat !1.3d-5
           fca = 17000
        elseif (isagri) then
           cos_P = COS_cro_production(sscolt%td(3)-273.15)
           kuncat = kuncat/5 !1.5d-5
           fca = 6500
        elseif (isdesert) then
           cos_P = COS_des_production(sscolt%td(3)-273.15)
           kuncat = kuncat !21d-5
           fca = 13000 !762
        endif
  
        call cos_flux_solver_Ogee(&
            cost%coscas, &           ! ambient COS concentration
            sscolt%td(1), &          ! T_s soil temperature (K)
            MAX(0.02, sscolt%eff_poros(2)), &   ! soil porosity (m3 m-3) ice corrected
            !MAX(0.01, sscolt%vol_liq(3))/poros, &
            MAX(0.01, sscolt%vol_liq(1)+sscolt%vol_ice(1)), &
            pressure, &                ! ambient pressure (Pa)
            cos_P*1d-12, &
            kuncat, &
            fca, &
            cost%cos_grnd_Ogee)                     ! soil COS flux (mol m-2 s-1)

        ! choice for the soil model is regulated by the soilogee_switch in
        ! namel_sibdrv
        if (soilogee_switch) then ! if the Ogee switch in namel_sibdrv is turned on,
          cost%cos_soil = cost%cos_grnd_Ogee
          !print*, 'Ogee soil flux'
        else
          cost%cos_soil = cost%cos_grnd
          !print*, 'Berry soil flux'
        endif 


        !print*, '  '
        !print*, 'cosm          ', cosm*1d12
        !print*, 'cos_grnd      ', cost%cos_grnd*1d12
        !print*, 'cos_grnd_Ogee ', cost%cos_grnd_Ogee*1d12
        !print*, 'cos_soil      ', cost%cos_soil*1d12
        !print*, 'poros         ', sscolt%eff_poros(2)
        !print*, 'swc           ', sscolt%vol_liq(1)
        !print*, '              ', MAX(0.02, sscolt%eff_poros(2))-MAX(0.01,sscolt%vol_liq(1))
        
        ! flux units of cos_assim, cos_flux, cos_grnd, cos_grnd_Ogee

        ! cos_assim is later converted to pmol m-2 s-1 when saved in output
        ! files (the other variables are not converted).

        ! Direction of fluxes:
        ! cos_assim: uptake is positive
        ! cos_flux: uptake is negative
        ! cos_grnd: uptake is positive
        ! cos_grnd_Ogee: uptake is positive
        ! cos_soil: uptake is positive
        
        !...prognostic equation for coscas: follow phosib
        cost%coscas = (cosaprev + dtsib/coscap &
                 * ((cost%cosm * gah2o) - cost%cos_soil - cost%cos_assim)) &
                 / (1.0 + (dtsib * gah2o)/coscap)

        !print*, 'coscas        ',   cost%coscas*1d12

        !...make sure COS stays positive
        if (cost%coscas .LT. dzero) then
           cost%badcos = cost%badcos + 1
           IF (cosbad_print) THEN
               print*,'Point/lon/lat: ',gref, lonsib, latsib
               print*,'PFT: ',pftref
               print*,'Bad COS (old/new): ',cosaprev,cost%coscas
               print*,'COS Fluxes (in/assim/grnd): ', &
                    cost%cosm*gah2o,cost%cos_assim,cost%cos_soil
               print*,'Vars input (coscap,gah2o,cosm,dtsib/coscap): ', &
                    coscap,gah2o,cost%cosm,dtsib/coscap

               print*,'Pieces: ',cost%cosm*gah2o,cost%cos_soil,cost%cos_assim, &
                    cost%cosm*gah2o-cost%cos_soil-cost%cos_assim

               print*,'Setting back to previous.'
            ENDIF
            IF (cosbad_stop) THEN
               print*,'Stopping In cos_calc.F90'
               STOP
            ENDIF
           cost%coscas = cosaprev
        endif
        cost%coscasp = cost%coscas * pressure

        !...new cos flux
        cost%cos_assim = cost%cosgt * cost%coscas

        !...back out concentrations at leaf surface and stomate
        cost%coss = cost%coscas - cost%cos_assim*(1.56/gbh2o)
        cost%cosi = cost%coss - cost%cos_assim * (1.94/gsh2onew)

        !...more diagnostics
        cost%cos_flux = gah2o * (cost%coscas - cost%cosm) !mol/m2/s

        !print*, 'cos_flux       ', cost%cos_flux
        
        cost%gsh2onew = gsh2onew

   endif !calculate COS

end subroutine cos_calc
