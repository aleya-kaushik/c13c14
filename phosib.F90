!==================SUBROUTINE PHOSIB===================================
!
!     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED   
!     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
!     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
!     BASE REFERENCE IS SE-92A
!
!                       CONVERSIONS
!                      -------------
!
!      1 MOL H2O           = 0.018 KG
!      1 MOL CO2           = 0.044 KG
!      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
!      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
!      GS  (CO2)           = GS (H2O) * 1./1.6
!      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO  
!               44.6 is the number of moles of air per cubic meter
!      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!      MM  (MOLAIR/MOLH2O) = 1.611
!
!                         OUTPUT
!                      -------------
!
!      ASSIM               = CANOPY GROSS ASSIMILATION RATE (MOL/M2/S)
!      ASSIMMFAC(4)        = CANOPY ASSIMILATION STRESS FACTORS
!          ASSIMFAC(1) :  RuBP (Rubisco)-limited rate 
!          ASSIMFAC(2) :  Light-limited rate
!          ASSIMFAC(3) :  Utlilzation/Export-limited rate
!          ASSIMFAC(4) :  Combined 1-3
!
!      PCO2CAS             = CAS CO2 CONCENTRATION (PA)
!      PCO2C               = CHLOROPLAST CO2 CONCENTRATION (PA)
!      PCO2I               = INTERNAL CO2 CONCENTRATION (PA)
!      PCO2S               = LEAF SURFACE CO2 CONCENTRATION (PA)
!
!      RST                 = CANOPY RESISTANCE (S/M)
!
!-----------------------------------------------------------------------
!
!                         NOTES
!                        -------
!
!...A co2 capacity is calculated
!      - represents mass of air under canopy top
!      - canopy top is MIN(4.0, z2)
!
!...All the carbon fluxes are expressed as mol C/m2/s
!
!----------------------------------------------------------------------

subroutine phosib(  &
    physcont, gprogt, gdiagt,  &
    tc, tcas, ra, rb, &
    resp_leaf, resp_auto, resp_het, &
    vegt, co2t)

    use kinds
    use module_param, only: &
        phys_param
    use module_phosib
    use module_pparams, only: &
        amagatwv, p0_sfc, tice,  &
        rstar, par_conv, po2m, &
        secs_per_day
    use module_sib, only: &
        gprog_type, gdiag_type, &
        veg_type, co2_type
    use module_sibconst, only: &
        nsoil, spinup
    use module_time!, only: &
        !dtsib, wt_clim, wt_daily

    implicit none

    !----------------------------------------------------------------------
    !...input variables
    type(phys_param), intent(in) :: physcont
    type(gprog_type), intent(inout) :: gprogt
    type(gdiag_type), intent(in) :: gdiagt
    real(r8), intent(in) :: tc, tcas, ra, rb
    real(r8), intent(in) :: resp_leaf, resp_auto, resp_het

    type(veg_type), intent(inout)  :: vegt
    type(co2_type), intent(inout)  :: co2t
    !type(pooll_type), intent(inout) :: poollt

    !...local variables
    integer :: icconv, ic, ic1, ic2 ! iteration loops
    !integer :: ictmp, icconv_final ! iteration loops
    logical :: cbad, cbad2
    real(r8) :: pco2cas_tmp, co2cas_tmp
    !----------------------------------
    !     SET LOCAL VARIABLES
    !----------------------------------
    !...VEGETATION VARIABLES
    c3 = dble(1.0-physcont%c4flag)
    c4 = dble(physcont%c4flag)
    atheta = dble(physcont%atheta)
    btheta = dble(physcont%btheta)
 
    !...ATMOSPHERIC VARIABLES
    pressure = dble(gprogt%ps) * 100.0D0
    !co2t%press_phosib = pressure
    !co2t%press_phosibps = dble(gprogt%ps)
    tprcor = tice*pressure/p0_sfc

    !...CO2 VARIABLES
    gprogt%pco2m = dble(gprogt%co2m)*pressure
    pco2cas = co2t%pco2cas
    pco2c = co2t%pco2c
    pco2i = co2t%pco2i
    pco2s = co2t%pco2s
    co2cas = pco2cas / pressure
    !co2m = gprogt%pco2m  / pressure 
    !if (spinup) then 
      !if spinning up with varying co2 off
      !very first timestep, read in from restart_read
      !in restart_read, bco2m is *p0_sfc 
     ! co2m = gprogt%pco2m / p0_sfc
    !else 
      !during the actual run with set_co2 always on,
      !at every subsequent timestep pco2m is calculated
      !within co2_module and updated by *pressure
    co2m = dble(gprogt%co2m)
    !print *,'co2m top of phosib', co2m
    !endif

    !print *,'pco2m being read into phosib :',gprogt%pco2m    
    !print *,'co2m being read into phosib :',co2m

    co2cap = co2t%casd * pressure/rstar/tcas
    resp_cas = resp_auto + resp_het

    !...DAMPING FACTORS
    pdamp = exp (-1.0 * zln2*(dtsib*dmin)/(dttin*ghalf))
    qdamp = 1.0 - pdamp

    !-----------------------------------------------------------------------
    !     CALCULATION OF RADIATION VALUES AND CANOPY PAR USE PARAMETER.
    !     APARKK      (PI)     : EQUATION (31) , SE-92A
    !-----------------------------------------------------------------------
    scatc = physcont%tran(1,1) +  physcont%ref(1,1)
    vegt%park = sqrt(1.-scatc) * vegt%gmudmu
    co2t%aparkk = vegt%fpar / vegt%park
    co2t%par = par_conv * (gdiagt%radvbc+gdiagt%radvdc)

    if (gdiagt%radvbc < minrad) then
        pfd=0.0
        toa_pfd = 0.0
    else
        pfd = par_conv * vegt%gmudmu * &
            ( gdiagt%radvbc + gdiagt%radvdc )
        toa_pfd = toatosfc * par_conv * vegt%gmudmu * &
            ( gdiagt%toa_radvbc + gdiagt%toa_radvdc)
    endif
    co2t%apar = pfd * physcont%effcon * (1. - scatc)
    toa_apar = toa_pfd * physcont%effcon * (1. - scatc)
    co2t%nspar = pfd * (1. - scatc)

    !--------------------------------------
    !  VMAX STRESS EFFECTS
    !  USES QT RELATIONSHIP : TABLE 2 IN SE-92A
    !---------------------------------------
    qt = 0.1 * ( tc - vmax_tref )

    !...temperature scaled Vmax
    vmaxts = vegt%vmax * vmax_q10**qt
    vmax_unscaled = vegt%vmax

    !...stress scaled Vmax (temperature and moisture)
    co2t%vmaxss = vmaxts*co2t%rstfac(2)*co2t%rstfac(3)
    !if (pawfrw .le. 1.E-6) co2t%vmaxss=dzero
  
    !-----------------------------------------------------------------------
    !     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
    !     COMPENSATION POINT       
    !
    !      ZKC          (KC)     : TABLE (2)     , SE-92A
    !      ZKO          (KO)     : TABLE (2)     , SE-92A
    !      SPFY         (S)      : TABLE (2)     , SE-92A
    !      GAMMA        (GAMMA-*): TABLE (2)     , SE-92A
    !-----------------------------------------------------------------------
     zkc     = 30. * 2.1**qt
     zko     = 30000. * 1.2**qt
     spfy    = 2600. * 0.75**qt
     co2t%gamma   = 0.5 * po2m/spfy * c3

    !-----------------------------------------
    ! UTILIZATION/EXPORT LIMITED ASSIMILATION
    !  OMSS (OMEGA-S): EQUATION (13) , SE-92A
    !-----------------------------------------
     rrkk = zkc*( 1. + po2m/zko ) * c3 &
              + vegt%vmax/5.* ( oms_q10**qt) * c4
     omss = ( vegt%vmax / 2.0 ) * ( oms_q10**qt ) &
             /co2t%rstfac(2) * c3 &
             + rrkk * co2t%rstfac(2) * c4

    !-----------------------------------------------------------------------
    !  OFFSET FOR BALL-BERRY RELATIONSHIP
    !   BINTC (B*ZLT)  : EQUATION (35) , SE-92A
    !-----------------------------------------------------------------------
    bintc = physcont%binter * vegt%lai * &
            co2t%rstfac(2) * co2t%soilfrz

    gminc = physcont%gmin * vegt%lai * &
            co2t%rstfac(2) * co2t%soilfrz

    !---------------------------------------
    ! CONVERT RESISTANCES TO CONDUCTANCES
    !---------------------------------------
    !...convert resistances to conductances (mol/m2/s)
    gsh2o  = 1.0/co2t%rst * amagatwv*tprcor/tc
    gbh2o  = 0.5/rb * amagatwv*tprcor/tc
    gah2o  = 1.0/ra * amagatwv*tprcor/gprogt%tm

    !...canopy biophysical-maximum conductance
    gxco2 = physcont%gmeso * vegt%vmax * co2t%aparkk * co2t%rstfac(2)

    !-------------------------------------------------------------------
    !                HERE IS THE ITERATION LOOP.
    !
    !...We iterate on PCO2C-sortin makes a 'first guess' at
    !...then orders PCO2C/Error pairs on increasing error size,
    !...then uses a combination of linear/quadratic fit to obtain 
    !...the 'next best guess' as iteration count increases.
    !...CYCALC uses that value of PCO2C to get the next value 
    !...of ASSIM. CO2CAS and CO2S follow.
    !-------------------------------------------------------------------

    !...first guess value of chloroplast CO2 is midway between
    !...compensation point and maximum
    range  = gprogt%pco2m * ( 1. - 1.6/physcont%gradm ) - co2t%gamma
    co2t%resp_casn = resp_cas * vegt%vmax * co2t%rstfac(2) * 2.0**qt
   
    eyy(:) = 0.
    pco2y(:) = 0.
    icconv=1
    co2cas_tmp=0.
    pco2cas_tmp=0.
    !icconv_final=1
    !ictmp = 1

    do ic = 1, numic
        call phosort( numic, ic, co2t%gamma, range, eyy, pco2y)
        call phocycalc( toa_apar, co2t%apar, co2t%aparkk,       &
                     vegt%vmax, co2t%vmaxss,     &
                     atheta, btheta, co2t%gamma,   &
                     rrkk, omss, c3, c4, pco2y(ic), &
                     assimc(ic), assime(ic), assims(ic), &
                     assimy(ic))

        !...prognose the new CAS CO2 according to flux divergence
        !...doing this in mol C / mol air (same as Pa C / Pa air)
        !assimn = assimy(ic) - resp_leaf
        !assimn = assimy(ic)

!        co2cas = (co2cas + (dtsib/co2cap) *  &
!                (co2t%resp_casn - assimn + co2m*gah2o) )  &
!                / (1+dtsib*gah2o/ co2cap )

        !co2cas = (co2cas + (dtsib/co2cap) *  & 
        !        (resp_cas - assimn + co2m*gah2o) )  &
        !        / (1+dtsib*gah2o/ co2cap ) 

        !co2cas_tmp = (co2cas + (dtsib/co2cap) *  & 
        !        (resp_cas - assimy(ic) + co2m*gah2o) )  &
        !        / (1+dtsib*gah2o/ co2cap ) 

        !.. test what happens if you turn above recycling off
        co2cas_tmp = co2m
        !co2cas = (co2m*gah2o)   &
        !        / (1+dtsib*gah2o/ co2cap ) 

        !co2cas = (co2cas + co2m*gah2o)   &
        !        / (1+dtsib*gah2o/ co2cap )         

        pco2cas_tmp = co2cas_tmp * pressure
        assimn = assimy(ic) - resp_leaf
        
        !...intermediate leaf surface CO2
        pco2s = pco2cas_tmp - (1.4/gbh2o * assimn * pressure)

        !...intermediate leaf internal CO2
        pco2i = pco2s - (1.6/gsh2o * assimn * pressure)

        !...intermediate leaf chloroplast CO2-this is what we iterate on 
        pco2c = pco2i - c3 * assimn * pressure * 1.0/gxco2
        eyy(ic) = pco2y(ic) - pco2c

        !... OLD phosib
        !if (ic>=2) then
        !    ic1 = ic-1
        !    if (abs(eyy(ic1))>=0.1)then
        !        icconv = ic
        !    else
        !        eyy(ic) = eyy(ic1)
        !        pco2y(ic) = pco2y(ic1)
        !    endif
        !endif

        !... NEW phosib
        !... check for convergence in the previous timestep
        ic1 = ic -1
        if (ic>=2) then
           if (abs(eyy(ic1)) .lt. 0.1) then
             icconv = ic1
             eyy(ic) = eyy(ic1)
             pco2y(ic) = pco2y(ic1)
           endif
           if (abs(eyy(ic1)) .ge. 0.1) then
             icconv = ic
             eyy(icconv) = eyy(ic)
             pco2y(icconv) = pco2y(ic)
           endif
        endif

    !if ( (ic .eq. numic) .and. (abs(eyy(ic)) .ge. 0.1) ) then
    !  eyy(:) = 0.
    !  pco2y(:) = 0.
    !  icconv=1
    !  do ic2=1,numic
    !    call phosort( numic, ic2, co2t%gamma, range, eyy, pco2y)
    !    call phocycalc( toa_apar, co2t%apar, co2t%aparkk,       &
    !              vegt%vmax, co2t%vmaxss,     &
    !              atheta, btheta, co2t%gamma,   &
    !              rrkk, omss, c3, c4, pco2y(ic2), &
    !              assimc(ic2), assime(ic2), assims(ic2), &
    !              assimy(ic2))

    !    !co2cas = co2m
    !    co2cas = co2m + (resp_cas - assimy(ic)) &
    !             /(co2cap/dtsib + gah2o)
    !    pco2cas = co2cas * pressure
    !    assimn = assimy(ic2) - resp_leaf

    !    !...intermediate leaf surface CO2
    !    pco2s = pco2cas - (1.4/gbh2o * assimn * pressure)

    !    !...intermediate leaf internal CO2
    !    pco2i = pco2s - (1.6/gsh2o * assimn * pressure)

    !    !...intermediate leaf chloroplast CO2-this is what we iterate on 
    !    pco2c = pco2i - c3 * assimn * pressure * 1.0/gxco2
    !    eyy(ic2) = pco2y(ic2) - pco2c

    !    ic1 = ic2 -1
    !    if ((ic2>=2) .and. (abs(eyy(ic1)) .lt. 0.1)) then
    !       icconv = ic1
    !       eyy(ic2) = eyy(ic1)
    !       pco2y(ic2) = pco2y(ic1)
    !    endif
    !    if ((ic2>=2) .and. (abs(eyy(ic1)) .ge. 0.1)) then
    !       icconv = ic2
    !       eyy(icconv) = eyy(ic2)
    !       pco2y(icconv) = pco2y(ic2)
    !    endif
    !  enddo !sub-iteration loop if not converged the first time
    !endif !if not converged the first time

!    enddo !main iteration loop

!cbad = 1.e6*dble(pco2c / pressure) > 1000
!cbad2 = 1.e6*dble(pco2cas_tmp / pressure) > 1000
!if ((cbad) .or. (cbad2)) then
! print*,' '
! print*,'co2c bad from phosib: ',1.e6*dble(pco2c / pressure)
! print*,'co2cas bad from phosib: ',1.e6*dble(pco2cas_tmp / pressure)
! print*,'eyy_ic :',eyy(ic)
! print*,'eyy_icconv :',eyy(icconv)
! print*,'year: ',year
! print*,'month: ',month
! print*,'doy: ',doy
! print*,'day: ',day
! print*,'hour: ',hour
! print*,'sec: ',sec_day
! print*,' '
! print*,'toa_apar: ',toa_apar
! print*,'co2t%apar: ',co2t%apar
! print*,'co2t%aparkk: ',co2t%aparkk
! print*,'vegt%vmax: ',vegt%vmax
! print*,'co2t%vmaxss: ',co2t%vmaxss
! print*,'rrkk: ',rrkk
! print*,'omss: ',omss
! print*,' '
!
!endif


    !... check if last set converged, if not then
    !... set co2cas to co2m and recalculate again
    !if ( (ic .eq. numic) .and. (abs(eyy(ic1)) .ge. 0.1) ) then
    !  eyy(:) = 0.
    !  pco2y(:) = 0.
    !  icconv=1
    !  do ic2=1,numic
    !    call phosort( numic, ic2, co2t%gamma, range, eyy, pco2y)
    !    call phocycalc( toa_apar, co2t%apar, co2t%aparkk,       &
    !              vegt%vmax, co2t%vmaxss,     &
    !              atheta, btheta, co2t%gamma,   &
    !              rrkk, omss, c3, c4, pco2y(ic2), &
    !              assimc(ic2), assime(ic2), assims(ic2), &
    !              assimy(ic2))

    !    co2cas = co2m
    !    pco2cas = co2cas * pressure
    !    assimn = assimy(ic2) - resp_leaf

    !    !...intermediate leaf surface CO2
    !    pco2s = pco2cas - (1.4/gbh2o * assimn * pressure)

    !    !...intermediate leaf internal CO2
    !    pco2i = pco2s - (1.6/gsh2o * assimn * pressure)

    !    !...intermediate leaf chloroplast CO2-this is what we iterate on 
    !    pco2c = pco2i - c3 * assimn * pressure * 1.0/gxco2
    !    eyy(ic2) = pco2y(ic2) - pco2c

    !    ic1 = ic2 -1
    !    if ((ic2>=2) .and. (abs(eyy(ic1)) .lt. 0.1)) then
    !       icconv = ic1
    !       eyy(ic2) = eyy(ic1)
    !       pco2y(ic2) = pco2y(ic1)
    !    endif
    !  enddo !sub-iteration loop if not converged the first time
    !endif !if not converged the first time

  enddo !main iteration loop

    co2t%icconv_phosib = icconv
    co2t%eyy_phosib = eyy(icconv)

    !...save variables
    assim_omc = assimc(icconv)
    assim_ome = assime(icconv)
    assim_oms = assims(icconv)
!    assimn = assimy(icconv) - resp_auto
    assimn = assimy(icconv) - resp_leaf
    
    co2cas = (co2cas + (dtsib/co2cap) *  &
                (resp_cas - assimy(icconv) + co2m*gah2o) )  &
                / (1+dtsib*gah2o/ co2cap )

    co2cas=co2m

    co2t%pco2c = pco2y(icconv)
    co2t%pco2i = co2t%pco2c + (c3 * assimn/gxco2*pressure)
    co2t%pco2s = co2t%pco2i + (1.6 * assimn/gsh2o*pressure)
    co2t%assim = MAX(0.0,assimy(icconv))
    co2t%assimd = co2t%assimd*(1. - wt_daily) &
                  + co2t%assim*wt_daily
    co2t%clim_assim = (1.-wt_clim)*co2t%clim_assim &
                  + wt_clim*co2t%assim


!cbad = 1.e6*dble(co2t%pco2c / pressure) > 2000
!if (cbad) then
! print*,'co2c bad from phosib: ',1.e6*dble(co2t%pco2c / pressure)
! print*,'eyy_ic :',eyy(ic)
! print*,'eyy_icconv :',eyy(icconv)
! print*,'year: ',year
! print*,'month: ',month
! print*,'doy: ',doy
! print*,'day: ',day
! print*,'hour: ',hour
! print*,'sec: ',sec_day
!endif

    !...calculate potential assimilations
    assimpot_omc = vegt%vmax*c4 + &
        vegt%vmax*(co2t%pco2i-co2t%gamma)/(co2t%pco2i+rrkk)*c3
    assimpot_ome = toa_apar*c4 + &
        toa_apar*(co2t%pco2i-co2t%gamma)/(co2t%pco2i+2.*co2t%gamma)*c3
    assimpot_oms = omss*(c3 + co2t%pco2i*c4)

    !...calculate assimilation rate stress factors
    assimfac(:) = 0.0
    IF (assimpot_omc .GT. minassim) THEN
        assimdiff = MAX(0.,assimpot_omc - assim_omc)
        assimfac(1) = MIN(1., MAX(0., &
                (assimdiff / assimpot_omc)))
    ENDIF
    assimfac(1) = 1. - assimfac(1)

    assimdiff = MAX(0.,assimpot_ome - assim_ome)
    IF (assimdiff .GT. minassimdiff) THEN
       assimfac(2) = MIN(1., MAX(0., &
              (assimdiff / assimpot_ome)))
    ENDIF
    assimfac(2) = 1. - assimfac(2)

    IF (assimpot_oms .GT. minassim) THEN
        assimdiff = MAX(0.,assimpot_oms - assim_oms)
        assimfac(3) = MIN(1., MAX(0., &
               (assimdiff / assimpot_oms)))
    ENDIF
    assimfac(3) = 1. - assimfac(3)

    assimfac(4) = assimfac(1) * &
        assimfac(2) * assimfac(3)

    !----------------------------------
    ! UPDATE STOMATAL CONDUCTANCE
    !----------------------------------
    !!!...Ball-Berry equation....!!!     
    co2s = MAX(co2t%pco2s,gprogt%pco2m*0.05) / pressure
    !print*,'pco2s ps :',co2t%pco2s
    !print*,'pco2m ps :',gprogt%pco2m
!    gsh2onew = (physcont%gradm *  &
!            MAX(1.0e-12_r8,co2t%assim)  &
!            * co2t%rstfac(1) * co2t%soilfrz / co2s) + bintc
    gsh2onew = (physcont%gradm *  &
            MAX(1.0e-12_r8,assimn)  &
            * co2t%rstfac(1) * co2t%soilfrz / co2s) + bintc

    !...calculate the change in stomatal resistance
    drst = co2t%rst * qdamp * ((gsh2o-gsh2onew)/  &
           (pdamp*gsh2o+qdamp*gsh2onew))

    co2t%rst = co2t%rst + drst

    !...update the smallest canopy stomatal resistance
    !...and convert it to s/m
    bintc = bintc * tc / ( amagatwv * tprcor)
    gminc = gminc * tc / ( amagatwv * tprcor)
    co2t%rst=MIN( 1./bintc, co2t%rst, 1./gminc )
    !co2t%rst=MIN( 1./bintc, co2t%rst )

    !if (co2t%rst > 1./bintc) then
    ! print*,'co2t%rst :', co2t%rst
    ! print*,'1./bintc :',1./bintc
    !endif

    !-----------------------------------
    ! SAVE CANOPY DIAGNOSTICS
    !-----------------------------------
    !...Carbon flux between CAS and reference level (mol C/m2/s)
    co2t%cflux = gsh2onew*(co2cas-co2m)

    !...CAN and CAS CO2 (Pa)
    co2t%pco2cas = co2cas * pressure
    co2t%co2cas = co2t%pco2cas / pressure    

    !...save the concentration variables
    co2t%co2s = dble(co2s)
    co2t%co2i = dble(co2t%pco2i / pressure)
    co2t%co2c = dble(co2t%pco2c / pressure)
    co2t%co2gamma = dble(co2t%gamma / pressure)
    co2t%co2m = dble(co2m)

    !print *,'co2cas out from phosib :',co2t%co2cas    
    !print *,'co2m out from phosib: ',co2t%co2m
    !print *,'co2s out from phosib: ',co2t%co2s
    !print *,'co2i out from phosib: ',co2t%co2i
    !print *,'co2c out from phosib: ',co2t%co2c
    !print *,'co2gamma out from phosib: ',co2t%co2gamma


    !...Potential top leaf photosynthesis
    ! Make assimn a top leaf, not the canopy
    assimnp = assimn / co2t%aparkk

    ! Bottom stopped assim
    antemp = MAX(dzero,assimnp)

    ! Potential intercellular CO2
    pco2ipot = pressure * (co2s-(1.6 * assimnp / &
               ((physcont%gradm * antemp / co2s) + bintc)))
    
    ! Potential rubisco limitation
    omcpot = vegt%vmax*2.1**qt*((pco2ipot-co2t%gamma) / &
               (pco2ipot+rrkk)*c3 + c4)

    ! Potential light limitation
    omepot = co2t%apar*((pco2ipot-co2t%gamma) / &
             (pco2ipot+2.*co2t%gamma)*c3 + c4)

    ! Potential sink or pep limitation
    omspot = (vegt%vmax / 2.)*(1.8**qt)*c3 &
              + rrkk*pco2ipot*c4

    ! Quad 1.
    sqrtin = MAX(dzero,((omepot+omcpot)**2 &
              - 4. * atheta * omepot * omcpot))

    ! Quad 1. Intermediate top leaf photosynthesis
    omppot = ((omepot+omcpot)-SQRT(sqrtin)) / &
              (2.*atheta)

    ! Quad 2.
    sqrtin = MAX(dzero,((omppot+omspot)**2 &
              -4.*btheta*omppot*omspot))

    ! Quad 2. Final potential top leaf photosynthesis
    if (omppot < 1.0E-14 ) omppot = dzero
    co2t%assimpot = ((omspot + omppot)-SQRT(sqrtin)) / &
                (2.*btheta)
    co2t%assimpot = MAX( co2t%assimpot, co2t%assim )

    !.. Call 13C fractionation code
    !IF (co2t%assim .GT. dzero) THEN
    !call cfrax_calc(physcont, gprogt, sibgl%co2t, &
    !     sibgl%fract, sibgl%fluxt, sibgl%co2t%co2cas, &
    !     sibgl%co2t%co2s, sibgl%co2t%co2i, sibgl%co2t%co2c, &
    !     sibgl%co2t%co2m, sibgl%co2t%co2gamma, sibgl%co2t%assim)

!print*,' '
!print*,'phosib ran, assim: ',co2t%assim
!print*,' '

    
end subroutine phosib
