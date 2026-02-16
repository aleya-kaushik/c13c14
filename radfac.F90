
!=================SUBROUTINE RADFAC======================================
subroutine radfac(chil, ref, tran, z1, z2,  &
       cosz, soilt, sscolt, vegt,  &
       hydrost, radt)

    !=======================================================================
    !
    !     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION
    !         ( DIRECT AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY
    !
    !     REFERENCE:
    !      Sellers, P.J., 1985: Canopy Reflectance, Photosynthesis and
    !                     Respiration. Int. J. Remote Sensing, 6(8)
    !                     1335-1372.
    !=======================================================================

    use kinds
    use module_pparams, only : &
        denh2o, denice, tice, &
        zlnd, gwlim
    use module_sibconst, only: &
        nsnow, nsoil, ntot
    use module_sib, only: &
        soil_type, sscol_type, &
        veg_type, &
        hydros_type, rad_type

    implicit none   

    !----------------------------------------------------------------------

    !...INPUT VARIABLES
    real(r4), intent(in) :: chil
    real(r4), dimension(2,2), intent(in) :: ref, tran
    real(r4), intent(in) :: z1, z2
    real(r8), intent(in) :: cosz
    type(soil_type), intent(in) :: soilt
    type(sscol_type), intent(in) :: sscolt
    type(veg_type), intent(in) :: vegt

    type(hydros_type), intent(inout) :: hydrost
    type(rad_type), intent(inout) :: radt

    !----------------------------------------------------------------------  

    !...LOCAL VARIABLES...
    integer(byte) :: irad  ! beam/diffuse loop variable
    integer(byte) :: iwave ! vis/nir loop variable
    real(r8), dimension(2) :: soref  ! local soref for loop
    real(r8), dimension(2,2) :: salb ! local total albedo                       

    real(r8) :: fff        ! minimum cosz value
    real(r8) :: facs       ! 1/20th of td(1) - Tice (0.0 - 0.4)
    real(r8) :: fmelt      ! 1-facs
    real(r8) :: scov       ! effective snow cover (0 - 0.5)
    real(r8) :: reff1, reff2  ! effective leaf reflectance (1=green, 2=brown)
    real(r8) :: tran1, tran2  ! effective leaf transmittance (1=green, 2=brown)
    real(r8) :: scat       ! average scattering coefficient (omega; eq. 1,2)
    real(r8) :: chiv       ! leaf angle distribution (-)
    real(r8) :: proj       ! leaf projection  (g(mu); eq. 13)
    real(r8) :: extkb      ! extinction coefficient (k, g(mu)/mu; eq. 1,2)
    real(r8) :: zmew       ! coefficient (INT(mu/g(mu); eq. 1,2)
    real(r8) :: acss       ! coefficient (A-S(mu); eq. 5)
    real(r8) :: upscat     ! leaf upwards scattering coefficient (omega*beta; eq. 3)
    real(r8), dimension(2) :: tranc1, tranc2, tranc3 ! canopy transmittances
    real(r8), dimension(2,2) :: albedoc, albedog  !canopy and ground albedo
                           !   (1,1) - visible, beam
                           !   (1,2) - visible, diffuse
                           !   (2,1) - nir, beam
                           !   (2,2) - nir, diffuse

    !...MISC COEFFICIENTS...
    !.....Intermediate variables identified in appendix
    real(r8) :: aa, bb
    real(r8) :: be, ce, de, fe, ge
    real(r8) :: bot, betao
    real(r8) :: hh1, hh2, hh3, hh4, hh5, &
                hh6, hh7, hh8, hh9, hh10
    real(r8) :: psi, zat
    real(r8) :: power1, power2
    real(r8) :: epsi, ek
    real(r8) :: f1, zp, den
    real(r8) :: zmk

    !--------------------------------------------------------------
    !...save local variables
    soref(1) = soilt%soref_vis
    soref(2) = soilt%soref_nir

    !...calculate minimum cosz
    fff = max(0.01746_r8,cosz)

    !...calculate snow vertical veg coverage (0-1)
    !.......the 5 right now is used to represent the average 
    !.......ratio of snow depth to water content 
    hydrost%snow_cvfc = ((hydrost%capacc_snow / denice)*5. - z1)/ &
                        (z2-z1)
    hydrost%snow_cvfc = MIN(done, MAX(dzero, hydrost%snow_cvfc))

    !...calculate snow vertical ground coverage (0-1)
    hydrost%snow_gvfc = min(done, hydrost%snow_gdepth / (zlnd*5.0))

    !...calculate saturation capacity depth (kg/m2)
    hydrost%satcapc = vegt%lait * 0.0001 * (done - hydrost%snow_cvfc) * denh2o
    hydrost%satcapg = gwlim

    !...calculate melting fraction
    facs = (sscolt%td(sscolt%nsl+1) - tice) * 0.04
    facs  = max( dzero, min( 0.4_r8, facs ) )
    fmelt = done - facs

    !-----------------------------------------------------------------------
    !...loop over 1=shortwave, 2=nir radiation
    do iwave = 1, 2
        scov =  min(0.5_r8, hydrost%wetfracc)
        reff1 = ( 1. - scov ) * ref(iwave,1) + &
                scov * ( 1.2 - iwave * 0.4 ) * fmelt
        reff2 = ( 1. - scov ) * ref(iwave,2) + &
                scov * ( 1.2 - iwave * 0.4 ) * fmelt
        tran1 = tran(iwave,1) * ( 1. - scov ) + &
                scov * ( 1.- ( 1.2 - iwave * 0.4 ) * fmelt )  &
                     * tran(iwave,1)
        tran2 = tran(iwave,2) * ( 1. - scov ) + &
                scov * ( 1.- ( 1.2 - iwave * 0.4 ) * fmelt ) * 0.9   &
                     * tran(iwave,2)

        !...calculate average scattering coefficient, leaf projection, 
        !.....and other coefficients
        scat = vegt%green*(tran1 + reff1) + &
               (1.- vegt%green) * (tran2 + reff2)
        chiv = chil

        if ( abs(chiv) .LE. 0.01 ) chiv = 0.01
        aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
        bb = 0.877 * ( 1. - 2. * aa )

        proj = aa + bb * fff
        extkb = ( aa + bb * fff ) / fff
        zmew = 1. / bb * ( 1. - aa / bb   &
            * log ( ( aa + bb ) / aa ) )
        acss = scat / 2. * proj / ( proj + fff * bb )
        acss = acss * ( 1. - fff * aa     &
            / ( proj + fff * bb ) * log ( ( proj   &
            +   fff * bb + fff * aa ) / ( fff * aa ) ) )

        upscat = vegt%green * tran1 + ( 1.- vegt%green ) * tran2
        upscat = 0.5 * ( scat + ( scat - 2. * upscat ) *   &
            (( 1. - chiv ) / 2. ) ** 2 )
        betao = ( 1. + zmew * extkb )   &
            / ( scat * zmew * extkb ) * acss

        !...calculate intermediate variables
        be = 1. - scat + upscat
        ce = upscat
        bot = ( zmew * extkb ) ** 2 + ( ce**2 - be**2 )

        if ( abs(bot) <= 1.e-10) then
            scat = scat* 0.98
            be = 1. - scat + upscat
            bot = ( zmew * extkb ) ** 2 + ( ce**2 - be**2 )
        endif

        de = scat * zmew * extkb * betao
        fe = scat * zmew * extkb * ( 1. - betao )
        hh1 = -de * be + zmew * de * extkb - ce * fe
        hh4 = -be * fe - zmew * fe * extkb - ce * de

        psi = sqrt(be**2 - ce**2)/zmew

        zat = vegt%lait/vegt%vcover*(1. - hydrost%snow_cvfc)

        power1 = min( psi*zat, 50.0_r8 )
        power2 = min( extkb*zat, 50.0_r8 )
        epsi = exp( - power1 )
        ek = exp ( - power2 )

        !...calculate ground albedos
        !......currently the diffuse and direct are the same
        !......making ge=1
        albedog(iwave,:) = &
            soref(iwave)*(1.-hydrost%snow_gvfc) + &
            ( 1.2-iwave*0.4 )* fmelt * hydrost%snow_gvfc
        ge = 1.0

        !...calculate canopy diffuse albedo
        f1 = be - ce / albedog(iwave,2)
        zp = zmew * psi

        den = ( be + zp ) * ( f1 - zp ) / epsi -   &
            ( be - zp ) * ( f1 + zp ) * epsi
        hh7 = ce * ( f1 - zp ) / epsi / den
        hh8 = -ce * ( f1 + zp ) * epsi / den
        f1 = be - ce * albedog(iwave,2)
        den = ( f1 + zp ) / epsi - ( f1 - zp ) * epsi

        hh9 = ( f1 + zp ) / epsi / den
        hh10 = - ( f1 - zp ) * epsi / den
        tranc2(iwave) = hh9 * epsi + hh10 / epsi
        albedoc(iwave,2) =  hh7 + hh8

        !...calculate canopy transmittance and direct albedo
        f1 = be - ce / albedog(iwave,2)
        zmk = zmew * extkb

        den = ( be + zp ) * ( f1 - zp ) / epsi -     &
            ( be - zp ) * ( f1 + zp ) * epsi
        hh2 = ( de - hh1 / bot * ( be + zmk ) )      &
            * ( f1 - zp ) / epsi -                   &
            ( be - zp ) * ( de - ce*ge - hh1 / bot   &
            * ( f1 + zmk ) ) * ek
        hh2 = hh2 / den
        hh3 = ( be + zp ) * (de - ce * ge -          &
            hh1 / bot * ( f1 + zmk )) * ek -         &
            ( de - hh1 / bot * ( be + zmk ) ) *      &
            ( f1 + zp ) * epsi
        hh3 = hh3 / den
        f1 = be - ce * albedog(iwave,2)
        den = ( f1 + zp ) / epsi - ( f1 - zp ) * epsi
        hh5 = - hh4 / bot * ( f1 + zp ) / epsi -     &
            ( fe + ce*ge*albedog(iwave,2) +         &
            hh4 / bot * ( zmk - f1 ) ) * ek
        hh5 = hh5 / den
        hh6 =   hh4 / bot * ( f1 - zp ) * epsi +     &
            ( fe + ce * ge * albedog(iwave,2) +     &
            hh4 / bot*( zmk - f1 ) ) * ek
        hh6 = hh6 / den
        tranc1(iwave) = ek
        tranc3(iwave) = hh4 / bot * ek + hh5 * epsi + hh6 / epsi
        albedoc(iwave,1) = hh1 / bot + hh2 + hh3

        !...calculate radfac, which multiplies incoming short wave fluxes
        !......to give absorption of radiation by canopy and ground
        radt%radfacg(iwave,1) = &
            (1. - vegt%vcover) * (1. - albedog(iwave,1)) +  & 
            vegt%vcover * (tranc1(iwave)*( 1.-albedog(iwave,1))  &
                               + tranc3(iwave)*(1.-albedog(iwave,2)))

        radt%radfacg(iwave,2) = &
            (1. - vegt%vcover) * (1. - albedog(iwave,2)) + &
            vegt%vcover * (tranc2(iwave) * (1.-albedog(iwave,2)))

        radt%radfacc(iwave,1) = vegt%vcover * &
              ( ( 1.-albedoc(iwave,1) )  &
                 - tranc1(iwave) * (1.-albedog(iwave,1))      &
                 - tranc3(iwave) * (1.-albedog(iwave,2)) )

        radt%radfacc(iwave,2) = vegt%vcover * &
             ( ( 1.-albedoc(iwave,2) ) &
                - tranc2(iwave) * ( 1.-albedog(iwave,2) ) )

    enddo   !iwave loop

   !...calculate total surface albedos
   do iwave = 1,2
      do irad = 1,2
          salb(iwave,irad) = (1.-vegt%vcover) * albedog(iwave,irad) + &
                              vegt%vcover * albedoc(iwave,irad)
      enddo
   enddo

   radt%albedo_visb = salb(1,1)
   radt%albedo_visd = salb(1,2)
   radt%albedo_nirb = salb(2,1)
   radt%albedo_nird = salb(2,2)

   !-----------------------------------------------------------------------
   !...calculate surface temperature based on ground/snow cover/temps
   radt%tsfc = min(tice,sscolt%td(sscolt%nsl+1)) * hydrost%snow_gvfc &
              + sscolt%td(1)*(1.-hydrost%snow_gvfc)

    !...the gap fraction/effective ground cover
    !.....reduced down to a linear relationship
    !.....between LAI and absorbed canopy thermal radiation
    !...calculate PFT gap fraction/effective ground cover 
    radt%effgc = MIN(done, 0.4*vegt%lait)

end subroutine radfac
