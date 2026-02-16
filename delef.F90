
!==================SUBROUTINE DELEF======================================

subroutine delef(psy, ros, em, ea, soilrh, &
                 rd, ra, &
                 ec, eg, es, fws)

    use kinds
    use module_local
    use module_pparams, only: &
        cp => spec_heat_cp, &
        snofac
    use module_time, only: dtsib

    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground latent
    !        heat fluxes with respect to Tc, Tg, Theta-m, and Qm.
    !     Calculation of initial latent heat fluxes.
    !
    !     --The ETC, ETG are vapor pressures at TC, TG
    !     --The BTC, BTG are derivatives of ETC, ETG with relation to TC, TG
    !
    !======================================================================== 

    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       EC             ECT + ECI
    !       EG             EGS + EGI
    !       ECDTC          dEC/dTC
    !       ECDTG          dEC/dTG
    !       ECDQM          dEC/dQM
    !       EGDTC          dEG/dTC
    !       EGDTG          dEG/dTG
    !       EGDQM          dEG/dQM
    !       BBC            dE/dTC
    !       BBG            dE/dTG
    !       BBM            dE/dQM
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    implicit none

    !----------------------------------------------------------------------

    !...input variables
    real(r8), intent(in) :: psy, ros, em, ea, soilrh
    real(r8), intent(in) :: rd, ra

    real(r8), intent(inout) :: ec, eg, es, fws

    !...Local variables
    real(r8) :: cpdpsy    ! cp/psy
    cpdpsy = cp / psy

    !-----------------------------------------------------------------------
    !                                                                       
    !     CALCULATION OF SURFACE RESISTANCE COMPONENTS, SEE EQUATIONS (64,66)
    !       OF SE-86
    !       gect(i)  =      (1. -wc(i)) /  rc(i)
    !       geci(i)  = epsc(i) * wc(i)  / (RB(I) + RB(I))
    !       gegs(i)  =        (1-wg(i)) / (rds(i))
    !       gegi(i)  = epsg(i) * wg(i)  /  rd(i)
    !                                                                       
    !-----------------------------------------------------------------------

    !... MODIFICATION FOR SOIL DRYNESS : REL. HUMIDITY
    cog1 = (gegi + gegs*soilrh )
    cog2 = (gegi + gegs        )

    !-----------------------------------------------------------------------
    !                                                                       
    !     FLUXES EXPRESSED IN JOULES M-2 
    !                                                                       
    !      ec         (EC)    : EQUATION (64) , SE-86
    !      eg         (EG)    : EQUATION (66) , SE-86
    !      es         (ES)    : EQUATION (66) , SE-86
    !      ea         (EA)    : EQUATION ????
    !-----------------------------------------------------------------------

    !...current fluxes in J/m2
    ec  = (etc - ea) * coc * ros  * cpdpsy * dtsib
    eg  = ( etg * cog1 - ea * cog2) &
            * ros * cpdpsy * dtsib
    es  = (ets - ea) / rd * ros * cpdpsy/snofac * dtsib

    fws = (ea  - em) / ra * ros * cpdpsy * dtsib

    !...fluxes for derivaties not accumulated, so W/m2
    !....for the canopy leaves vapor pressure: W/ (m2* K)
    ecdtc = getc * coc * ros * cpdpsy
    ecdea = - coc * ros * cpdpsy

    !....for ground latent heat fluxes: W/ (m2* K)
    egdtg = getg * cog1 * ros * cpdpsy
    egdea = - cog2 * ros * cpdpsy               

    !....for snow latent heat fluxes: W/ (m2* K)
    esdts =   gets * ros * cpdpsy / rd / snofac

    !....for snow latent heat fluxes: W/ (m2 * Pa)
    esdea = - ros * cpdpsy / rd / snofac        

    !pl...for CAS latent heat fluxes: W/ (m2* Pa)
    eadea = ros * cpdpsy / ra 
    eadem = - eadea            

end subroutine delef                                              
