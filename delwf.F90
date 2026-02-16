
!===================SUBROUTINE DELWF=====================================

subroutine delwf()

    use kinds
    use module_local

    !========================================================================
    !
    !     Calculation of partial derivatives of canopy and ground radiative
    !        heat fluxes with respect to Tc, Tg
    !     Here we are doing only the long wave radiative loss, which is the
    !     only radiative quantity we are trying to bring to the next time step.
    !
    !======================================================================== 

    !------------------------------INPUT is coming from Netrad-------------
    !
    !       dtc4, dtg4, dts4, which are the derivatives of the LW loss
    !
    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
    !
    !       LCDTC          dLC/dTC 
    !       LCDTG          dLC/dTG
    !       LCDTS          dLC/dTS
    !       LGDTG          dLG/dTG
    !       LGDTC          dLG/dTC
    !       LSDTS          dLS/dTS
    !       LSDTC          dLS/dTC
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT none

    !----------------------------------------------------------------------

    !   canopy leaves:
    lcdtc = 2 * dtc4 * fac1
    lcdtg =     - dtg4 * fac1
    lcdts =     - dts4 * fac1
    !
    !   ground:
    !
    lgdtg =   dtg4
    lgdtc = - dtc4 * fac1
    !
    !   snow:
    !
    lsdts =   dts4
    lsdtc = - dtc4 * fac1

end subroutine delwf
