!==================SUBROUTINE NETRAD====================================
!        CALCULATES NET RADIATION, RADT, USING RADIATION FROM 
!         PHYSICS AND CURRENT LOSSES FROM CANOPY AND GROUND          
!=======================================================================
subroutine radnet ( nsl, &
    tsfc, tc, radt)

    use kinds
    use module_local
    use module_pparams, only: &
        stefan, tice
    use module_sib, only: rad_type

    implicit none

    !...input variables
    integer(byte), intent(in) :: nsl
    real(r8), intent(in) :: tsfc, tc
    type(rad_type), intent(inout) :: radt

    !...local variables
    real(r8) :: &
        tg,      & ! ground temperature (K)
        ts,      & ! snow temperature (K)
        tc4,     & ! canopy temp **4
        tg4,     & ! ground temp **4
        ts4        ! snow temp **4

    !------------------------------------------
    !...save effective ground cover in local module
    fac1 = radt%effgc

    !...calculate temperatures
    tg = tsfc
    ts = tsfc

    tc4 = tc**4
    tg4 = tg**4
    ts4 = ts**4

    !...derivatives
    dtc4 = 4*stefan * tc**3
    dtg4 = 4*stefan * tg**3
    dts4 = 4*stefan * ts**3

    !...canopy leaves thermal radiation contributions
    closs =  2 * fac1 * stefan * tc4

    if (nsl == 0) then
       closs = closs - fac1 * stefan * tg4
    else
       closs = closs - fac1 * stefan * ts4
    endif 

    !...ground thermal radiation contributions
    gloss =  stefan * tg4 - fac1 * stefan * tc4

    !...snow thermal radiation contributions
    sloss = stefan * ts4 - fac1 * stefan * tc4

    !...net radiation terms
    radt%radtc = radt%radc3c - closs
    radt%radtg = radt%radc3g - gloss
    radt%radts = radt%radc3g - sloss

end subroutine radnet
