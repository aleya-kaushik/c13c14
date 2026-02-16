 
!===================SUBROUTINE DELHF25=====================================
 
subroutine delhf(tm, bps, ros, &
                 ta, nsl, tsfc, &
                 fss, hc, hg, hs, &
                 tc, ra, rb, rd)

use kinds
use module_local
use module_pparams, only: cp => spec_heat_cp
use module_time, only: dtsib

!========================================================================
!
!     Calculation of partial derivatives of canopy and ground sensible
!        heat fluxes with respect to Tc, Tg, and Theta-m.
!     Calculation of initial sensible heat fluxes.
!
!======================================================================== 


!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       HC             CANOPY SENSIBLE HEAT FLUX (J M-2)
!       HG             GROUND SENSIBLE HEAT FLUX (J M-2)
!       HS             SNOW   SENSIBLE HEAT FLUX (J M-2)
!       HA             CAS    SENSIBLE HEAT FLUX (J M-2)
!       HCDTC          dHC/dTC 
!       HCDTA          dHC/dTA
!       HGDTG          dHG/dTG
!       HGDTA          dHG/dTA
!       HSDTS          dHS/dTS
!       HSDTA          dHS/dTA
!       HADTA          dHA/dTA
!       HADTH          dHA/dTH
!       AAC            dH/dTC
!       AAG            dH/dTG
!       AAM            dH/dTH
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

!----------------------------------------------------------------------

!...input variables
real(r8), intent(in) :: tm, bps(2), ros
integer(byte), intent(in) :: nsl
real(r8), intent(in) :: ta, tsfc

real(r8), intent(inout) :: fss, hc, hg, hs, tc, ra, rb, rd

!...local variables
real(r8) :: rai    ! 1/ra
real(r8) :: rbi    ! 1/rb
real(r8) :: rdi    ! 1/rd

!-----------------------------------------------------------------------
!                                                                       
!     FLUXES EXPRESSED IN JOULES/M2, although in SIBSLV WE THEN WANT W/m2
!
!      HC          (HC)    : EQUATION (63) , SE-86
!      HG          (HG)    : EQUATION (65) , SE-86
!      HS          (HS)    : EQUATION (65) , SE-86
!-----------------------------------------------------------------------

    rai = 1.0 / ra
    rbi = 1.0 / rb  
    rdi = 1.0 / rd              

    !    these are the current time step fluxes in J/m2
    hc   = cp * ros * (tc - ta) * rbi * dtsib

    if (nsl == 0) then !no snow case
        hg   = cp * ros * (tsfc - ta) * rdi * dtsib 
        hs   = 0.0
    else
        hg = 0.0
        hs   = cp * ros * (tsfc - ta) * rdi * dtsib 
    endif

    fss  = cp * ros * (ta - tm) * rai * dtsib

    !    now we do the partial derivatives
    !    these are done assuming the fluxes in W/m2        
    !    for canopy leaves sensible heat flux: W/(m2 * K)
    ! 
    hcdtc =   cp * ros * rbi
    hcdta = - hcdtc
    !
    !    for ground and snow sensible heat fluxes: W/(m2 * K)
    !
    hgdtg =   cp * ros * rdi  
    hsdts =   hgdtg
    hgdta = - hgdtg
    hsdta = - hgdtg
    !
    !    for the canopy air space (CAS) sensible heat flux: W/(m2 * K)
    !
    hadta =   cp * ros * rai
    hadth = - hadta/bps(1)

end subroutine delhf                                              
