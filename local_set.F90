!--------------------------------------------------------
subroutine local_set(ps, tc, &
            capacc_liq, capacc_snow, &
            capacg, snow_gmass, &
            nsl, www_liq, www_ice, td)

!--------------------------------------------------------

!Sets the local sib variables to values
!   ready for the next timestep.

use kinds
use module_local
use module_pparams, only: &
    tice
use module_sibconst, only: &
    nsoil, nsnow


implicit none

!...input variables
real(r8), intent(in) :: ps, tc
real(r8), intent(in) :: capacc_liq, capacc_snow
real(r8), intent(in) :: capacg, snow_gmass
integer(byte), intent(in) :: nsl
real(r8), dimension(-nsnow+1:nsoil), intent(in) :: &
     www_liq, www_ice, td

!...local variables
integer(i4) :: j
real(r8) :: toth2o_old

!---------------------------------------------------------

     !...previous time-step values
     capaccl_old = capacc_liq
     capaccs_old = capacc_snow
     capacg_old = capacg
     snow_gmass_old =snow_gmass
     nsl_old = nsl

     td_old(:) = td(:)
     wwwliq_old(:) = www_liq(:)
     wwwice_old(:) = www_ice(:)

     do j=nsl+1,nsoil
        toth2o_old = (wwwice_old(j) + wwwliq_old(j))
        if (toth2o_old > 0.0) then
            frac_iceold(j) = wwwice_old(j) / toth2o_old
        endif
     enddo

     !...saturation vapor pressures for canopy, soil, and snow
     !......(vapor pressures converted from Pa to hPa/mb)
     ppl(1) = ps*100.0
     ttl(1) = tc
     call dtess_eau(1,ppl,ttl,esst,dtesst)

     etc = esst(1)/100.0
     getc = dtesst(1)/100.0

     ttl(1) = td(1)
     call dtess_eau(1,ppl,ttl,esst,dtesst)
     etg = esst(1)/100.0
     getg = dtesst(1)/100.0

     ttl(1) = td(nsl+1)  !min(tice, td(nsl+1))
     call dtess_eau(1,ppl,ttl,esst,dtesst)
     ets = esst(1)/100.0
     gets = dtesst(1)/100.0

     !...canopy air space (CAS) and canopy
     dtg = dzero
     dts = dzero
     dtc = dzero
     dth = dzero
     dqm = dzero
     dta = dzero
     dea = dzero

     gect = dzero
     geci = dzero
     gegs = dzero
     gegi = dzero
     coc = dzero
     cog1 = dzero
     cog2 = dzero

     !...radiation
     dtc4 = dzero
     dtg4 = dzero
     dts4 = dzero
     lcdtc = dzero
     lcdtg = dzero
     lcdts = dzero
     lgdtc = dzero
     lgdtg = dzero
     lsdts = dzero
     lsdtc = dzero
     hcdtc = dzero
     hcdta = dzero
     hgdta = dzero
     hgdtg = dzero
     hsdta = dzero
     hsdts = dzero
     hadta = dzero
     hadth = dzero
     ecdtc = dzero
     ecdea = dzero
     egdtg = dzero
     egdea = dzero
     esdts = dzero
     esdea = dzero
     eadea = dzero
     eadem = dzero
     closs = dzero
     gloss = dzero
     sloss = dzero
     fac1 = dzero

     !...soil/snow
     dtd(:) = dzero
     imelt(:) = dzero

     !...energy balance
     radttc = dzero
     radttg = dzero
     radtts = dzero
     cas_e_storage = dzero
     cas_w_storage = dzero


end subroutine local_set
