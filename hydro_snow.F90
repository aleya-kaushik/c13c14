subroutine hydro_snow(poros, hydrost, sscolt)

!----------------------------------------------------------------------
!
!     Evaluate and update snow mass/snow water.
!                 based on code in CLM_HYDRO_SNOW.F90
!
!     The following is taken verbatim from the comments in CLM_HYDRO_SNOW
!
!     Description:
!       Evaluate the change of snow mass and the snow water onto soil. 
!     Water flow within snow is computed by an expicit and non-physical
!     capacity (a tentative value is used, i.e. equal to 0.33*porosity)
!     to percolate into the underlying layer.  Except for cases where 
!     the porosity of one of the two neighboring layers is less than 0.05,
!     zero flow is assumed. The water flow out of the bottom of the snow
!     pack will participate as the input of the soil water and runoff.
!
!     Revision History:
!      15 September 1999: Yonjiu Dai; original code
!      15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!      15 November  2000: Mariana Vertenstein
!      22 January   2002: Ian Baker - integration into SiB
!      18 October   2015: Kathy Haynes - modified for SiB4
!----------------------------------------------------------------------

use kinds

use module_pparams, only: &
    denice, denh2o, &
    gwctog, ssi
use module_sib, only: &
    hydros_type, &
    sscol_type
use module_time, only: dtsib, dtisib


implicit none

!----------------------------------------------------------------------
!...input variables
real(r8), intent(in) :: poros
type(hydros_type), intent(inout) :: hydrost
type(sscol_type), intent(inout) :: sscolt

!...LOCAL VARIABLES...
integer(i4) :: j
real(r8)    :: qin, qout, qout_snow

!----------------------------------------------------------------------

!...no snow, save water and return
if (sscolt%nsl == 0) then  !no snow case
   hydrost%www_inflow = hydrost%pcpg_rain  ! units kg/m^2/sec 

!...snow exists, proceed through routine
else  !snow present

   !...Note: Not adjusting top snow layer liquid by fluxes.
   !.......not adding dew not subdracting sublimation
   !...Note: Snow was added to top snow layer in 
   !.......hydro_canopy.F90

   !...add rain to top snow layer
   j=sscolt%nsl+1
   sscolt%www_liq(j) = &
         sscolt%www_liq(j) + hydrost%pcpg_rain*dtsib

   !...porosity and partial volume - recalculated
   do j=sscolt%nsl+1,0
       sscolt%vol_ice(j) = min(poros, &
             sscolt%www_ice(j)/(sscolt%dz(j)*denice))
       sscolt%eff_poros(j) = 1.0 - sscolt%vol_ice(j)
       sscolt%vol_liq(j) = min(sscolt%eff_poros(j), &
             sscolt%www_liq(j)/(sscolt%dz(j)*denh2o))
   enddo

   !---this directly from the comments in CLM_HYDRO_SNOW---
   !...Capillary forces within snow are usually two or more orders of 
   !...magnitude less than those of gravity. Only gravity terms are
   !...considered.  The general expression for water flow is "K" * ss**3",
   !...however, no effective parameterization for "K".  Thus, a very
   !...simple consideration (not physically based) is introduced: when
   !...the liquid water of layer exceeds the layer's holding capacity,
   !...the excess meltwater adds to the underlying neighbor water.
   !...CLM sets 0.033 as irreducible water saturation for snow (ssi)
   qin = 0.0
   do j = sscolt%nsl+1,0
      sscolt%www_liq(j) = sscolt%www_liq(j) + qin
      if (j <= -1 ) then
         !...no runoff over snow surface, just ponding on surface
         if (sscolt%eff_poros(j) < 0.05 .or. &
             sscolt%eff_poros(j+1) < 0.05) then
             qout = 0.0
         else
             qout = max(0.0_r8,(sscolt%vol_liq(j)- &
                    ssi * sscolt%eff_poros(j)) * sscolt%dz(j))
             qout = min(qout,(1.0-sscolt%vol_ice(j+1) - &
                    sscolt%vol_liq(j+1))*sscolt%dz(j+1))
          endif
       else
          qout = max(0.0_r8,(sscolt%vol_liq(j) - &
                   ssi*sscolt%eff_poros(j))*sscolt%dz(j))
       endif

       qout           = qout * 1000.0
       sscolt%www_liq(j) = sscolt%www_liq(j) - qout
       qin            = qout
   enddo

   !...liquid out of the bottom of the snow into the top soil layer...
   qout_snow = qout*dtisib
   hydrost%www_inflow = qout_snow

endif   ! snow present/not present condition

!...put www_inflow into capac(2) - ground interception
!...before it is infiltrated.
hydrost%capacg = hydrost%capacg + hydrost%www_inflow * dtsib

!...canopy/ground wetness fractions
hydrost%wetfracg = MAX(dzero, MIN( done, &
         gwctog * (hydrost%capacg/hydrost%satcapg)))

end subroutine hydro_snow
