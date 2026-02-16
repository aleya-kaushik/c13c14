!----------------------------------------------------------------------

subroutine snow_compact(nsl, td, &
           www_liq, www_ice, dz)

!----------------------------------------------------------------------
!
!   Based on CLM subroutine clm_compact
!
!   CLM web info : http://clm.gsfc.nasa.gov
!
!   Description:
!   Three metamorphisms of changing snow characteristics are 
!   implemented, i.e., destructive, overburden and melt. The
!   treatments of the former two are from SNTHERM.89 and SNTHERM.99 
!   (1991, 1999). The contribution due to melt metamorphism is 
!   simply taken as a ratio of snow ice fraction after the melting 
!   versus before the melting.
!
!   Revision History:
!   15 September 1999: Yongjiu Dai; initial code
!   15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!   30 January   2002: Ian Baker, SiB integration

use kinds
use module_local, only: &
    frac_iceold, imelt
use module_oparams, only: &
    snow_c2, snow_c3, snow_c4, &
    snow_c5, snow_dm, snow_eta0
use module_pparams, only: &
    denice, denh2o, tice
use module_time, only: &
    dtsib, dtisib

implicit none

!----------------------------------------------------------------------
!...input variables
integer(byte), intent(in) :: nsl
real(r8), dimension(nsl+1:0), intent(in) :: td, www_liq, www_ice
real(r8), dimension(nsl+1:0), intent(inout) :: dz

!...local variables
integer(i4)        :: j

real(r8) :: burden    ! pressure of overlying snow (kg m^-2)
real(r8) :: wx        ! water mass (ice + liquid) (kg m^-2)
real(r8) :: void      ! void = 1 - vol_ice - vol_liquid
real(r8) :: bi        ! partial density of ice (kg m^-3)
real(r8) :: fi        ! fraction of ice relative to total
!  water content
real(r8) :: delt      ! snow sib%prog%td - tice (K)
real(r8) :: dexpf     ! expf = exp(-c4*(tice-sib%prog%td))
real(r8) :: ddz1      ! rate of settling snowpack due to 
!  destructive metamorphism
real(r8) :: ddz2      ! rate of compaction of snowpack due
!  to overburden
real(r8) :: ddz3      ! rate of compaction of snowpack due
!  to melt
real(r8) :: pdzdtc    ! nodal rate of change in fractonal
!  thickness due to compaction (fraction sec^-1)


!----------------------------------------------------------------------

    burden = 0.0

    do j=nsl+1,0

        wx = www_ice(j) + www_liq(j)
        void = 1.0 - (www_ice(j)/denice + www_liq(j)/ denh2o)/dz(j)

        !...disallow compaction for water saturated node and lower ice lens node
        if(void <= 0.001  .or.  www_ice(j) <= 0.1) then
            burden = burden + wx
            cycle
        endif

        bi = www_ice(j)/dz(j)
        fi = www_ice(j)/wx
        delt = tice - td(j)
        dexpf = exp(-snow_c4*delt)

        !...settling as a result of desctructive metamorphism
        ddz1 = -snow_c3*dexpf
        if(bi > snow_dm) ddz1 = ddz1*exp(-46.e-3*(bi-snow_dm))

        !...liquid water term
        if(www_liq(j) > 0.01*dz(j)) ddz1 = ddz1*snow_c5

        !...compaction due to overburden
        ddz2 = -burden*exp(-0.08*delt-snow_c2*bi)/snow_eta0

        !...compaction occurring during melt
        if(imelt(j) == 1 )then
            ddz3 = -1.*dtisib * max(0.0_r8, &
                (frac_iceold(j) - fi)/frac_iceold(j))
        else
            ddz3 = 0.0
        endif

        !...time rate of fractional change in dz (units of sec-1)
        pdzdtc = ddz1 + ddz2 + ddz3

        !...change in dz due to compaction
        dz(j) = dz(j) * (1.0+pdzdtc*dtsib)

        !...pressure of overlying snow
        burden = burden + wx

    enddo ! nsnow loop


end subroutine snow_compact
