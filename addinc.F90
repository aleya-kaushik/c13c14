!-SUBROUTINE: addinc-------------------------------------------------
! Add prognostic variable increments to prognostic variables and
! diagnose heat fluxes and mixing ratios.
!----------------------------------------------------------------
!
!                           OUTPUT
!
!   tcas    CAS Temperature (K)
!   tc      Canopy Temperature (K)
!   eacas   CAS Vapor Preassure (hPa or mb)
!   td      Deep Soil Temperature (K)
!   fss     CAS to Mixed Layer Sensible Heast Flux (W m^-2)
!   fws     CAS to Mixed Layer Latent Heat Flux (W m^-2)
!   sh      Mixed Layer Mixing Ratio (kg/kg)
!   shcas   CAS Mixing Ratio (kg/kg)
!
!----------------------------------------------------------------
subroutine addinc (gref, pftref, &
                   lonsib, latsib,  &
                   gdiagt, gprogt, &
                   cast, fluxt, sscolt)

use kinds
use module_local, only: &
   dea, dta, dtc, dtd
use module_pparams, only: &
   spec_heat_cp
use module_sib, only: &
   gdiag_type, gprog_type, &
   cas_type, flux_type, sscol_type
use module_sibconst, only: nsoil

implicit none

!-Input--------------
integer(i4), intent(in) :: gref, pftref
real(r4), intent(in) :: lonsib, latsib
type(gdiag_type), intent(inout) :: gdiagt
type(gprog_type), intent(inout) :: gprogt
type(cas_type), intent(inout)   :: cast
type(flux_type), intent(inout)  :: fluxt
type(sscol_type), intent(inout) :: sscolt

!-Local Variables----
integer(i4) :: j

!----------------------------------------------------------------

    cast%eacas = cast%eacas + dea
    cast%tcas = cast%tcas + dta
    cast%tc   = cast%tc + dtc
    
    do j = sscolt%nsl+1, nsoil
        sscolt%td(j) = sscolt%td(j) + dtd(j)
        if (sscolt%td(j) < 100. .or. sscolt%td(j) > 400.) then
            print*, 'Point/lon/lat: ',gref, lonsib, latsib
            print*, 'PFT: ',pftref
            print*, 'BAD td VALUE AT LEVEL:',j
            print*, 'td: ', sscolt%td(j)-dtd(j),sscolt%td(j)
            stop
        endif
    enddo

    if ( gdiagt%em <= 0.) gdiagt%em=1.e-3
    if ( cast%eacas <= 0.) cast%eacas = gdiagt%em
    if ( cast%tc < 100.0 .or. cast%tcas < 100.0 .or. &
         cast%tc > 500. .or. cast%tcas > 400.) then

        print*, 'Point/lon/lat: ',gref, lonsib, latsib
        print*, 'PFT: ',pftref
        print *, 'BAD ta or tc VALUE:'
        print *, 'ta: ', cast%tcas-dta, cast%tcas  
        print *, 'tc: ', cast%tc-dtc, cast%tc
        print *, 'tsnow:  ', &
              sscolt%td(sscolt%nsl+1) - dtd(sscolt%nsl+1), &
              sscolt%td(sscolt%nsl+1)
        print*,  'tsoil1: ', &
              sscolt%td(1)-dtd(1),sscolt%td(1)
        print *, 'tnsoil: ', sscolt%td(nsoil)-dtd(nsoil), &
              sscolt%td(nsoil)
        print *, 'ea: ', cast%eacas-dea, cast%eacas
        print *, 'em: ', gdiagt%em
        print *, ''

        stop
    endif

    ! Calucluate latent and sunsible fluxes between CAS and
    ! mixed (boundry) layer
    fluxt%fss = gdiagt%ros * spec_heat_cp &
                  * (cast%tcas - gprogt%tm) &
                  / fluxt%ra
    fluxt%fws = (cast%eacas - gdiagt%em) &
                 / fluxt%ra * spec_heat_cp &
                  * gdiagt%ros / gdiagt%psy

    ! Recalculate mixing ratios
    gprogt%sh  = 0.622 / (gprogt%ps / gdiagt%em - 1.)
    cast%shcas = 0.622 / (gprogt%ps / cast%eacas - 1.)
    
end subroutine addinc

