!Routine to update CAS variables
!
subroutine cas_update(bps1, psy, ros, &
     casd, capacc_liq, capacc_snow, lai, &
     cast)

use kinds
use module_pparams, only: &
    denh2o, denice, &
    leafhc, h2ohc, &
    spec_heat_cp, cv
use module_sib, only: &
    cas_type

implicit none

!...input variables
real(r8), intent(in) :: bps1, psy, ros
real(r8), intent(in) :: capacc_liq, capacc_snow
real(r8), intent(in) :: casd, lai
type(cas_type), intent(inout) :: cast


    !...canopy heat capacity
    cast%hcapc = lai*leafhc + &
                 (capacc_snow/denice + &
                  capacc_liq/denh2o) * h2ohc

    !...canopy potential temperature and water vapor
    cast%thcas = cast%tcas / bps1

    !...canopy storage inertia terms
    cast%hcapcas = ros * spec_heat_cp * casd
    cast%vcapcas = ros * cv * casd / psy

end subroutine cas_update
