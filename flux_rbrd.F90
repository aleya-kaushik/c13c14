!==================SUBROUTINE RBRD=======================================
subroutine flux_rbrd(z2, u2, &
                tc, ta, td, &
                snow_cvfc, &
                lai, ccc1, ccc2, &
                rbc, rdc, rb, rd)

use kinds
use module_pparams, only: grav   
use module_sibconst, only: nsnow, nsoil

implicit none

!----------------------------------------------------------------------

!...input variables
real(r4), intent(in) :: z2
real(r8), intent(in) :: u2
real(r8), intent(in) :: tc, ta, td
real(r8), intent(in) :: snow_cvfc
real(r8), intent(in) :: lai, ccc1, ccc2
real(r8), intent(inout) :: rbc, rdc
real(r8), intent(inout) :: rb, rd

!----------------------------------------------------------------------  

!      Reference

!      Sellers, P.J. and Mintz, Y., Y.C. Sud, A. Dalcher, 1986: A Simple 
!                     Biospher Model (SiB) for use Within General 
!                     Circulation Models. JAS, 43(6),505-531.

!========================================================================
!
!      CALCULATION OF RB AND RD AS FUNCTIONS OF U2 AND TEMPERATURES
!
!======================================================================== 


!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       RB (GRB)       CANOPY TO CAS AERODYNAMIC RESISTANCE (Sec M-1)
!       RD (GRD)       GROUND TO CAS AERODYNAMIC RESISTANCE (Sec M-1)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!Bio...LOCAL VARIABLES
real(r8) :: temdif ! vegetation-CAS temperature difference (K)
real(r8) :: fac    ! factor for vegetation-CAS resistance calc
real(r8) :: fih    ! factor for ground-CAS resistance calc

    !-----------------------------------------------------------------------
    !      RBC and RDC         : UPDATE FOR SNOW
    !-----------------------------------------------------------------------
    rbc = ccc1 / (1. - snow_cvfc)
    rdc = ccc2 * (1. - snow_cvfc)

    !-----------------------------------------------------------------------
    !      RB       (RB)       : EQUATION (A9), SE-86
    !-----------------------------------------------------------------------

    temdif  = MAX( 0.01_r8,  tc - ta)
    fac     = lai / 890.* (temdif * 20.0)**0.25
    rb  = 1.0 / ( (SQRT(u2) / rbc) + fac )

    !-----------------------------------------------------------------------
    !      RD       (RD)       : EQUATION (A15), SE-86
    !-----------------------------------------------------------------------

    !...for rd, we use the composite soil/snow temperature
    temdif = MAX( 0.1_r8, td - ta )
    fih = SQRT( 1.+9.* grav * temdif * z2 / (td*u2*u2) )
    rd  = rdc / (u2 * fih) 

end subroutine flux_rbrd
