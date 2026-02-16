!!!Calculation of fluorescence (SIF).
!!!
!!!These equations follow: 
!!!van der Tol et al. (2014)
!!!
subroutine sif_calc(nspar, fpar, &
      assim, assimpot, tc, vmax, sift)

use kinds
use module_oparams, only: &
    sif_a1, sif_a2, &
    sif_kf0, sif_kn0, sif_kp0
use module_pparams, only: &
    mol_to_umol
use module_sib, only: sif_type

implicit none

!...source

!...input variables
real(r8), intent(in) :: nspar    !non-scattered PAR (mol/m2/s)
real(r8), intent(in) :: fpar     !absorbed fraction of PAR (-)
real(r8), intent(in) :: assim    !asimilation (mol C/m2/s)
real(r8), intent(in) :: assimpot !potential assimilation (mol C/m2/s)
real(r8), intent(in) :: tc       !canopy temperature (K)
real(r8), intent(in) :: vmax     !rubisco velocity (mol/m2/s)
type(sif_type), intent(inout) :: sift

!...local variables
real(r8) :: a3, po0, ps, totk, toty
real(r8) :: fo, fm, kappa

if (nspar .le. dzero) then
   sift%sif_je = dzero
   sift%sif_jo = dzero
   sift%sif_jejo = dzero

   sift%sif_kd = dzero
   sift%sif_kn = dzero
   sift%sif_kp = dzero

   sift%phi_d = dzero
   sift%phi_f = dzero
   sift%phi_n = dzero
   sift%phi_p = dzero

   sift%sif = dzero
else

   ! Potential and actual electron transport rate
   sift%sif_je = assim
   sift%sif_jo = assimpot
   IF (sift%sif_jo .GT. dzero) THEN
       sift%sif_jejo = sift%sif_je / sift%sif_jo
   ELSE
      sift%sif_jejo = dzero
   ENDIF
       
   !----Probabilities of excitons to follow pathways
   ! Heat dissipation probability
   !...vdt14 page 3, with temperature adjustment from
   !...Joe Berry personal communication
   if (tc < 300.) then
      sift%sif_kd = 0.95
   else
      sift%sif_kd = 0.95 + ((tc - 300.) * 0.0236)
   endif

   ! Non-photochemical quenching (NPQ) probability
   !...x-factor (vdt14, Eq. 16)
   !.....ZERO if GPP = potential GPP
   !.....ONE  if GPP = 0
   po0 = sif_kp0 / (sif_kf0 + sift%sif_kd + sif_kp0)
   ps = po0 * sift%sif_jejo
   sift%sif_x = 1.0 - ps / po0
   sift%sif_x = MIN(MAX(sift%sif_x,0.0), 1.0)

   !...Empirical part of formulation and may be PFT-dependent
   a3 = ((1.0 + sif_a2)*sift%sif_x**sif_a1)/(sif_a2 + sift%sif_x**sif_a1)
   sift%sif_kn = sif_kn0 * a3

   ! Photosynthesis probability (backed out using vdt14, Eq. 2)
   sift%sif_kp = (ps * &
        (sift%sif_kn + sift%sif_kd + sif_kf0))/ &
        (1.0 - ps)

   !----Yields
   totk = sift%sif_kd + sif_kf0 + sift%sif_kn + sift%sif_kp

   !...dissipation
   sift%phi_d = sift%sif_kd / totk

   !...SIF
   fo = sif_kf0 / totk
   fm = sif_kf0 / (sift%sif_kd + sif_kf0 + sift%sif_kn)
   sift%phi_f = fm * (1.0 - ps)

   !...NPQ
   sift%phi_n = sift%sif_kn / totk

   !...photosynthesis
   sift%phi_p = sift%sif_kp / totk

   !...sanity check
   toty = sift%phi_d +  sif_kf0/totk + sift%phi_n + sift%phi_p
   if ((toty .gt. 1.01) .or. (toty .lt. 0.99)) then
       print*,'SIF Yields Do Not Add To 1.0!!'
       stop
    endif

   !----SIF
   !...Follows J.-E. Lee et al. (2015)
   kappa = 0.04 * vmax * mol_to_umol + 8.1
   sift%sif = (nspar * mol_to_umol * fpar) * sift%phi_f / kappa

endif

end subroutine sif_calc
