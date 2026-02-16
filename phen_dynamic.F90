!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Dynamic Phenology Stages
!   Number Of Stages = 5
!    1: Leaf-Out
!    2: Growth
!    3: Mature
!    4: Stress 
!    5: Senescence/Dormant
!
!   Stages are determined by calculating the 
!   phenology stage index by combining:
!    - Growth Potential
!    - Day Length Potential
!    - Weather Potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine phen_dynamic(phencont, &
     dlenmax, dlen, dlendt, lai, &
     phent)

use kinds
use module_param, only: phen_param
use module_sib, only: phen_type
use module_time, only: wt_seas

implicit none

!...input variables
type(phen_param), intent(in) :: phencont
real(r4), intent(in) :: dlenmax
real(r8), intent(in) :: dlen, dlendt, lai
type(phen_type), intent(inout) :: phent

!...local variables
real(r8) :: daylenadd, daylenoff, slope, yint

!Calculate phenology stage determinants
!...Day Length Potential
phent%phens_dayl = done
IF (phencont%psdayl_ref .ge. rzero) THEN
   daylenoff = dlenmax - phencont%psdayl_ref
   if ((dlendt .le. rzero) .and. (dlen .le. daylenoff)) then
      phent%phens_dayl = rone - (daylenoff - dlen)*phencont%psdayl_mul
   endif
ELSE
   daylenoff = dlenmax + phencont%psdayl_ref
   if (dlendt .gt. rzero) then
      if (dlen .gt. daylenoff) then
         phent%phens_dayl = rone - (dlen - daylenoff)*phencont%psdayl_mul
      endif
   else
      daylenadd = (dlenmax - daylenoff)*phencont%psdayl_mul
      phent%phens_dayl = rone - daylenadd - (dlenmax - dlen)*phencont%psdayl_mul
   endif
ENDIF
phent%phens_dayl = MIN(done, MAX(phencont%psdayl_min, &
                    phent%phens_dayl))

!...Growth Potential
IF (lai .lt. phent%phenc_laimin) THEN
   phent%phens_grw = done
ELSEIF (lai .lt. phent%phenc_laimax) THEN
   slope = (done - phencont%psg_min) &
         / (phent%phenc_laimin - phent%phenc_laimax)
   yint = phencont%psg_min - slope*phent%phenc_laimax
   phent%phens_grw = (1.-wt_seas)*(slope*lai + yint) &
                   + wt_seas*phent%phens_grw
ELSE
   phent%phens_grw = phencont%psg_min
ENDIF

!...Weather Potential
IF (phent%phenave_wacsm .gt. dzero) THEN
   phent%phens_wx = phent%phenave_wac &
            / phent%phenave_wacsm
ELSE
   phent%phens_wx = done
ENDIF


!...Combined Index
phent%phen_pi = MIN(phent%phens_dayl,phent%phens_grw,phent%phens_wx)


!...Set Phenology Stage
IF (phent%phen_pi .gt. phencont%threshp(1)) THEN
   phent%phen_istage = 1
ELSEIF (phent%phen_pi .gt. phencont%threshp(2)) THEN
   phent%phen_istage = 2
ELSEIF (phent%phen_pi .gt. phencont%threshp(3)) THEN
   phent%phen_istage = 3
ELSEIF (phent%phen_pi .gt. phencont%threshp(4)) THEN
   phent%phen_istage = 4
ELSE
   phent%phen_istage = 5
ENDIF

end subroutine phen_dynamic
