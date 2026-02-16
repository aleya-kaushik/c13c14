subroutine phen_gss(phencont, daylenmax, daylen, daylendt, phent)

! Determines if meteorological/climatological conditions are met
! for the growing season to start

use kinds
use module_param, only: phen_param
use module_sib, only: phen_type
use module_time, only: doy, curday_per_year

implicit none

!...parameters
real(r8), parameter :: cyears = 10. !mean number of years
real(r8), parameter :: wt_climyr = 1./cyears !calculation weight

! Input
type(phen_param), intent(in) :: phencont
real(r4), intent(in) :: daylenmax
real(r8), intent(in) :: daylen, daylendt
type(phen_type), intent(inout) :: phent

!------------------------------
!Growing Season Start Flag
!...Daylength Trigger
phent%phenflag_daylen = .false.
IF (daylendt .GE. 0.) THEN
   IF (daylen .GT. phencont%daylen_mini) phent%phenflag_daylen = .true.
ELSE
   IF (phencont%daylen_offd .GT. rzero) THEN
      IF ((daylenmax - daylen) .LT. phencont%daylen_offd) &
           phent%phenflag_daylen = .true.
   ELSE
      IF (daylen .LT. abs(phencont%daylen_offd)) &
           phent%phenflag_daylen = .true.
   ENDIF
ENDIF

!...Soil Moisture Trigger
IF (phent%phenave_tawftop .GT. phencont%tawftop_min) THEN
    phent%phenflag_moist=.true.
ELSE
    phent%phenflag_moist=.false.
ENDIF

!...Temperature Trigger
phent%phenflag_temp=.false.
IF ((phent%phenave_tm .GT. dble(phencont%tm_min)) .AND. &
    (phent%phenave_tm .LT. dble(phencont%tm_max))) THEN
   phent%phenflag_temp=.true.
ENDIF

!...Combined Flag
phent%phenflag_gsspass = .false.
IF (phent%phenflag_daylen .AND. &
    phent%phenflag_moist .AND. &
    phent%phenflag_temp) THEN
   phent%phenflag_gsspass = .true.
ENDIF


!-------------------------
!Growing Season Reset Flag
IF (phent%phenave_assimsm .GT. dzero) THEN
    phent%phenave_assimpot = phent%phenave_assim / phent%phenave_assimsm
ELSE
    phent%phenave_assimpot = dzero
ENDIF
IF (phent%phenave_assimpot .lt. phencont%assim_resetv) THEN
   phent%phenflag_assimlow = .true.
ENDIF

!--------------------------
!Precipitation Window Flag
IF ((phent%phenave_pr .GT. phent%phenave_prsm) .AND. &
    (phent%phenflag_gsspass)) THEN
    phent%phenave_prsm = phent%phenave_pr
    phent%phenave_prsdoy = doy
ENDIF

IF ((doy .LT. phent%phenave_prcdoy + phencont%precip_aft) .AND. &
    (doy .GT. phent%phenave_prcdoy - phencont%precip_bef)) THEN
     phent%phenflag_precip = .true.
ELSE
     phent%phenflag_precip = .false.
ENDIF
     
IF (doy .EQ. curday_per_year) THEN
   IF (phent%phenave_prcdoy .ge. done) THEN
      phent%phenave_prcdoy = phent%phenave_prsdoy*wt_climyr &
           + (1.-wt_climyr)*phent%phenave_prcdoy
   ELSE
      phent%phenave_prcdoy = phent%phenave_prsdoy
   ENDIF

   phent%phenave_prsm = dzero
   phent%phenave_prsdoy = dzero
ENDIF

!---------------------------------
!Reset growing-season variables
IF ((phent%phenflag_assimlow) .and. &
    (phent%phenflag_gsspass)) THEN
    phent%phenave_wacsm = phent%phenave_wac
    phent%phenave_assimsm = phent%phenave_assim
    phent%phenflag_assimlow = .false.

    phent%nd_dormant = izero
    phent%nd_gs = izero
    phent%nd_stg(:) = izero
ENDIF

end subroutine
