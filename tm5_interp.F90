!--------------------------------------------------------------
subroutine tm5_interp(indx, lon, lat, gdiagt, gprogt)
!--------------------------------------------------------------
!
! This subroutine interpolates the sibdrv forcing surface meteorological
! variables between their read times
!
! Written by Ara Cho (2020)
!--------------------------------------------------------------

use kinds
use module_io
use module_oparams, only: &
   mu_1, perih
use module_pparams
use module_sib, only: &
   gdiag_type, gprog_type
use module_sibconst, only:  &
   cos_dec, sin_dec
use module_time, only: &
   dtsib, sec_day, sec_tot, &
   steps_per_day, &
   wt_clim, wt_daily, wt_seas

implicit none

!...input variables
integer(i4), intent(in) :: indx
real(r4), intent(in) :: lon, lat
type (gdiag_type) :: gdiagt
type (gprog_type) :: gprogt

!...sif variables
real(r8) :: dbarod, coslon

!...local variables
real(r8) :: localtime  ! local time, fraction of day
real(r8) :: localtime_tm5step1, localtime_tm5step2  ! local time of prev/next tm5 step
real(r8) :: cosz_old, cosz_tm5step1, cosz_tm5step2  ! cosz of prev/next tm5 step
real(r8) :: facsibdrv  ! scaling factor between tm5 data points
real(r8) :: facswdwn   ! scaling factor based on cosz

    ! reset variables
    gdiagt%sif_flag(:) = .false.

    ! calculate cosine zenith angle
    localtime = dble(sec_day) / dble(secs_per_day) + (dble(lon) / dble(360.0))
    call zenith_angle( lat, cos_dec, sin_dec, localtime, gdiagt%cosz )

    ! get scaling factors
    facsibdrv = MAX(dzero, dble(tm5_seccur-sec_tot) / dble(tm5_step))

    cosz_old = gdiagt%cosz
    localtime_tm5step1 = dble(sec_day-tm5_step*(1.0-facsibdrv)) / dble(secs_per_day) + &
                               (dble(lon) / dble(360.0))
    if (localtime_tm5step1 > done) localtime_tm5step1 = localtime_tm5step1 - 1.0
    if (localtime_tm5step1 < dzero) localtime_tm5step1 = localtime_tm5step1 + 1.0
    call zenith_angle( lat, cos_dec, sin_dec, localtime_tm5step1, cosz_tm5step1 )

    if (localtime > done) localtime = localtime - 1.0
    if (localtime < dzero) localtime = localtime + 1.0

    if (localtime_tm5step2 > done) localtime_tm5step2 = localtime_tm5step2 - 1.0
    if (localtime_tm5step2 < dzero) localtime_tm5step2 = localtime_tm5step2 + 1.0
    call zenith_angle( lat, cos_dec, sin_dec, localtime_tm5step2, cosz_tm5step2 )

    if (abs(cosz_tm5step1-cosz_tm5step2) .gt. 1.E-12) then
       facswdwn = abs(gdiagt%cosz-cosz_tm5step2)/abs(cosz_tm5step1-cosz_tm5step2)
    else
       facswdwn = facsibdrv
    endif

    if (((gdiagt%cosz .gt. cosz_tm5step1) .and. (gdiagt%cosz .gt. cosz_tm5step2)) .or. &
        ((gdiagt%cosz .lt. cosz_tm5step1) .and. (gdiagt%cosz .lt. cosz_tm5step2)) .or. &
        ((facswdwn .lt. 0.2) .and. (facsibdrv .gt. 0.8))) then
        facswdwn=facsibdrv
    endif

    !...Interpolate tm5 data
    ! calculate tm5 data interpolation factor
    facsibdrv = MAX(0., dble(tm5_seccur-sec_tot) / dble(tm5_step))

    ! interpolate tm5 temperature
    !gprogt%pressure_tm5 = facsibdrv*gprogt%pressure_tm51 &
    ! + (1.-facsibdrv) * gprogt%pressure_tm52

    ! interpolate tm5 temperature
    gprogt%cosm_tm5 = facsibdrv*gprogt%cosm_tm51 &
     + (1.-facsibdrv) * gprogt%cosm_tm52
end subroutine tm5_interp
