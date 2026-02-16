!--------------------------------------------------------------
subroutine driver_interp(indx, lon, lat, gdiagt, gprogt)
!--------------------------------------------------------------
!
! This subroutine interpolates the sibdrv forcing surface meteorological
! variables between their read times
!
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
real(r8) :: localtime_driverstep1, localtime_driverstep2  ! local time of prev/next driver step
real(r8) :: cosz_old, cosz_driverstep1, cosz_driverstep2  ! cosz of prev/next driver step
real(r8) :: facsibdrv  ! scaling factor between driver data points
real(r8) :: facswdwn   ! scaling factor based on cosz

    ! reset variables
    gdiagt%sif_flag(:) = .false.

    ! calculate cosine zenith angle
    localtime = dble(sec_day) / dble(secs_per_day) + (dble(lon) / dble(360.0))
    call zenith_angle( lat, cos_dec, sin_dec, localtime, gdiagt%cosz )

    ! get scaling factors
    facsibdrv = MAX(dzero, dble(driver_seccur-sec_tot) / dble(driver_step))

    cosz_old = gdiagt%cosz
    localtime_driverstep1 = dble(sec_day-driver_step*(1.0-facsibdrv)) / dble(secs_per_day) + &
                               (dble(lon) / dble(360.0))
    if (localtime_driverstep1 > done) localtime_driverstep1 = localtime_driverstep1 - 1.0
    if (localtime_driverstep1 < dzero) localtime_driverstep1 = localtime_driverstep1 + 1.0
    call zenith_angle( lat, cos_dec, sin_dec, localtime_driverstep1, cosz_driverstep1 )

    if (localtime > done) localtime = localtime - 1.0
    if (localtime < dzero) localtime = localtime + 1.0

    localtime_driverstep2 = dble(sec_day+driver_step*facsibdrv) / dble(secs_per_day) + &
                               (dble(lon) / dble(360.0))
    if (localtime_driverstep2 > done) localtime_driverstep2 = localtime_driverstep2 - 1.0
    if (localtime_driverstep2 < dzero) localtime_driverstep2 = localtime_driverstep2 + 1.0
    call zenith_angle( lat, cos_dec, sin_dec, localtime_driverstep2, cosz_driverstep2 )

    if (abs(cosz_driverstep1-cosz_driverstep2) .gt. 1.E-12) then
       facswdwn = abs(gdiagt%cosz-cosz_driverstep2)/abs(cosz_driverstep1-cosz_driverstep2)
    else
       facswdwn = facsibdrv
    endif

    if (((gdiagt%cosz .gt. cosz_driverstep1) .and. (gdiagt%cosz .gt. cosz_driverstep2)) .or. &
        ((gdiagt%cosz .lt. cosz_driverstep1) .and. (gdiagt%cosz .lt. cosz_driverstep2)) .or. &
        ((facswdwn .lt. 0.2) .and. (facsibdrv .gt. 0.8))) then
        facswdwn=facsibdrv
    endif

    !.....TOA calculations
    ! calculate TOA solar downwelling
    if (gdiagt%cosz > 0.0) then
        coslon = cos( pi180 * lon)
        dbarod = 1.0 + eccn * (coslon - perih)
        gdiagt%toa_solar = solar_const * dbarod * dbarod * gdiagt%cosz
    else
        gdiagt%toa_solar = dzero
    endif
    ! split downwelling TOA radiation to visible and NIR (direct and diffuse)
     call raddrv(gdiagt%toa_solar,gdiagt%cosz,gdiagt%toa_radvbc,  &
                 gdiagt%toa_radvdc,gdiagt%toa_radnbc,gdiagt%toa_radndc)
    ! calculate TOA PAR
    gdiagt%toa_par = par_conv * (gdiagt%toa_radvbc + gdiagt%toa_radvdc)

    !...SW radiation calculations
    ! scale SW radiation using cosine of zenith angle
    if (gdiagt%cosz > dzero) then
       gprogt%sw_dwn = gprogt%sw_dwn1*facswdwn + gprogt%sw_dwn2*(1-facswdwn)
    else
       gprogt%sw_dwn = dzero
    endif
    if (gprogt%sw_dwn < dzero) gprogt%sw_dwn = dzero

    ! split downwelling SW radiation to visible and NIR (direct and diffuse)
     call raddrv(gprogt%sw_dwn,gdiagt%cosz,gdiagt%radvbc,  &
                 gdiagt%radvdc,gdiagt%radnbc,gdiagt%radndc)

    !...SIF calculations
    ! calculate aerosol + cloud optical depth
    !     based on surface vs TOA insolation
    if ((gdiagt%toa_solar > 20.0) .AND. &
        ((gdiagt%radvbc + gdiagt%radnbc) > 0.0)) then
         gdiagt%aod = mu_1 * LOG(ABS((gdiagt%radvbc + gdiagt%radnbc) &
                      / gdiagt%toa_solar))

         ! determine SIF attenuation
         !    based on Frankenberg et al. (2012; Figure 12):
         !    (AMT 5, 2081-2094 doi:10.5194/amt-5-2081-2012)
         gdiagt%sif_atten = MIN(done, MAX(dzero, &
                 EXP((-0.05 / gdiagt%cosz) * gdiagt%aod)))
    else
         gdiagt%aod = dzero
         gdiagt%sif_atten = dzero
    endif

    ! determine if conditions are suitable to take SIF satellite measurements
    if (qp_savegf .OR. pbp_savegf .OR. hr_savegf) then
       if (gdiagt%sif_atten > 0.75) then
          if ((localtime >= sif_gome2_tstart) .AND. &
              (localtime < sif_gome2_tstop)) then !GOME-2
              gdiagt%sif_flag(1) = .true.
              IF (qp_savegf) &
                  qp_sif_satcount(indx,qp_countsave,1) = &
                     qp_sif_satcount(indx,qp_countsave,1) + 1
              IF (hr_savegf) &
                  hr_sif_satcount(indx,hr_countsave,1) = &
                      hr_sif_satcount(indx,hr_countsave,1) + 1
              IF (pbp_savegf) &
                  pbp_sif_satcount(indx,pbp_countsave,1) = &
                      pbp_sif_satcount(indx,pbp_countsave,1) + 1
          endif

          if ((localtime >= sif_oco2_tstart) .AND. &
              (localtime < sif_oco2_tstop)) then !GOSAT/OCO-2
              gdiagt%sif_flag(2) = .true.
              IF (qp_savegf) &
                  qp_sif_satcount(indx,qp_countsave,2) = &
                     qp_sif_satcount(indx,qp_countsave,2) + 1
              IF (hr_savegf) &
                   hr_sif_satcount(indx,hr_countsave,2) = &
                      hr_sif_satcount(indx,hr_countsave,2) + 1
              IF (pbp_savegf) &
                   pbp_sif_satcount(indx,pbp_countsave,2) = &
                      pbp_sif_satcount(indx,pbp_countsave,2) + 1
           endif
       endif
    endif

    !...Interpolate driver data
    ! calculate driver data interpolation factor
    facsibdrv = MAX(0., dble(driver_seccur-sec_tot) / dble(driver_step))

    ! interpolate driver temperature
    gprogt%tm = facsibdrv*gprogt%tm1 + (1.-facsibdrv) * gprogt%tm2

    !  set daily value
    gprogt%tmd = gprogt%tmd*(done-wt_daily) + gprogt%tm*wt_daily
    gdiagt%tmdf = ((gprogt%tmd - tref)*1.8) + 32.0
    
    ! interpolate driver humidity
    gprogt%sh = facsibdrv*gprogt%sh1 + (1.-facsibdrv) * gprogt%sh2

    ! interpolate driver surface pressure
    gprogt%ps = facsibdrv*gprogt%ps1 + (1.-facsibdrv) * gprogt%ps2

    ! interpolate driver wind speed
    gprogt%spdm = facsibdrv*gprogt%spdm1 + (1.-facsibdrv) * gprogt%spdm2
    gprogt%spdm = MAX(gprogt%spdm,1.0_r8)

    ! interpolate driver large scale precipitation
    !gprogt%lspr =  facsibdrv*gprogt%lspr1 + (1.-facsibdrv) * gprogt%lspr2
    gprogt%lspr = gprogt%lspr1             !(mm/s)
    gprogt%lsprt = gprogt%lspr * 0.001     !(m/s)

    ! interpolate driver cumulus or convective precipitation
    !gprogt%cupr =  facsibdrv*gprogt%cupr1 + (1.-facsibdrv) * gprogt%cupr2
    gprogt%cupr = gprogt%cupr1           !(mm/s)
    gprogt%cuprt = gprogt%cupr * 0.001   !(m/s)

    ! check for ridiculously low precipitation
    if ((gprogt%cuprt + gprogt%lsprt) < 1.0E-10) then
        gprogt%cuprt = dzero
        gprogt%lsprt = dzero
    endif

    ! interpolate driver longwave radiation
    gprogt%dlwbot = facsibdrv*gprogt%dlwbot1 + (1.-facsibdrv) * gprogt%dlwbot2

    ! calculate reference level potential temperature and 
    ! some related quantities
    gdiagt%bps(1) = (0.001*gprogt%ps)**kappa

    gdiagt%bps(2) = (0.001*(gprogt%ps-psb))**kappa

    gdiagt%thm = gprogt%tm / gdiagt%bps(1)
    gdiagt%em = gprogt%sh * gprogt%ps / (0.622 + gprogt%sh)

    ! calculate reference level air density
    gdiagt%ros = rgfac * gprogt%ps / gprogt%tm

    ! calculate psycrometric constant
    gdiagt%psy = spec_heat_cp / lvap * gprogt%ps / 0.622

    ! calculate climatological variables
    gprogt%seas_precip = (1-wt_seas)*gprogt%seas_precip &
         + wt_seas*(gprogt%cupr + gprogt%lspr) &
         * dtsib * steps_per_day !convert to mm/day
    gprogt%seas_tm = (1-wt_seas)*gprogt%seas_tm &
         + wt_seas*gprogt%tm
    
    gprogt%clim_cupr = (1 - wt_clim)*gprogt%clim_cupr &
         + wt_clim*gprogt%cupr*dtsib*steps_per_day !convert to mm/day
    gprogt%clim_precip = (1 - wt_clim)*gprogt%clim_precip &
         + wt_clim*(gprogt%cupr + gprogt%lspr) &
         * dtsib * steps_per_day !convert to mm/day
    gprogt%clim_tm = (1 - wt_clim)*gprogt%clim_tm &
         + wt_clim*gprogt%tm

end subroutine driver_interp

!---------------------------------------------------------------------
subroutine raddrv(swdown,sunang,radvbc,radvdc,radnbc,radndc)
!---------------------------------------------------------------------
!
!  Radiation code to use the downward short-wave radiation
!  and the sun angle to estimate radvbc,radvdc, radndc, radndc
!
!---------------------------------------------------------------------

use kinds
use module_oparams, only: &
    rad_c1, rad_c2, rad_c3, &
    rad_c4, rad_c5

implicit none

!...input variables
real(r8), intent(in) :: swdown, sunang
real(r8), intent(inout) :: radvbc, radvdc, radnbc, radndc

!...local variables
real(r8) ::  &
    cloud, difrat, vnrat
real(r8) ::  &
    stemp, localcosz

     localcosz = max( 0.001_r8, sunang )
     stemp = swdown
     stemp = MAX(stemp,0.01_r8)

     cloud = (rad_c5 * localcosz - stemp) / (rad_c4 * localcosz)                   
     cloud = max(cloud,0.)                                                
     cloud = min(cloud,1.)                                                  

     difrat = 0.0604 / ( sunang-0.0223 + 1.0e-10 ) + 0.0683
     if ( difrat < 0. ) difrat = 0.
     if ( difrat > 1. ) difrat = 1.
     difrat = difrat + ( 1. - difrat ) * cloud

     vnrat = ( rad_c1 - cloud*rad_c2 ) / ( ( rad_c1 - cloud*rad_c3 ) + ( rad_c1 - cloud*rad_c2 ) )

     radvbc = (1.-difrat)*vnrat*stemp
     radvdc = difrat*vnrat*stemp
     radnbc = (1.-difrat)*(1.-vnrat)*stemp
     radndc = difrat*(1.-vnrat)*stemp

end subroutine raddrv
