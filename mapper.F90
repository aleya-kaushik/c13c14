subroutine aerointerpolate(gref, glon, glat, &
    pnum, pref, lai, fvcover, &
    z0, zp_disp, rbc, rdc)

    !----------------------------------------------------------------
    ! This subroutine calculates the aerodynamic properties by
    ! bi-linear interpolation from a lookup table of previously
    ! calculated values.
    ! The interporation table is a numpts x numpts LAI grid
    ! with LAI ranging from 0.02 to 10.
    !----------------------------------------------------------------

    use kinds
    use module_param, only: &
         ngrid, aerovar, &
         fVCovergrid, LAIgrid

    IMPLICIT NONE

    ! parameters
    logical, parameter :: high_lai = .true.
    
    ! input variables
    integer(i4), intent(in) :: gref, pnum, pref
    real(r4), intent(in) :: glon, glat
    real(r8), intent(in) :: lai, fvcover
    real(r8), intent(inout) :: z0, zp_disp, rbc, rdc

    ! locals
    integer(i4)  :: i      ! index for LAI grid
    integer(i4)  :: j      ! index for fVCover grid

    real(r8) :: loclai     ! local LAI
    real(r8) :: dlai       ! grid spacing between LAI values
    real(r8) :: laiwt1, laiwt2     ! linear weights for interpolation between LAI

    real(r8) :: locfvcover ! local fVCover
    real(r8) :: dfvcover   ! grid spacing between fVCover values
    real(r8) :: fvcovwt1, fvcovwt2 ! linear weights for interpolation between fVCover

    ! Calculate grid spacing
    dfvcover = fvcovergrid(2) - fvcovergrid(1)
    dlai = laigrid(2) - laigrid(1)

    ! Assign input LAI to local variables and make sure 
    ! they lie within the limits of the interpolation tables.
    locfvcover = max(fvcover, 0.01_r8)
    loclai = max(lai, 0.02_r8)

    ! Determine the nearest array location for the desired LAI
    i = int(loclai / dlai) + 1
    if (i+1 > ngrid) then
       if (high_lai) THEN
          i = ngrid-1
       else
          print*,'LAI higher than interpolation tables: ',loclai
          print'(a,2F10.2,I8,I4)','Point (lon/lat/sibpt/pft): ',glon,glat,gref,pref
          stop
       endif
    endif
    laiwt1 = (loclai - laigrid(i)) / dlai
    laiwt2 = (laigrid(i+1) - loclai) / dlai

    ! Determine the nearest array location for the desired fVCover
    j = int(locfvcover / dfvcover) + 1

    if (j+1 > ngrid) then
       if (high_lai) then
          j = ngrid-1
       else
          print*,'fVCover higher than interpolation tables: ',locfvcover
          print'(a,2F10.2,I8,I4)','Point (lon/lat/sibpt/pft): ',glon,glat,gref,pref
          stop
       endif
    endif
    fvcovwt1 = (locfvcover - fvcovergrid(j)) / dfvcover
    fvcovwt2 = (fvcovergrid(j+1) - locfvcover) / dfvcover

     ! Linearly interpolate RbC, RdC, roughness length, zero plane displacement
     rbc = (laiwt2*aerovar(pnum,i,j)%rbc + laiwt1*aerovar(pnum,i+1,j)%rbc)*0.5 + &
           (fvcovwt2*aerovar(pnum,i,j)%rbc + fvcovwt1*aerovar(pnum,i,j+1)%rbc)*0.5
     rdc = (laiwt2*aerovar(pnum,i,j)%rdc + laiwt1*aerovar(pnum,i+1,j)%rdc)*0.5 + &
           (fvcovwt2*aerovar(pnum,i,j)%rdc + fvcovwt1*aerovar(pnum,i,j+1)%rdc)*0.5
     z0  = (laiwt2*aerovar(pnum,i,j)%zo  + laiwt1*aerovar(pnum,i+1,j)%zo)*0.5 + &
           (fvcovwt2*aerovar(pnum,i,j)%zo + fvcovwt1*aerovar(pnum,i,j+1)%zo)*0.5
     zp_disp = (laiwt2*aerovar(pnum,i,j)%zp_disp + laiwt1*aerovar(pnum,i+1,j)%zp_disp)*0.5 + &
               (fvcovwt2*aerovar(pnum,i,j)%zp_disp + fvcovwt1*aerovar(pnum,i,j+1)%zp_disp)*0.5

end subroutine aerointerpolate


! calculates time mean leaf projection relative to the Sun.
subroutine gmuder(Lat, DOY, ChiL, gmudmu)
    use kinds
    use module_pparams, only: pi, pi180

    implicit none

    ! input variables
    real(r4), intent(in)    :: Lat  ! latitude in degrees
    integer(i4), intent(in) :: DOY  ! day of year
    real(r4), intent(in)  :: ChiL   ! leaf angle distribution factor
    real(r8), intent(out) :: gmudmu ! time mean projected leaf area normal to Sun

    ! locals
    integer :: itime ! time counter
    real(r8) :: gtime    ! time from 0:00 Greenwhich Mean Time (GMT)
    real(r8) :: coshr    ! cosine of the Greenwhich Meridian (GM) Hour angle
    real(r8) :: mu       ! cosine of the Solar zenith angle
    real(r8) :: chiv     ! dummy variable for leaf angle distribution factor
    real(r8) :: dec      ! declination of the Sun (Solar Declination)
    real(r8) :: sin_dec  ! sine of the solar declination
    real(r8) :: cos_dec  ! cosine of the solar declination
    real(r8) :: aa       ! minimum possible LAI projection vs. cosine Sun angle
    real(r8) :: bb       ! slope leaf area projection vs. cosine Sun angle
    real(r8) :: cloud    ! (?) Cloud cover fraction
    real(r8) :: fb       ! (?) mean cosine of Sun Angle
    real(r8) :: swdown   ! (?) downward shortwave flux
    real(r8) :: pardif   ! (?) PAR flux difracted into canopy by clouds
    real(r8) :: pardir   ! (?) PAR flux directly onto canopy
    real(r8) :: difrat   ! (?) fraction of shortwave flux defracted by clouds
    real(r8) :: vnrat    ! (?) shortwave flux ratio of some sort
    real(r8) :: tor      ! (?) TBD
    real(r8) :: topint   ! (?) Total flux onto canopy adjusted for sun angle
    real(r8) :: botint   ! (?) total PAR flux onto canopy during 24 hour period

    cloud=0.5

    ! Calculate solar declination in radians
    dec=pi180*23.5*sin(1.72e-2*(DOY-80))

    ! Calculate sine and cosine of solar declination
    sin_dec=sin(dec)                                                         
    cos_dec=cos(dec)

    ! Begin time loop to average leaf projection over 24 hour period
    topint=0.
    botint=0.

    do itime=1, 48, 1
        ! Calculate time from zero Greenwhich Mean Time (GMT)
        gtime=0.5*real(itime) 

        ! Calculate cosine of hour angle of Grenwhich Meridion (GM)
        coshr=cos(-pi+gtime/24.*2.*pi)

        ! Calculate cosine of the Sun angle (mu)
        !     longitude=GM=0 degrees, latitude=Lat
        mu=sin(Lat*pi180)*sin_dec+cos(Lat*pi180)*cos_dec*coshr

        ! Ensure the cosine of Sun angle is positive, but not zero
        !     e.g., daylight, sun angle<=89.4 degrees (about start disc set/rise)
        mu=dble(amax1(0.01, real(mu))) 

        ! It looks like this code calculates the direct and difracted PAR based
        ! on the solar constant and a cloud fraction of 0.5 at the top and
        ! bottom of the atmosphere.  During night, mu=0.01, a constant.  These
        ! calculations do not match the definition of G(mu)/mu described in 
        ! Bonan (1996) and Sellers (1985).
        tor    = 0.7**(1./mu)
        swdown = 1375.*mu*(tor+0.271-0.294*tor)
        difrat = 0.0604/(mu-0.0223)+0.0683
        difrat = max(difrat,0.)
        difrat = min(difrat,1.)
        difrat = difrat+(1.-difrat)*cloud
        vnrat  = (580.-cloud*464.)/((580.-cloud*499.) + (580.-cloud*464.))
        pardir = (1.-difrat)*vnrat*swdown
        pardif = difrat*vnrat*swdown
        topint = topint+pardir*mu+pardif*0.5
        botint = botint+pardir+pardif
    enddo

    ! Calculate what looks like a mean value of something
    fb=topint/botint

    ! Calculate min and slope of LAI projection 
    chiv=ChiL                                                               
    if (abs(chiv) <= 0.01) chiv=0.01
    !   calculate minimum value of projected leaf area
    aa=0.5-0.633*chiv-0.33*chiv*chiv
    !   calculate slope of projected leaf area wrt cosine sun angle
    bb=0.877*(1.-2.*aa)                                             

    ! Calculate mean projected LAI in Sun direction assuming fb approximates
    ! the mean cosine of the sun angle
    gmudmu=(aa+bb*fb)/fb

end subroutine gmuder

