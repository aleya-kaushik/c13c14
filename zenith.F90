!-------------------------------------------------------------------------------
subroutine init_solar_dec( doy, cureqnx, curday_per_year, lonearth )
!-------------------------------------------------------------------------------
! Initializes the declination of the Sun
!
! FUNCTIONS CALLED:  rt_asc
!-----------------------------------------------------------------------

use kinds
use module_pparams, only: pidaypy

implicit none

!...input variables
integer(i4), intent(in) :: doy
integer(i4), intent(in) :: cureqnx
integer(i4), intent(in) :: curday_per_year
real(r8), intent(inout) :: lonearth

! local variables
real(r8) :: t1       ! 1st factor to determine longitude of earth (lonearth)
real(r8) :: t2       ! 2nd factor to determine longitude of earth (lonearth)
real(r8) :: t3       ! 3rd factor to determine longitude of earth (lonearth)
real(r8) :: t4       ! 4th factor to determine longitude of earth (lonearth)
integer(i4) :: iday    ! day of year since vernal equinox variable
real(r8) :: rt_asc

    ! lon of Earth from equinox at start of simulation
    iday = doy - cureqnx
    if ( iday < 0 ) iday = iday + curday_per_year
    lonearth=0.0
    if ( iday /= 0 ) then
      do while ( iday > 1 )
        iday = iday - 1
        t1 = rt_asc( lonearth ) * pidaypy
        t2 = rt_asc( lonearth+t1*.5 ) * pidaypy
        t3 = rt_asc( lonearth+t2*.5 ) * pidaypy
        t4 = rt_asc( lonearth+t3 ) * pidaypy
        lonearth = lonearth + (t1+2.*(t2+t3)+t4) / 6.
      enddo
    endif

end subroutine init_solar_dec


subroutine solar_dec( doy, cureqnx, lonearth, &
                      sin_dec, cos_dec, tan_dec  )
    !-------------------------------------------------------------------
    ! Calculates the declination of the Sun
    !-------------------------------------------------------------------

    use kinds
    use module_pparams, only: &
          pi180, pidaypy, decmax

    implicit none

    !...input variables
    integer(i4), intent(in) :: doy
    integer(i4), intent(in) :: cureqnx

    real(r8), intent(inout) :: lonearth
    real(r8), intent(inout) :: sin_dec, cos_dec, tan_dec

    ! local variables
    real(r8) :: t1  ! 1st factor to determine longitude of earth (lonearth)
    real(r8) :: t2  ! 2nd factor to determine longitude of earth (lonearth)
    real(r8) :: t3  ! 3rd factor to determine longitude of earth (lonearth)
    real(r8) :: t4  ! 4th factor to determine longitude of earth (lonearth)
    real(r8) :: rt_asc

    ! reset lon of Earth from equinox
    if ( doy == cureqnx ) lonearth = 0.0

    ! Increment Longitude of Earth
    t1 = rt_asc( lonearth ) * pidaypy
    t2 = rt_asc( lonearth+t1*.5 ) * pidaypy
    t3 = rt_asc( lonearth+t2*.5 ) * pidaypy
    t4 = rt_asc( lonearth+t3 ) * pidaypy
    lonearth = lonearth + (t1+2.*(t2+t3)+t4) / 6.

    ! Calculate the sine and cosine of Solar declination
    sin_dec = sin( decmax*pi180 ) * sin( lonearth )
    cos_dec = sqrt( 1. - sin_dec * sin_dec )
    if (cos_dec .ne. 0.) then
        tan_dec = sin_dec / cos_dec
    else
        tan_dec = 0.
    endif

end subroutine solar_dec

!=======================================================================
subroutine zenith_angle ( latsib, cos_dec, sin_dec, &
                           loctime, cosz )
!=======================================================================      
! calculates the zenith angle 
!----------------------------------------------------------------------

use kinds
use module_pparams, only: pi180

implicit none

!...input variables
real(r4), intent(in) :: latsib
real(r8), intent(in) :: cos_dec, sin_dec
real(r8), intent(in) :: loctime
real(r8), intent(inout) :: cosz

! internal variables
real(r8) :: sinlat      ! sine of latitude
real(r8) :: coslat      ! cosine of latitude
real(r8) :: hrang       ! hour angle; longitude of Sun from Greenwhich meridian
real(r8) :: cos_hour    ! cosine delta longitude between Sun & SiB point

    ! Calculate hour angle (longitude of Sun from Greenwhich meridian)
    hrang = 360.*loctime-180.

    ! Calculate cosine of solar zenith angle
    cos_hour = cos( pi180 * hrang )
    sinlat   = sin( pi180 * latsib )
    coslat   = cos( pi180 * latsib )
    cosz = coslat * cos_dec * cos_hour + sinlat * sin_dec

end subroutine zenith_angle


!-------------------------------------------------------------------------------
function rt_asc( ang )
!-------------------------------------------------------------------------------
! calculates correction for longitude (right ascension) of the Earth
! from vernal equinox based on the angle around Sun traversed
! since beginning of the year

    use kinds
    use module_pparams, only: &
        eccn, pi, perhl

    implicit none

real(r8) :: rt_asc ! right ascension correction
real(r8) :: ang    ! angle from beginning of year (radians)

! local variables
real(r8) :: reccn
real(r8) :: perhlr

    ! rt ascension correction based on Earth's orbit
    reccn  = 1. / (1. - eccn * eccn ) **1.5
    perhlr = perhl * (pi/180.)
    rt_asc = reccn * (1. - eccn * cos( ang - perhlr ) ) **2

end function rt_asc

!-------------------------------------------------------------------------------
subroutine day_length(subcount, pi180, tan_dec, lat, &
                      dlength, dlengthdt)
!-------------------------------------------------------------------------------
! calculates the length of day

use kinds

implicit none

!...input variables
integer(i4), intent(in) :: subcount !number of points
real(r8), intent(in) :: pi180 !conversion from degrees to radians
real(r8), intent(in) :: tan_dec  !tangent of solar declination angle (-)
real(r4), dimension(subcount), intent(in) :: lat   !latitude (degrees)
real(r8), dimension(subcount), intent(inout) :: &
     dlength   !length of daylight (hr)
real(r8), dimension(subcount), intent(inout) :: &
     dlengthdt !change in daylight length (hr)

!...local variables
integer(i4) :: i
real(r8) :: tan_lat, cos_hr, pdlength

do i=1,subcount
   if ((lat(i) .gt. -90.) .and. (lat(i) .lt. 90.)) then
        tan_lat = tan(pi180*lat(i))
   else
        tan_lat = 0.
   endif

   cos_hr = -tan_lat * tan_dec
   if (cos_hr .lt. -1.0) then
      dlength(i) = 24.
      dlengthdt(i) = 0.
   elseif (cos_hr .gt. 1.0) then
      dlength(i) = 0.
      dlengthdt(i) = 0.
   else
      pdlength = dlength(i)
      dlength(i) = (acos(cos_hr)/pi180*24.)/180.
      dlengthdt(i) = dlength(i) - pdlength
   endif
enddo

end subroutine day_length

!-------------------------------------------------------------------------------
subroutine day_lengthpt(pi180, tan_dec, lat, &
                      dlength, dlengthdt)
!-------------------------------------------------------------------------------
! calculates the length of day

use kinds

implicit none

!...input variables
real(r8), intent(in) :: pi180 !conversion from degrees to radians
real(r8), intent(in) :: tan_dec  !tangent of solar declination angle (-)
real(r4), intent(in) :: lat !latitude (degrees)
real(r8), intent(inout) :: dlength   !length of daylight (hr)
real(r8), intent(inout) :: dlengthdt !change in daylight length (hr)

!...local variables
real(r8) :: tan_lat, cos_hr, pdlength

if ((lat .gt. -90.) .and. (lat .lt. 90.)) then
     tan_lat = tan(pi180*lat)
else
     tan_lat = 0.
endif

cos_hr = -tan_lat * tan_dec
if (cos_hr .lt. -1.0) then
   dlength = 24.
   dlengthdt = 0.
elseif (cos_hr .gt. 1.0) then
   dlength = 0.
   dlengthdt = 0.
else
   pdlength = dlength
   dlength = (acos(cos_hr)/pi180*24.)/180.
   dlengthdt = dlength - pdlength
endif


end subroutine day_lengthpt

!-----------------------------------------------------
subroutine max_day_length(subcount, lat, max_dlength) 
!-----------------------------------------------------
! calculates the minimum and maximum day length of day

use kinds
use module_pparams, only: &
   pi180, days_per_year, &
   eqnx, nhsolstice, shsolstice

implicit none

!...input variables
integer(i4), intent(in) :: subcount !number of points
real(r4), dimension(subcount), intent(in) :: &
     lat   !latitude (degrees)

!real(r4), dimension(subcount), intent(inout) :: &
!    min_dlength   !minimum day length (hr)

real(r4), dimension(subcount), intent(inout) :: &
     max_dlength   !maximum day length (hr)


!...local variables
integer(i4) :: i, mydoymin, mydoymax
real(r8) :: lonearth, cos_dec, sin_dec, tan_dec
real(r8) :: tan_lat, cos_hr

!...initialize local variables:
lonearth = 0.
tan_dec = 0.

do i=1,subcount
   if ((lat(i) .gt. -90.) .and. (lat(i) .lt. 90.)) then
        tan_lat = tan(pi180*lat(i))
   else
        tan_lat = 0.
   endif

   if (lat(i) .gt. 0.) then
      mydoymax = nhsolstice
      mydoymin = shsolstice
   else
      mydoymax = shsolstice
      mydoymin = nhsolstice
   endif
   
   !calculate daylength of longest day of year
   call init_solar_dec( mydoymax, eqnx, &
                  days_per_year, lonearth)
   call solar_dec(mydoymax, eqnx, lonearth, &
                  sin_dec, cos_dec, tan_dec  )

   cos_hr = -tan_lat * tan_dec
   if (cos_hr .lt. -1.0) then
      max_dlength(i) = 24.
   elseif (cos_hr .gt. 1.0) then
      max_dlength(i) = 0.
   else
      max_dlength(i) = real((acos(cos_hr)/pi180*24.)/180.)
   endif

   !calculate daylength of shortest day of year
   !call init_solar_dec( mydoymin, eqnx, &
   !               days_per_year, lonearth)
   !call solar_dec(mydoymin, eqnx, lonearth, &
   !               sin_dec, cos_dec, tan_dec  )

   !cos_hr = -tan_lat * tan_dec
   !if (cos_hr .lt. -1.0) then
   !   min_dlength(i) = 24.
   !elseif (cos_hr .gt. 1.0) then
   !   min_dlength(i) = 0.
   !else
   !   min_dlength(i) = real((acos(cos_hr)/pi180*24.)/180.)
   !endif

enddo


end subroutine max_day_length


!-------------------------------------------------------------------------------
subroutine calc_lst(dayfrac, &
                    subcount, lonsib, &
                    dayfrac_lst)
!-------------------------------------------------------------------------------
!
! calculates the local standard time hour fraction
!   set to the nearest hour
!

use module_sibconst, only: updatelst_switch
use kinds
implicit none

!...input variables
real(r8), intent(in) :: dayfrac
integer(i4), intent(in) :: subcount
real(r4), dimension(subcount), intent(in) :: lonsib
real(r8), dimension(subcount), intent(inout) :: dayfrac_lst
real(r8), dimension(25) :: dayfrac_hr

!...local variables
integer(i4) :: hr, pt
integer(i4) :: nearhr
real(r8) :: lstfrac

dayfrac_hr(:) = 0
do hr=0,23
   dayfrac_hr(hr+1) = hr / dble(24.)
enddo

do pt=1, subcount
   if (updatelst_switch) then
       lstfrac = dble(lonsib(pt)) / dble(360.)
       if (lstfrac > 1.0) lstfrac = lstfrac - 1.0
       if (lstfrac < 0.0) lstfrac = lstfrac + 1.0

       nearhr = nint(lstfrac*24)
       lstfrac = dayfrac_hr(nearhr+1)
   else
       lstfrac = 0.
   endif

   dayfrac_lst(pt) = dayfrac + lstfrac
enddo

end subroutine calc_lst

!-------------------------------------------------------------------------------
subroutine update_lst(subcount, dtdayfrac, dayfrac_lst, &
                       new_day_lst)
!-------------------------------------------------------------------------------
!
! updates the local standard time hour fraction
!

use kinds
implicit none


!...input variables
integer(i4), intent(in) :: subcount
real(r8), intent(in) :: dtdayfrac
real(r8), dimension(subcount), intent(inout) :: dayfrac_lst
logical, dimension(subcount), intent(inout) :: new_day_lst

where(new_day_lst)
   new_day_lst = .false.
endwhere

dayfrac_lst(:) = dayfrac_lst(:) + dtdayfrac
where(dayfrac_lst .gt. 0.9999) 
    new_day_lst = .true.
    dayfrac_lst = dzero
end where

end subroutine update_lst
