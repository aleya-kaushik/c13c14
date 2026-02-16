!================SUBROUTINE CO2=========================
!
!    DERIVATION OF ATMOSPHERIC CO2 LEVEL
!    BASED ON PARAMETRIZATION MADE BY
!    I. VAN DER VELDE IN SIBCASA
!
!-------------------------------------------------------

subroutine set_co2( year, mon, gprogt )

    use module_pparams, only: &
        p0_sfc, bco2m
    use module_sib, only: &
        sib, gprog_type
    use module_sibconst, only: &
        subcount, varco2_switch
    !use module_phosib, only: &
    !    pressure 

    implicit none

    !---------------------------------------------------
    !...input variables
    integer :: year   ! year
    integer :: mon ! month of year

    type(gprog_type), intent(inout) :: gprogt

    !---------------------------------------------------
    !...local variables
    integer :: i, n    ! iteration 
    real, dimension(subcount) :: latarg ! matrix lat argument
    real, dimension(5) ::    flask_lat  ! flask latitudes
    real, dimension(5,12) :: flask_co2  ! flask seasonal cycle
    real :: seas ! seasonal variation adjustment
    integer :: num1 ! southern flask number index for interpolation
    integer :: num2 ! northern flask number index for interpolation

    !..Set CO2 global base value in ppmv
    !pressure = dble(gprogt%ps) * 100.0D0

    
    if (varco2_switch) then

      sib%g%gprogt%co2m = dble(280.+0.27*exp(0.019325329*(year-1700.)))
      !
      ! assign latitudes of flask data
      flask_lat(1)=90.     ! north pole
      flask_lat(2)=71.32   ! BRW (Barrow)
      flask_lat(3)=19.530  ! MLO (Mauna loa)
      flask_lat(4)=-89.980 ! SPO (South Pole Station
      flask_lat(5)=-90.    ! south pole
      !
      ! relative average seasonal cycle in co2 concentration (ppm)
      ! derived from flask observations trends and mean values removed
      ! BRW Season
      flask_co2(2,1)=0.32268710E+01
      flask_co2(2,2)=0.39178782E+01
      flask_co2(2,3)=0.42880478E+01
      flask_co2(2,4)=0.45407581E+01
      flask_co2(2,5)=0.48192682E+01
      flask_co2(2,6)=0.25147221E+01
      flask_co2(2,7)=-0.45384240E+01
      flask_co2(2,8)=-0.10298231E+02
      flask_co2(2,9)=-0.87144508E+01
      flask_co2(2,10)=-0.35581274E+01
      flask_co2(2,11)=0.55232322E+00
      flask_co2(2,12)=0.32477596E+01
      !
      ! MLO Season
      flask_co2(3,1)=-0.21220398E+00
      flask_co2(3,2)=0.57764530E+00
      flask_co2(3,3)=0.16015472E+01
      flask_co2(3,4)=0.27818632E+01
      flask_co2(3,5)=0.33004551E+01
      flask_co2(3,6)=0.26449833E+01
      flask_co2(3,7)=0.18991274E+00
      flask_co2(3,8)=-0.19132940E+01
      flask_co2(3,9)=-0.33598199E+01
      flask_co2(3,10)=-0.31779656E+01
      flask_co2(3,11)=-0.17802587E+01
      flask_co2(3,12)=-0.33074304E+00
      !
      ! SPO Season     
      flask_co2(4,1)=-0.40703124E+00
      flask_co2(4,2)=-0.10044072E+01
      flask_co2(4,3)=-0.10342996E+01
      flask_co2(4,4)=-0.86555952E+00
      flask_co2(4,5)=-0.61372006E+00
      flask_co2(4,6)=-0.31292197E+00
      flask_co2(4,7)=0.49972042E-01
      flask_co2(4,8)=0.52145582E+00
      flask_co2(4,9)=0.80552232E+00
      flask_co2(4,10)=0.91416389E+00
      flask_co2(4,11)=0.89274549E+00
      flask_co2(4,12)=0.79357713E+00
      
      !
      ! Assume south pole has same values as SPO; north pole same as BRW
      do i=1,12
        flask_co2(1,i)=flask_co2(2,i)
        flask_co2(5,i)=flask_co2(4,i)
      enddo

      !print*,sib%g%gprogt%pco2m

      ! calculate co2 concentration
      do n = 1,subcount
        do i=5,1, -1
          if (flask_lat(i)<sib%g(n)%lat) then
            num1=i-1
            num2=i
          endif
        enddo

      ! calculate latitude matrix argument
        latarg(n) = (sib%g(n)%lat - 0.25) + 90.

      !
      ! adjust for latitude dependent seasonal cycle interpolating between flas stations
        seas=(flask_co2(num1,mon)-flask_co2(num2,mon)) &
            / (flask_lat(num1)-flask_lat(num2)) &
            * (sib%g(n)%lat-flask_lat(num2)) &
            + flask_co2(num2,mon)
        sib%g(n)%gprogt%co2m = (sib%g(n)%gprogt%co2m+seas)/1.0D6
        !print *,'ppm gprogt%co2m from co2_mod:', sib%g(n)%gprogt%pco2m
      !
      ! convert from ppmv to pa
        !sib%g(n)%gprogt%pco2m = (sib%g(n)%gprogt%pco2m*p0_sfc)/1.E6
        !sib%g(n)%gprogt%pco2m = dble(sib%g(n)%gprogt%pco2m*pressure)/1.E6    
        !print *,'Pa gprogt%pco2m from co2_mod:', sib%g(n)%gprogt%pco2m
    enddo

  else ! for flat-CO2 case
      do n = 1,subcount
       !print*,'ppm gprogt%pco2m from co2_mod:',bco2m
       !sib%g(n)%gprogt%pco2m = dble(bco2m*pressure)/1.E6
       sib%g(n)%gprogt%co2m = dble(bco2m)/1.0D6
       !print *,'Pa gprogt%pco2m from co2_mod:', sib%g(n)%gprogt%pco2m
      enddo
  endif

end subroutine set_co2
