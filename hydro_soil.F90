
subroutine hydro_soil(ect, egs, egsmax, &
                      rootf, wp_eff,    &
                      soilt, hydrost, sscolt)

!----------------------------------------------------------------------
!
!     soil water - based on CLM_HYDRO_SOIL.F90 code...
!
!----------------------------------------------------------------------
!
!     following taken (almost) verbatim from CLM_HYDRO_SOIL comments...
!
!     main subroutine used to execute the calculation of water 
!     processes over soil.
!
!     1) water within soil (subroutine soil_water)
! 
!     2) runoff
!        The original code was provided by Robert E. Dickinson based on
!        following clues: exponential decrease of Ksat, a water table
!        level determination level including highland and lowland levels
!        and fractional area of wetland (water table above the surface).
!        Runoff is parameterized from the lowlands in terms of precip
!        incident on wet areas and a base flow, where these are estimated
!        using ideas from TOPMODEL.
!
!     The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
!     *  using a new method to determine water table depth and the 
!        fractional wet area (fwsoil).
!     *  computing runoff (surface and subsurface) from this fraction and
!        the remaining fraction (i.e. 1-fwsoil)
!     *  for the 1-fwsoil part, using BATS1e method to compute surface and
!        subsurface runoff.
!
!     The original code on soil moisture and runoff was provided by R.E.
!     Dickinson in July 1996.
!
!     Revision History:
!     15 September 1999: Yongjiu Dai; Initial code
!     12 November  1999: Z.-L. Yang and G.-Y. Niu
!     15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!     25 January   2002: Ian Baker; Integration into SiB
!     09 September 2015: Kathy Haynes; Set runoff as an option due
!                            to equation descrepancies and substantial
!                            differences from the CLM Tech Manual
!----------------------------------------------------------------------

use kinds

use module_pparams, only: &
    denice, denh2o, &
    wtfact 
use module_sibconst, only: &
    nsoil
use module_sib, only: &
    soil_type, hydros_type, &
    sscol_type
use module_time, only: dtsib, dtisib

implicit none

!----------------------------------------------------------------------
!...input variables
real(r8), intent(in) :: ect, egs, egsmax
real(r8), dimension(nsoil), intent(in) :: rootf, wp_eff
type(soil_type), intent(in) :: soilt
type(hydros_type), intent(inout) :: hydrost
type(sscol_type), intent(inout) :: sscolt

!...local variables
integer :: j      ! loop variable
real(r8)    :: zwice    ! total ice mass in soil column !(kg/m2)
real(r8)    :: www_tot(1:nsoil)  ! total pore space occupied by ice+liquid (-)
real(r8)    :: dw_liq(1:nsoil)   ! change in layer liquid-volumetric (m3/m3)
real(r8)    :: dzmm(1:nsoil)     ! soil layer thickness (mm)
real(r8)    :: zmm(1:nsoil)      ! node depth (mm)
real(r8)    :: wmean    ! averaged soil wetness in top layers
real(r8)    :: zmean    ! top soil layers contributing to runoff
real(r8)    :: fz       ! factor (-)
real(r8)    :: zwt      ! water table depth (m)
real(r8)    :: infil    ! infiltration into soil (kg/m2/s)
real(r8)    :: fwsoil   ! saturation fraction
real(r8)    :: watmin   ! minimum soil moisture 
real(r8)    :: xs       ! excess soil moisture (kg/m2) 
real(r8)    :: qqq       ! placeholder for runoff (kg/m2)
real(r8)    :: roff_t    ! placeholder for subsurface runoff (kg/m2)
real(r8)    :: roffo_t   ! placeholder for overland runoff (kg/m2)

!---------------------------------------------------------------------
!...set variables
zwice = 0.0  ! sum of ice mass of soil

do j=1,nsoil
   zwice = zwice + sscolt%www_ice(j)

   sscolt%vol_ice(j) = min(soilt%poros, &
                sscolt%www_ice(j)/(sscolt%dz(j)*denice))
   sscolt%eff_poros(j) = soilt%poros - sscolt%vol_ice(j)
   sscolt%vol_liq(j) = min(sscolt%eff_poros(j), &
                 sscolt%www_liq(j)/(sscolt%dz(j)*denh2o))
   if (sscolt%vol_liq(j) == 0.0 .and. sscolt%www_liq(j) > 0.0 ) then
       hydrost%roff = hydrost%roff + sscolt%www_liq(j) 
       sscolt%www_liq(j) = 0.0
   endif

   www_tot(j) = min(1.0_r8,(sscolt%vol_ice(j)+sscolt%vol_liq(j))/ &
                           soilt%poros)
enddo

   !...determine water table
   wmean = 0.0
   fz    = 1.0
   do j=1,nsoil
       wmean = wmean + www_tot(j)*sscolt%dz(j)
   enddo
   zwt = fz * (abs(sscolt%layer_z(nsoil)) - wmean)

   !...saturation fraction
   fwsoil = wtfact * min(1.0_r8,exp(-zwt))
   !print'(2(a,g16.6))','Water Table Depth, m=',zwt,' Sat Fraction=',fwsoil

   !...these soil calculations are hardwired for a 10-layer soil model
   wmean = 0.0
   zmean = 0.0

   do j=1,3
      zmean = zmean + sscolt%dz(j)
      wmean = wmean + www_tot(j)*sscolt%dz(j)
   enddo

   wmean = wmean/zmean

   !...Modifying infiltration from CLM. 
   !...We put precipitation/snowmelt into surface interception store (capacg).
   !...Any infiltration will come from that reservoir. 
   !...roffo_t represents everything that cannot be infiltrated this timestep.
   roffo_t  = max(0.0_r8,fwsoil*hydrost%capacg*dtisib) +  &
              max(0.0_r8,(1.0 - fwsoil) * &
              min(1.0_r8,wmean**4*hydrost%capacg*dtisib))
   roffo_t = roffo_t * dtsib !kg/m2

   !...infiltration into surface soil layer
   infil = hydrost%capacg - roffo_t  !kg/m2

   !...Capacg is constrained by satcapg, which is set in subroutine radfac.F90.
   hydrost%capacg = hydrost%capacg - infil 
   roffo_t = hydrost%capacg - hydrost%satcapg
   if (roffo_t .gt. 0.) then
      hydrost%roffo = hydrost%roffo + roffo_t  ! kg/m^2
      hydrost%capacg = hydrost%capacg - roffo_t
   endif

!...convert infil from kg/m2 to kg/m2/s
if (infil > 0.0) then
    infil = infil * dtisib
endif

!...set effective root fractions for soil
!...moisture transport
call soil_rootr( soilt%fieldcap, &
                 ect, egs, egsmax,  &
                 rootf, wp_eff,  &
                 sscolt%td(1:nsoil), &
                 sscolt%vol_liq(1:nsoil), &
                 sscolt%www_liq(1:nsoil), &
                 sscolt%rootr)

!...set up r, a, b, and c vectors (following Bonan's (1996) soil)
!...for tridiagonal matrix.
!...(length units will be millimeters)
do j = 1,nsoil
    zmm(j)  = abs(sscolt%node_z(j) *1000.0)
    dzmm(j) = sscolt%dz(j) * 1000.0
enddo

call soil_water( infil, ect, egs, &
                soilt%poros, soilt%satco, &
                soilt%bee, soilt%phsat, &
                dzmm, zmm, &
                sscolt%rootr, &
                sscolt%td(1:nsoil), sscolt%vol_liq(1:nsoil), &
                sscolt%eff_poros(1:nsoil), &
                hydrost%roff, dw_liq )

!...update the liquid water mass (kg/m^2)
do j=1,nsoil
   sscolt%www_liq(j) = sscolt%www_liq(j) + &
                       dw_liq(j)*sscolt%dz(j)*denh2o
enddo

!...limit www_liq to be greater than or equal to watmin
!...get water needed to bring www_liq equal to watmin from lower level
qqq = 0.0
watmin = 0.0
do j=1,nsoil-1
   if (sscolt%www_liq(j) < watmin*sscolt%dz(j)*denh2o ) then
       xs = watmin * sscolt%dz(j) * denh2o - sscolt%www_liq(j)
   else
       xs = 0.0
   endif
  
   if (xs > 0.0) then
      sscolt%www_liq(j)   = sscolt%www_liq(j)   + xs
      sscolt%www_liq(j+1) = sscolt%www_liq(j+1) - xs
   endif
enddo

j = nsoil
if (sscolt%www_liq(j) < watmin) then
    xs = watmin * sscolt%dz(j) * denh2o - sscolt%www_liq(j)
else
    xs = 0.0
endif
sscolt%www_liq(j) = sscolt%www_liq(j) + xs
qqq = -xs


!...determine water in excess of saturation
roff_t = sscolt%www_liq(1) - sscolt%eff_poros(1) * &
          sscolt%dz(1)*denh2o
xs = max(0.0_r8,roff_t)
if (xs > 0.0) then 
    sscolt%www_liq(1) = sscolt%eff_poros(1) * &
       sscolt%dz(1)*denh2o
endif

do j=2,nsoil
    roff_t =  sscolt%www_liq(j) - sscolt%eff_poros(j) * &
          sscolt%dz(j)*denh2o
    xs = xs + max(0.0_r8,roff_t)

    sscolt%www_liq(j) = min(sscolt%eff_poros(j) * &
        sscolt%dz(j)*denh2o,sscolt%www_liq(j))
enddo

!...sub-surface runoff and drainage
qqq = qqq + xs                    
hydrost%roff = hydrost%roff + qqq
hydrost%infil = infil * dtsib !units of kg/m2

end subroutine hydro_soil


!=========================================================================
!=========================================================================
subroutine soil_rootr( fieldcap, &
                       ect, egs, egsmax, &
                       rootf, wp_eff,    &
                       td, vol_liq, www_liq, &
                       rootr)
!=========================================================================

!     subroutine to calculate water-adjusted root fraction
!
!----------------------------------------------------------------------

use kinds
use module_oparams, only:  &
    near_zero
use module_pparams, only:  &
    lvap, tice
use module_sibconst, only: &
    nsoil

implicit none

!...input variables
real(r8), intent(in) :: fieldcap
real(r8), intent(in) :: ect, egs, egsmax
real(r8), dimension(nsoil), intent(in) :: &
     rootf, wp_eff, &
     td, vol_liq, www_liq
real(r8), dimension(nsoil), intent(inout) :: rootr

!...local variables
real(r8) :: btran

!...misc variables
integer(i4) :: i

!----------------------------------------------------------------------

     btran = near_zero

     !.....surface soil layer
     if (egs > dzero ) then
         if (www_liq(1) > near_zero) then
            if ((egs + ect*rootf(1)) > egsmax ) then
                 rootr(1) = (0.5*www_liq(1)*lvap &
                            - egs ) / ect
                 rootr(1) = MAX( MIN(rootr(1), done), dzero)
                 rootr(1) = rootr(1) * rootf(1)
                 if (www_liq(1) <= wp_eff(1)) then
                     rootr(1) = dzero
                 endif
             else
                 rootr(1) = &
                      (1.0 - wp_eff(1)/vol_liq(1)) /  &
                      (1.0 - wp_eff(1)/fieldcap)
                 rootr(1) = MAX( MIN(rootr(1), done), dzero)
                 rootr(1) = rootr(1) * rootf(1)
             endif  !(egs + ect) > egsmax
         else
            rootr(1) = dzero
        endif  !www_liq(1) > near_zero
     else  !egs <= zero
          if (www_liq(1) > near_zero) then
              rootr(1) = &
                   (1.0 - wp_eff(1)/vol_liq(1)) /  &
                   (1.0 - wp_eff(1)/fieldcap)
              rootr(1) = MAX( MIN(rootr(1), done), dzero)
              rootr(1) = rootr(1) * rootf(1)
         else
              rootr(1) = dzero
         endif
     endif

     !.......no transpiration loss from frozen soil.......
     if(td(1) < tice) rootr(1) = dzero

     !.....soil column
     btran = btran + rootr(1)

     do i=2,nsoil
        if (td(i) >= tice .and. vol_liq(i) > 0.0 ) then
            rootr(i) = &
                   (1.0 - wp_eff(i)/vol_liq(i)) /  &
                   (1.0 - wp_eff(i)/fieldcap)
            rootr(i) = MAX( MIN(rootr(i), done), dzero)
            rootr(i) = rootr(i) * rootf(i)
            btran = btran + rootr(i)           
         else  !td >= tice
            rootr(i) = dzero
         endif
     enddo !i=2,nsoil

     !.......normalize to get layer contribution
     if (btran .gt. dzero) then
         do i=1,nsoil
             rootr(i) = rootr(i) / btran
             rootr(i) = max(rootr(i), dzero)
         enddo
     else
         rootr(:) = dzero
     endif

end subroutine soil_rootr


!=========================================================================
!=========================================================================
subroutine soil_water( infil, ect, egs, &
                       poros, satco, bee, phsat, &
                       dzmm, zmm, rootr, &
                       td, vol_liq, eff_poros, &
                       roff, dw_liq)
!=========================================================================
!
!     subroutine based on  CLM_SOILWATER
!
!     CLM web info: http://clm.gsfc.nasa.gov
!
!     Description: (taken directly from CLM_SOILWATER.F90)
!     Soil moisture is predicted from a 10-layer model (as with soil
!     temperature), in which the vertical soil moisture transport is 
!     governed by infiltration, runoff, gradient diffusion, gravity 
!     and root extraction through canopy transpiration.  The net water
!     applied to the surface layer is the snowmelt plus precipitation 
!     plus the throughfall of canopy dew minus surface runoff and 
!     evaporation. 
!
!     The vertical water flow in an unsaturated porous media is 
!     described by Darcy's Law, and the hydraulic conductivity and the
!     soil negative potential vary with soil water content and soil 
!     texture based on the work of Clapp and Hornberger (1978) and Cosby 
!     et al. (1984). The equation is integrated over the layer thickness,
!     in which the time rate of change in water mass must equal the net 
!     flow across the bounding interface, plus the rate of internal source
!     or sink.  The terms of water flow across the layer interfaces are 
!     linearly expanded by using first-order Taylor expansion.  The 
!     equations result in a tridiagonal system of equations.

!     Note: lentgh units here are all millimeter
!
!     Richards Equation
!
!     d wat     d      d wat   d psi
!     ----- = - -- [ k(-----   ----- - 1) ] + S
!      dt       dz       dz    d wat
!
!     where:
!     wat = volumetric water content (m^3/m^3)
!     psi = soil matric potential (mm)
!     dt  = time step (sec)
!     z   = depth (mm)
!     qin = inflow at top of layer (mm H2O/sec)
!     qout= outflow at bottom of layer (mm H2O/sec)
!     s   = source/sink flux (mm H2O/sec)
!     k   = hydraulic conductivity (mm H2O/sec)
!
!     Solution:
!     linearize k and psi about d wat and use tridiagonal system of 
!     equations to solve for d wat, where for layer j
!
!     r_j = a_j[d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
!     Revision History
!     15 September 1999: Yongjiu Dai; Initial code
!     15 December  1999: Paul Houser and Jon Radakovich; F90 Revision
!     25 January   2002: Ian Baker; SiB integration
!----------------------------------------------------------------------

use kinds
use module_pparams, only:  &
    grav,  tice, lvap, wimp, phmin
use module_sibconst, only:  &
    nsnow, nsoil
use module_time, only:  &
    dtsib, dtisib

implicit none

!----------------------------------------------------------------------

!...INPUT VARIABLES
real(r8),intent(in)  :: infil   ! infiltration into soil
                                             !  (kg m^-2 sec^-1) 
real(r8),intent(in)  :: ect, egs
real(r8), intent(in) :: poros, satco, bee, phsat
real(r8), dimension(nsoil), intent(in)  :: dzmm, zmm, rootr
real(r8), dimension(nsoil), intent(in) :: td, vol_liq, eff_poros

!...OUTPUT VARIABLES
real(r8), intent(inout) :: roff
real(r8),intent(out) :: dw_liq(1:nsoil) ! change in layer liquid
                                        !  (m^3/m^3)(volumetric)

!...LOCAL VARIABLES
integer ::   j     ! loop variables
real(r8)  :: lvapi     ! 1/lvap
real(r8)  :: s_node    ! volumetric wetness of node
real(r8)  :: s1        ! wetness at interface of layer
real(r8)  :: s2        ! conductivity*wetness**(2b+2)
real(r8)  :: hk(1:nsoil)     ! hydraulic conductivity
                             ! (mm H2O sec^1)
real(r8)  :: dhkdw(1:nsoil)  ! d(hk)/d(water)
real(r8)  :: smp(1:nsoil)
                                    ! soil matric potential (mm)
real(r8)  :: dsmpdw(1:nsoil)
                                    ! d(smp)/d(wetness)
real(r8)  :: qin       ! flux of water into layer 
                                    !  (mm H2O sec^-1)
real(r8)  :: qout      ! flux of water out of layer 
                                    !  (mm H2O sec^-1)  
real(r8)  :: den       ! used in calculating qin,qout   
real(r8)  :: num       ! used in calculating qin,qout
real(r8)  :: dqidw0    ! d(qin)/d(vol_liq(j-1)) 
real(r8)  :: dqidw1    ! d(qin)/d(vol_liq(j))
real(r8)  :: dqodw1    ! d(qout)/d(vol_liq(j))
real(r8)  :: dqodw2    ! d(qout)/d(vol_liq(j+1))
real(r8)  :: rmx(1:nsoil)
                                    ! "r" forcing term of tridiag matrix
real(r8)  :: amx(1:nsoil)
                                    ! "a" left off diagonal of tridiag mat
real(r8)  :: bmx(1:nsoil)
                                    ! "b" diagonal column for tridiag mat
real(r8)  :: cmx(1:nsoil)
                                    ! "c" right off diagonal of tridiag mat

!----------------------------------------------------------------------

    lvapi = 1.0/lvap
    !...initialize local variables to zero.
    dw_liq(:) = dzero
    hk(:) = dzero
    dhkdw(:) = dzero
    smp(:) = dzero
    dsmpdw(:) = dzero
    rmx(:) = dzero
    amx(:) = dzero
    bmx(:) = dzero
    cmx(:) = dzero

    !...set hydraulic conductivity to zero if effective porosity 5% in 
    !...any two adjoining layers, or if volumetric water (liquid) content
    !...less than 0.001
    do j=1,nsoil
        if( ((eff_poros(j) < wimp)            .or.   &
            (eff_poros(min(nsoil,j+1)) < wimp)) .or. &
            (vol_liq(j) <= 1.E-3)) then
            hk(j)    = 0.0
            dhkdw(j) = 0.0
        else
            s1 = 0.5*(vol_liq(j)+vol_liq(min(nsoil,j+1)))/ &
                 poros

            s2 = satco*1000.0*s1**(2.0*bee+2.0)
            hk(j) = s1*s2

            dhkdw(j) = (2.0*bee+3.0)*s2*0.5/poros
            if (j==nsoil) dhkdw(j) = dhkdw(j)*2.0
        endif

        !...evaluate hydraulic conductivity, soil matric potential,
        !...d(smp)/d(vol_liq) and d(hk)/d(vol_liq)
        if(td(j) > tice) then

            s_node = max(vol_liq(j)/poros,0.01_r8)
            s_node = min(1.0_r8,s_node)

            smp(j) = phsat*1000.0*s_node**(-bee)
            smp(j) = max(phmin,smp(j))
            dsmpdw(j) = -bee*smp(j)/(s_node*poros)

        else

        !...when ice is present, the matric potential is only related to
        !...temperature by (Fuchs et. al., 1978; Soil Sci. Soc. of Amer. J.
        !...42(3); 379-385) 
        !...Unit 1 joule = 1kg m2/sec2 j/kg/(m/s2) ==> m ==>1e3mm

            smp(j) = 1.e3_r8 * 0.3336e6_r8/grav *  &
                (td(j) - tice)/td(j)
            smp(j) = max(phmin,smp(j))
            dsmpdw(j) = 0.0

        endif
    enddo

    !...set up a, b, c and r vectors for tridiagonal solver
    !...node j=1

    j = 1
    qin    = infil
    den    = zmm(j+1) - zmm(j)
    num    = (smp(j+1)-smp(j)) - den
    qout   = -hk(j)*num/den
    dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
    dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
    rmx(j) = qin - qout  - (((ect * rootr(j)) + egs) * dtisib * lvapi)
    amx(j) = 0.0
    bmx(j) = dzmm(j) * (dtisib) + dqodw1
    cmx(j) = dqodw2

    !...nodes 2 through nsoil-1
    do j = 2,nsoil-1 
        den    = zmm(j) - zmm(j-1)
        num    = (smp(j)-smp(j-1)) - den
        qin    = -hk(j-1)*num/den
        dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
        dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den  
        den    = zmm(j+1)-zmm(j)
        num    = smp(j+1)-smp(j) - den
        qout   = -hk(j)*num/den
        dqodw1 = -(-hk(j)*dsmpdw(j)  +  num*dhkdw(j))/den
        dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den
        rmx(j) = qin - qout  - (ect * dtisib * rootr(j) * lvapi)
        amx(j) = -dqidw0
        bmx(j) = dzmm(j)*dtisib - dqidw1 + dqodw1
        cmx(j) = dqodw2
    enddo

    !...node j=nsoil
    j = nsoil
    den    = zmm(j) - zmm(j-1)
    num    = smp(j) - smp(j-1) - den
    qin    = -hk(j-1) * num/den
    dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
    dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den
    qout   = hk(j)
    dqodw1 = dhkdw(j)
    rmx(j) = qin - qout  - (ect * dtisib * rootr(j) * lvapi)
    amx(j) = -dqidw0
    bmx(j) = dzmm(j)*dtisib - dqidw1 + dqodw1
    cmx(j) = 0.0

    ! Add qout out of the bottom layer to runoff
    roff = roff + qout*dtsib

    !...solve
    call  clm_tridia (nsoil, amx, bmx, cmx, rmx, dw_liq)

end subroutine soil_water


!=========================================================================
!=========================================================================
subroutine clm_tridia (n, a, b, c, r, u )

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!
! REVISION HISTORY:
!  15 September 1999: Yongjiu Dai; Initial code
!  15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision 
!  15 January 2002: Ian Baker revision to work in SiB
!=========================================================================
! $Id: tridiag_solver.F90,v 1.1 2005/04/06 22:24:04 chorak Exp $
!=========================================================================

    use kinds
    implicit none

    !=== Arguments ===========================================================

    integer , intent(in)  :: n
    real(r8), intent(in)  :: a(1:n),b(1:n),c(1:n),r(1:n)
    real(r8), intent(out) :: u(1:n)

    !=== Local Variables =====================================================

    integer j
    real(r8) gam(1:n),bet

    !=== End Variable List ===================================================


    bet  = b(1)
    u(1) = r(1) / bet
    do j = 2, n
        gam(j) = c(j-1) / bet
        bet = b(j) - a(j) * gam(j)
        u(j) = (r(j) - a(j)*u(j-1)) / bet
    enddo

    do j = n-1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
    enddo

end subroutine clm_tridia

