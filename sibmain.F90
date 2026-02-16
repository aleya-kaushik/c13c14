!----------------------------------------------------------------------
subroutine sib_main ( gref, lonsib, latsib, &
       pref, pnum, physcont,  &
       gprogt, gdiagt, sibgl )

    use kinds
    use module_local
    use module_oparams, only: &
       canlai_min, cos_i3
    use module_param, only: &
       phys_param
    use module_pftinfo
    use module_sib, only: &
       gprog_type, gdiag_type, lu_type

    implicit none

!-----------------------------------------------------------------------

!     REFERENCES: 
!     Baker, I. T., A. B. Harper, H. R. da Rocha, A. S. Denning, 
!          A. C. Araujo, L. S. Borma, and S. C. Wofsy (2013), Surface
!          ecophysiological behavior across vegetation and moisture
!          gradients in tropical South America.  Ag. Forest Met.,
!          doi:10.1016/j.agformet.2012.11.015.
!
!     Lokupitiya, E. Y., A. S. Denning, K. Paustian, I. T. Baker, 
!          K. Schaefer, S. Verma, T. Meyers, C. Bernacchi, A Suyker,
!          and M. Fischer (2009), Incorporation of crop phenology in the
!          Simple Bioshphere Model (SiBcrop) to improve land-atmopshere
!          carbon exchanges from croplands.  Biogeosciences, 6, 969-986.
!
!     Schaefer, K., P. Tans, S. Denning, I. Baker, J. Berry, L. Prihodko,
!          N. Suits, and A. Philpott (2008), The combined Simple 
!          Biosphere/Carnegie-Ames-Stanford Approach (SiBCASA) model.
!          J. Geophys. Res., 113, G03034, doi:10.1029/2007JG000603.
!
!     Sellers, P. J., D. A. Randall, C. J. Collatz, J. A. Berry,
!          C. B. Field, D. A. Dazlich, C. Zhang, G. Collelo (1996), A revised 
!          land-surface parameterization (SiB2) for atmospheric GCMs. Part 1:
!          Model formulation.  J. Clim, 9, 676-705.
!
!     Sellers, P.J., S. O. Los, C. J. Tucker, C. O. Justics, D. A. Dazlich,
!          G. J. Collatz, and D. A. Randall (1996), A revised land
!          surface parameterization (SiB2) for atmopsheric GCMs Part II:
!          The generation of global fields of terrestrial biospherical
!          parameters from satellite data.  J. Clim, 9, 706-737.
!
!-----------------------------------------------------------------------

!---INPUT VARIABLES
integer(i4), intent(in) :: gref, pref, pnum
real(r4), intent(in) :: lonsib, latsib
type (phys_param), intent(in) :: physcont
type (gdiag_type), intent(inout) :: gdiagt
type (gprog_type), intent(inout) :: gprogt
type (lu_type), intent(inout) :: sibgl

!---LOCAL VARIABLES
integer(byte) :: nsl  !local number of snow layers

!-------------------------------------------------------------------------

!--SET LOCAL VARIABLES
call local_set(gprogt%ps, sibgl%cast%tc, &
     sibgl%hydrost%capacc_liq, sibgl%hydrost%capacc_snow, &
     sibgl%hydrost%capacg, sibgl%hydrost%snow_gmass, &
     sibgl%sscolt%nsl, sibgl%sscolt%www_liq, &
     sibgl%sscolt%www_ice, sibgl%sscolt%td)

!--SET HYDROLOGY DIAGNOSTICS
!...soil properties
call hydro_sets(sibgl%soilt, sibgl%hydrost, sibgl%sscolt)

!...vegetation root zone stresses
call hydro_setv(sibgl%soilt%fc_eff, sibgl%soilt%wp_eff, &
     sibgl%vegt%rootf, sibgl%sscolt, sibgl%hydrovt)    

!--UPDATE RADIATION
!...calculate albedo/reflectance/transmissivity
call radfac(physcont%chil, physcont%ref, &
     physcont%tran, physcont%z1, physcont%z2, &
     gdiagt%cosz, sibgl%soilt, sibgl%sscolt, sibgl%vegt,   &
     sibgl%hydrost, sibgl%radt)

!...calculate the absorbed radiation (radc3)
call radabs(gprogt%dlwbot, gdiagt, sibgl%radt)

!...calculate the net radiation (radt)
call radnet(sibgl%sscolt%nsl, sibgl%radt%tsfc, &
     sibgl%cast%tc, sibgl%radt)

!--UPDATE CAS CAPACITIES
call cas_update(gdiagt%bps(1), gdiagt%psy,  &
     gdiagt%ros, sibgl%co2t%casd, &
     sibgl%hydrost%capacc_liq, sibgl%hydrost%capacc_snow, &
     sibgl%vegt%lai, sibgl%cast)

!--UPDATE CANOPY CONDUCTANCES/RESISTANCES
!.....calculate canopy and ground resistances
nsl = sibgl%sscolt%nsl
call flux_vrbrd(pnum, physcont%z2, nsl, &
     sibgl%co2t%rst, sibgl%sscolt%td(nsl+1), &
     sibgl%cast, sibgl%vegt, gprogt, gdiagt,  &
     sibgl%fluxt, sibgl%hydrost)

!--UPDATE CARBON CYCLE
!..calculate photosynthesis stress factors
call phostress(pnum, physcont, &
     gprogt%ps, etc, sibgl%cast%tc, sibgl%cast%eacas, &
     sibgl%fluxt%rb, sibgl%hydrost%ecmass,  &
     sibgl%sscolt%td(1), sibgl%sscolt%td(2), &
     sibgl%hydrovt%pawfzw, sibgl%hydrovt%tawfrw, &
     sibgl%cast%tcmin, sibgl%co2t, sibgl%cast%vpd, sibgl%fluxt%ect)

if (sibgl%vegt%lai .gt. canlai_min) then
   !..calculate canopy conductance and photosynthesize
   call phosib(physcont, gprogt, gdiagt, &
        sibgl%cast%tc, sibgl%cast%tcas, sibgl%fluxt%ra, &
        sibgl%fluxt%rb, sibgl%poollt%resp_leaf, &
        sibgl%poollt%resp_auto, sibgl%pooldt%resp_het, &
        sibgl%vegt, sibgl%co2t)

   !..calculate 13CO2 fractionation
   call cfrax_calc(physcont, gprogt, sibgl%co2t, &
         sibgl%fract, sibgl%fluxt, &! sibgl%co2t%co2cas, &
         !sibgl%co2t%co2s, sibgl%co2t%co2i, sibgl%co2t%co2c, &
         sibgl%co2t%co2m, sibgl%co2t%co2gamma, sibgl%co2t%assim)

   !..calculate COS
   call cos_calc(gref, pref, lonsib, latsib, &
         cos_i3, gprogt%pcosm, sibgl%cast%tcas, &
         sibgl%co2t%aparkk, sibgl%co2t%assim, &
         sibgl%co2t%rstfac(2), sibgl%soilt%poros, &
         sibgl%soilt%zm, sibgl%soilt%woptzm, &
         sibgl%soilt%wsat, sibgl%vegt%rootf(1:cos_i3), &
         sibgl%sscolt%dz(1:cos_i3), sibgl%sscolt%www_ice(1:cos_i3), &
         sibgl%sscolt%www_liq(1:cos_i3), sibgl%pooldt%resp_soilnr_lay, &
         sibgl%pooldt%mhrt_soil_hot, sibgl%cast%tc, &
         sibgl%co2t, sibgl%cost, pref, sibgl%sscolt)

   !..calculate SIF
   call sif_calc(sibgl%co2t%nspar, sibgl%vegt%fpar, &
         sibgl%co2t%assim, sibgl%co2t%assimpot, &
         sibgl%cast%tc, sibgl%vegt%vmax, sibgl%sift)

else
   !call phonveg(gdiagt%radvbc, gdiagt%radvdc, &
   !             gprogt%co2m, sibgl%co2t, gprogt%ps)
   call phonveg(gdiagt%radvbc, gdiagt%radvdc, &
                gprogt, sibgl%co2t, gprogt%ps)
endif

!----UPDATE / INCREMENT PROGNOSTIC VARIABLES-----
!.....calculate partial derivatives of the various heat fluxes
!.....with respect to ground/canopy/snow temp
call delwf()

call delef(gdiagt%psy, gdiagt%ros, gdiagt%em, &
     sibgl%cast%eacas, sibgl%hydrost%rhsoil, &
     sibgl%fluxt%rd, sibgl%fluxt%ra, sibgl%fluxt%ec, &
     sibgl%fluxt%eg, sibgl%fluxt%es, sibgl%fluxt%fws)

call delhf(gprogt%tm, gdiagt%bps, gdiagt%ros, &
     sibgl%cast%tcas, sibgl%sscolt%nsl, sibgl%radt%tsfc, &
     sibgl%fluxt%fss, sibgl%fluxt%hc, sibgl%fluxt%hg, &
     sibgl%fluxt%hs, sibgl%cast%tc, sibgl%fluxt%ra, &
     sibgl%fluxt%rb, sibgl%fluxt%rd)

!...solve matrix of prognostic variables
call sibslv(sibgl%cast, sibgl%fluxt, sibgl%radt, sibgl%sscolt)

!...update prognostic variables to get total latent and sensible fluxes
call addinc(gref, pref, lonsib, latsib, gdiagt, &
     gprogt, sibgl%cast, sibgl%fluxt, sibgl%sscolt)

!--CALCULATE FLUXES
!...calculates fluxes and updates 
!...canopy and ground surface water stores
call flux_update(gdiagt%psy, gdiagt%ros, sibgl%radt%radtc, &
     sibgl%radt%radtg, sibgl%radt%radts, sibgl%radt%radc3c, &
     sibgl%radt%radc3g, sibgl%soilt%poros, &
     sibgl%hydrovt%paw_lay, sibgl%cast, sibgl%fluxt, &
     sibgl%hydrost, sibgl%sscolt, gprogt%ps, gprogt%spdm, &
     gprogt%dlwbot)

!--UPDATE SOIL/SNOW HYDROLOGY STORES
!...update water storage/interception on canopy,
!...and calculate throughfall
if (sibgl%vegt%lai .gt. canlai_min) then
   call hydro_canopy(gref, lonsib, latsib, pref, &
         physcont%chil, gprogt%tm,  &
         sibgl%cast%tcas, sibgl%vegt%lait, sibgl%vegt%vcover, &
         gprogt%cuprt, gprogt%lsprt, &
         sibgl%cast%tc, sibgl%hydrost, sibgl%sscolt)
else
    call hydro_nveg(gprogt%tm, sibgl%cast%tcas, &
         gprogt%cuprt, gprogt%lsprt, sibgl%hydrost, sibgl%sscolt)
endif

!...update precipitation onto snowcover
call hydro_snow(sibgl%soilt%poros, sibgl%hydrost, sibgl%sscolt)

!...update infiltration and soil water
call hydro_soil(sibgl%fluxt%ect, sibgl%fluxt%egs, &
     sibgl%fluxt%egsmax, sibgl%vegt%rootf, sibgl%soilt%wp_eff, &
     sibgl%soilt, sibgl%hydrost, sibgl%sscolt)

!...update snow layers, if present
nsl = sibgl%sscolt%nsl
if (nsl < 0) then
   !.....compact snow levels
   call snow_compact(nsl, sibgl%sscolt%td(nsl+1:0), &
        sibgl%sscolt%www_liq(nsl+1:0), &
        sibgl%sscolt%www_ice(nsl+1:0), &
        sibgl%sscolt%dz(nsl+1:0))

   !.....combine snow levels when necessary
   call snow_combine(gref, lonsib, latsib, pref, &
        sibgl%vegt%lai, sibgl%soilt%poros, &
        sibgl%hydrost, sibgl%sscolt)

   !.....subdivide snow levels when necessary
   call snow_subdivide(sibgl%sscolt)

endif  !snow layers present

!--ENERGY AND WATER BALANCE CHECK
call balan_eh2o(gref, lonsib, latsib, pref, &
     gprogt%lspr, gprogt%cupr, sibgl%vegt%lai, &
     sibgl%fluxt, sibgl%sscolt, sibgl%hydrost, &
     sibgl%fluxt%ebalnum, sibgl%fluxt%wbalnum)

end subroutine sib_main
