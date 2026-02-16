! Sets up the SiB4 diagnostic variables for a simulation
subroutine sibtype_setup()

use kinds

use module_oparams, only: &
     cos_casd_min
use module_pparams, only: &
     pi180, tref
use module_param
use module_pftinfo, only: &
     pft_num, pft_type, &
     type_bare, type_crop, type_grass
use module_phosib, only: co2_casd_min
use module_sib, only: sib
use module_sibconst, only: &
     nsoil, subcount, subset, &
     sublonsib, sublatsib, &
     tan_dec
use module_sibvs, only: sibvs
use module_time, only: &
     daylenmax, starttime

implicit none

!local variables
integer(i4) :: gref
real(r4) :: glon, glat
integer(i4) :: g,l

integer(byte) :: ptyperef
integer(i4) :: dnum, pref, pnum
logical :: isbare, iscrop, isgrass


!------------------------------------------
print*,'Setting SiB4 Variables'

! Set the diagnostic variables
do g=1, subcount
   gref = subset(g)
   glon = sublonsib(g)
   glat = sublatsib(g)

   !...Daylength
   call day_lengthpt(pi180, tan_dec, glat, &
        sib%g(g)%gdiagt%daylen, sib%g(g)%gdiagt%daylendt)

   IF (abs(sib%g(g)%gdiagt%daylendt) .lt. 1.E-6) THEN
       call day_lengthpt(pi180, tan_dec, glat, &
            sib%g(g)%gdiagt%daylen, sib%g(g)%gdiagt%daylendt)
   ENDIF

   do l=1, sib%g(g)%g_nlu

      pref = sib%g(g)%l(l)%ipft
      pnum = pft_num(pref)
      ptyperef = pft_type(pnum)
      isbare = (ptyperef .eq. type_bare)
      iscrop = (ptyperef .eq. type_crop)
      isgrass = (ptyperef .eq. type_grass)

      !...Time Invariant
      !Set soil properties
      call setup_soilt( &
          sibvs(gref)%clayfrac(l),  &
          sibvs(gref)%sandfrac(l),  &
          sibvs(gref)%soref_vis(l), &
          sibvs(gref)%soref_nir(l), &
          physcon(pnum)%fc_min, &
          physcon(pnum)%wp_min, &
          sib%g(g)%l(l)%soilt)

      !Set soil column variables
      call setup_sscolt( &
          sib%g(g)%l(l)%soilt%poros, &
          sib%g(g)%l(l)%sscolt)

      !Set vegetation rooting fraction
      call setup_rootf( &
          physcon(pnum)%kroot, &
          physcon(pnum)%rootd, &
          sib%g(g)%l(l)%sscolt%dz(1:nsoil), &
          sib%g(g)%l(l)%sscolt%layer_z(1:nsoil), &
          sib%g(g)%l(l)%sscolt%node_z(1:nsoil),  &
          sib%g(g)%l(l)%vegt%rootf)    


      !...Time-Step Varying
      !Set co2 variables
      sib%g(g)%l(l)%co2t%casd = MAX(dble(co2_casd_min), &
                                    dble(physcon(pnum)%z2))

      !Set cos variables
      sib%g(g)%l(l)%cost%cos_casd = MAX(dble(cos_casd_min), &
                                        dble(physcon(pnum)%z2))

      !Set initial cos soil profile
      sib%g(g)%l(l)%cost%cos_s = 2d-8 ! ~500 ppt

      !Set pool variables
      call setup_poolt( poolcon(pnum), &
           sib%g(g)%l(l)%vegt%rootf, &
           sib%g(g)%l(l)%equibdt, sib%g(g)%l(l)%pooldt, &
           sib%g(g)%l(l)%equiblt, sib%g(g)%l(l)%poollt, &
           sib%g(g)%l(l)%fract)

      !Set hydrological snow/soil variables
      call setup_hydrost( physcon(pnum)%z1, &
           physcon(pnum)%z2, &
           sib%g(g)%l(l)%sscolt, &
           sib%g(g)%l(l)%hydrost)
      call hydro_sets(sib%g(g)%l(l)%soilt, &
           sib%g(g)%l(l)%hydrost, sib%g(g)%l(l)%sscolt)

      !Set hydrological vegetation variables
      call hydro_setv(  &
           sib%g(g)%l(l)%soilt%fc_eff, &
           sib%g(g)%l(l)%soilt%wp_eff, &
           sib%g(g)%l(l)%vegt%rootf, &
           sib%g(g)%l(l)%sscolt, sib%g(g)%l(l)%hydrovt)

      !Set vegetation variables
      dnum = MAX(starttime-1, 0)
      call veg_update( dnum, &
          gref, glon, glat, &
          pnum, pref, iscrop, isgrass, &
          physcon(pnum), &
          sib%g(g)%l(l)%pooldt%poollu, &
          sib%g(g)%l(l)%hydrost%snow_cvfc, &
          sib%g(g)%l(l)%sscolt%node_z(1:nsoil), &
          sib%g(g)%l(l)%poollt, sib%g(g)%l(l)%vegt)
      !Set phenology variables
      if (.not. isbare) then
          call setup_phent( &
               g, l, pnum, pref, &
               phencon(pnum), poolcon(pnum),  &
               sib%g(g)%gdiagt%daylen, sib%g(g)%gdiagt%daylendt, &
               daylenmax(g), &
               sib%g(g)%l(l)%vegt%lai, sib%g(g)%gprogt%tm, &
               sib%g(g)%l(l)%vegt%vmax, sib%g(g)%l(l)%pooldt, &
               sib%g(g)%l(l)%poollt, sib%g(g)%l(l)%phent, &
               physcon(pnum))
       endif
   enddo !l=1, g_nlu

   !Set grid cell diagnostic variables
   sib%g(g)%gdiagt%tmdf = 32.0 &
       + (sib%g(g)%gprogt%tmd-tref * 1.8)

enddo !g=1, subcount


end subroutine sibtype_setup


!------------------------------------------
! SET-UP TIME-INVARIANT VARIABLES
!------------------------------------------

!----------------------------------------
!Routine to set the soil variables
!
!...These are set using %sand and %clay,
!...using relationships from:
!   Clapp and Hornberger, 1978; 
!   Cosby et al., 1984; and
!   Lawrence and Slater, 2008
subroutine setup_soilt( &
     clayfrac, sandfrac, &
     soref_vis, soref_nir, &
     fc_min, wp_min, &
     soilt)

use kinds
use module_pparams, only: &
    wpotfc, wpotwp
use module_sib, only: &
    soil_type
use module_sibconst, only: nsoil, &
    print_soil, print_stop

implicit none

!...input variables
real(r8), intent(in) :: clayfrac, sandfrac
real(r8), intent(in) :: soref_vis, soref_nir
real(r4), dimension(nsoil/2), intent(in) :: &
    fc_min, wp_min 
type(soil_type), intent(inout) :: soilt

!...local variables
integer(i4) :: count, k
real(r8) :: cfrac, sfrac
real(r8) :: tkmsoil  ! soil mineral conductivity (W/m K)
real(r8) :: bd       ! bulk density of dry soil material (kg/m3)

     !...basic soil properties
     soilt%clayfrac = clayfrac
     soilt%sandfrac = sandfrac
     soilt%soref_vis = soref_vis
     soilt%soref_nir = soref_nir

     cfrac = clayfrac*100.
     sfrac = sandfrac*100.

     soilt%poros = 0.489-0.00126*sfrac
     soilt%satco = 0.0070556*10**(-0.884+0.0153*sfrac)/1000.

     !...misc soil variables with depth
     soilt%csolid = (2.128*sandfrac + 2.385*clayfrac) / &
                      (sandfrac+clayfrac)*1.0E6

     bd = (1.0 - soilt%poros) * 2.7E3
     soilt%tkdry = (0.135*bd + 64.7) / (2.7E3 - 0.947*bd)

     tkmsoil = (8.8*sandfrac + 2.92*clayfrac ) / &
               (sandfrac + clayfrac)
     soilt%tkmg = tkmsoil**(1.0-soilt%poros) 
     soilt%tksat = soilt%tkmg*0.57**soilt%poros

     soilt%bee = 2.91+0.159*cfrac
     soilt%phsat = -10.*10.**(1.88-0.0131*sfrac)/1000.
     soilt%fieldcap = soilt%poros * &
            ((wpotfc/9.8)/soilt%phsat) ** (-1.0/soilt%bee)
     soilt%vwcmin = soilt%poros * &
            ((wpotwp/9.8)/soilt%phsat) ** (-1.0/soilt%bee)

     soilt%wopt = (-0.08*clayfrac**2 + &
                    0.22*clayfrac+0.59)*100.
     soilt%wsat = 0.25*clayfrac+0.5
     soilt%zm = -2*clayfrac**3 - 0.4491*clayfrac**2 + &
                    0.2101*clayfrac+0.3478
     soilt%woptzm = (soilt%wopt/100.)**soilt%zm

     !...effective field capacity and wilting point
     count=1
     do k=1,nsoil
        soilt%fc_eff(k) = MIN(soilt%fieldcap, &
                              fc_min(count))
        soilt%wp_eff(k) = MIN(soilt%vwcmin, &
                              wp_min(count))
        if (mod(k,2) .eq. 0) count=count+1
     enddo

     !...print if requested
     if (print_soil) then 
         print*,''
         print('(a,2F6.3)'),' Clay/Sand Fracs: ', clayfrac, sandfrac
         print('(a,2F6.3)'),' Soil Reflectivity (vis/nir): ', &
              soref_vis, soref_nir
         print('(a,F6.3)'),' Soil Tension at Sat: ', soilt%phsat
         print('(a,F6.3)'),' Porosity: ',soilt%poros
         print('(a,2F6.3)'),' Field Capacity/Wilting Point: ', &
              soilt%fieldcap, soilt%vwcmin
         do k=1,nsoil
             print('(a,i2,2F6.3)'),' Effective FC/WP: ', &
                   k,soilt%fc_eff(k), soilt%wp_eff(k)
         enddo
         print*,''

         if (print_stop) stop
     endif

end subroutine setup_soilt


!----------------------------------------
!Routine to set the soil column variables
subroutine setup_sscolt(poros, sscolt)

use kinds
use module_oparams, only: wsat_default
use module_pparams, only: denh2o
use module_sib, only: sscol_type
use module_sibconst, only: &
    nsnow, nsoil, &
    print_sscol, print_stop

implicit none

!...input variables
real(r8), intent(in) :: poros
type(sscol_type), intent(inout) :: sscolt

!...local variables
integer(i4) :: j,k
real(r8) :: scalez
real(r8) :: www_tot


      !...chose soil layer vertical scaling factor
      !scalez = 0.025 !...used for SiBCrop
      scalez = 0.073

      !...soil layers
      do k=1, nsoil
          sscolt%node_z(k) = (scalez*(exp(0.5*(k-0.5))-1.0))
      enddo

      sscolt%dz(1) = 0.5*(sscolt%node_z(1) + sscolt%node_z(2))
      do k=2,nsoil-1
         sscolt%dz(k) = 0.5*(sscolt%node_z(k+1) - sscolt%node_z(k-1))
      enddo
      sscolt%dz(nsoil) = sscolt%node_z(nsoil) - sscolt%node_z(nsoil-1)

      do k=1,nsoil-1
          sscolt%layer_z(k) = 0.5*(sscolt%node_z(k) + sscolt%node_z(k+1))
      enddo
      sscolt%layer_z(nsoil) = sscolt%node_z(nsoil) + 0.5*sscolt%dz(nsoil)

      !...snow layers
      do j=0,sscolt%nsl+1,-1
         sscolt%node_z(j) = sscolt%layer_z(j) - 0.5*sscolt%dz(j)
         sscolt%layer_z(j-1) = sscolt%node_z(j) - 0.5*sscolt%dz(j)
      enddo

      !...snow/soil column moisture
      www_tot = sum(sscolt%www_liq(1:nsoil)) + sum(sscolt%www_ice(1:nsoil))
      if (www_tot .le. 1.E-12) then
          !print*,'   Saturated Soil Moisture'
          do k=sscolt%nsl+1,nsoil
              sscolt%www_liq(k) = &
                  sscolt%dz(k) * denh2o * poros * wsat_default
          enddo
      endif

      !...print soil/snow layer info
      IF (print_sscol) THEN
          print*, '     Number of Snow Layers: ', abs(sscolt%nsl)
          print'(2a)','        lev   layer (m)  node (m)  dz (m)      ', &
                      'td (K)      liq (kg/m2)    ice (kg/m2)'
          do k=-nsnow+1,nsoil
              print'(a,i6,3f10.5,3f14.5)','     ',k, &
                  sscolt%layer_z(k), sscolt%node_z(k), &
                  sscolt%dz(k), sscolt%td(k), &
                  sscolt%www_liq(k), sscolt%www_ice(k)
          enddo
          if (print_stop) stop
      ENDIF

end subroutine setup_sscolt


!----------------------------------------
!Routine to set the vegetation rooting fraction
subroutine setup_rootf( &
      kroot, rootd, &
      dz, layer_z, node_z, rootf)

use kinds
use module_sibconst, only: nsoil

implicit none

!...input variables
real(r4), intent(in) :: kroot, rootd
real(r8), dimension(nsoil), intent(in) :: &
      dz, layer_z, node_z
real(r8), dimension(nsoil), intent(inout) :: &
      rootf

!...local variables
integer(i4) :: k
real(r8) :: rtot
real(r8) :: zbot, ztop
real(r8) :: totalroot

      !...calculate root distribution
      ztop = dzero
      if (kroot > 0.) then
          totalroot = (1.0 - exp(-kroot*layer_z(nsoil)))/kroot

          do k=1,nsoil
              zbot = ztop + dz(k)
              rootf(k) = (exp(-kroot*ztop) - exp(-kroot*zbot)) &
                          / (kroot*totalroot)
              ztop = zbot
          enddo
      endif

      !...modify root distribution
      if (rootd < layer_z(nsoil)) then
          rtot = dzero

          do k=1,nsoil
              if (node_z(k) > rootd) then
                  rootf(k) = dzero
              else
                  rtot = rtot + rootf(k)
              endif
          enddo

          if (rtot > 0.) then
              do k=1,nsoil
                  rootf(k) = rootf(k)/rtot
              enddo
          endif
      endif

end subroutine setup_rootf


!------------------------------------------
! SET-UP TIME-STEP VARYING VARIABLES
!------------------------------------------

!-----------------------------------------------
!Routine to set hydrological snow/soil variables
subroutine setup_hydrost( z1, z2, &
      sscolt, hydrost)
 
use kinds
use module_pparams, only:  &
     cwlim, gwlim, gwctog, &
     denice
use module_sib, only: &
     soil_type, sscol_type, &
     hydros_type
use module_sibconst, only: nsnow

implicit none

!...input variables
real(r4), intent(in) :: z1, z2
type(sscol_type), intent(in) :: sscolt
type(hydros_type), intent(inout) :: hydrost

!...local variables
integer(i4) :: k

      hydrost%satcapc = cwlim
      hydrost%satcapg = gwlim
      hydrost%wetfracc = &
          MAX(dzero, MIN(done, &
              (hydrost%capacc_liq + &
               hydrost%capacc_snow) / &
               hydrost%satcapc))
      hydrost%wetfracg = &
          MAX(dzero, MIN(done, gwctog * &
              (hydrost%capacg / &
               hydrost%satcapg)))

      do k=-nsnow+1,0
          hydrost%snow_gmass = hydrost%snow_gmass + &
              sscolt%www_liq(k) + sscolt%www_ice(k)
          hydrost%snow_gdepth = hydrost%snow_gdepth + &
              sscolt%dz(k)
      enddo

      hydrost%www_tot = sum(sscolt%www_liq) + &
          sum(sscolt%www_ice)

      hydrost%snow_cvfc = &
          ((hydrost%capacc_snow/denice)*5. - z1) &
            / (z2-z1)
      hydrost%snow_cvfc = MIN(done, &
          MAX(dzero, hydrost%snow_cvfc))

end subroutine setup_hydrost


!-------------------------------------
!Routine to set the pool variables
subroutine setup_poolt( poolcont, &
      rootf, equibdt, pooldt, &
      equiblt, poollt, fract)
 
use kinds
use module_poolinfo, only: &
    pool_indx_lay
use module_sib, only: &
    equibd_type, equibl_type, &
    poold_type, pooll_type, &
    fract_type
use module_sibconst, only: &
    npoollu, npoolpft, nsoil, &
    spinup, spinup_default, &
    spinup_continue
use module_param, only: &
    pool_param

implicit none

!...input variables
!real(r8), dimension(npoolpft), intent(in) :: poolpft_min
type(pool_param), intent(in) :: poolcont
real(r8), dimension(nsoil), intent(in) :: rootf
type(equibd_type), intent(inout) :: equibdt
type(poold_type), intent(inout) :: pooldt
type(equibl_type), intent(inout) :: equiblt
type(pooll_type), intent(inout) :: poollt
type(fract_type), intent(in) :: fract

!...local variables
integer(i4) :: p,dlay,tcref

!----------------------------------------
!...reset pools if spinning up
if ((spinup) .and. (spinup_default) .and. &
    (.not. spinup_continue)) then
    pooldt%poollu(:) = dzero
    poollt%poolpft(:) = poolcont%poolpft_min
    poollt%poolpftp(:) = poolcont%poolpft_min
endif

!...set vertical distribution information
do p=1,npoollu/2 !1,6
   dlay = pool_indx_lay(p+npoolpft/2) !6,12
   if (dlay .eq. 1) then
       pooldt%poollu_flay(p,1) = done
   else
       pooldt%poollu_flay(p,:) = rootf(:)
   endif
enddo

do p=1,npoolpft/2 !1,5
   dlay = pool_indx_lay(p) !1,5
   if (dlay .eq. 1) then
      poollt%poolpft_flay(p,1) = done
   else
      poollt%poolpft_flay(p,:) = rootf(:)
   endif
enddo

!...same as above but for C13 pools
do p=npoollu/2+1,npoollu !7,12
   dlay = pool_indx_lay(p+npoolpft) !17,22
   if (dlay .eq. 1) then
       pooldt%poollu_flay(p,1) = done
   else
       pooldt%poollu_flay(p,:) = rootf(:)
   endif
enddo

do p=npoolpft/2+1,npoolpft !6,10
   dlay = pool_indx_lay(p+npoolpft/2+1) !12,16
   if (dlay .eq. 1) then
      poollt%poolpft_flay(p,1) = done
   else
      poollt%poolpft_flay(p,:) = rootf(:)
   endif
enddo


!...set pool vertical distribution and
!...equilibrium information
do p=1,npoollu/2
   pooldt%poollup(p) = pooldt%poollu(p)
   if (sum(pooldt%poollu_lay(p,:)) .ne. pooldt%poollu(p)) then
       pooldt%poollu_lay(p,:) = pooldt%poollu(p) &
           * pooldt%poollu_flay(p,:)
   endif

!   if (.not. spinup_continue) then
      equibdt%poollu_init(p) = pooldt%poollu(p)
      equibdt%poollu_min(p) = pooldt%poollu(p)
      equibdt%poollu_max(p) = pooldt%poollu(p)
!   else
!      equibdt%poollu_init(p) = equibdt%poollu_equib(p)
!      equibdt%poollu_min(p) = equibdt%poollu_min(p)
!      equibdt%poollu_max(p) = equibdt%poollu_max(p)
!   endif
enddo

do p=npoollu/2+1,npoollu
   tcref=p-6
   pooldt%poollup(p) = pooldt%poollu(p)
   if (sum(pooldt%poollu_lay(p,:)) .ne. pooldt%poollu(p)) then
       pooldt%poollu_lay(p,:) = pooldt%poollu(p) &
           * pooldt%poollu_flay(p,:)
   endif

!   if (.not. spinup_continue) then
      equibdt%poollu_init(p) = pooldt%poollu(p)
      equibdt%poollu_min(p) = pooldt%poollu(p)
      equibdt%poollu_max(p) = pooldt%poollu(p)
!   else
!      equibdt%poollu_init(p) = equibdt%poollu_equib(p)
!      equibdt%poollu_min(p) = equibdt%poollu_min(p)
!      equibdt%poollu_max(p) = equibdt%poollu_max(p)
!   endif
enddo


do p=1,npoolpft/2
   poollt%poolpftp(p) = poollt%poolpft(p)
   if (sum(poollt%poolpft_lay(p,:)) .ne. poollt%poolpft(p)) then
       poollt%poolpft_lay(p,:) = poollt%poolpft(p) &
           * poollt%poolpft_flay(p,:)
   endif

!   if (.not. spinup_continue) then
      equiblt%poolpft_init(p) = poollt%poolpft(p)
      equiblt%poolpft_min(p) = poollt%poolpft(p)
      equiblt%poolpft_max(p) = poollt%poolpft(p)
!   else
!      equiblt%poolpft_init(p) = equiblt%poolpft_equib(p)
!      equiblt%poolpft_min(p) = equiblt%poolpft_min(p)
!      equiblt%poolpft_max(p) = equiblt%poolpft_max(p)
!   endif
enddo

do p=npoolpft/2+1,npoolpft
   tcref=p-5
   poollt%poolpftp(p) = poollt%poolpft(p)
   if (sum(poollt%poolpft_lay(p,:)) .ne. poollt%poolpft(p)) then
       poollt%poolpft_lay(p,:) = poollt%poolpft(p) &
           * poollt%poolpft_flay(p,:)
   endif

!   if (.not. spinup_continue) then
      equiblt%poolpft_init(p) = poollt%poolpft(p)
      equiblt%poolpft_min(p) = poollt%poolpft(p)
      equiblt%poolpft_max(p) = poollt%poolpft(p)
!   else
!      equiblt%poolpft_init(p) = equiblt%poolpft_equib(p)
!      equiblt%poolpft_min(p) = equiblt%poolpft_min(p)
!      equiblt%poolpft_max(p) = equiblt%poolpft_max(p)
!   endif
enddo

end subroutine setup_poolt


!-------------------------------------
!Routine to set the phenology variables
subroutine setup_phent( &
      gnum, lnum, pnum, pref, &
      phencont, poolcont, &
      daylen, daylendt, daylenmax, &
      lai, tm, vmax, &
      pooldt, poollt, phent, physcont)

use kinds
use module_param, only: &
   phen_param, pool_param, phys_param
use module_pftinfo, only: &
   pft_pmeth, pmeth_gdd, pmeth_stg, pmeth_nvg
use module_pparams, only: tref
use module_sib, only: &
   poold_type, pooll_type, phen_type

implicit none

!...input variables
integer(i4), intent(in) :: gnum, lnum, pnum, pref
type(phen_param), intent(in) :: phencont
type(pool_param), intent(in) :: poolcont
type(phys_param), intent(in) :: physcont
real(r4), intent(in) :: daylenmax
real(r8), intent(in) :: daylen, daylendt
real(r8), intent(in) :: lai, tm

real(r8), intent(inout) :: vmax
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(phen_type), intent(inout) :: phent

!...local variables
integer(i4) :: ipft, ips
real(r8) :: tmdf

!...setup local variables
ipft = pref
tmdf = ((tm - tref)*1.8) + 32.0
phent%phenave_wac = 0.5 &
    * (phent%phenave_wac + phent%phenave_env)

!...setup growing season start variables
call phen_gss(phencont, daylenmax, daylen, daylendt, phent)

!...setup phenology stage index and phenology stage
if (phent%phen_istage .le. 0) then
   if (pft_pmeth(pnum) .eq. pmeth_stg) then
      !print('(a)'),'   Setting Up Dynamic Phenology'
      call phen_dynamic( phencont, &
           daylenmax, daylen, daylendt, lai, &
           phent)
   elseif (pft_pmeth(pnum) .eq. pmeth_gdd) then
      !print('(a)'),'   Setting Up Defined Phenology'
      call phen_defined( &
           gnum, lnum, ipft, pnum, &
           phencont, poolcont, &
           lai, tmdf, phent, pooldt, poollt, physcont)
   elseif (pft_pmeth(pnum) .ne. pmeth_nvg) then
          print('(a)'),'Unknown phenology method.'
          print('(a,i3)'),'Method Selected: ',pft_pmeth(pnum)
          print('(a,i3)'),'Dynamic Phenology Method == ',pmeth_stg
          print('(a,i3)'),'Defined Phenology Method == ',pmeth_gdd
          print('(a)'),'Stopping.'
          stop
   endif
endif

!...setup associated variables
ips = phent%phen_istage

!.....Photosynthetic Allocation
poollt%alloc_phen(:) = phencont%allocp(:,ips)
poollt%alloc(:) = poollt%alloc_phen(:)

!.....Dynamical Allocation
!.....not allowed during leaf-out or senescence
poollt%aadj_moist = phencont%adj_moist
poollt%aadj_temp = phencont%adj_temp
IF ((ips .le. 1) .or. (ips .eq. phencont%npstg)) THEN
    poollt%aadj_moist = .false.
    poollt%aadj_temp = .false.
ENDIF

!.....Canopy pool transfer rate
poollt%tfl_pstage = phencont%lptransfer(ips)

!.....Vmax
vmax = phencont%vmax0(ips)


end subroutine setup_phent
