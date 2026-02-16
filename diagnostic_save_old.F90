 !============================================================================
!
! Saves the variables and calculates time averages for output for:
!  -- single point (pbp)
!  -- entire domain, monthly (qp)
!  -- entire domain, daily/hourly (hr)
!
!============================================================================

!-----------------------------------------------------------
subroutine diagnostic_save_lall( &
     gnum, lnum, pnum, &
     soilt, cast, co2t, cost, &
     fract, fluxt, hydrost, hydrovt, &
     phent, pooldt, poollt, radt, &
     sift, sscolt, vegt)

use kinds
use module_io
use module_sib


implicit none

!...input variables
integer(i4), intent(in) :: gnum, lnum, pnum
type(soil_type), intent(in) :: soilt
type(cas_type), intent(in) :: cast
type(co2_type), intent(in) :: co2t
type(cos_type), intent(in) :: cost
type(fract_type), intent(in) :: fract
type(flux_type), intent(in) :: fluxt
type(hydros_type), intent(in) :: hydrost
type(hydrov_type), intent(in) :: hydrovt
type(phen_type), intent(in) :: phent
type(poold_type), intent(in) :: pooldt
type(pooll_type), intent(in) :: poollt
type(rad_type), intent(in) :: radt
type(sif_type), intent(in) :: sift
type(sscol_type), intent(in) :: sscolt
type(veg_type), intent(in) :: vegt
!type(pool_param), intent(in) :: poolcont

!...local variables
integer(i4) :: pbpref

!--------------------------------------------
IF (hr_saveluf) THEN
   call diagnostic_savelu(pnum, &
        hr_nvarlu, hr_nsaveperout, &
        hr_countsave, hr_vreflu, &
        soilt, cast, co2t, cost, &
        fract, fluxt, hydrost, hydrovt, &
        phent, pooldt, poollt, radt, &
        sift, sscolt, vegt, hr_lu(gnum,lnum,:,:))
ENDIF

IF (pbp_saveluf) THEN
    pbpref = pbp_outref(gnum)
    IF (pbpref .gt. izero) THEN
       call diagnostic_savelu(pnum, &
            pbp_nvarlu, pbp_nsaveperout, &
            pbp_countsave, pbp_vreflu, &
            soilt, cast, co2t, cost, &
            fract, fluxt, hydrost, hydrovt, &
            phent, pooldt, poollt, radt, &
            sift, sscolt, vegt, pbp_lu(pbpref,lnum,:,:))
    ENDIF
ENDIF

IF (qp_saveluf) THEN
    call diagnostic_savelu(pnum, &
         qp_nvarlu, qp_nsaveperout, &
         qp_countsave, qp_vreflu, &
         soilt, cast, co2t, cost, &
         fract, fluxt, hydrost, hydrovt, &
         phent, pooldt, poollt, radt, &
         sift, sscolt, vegt, qp_lu(gnum,lnum,:,:))
ENDIF

end subroutine diagnostic_save_lall

!-----------------------------------------------------------
subroutine diagnostic_save_gall( &
     gnum, sibg, gprogt)

use kinds
use module_io
use module_sib, only: &
   gridcell_type, gprog_type

implicit none

!...input variables
integer(i4), intent(in) :: gnum
type (gridcell_type), intent(in) :: sibg
type(gprog_type), intent(in) :: gprogt

!...local variables
integer(i4) :: pbpref

IF (hr_savegf) THEN
    call diagnostic_saveg( &
         hr_nvarg, hr_nsaveperout, &
         hr_vrefg, hr_countsave, &
         sibg, hr_g(gnum,:,:), gprogt)
ENDIF

IF (pbp_savegf) THEN
    pbpref = pbp_outref(gnum)
    IF (pbpref .gt. izero) THEN
        call diagnostic_saveg( &
              pbp_nvarg, pbp_nsaveperout, &
              pbp_vrefg, pbp_countsave, &
              sibg, pbp_g(pbpref,:,:), gprogt)
     ENDIF
ENDIF

IF (qp_savegf) THEN
    call diagnostic_saveg( &
         qp_nvarg, qp_nsaveperout, &
         qp_vrefg, qp_countsave, &
         sibg, qp_g(gnum,:,:), gprogt)
ENDIF

end subroutine diagnostic_save_gall

!-----------------------------------------------------------
subroutine diagnostic_savelu(pnum, &
     nvars, nsave, refsave, refvars, &
     soilt, cast, co2t, cost, fract, &
     fluxt, hydrost, hydrovt, &
     phent, pooldt, poollt, &
     radt, sift, sscolt, vegt, &
     outsave)

use kinds
use module_pparams, only: &
     mol_to_bu_mze, mol_to_bu_soy, &
     mol_to_bu_wwt, mol_to_dw, &
     mol_to_mg, mol_to_pmol, &
     mol_to_umol, mwc, &
     molc13_to_mg, mwc13
use module_pftinfo, only: &
     pft_mze, pft_soy, pft_wwt
use module_poolinfo
use module_sibconst, only: &
     npoolpft, nsnow, nsoil, ntot
use module_phosib
use module_fractsib
use module_sib
use module_time, only: dtisib


implicit none

!...input variables
integer(i4), intent(in) :: pnum
integer(i4), intent(in) :: nvars, nsave, refsave
integer(i4), dimension(nvars), intent(in) :: refvars
type(soil_type), intent(in) :: soilt
type(cas_type), intent(in) :: cast
type(co2_type), intent(in) :: co2t
type(cos_type), intent(in) :: cost
type(fract_type), intent(in) :: fract
type(flux_type), intent(in) :: fluxt
type(hydros_type), intent(in) :: hydrost
type(hydrov_type), intent(in) :: hydrovt
type(phen_type), intent(in)  :: phent
type(poold_type), intent(in) :: pooldt
type(pooll_type), intent(in) :: poollt
type(rad_type), intent(in) :: radt
type(sif_type), intent(in) :: sift
type(sscol_type), intent(in) :: sscolt
type(veg_type), intent(in) :: vegt
real(r8), intent(inout) :: outsave(nvars,nsave)
!type(pool_param), intent(in) :: poolcont

!...local variables
integer(i4) :: i, k
integer(i4) :: lp, frp, crp, wp, pp 
integer(i4) :: lpc13, frpc13, crpc13, wpc13, ppc13
integer(i4) :: cdb, l1p, l2p, slp, slow, arm
integer(i4) :: cdbc13, l1pc13, l2pc13, slpc13, slowc13, armc13

!--------------------------------
!...set local variables
!... matches set-up in equipools_control.F90
lp = pool_indx_leaf
frp  = pool_indx_froot
crp  = pool_indx_croot
wp   = pool_indx_stwd
pp   = pool_indx_prod
cdb  = pool_indx_cdb - npoolpft/2
l1p  = pool_indx_metl - npoolpft/2
l2p  = pool_indx_strl - npoolpft/2
slp  = pool_indx_slit - npoolpft/2
slow = pool_indx_slow - npoolpft/2
arm  = pool_indx_arm - npoolpft/2
lpc13 = pool_indx_leaf_c13-6
frpc13  = pool_indx_froot_c13-6
crpc13  = pool_indx_croot_c13-6
wpc13   = pool_indx_stwd_c13-6
ppc13   = pool_indx_prod_c13-6
cdbc13  = pool_indx_cdb_c13 - npoolpft
l1pc13  = pool_indx_metl_c13 - npoolpft
l2pc13  = pool_indx_strl_c13 - npoolpft
slpc13  = pool_indx_slit_c13 - npoolpft
slowc13 = pool_indx_slow_c13 - npoolpft
armc13  = pool_indx_arm_c13 - npoolpft

!--------------------------------
!...loop through output
k=nsoil
do i=1,nvars

  select case (refvars(i))

   !!!------------------------!!!
   !!!TIME-INVARIANT VARIABLES!!!
   !!!------------------------!!!

   !!!Soil Parameters!!!
   case (1)
      outsave(i,refsave) = outsave(i,refsave) + soilt%sandfrac
   case (2)
      outsave(i,refsave) = outsave(i,refsave) + soilt%clayfrac
   case (3)
      outsave(i,refsave) = outsave(i,refsave) + soilt%soref_vis
   case (4)
     outsave(i,refsave) = outsave(i,refsave) + soilt%soref_nir
   case (5)
      outsave(i,refsave) = outsave(i,refsave) + soilt%poros
   case (6)
      outsave(i,refsave) = outsave(i,refsave) + soilt%satco
   case (7)
      outsave(i,refsave) = outsave(i,refsave) + soilt%csolid
   case (8)
      outsave(i,refsave) = outsave(i,refsave) + soilt%tkdry
   case (9)
      outsave(i,refsave) = outsave(i,refsave) + soilt%tkmg
   case (10)
      outsave(i,refsave) = outsave(i,refsave) + soilt%tksat
   case (11)
      outsave(i,refsave) = outsave(i,refsave) + soilt%bee
   case (12)
      outsave(i,refsave) = outsave(i,refsave) + soilt%phsat
   case (13)
      outsave(i,refsave) = outsave(i,refsave) + soilt%fieldcap
   case (14)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + soilt%fc_eff(k)
   case (15)
      outsave(i,refsave) = outsave(i,refsave) + soilt%vwcmin
   case (16)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + soilt%wp_eff(k)
   case (17)
      outsave(i,refsave) = outsave(i,refsave) + soilt%wopt
   case (18)
      outsave(i,refsave) = outsave(i,refsave) + soilt%woptzm
   case (19)
      outsave(i,refsave) = outsave(i,refsave) + soilt%wsat
   case (20)
      outsave(i,refsave) = outsave(i,refsave) + soilt%zm

   !!!----------------------!!!
   !!!TIME-VARYING VARIABLES!!!
   !!!----------------------!!!

   !!!Climatological Variables!!!
   case (30)
       outsave(i,refsave) = outsave(i,refsave) + co2t%clim_assim*mol_to_umol
   case (31)
       outsave(i,refsave) = outsave(i,refsave) + vegt%clim_lai
   case (32)
       outsave(i,refsave) = outsave(i,refsave) + hydrovt%clim_pawfrw
   case (33)
       outsave(i,refsave) = outsave(i,refsave) + hydrovt%clim_tawfrw

   !!!Canopy and CAS Variables!!!
   case (48)
      outsave(i,refsave) = outsave(i,refsave) + cast%eacas
   case (49)
      outsave(i,refsave) = outsave(i,refsave) + cast%hcapc * dtisib
   case (50)
      outsave(i,refsave) = outsave(i,refsave) + cast%hcapcas * dtisib
   case (51)
      outsave(i,refsave) = outsave(i,refsave) + cast%shcas
   case (52)
      outsave(i,refsave) = outsave(i,refsave) + cast%tc
   case (53)
      outsave(i,refsave) = outsave(i,refsave) + cast%tcmin
   case (54)
      outsave(i,refsave) = outsave(i,refsave) + cast%tcas
   case (55)
      outsave(i,refsave) = outsave(i,refsave) + cast%tkecas
   case (56)
      outsave(i,refsave) = outsave(i,refsave) + cast%vcapcas * dtisib
   case (57)
      outsave(i,refsave) = outsave(i,refsave) + cast%vpd

   !!!CO2 Variables!!!
   case (60)
      outsave(i,refsave) = outsave(i,refsave) + co2t%assim*mol_to_umol
   case (61)
      outsave(i,refsave) = outsave(i,refsave) + co2t%assimd*mol_to_umol
   case (64)
      outsave(i,refsave) = outsave(i,refsave) + assim_omc*mol_to_umol
   case (65)
      outsave(i,refsave) = outsave(i,refsave) + assim_ome*mol_to_umol
   case (66)
      outsave(i,refsave) = outsave(i,refsave) + assim_oms*mol_to_umol
   case (67)
      outsave(i,refsave) = outsave(i,refsave) + co2t%assimpot*mol_to_umol
   case (68)
      outsave(i,refsave) = outsave(i,refsave) + assimpot_omc*mol_to_umol
   case (69)
      outsave(i,refsave) = outsave(i,refsave) + assimpot_ome*mol_to_umol
   case (70)
      outsave(i,refsave) = outsave(i,refsave) + assimpot_oms*mol_to_umol
   case (71)
      outsave(i,refsave) = outsave(i,refsave) + assimfac(1)
   case (72)
      outsave(i,refsave) = outsave(i,refsave) + assimfac(2)
   case (73)
      outsave(i,refsave) = outsave(i,refsave) + assimfac(3)
   case (74)
      outsave(i,refsave) = outsave(i,refsave) + assimfac(4)
   case (75)
      outsave(i,refsave) = outsave(i,refsave) + co2t%vmaxss*mol_to_umol
   case (76)
      outsave(i,refsave) = outsave(i,refsave) + co2t%apar*mol_to_umol
   case (77)
      outsave(i,refsave) = outsave(i,refsave) + co2t%aparkk
   case (78)
      outsave(i,refsave) = outsave(i,refsave) + co2t%gamma
   case (79)
      outsave(i,refsave) = outsave(i,refsave) + co2t%par*mol_to_umol
   case (80)
      outsave(i,refsave) = outsave(i,refsave) + co2t%nspar*mol_to_umol

   case (82)
      outsave(i,refsave) = outsave(i,refsave) + co2t%casd
   case (83)
      outsave(i,refsave) = outsave(i,refsave) + co2t%cflux*mol_to_umol
   case (84)
      outsave(i,refsave) = outsave(i,refsave) + co2t%pco2cas
   case (85)
      outsave(i,refsave) = outsave(i,refsave) + co2t%pco2c
   case (86)
      outsave(i,refsave) = outsave(i,refsave) + co2t%pco2i
   case (87)
      outsave(i,refsave) = outsave(i,refsave) + co2t%pco2s
   case (88)
      outsave(i,refsave) = outsave(i,refsave) + gah2o
   case (89)
      outsave(i,refsave) = outsave(i,refsave) + gbh2o
   case (90)
      outsave(i,refsave) = outsave(i,refsave) + gsh2o
   case (91)
      outsave(i,refsave) = outsave(i,refsave) + co2t%rst
   case (92)
      outsave(i,refsave) = outsave(i,refsave) + co2t%soilfrz
   case (93)
      outsave(i,refsave) = outsave(i,refsave) + co2t%soilfrztg
   case (94)
      outsave(i,refsave) = outsave(i,refsave) + co2t%soilfrztd
   case (95)
      outsave(i,refsave) = outsave(i,refsave) + co2t%rstfac(1)
   case (96)
      outsave(i,refsave) = outsave(i,refsave) + co2t%rstfac(2)
   case (97)
      outsave(i,refsave) = outsave(i,refsave) + co2t%rstfac(3)
   case (98)
      outsave(i,refsave) = outsave(i,refsave) + co2t%rstfac(4)

   !!!Carbonyl Sulfide (COS) Values!!!
   case (101)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_assim*mol_to_pmol
   case (102)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_flux*mol_to_pmol
   case (103)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_grnd*mol_to_pmol
   case (104)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_lru
   case (105)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_lru2
   case (106)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_lru3
   case (107)
      outsave(i,refsave) = outsave(i,refsave) + cost%cos_lru4
   case (108)
      outsave(i,refsave) = outsave(i,refsave) + cost%coscas*mol_to_pmol
   case (109)
      outsave(i,refsave) = outsave(i,refsave) + cost%cosi*mol_to_pmol
   case (110)
      outsave(i,refsave) = outsave(i,refsave) + cost%coss*mol_to_pmol
   case (111)
      outsave(i,refsave) = outsave(i,refsave) + cost%coscasp
   case (112)
      outsave(i,refsave) = outsave(i,refsave) + cost%cosgm
   case (113)
      outsave(i,refsave) = outsave(i,refsave) + cost%cosgt

   !!!Flux Variables!!!
   case (120)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%ct
   case (121)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%cu
   case (122)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%drag
   case (124)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%ustar
   case (125)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%ventmf
   case (126)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%ec * dtisib
   case (127)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%eci * dtisib
   case (128)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%ect * dtisib
   case (129)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%eg * dtisib
   case (130)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%egi * dtisib
   case (131)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%egs * dtisib
   case (132)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%egsmax * dtisib
   case (133)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%es * dtisib
   case (134)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%fws
   case (135)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%hc * dtisib
   case (136)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%hg * dtisib
   case (137)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%hs * dtisib
   case (138)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%fss
   case (139)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%storhc
   case (140)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%storhg
   case (141)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%ra
   case (142)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%rb
   case (143)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%rbc
   case (144)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%rc
   case (145)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%rd
   case (146)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%rdc
   case (147)
      outsave(i,refsave) = outsave(i,refsave) + fluxt%rds

   !!!Hydrological Soil Variables!!!
   case (151)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%rhsoil
   case (152)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%rsoil
   case (153)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%ecmass
   case (154)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%egmass
   case (155)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%infil
   case (156)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%p0
   case (157)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%pcpg_rain
   case (158)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%pcpg_snow
   case (159)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%roffo
   case (160)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%roff
   case (161)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%snow_gdepth
   case (162)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%snow_gmass
   case (163)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%snow_gvfc
   case (164)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%www_tot
   case (165)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%www_inflow
   case (166)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%satfrac
   case (167)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%capacc_liq
   case (168)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%capacc_snow
   case (169)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%capacg
   case (170)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%satcapc
   case (171)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%satcapg
   case (172)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%snow_cvfc
   case (173)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%wetfracc
   case (174)
      outsave(i,refsave) = outsave(i,refsave) + hydrost%wetfracg

   !!!Hydrological Vegetation Variables!!!
   case (180)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%paw_lay(k)
   case (181)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%pawmax_lay(k)
   case (182)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%pawfrac_lay(k)
   case (183)
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%pawfrw
   case (184)
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%pawftop
   case (185)
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%pawfzw
   case (186)
         if (k >= nsoil) k=0
         k=k+1
         outsave(i,refsave) = outsave(i,refsave) + hydrovt%taw_lay(k)
   case (187)
         if (k >= nsoil) k=0
         k=k+1
         outsave(i,refsave) = outsave(i,refsave) + hydrovt%tawfrac_lay(k)
   case (188)
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%tawfrw
   case (189)
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%tawftop
   case (190)
      outsave(i,refsave) = outsave(i,refsave) + hydrovt%tawfzw

   !!!Phenology Variables!!!
   case (201)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_assim*mol_to_umol
   case (202)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_assimsm*mol_to_umol
   case (203)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_assimpot
   case (204)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_pr
   case (205)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_prsm
   case (206)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_prsdoy
   case (207)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_prcdoy
   case (208)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_prpot
   case (209)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_tawftop
   case (210)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_tm
   case (211)
      if (phent%phenflag_assimlow) then
          outsave(i,refsave) = outsave(i,refsave) + done
      else
          outsave(i,refsave) = outsave(i,refsave) + dzero
       endif
   case (212)
      if (phent%phenflag_precip) then
          outsave(i,refsave) = outsave(i,refsave) + done
      else
          outsave(i,refsave) = outsave(i,refsave) + dzero
      endif
   case (213)
      if (phent%phenflag_gsspass) then
          outsave(i,refsave) = outsave(i,refsave) + done
      else
          outsave(i,refsave) = outsave(i,refsave) + dzero
      endif
   case (214)
      if (phent%phenflag_daylen) then
          outsave(i,refsave) = outsave(i,refsave) + done
      else
          outsave(i,refsave) = outsave(i,refsave) + dzero
      endif
   case (215)
      if (phent%phenflag_moist) then
          outsave(i,refsave) = outsave(i,refsave) + done
      else
          outsave(i,refsave) = outsave(i,refsave) + dzero
      endif
   case (216)
      if (phent%phenflag_temp) then
          outsave(i,refsave) = outsave(i,refsave) + done
      else
          outsave(i,refsave) = outsave(i,refsave) + dzero
      endif

   case (220)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%phen_istage)
   case (221)
      outsave(i,refsave) = outsave(i,refsave) + phent%phen_pi
   case (222)
      outsave(i,refsave) = outsave(i,refsave) + phent%phens_dayl
   case (223)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenc_climp
   case (224)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenc_laimax
   case (225)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenc_laimin
   case (226)
      outsave(i,refsave) = outsave(i,refsave) + phent%phens_grw
   case (227)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_env
   case (228)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_wa
   case (229)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_wac
   case (230)
      outsave(i,refsave) = outsave(i,refsave) + phent%phenave_wacsm
   case (231)
      outsave(i,refsave) = outsave(i,refsave) + phent%phens_wx

   case (232)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%nd_gs)
   case (233)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%nd_stg(1))
   case (234)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%nd_stg(2))
   case (235)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%nd_stg(3))
   case (236)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%nd_stg(4))
   case (237)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%nd_stg(5))
   case (238)
      outsave(i,refsave) = outsave(i,refsave) + poollt%nd_fire
   case (239)
      outsave(i,refsave) = outsave(i,refsave) + poollt%nd_grz
   case (240)
      outsave(i,refsave) = outsave(i,refsave) + phent%gdd
   case (241)
      outsave(i,refsave) = outsave(i,refsave) + dble(phent%ipd)
   case (242)
      outsave(i,refsave) = outsave(i,refsave) + phent%dapd
   case (243)
      outsave(i,refsave) = outsave(i,refsave) + phent%dapdaf
   case (244)  !above ground biomass
      outsave(i,refsave) = outsave(i,refsave) &
             + (poollt%poolpft(lp) + poollt%poolpft(wp) + &
                poollt%poolpft(pp)) * mwc
   case (245)  !brown carbon
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(cdb)*mwc
   case (246)  !seed carbon
      outsave(i,refsave) = outsave(i,refsave) + phent%seed_pool*mwc
   case (247)  !leaf carbon
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(lp)*mwc
   case (248)  !fine root carbon
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(frp)*mwc
   case (249)  !stem carbon
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(wp)*mwc
   case (250)  !product carbon
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(pp)*mwc
   case (251)  !brown dry weight
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(cdb)*mol_to_dw
   case (252)  !seed pool dry weight
      outsave(i,refsave) = outsave(i,refsave) + phent%seed_pool*mol_to_dw
   case (253)  !leaf dry weight
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(lp)*mol_to_dw
   case (254)  !fine root dry weight
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(frp)*mol_to_dw
   case (255)  !stem dry weight
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(wp)*mol_to_dw
   case (256)  !product dry weight
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(pp)*mol_to_dw
   case (257)
      outsave(i,refsave) = outsave(i,refsave) + phent%seed_pool*mol_to_dw
   case (258)
      outsave(i,refsave) = outsave(i,refsave) &
              + sum(poollt%gain_seed)*mol_to_umol
   case (259)  !yield
      if (pnum .eq. pft_mze) then
         outsave(i,refsave) = outsave(i,refsave) &
              + poollt%poolpft(pp)*mol_to_bu_mze
      elseif (pnum .eq. pft_soy) then
         outsave(i,refsave) = outsave(i,refsave) &
              + poollt%poolpft(pp)*mol_to_bu_soy
      elseif (pnum .eq. pft_wwt) then
         outsave(i,refsave) = outsave(i,refsave) &
              + poollt%poolpft(pp)*mol_to_bu_wwt
      else
         outsave(i,refsave) = outsave(i,refsave) &
              + poollt%poolpft(pp)*4.0
      endif

   case (270)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc(lp)
   case (271)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc(frp)
   case (272)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc(crp)
   case (273)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc(wp)
   case (274)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc(pp)
   case (275)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_phen(lp)
   case (276)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_phen(frp)
   case (277)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_phen(crp)
   case (278)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_phen(wp)
   case (279)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_phen(pp)
   case (280)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_moist(lp)
   case (281)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_moist(frp)
   case (282)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_moist(crp)
   case (283)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_moist(wp)
   case (284)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_moist(pp)
   case (285)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_temp(lp)
   case (286)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_temp(frp)
   case (287)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_temp(crp)
   case (288)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_temp(wp)
   case (289)
      outsave(i,refsave) = outsave(i,refsave) + poollt%alloc_temp(pp)

   !!!Pool Variables!!!
   !..Autotrophic Respiration/Transfer Scalars...
   case (299)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mcr_assim
   case (300)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mcr_freeze
   case (301)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mcr_hot
   case (302)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mcr_scale
   case (303)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_freeze_lay(k)
   case (304)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_hot_lay(k)
   case (305)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_scale_lay(k)
   case (306)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_assim
   case (307)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_freeze
   case (308)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_hot
   case (309)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_lai
   case (310)
      outsave(i,refsave) = outsave(i,refsave) + poollt%mrr_scale
   case (315)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tf_turnover(lp)
   case (316)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tf_turnover(frp)
   case (317)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tf_turnover(crp)
   case (318)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tf_turnover(wp)
   case (319)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tf_turnover(pp)
   case (320)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfl_daylen
   case (321)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfl_freeze
   case (322)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfl_dry
   case (323)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfl_pstage
   case (324)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfl_total

   !..Heterotrophic Respiration/Transfer Scalars...
   case (329)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_sfc_assim
   case (330)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_sfc_freeze
   case (331)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_sfc_hot
   case (332)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_sfc_precip
   case (333)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_sfc_scale
   case (336)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_freeze_lay(k)
   case (337)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_hot_lay(k)
   case (338)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_moist_lay(k)
   case (339)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_pawf_lay(k)
   case (340)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_scale_lay(k)
   case (341)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_assim
   case (342)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_freeze
   case (343)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_hot
   case (344)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_moist
   case (345)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_pawfrw
   case (346)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%mhrt_soil_scale


   !...Net Ecosystem Exchange
   case (350)  !nee
      outsave(i,refsave) = outsave(i,refsave) + &
              (poollt%resp_auto + pooldt%resp_het + poollt%resp_nveg &
             + poollt%resp_grz + poollt%resp_hrvst - co2t%assim)*mol_to_umol

   !...Respirations
   case (351)
     outsave(i,refsave) = outsave(i,refsave) + poollt%resp_auto*mol_to_umol
   case (352)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_grow*mol_to_umol
   case (353)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_leaf*mol_to_umol
   case (354)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_mntn*mol_to_umol
   case (355)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_root*mol_to_umol
   case (356)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%resp_het*mol_to_umol
   case (357)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_fire*mol_to_umol
   case (358)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_grz*mol_to_umol
   case (359)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_hrvst*mol_to_umol
   case (360)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_nveg*mol_to_umol
   case (361)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%resp_soil*mol_to_umol
   case (362)
      outsave(i,refsave) = outsave(i,refsave) &
           + (poollt%resp_auto + pooldt%resp_het + poollt%resp_nveg &
           + poollt%resp_grz + poollt%resp_hrvst)*mol_to_umol
   case (363)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rmmd_fire*mol_to_umol
   case (364)
      outsave(i,refsave) = outsave(i,refsave) + sum(poollt%loss_grz)*mol_to_umol
   case (365)
      outsave(i,refsave) = outsave(i,refsave) &
                            + sum(poollt%loss_hrvst_lay(1:5,:))*mol_to_umol
   case (366)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rmvd_hrvst*mol_to_umol


   !...Live Biomass...
   case (370)  !leaf pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(lp)*mol_to_mg
   case (371)  !leaf_npp
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(lp) + poollt%gain_seed(lp) &
              - poollt%loss_gresp(lp) - sum(poollt%loss_mresp_lay(lp,:)))
   case (372)  !leaf_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(lp) + poollt%gain_seed(lp))
   case (373)  !leaf_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%loss_gresp(lp) + sum(poollt%loss_mresp_lay(lp,:)) &
              + sum(poollt%loss_trans_lay(lp,:)) + poollt%loss_grz(1) &
              + sum(poollt%loss_hrvst_lay(lp,:)) + sum(poollt%loss_fire_lay(lp,:)))
   case (374)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_fire_lay(lp,1)*mol_to_umol
   case (375)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_grz(1)*mol_to_umol
   case (376)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_hrvst_lay(lp,:))
   case (377)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (sum(poollt%loss_mresp_lay(lp,:)) + poollt%loss_gresp(lp))
   case (378)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * poollt%loss_gresp(lp)
   case (379)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(lp,:))
   case (380)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_trans_lay(lp,:))

   case (381)  !fine root pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(frp)*mol_to_mg
   case (382)  !froot_npp
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(frp) + poollt%gain_seed(frp) &
              - poollt%loss_gresp(frp) - sum(poollt%loss_mresp_lay(frp,:)))
   case (383)  !froot_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(frp) + poollt%gain_seed(frp))
   case (384)  !froot_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%loss_gresp(frp) + sum(poollt%loss_mresp_lay(frp,:)) &
              + sum(poollt%loss_trans_lay(frp,:)) &
              + sum(poollt%loss_hrvst_lay(frp,:)) &
              + sum(poollt%loss_fire_lay(frp,:)))
   case (385)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(poollt%loss_fire_lay(frp,:))
   case (386)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_hrvst_lay(frp,:))
   case (387)
       outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (sum(poollt%loss_mresp_lay(frp,:)) + poollt%loss_gresp(frp))
   case (388)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * poollt%loss_gresp(frp)
   case (389)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(frp,:))
   case (390)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_trans_lay(frp,:))

   case (391)  !coarse root pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(crp)*mol_to_mg
   case (392)  !croot_npp
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(crp) + poollt%gain_seed(crp) &
              - poollt%loss_gresp(crp) - sum(poollt%loss_mresp_lay(crp,:)))
   case (393)  !croot_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(crp) + poollt%gain_seed(crp))
   case (394)  !croot_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%loss_gresp(crp) + sum(poollt%loss_mresp_lay(crp,:)) &
              + sum(poollt%loss_trans_lay(crp,:)) &
              + sum(poollt%loss_hrvst_lay(crp,:)) &
              + sum(poollt%loss_fire_lay(crp,:)))
   case (395)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(poollt%loss_fire_lay(crp,:))
   case (396)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_hrvst_lay(crp,:))
   case (397)
       outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (sum(poollt%loss_mresp_lay(crp,:)) + poollt%loss_gresp(crp))
   case (398)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * poollt%loss_gresp(crp)
   case (399)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(crp,:))
   case (400)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_trans_lay(crp,:))

   case (401)  !stem/wood pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(wp)*mol_to_mg
   case (402)  !stwd_npp
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(wp) + poollt%gain_seed(wp) &
              - poollt%loss_gresp(wp) - sum(poollt%loss_mresp_lay(wp,:)))
   case (403)  !stwd_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(wp) + poollt%gain_seed(wp))
   case (404)  !stwd_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%loss_gresp(wp) + sum(poollt%loss_mresp_lay(wp,:)) &
              + sum(poollt%loss_trans_lay(wp,:)) + poollt%loss_grz(2) &
              + sum(poollt%loss_hrvst_lay(wp,:)) + poollt%loss_fire_lay(wp,1))
   case (405)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_fire_lay(wp,1)*mol_to_umol
   case (406)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_grz(2)*mol_to_umol
   case (407)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_hrvst_lay(wp,:))
   case (408)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (sum(poollt%loss_mresp_lay(wp,:)) + poollt%loss_gresp(wp))
   case (409)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_gresp(wp)*mol_to_umol
   case (410)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(wp,:))
   case (411)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_trans_lay(wp,:))

   case (412)  !product pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(pp)*mol_to_mg
   case (413)  !prod_npp
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(pp) + poollt%gain_seed(pp) &
              - poollt%loss_gresp(pp) - sum(poollt%loss_mresp_lay(pp,:)))
   case (414)  !prod_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(pp) + poollt%gain_seed(pp))
   case (415)  !prod_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%loss_gresp(pp) + sum(poollt%loss_mresp_lay(pp,:)) &
              + sum(poollt%loss_trans_lay(pp,:)) + poollt%loss_grz(3) &
              + sum(poollt%loss_hrvst_lay(pp,:)) + poollt%loss_fire_lay(pp,1))
   case (416)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_fire_lay(pp,1)*mol_to_umol
   case (417)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_grz(3)*mol_to_umol
   case (418)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_hrvst_lay(pp,:))
   case (419)
       outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (sum(poollt%loss_mresp_lay(pp,:)) + poollt%loss_gresp(pp))
   case (420)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_gresp(pp)*mol_to_umol
   case (421)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(pp,:))
   case (422)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_trans_lay(pp,:))


   !!!Pool Dead Stores!!!
   case (430)  !coarse dead biomass pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(cdb)*mol_to_mg
   case (431)  !cdb_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
         * (sum(pooldt%gain_grz_lay(cdb,:)) + sum(pooldt%gain_hrvst_lay(cdb,:)) &
            + sum(pooldt%gain_transd_lay(cdb,:)) + sum(pooldt%gain_transl_lay(cdb,:)))
   case (432)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_grz_lay(cdb,:))
   case (433)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_hrvst_lay(cdb,:))
   case (434)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(cdb,:))
   case (435)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(cdb,:))
   case (436)  !cdb_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (sum(pooldt%loss_resp_lay(cdb,:)) &
           + sum(pooldt%loss_trans_lay(cdb,:)) &
           + sum(pooldt%loss_fire_lay(cdb,:)))
   case (437)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(pooldt%loss_fire_lay(cdb,:))
   case (438)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(cdb,:))
   case (439)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(cdb,:))

   case (442)  !metabolic litter pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(l1p)*mol_to_mg
   case (443)  !lmet_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
         * (sum(pooldt%gain_grz_lay(l1p,:)) + sum(pooldt%gain_hrvst_lay(l1p,:)) &
            + sum(pooldt%gain_transd_lay(l1p,:)) + sum(pooldt%gain_transl_lay(l1p,:)))
   case (444)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_grz_lay(l1p,:))
   case (445)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_hrvst_lay(l1p,:))
   case (446)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(l1p,:))
   case (447)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(l1p,:))
   case (448)  !lmet_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (sum(pooldt%loss_resp_lay(l1p,:)) &
           + sum(pooldt%loss_trans_lay(l1p,:)) &
           + sum(pooldt%loss_fire_lay(l1p,:)))
   case (449)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(pooldt%loss_fire_lay(l1p,:))
   case (450)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(l1p,:))
   case (451)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(l1p,:))

   case (453)  !structural litter pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(l2p)*mol_to_mg
   case (454)  !lstr_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
         * (sum(pooldt%gain_grz_lay(l2p,:)) + sum(pooldt%gain_hrvst_lay(l2p,:)) &
            + sum(pooldt%gain_transd_lay(l2p,:)) + sum(pooldt%gain_transl_lay(l2p,:)))
   case (455)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_grz_lay(l2p,:))
   case (456)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_hrvst_lay(l2p,:))
   case (457)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(l2p,:))
   case (458)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(l2p,:))
   case (459)  !lstr_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (sum(pooldt%loss_resp_lay(l2p,:)) &
           + sum(pooldt%loss_trans_lay(l2p,:))&
           + sum(pooldt%loss_fire_lay(l2p,:)))
   case (460)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(pooldt%loss_fire_lay(l2p,:))
   case (461)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(l2p,:))
   case (462)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(l2p,:))

   case (464)  !soil litter pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(slp)*mol_to_mg
   case (465)  !slit_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
         * (sum(pooldt%gain_grz_lay(slp,:)) + sum(pooldt%gain_hrvst_lay(slp,:)) &
            + sum(pooldt%gain_transd_lay(slp,:)) + sum(pooldt%gain_transl_lay(slp,:)))
   case (466)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_grz_lay(slp,:))
   case (467)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_hrvst_lay(slp,:))
   case (468)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(slp,:))
   case (469)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(slp,:))
   case (470)  !slit_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (sum(pooldt%loss_resp_lay(slp,:)) &
           + sum(pooldt%loss_trans_lay(slp,:)) &
           + sum(pooldt%loss_fire_lay(slp,:)))
   case (471)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(pooldt%loss_fire_lay(slp,:))
   case (472)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(slp,:))
   case (473)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(slp,:))

   case (475)  !slow pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(slow)*mol_to_mg
   case (476)  !slow_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
         * (sum(pooldt%gain_grz_lay(slow,:)) + sum(pooldt%gain_hrvst_lay(slow,:)) &
            + sum(pooldt%gain_transd_lay(slow,:)) + sum(pooldt%gain_transl_lay(slow,:)))
   case (477)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_grz_lay(slow,:))
   case (478)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_hrvst_lay(slow,:))
   case (479)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(slow,:))
   case (480)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(slow,:))
   case (481)  !slow_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (sum(pooldt%loss_resp_lay(slow,:)) &
           + sum(pooldt%loss_trans_lay(slow,:)) &
           + sum(pooldt%loss_fire_lay(slow,:)))
   case (482)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(pooldt%loss_fire_lay(slow,:))
   case (483)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(slow,:))
   case (484)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(slow,:))

   case (486)  !armored/passive pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(arm)*mol_to_mg
   case (487)  !arm_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
         * (sum(pooldt%gain_grz_lay(arm,:)) + sum(pooldt%gain_hrvst_lay(arm,:)) &
            + sum(pooldt%gain_transd_lay(arm,:)) + sum(pooldt%gain_transl_lay(arm,:)))
   case (488)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_grz_lay(arm,:))
   case (489)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_hrvst_lay(arm,:))
   case (490)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(arm,:))
   case (491)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(arm,:))
   case (492)  !arm_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (sum(pooldt%loss_resp_lay(arm,:)) &
           + sum(pooldt%loss_trans_lay(arm,:)) &
           + sum(pooldt%loss_fire_lay(arm,:)))
   case (493)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * sum(pooldt%loss_fire_lay(arm,:))
   case (494)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(arm,:))
   case (495)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(arm,:))

   !!!Radiation Values!!!
   case (510)
     outsave(i,refsave) = outsave(i,refsave) + radt%albedo_visb
   case (511)
     outsave(i,refsave) = outsave(i,refsave) + radt%albedo_visd
   case (512)
     outsave(i,refsave) = outsave(i,refsave) + radt%albedo_nirb
   case (513)
     outsave(i,refsave) = outsave(i,refsave) + radt%albedo_nird
   case (514)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacc(1,1)
   case (515)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacc(1,2)
   case (516)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacc(2,1)
   case (517)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacc(2,2)
   case (518)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacg(1,1)
   case (519)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacg(1,2)
   case (520)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacg(2,1)
   case (521)
     outsave(i,refsave) = outsave(i,refsave) + radt%radfacg(2,2)
   case (522)
      outsave(i,refsave) = outsave(i,refsave) + radt%radc3c
   case (523)
      outsave(i,refsave) = outsave(i,refsave) + radt%radc3g
   case (524)
      outsave(i,refsave) = outsave(i,refsave) + radt%radtc
   case (525)
      outsave(i,refsave) = outsave(i,refsave) + radt%radtg
   case (526)
      outsave(i,refsave) = outsave(i,refsave) + radt%radts
   case (527)
      outsave(i,refsave) = outsave(i,refsave) + radt%effgc
   case (528)
      outsave(i,refsave) = outsave(i,refsave) + radt%tsfc

   !!!Fluorescence (SIF) Values!!!
   case (550)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_je
   case (551)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_jo
   case (552)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_jejo
   case (553)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_x
   case (554)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_kd
   case (555)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_kn
   case (556)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif_kp
   case (557)
      outsave(i,refsave) = outsave(i,refsave) + sift%phi_d
   case (558)
      outsave(i,refsave) = outsave(i,refsave) + sift%phi_f
   case (559)
      outsave(i,refsave) = outsave(i,refsave) + sift%phi_n
   case (560)
      outsave(i,refsave) = outsave(i,refsave) + sift%phi_p
   case (561)
      outsave(i,refsave) = outsave(i,refsave) + sift%sif


   !!!Soil Column Variables!!!
   case (599)
      outsave(i,refsave) = outsave(i,refsave) + dble(sscolt%nsl)
   case (600)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%rootr(k)
   case (601)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(k)
   case (602)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%eff_poros(k)
   case (603)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%layer_z(k)
   case (604)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%node_z(k)
   case (605)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%shcap(k) * dtisib
   case (606)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%slamda(k)
   case (607)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%tksoil(k)
   case (608)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%vol_ice(k)
   case (609)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%vol_liq(k)
   case (610)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%dz(k)
   case (611)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(k)
   case (612)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_ice(k)
   case (613)
      if (k >= nsoil) k=-nsnow
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(k)

   case (620)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(1)
   case (621)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(2)
   case (622)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(3)
   case (623)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(4)
   case (624)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(5)
   case (625)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(6)
   case (626)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(7)
   case (627)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(8)
   case (628)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(9)
   case (629)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%td(10)
   case (630)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(1)
   case (631)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(2)
   case (632)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(3)
   case (633)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(4)
   case (634)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(5)
   case (635)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(6)
   case (636)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(7)
   case (637)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(8)
   case (638)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(9)
   case (639)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%www_liq(10)
   case (640)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(1)
   case (641)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(2)
   case (642)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(3)
   case (643)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(4)
   case (644)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(5)
   case (645)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(6)
   case (646)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(7)
   case (647)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(8)
   case (648)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(9)
   case (649)
      outsave(i,refsave) = outsave(i,refsave) + sscolt%satfrac_lay(10)

   !!!Vegetation Values !!!
   case (700)
      outsave(i,refsave) = outsave(i,refsave) + vegt%z0
   case (701)
      outsave(i,refsave) = outsave(i,refsave) + vegt%z0d
   case (702)
      outsave(i,refsave) = outsave(i,refsave) + vegt%zp_dispd
   case (703)
      outsave(i,refsave) = outsave(i,refsave) + vegt%zpd_adj
   case (704)
      outsave(i,refsave) = outsave(i,refsave) + vegt%zztemp
   case (705)
      outsave(i,refsave) = outsave(i,refsave) + vegt%zzwind
   case (706)
      outsave(i,refsave) = outsave(i,refsave) + vegt%cc1
   case (707)
      outsave(i,refsave) = outsave(i,refsave) + vegt%cc2
   case (708)
      if (k >= nsoil) k=0
      k=k+1
      outsave(i,refsave) = outsave(i,refsave) + vegt%rootf(k)
   case (709)
      outsave(i,refsave) = outsave(i,refsave) + vegt%fpar
   case (710)
      outsave(i,refsave) = outsave(i,refsave) + vegt%green
   case (711)
      outsave(i,refsave) = outsave(i,refsave) + vegt%lai
   case (712)
      outsave(i,refsave) = outsave(i,refsave) + vegt%lait
   case (713)
      outsave(i,refsave) = outsave(i,refsave) + vegt%vcover
   case (714)
      outsave(i,refsave) = outsave(i,refsave) + vegt%gmudmu
   case (715)
      outsave(i,refsave) = outsave(i,refsave) + vegt%park
   case (716)
      outsave(i,refsave) = outsave(i,refsave) + vegt%vmax*mol_to_umol
   case (717)
      outsave(i,refsave) = outsave(i,refsave) + vmaxts*mol_to_umol

    case(915)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_s(1)
    case(916)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_s(2)
    case(917)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_s(3)
    case(918)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_s(4)
    case(919)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_s(5)
    case(920)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_s(10)

    case(921)
       outsave(i,refsave) = outsave(i,refsave) + cost%gsh2onew
    case(922)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_grnd_Ogee
    case(923)
       outsave(i,refsave) = outsave(i,refsave) + cost%cos_soil
    case(924)
       outsave(i,refsave) = outsave(i,refsave) + cost%cosm
    case(925)
       outsave(i,refsave) = outsave(i,refsave) + cost%badcos

    !!!C-13 Values !!!
    case(930)
       outsave(i,refsave) = outsave(i,refsave) + fract%d13cca
    case(931)
       outsave(i,refsave) = outsave(i,refsave) + fract%d13cm
    case(932)
       outsave(i,refsave) = outsave(i,refsave) + fract%c13ca*mol_to_umol
    case(933)
       outsave(i,refsave) = outsave(i,refsave) + fract%c12ca*mol_to_umol
    case(934)
       outsave(i,refsave) = outsave(i,refsave) + fract%c13cm*mol_to_umol
    case(935)
       outsave(i,refsave) = outsave(i,refsave) + fract%c12cm*mol_to_umol
    case(936)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps
    case(937)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps_nog
    case(938)
       outsave(i,refsave) = outsave(i,refsave) + fract%d13cassim
    case(939)
       outsave(i,refsave) = outsave(i,refsave) + fract%d13cassim_nog
    case(940)
       outsave(i,refsave) = outsave(i,refsave) + fract%c13assim*mol_to_umol
    case(941)
       outsave(i,refsave) = outsave(i,refsave) + fract%c13assimd*mol_to_umol
    case(942)
       outsave(i,refsave) = outsave(i,refsave) + fract%c12assim*mol_to_umol

    case(943)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps_k1
    case(944)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps_k2
    case(945)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps_k3
    case(946)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps_k4
    case(947)
       outsave(i,refsave) = outsave(i,refsave) + fract%kiecps_k5
    
    case(948)
       outsave(i,refsave) = outsave(i,refsave) + fract%rcassim

   case (950)
      outsave(i,refsave) = outsave(i,refsave) + co2t%co2cas*mol_to_umol
   case (951)
      outsave(i,refsave) = outsave(i,refsave) + co2t%co2s*mol_to_umol
   case (952)
      outsave(i,refsave) = outsave(i,refsave) + co2t%co2i*mol_to_umol
   case (953)
      outsave(i,refsave) = outsave(i,refsave) + co2t%co2c*mol_to_umol
   case (954)
      outsave(i,refsave) = outsave(i,refsave) + co2t%co2m*mol_to_umol
   case (955)
      outsave(i,refsave) = outsave(i,refsave) + co2t%co2gamma*mol_to_umol


   !...Live C13 Biomass...
   case (960)  !leaf pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(lpc13)*molc13_to_mg
   case (961)  !fine root pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(frpc13)*molc13_to_mg
   case (962)  !coarse root pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(crpc13)*molc13_to_mg
   case (963)  !stem/wood pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(wpc13)*molc13_to_mg
   case (964)  !product pool
      outsave(i,refsave) = outsave(i,refsave) + poollt%poolpft(ppc13)*molc13_to_mg
   !...Dead C13 Biomass...
   case (965)  !coarse dead biomass pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(cdbc13)*molc13_to_mg
   case (966)  !metabolic litter pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(l1pc13)*molc13_to_mg
   case (967)  !structural litter pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(l2pc13)*molc13_to_mg
   case (968)  !soil litter pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(slpc13)*molc13_to_mg
   case (969)  !slow soil pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(slowc13)*molc13_to_mg
   case (970)  !armored/passive pool
      outsave(i,refsave) = outsave(i,refsave) + pooldt%poollu(armc13)*molc13_to_mg


   !...C13 Respirations
   case (971)
     outsave(i,refsave) = outsave(i,refsave) + poollt%resp_autoc13*mol_to_umol
   case (972)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_growc13*mol_to_umol
   case (973)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_leafc13*mol_to_umol
   case (974)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_mntnc13*mol_to_umol
   case (975)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_rootc13*mol_to_umol
   case (976)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%resp_hetc13*mol_to_umol
   case (977)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_firec13*mol_to_umol
   case (978)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_grzc13*mol_to_umol
   case (979)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_hrvstc13*mol_to_umol
   case (980)
      outsave(i,refsave) = outsave(i,refsave) + poollt%resp_nvegc13*mol_to_umol
   case (981)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%resp_soilc13*mol_to_umol
   case (982)
      outsave(i,refsave) = outsave(i,refsave) &
           + (poollt%resp_autoc13 + pooldt%resp_hetc13 + poollt%resp_nvegc13 &
           + poollt%resp_grzc13 + poollt%resp_hrvstc13)*mol_to_umol
   case (983)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rmmd_firec13*mol_to_umol
   case (984)
      outsave(i,refsave) = outsave(i,refsave) + sum(poollt%loss_grzc13)*mol_to_umol
   case (985)
      outsave(i,refsave) = outsave(i,refsave) &
                           + sum(poollt%loss_hrvst_lay(6:10,:))*mol_to_umol
   case (986)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rmvd_hrvstc13*mol_to_umol

   case (987)
      outsave(i,refsave) = outsave(i,refsave) + fract%c13resptot*mol_to_umol
   case (988)
      outsave(i,refsave) = outsave(i,refsave) + fract%c12resptot*mol_to_umol

   !...Isotope signatures: pools
   case (990)
       outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_leafc13
   case (991)
       outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_frootc13
   case (992)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_crootc13
   case (993)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_stwdc13
   case (994)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_prodc13
   case (995)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_cdbc13
   case (996)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_metlc13
   case (997)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_strlc13
   case (998)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_slitc13
   case (999)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_slowc13
   case (1000)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cpool_armc13

   !...Isotope signatures: respirations
   case (1001)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_autoc13
   case (1002)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_growc13
   case (1003)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_leafc13
   case (1004)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_mntnc13
   case (1005)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_rootc13
   case (1006)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_hetc13
   case (1007)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_firec13
   case (1008)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_grzc13
   case (1009)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_hrvstc13
   case (1010)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_nvegc13
   case (1011)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_soilc13
   case (1012)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cresp_totc13
   case (1013)
      outsave(i,refsave) = outsave(i,refsave) + fract%d13cemisfire_pool
   case (1014)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol* &
         (poollt%loss_fire_lay(lp,1) + poollt%loss_fire_lay(wp,1) &
          + poollt%loss_fire_lay(pp,1) &
          + sum(poollt%loss_fire_lay(crp,:)) + sum(poollt%loss_fire_lay(frp,:)) &
          + sum(pooldt%loss_fire_lay(cdb,:)) + sum(pooldt%loss_fire_lay(l1p,:)) &
          + sum(pooldt%loss_fire_lay(l2p,:)) + sum(pooldt%loss_fire_lay(slp,:)) &
          + sum(pooldt%loss_fire_lay(slow,:)) + sum(pooldt%loss_fire_lay(arm,:)))

   !...Additional diagnostics for pressure trace
   ! case(1013)
   !    outsave(i,refsave) = outsave(i,refsave) + fract%co2m_cfrax*mol_to_umol
   ! case(1014)
   !    outsave(i,refsave) = outsave(i,refsave) + cost%press_coscalc
   case(1015)
      outsave(i,refsave) = outsave(i,refsave) + co2t%eyy_phosib
   case(1016)
      outsave(i,refsave) = outsave(i,refsave) + gxco2
   case(1017)
      outsave(i,refsave) = outsave(i,refsave) + co2t%icconv_phosib
   case(1018)
      outsave(i,refsave) = outsave(i,refsave) + co2t%resp_casn*mol_to_umol
   !case(1019)
   !   outsave(i,refsave) = outsave(i,refsave) + cast%press_flxupdate
   !case(1020)
   !   outsave(i,refsave) = outsave(i,refsave) + cast%press_flxupdateps
   ! case(1021)
   !    outsave(i,refsave) = outsave(i,refsave) + co2t%press_phosibps
   ! case(1022)
   !    outsave(i,refsave) = outsave(i,refsave) + fract%press_cfraxps

   case (1020)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin(frp))
   case (1021)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin_updated(frpc13))
   case (1022)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin(crp))
   case (1023)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin_updated(crpc13))
   case (1024)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin(wp))
   case (1025)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin_updated(wpc13))

   !case (1019)
   !   outsave(i,refsave) = outsave(i,refsave) + poollt%isofactorp(lpc13)
   !case (1020)
   !   outsave(i,refsave) = outsave(i,refsave) + poollt%isofactorp(frpc13)
   !case (1021)
   !   outsave(i,refsave) = outsave(i,refsave) + poollt%isofactorp(crpc13)
   !case (1022)
   !   outsave(i,refsave) = outsave(i,refsave) + poollt%isofactorp(wpc13)

   !...Additional diagnostics for C13 pool trace
   case (1030)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * poollt%loss_gresp(lpc13)
   case (1031)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(lpc13,:))
   case (1032)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * poollt%loss_gresp(frpc13)
   case (1033)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(frpc13,:))
   case (1034)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * poollt%loss_gresp(crpc13)
   case (1035)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(crpc13,:))
   case (1036)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_gresp(wpc13)*mol_to_umol
   case (1037)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(wpc13,:))
   case (1038)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_gresp(ppc13)*mol_to_umol
   case (1039)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_mresp_lay(ppc13,:))
   case (1040)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%poolpft_lay(lp,:))
   case (1041)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%poolpft_lay(lpc13,:))

   case (1042)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin(lp))
   case (1043)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin(lpc13))
   case (1044)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%poolpft_dgain(lp,:))
   case (1045)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%poolpft_dgain(lpc13,:))
   case (1046)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%poolpft_dloss(lp,:))
   case (1047)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%poolpft_dloss(lpc13,:))
   case (1048)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpft(lp))
   case (1049)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpft(lpc13))
   case (1050)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_fire_lay(lpc13,1)*mol_to_umol
   case (1051)
      outsave(i,refsave) = outsave(i,refsave) + poollt%loss_grzc13(1)*mol_to_umol

   case (1052)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%loss_trans_lay(lpc13,:))
   case (1053)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfrac_lp
   case (1054)
      outsave(i,refsave) = outsave(i,refsave) + poollt%tfrac_lpc13
   case (1055)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%pftavail_lay(lp,:))
   case (1056)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%pftavail_lay(lpc13,:))

   case (1060)  !leaf_gain
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%gain_assim(lpc13) + poollt%gain_seed(lpc13))
   case (1061)  !leaf_loss
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
           * (poollt%loss_gresp(lpc13) + sum(poollt%loss_mresp_lay(lpc13,:)) &
              + sum(poollt%loss_trans_lay(lpc13,:)) + poollt%loss_grzc13(1) &
              + sum(poollt%loss_hrvst_lay(lpc13,:)) + sum(poollt%loss_fire_lay(lpc13,:)))

   case (1062)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%assim_dgain(lp,:))
   case (1063)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%assim_dloss(lp,:))
   case (1064)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%assim_dgain(lpc13,:))
   case (1065)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%assim_dloss(lpc13,:))
   case (1066)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%autoresp_dgain(lp,:))
   case (1067)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%autoresp_dloss(lp,:))
   case (1068)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            *  sum(poollt%autoresp_dgain(lpc13,:))
   case (1069)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%autoresp_dloss(lpc13,:))
   case (1070)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%autotran_dloss(lp,:))
   case (1071)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * sum(poollt%autotran_dloss(lpc13,:))
   case (1072)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
            * (poollt%poolpftmin_updated(lpc13))

   case (1073)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(cdbc13,:))
   case (1074)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(cdbc13,:))
   case (1075)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(cdbc13,:))
   case (1076)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(cdbc13,:))
   case (1077)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(l1pc13,:))
   case (1078)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(l1pc13,:))
   case (1079)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(l1pc13,:))
   case (1080)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(l1pc13,:))
   case (1081)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(l2pc13,:))
   case (1082)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(l2pc13,:))
   case (1083)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(l2pc13,:))
   case (1084)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(l2pc13,:))
   case (1085)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(slpc13,:))
   case (1086)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(slpc13,:))
   case (1087)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(slpc13,:))
   case (1088)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(slpc13,:))
   case (1089)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(slowc13,:))
   case (1090)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(slowc13,:))
   case (1091)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(slowc13,:))
   case (1092)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(slowc13,:))
   case (1093)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transl_lay(armc13,:))
   case (1094)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%gain_transd_lay(armc13,:))
   case (1095)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_resp_lay(armc13,:))
   case (1096)
      outsave(i,refsave) = outsave(i,refsave) + mol_to_umol &
          * sum(pooldt%loss_trans_lay(armc13,:))

   case (1100)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rcpoolpft(lpc13)
   case (1101)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rcpoolpft(wpc13)
   case (1102)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rcpoolpft(frpc13)
   case (1103)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rcpoolpft(crpc13)
   case (1104)
      outsave(i,refsave) = outsave(i,refsave) + poollt%rcpoolpft(ppc13)
   case (1105)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%rcpoollu(cdbc13)
   case (1106)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%rcpoollu(l1pc13)
   case (1107)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%rcpoollu(l2pc13)
   case (1108)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%rcpoollu(slpc13)
   case (1109)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%rcpoollu(slowc13)
   case (1110)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%rcpoollu(armc13)

   !case (1111)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%kiecpsxassim*mol_to_umol
   !case (1112)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%kiecpsk1xassim*mol_to_umol
   !case (1113)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%kiecpsk2xassim*mol_to_umol
   !case (1114)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%kiecpsk3xassim*mol_to_umol
   !case (1115)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%kiecpsk4xassim*mol_to_umol
   !case (1116)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%kiecpsk5xassim*mol_to_umol

   !case (1117)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cassimxassim*mol_to_umol 
   !case (1118)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crgrowxrgrow*mol_to_umol
   !case (1119)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crleafxrleaf*mol_to_umol
   !case (1120)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crmntnxrmntn*mol_to_umol
   !case (1121)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crrootxrroot*mol_to_umol
   !case (1122)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crfirexrfire*mol_to_umol
   !case (1123)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crgrzxrgrz*mol_to_umol
   !case (1124)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crhrvstxrhrvst*mol_to_umol
   !case (1125)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crnvegxrnveg*mol_to_umol
   !case (1126)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crhetxrhet*mol_to_umol
   !case (1127)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crsoilxrsoil*mol_to_umol
   !case (1128)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crtotxrtot*mol_to_umol
   !case (1129)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13crautoxrauto*mol_to_umol

   case (1130)
      outsave(i,refsave) = outsave(i,refsave) + poollt%curpoolpft(lp)
   case (1131)
      outsave(i,refsave) = outsave(i,refsave) + poollt%curpoolpft(wp)
   case (1132)
      outsave(i,refsave) = outsave(i,refsave) + poollt%curpoolpft(frp)
   case (1133)
      outsave(i,refsave) = outsave(i,refsave) + poollt%curpoolpft(crp)
   case (1134)
      outsave(i,refsave) = outsave(i,refsave) + poollt%curpoolpft(pp)
   case (1135)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%curpoollu(cdb)
   case (1136)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%curpoollu(l1p)
   case (1137)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%curpoollu(l2p)
   case (1138)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%curpoollu(slp)
   case (1139)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%curpoollu(slow)
   case (1140)
      outsave(i,refsave) = outsave(i,refsave) + pooldt%curpoollu(arm)

   !case (1141)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpleafxpleaf
   !case (1142)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpfrootxpfroot
   !case (1143)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpcrootxpcroot
   !case (1144)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpstwdxpstwd
   !case (1145)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpprodxpprod
   !case (1146)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpcdbxpcdb
   !case (1147)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpmetlxpmetl
   !case (1148)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpstrlxpstrl
   !case (1149)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpslitxpslit
   !case (1150)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cpslowxpslow
   !case (1151)
   !   outsave(i,refsave) = outsave(i,refsave) + fract%d13cparmxparm

  end select

enddo  !i=1,nvars

end subroutine diagnostic_savelu


!------------------------------------------------------------------
subroutine diagnostic_saveg( &
     nvars, nsave, refvars, refsave, &
     sibg, outsave, gprogt)

use kinds
use module_pftinfo, only: &
    pft_num, pft_type, type_crop, &
    pft_group, group_grass, group_shrub, &
    group_ndlfor, group_bdlfor
use module_pparams, only: &
    mol_to_pmol, mol_to_umol, &
    p0_sfc, tref, mol_to_mg
use module_sib, only: &
    gridcell_type, gprog_type
use module_phosib, only: &
    pressure
use module_poolinfo
use module_sibconst, only: &
    npoolpft

implicit none

!...input variables
integer(i4), intent(in) :: nvars, nsave, refsave
integer(i4), dimension(nvars), intent(in) :: refvars
type(gridcell_type), intent(in) :: sibg
real(r8), dimension(nvars,nsave), intent(inout) :: outsave

type(gprog_type), intent(in) :: gprogt

!...local variables
integer(byte) :: ptype
integer(byte) :: groupref
integer(i4) :: i,l
integer(i4) :: pnum
real(r8) :: assim, resp
real(r8) :: assimc13, d13cassimxassim
real(r8) :: femis, d13cfemisxemis, firec13, firetotc
real(r8) :: d13cresptotxresptot, resptotc13
real(r8) :: kiecpsxassim, kiecpsk1xassim, &
            kiecpsk2xassim, kiecpsk3xassim, &
            kiecpsk4xassim, kiecpsk5xassim
real(r8) :: totp, totlp, totdp
integer(i4) :: lp, frp, crp, wp, pp
integer(i4) :: cdb, l1p, l2p, slp, slow, arm
real(r8) :: nzero=1.E-14
logical :: isgrass, isshrub, isfor
real(r8) :: totpc13, totlpc13, totdpc13
integer(i4) :: lpc13, frpc13, crpc13, wpc13, ppc13
integer(i4) :: cdbc13, l1pc13, l2pc13, slpc13, slowc13, armc13
real(r8) :: fgra_agb,fgra_bgb,fshb_agb,fshb_bgb,ffor_agb,ffor_bgb

!--------------------------------
!...set local variables
!... matches set-up in equipools_control.F90
lp = pool_indx_leaf
frp  = pool_indx_froot
crp  = pool_indx_croot
wp   = pool_indx_stwd
pp   = pool_indx_prod
cdb  = pool_indx_cdb - npoolpft/2
l1p  = pool_indx_metl - npoolpft/2
l2p  = pool_indx_strl - npoolpft/2
slp  = pool_indx_slit - npoolpft/2
slow = pool_indx_slow - npoolpft/2
arm  = pool_indx_arm - npoolpft/2
lpc13 = pool_indx_leaf_c13-6
frpc13  = pool_indx_froot_c13-6
crpc13  = pool_indx_croot_c13-6
wpc13   = pool_indx_stwd_c13-6
ppc13   = pool_indx_prod_c13-6
cdbc13  = pool_indx_cdb_c13 - npoolpft
l1pc13  = pool_indx_metl_c13 - npoolpft
l2pc13  = pool_indx_strl_c13 - npoolpft
slpc13  = pool_indx_slit_c13 - npoolpft
slowc13 = pool_indx_slow_c13 - npoolpft
armc13  = pool_indx_arm_c13 - npoolpft

pressure = dble(gprogt%ps) * 100.0D0
!!!!!!!!!!!!!
do i=1,nvars

   select case (refvars(i))

   !!!Grid Cell Parameter Variables!!!
   case (800)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%cosz
   case (801)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%daylen
   case (802)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%daylendt
   case (803)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%tmd
   case (804)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%tmd - tref
   case (805)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%tmdf
   case (806)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%thm
   case (807)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%em
   case (808)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%ros
   case (809)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%psy
   case (810)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%radvbc
   case (811)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%radvdc
   case (812)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%radnbc
   case (813)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%radndc
   case (814)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%toa_solar
   case (815)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%toa_radvbc
   case (816)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%toa_radvdc
   case (817)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%toa_radnbc
   case (818)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%toa_radndc
   case (819)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%toa_par*mol_to_umol
   case (820)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%aod
   case (821)
       outsave(i,refsave) = outsave(i,refsave) + sibg%gdiagt%sif_atten


   !!!Grid Cell Prognostic Variables!!!
   case (850)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%sw_dwn
   case (851)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%dlwbot
   case (852)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%tm
   case (853)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%sh
   case (854)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%ps
   case (855)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%cupr*3600.
   case (856)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%lspr*3600.
   case (857)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%spdm
!   case (858)
!      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%pco2m*1.E6/p0_sfc
!   case (858)
!      outsave(i,refsave) = outsave(i,refsave) + dble(sibg%gprogt%pco2m)*1.0D6/pressure
   case (858)
      outsave(i,refsave) = outsave(i,refsave) + dble(sibg%gprogt%co2m)*mol_to_umol
   case (859)
      outsave(i,refsave) = outsave(i,refsave) + dble(sibg%gprogt%pco2m)
   case (860)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%pcosm*1.E12/pressure
   case (861)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%pcosm
   case (862)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%firec*mol_to_umol
   case (863)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%fireco2*mol_to_umol
   case (864)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%seas_precip
   case (866)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%seas_tm
   case (867)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%clim_cupr
   case (868)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%clim_precip
   case (869)
      outsave(i,refsave) = outsave(i,refsave) + sibg%gprogt%clim_tm


   !!!Grid-Cell Specialty Cases for Isotope code!!!
   case (870)  !kiecps
      assim = dzero
      kiecpsxassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         kiecpsxassim = kiecpsxassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * & 
                (sibg%l(l)%fract%kiecps)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + kiecpsxassim/assim   
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (871)  !kiecps_k1
      assim = dzero
      kiecpsk1xassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         kiecpsk1xassim = kiecpsk1xassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * &
                (sibg%l(l)%fract%kiecps_k1)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + kiecpsk1xassim/assim
 !     else
 !         outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (872)  !kiecps_k2
      assim = dzero
      kiecpsk2xassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         kiecpsk2xassim = kiecpsk2xassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * &
                (sibg%l(l)%fract%kiecps_k2)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + kiecpsk2xassim/assim
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (873)  !kiecps_k3
      assim = dzero
      kiecpsk3xassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         kiecpsk3xassim = kiecpsk3xassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * &
                (sibg%l(l)%fract%kiecps_k3)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + kiecpsk3xassim/assim
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (874)  !kiecps_k4
      assim = dzero
      kiecpsk4xassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         kiecpsk4xassim = kiecpsk4xassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * &
                (sibg%l(l)%fract%kiecps_k4)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + kiecpsk4xassim/assim
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (875)  !kiecps_k5
      assim = dzero
      kiecpsk5xassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         kiecpsk5xassim = kiecpsk5xassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * &
                (sibg%l(l)%fract%kiecps_k5)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + kiecpsk5xassim/assim
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (876)  !d13cassim
      assim = dzero
      d13cassimxassim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
         d13cassimxassim = d13cassimxassim + (mol_to_umol * &
                sibg%l(l)%co2t%assim * sibg%l(l)%larea) * &
                (sibg%l(l)%fract%d13cassim)
      enddo
      if (assim .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + d13cassimxassim/assim
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (877)  !d13cresptot
      resp = dzero
      d13cresptotxresptot = dzero
      do l=1,sibg%g_nlu
         resp = resp + mol_to_umol*sibg%l(l)%larea &
                  * (sibg%l(l)%pooldt%resp_het &
                   + sibg%l(l)%poollt%resp_auto &
                   + sibg%l(l)%poollt%resp_grz &
                   + sibg%l(l)%poollt%resp_hrvst &
                   + sibg%l(l)%poollt%resp_nveg )
         d13cresptotxresptot = d13cresptotxresptot + &
                   (mol_to_umol*sibg%l(l)%larea &
                   * (sibg%l(l)%pooldt%resp_het &
                   + sibg%l(l)%poollt%resp_auto &
                   + sibg%l(l)%poollt%resp_grz &
                   + sibg%l(l)%poollt%resp_hrvst &
                   + sibg%l(l)%poollt%resp_nveg )) &
                   * (sibg%l(l)%fract%d13cresp_totc13)
      enddo
      if (resp .gt. nzero) then
          outsave(i,refsave) = outsave(i,refsave) + d13cresptotxresptot/resp
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (878)  !gppc13
      assimc13 = dzero
      do l=1,sibg%g_nlu
         assimc13 = assimc13 + sibg%l(l)%fract%c13assim * mol_to_umol &
                 * sibg%l(l)%larea
      enddo
      outsave(i,refsave) = outsave(i,refsave) + assimc13

   case (879)  !recoc13
      resptotc13 = dzero
      do l=1,sibg%g_nlu
         resptotc13 = resptotc13 + sibg%l(l)%fract%c13resptot * mol_to_umol &
                 * sibg%l(l)%larea
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resptotc13

   case (880)  !firec13
      firec13 = dzero
      do l=1,sibg%g_nlu
         firec13 = firec13 + sibg%l(l)%poollt%resp_firec13 * mol_to_umol !&
                 !* sibg%l(l)%larea
      enddo
      outsave(i,refsave) = outsave(i,refsave) + firec13

   case (881)  !firetotc
      firetotc = dzero
      do l=1,sibg%g_nlu
         firetotc = firetotc + sibg%l(l)%poollt%resp_fire * mol_to_umol !&
                 !* sibg%l(l)%larea
      enddo
      outsave(i,refsave) = outsave(i,refsave) + firetotc

   case (882)  !d13cemis_fire
      femis = dzero
      d13cfemisxemis = dzero
      do l=1,sibg%g_nlu
        femis = femis + mol_to_umol* &
          ! sibg%l(l)%larea 
          (sibg%l(l)%poollt%loss_fire_lay(lp,1) &
          + sibg%l(l)%poollt%loss_fire_lay(wp,1) &
          + sibg%l(l)%poollt%loss_fire_lay(pp,1) &
          + sum(sibg%l(l)%poollt%loss_fire_lay(crp,:)) &
          + sum(sibg%l(l)%poollt%loss_fire_lay(frp,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(cdb,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(l1p,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(l2p,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(slp,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(slow,:)) & 
          + sum(sibg%l(l)%pooldt%loss_fire_lay(arm,:)))

        d13cfemisxemis = d13cfemisxemis + (mol_to_umol* &
          !sibg%l(l)%larea * 
          (sibg%l(l)%poollt%loss_fire_lay(lp,1) &
          + sibg%l(l)%poollt%loss_fire_lay(wp,1) &
          + sibg%l(l)%poollt%loss_fire_lay(pp,1) &
          + sum(sibg%l(l)%poollt%loss_fire_lay(crp,:)) &
          + sum(sibg%l(l)%poollt%loss_fire_lay(frp,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(cdb,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(l1p,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(l2p,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(slp,:)) &
          + sum(sibg%l(l)%pooldt%loss_fire_lay(slow,:)) & 
          + sum(sibg%l(l)%pooldt%loss_fire_lay(arm,:)))) &
          * (sibg%l(l)%fract%d13cemisfire_pool)
      enddo
      if (femis .gt. dzero) then
          outsave(i,refsave) = outsave(i,refsave) + d13cfemisxemis/femis
!      else
!          outsave(i,refsave) = outsave(i,refsave) + (-9999.0)
      endif

   case (883)  !total pools
      totp = dzero
      do l=1,sibg%g_nlu
      totp = totp + mol_to_mg * &
             sibg%l(l)%larea * (sibg%l(l)%poollt%poolpft(lp) &
             + sibg%l(l)%poollt%poolpft(wp) + sibg%l(l)%poollt%poolpft(pp) &
             + sibg%l(l)%poollt%poolpft(crp) + sibg%l(l)%poollt%poolpft(frp) &
             + sibg%l(l)%pooldt%poollu(cdb) + sibg%l(l)%pooldt%poollu(l1p) &
             + sibg%l(l)%pooldt%poollu(l2p) + sibg%l(l)%pooldt%poollu(slp) &
             + sibg%l(l)%pooldt%poollu(slow) + sibg%l(l)%pooldt%poollu(arm))
      enddo
      outsave(i,refsave) = outsave(i,refsave) + totp              

   case (884)  !live pools
      totlp = dzero
      do l=1,sibg%g_nlu
      totlp = totlp + mol_to_mg * &
             sibg%l(l)%larea * (sibg%l(l)%poollt%poolpft(lp) &
             + sibg%l(l)%poollt%poolpft(wp) + sibg%l(l)%poollt%poolpft(pp) &
             + sibg%l(l)%poollt%poolpft(crp) + sibg%l(l)%poollt%poolpft(frp))
      enddo
      outsave(i,refsave) = outsave(i,refsave) + totlp 

   case (885)  !dead pools
      totdp = dzero
      do l=1,sibg%g_nlu
      totdp = totdp + mol_to_mg * &
             sibg%l(l)%larea * & 
             (sibg%l(l)%pooldt%poollu(cdb) + sibg%l(l)%pooldt%poollu(l1p) &
             + sibg%l(l)%pooldt%poollu(l2p) + sibg%l(l)%pooldt%poollu(slp) &
             + sibg%l(l)%pooldt%poollu(slow) + sibg%l(l)%pooldt%poollu(arm))
      enddo
      outsave(i,refsave) = outsave(i,refsave) + totdp

   case (886)  !total C13 pools
      totpc13 = dzero
      do l=1,sibg%g_nlu
      totpc13 = totpc13 + mol_to_mg * &
             sibg%l(l)%larea * (sibg%l(l)%poollt%poolpft(lpc13) &
             + sibg%l(l)%poollt%poolpft(wpc13) + sibg%l(l)%poollt%poolpft(ppc13) &
             + sibg%l(l)%poollt%poolpft(crpc13) + sibg%l(l)%poollt%poolpft(frpc13) &
             + sibg%l(l)%pooldt%poollu(cdbc13) + sibg%l(l)%pooldt%poollu(l1pc13) &
             + sibg%l(l)%pooldt%poollu(l2pc13) + sibg%l(l)%pooldt%poollu(slpc13) &
             + sibg%l(l)%pooldt%poollu(slowc13) + sibg%l(l)%pooldt%poollu(armc13))
      enddo
      outsave(i,refsave) = outsave(i,refsave) + totpc13

   case (887)  !live C13 pools
      totlpc13 = dzero
      do l=1,sibg%g_nlu
      totlpc13 = totlpc13 + mol_to_mg * &
             sibg%l(l)%larea * (sibg%l(l)%poollt%poolpft(lpc13) &
             + sibg%l(l)%poollt%poolpft(wpc13) + sibg%l(l)%poollt%poolpft(ppc13) &
             + sibg%l(l)%poollt%poolpft(crpc13) + sibg%l(l)%poollt%poolpft(frpc13))
      enddo
      outsave(i,refsave) = outsave(i,refsave) + totlpc13

   case (888)  !dead C13 pools
      totdpc13 = dzero
      do l=1,sibg%g_nlu
      totdpc13 = totdpc13 + mol_to_mg * &
             sibg%l(l)%larea * &
             (sibg%l(l)%pooldt%poollu(cdbc13) + sibg%l(l)%pooldt%poollu(l1pc13) &
             + sibg%l(l)%pooldt%poollu(l2pc13) + sibg%l(l)%pooldt%poollu(slpc13) &
             + sibg%l(l)%pooldt%poollu(slowc13) + sibg%l(l)%pooldt%poollu(armc13))
      enddo
      outsave(i,refsave) = outsave(i,refsave) + totdpc13

   case (889)  !burned grass agb
      fgra_agb = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isgrass = (pft_group(pnum) .eq. group_grass)
         if (isgrass) then
          fgra_agb = fgra_agb + mol_to_umol* &
             !sibg%l(l)%larea * 
             (sibg%l(l)%poollt%loss_fire_lay(lp,1) &
             + sibg%l(l)%poollt%loss_fire_lay(wp,1) &
             + sibg%l(l)%poollt%loss_fire_lay(pp,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(cdb,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(l1p,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(l2p,1))
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + fgra_agb

   case (890)  !burned grass bgb
      fgra_bgb = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isgrass = (pft_group(pnum) .eq. group_grass)
         if (isgrass) then
          fgra_bgb = fgra_bgb + mol_to_umol* &
             !sibg%l(l)%larea * &
             (sum(sibg%l(l)%poollt%loss_fire_lay(crp,:)) &
             + sum(sibg%l(l)%poollt%loss_fire_lay(frp,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(slp,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(slow,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(arm,:)))
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + fgra_bgb

   case (891)  !burned shrub agb
      fshb_agb = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isshrub = (pft_group(pnum) .eq. group_shrub)
         if (isshrub) then
          fshb_agb = fshb_agb + mol_to_umol* &
             !sibg%l(l)%larea * 
             (sibg%l(l)%poollt%loss_fire_lay(lp,1) &
             + sibg%l(l)%poollt%loss_fire_lay(wp,1) &
             + sibg%l(l)%poollt%loss_fire_lay(pp,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(cdb,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(l1p,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(l2p,1))
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + fshb_agb

   case (892)  !burned shrub bgb
      fshb_bgb = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isshrub = (pft_group(pnum) .eq. group_shrub)
         if (isshrub) then
          fshb_bgb = fshb_bgb + mol_to_umol* &
             !sibg%l(l)%larea * &
             (sum(sibg%l(l)%poollt%loss_fire_lay(crp,:)) &
             + sum(sibg%l(l)%poollt%loss_fire_lay(frp,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(slp,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(slow,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(arm,:)))
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + fshb_bgb

   case (893)  !burned forest agb
      ffor_agb = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isfor = ((pft_group(pnum) .eq. group_ndlfor) .or. &
                  (pft_group(pnum) .eq. group_bdlfor))
         if (isfor) then
          ffor_agb = ffor_agb + mol_to_umol* &
             !sibg%l(l)%larea * 
             (sibg%l(l)%poollt%loss_fire_lay(lp,1) &
             + sibg%l(l)%poollt%loss_fire_lay(wp,1) &
             + sibg%l(l)%poollt%loss_fire_lay(pp,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(cdb,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(l1p,1) &
             + sibg%l(l)%pooldt%loss_fire_lay(l2p,1))
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + ffor_agb

   case (894)  !burned forest bgb
      ffor_bgb = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isfor = ((pft_group(pnum) .eq. group_ndlfor) .or. &
                  (pft_group(pnum) .eq. group_bdlfor))
         if (isfor) then
          ffor_bgb = ffor_bgb + mol_to_umol* &
             !sibg%l(l)%larea * &
             (sum(sibg%l(l)%poollt%loss_fire_lay(crp,:)) &
             + sum(sibg%l(l)%poollt%loss_fire_lay(frp,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(slp,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(slow,:)) &
             + sum(sibg%l(l)%pooldt%loss_fire_lay(arm,:)))
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + ffor_bgb

   case (895)  !gpp grasses
      assim = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isgrass = (pft_group(pnum) .eq. group_grass)
         if (isgrass) then
            assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                    * sibg%l(l)%larea
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + assim

   case (896)  !gpp shrubs
      assim = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isshrub = (pft_group(pnum) .eq. group_shrub)
         if (isshrub) then
            assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                    * sibg%l(l)%larea
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + assim

   case (897)  !gpp forests
      assim = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isfor = ((pft_group(pnum) .eq. group_ndlfor) .or. &
                  (pft_group(pnum) .eq. group_bdlfor))
         if (isfor) then
            assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                    * sibg%l(l)%larea
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + assim

   case (898)  !reco grasses
      resp = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isgrass = (pft_group(pnum) .eq. group_grass)
         if (isgrass) then
             resp = resp + mol_to_umol*sibg%l(l)%larea &
                    * (sibg%l(l)%pooldt%resp_het &
                     + sibg%l(l)%poollt%resp_auto &
                     + sibg%l(l)%poollt%resp_grz &
                     + sibg%l(l)%poollt%resp_hrvst &
                     + sibg%l(l)%poollt%resp_nveg)
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resp

   case (899)  !reco shrubs
      resp = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isshrub = (pft_group(pnum) .eq. group_shrub)
         if (isshrub) then
             resp = resp + mol_to_umol*sibg%l(l)%larea &
                    * (sibg%l(l)%pooldt%resp_het &
                     + sibg%l(l)%poollt%resp_auto &
                     + sibg%l(l)%poollt%resp_grz &
                     + sibg%l(l)%poollt%resp_hrvst &
                     + sibg%l(l)%poollt%resp_nveg)
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resp

   case (900)  !reco forests
      resp = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         isfor = ((pft_group(pnum) .eq. group_ndlfor) .or. &
                  (pft_group(pnum) .eq. group_bdlfor))
         if (isfor) then
             resp = resp + mol_to_umol*sibg%l(l)%larea &
                    * (sibg%l(l)%pooldt%resp_het &
                     + sibg%l(l)%poollt%resp_auto &
                     + sibg%l(l)%poollt%resp_grz &
                     + sibg%l(l)%poollt%resp_hrvst &
                     + sibg%l(l)%poollt%resp_nveg)
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resp

   !!!Grid-Cell Specialty Cases!!!
   case (901)  !gpp
      assim = dzero
      do l=1,sibg%g_nlu
         assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                 * sibg%l(l)%larea
      enddo
      outsave(i,refsave) = outsave(i,refsave) + assim

   case (902)  !gpp not crop
      assim = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         ptype = pft_type(pnum)
         if (ptype .ne. type_crop) then
            assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                    * sibg%l(l)%larea
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + assim

   case (903)  !gpp crop
      assim = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         ptype = pft_type(pnum)
         if (ptype .eq. type_crop) then
            assim = assim + sibg%l(l)%co2t%assim * mol_to_umol &
                    * sibg%l(l)%larea
         endif
      enddo

   case (904)  !reco
      resp = dzero
      do l=1,sibg%g_nlu
         resp = resp + mol_to_umol*sibg%l(l)%larea &
                  * (sibg%l(l)%pooldt%resp_het &
                   + sibg%l(l)%poollt%resp_auto &
                   + sibg%l(l)%poollt%resp_grz &
                   + sibg%l(l)%poollt%resp_hrvst &
                   + sibg%l(l)%poollt%resp_nveg )
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resp

   case (905)  !reco no crops
      resp = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         ptype = pft_type(pnum)
         if (ptype .ne. type_crop) then
             resp = resp + mol_to_umol*sibg%l(l)%larea &
                    * (sibg%l(l)%pooldt%resp_het &
                     + sibg%l(l)%poollt%resp_auto &
                     + sibg%l(l)%poollt%resp_grz &
                     + sibg%l(l)%poollt%resp_hrvst &
                     + sibg%l(l)%poollt%resp_nveg)
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resp

   case (906)  !reco crops
      resp = dzero
      do l=1,sibg%g_nlu
         pnum = pft_num(sibg%l(l)%ipft)
         ptype = pft_type(pnum)
         if (ptype .eq. type_crop) then
             resp = resp + mol_to_umol*sibg%l(l)%larea &
                    * (sibg%l(l)%pooldt%resp_het &
                     + sibg%l(l)%poollt%resp_auto &
                     + sibg%l(l)%poollt%resp_grz &
                     + sibg%l(l)%poollt%resp_hrvst &
                     + sibg%l(l)%poollt%resp_nveg)
         endif
      enddo
      outsave(i,refsave) = outsave(i,refsave) + resp

   case (907)  !nee
      resp = dzero
      assim = dzero
      do l=1,sibg%g_nlu
         resp = resp + mol_to_umol*sibg%l(l)%larea &
                * (sibg%l(l)%pooldt%resp_het &
                 + sibg%l(l)%poollt%resp_auto &
                 + sibg%l(l)%poollt%resp_grz &
                 + sibg%l(l)%poollt%resp_hrvst &
                 + sibg%l(l)%poollt%resp_nveg)
         assim = assim + sibg%l(l)%co2t%assim*mol_to_umol*sibg%l(l)%larea
      enddo
      outsave(i,refsave) = outsave(i,refsave) + (resp-assim)

   case (908)  !cos
      do l=1,sibg%g_nlu
         outsave(i,refsave) = outsave(i,refsave) &
                    + sibg%l(l)%cost%cos_flux &
                    * mol_to_pmol * sibg%l(l)%larea
      enddo

   case (909)  !sif_gome2
      if (sibg%gdiagt%sif_flag(1)) then
         do l=1,sibg%g_nlu
            outsave(i,refsave) = outsave(i,refsave) &
                    + sibg%l(l)%sift%sif * sibg%l(l)%larea
         enddo
      endif

   case (910)  !sif_oco2
      if (sibg%gdiagt%sif_flag(2)) then
         do l=1,sibg%g_nlu
            outsave(i,refsave) = outsave(i,refsave) &
                    + sibg%l(l)%sift%sif * sibg%l(l)%larea
         enddo
      endif

    case (911)  !Fcos_assim
       do l=1,sibg%g_nlu
          outsave(i,refsave) = outsave(i,refsave) &
                     + sibg%l(l)%cost%cos_assim &
                     * mol_to_pmol * sibg%l(l)%larea
       enddo

     case (912)  !Fcos_soil
        do l=1,sibg%g_nlu
           outsave(i,refsave) = outsave(i,refsave) &
                      + sibg%l(l)%cost%cos_soil &
                      * mol_to_pmol * sibg%l(l)%larea
        enddo

     case (913)  !Fcos_soil_Ogee
        do l=1,sibg%g_nlu
           outsave(i,refsave) = outsave(i,refsave) &
                      + sibg%l(l)%cost%cos_grnd_Ogee &
                      * mol_to_pmol * sibg%l(l)%larea
        enddo

   end select

enddo !i=1,nvars

end subroutine diagnostic_saveg
