!---------------------------------------------------------------------
subroutine sibtype_init()
!---------------------------------------------------------------------

!Initializes all variables in sibtype

use kinds
use module_param, only: &
    phencon
use module_pftinfo, only: &
    pft_num
use module_sib
use module_sibconst
use module_sibvs
use module_time, only: &
    startyear
use module_isodata, only: &
    isoyr, globc13, globc14

implicit none


! local variables
integer(i4) :: g,l,countl
integer(i4) :: gref,pnum,pref

!integer(i4), intent(in) :: cantot

   !...Initialize structure variables
   print*,''
   print*,'Initializing SiB Variables'

   !!!...Gridcell Variables
   if (allocated(sib%g)) deallocate(sib%g)
   allocate(sib%g(subcount))

   do g=1, subcount

       gref = subset(g)
       sib%g(g)%lat = latsib(gref)
       sib%g(g)%lon = lonsib(gref)

       call init_gdiagt(sib%g(g)%gdiagt)
       call init_gprogt(sib%g(g)%gprogt)

       !!!...Landunit/PFT Variables
       countl = sibvs(gref)%gnlu
       sib%g(g)%g_nlu = countl
       allocate(sib%g(g)%l(countl))
       sib%g(g)%l(:)%ipft = izero
       sib%g(g)%l(:)%larea = rzero

       do l=1,countl
!          cantot=npoolcan+npoolcanc13
          pref = sibvs(gref)%pftref(l)
          pnum = pft_num(pref)

          sib%g(g)%l(l)%ipft = pref
          sib%g(g)%l(l)%larea = sibvs(gref)%larea(l)

          call init_soilt(nsoil, sib%g(g)%l(l)%soilt)
          call init_cast(sib%g(g)%l(l)%cast)
          call init_co2t(sib%g(g)%l(l)%co2t)
          call init_cost(sib%g(g)%l(l)%cost, nsoil)
          call init_fract(sib%g(g)%l(l)%fract)
          call init_equibdt(npoollu,sib%g(g)%l(l)%equibdt)
          call init_equiblt(npoolpft,sib%g(g)%l(l)%equiblt)
          call init_fluxt(sib%g(g)%l(l)%fluxt)
          call init_hydrost(sib%g(g)%l(l)%hydrost)
          call init_hydrovt(nsoil, sib%g(g)%l(l)%hydrovt)
          call init_phent(phencon(pnum)%npstg, sib%g(g)%l(l)%phent)
          call init_pooldt(nsoil, npoollu, sib%g(g)%l(l)%pooldt)
          call init_poollt(nsoil, npoolcan, npoolcanc13, npoolcanc14, npoolpft, sib%g(g)%l(l)%poollt)
          call init_radt(sib%g(g)%l(l)%radt)
          call init_sift(sib%g(g)%l(l)%sift)
          call init_sscolt(nsnow, nsoil, sib%g(g)%l(l)%sscolt)
          call init_vegt(nsoil, sib%g(g)%l(l)%vegt)

       enddo  !l=1,countl
    enddo  !g=1,subcount

end subroutine sibtype_init

!--------------------------------------------------------
!  GRID CELL VARIABLES
!--------------------------------------------------------

!Routine to initialize the gridcell diagnostic variables
subroutine init_gdiagt(gpt)

use kinds
use module_sib, only: &
   gdiag_type

implicit none

!...input variables
type (gdiag_type), intent(inout) :: gpt

     gpt%cosz = dzero
     gpt%daylen = dzero
     gpt%daylendt = dzero

     gpt%tmdf = dzero
     gpt%thm = dzero
     gpt%bps(1:2) = dzero
     gpt%em = dzero
     gpt%ros = dzero
     gpt%psy = dzero

     gpt%radvbc = dzero
     gpt%radvdc = dzero
     gpt%radnbc = dzero
     gpt%radndc = dzero

     gpt%toa_solar = dzero
     gpt%toa_radvbc = dzero
     gpt%toa_radvdc = dzero
     gpt%toa_radnbc = dzero
     gpt%toa_radndc = dzero
     gpt%toa_par = dzero
     gpt%aod = dzero

     gpt%sif_atten = dzero
     gpt%sif_flag(:) = .false.

     gpt%gridcell_spunup = .false.

end subroutine init_gdiagt

!--------------------------------------------------------
!Routine to initialize the gridcell prognostic variables
subroutine init_gprogt(gpt)

use kinds
use module_sib, only: gprog_type

implicit none

!...input variables
type (gprog_type), intent(inout) :: gpt

     gpt%cupr = dzero
     gpt%cupr1 = dzero
     gpt%cupr2 = dzero
     gpt%cuprt = dzero
     gpt%dlwbot = dzero
     gpt%dlwbot1 = dzero
     gpt%dlwbot2 = dzero
     gpt%lspr = dzero
     gpt%lspr1 = dzero
     gpt%lspr2 = dzero
     gpt%lsprt = dzero
     gpt%ps = dzero
     gpt%ps1 = dzero
     gpt%ps2 = dzero
     gpt%sh = dzero
     gpt%sh1 = dzero
     gpt%sh2 = dzero
     gpt%spdm = dzero
     gpt%spdm1 = dzero
     gpt%spdm2 = dzero
     gpt%sw_dwn = dzero
     gpt%sw_dwn1 = dzero
     gpt%sw_dwn2 = dzero
     gpt%tm = dzero
     gpt%tm1 = dzero
     gpt%tm2 = dzero

     gpt%firec = dzero
     gpt%firec1 = dzero
     gpt%firec2 = dzero
     gpt%fireco2 = dzero
     gpt%fireco21 = dzero
     gpt%fireco22 = dzero
     
     gpt%pco2m = dzero
     gpt%pcosm = dzero
     gpt%co2m = dzero
     gpt%cosm = dzero
     !gpt%pressure_tm5 = dzero !ARACHO; to consider tm5 pressure

     gpt%tmd = dzero
     gpt%seas_precip = dzero
     gpt%seas_tm = dzero
     gpt%clim_cupr = dzero
     gpt%clim_precip = dzero
     gpt%clim_tm = dzero
     
end subroutine init_gprogt

!--------------------------------------------------------------
!  LAND UNIT TIME-INVARIANT VARIABLES
!--------------------------------------------------------------

!--------------------------------------------------------------
!Routine to initialize the soil characteristics/properties
subroutine init_soilt(nsoil, soilt)

use kinds
use module_sib, only: soil_type

implicit none

!...input variables
integer, intent(in) :: nsoil
type (soil_type), intent(inout) :: soilt

     soilt%sandfrac = dzero
     soilt%clayfrac = dzero
     soilt%soref_vis = dzero
     soilt%soref_nir = dzero

     soilt%poros = dzero
     soilt%satco = dzero

     soilt%csolid = dzero
     soilt%tkdry = dzero
     soilt%tkmg = dzero
     soilt%tksat = dzero

     soilt%bee = dzero
     soilt%phsat = dzero
     soilt%fieldcap = dzero
     soilt%vwcmin = dzero

     soilt%wopt = dzero
     soilt%woptzm = dzero
     soilt%wsat = dzero
     soilt%zm = dzero

     allocate(soilt%fc_eff(nsoil))
     soilt%fc_eff(:) = dzero
     allocate(soilt%wp_eff(nsoil))
     soilt%wp_eff(:) = dzero

end subroutine init_soilt


!--------------------------------------------------------------
!  LAND UNIT TIME-VARYING VARIABLES
!--------------------------------------------------------------

!--------------------------------------------------------------
!Routine to initialize the canopy air space (CAS) variables
subroutine init_cast(cast)

use kinds
use module_sib, only: cas_type

implicit none

!...input variables
type (cas_type), intent(inout) :: cast

     cast%tc = dzero

     cast%eacas = dzero
     cast%shcas = dzero
     cast%tcas = dzero
     cast%tkecas = dzero

     cast%hcapc   = dzero
     cast%tcmin   = dzero

     cast%hcapcas = dzero
     cast%thcas   = dzero
     cast%vcapcas = dzero

     cast%press_addinc = dzero
     cast%press_flxvrbrd = dzero
     cast%press_flxupdate = dzero
     cast%press_flxupdateps = dzero

end subroutine init_cast


!------------------------------------------------------------
!Routine to initialize the CO2 variables
subroutine init_co2t(co2t)

use kinds
use module_sib, only: co2_type

implicit none

!...input variables
type(co2_type), intent(inout) :: co2t

     co2t%assim = dzero
     co2t%assimd = dzero
     co2t%clim_assim = dzero

     co2t%assimpot = dzero
     co2t%apar = dzero
     co2t%aparkk = dzero
     co2t%gamma = dzero
     co2t%par = dzero
     co2t%nspar = dzero

     co2t%casd = dzero
     co2t%cflux = dzero

     co2t%pco2cas = dzero
     co2t%pco2c = dzero
     co2t%pco2i = dzero
     co2t%pco2s = dzero
     co2t%pco2m = dzero

     co2t%co2cas = dzero
     co2t%co2c = dzero
     co2t%co2i = dzero
     co2t%co2s = dzero
     co2t%co2m = dzero
     co2t%co2gamma = dzero

     co2t%rst = dzero

     co2t%soilfrz = dzero
     co2t%soilfrztg = dzero
     co2t%soilfrztd = dzero

     co2t%rstfac(:) = dzero
     co2t%vmaxss = dzero
     co2t%gs_stress(:) = dzero

     !co2t%press_phosib = dzero
     !co2t%press_phosibps = dzero

     co2t%icconv_phosib = izero
     co2t%eyy_phosib = dzero
     !co2t%icconv_phosib = izero
     co2t%resp_casn = dzero

end subroutine init_co2t

!--------------------------------------------------------
!Routine to initialize C13 variables
subroutine init_fract(fract)

use kinds
use module_sib, only: fract_type

implicit none

!...input variables
type(fract_type), intent(inout) :: fract

!.. Local variables
integer(i4) :: yrnow, i, loc

    do i=1,nisodatayr
      if (floor(isoyr(i)) .eq. startyear) then
       loc=i
       exit
      endif
    enddo
    fract%d13cm = dble(globc13(loc))
    fract%d14cm = dble(globc14(loc))

    fract%d13cca = dzero
    fract%c13ca = dzero
    fract%c12ca = dzero
    fract%c13cca = dzero
    fract%c12cca = dzero
    fract%c13cm = dzero
    fract%c12cm = dzero
    fract%kiecps = dzero
    fract%rcassim = dzero
    fract%rcassimfac = dzero
    fract%d13cassim = dzero
    fract%c13assim = dzero
    fract%kiecps_nog = dzero
    fract%rcassim_nog = dzero
    fract%d13cassim_nog = dzero
    fract%c13assim_nog = dzero
    fract%c13assimd = dzero
    fract%c12assim = dzero
    fract%c13resptot = dzero
    fract%c12resptot = dzero

    fract%c12assimn = dzero
    fract%c13assimn = dzero

    fract%d14cca = dzero
    fract%c14ca = dzero
    fract%c14cm = dzero
    fract%rcassimc14 = dzero
    fract%rcassimfacc14 = dzero
    fract%d14cassim = dzero
    fract%c14assim = dzero
    fract%c14assimd = dzero
    fract%c14resptot = dzero

    fract%c14alpha = dzero

    fract%kiecps_k1 = dzero
    fract%kiecps_k2 = dzero
    fract%kiecps_k3 = dzero
    fract%kiecps_k4 = dzero
    fract%kiecps_k5 = dzero

    fract%d13cpool_leafc13 = dzero
    fract%d13cpool_frootc13 = dzero
    fract%d13cpool_crootc13 = dzero
    fract%d13cpool_stwdc13 = dzero
    fract%d13cpool_prodc13 = dzero
    fract%d13cpool_cdbc13 = dzero
    fract%d13cpool_metlc13 = dzero
    fract%d13cpool_strlc13 = dzero
    fract%d13cpool_slitc13 = dzero
    fract%d13cpool_slowc13 = dzero
    fract%d13cpool_armc13 = dzero

    fract%d13cresp_autoc13 = dzero
    fract%d13cresp_growc13 = dzero
    fract%d13cresp_leafc13 = dzero
    fract%d13cresp_mntnc13 = dzero
    fract%d13cresp_rootc13 = dzero
    fract%d13cresp_hetc13 = dzero
    fract%d13cresp_firec13 = dzero
    fract%d13cemisfire_pool = dzero
    fract%d13cresp_grzc13 = dzero
    fract%d13cresp_hrvstc13 = dzero
    fract%d13cresp_nvegc13 = dzero
    fract%d13cresp_soilc13 = dzero
    fract%d13cresp_totc13 = dzero

    fract%d14cpool_leafc14 = dzero
    fract%d14cpool_frootc14 = dzero
    fract%d14cpool_crootc14 = dzero
    fract%d14cpool_stwdc14 = dzero
    fract%d14cpool_prodc14 = dzero
    fract%d14cpool_cdbc14 = dzero
    fract%d14cpool_metlc14 = dzero
    fract%d14cpool_strlc14 = dzero
    fract%d14cpool_slitc14 = dzero
    fract%d14cpool_slowc14 = dzero
    fract%d14cpool_armc14 = dzero

    fract%d14cresp_autoc14 = dzero
    fract%d14cresp_growc14 = dzero
    fract%d14cresp_leafc14 = dzero
    fract%d14cresp_mntnc14 = dzero
    fract%d14cresp_rootc14 = dzero
    fract%d14cresp_hetc14 = dzero
    fract%d14cresp_firec14 = dzero
    fract%d14cemisfire_pool = dzero
    fract%d14cresp_grzc14 = dzero
    fract%d14cresp_hrvstc14 = dzero
    fract%d14cresp_nvegc14 = dzero
    fract%d14cresp_soilc14 = dzero
    fract%d14cresp_totc14 = dzero

    fract%co2m_cfrax = dzero
    fract%press_cfrax = dzero
    fract%press_cfraxps = dzero

    fract%rcpoolfire = dzero
    fract%rcpoolfirec14 = dzero
    fract%poolemistotC = dzero
    fract%poolemisc13 = dzero
    fract%poolemisc14 = dzero
    !fract%d13cassimxassim = dzero
    !fract%kiecpsxassim = dzero
    !fract%kiecpsk1xassim = dzero
    !fract%kiecpsk2xassim = dzero
    !fract%kiecpsk3xassim = dzero
    !fract%kiecpsk4xassim = dzero
    !fract%kiecpsk5xassim = dzero

    !fract%d13crautoxrauto = dzero
    !fract%d13crgrowxrgrow = dzero
    !fract%d13crleafxrleaf = dzero
    !fract%d13crmntnxrmntn = dzero
    !fract%d13crrootxrroot = dzero
    !fract%d13crfirexrfire = dzero
    !fract%d13crgrzxrgrz = dzero
    !fract%d13crhrvstxrhrvst = dzero
    !fract%d13crnvegxrnveg = dzero
    !fract%d13crhetxrhet = dzero
    !fract%d13crsoilxrsoil = dzero
    !fract%d13crtotxrtot = dzero    

    !fract%d13cpleafxpleaf = dzero
    !fract%d13cpfrootxpfroot = dzero
    !fract%d13cpcrootxpcroot = dzero
    !fract%d13cpstwdxpstwd = dzero
    !fract%d13cpprodxpprod = dzero
    !fract%d13cpcdbxpcdb = dzero
    !fract%d13cpmetlxpmetl = dzero
    !fract%d13cpstrlxpstrl = dzero
    !fract%d13cpslitxpslit = dzero
    !fract%d13cpslowxpslow = dzero
    !fract%d13cparmxparm = dzero

!    allocate(fract%isoyr(nisodatayr))
!    fract%isoyr = dzero
!    allocate(fract%nhc14(nisodatayr))
!    fract%nhc14 = dzero
!    allocate(fract%shc14(nisodatayr))
!    fract%shc14 = dzero
!    allocate(fract%tropc14(nisodatayr))
!    fract%tropc14 = dzero
!    allocate(fract%globc13(nisodatayr))
!    fract%globc13 = dzero

end subroutine init_fract


!--------------------------------------------------------
!Routine to initialize the carbonyl sulfide (COS) variables
subroutine init_cost(cost, nsoil)

use kinds
use module_sib, only: cos_type

implicit none

!...input variables
type(cos_type), intent(inout) :: cost
integer(i4), intent(in) :: nsoil

     cost%cos_casd = dzero
     cost%cos_flux = dzero

     cost%coscas = dzero
     cost%coss = dzero
     cost%cosi = dzero
     cost%coscasp = dzero

     cost%cos_assim = dzero
     cost%cos_lru = dzero
     cost%cos_lru2 = dzero
     cost%cos_lru3 = dzero
     cost%cos_lru4 = dzero
     cost%cos_grnd = dzero
     cost%cosgm = dzero
     cost%cosgt = dzero

     cost%cos_grnd_Ogee = dzero
     cost%cos_soil = dzero

     cost%gsh2onew = dzero
     cost%cosm = dzero
     cost%badcos = dzero

     allocate(cost%cos_s(nsoil))
     cost%cos_s(:) = dzero

     cost%press_coscalc = dzero
     
end subroutine init_cost


!---------------------------------------------------------
!Routine to initialize the dead equilibrium variables
subroutine init_equibdt(npoollu, equibdt)

use kinds
use module_sib, only: equibd_type

implicit none

!...input variables
integer(i4), intent(in) :: npoollu
!integer(i4), intent(in) :: npoollu_c13
type(equibd_type), intent(inout) :: equibdt

     allocate(equibdt%poollu_totgain(npoollu))
     equibdt%poollu_totgain(:) = dzero
     allocate(equibdt%poollu_totloss(npoollu))
     equibdt%poollu_totloss(:) = dzero

!     allocate(equibdt%poollu_totgain(npoollu_c13))
!     equibdt%poollu_totgain(:) = dzero
!     allocate(equibdt%poollu_totloss(npoollu_c13))
!     equibdt%poollu_totloss(:) = dzero

     allocate(equibdt%poollu_init(npoollu))
     equibdt%poollu_init(:) = dzero
     allocate(equibdt%poollu_end(npoollu))
     equibdt%poollu_end(:) = dzero
     allocate(equibdt%poollu_min(npoollu))
     equibdt%poollu_min(:) = dzero
     allocate(equibdt%poollu_max(npoollu))
     equibdt%poollu_max(:) = dzero
     allocate(equibdt%poollu_gain(npoollu))
     equibdt%poollu_gain(:) = dzero
     allocate(equibdt%poollu_loss(npoollu))
     equibdt%poollu_loss(:) = dzero
     allocate(equibdt%poollu_ratio(npoollu))
     equibdt%poollu_ratio(:) = dzero
     allocate(equibdt%poollu_equib(npoollu))
     equibdt%poollu_equib(:) = dzero
     allocate(equibdt%poollu_notdone(npoollu))
     equibdt%poollu_notdone(:) = .true. 

     !.....variables for surface pools
     equibdt%deadsfc_init = dzero
     equibdt%deadsfc_end = dzero
     equibdt%deadsfc_gain = dzero
     equibdt%deadsfc_loss = dzero
     equibdt%deadsfc_ratio = dzero
     equibdt%deadsfc_notdone = .true.

     !.....variables for soil pools
     equibdt%deadsoil_init = dzero
     equibdt%deadsoil_end = dzero
     equibdt%deadsoil_gain = dzero
     equibdt%deadsoil_loss = dzero
     equibdt%deadsoil_ratio = dzero
     equibdt%deadsoil_notdone = .true.

     !.....variables for surface C13 pools
     equibdt%deadsfc_c13_init = dzero
     equibdt%deadsfc_c13_end = dzero
     equibdt%deadsfc_c13_gain = dzero
     equibdt%deadsfc_c13_loss = dzero
     equibdt%deadsfc_c13_ratio = dzero
     equibdt%deadsfc_c13_notdone = .true.

     !.....variables for soil C13 pools
     equibdt%deadsoil_c13_init = dzero
     equibdt%deadsoil_c13_end = dzero
     equibdt%deadsoil_c13_gain = dzero
     equibdt%deadsoil_c13_loss = dzero
     equibdt%deadsoil_c13_ratio = dzero
     equibdt%deadsoil_c13_notdone = .true.

     !.....variables for surface C14 pools
     equibdt%deadsfc_c14_init = dzero
     equibdt%deadsfc_c14_end = dzero
     equibdt%deadsfc_c14_gain = dzero
     equibdt%deadsfc_c14_loss = dzero
     equibdt%deadsfc_c14_ratio = dzero
     equibdt%deadsfc_c14_notdone = .true.

     !.....variables for soil C14 pools
     equibdt%deadsoil_c14_init = dzero
     equibdt%deadsoil_c14_end = dzero
     equibdt%deadsoil_c14_gain = dzero
     equibdt%deadsoil_c14_loss = dzero
     equibdt%deadsoil_c14_ratio = dzero
     equibdt%deadsoil_c14_notdone = .true.

     !.....variables for spin-up
     equibdt%lupft_spunup = .false.

end subroutine init_equibdt


!---------------------------------------------------------
!Routine to initialize the live equilibrium variables
subroutine init_equiblt(npoolpft, equiblt)

use kinds
use module_sib, only: equibl_type

implicit none

!....input variables
integer(i4), intent(in) :: npoolpft
type(equibl_type), intent(inout) :: equiblt

     allocate(equiblt%poolpft_totgain(npoolpft))
     equiblt%poolpft_totgain(:) = dzero
     allocate(equiblt%poolpft_totloss(npoolpft))
     equiblt%poolpft_totloss(:) = dzero

     allocate(equiblt%poolpft_init(npoolpft))
     equiblt%poolpft_init(:) = dzero
     allocate(equiblt%poolpft_end(npoolpft))
     equiblt%poolpft_end(:) = dzero
     allocate(equiblt%poolpft_min(npoolpft))
     equiblt%poolpft_min(:) = dzero
     allocate(equiblt%poolpft_max(npoolpft))
     equiblt%poolpft_max(:) = dzero
     allocate(equiblt%poolpft_gain(npoolpft))
     equiblt%poolpft_gain(:) = dzero
     allocate(equiblt%poolpft_loss(npoolpft))
     equiblt%poolpft_loss(:) = dzero
     allocate(equiblt%poolpft_ratio(npoolpft))
     equiblt%poolpft_ratio(:) = dzero
     allocate(equiblt%poolpft_equib(npoolpft))
     equiblt%poolpft_equib(:) = dzero
    
     allocate(equiblt%poolpft_notdone(npoolpft))
     equiblt%poolpft_notdone(:) = .true.

     !.....variables for live pools
     equiblt%live_init = dzero
     equiblt%live_end = dzero
     equiblt%live_gain = dzero
     equiblt%live_loss = dzero
     equiblt%live_ratio = dzero
     equiblt%live_notdone = .true.

     !.....variables for live pools (C13)
     equiblt%live_c13_init = dzero
     equiblt%live_c13_end = dzero
     equiblt%live_c13_gain = dzero
     equiblt%live_c13_loss = dzero
     equiblt%live_c13_ratio = dzero
     equiblt%live_c13_notdone = .true.

     !.....variables for live pools (C14)
     equiblt%live_c14_init = dzero
     equiblt%live_c14_end = dzero
     equiblt%live_c14_gain = dzero
     equiblt%live_c14_loss = dzero
     equiblt%live_c14_ratio = dzero
     equiblt%live_c14_notdone = .true.

end subroutine init_equiblt


!--------------------------------------------------------------
!Routine to initialize the flux variables
subroutine init_fluxt(fluxt)

use kinds
use module_sib, only: flux_type

implicit none

!...input variables
type (flux_type), intent(inout) :: fluxt

    fluxt%ct = dzero
    fluxt%cu = dzero
    fluxt%drag = dzero
    fluxt%ustar = dzero
    fluxt%ventmf = dzero

    fluxt%ec = dzero
    fluxt%eci = dzero
    fluxt%ect = dzero
    fluxt%eg = dzero
    fluxt%egi = dzero
    fluxt%egs = dzero
    fluxt%egsmax = dzero
    fluxt%es = dzero
    fluxt%fws = dzero

    fluxt%hc = dzero
    fluxt%hg = dzero
    fluxt%hs = dzero
    fluxt%fss = dzero
    fluxt%storhc = dzero
    fluxt%storhg = dzero

    fluxt%ra = dzero
    fluxt%rb = dzero
    fluxt%rbc = dzero
    fluxt%rc = dzero
    fluxt%rd = dzero
    fluxt%rdc = dzero
    fluxt%rds = dzero

    fluxt%ebalnum = izero
    fluxt%wbalnum = izero

    fluxt%et0 = dzero
    fluxt%et1 = dzero
    fluxt%et0a = dzero
    fluxt%et1a = dzero

end subroutine init_fluxt


!--------------------------------------------------------------
!Routine to initialize the soil hydrological variables
subroutine init_hydrost(hydrost)

use kinds
use module_sib, only: hydros_type

implicit none

!...input variables
type (hydros_type), intent(inout) :: hydrost

     hydrost%rhsoil = dzero
     hydrost%rsoil = dzero

     hydrost%ecmass = dzero
     hydrost%egmass = dzero
     hydrost%infil = dzero
     hydrost%p0 = dzero
     hydrost%pcpg_rain = dzero
     hydrost%pcpg_snow = dzero
     hydrost%roff = dzero
     hydrost%roffo = dzero
     hydrost%snow_gdepth = dzero
     hydrost%snow_gmass = dzero
     hydrost%snow_gvfc  = dzero

     hydrost%www_tot = dzero
     hydrost%www_inflow = dzero
     hydrost%satfrac = dzero

     hydrost%capacc_liq = dzero
     hydrost%capacc_snow = dzero
     hydrost%capacg = dzero
     hydrost%satcapc = dzero
     hydrost%satcapg = dzero
     hydrost%snow_cvfc = dzero
     hydrost%wetfracc = dzero
     hydrost%wetfracg = dzero

end subroutine init_hydrost


!--------------------------------------------------------------
!Routine to initialize the vegetation hydrological variables
subroutine init_hydrovt(nsoil, hydrovt)

use kinds
use module_sib, only: hydrov_type

implicit none

!...input variables
integer(i4), intent(in) :: nsoil
type (hydrov_type), intent(inout) :: hydrovt

     allocate(hydrovt%paw_lay(nsoil))
     hydrovt%paw_lay(:) = dzero
     allocate(hydrovt%pawmax_lay(nsoil))
     hydrovt%pawmax_lay(:) = dzero
     allocate(hydrovt%pawfrac_lay(nsoil))
     hydrovt%pawfrac_lay(:) = dzero

     hydrovt%pawfrw = dzero
     hydrovt%pawftop = dzero
     hydrovt%pawfzw = dzero

     allocate(hydrovt%taw_lay(nsoil))
     hydrovt%taw_lay(:) = dzero
     allocate(hydrovt%tawfrac_lay(nsoil))
     hydrovt%tawfrac_lay(:) = dzero

     hydrovt%tawfrw = dzero
     hydrovt%tawftop = dzero
     hydrovt%tawfzw = dzero

     hydrovt%clim_pawfrw = dzero
     hydrovt%clim_tawfrw = dzero

end subroutine init_hydrovt

!--------------------------------------------------------
!Routine to initialize the phenology variables
subroutine init_phent(npstg, phent)

use kinds
use module_sib, only: phen_type

implicit none

!...input variables
integer(i4), intent(in) :: npstg
type(phen_type), intent(inout) :: phent

     phent%phenave_assim = dzero
     phent%phenave_assimsm = dzero
     phent%phenave_assimpot = dzero
     phent%phenflag_assimlow = .false.

     phent%phenave_pr = dzero
     phent%phenave_prsm = dzero
     phent%phenave_prsdoy = dzero
     phent%phenave_prcdoy = dzero
     phent%phenave_prpot = dzero
     phent%phenflag_precip = .false.

     phent%phenave_tawftop = dzero
     phent%phenflag_moist = .false.
     
     phent%phenave_tm = dzero
     phent%phenflag_temp = .false.
     
     phent%phenflag_daylen = .false.
     phent%phenflag_gsspass = .false.

     phent%nd_dormant = izero
     phent%nd_gs = izero
     allocate(phent%nd_stg(npstg))
     phent%nd_stg(:) = izero

     phent%phen_istage = ione
     phent%phen_pi = done
     phent%phens_dayl = done

     phent%phenc_climp = dzero
     phent%phenc_laimax = dzero
     phent%phenc_laimin = dzero
     phent%phens_grw = dzero

     phent%phenave_env = dzero
     phent%phenave_wa = dzero
     phent%phenave_wac = dzero
     phent%phenave_wacsm = dzero
     phent%phens_wx = dzero

     phent%ipd = izero
     phent%dapd = izero
     phent%dapdaf = izero
     phent%gdd = dzero

     phent%seed_pool = dzero

end subroutine init_phent


!--------------------------------------------------------------
!Routine to initialize the dead pool variables
subroutine init_pooldt(nsoil, npoollu, pooldt)

use kinds
use module_sib, only: poold_type

implicit none

!...input variables
integer(i4), intent(in) :: nsoil, npoollu
type (poold_type), intent(inout) :: pooldt

     allocate(pooldt%gain_grz_lay(npoollu,nsoil))
     pooldt%gain_grz_lay(:,:) = dzero

     allocate(pooldt%gain_hrvst_lay(npoollu,nsoil))
     pooldt%gain_hrvst_lay(:,:) = dzero

     allocate(pooldt%gain_transl_lay(npoollu,nsoil))
     pooldt%gain_transl_lay(:,:) = dzero

     allocate(pooldt%gain_transd_lay(npoollu,nsoil))
     pooldt%gain_transd_lay(:,:) = dzero

     allocate(pooldt%loss_fire_lay(npoollu,nsoil))
     pooldt%loss_fire_lay(:,:) = dzero
     
     pooldt%mhrt_sfc_assim = done
     pooldt%mhrt_sfc_freeze = done
     pooldt%mhrt_sfc_hot = done
     pooldt%mhrt_sfc_precip = done
     pooldt%mhrt_sfc_scale = done

     allocate(pooldt%mhrt_soil_freeze_lay(nsoil))
     allocate(pooldt%mhrt_soil_hot_lay(nsoil))
     allocate(pooldt%mhrt_soil_moist_lay(nsoil))
     allocate(pooldt%mhrt_soil_pawf_lay(nsoil))
     allocate(pooldt%mhrt_soil_scale_lay(nsoil))
     pooldt%mhrt_soil_freeze_lay(:) = done
     pooldt%mhrt_soil_hot_lay(:) = done
     pooldt%mhrt_soil_moist_lay(:) = done
     pooldt%mhrt_soil_pawf_lay(:) = done
     pooldt%mhrt_soil_scale_lay(:) = done
     pooldt%mhrt_soil_assim = done
     pooldt%mhrt_soil_freeze = done
     pooldt%mhrt_soil_hot = done
     pooldt%mhrt_soil_moist = done
     pooldt%mhrt_soil_pawfrw = done
     pooldt%mhrt_soil_scale = done

     allocate(pooldt%kratert_lay(npoollu,nsoil))
     pooldt%kratert_lay(:,:) = dzero
     allocate(pooldt%loss_resp_lay(npoollu,nsoil))
     pooldt%loss_resp_lay(:,:) = dzero
     allocate(pooldt%loss_trans_lay(npoollu,nsoil))
     pooldt%loss_trans_lay(:,:) = dzero

     pooldt%resp_het = dzero
     pooldt%resp_soil = dzero
     allocate(pooldt%resp_soil_lay(nsoil))
     pooldt%resp_soil_lay(:) = dzero
     pooldt%resp_soilnr = dzero
     allocate(pooldt%resp_soilnr_lay(nsoil))
     pooldt%resp_soilnr_lay(:) = dzero

     pooldt%resp_hetc13 = dzero
     pooldt%resp_soilc13 = dzero
     allocate(pooldt%resp_soilc13_lay(nsoil))
     pooldt%resp_soilc13_lay(:) = dzero
     pooldt%resp_soilnrc13 = dzero
     allocate(pooldt%resp_soilnrc13_lay(nsoil))
     pooldt%resp_soilnrc13_lay(:) = dzero

     pooldt%resp_hetc14 = dzero
     pooldt%resp_soilc14 = dzero
     allocate(pooldt%resp_soilc14_lay(nsoil))
     pooldt%resp_soilc14_lay(:) = dzero
     pooldt%resp_soilnrc14 = dzero
     allocate(pooldt%resp_soilnrc14_lay(nsoil))
     pooldt%resp_soilnrc14_lay(:) = dzero

     allocate(pooldt%poollu_dgain(npoollu,nsoil))
     pooldt%poollu_dgain(:,:) = dzero
     allocate(pooldt%poollu_dloss(npoollu,nsoil))
     pooldt%poollu_dloss(:,:) = dzero
     allocate(pooldt%poollu(npoollu))
     pooldt%poollu(:) = dzero
     allocate(pooldt%poollu_lay(npoollu,nsoil))
     pooldt%poollu_lay(:,:) = dzero
     allocate(pooldt%poollu_flay(npoollu,nsoil))
     pooldt%poollu_flay(:,:) = dzero

     allocate(pooldt%poollup(npoollu))
     pooldt%poollup(:) = dzero

     allocate(pooldt%rcpoollu(npoollu))
     pooldt%rcpoollu(:) = dzero
     allocate(pooldt%rcpoollu_lay(npoollu,nsoil))
     pooldt%rcpoollu_lay(:,:) = dzero
     allocate(pooldt%curpoollu(npoollu))
     pooldt%curpoollu(:) = dzero

end subroutine init_pooldt


!--------------------------------------------------------
!Routine to initialize the live pool variables
subroutine init_poollt(nsoil, npoolcan, npoolcanc13, npoolcanc14, npoolpft, poollt)

use kinds
use module_sib, only: pooll_type

implicit none

!...input variables
integer(i4), intent(in) :: nsoil, npoolcan, npoolcanc13, npoolcanc14, npoolpft
type(pooll_type), intent(inout) :: poollt

     allocate(poollt%alloc(npoolpft))
     poollt%alloc(:) = dzero
     poollt%aadj_moist = .false.
     poollt%aadj_temp = .false.

     allocate(poollt%alloc_phen(npoolpft))
     poollt%alloc_phen(:) = dzero
     allocate(poollt%alloc_moist(npoolpft))
     poollt%alloc_moist(:) = dzero
     allocate(poollt%alloc_temp(npoolpft))
     poollt%alloc_temp(:) = dzero

     allocate(poollt%gain_assim(npoolpft))
     poollt%gain_assim(:) = dzero
     allocate(poollt%gain_seed(npoolpft))
     poollt%gain_seed(:) = dzero

     allocate(poollt%loss_gresp(npoolpft))
     poollt%loss_gresp(:) = dzero

     poollt%mcr_assim  = done
     poollt%mcr_freeze = done
     poollt%mcr_hot    = done
     poollt%mcr_scale  = done

     allocate(poollt%mrr_freeze_lay(nsoil))
     allocate(poollt%mrr_hot_lay(nsoil))
     allocate(poollt%mrr_scale_lay(nsoil))
     poollt%mrr_freeze_lay(:) = done
     poollt%mrr_hot_lay(:) = done
     poollt%mrr_scale_lay(:) = done
     poollt%mrr_assim = done
     poollt%mrr_freeze = done
     poollt%mrr_hot = done
     poollt%mrr_lai = done
     poollt%mrr_scale = done

     allocate(poollt%krater_lay(npoolpft,nsoil))
     poollt%krater_lay(:,:) = dzero
     allocate(poollt%loss_mresp_lay(npoolpft,nsoil))
     poollt%loss_mresp_lay(:,:) = dzero

     poollt%resp_auto = dzero
     poollt%resp_grow = dzero
     poollt%resp_leaf = dzero
     poollt%resp_mntn = dzero
     poollt%resp_nveg = dzero
     poollt%resp_root = dzero
     allocate(poollt%resp_nveg_tmp(npoolpft/3)) !5 lives pools for each C, based on npoolpft=15
     poollt%resp_nveg_tmp(:) = dzero
     allocate(poollt%resp_nvegc13_tmp(npoolpft/3))
     poollt%resp_nvegc13_tmp(:) = dzero
     allocate(poollt%resp_nvegc14_tmp(npoolpft/3))
     poollt%resp_nvegc14_tmp(:) = dzero

     poollt%resp_autoc13 = dzero
     poollt%resp_growc13 = dzero
     poollt%resp_leafc13 = dzero
     poollt%resp_mntnc13 = dzero
     poollt%resp_nvegc13 = dzero
     poollt%resp_rootc13 = dzero

     poollt%resp_autoc14 = dzero
     poollt%resp_growc14 = dzero
     poollt%resp_leafc14 = dzero
     poollt%resp_mntnc14 = dzero
     poollt%resp_nvegc14 = dzero
     poollt%resp_rootc14 = dzero

     poollt%resp_leafc12 = dzero
     !poollt%resp_casn = dzero

     poollt%tfl_daylen = dzero
     poollt%tfl_freeze = dzero
     poollt%tfl_dry = dzero
     poollt%tfl_pstage = dzero
     poollt%tfl_total = dzero
     poollt%tfl_totalc13 = dzero
     poollt%tfl_totalc14 = dzero

     allocate(poollt%tf_turnover(npoolpft))
     poollt%tf_turnover(:) = dzero
     allocate(poollt%loss_trans_lay(npoolpft,nsoil))
     poollt%loss_trans_lay(:,:) = dzero

     poollt%nd_fire = dzero
     poollt%rmmd_fire = dzero
     poollt%resp_fire = dzero
     allocate(poollt%loss_fire_lay(npoolpft,nsoil))
     poollt%loss_fire_lay(:,:) = dzero
     poollt%rmmd_firec13 = dzero
     poollt%resp_firec13 = dzero
     poollt%rmmd_firec14 = dzero
     poollt%resp_firec14 = dzero
     !allocate(poollt%loss_firec13_lay(npoolpft,nsoil))
     !poollt%loss_firec13_lay(:,:) = dzero
     
     poollt%nd_grz = dzero
     poollt%resp_grz = dzero
     poollt%resp_grzc13 = dzero
     allocate(poollt%loss_grz(npoolcan))
     poollt%loss_grz = dzero
     allocate(poollt%loss_grzc13(npoolcanc13))
     poollt%loss_grzc13 = dzero
     allocate(poollt%loss_grzc14(npoolcanc14))
     poollt%loss_grzc14 = dzero

     poollt%rmvd_hrvst = dzero
     poollt%rmvd_hrvstc13 = dzero
     poollt%resp_hrvst = dzero
     poollt%resp_hrvstc13 = dzero
     poollt%resp_hrvst = dzero
     poollt%resp_hrvstc14 = dzero
     allocate(poollt%loss_hrvst_lay(npoolpft,nsoil))
     poollt%loss_hrvst_lay(:,:) = dzero

     allocate(poollt%poolpft_dloss(npoolpft,nsoil))
     poollt%poolpft_dloss(:,:) = dzero
     allocate(poollt%poolpft_dgain(npoolpft,nsoil))
     poollt%poolpft_dgain(:,:) = dzero

     allocate(poollt%poolpft(npoolpft))
     poollt%poolpft(:) = dzero
     allocate(poollt%poolpftmin(npoolpft))
     poollt%poolpftmin(:) = dzero
     allocate(poollt%poolpft_lay(npoolpft,nsoil))
     poollt%poolpft_lay(:,:) = dzero
     allocate(poollt%poolpft_flay(npoolpft,nsoil))
     poollt%poolpft_flay(:,:) = dzero

     allocate(poollt%poolpftp(npoolpft))
     poollt%poolpftp(:) = dzero

     !... Additional diagnostics for C13
     allocate(poollt%mrespr(npoolpft))
     poollt%mrespr(:) = dzero
     allocate(poollt%pftavail_lay(npoolpft,nsoil))
     poollt%pftavail_lay(:,:) = dzero
     poollt%tfrac_lp = dzero
     poollt%tfrac_lpc13 = dzero
     poollt%tfrac_lpc14 = dzero
     allocate(poollt%assim_dgain(npoolpft,nsoil))
     poollt%assim_dgain(:,:) = dzero
     allocate(poollt%assim_dloss(npoolpft,nsoil))
     poollt%assim_dloss(:,:) = dzero
     allocate(poollt%autoresp_dgain(npoolpft,nsoil))
     poollt%autoresp_dgain(:,:) = dzero
     allocate(poollt%autoresp_dloss(npoolpft,nsoil))
     poollt%autoresp_dloss(:,:) = dzero
     allocate(poollt%autotran_dloss(npoolpft,nsoil))
     poollt%autotran_dloss(:,:) = dzero
     allocate(poollt%poolpftmin_updated(npoolpft))
     poollt%poolpftmin_updated(:) = dzero
     allocate(poollt%isofactorp(npoolpft))
     poollt%isofactorp(:) = dzero
     allocate(poollt%rcpoolpft_lay(npoolpft,nsoil))
     poollt%rcpoolpft_lay(:,:) = dzero
     allocate(poollt%rcpoolpft(npoolpft))
     poollt%rcpoolpft(:) = dzero
     allocate(poollt%curpoolpft(npoolpft))
     poollt%curpoolpft(:) = dzero

     !poollt%rcpoolfire = dzero
     !poollt%poolemistotC = dzero
     !poollt%poolemisc13 = dzero

end subroutine init_poollt


!--------------------------------------------------------------
!Routine to initialize the radiation variables
subroutine init_radt(radt)

use kinds
use module_sib, only: rad_type

implicit none

!...input variables
type (rad_type), intent(inout) :: radt

     radt%albedo_visb = dzero
     radt%albedo_visd = dzero
     radt%albedo_nirb = dzero
     radt%albedo_nird = dzero

     radt%radfacc(:,:) = dzero
     radt%radfacg(:,:) = dzero

     radt%radc3c = dzero
     radt%radc3g = dzero
     radt%radtc = dzero
     radt%radtg = dzero
     radt%radts = dzero
     radt%effgc = dzero

     radt%tsfc = dzero

end subroutine init_radt


!--------------------------------------------------------
!Routine to initialize the fluorescence (SIF) variables
subroutine init_sift(sift)

use kinds
use module_sib, only: sif_type

implicit none

!...input variables
type(sif_type), intent(inout) :: sift

     sift%sif_je = dzero
     sift%sif_jo = dzero
     sift%sif_jejo = dzero
     sift%sif_x = dzero

     sift%sif_kd = dzero
     sift%sif_kn = dzero
     sift%sif_kp = dzero

     sift%phi_d = dzero
     sift%phi_f = dzero
     sift%phi_n = dzero
     sift%phi_p = dzero

     sift%sif = dzero

end subroutine init_sift


!--------------------------------------------------------
!Routine to initialize the soil/snow column variables
subroutine init_sscolt(nsnow, nsoil, sscolt)

use kinds
use module_sib, only: sscol_type

implicit none

!...input variables
integer(i4), intent(in) :: nsnow, nsoil
type (sscol_type), intent(inout) :: sscolt

     sscolt%nsl = bzero

     allocate(sscolt%rootr(nsoil))
     allocate(sscolt%satfrac_lay(nsoil))
     allocate(sscolt%eff_poros(-nsnow+1:nsoil))
     allocate(sscolt%layer_z(-nsnow:nsoil))
     allocate(sscolt%node_z(-nsnow+1:nsoil))
     allocate(sscolt%shcap(-nsnow+1:nsoil))
     allocate(sscolt%slamda(-nsnow+1:nsoil))
     allocate(sscolt%tksoil(-nsnow+1:nsoil))
     allocate(sscolt%vol_liq(-nsnow+1:nsoil))
     allocate(sscolt%vol_ice(-nsnow+1:nsoil))

     allocate(sscolt%td(-nsnow+1:nsoil))
     allocate(sscolt%www_liq(-nsnow+1:nsoil))
     allocate(sscolt%www_ice(-nsnow+1:nsoil))
     allocate(sscolt%dz(-nsnow+1:nsoil))

     sscolt%rootr(:) = dzero
     sscolt%satfrac_lay(:) = dzero
     sscolt%eff_poros(:) = dzero
     sscolt%layer_z(:) = dzero
     sscolt%node_z(:) = dzero
     sscolt%shcap(:) = dzero
     sscolt%slamda(:) = dzero
     sscolt%tksoil(:) = dzero
     sscolt%vol_liq(:) = dzero
     sscolt%vol_ice(:) = dzero

     sscolt%dz(:) = dzero
     sscolt%td(:) = dzero
     sscolt%www_liq(:) = dzero
     sscolt%www_ice(:) = dzero

end subroutine init_sscolt


!--------------------------------------------------------
!Routine to initialize the vegetation variables
subroutine init_vegt(nsoil, vegt)

use kinds
use module_sib, only: veg_type

implicit none

!...input variables
integer(i4), intent(in) :: nsoil
type (veg_type), intent(inout) :: vegt

     vegt%z0d = dzero
     vegt%zp_dispd = dzero
     vegt%zpd_adj = dzero
     vegt%zztemp = dzero
     vegt%zzwind = dzero
     vegt%cc1 = dzero
     vegt%cc2 = dzero
    
     allocate(vegt%rootf(nsoil))
     vegt%rootf(:) = dzero

     vegt%fpar = dzero
     vegt%green = dzero
     vegt%lai = dzero
     vegt%lait = dzero
     vegt%vcover = dzero

     vegt%gmudmu = dzero
     vegt%park = dzero
     vegt%vmax = dzero

     vegt%clim_lai = dzero

end subroutine init_vegt

