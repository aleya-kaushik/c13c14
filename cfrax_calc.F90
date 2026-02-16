!  CFRAX calculates 13C and 12C fluxes and concentrations in the canopy,
!  mixed layer, carbon assimilation (photosynthate), respired soil
!  carbon,
!  assuming that discrimination against 13C during photosynthesis is a 
!  function of discrimination during diffusion and the pCO2i/pCO2c
!  ratio.
!  C4 discrimination against 13C only results from diffusion. 
!
subroutine cfrax_calc( &
           physcont, gprogt, &
           co2t, fract, fluxt, &
           !co2cas, co2s, co2i, co2c, &
           co2m, co2gamma, co2_assim)!, &
!           isoyrtmp, globc13tmp)

use kinds
use module_fractsib
use module_param, only: &
    phys_param
use module_pparams, only: &
    p0_sfc, pdb, &     
    kieclfbl, kiecstom, &
    kieclphas, kiecdis, &
    kiecrbsco, kiecphtrsp
use module_phosib, only: &
    c4, gah2o, co2cap!, &
    !co2s, co2c, co2i, &
    !co2cas, co2m, co2gamma
use module_sib, only: &
    gprog_type, fract_type, &
    co2_type, flux_type
use module_time!, only: &
    !dtsib, wt_daily, year, &
    !month, hour, day
use module_sibconst, only: &
    nisodatayr, varc13m_switch, &
    varc13_switch, varco2_switch
use module_isodata, only: &
    isoyr, globc13

implicit none

!.. Input Variables
type(phys_param), intent(in) :: physcont
type(gprog_type), intent(in) :: gprogt
type(flux_type), intent(in) :: fluxt
type(fract_type), intent(inout) :: fract
type(co2_type), intent(in) :: co2t
!real(r8), intent(in) :: co2cas, co2s, co2i, co2c
real(r8), intent(in) :: co2m, co2gamma, co2_assim

!.. Local variables
integer(i4) :: yrnow, i, loc
real(r8) :: xx
logical :: cbad
!real(r8), intent(in) :: isoyr, globc13
!intrinsic :: findloc
    !Local Variables
    !real(r8), intent(in)  :: pco2s
    !real(r8), intent(in)  :: pco2i
    !real(r8), intent(in)  :: pco2c
    !real(r8), intent(in) :: gamma
    !real(r8) :: pco2m
    !real(r8), intent(in) :: ps, pco2m
    !real(r8) :: xco2cas
    !real(r8) :: xco2s
    !real(r8) :: xco2i
    !real(r8) :: xco2c
    !real(r8) :: xco2m
    !real(r8) :: xco2gamma
real(r8) :: rca
real(r8) :: rcm
real(r8) :: nzero=1.E-14
    
    !real(r8) :: d13cca
    !real(r8) :: d13cm

fract%press_cfrax = dble(gprogt%ps) * 100.0D0
fract%press_cfraxps = dble(gprogt%ps)


!..reset this value
!fract%d13cassimxassim = dzero

if (varc13m_switch) then
   !..Update d13cm from isodata input
   yrnow=year
   !idx=findloc(isoyrtmp,yrnow+0.5)
   do i=1,nisodatayr
     if (floor(isoyr(i)) .eq. yrnow) then
      loc=i 
      exit
     endif
   enddo
  fract%d13cm = dble(globc13(loc))
else
  fract%d13cm = -7.6
endif

!d13cca=-8.0
if (varc13_switch .or. varco2_switch) then
   !..Update d13cca from c13_iso_calc
!  if ((fract%d13cca_updated .gt. -100.) .and. & 
!      (fract%d13cca_updated .lt. 0.)) then
!   fract%d13cca = fract%d13cca_updated
!  else
!   fract%d13cca = fract%d13cm
!  endif
   !print*,'d13cca: ',fract%d13cca
   ! below to turn off recyling for isotopes as well
   fract%d13cca = fract%d13cm
else
   fract%d13cca = -7.6
endif

!fract%d13cca = d13cca
!fract%d13cm  = d13cm
fract%co2m_cfrax=co2m

!print*,'d13cca :',fract%d13cca

!!!!!!!!!!!!!!!!!!!!!!!!!!

    !ps = gprogt%ps * 100.0

    ! Convert from pa to ppm
    !xco2cas = co2t%pco2cas / (gprogt%ps * 100.0)
    !xco2s   = co2t%pco2s / (gprogt%ps * 100.0)
    !xco2i   = co2t%pco2i / (gprogt%ps * 100.0)
    !xco2c   = co2t%pco2c / (gprogt%ps * 100.0)
    !xco2m = co2m
    !xco2m = gprogt%pco2m / (gprogt%ps * 100.0) 
    !xco2cas = co2t%co2cas
    !xco2s = co2t%co2s 
    !xco2i = co2t%co2i 
    !xco2c = co2t%co2c 
    !xco2m = co2t%co2m  
    ! Get the CO2 compensation point (gamma) term   
    !xco2_gamma = co2t%gamma
    !xco2gamma = co2t%co2gamma

    ! d13Cca and d13Cm are converted to concentrations (moles/m3) of 13C 
    ! and 12C by first calculating isotope ratios (13C/12C) of the
    ! canopy (ca) and mixed layer (m). 

rca   = dble((dble(fract%d13cca * pdb) / 1000.0D0) + pdb) !0.0111472024
fract%c13ca = (dble(rca) * dble(co2t%co2cas)) / (1.0D0 + dble(rca))
fract%c12ca = dble(co2t%co2cas) / (1.0D0 + dble(rca))
!xx=co2cas-c13ca-c12ca
!print *,'co2cas into cfrax :',co2cas
!print *,'c13ca in cfrax :',c13ca
!print *,'c12ca in cfrax :',c12ca
!print*,'co2cas-c12ca-c13ca:',xx

!fract%c13ca = dble(c13ca)
!fract%c12ca = dble(c12ca)
rcm   = dble((dble(fract%d13cm * pdb) / 1000.0D0) + pdb)
fract%c13cm = (dble(rcm) * dble(co2m)) / (1.0D0 + dble(rcm))
fract%c12cm = dble(co2m) / (1.0D0 + dble(rcm))
!fract%c13cm = dble(c13cm)
!fract%c12cm = dble(c12cm)

! C13/C12 discrimination for C3 plants.  The isotope effect during
! C3 photosynthesis is a function of a combination of the isotope
! effects associated with molecular transport of CO2 across the leaf
! boundary layer (lfbl), into the stoma (stom), dissolution to in mesophyll
! H2O (dis), and transport in the liquid phase (lphas).  The isotope
! effect during C4 photosynthesis is only a function (for now) of transport
! into the stoma.    
! note: IECpsC3 is the isotope effect for carbon isotopic  
! discrimination during photosynthesis of C3 plants.  Similarly for
! KIECpsC4, but for C4 plants. 

!print*,'co2_assim into cfrax: ',co2_assim

!... ADD a check on pco2c until phosib is figured out  
cbad = co2t%co2c*1.e6 > 2000
!if (cbad) then
! print*,'co2c bad from cfrax: ',co2t%co2c*1.e6
! print*,'year: ',year
! print*,'month: ',month
! print*,'doy: ',doy
! print*,'day: ',day
! print*,'hour: ',hour
! print*,'sec: ',sec_day
!endif

!if (yrnow .eq. 1941) then
!  print*,'co2c :',co2c
!endif

!if ((co2t%assim .GT. dzero) .and. (.not. cbad)) then
if (co2t%assim .GT. nzero) then
    if (c4 .EQ. dzero) then
!            fract%kiecps = ( kieclfbl * co2cas + &
!                     (kiecstom-kieclfbl) * co2s + &
!                     (kiecdis + kieclphas - kiecstom) * co2i + &
!                     (kiecrbsco-kiecdis-kieclphas) * co2c ) &
!                     / co2cas
        ! uncheck below to restore
        if (varc13_switch .or. varco2_switch) then 
          fract%kiecps = dble(( kieclfbl * co2t%co2cas + &
                   (kiecstom-kieclfbl) * co2t%co2s + &
                   (kiecdis + kieclphas - kiecstom) * co2t%co2i + &
                   (kiecrbsco-kiecdis-kieclphas) * co2t%co2c - &
                   kiecphtrsp * co2gamma  ) &
                   / co2t%co2cas)
          fract%kiecps_k1 = dble(kieclfbl * co2t%co2cas/co2t%co2cas)
          fract%kiecps_k2 = dble((kiecstom-kieclfbl) * co2t%co2s/co2t%co2cas)
          fract%kiecps_k3 = dble((kiecdis + kieclphas - kiecstom) * &
                                 co2t%co2i/co2t%co2cas)
          fract%kiecps_k4 = dble((kiecrbsco-kiecdis-kieclphas) * co2t%co2c/co2t%co2cas)
          fract%kiecps_k5 = dble(kiecphtrsp * co2gamma/co2t%co2cas)

          fract%kiecps_nog = dble(( kieclfbl * co2t%co2cas + &
                   (kiecstom-kieclfbl) * co2t%co2s + &
                   (kiecdis + kieclphas - kiecstom) * co2t%co2i + &
                   (kiecrbsco-kiecdis-kieclphas) * co2t%co2c ) &
                   / co2t%co2cas)
         
         ! print*,'kiecps: ',fract%kiecps  
         ! fract%kiecpsxassim = fract%kiecps * co2t%assim
         ! fract%kiecpsk1xassim = fract%kiecps_k1 * co2t%assim
         ! fract%kiecpsk2xassim = fract%kiecps_k2 * co2t%assim
         ! fract%kiecpsk3xassim = fract%kiecps_k3 * co2t%assim
         ! fract%kiecpsk4xassim = fract%kiecps_k4 * co2t%assim
         ! fract%kiecpsk5xassim = fract%kiecps_k5 * co2t%assim

        !else below is purely for testing to make sure responses are flat
        else ! const frac for flatc13 or flatco2 case      
          fract%kiecps = -18.0
          fract%kiecps_nog = -18.0
        !  fract%kiecpsxassim = fract%kiecps * co2t%assim
        endif
    else      !C4 plants given constant KIE = 4.4per mil  
        fract%kiecps = dble(kiecstom)
        fract%kiecps_nog = dble(kiecstom)
        !fract%kiecpsxassim = fract%kiecps * co2t%assim
    endif !plant type selection based on c4
    !!!else
    !
    ! We need values for when Assimn < 0.0, i.e. nighttime.  We revert
    ! to the del13C of respiration both for convenience and because the
    ! nighttime flux to the atmosphere reflects the 13C/12C ratio of
    ! respiration
    !!!  if(c4flag == 0.0) then
    !!!     if(pdiagt%kiecps_7d == 0.) then
    !!!        pdiagt%kiecps_7d = -19.2
    !!!     endif
    !!!     pdiagt%kiecps = -19.2 !pdiagt%kiecps_7d
    !!!  else
    !!!     if(pdiagt%kiecps_7d == 0.) then
    !!!        pdiagt%kiecps_7d =  kiecstom
    !!!     endif
    !!!     pdiagt%kiecps = kiecstom !pdiagt%kiecps_7d
    !!!  endif
!endif
    !fract%kiecps = kiecps
    !fract%kiecps_nog = kiecps_nog
    !fract%kiecps_k1 = kiecps_k1    
    !fract%kiecps_k2 = kiecps_k2
    !fract%kiecps_k3 = kiecps_k3
    !fract%kiecps_k4 = kiecps_k4
    !fract%kiecps_k5 = kiecps_k5
!
! calculates d13C of carbon and fluxes of 13C and 12C assimilated 
! IV First calculate during daytime conditions using the actual
! kiecps
!if(co2t%assim .GT. dzero) then
    !kiecps=-18.1451612903
!    rcassim=0.0109450328
    fract%rcassim   =  rca*((fract%kiecps / &
                       1000.0D0) + 1.0D0)
    fract%rcassim_nog   =  rca*((fract%kiecps_nog / &
                       1000.0D0) + 1.0D0)
    !!!    pdiagt%rcrespcan   =   rca*((pdiagt%kiecps_7d / &
    !!!        1000.0) + 1.)
    fract%rcassimfac = fract%rcassim/(fract%rcassim+1.0D0)

    !put in rcassim = hard-coded # corresponding to d13cassim=-26
    fract%d13cassim =  ((fract%rcassim - pdb) / pdb ) * &
                        1000.0D0
    !fract%d13cassimxassim =  fract%d13cassim * co2t%assim

    fract%d13cassim_nog =  ( ((rca*((fract%kiecps_nog / 1000.0D0) &
                           + 1.D0))-pdb)/ pdb) * 1000.0D0
    fract%c13assim_nog  =  (fract%rcassim_nog * co2t%assim)/ &
                        (fract%rcassim_nog+1.0D0)
    fract%c13assim  =  (fract%rcassim * co2t%assim)/ &
                        (fract%rcassim+1.0D0)
    fract%c12assim  =  co2t%assim / (fract%rcassim +1.0D0)

    if (fract%d13cca .lt. -1000. .or. fract%d13cassim .lt. -400.) then ! &
       ! fract%d13cca .gt. 10. .or. fract%d13cassim .lt. -400.) then
     print*,'FROM CFRAX : '
     print*,'d13cassim: ',fract%d13cassim
     print*,'assim: ',co2t%assim*1.0D6
     print*,'kiecps: ',fract%kiecps
     print*,'d13cca cfrax: ',fract%d13cca
     print*,'co2c :',co2t%co2c*1.0D6
     print*,'co2cas :',co2t%co2cas*1.0D6
     print*,'co2s :',co2t%co2s*1.0D6
     print*,'co2i :',co2t%co2i*1.0D6
    endif
 


    !!!    pdiagt%c13assimn  = ( pdiagt%rcassimn * &
    !!!        pdiagt%assimn)/(1. + pdiagt%rcassimn)
    !!!    pdiagt%c12assimn  =  pdiagt%assimn / (1.+ &
    !!!        pdiagt%rcassimn)

        ! IV Summation of fluxes to calculate 7 day mean kiecps in
        ! SiBDRV.F90
    !!!    pdiagt%c13assim_s = pdiagt%c13assim_s +pdiagt%c13assim
    !!!    pdiagt%c12assim_s = pdiagt%c12assim_s + pdiagt%c12assim

    !!!    pdiagt%c13assimn_s = pdiagt%c13assimn_s + pdiagt%c13assimn
    !!!    pdiagt%c12assimn_s = pdiagt%c12assimn_s + pdiagt%c12assimn
    !!!else
     !fract%c13assim=c13assim
     !fract%c12assim=c12assim
else !co2_assim = dzero
   fract%c13assim_nog=dzero
   fract%c13assim=dzero
   fract%c12assim=dzero
   
endif

!fract%rcassim = rcassim
!fract%d13cassim = d13cassim
!fract%d13cassim_nog = d13cassim_nog

fract%c13assimd = dble(fract%c13assim*(1.0D0 - wt_daily) &
                   + fract%c13assim*wt_daily)

    !
    ! Canopy concentrations at time n+1 is calculated using an implicit
    ! scheme.

    ! Note: physcon%z2 in SiB4 is canopy top and PFT-dependent
    ! needs to be read in like physcon(i)%z2

!!! THIS WAS MOVED TO c13_iso_calc
!!!! uncheck below to activate the c13ca/c12ca calculation
!!.. following similar calculation in phosib for CO2
!fract%c13ca = dble(fract%c13ca + dble(dtsib / co2cap) * dble(fract%c13resptot -
!&
!              fract%c13assim + (fract%c13cm*gah2o ))) &
!              / (1.0D0 + dble(dtsib*gah2o / co2cap))
!fract%c12ca = dble(fract%c12ca + dble(dtsib / co2cap) * dble(fract%c12resptot -
!&
!              fract%c12assim + (fract%c12cm*gah2o ))) &
!              / (1.0D0 +dble(dtsib*gah2o / co2cap))
!fract%d13cca = ((fract%c13ca/fract%c12ca - pdb) / pdb) *1000.D0




!!... original from old cfrax code
!fract%c13ca = (fract%c13ca + (dtsib / physcont%z2) * (fract%c13resptot - &
!              fract%c13assim + (fract%c13cm / fluxt%ra ))) &
!              / (1. + (dtsib / fluxt%ra ) / physcont%z2)
!fract%c12ca = (fract%c12ca + (dtsib / physcont%z2) * (fract%c12resptot - &
!              fract%c12assim + (fract%c12cm / fluxt%ra ))) &
!              / (1. + (dtsib / fluxt%ra ) / physcont%z2)

    !
    ! del13C of canopy is recalculated using concentrations of 13Cca and 
    ! 12Cca.  The fluxes (moles/m2/sec) of 13C and 12C out of the canopy 
    ! (the turbulent flux), and the del13C value (per mil vs PDB) of
    ! this 
    ! flux are calculated. 

    !
    ! Use the following if you want the net flux from the canopy to 
    ! be based on differences in 12C and 13C net fluxes from respiration
    ! and photosynthesis
    !           sib%diag%flux13c = sib%diag%resp_grnd * rcresp / (1. +
    !           rcresp) -
    !     &  (sib%diag%assimn(6) * sib%diag%rcassimn / (1. +
    !     sib%diag%rcassimn))
    
    !           sib%diag%flux12c = sib%diag%resp_grnd / (1. + rcresp) -
    !     &  (sib%diag%assimn(6) / (1. + sib%diag%rcassimn))

    ! Use the following if you want the net flux from the canopy to 
    ! be based on differences  in concentration gradients between the 
    ! canopy and overlying atmosphere.
    !!!pdiagt%flux13c = ( pprogt%c13ca - pprogt%c13cm) / pdiagt%ra
    !!!pdiagt%flux12c = ( pprogt%c12ca - pprogt%c12cm) / pdiagt%ra
    !!!pdiagt%flux_turb = pdiagt%flux13c + pdiagt%flux12c


!print*,' '
!print*,'d13cassim: ',fract%d13cassim
!print*,' '

end subroutine cfrax_calc


