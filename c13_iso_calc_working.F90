!  c13_iso_calc calculates 13C and 12C 
!  for respiration fluxes and C13 pools
!
subroutine c13_iso_calc( &
           poolcont, &
           poollt, pooldt, &
           fract)

use kinds
use module_fractsib
use module_sib, only: &
    pooll_type, poold_type, &
    fract_type
use module_pparams,only: &
    pdb
use module_time, only: &
    dtsib, wt_daily
use module_poolinfo
use module_sibconst, only: &
    npoolpft, nsoil, npoollu, &
    varc13_switch, varco2_switch
!    npoollu, &
!    npoolsfc, npoolsoil, &
!    npoolsfcc13, npoolsoilc13
use module_phosib, only: &
    c4, gah2o, co2cap!, &
use module_param, only: &
   pool_param

implicit none

!Input Variables
type(fract_type), intent(inout) :: fract
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(pool_param), intent(inout) :: poolcont
!real(r8), intent(in) :: &
!    resp_autoc13, resp_growc13, &
!    resp_leafc13, resp_mntnc13, &
!    resp_rootc13, resp_hetc13, &
!    resp_firec13, resp_grzc13, &
!    resp_hrvstc13, resp_nvegc13, &
!    resp_soilc13

integer(i4) :: n,s,sref,tcref

integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13

real(r8) :: tottemp, totc13temp
real(r8) :: tmpmin
real(r8) :: tmpval

lp =  pool_indx_leaf !ntpool index 1
frp = pool_indx_froot !ntpool index 2
crp = pool_indx_croot !ntpool index 3
wp =  pool_indx_stwd !ntpool index 4
pp =  pool_indx_prod !ntpool index 5

lpc13 =  pool_indx_leaf_c13-6 !ntpool index 12, npoolpft index 6
frpc13 = pool_indx_froot_c13-6 !ntpool index 13, npoolpft index 7
crpc13 = pool_indx_croot_c13-6 !ntpool index 14, npoolpft index 8
wpc13 =  pool_indx_stwd_c13-6 !ntpool index 15, npoolpft index 9
ppc13 =  pool_indx_prod_c13-6 !ntpool index 16, npoolpft index 10

cdbp  = pool_indx_cdb-npoolpft/2 !ntpool index 6, npoollu index 1
metlp = pool_indx_metl-npoolpft/2 !ntpool index 7, npoollu index 2
strlp = pool_indx_strl-npoolpft/2 !ntpool index 8, npoollu index 3
slitp = pool_indx_slit-npoolpft/2 !ntpool index 9, npoollu index 4
slowp = pool_indx_slow-npoolpft/2 !ntpool index 10, npoollu index 5
armp  = pool_indx_arm-npoolpft/2 !ntpool index 11, npoollu index 6

cdbpc13  = pool_indx_cdb_c13 - npoolpft !ntpool index 17, npoollu index 7
metlpc13  = pool_indx_metl_c13 - npoolpft !ntpool index 18, npoollu index 8
strlpc13  = pool_indx_strl_c13 - npoolpft !ntpool index 19, npoollu index 9
slitpc13  = pool_indx_slit_c13 - npoolpft !ntpool index 20, npoollu index 10
slowpc13 = pool_indx_slow_c13 - npoolpft !ntpool index 21, npoollu index 11
armpc13  = pool_indx_arm_c13 - npoolpft !ntpool index 22, npoollu index 12


!...Initialize values
!fract%d13cpool_leafc13 = dzero
!fract%d13cpool_frootc13 = dzero
!fract%d13cpool_crootc13 = dzero
!fract%d13cpool_stwdc13 = dzero
!fract%d13cpool_prodc13 = dzero
!fract%d13cpool_cdbc13 = dzero
!fract%d13cpool_metlc13 = dzero
!fract%d13cpool_strlc13 = dzero
!fract%d13cpool_slitc13 = dzero
!fract%d13cpool_slowc13 = dzero
!fract%d13cpool_armc13 = dzero

!-------Calculate isotope signatures of pools and respiration fluxes--------

!...Isotope signatures of pools
if (poollt%poolpft(lpc13) .GT. dzero) then
  fract%d13cpool_leafc13 = dble(( ( (poollt%poolpft(lpc13)/ &
                              (poollt%poolpft(lp)-poollt%poolpft(lpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif

if (poollt%poolpft(frpc13) .GT. dzero) then
  fract%d13cpool_frootc13 = dble(( ( (poollt%poolpft(frpc13)/ &
                              (poollt%poolpft(frp)-poollt%poolpft(frpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif

!print*,' '
!print*,'d13cpool_leafc13: ',fract%d13cpool_leafc13
!print*,' '
!print*,' '
!print*,'d13cpool_frootc13: ',fract%d13cpool_frootc13
!print*,' '

if (poollt%poolpft(crpc13) .GT. dzero) then
  fract%d13cpool_crootc13 = dble(( ( (poollt%poolpft(crpc13)/ &
                              (poollt%poolpft(crp)-poollt%poolpft(crpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif
if (poollt%poolpft(wpc13) .GT. dzero) then
  fract%d13cpool_stwdc13 = dble(( ( (poollt%poolpft(wpc13)/ &
                              (poollt%poolpft(wp)-poollt%poolpft(wpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif
if (poollt%poolpft(ppc13) .GT. dzero) then
  fract%d13cpool_prodc13 = dble(( ( (poollt%poolpft(ppc13)/ &
                              (poollt%poolpft(pp)-poollt%poolpft(ppc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif

if (pooldt%poollu(cdbpc13) .GT. dzero) then
  fract%d13cpool_cdbc13 = dble(( ( (pooldt%poollu(cdbpc13)/ &
                              (pooldt%poollu(cdbp)-pooldt%poollu(cdbpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif
if (pooldt%poollu(metlpc13) .GT. dzero) then
  fract%d13cpool_metlc13 = dble(( ( (pooldt%poollu(metlpc13)/ &
                              (pooldt%poollu(metlp)-pooldt%poollu(metlpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif 
if (pooldt%poollu(strlpc13) .GT. dzero) then
  fract%d13cpool_strlc13 = dble(( ( (pooldt%poollu(strlpc13)/ &
                              (pooldt%poollu(strlp)-pooldt%poollu(strlpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif
if (pooldt%poollu(slitpc13) .GT. dzero) then
  fract%d13cpool_slitc13 = dble(( ( (pooldt%poollu(slitpc13)/ &
                              (pooldt%poollu(slitp)-pooldt%poollu(slitpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif
if (pooldt%poollu(slowpc13) .GT. dzero) then
  fract%d13cpool_slowc13 = dble(( ( (pooldt%poollu(slowpc13)/ &
                              (pooldt%poollu(slowp)-pooldt%poollu(slowpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif
if (pooldt%poollu(armpc13) .GT. dzero) then
  fract%d13cpool_armc13 = dble(( ( (pooldt%poollu(armpc13)/ &
                              (pooldt%poollu(armp)-pooldt%poollu(armpc13)))/pdb) &
                               - 1.0D0)*1000.0D0)
endif

!...Isotope signatures of respiration fluxes

if (poollt%resp_auto .GT. dzero) then
  fract%d13cresp_autoc13 = &
      (((poollt%resp_autoc13/(poollt%resp_auto-poollt%resp_autoc13))/pdb)-1.0D0)*1000.0D0  
endif

if (poollt%resp_grow .GT. dzero) then
  fract%d13cresp_growc13 = &
      (((poollt%resp_growc13/(poollt%resp_grow-poollt%resp_growc13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_leaf .GT. dzero) then
  fract%d13cresp_leafc13 = &
      (((poollt%resp_leafc13/(poollt%resp_leaf-poollt%resp_leafc13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_mntn .GT. dzero) then
  fract%d13cresp_mntnc13 = &
      (((poollt%resp_mntnc13/(poollt%resp_mntn-poollt%resp_mntnc13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_root .GT. dzero) then
  fract%d13cresp_rootc13 = &
      (((poollt%resp_rootc13/(poollt%resp_root-poollt%resp_rootc13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_fire .GT. dzero) then
  fract%d13cresp_firec13 = &
      (((poollt%resp_firec13/(poollt%resp_fire-poollt%resp_firec13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_grz .GT. dzero) then
  fract%d13cresp_grzc13 = &
      (((poollt%resp_grzc13/(poollt%resp_grz-poollt%resp_grzc13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_hrvst .GT. dzero) then
  fract%d13cresp_hrvstc13 = &
      (((poollt%resp_hrvstc13/(poollt%resp_hrvst-poollt%resp_hrvstc13))/pdb)-1.0D0)*1000.0D0
endif

if (poollt%resp_nveg .GT. dzero) then
  fract%d13cresp_nvegc13 = &
      (((poollt%resp_nvegc13/(poollt%resp_nveg-poollt%resp_nvegc13))/pdb)-1.0D0)*1000.0D0
endif


if (pooldt%resp_het .GT. dzero) then
  fract%d13cresp_hetc13 = &
      (((pooldt%resp_hetc13/(pooldt%resp_het-pooldt%resp_hetc13))/pdb)-1.0D0)*1000.0D0
endif

if (pooldt%resp_soil .GT. dzero) then
  fract%d13cresp_soilc13 = &
      (((pooldt%resp_soilc13/(pooldt%resp_soil-pooldt%resp_soilc13))/pdb)-1.0D0)*1000.0D0
endif

tottemp = (poollt%resp_auto + pooldt%resp_het + poollt%resp_nveg &
           + poollt%resp_grz + poollt%resp_hrvst)

totc13temp = (poollt%resp_autoc13 + pooldt%resp_hetc13 + poollt%resp_nvegc13 &
           + poollt%resp_grzc13 + poollt%resp_hrvstc13)

fract%c13resptot=totc13temp
fract%c12resptot=(tottemp-totc13temp)

if (tottemp .GT. dzero) then
  fract%d13cresp_totc13 = &
            (((totc13temp/(tottemp-totc13temp))/pdb)-1.0D0)*1000.0D0
endif

!!.. Update d13cca
!!.. following similar calculation in phosib for CO2

if (varc13_switch .or. varco2_switch) then
  fract%c13ca = dble(fract%c13ca + dble(dtsib / co2cap) * dble(fract%c13resptot - &
                    fract%c13assim + (fract%c13cm*gah2o ))) &
                    / (1.0D0 + dble(dtsib*gah2o / co2cap))
  fract%c12ca = dble(fract%c12ca + dble(dtsib / co2cap) * dble(fract%c12resptot - &
                     fract%c12assim + (fract%c12cm*gah2o ))) &
                     / (1.0D0 +dble(dtsib*gah2o / co2cap))
  fract%d13cca_updated = ((fract%c13ca/fract%c12ca - pdb) / pdb) *1000.D0
endif


!...Calculate rcpool for each pool/layer
!...uses tmpval to check for mass in total pools first

do n=npoolpft/2+1,npoolpft !(6,10) C13 pools
   tcref=n-5 !(1,5) totC live pools
   !tsref=n-5 !ntpool for totC live pools
   sref=n+npoolpft/2+1 !12,16 ntpool

   !if (poollt%poolpft(tcref) .gt. dzero) then
   if (pool_indx_lay(sref) .eq. 1) then
      tmpval=( poollt%poolpft(tcref) &
              + poollt%poolpft_dgain(tcref,1) &
              - poollt%poolpft_dloss(tcref,1) )
     if (tmpval .gt. dzero) then
         poollt%rcpoolpft_lay(n,1) = &
            ( poollt%poolpft(n) &
              + poollt%poolpft_dgain(n,1) &
              - poollt%poolpft_dloss(n,1) )/ &
            ( poollt%poolpft(tcref) &
              + poollt%poolpft_dgain(tcref,1) &
              - poollt%poolpft_dloss(tcref,1) )
         poollt%rcpoolpft(n) = poollt%rcpoolpft_lay(n,1)
     endif !tmpval >0
   else
       do s=1,pool_indx_lay(sref) 
         tmpval=( poollt%poolpft(tcref) &
                + poollt%poolpft_dgain(tcref,s) &
                - poollt%poolpft_dloss(tcref,s) )
         if (tmpval .gt. dzero) then      
           !pool_indx_lay the same sref for live pools in sequence
           poollt%rcpoolpft_lay(n,s) = &
              ( poollt%poolpft(n) &
                + poollt%poolpft_dgain(n,s) &
                - poollt%poolpft_dloss(n,s) )/ &
              ( poollt%poolpft(tcref) &
                + poollt%poolpft_dgain(tcref,s) &
                - poollt%poolpft_dloss(tcref,s) )
         endif !tmpval >0
       enddo !soil layers
       tmpval=( poollt%poolpft(tcref) &
              + sum(poollt%poolpft_dgain(tcref,:)) &
              - sum(poollt%poolpft_dloss(tcref,:)) )
       if (tmpval .gt. dzero) then
         poollt%rcpoolpft(n) = &
            ( poollt%poolpft(n) &
              + sum(poollt%poolpft_dgain(n,:)) &
              - sum(poollt%poolpft_dloss(n,:)) )/ &
            ( poollt%poolpft(tcref) &
              + sum(poollt%poolpft_dgain(tcref,:)) &
              - sum(poollt%poolpft_dloss(tcref,:)) )
       endif !tmpval >0
   endif ! pool_indx_lay
enddo !poolpft loop

do n=npoollu/2+1,npoollu !(7,12) C13 pools
   tcref=n-6 !(1,6) totC dead pools
   !tsref=n !ntpool for totC dead pools
   sref=n+npoolpft !17,22 ntpool
   !if (pooldt%poollu(tcref) .gt. dzero) then
   if (pool_indx_lay(sref) .eq. 1) then
     tmpval=( pooldt%poollu(tcref) &
              + pooldt%poollu_dgain(tcref,1) &
              - pooldt%poollu_dloss(tcref,1) )
     if (tmpval .gt. dzero) then
       pooldt%rcpoollu_lay(n,1) = &
          ( pooldt%poollu(n) &
            + pooldt%poollu_dgain(n,1) &
            - pooldt%poollu_dloss(n,1) )/ &
          ( pooldt%poollu(tcref) &
            + pooldt%poollu_dgain(tcref,1) &
            - pooldt%poollu_dloss(tcref,1) )
       pooldt%rcpoollu(n) = pooldt%rcpoollu_lay(n,1)
     endif !tmpval >0
   else
       do s=1,pool_indx_lay(sref)
       !pool_indx_lay the same sref for dead pools in sequence
         tmpval=( pooldt%poollu(tcref) &
              + pooldt%poollu_dgain(tcref,s) &
              - pooldt%poollu_dloss(tcref,s) )
         if (tmpval .gt. dzero) then
           pooldt%rcpoollu_lay(n,s) = &
            ( pooldt%poollu(n) &
              + pooldt%poollu_dgain(n,s) &
              - pooldt%poollu_dloss(n,s) )/ &
            ( pooldt%poollu(tcref) &
              + pooldt%poollu_dgain(tcref,s) &
              - pooldt%poollu_dloss(tcref,s) )
         endif !tmpval >0
       enddo !soil layers
       tmpval=( pooldt%poollu(tcref) &
              + sum(pooldt%poollu_dgain(tcref,:)) &
              - sum(pooldt%poollu_dloss(tcref,:)) )
       if (tmpval .gt. dzero) then
         pooldt%rcpoollu(n) = &
            ( pooldt%poollu(n) &
              + sum(pooldt%poollu_dgain(n,:)) &
              - sum(pooldt%poollu_dloss(n,:)) )/ &
            ( pooldt%poollu(tcref) &
              + sum(pooldt%poollu_dgain(tcref,:)) &
              - sum(pooldt%poollu_dloss(tcref,:)) )
       endif !tmpval>0
   endif !pool_indx_lay
enddo !poollu loop

!...Update poolpft_min for the isotope pools
tmpmin = dzero
!...Update poolpft_min for isotope pools based on current pool sizes
do n=npoolpft/2+1,npoolpft !(6,10)
   tcref=n-5
   if (poollt%rcpoolpft(n) .gt. dzero) then
   !tmpmin = (fract%rcassim/(fract%rcassim+1.0D0)) &
   !                                * poolcont%poolpft_min(tcref)
     tmpmin = (poollt%rcpoolpft(n)) * poolcont%poolpft_min(tcref)
     poolcont%poolpft_min(n) = tmpmin
   endif
!print*,'tmpmin: ',tmpmin

enddo


!.. Get grid cell pool ratio for (lp+wp+cdb+metl+strl) for fire

tmpval = ( (poollt%poolpft(lp) &
            + poollt%poolpft(wp) &
            + pooldt%poollu(cdbp) &
            + pooldt%poollu(metlp) &
            + pooldt%poollu(strlp)) &
         + (poollt%poolpft_dgain(lp,1) &
            + poollt%poolpft_dgain(wp,1) &
            + pooldt%poollu_dgain(cdbp,1) &
            + pooldt%poollu_dgain(metlp,1) &
            + pooldt%poollu_dgain(strlp,1)) &
         - (poollt%poolpft_dloss(lp,1) &
            + poollt%poolpft_dloss(wp,1) &
            + pooldt%poollu_dloss(cdbp,1) &
            + pooldt%poollu_dloss(metlp,1) &
            + pooldt%poollu_dloss(strlp,1)) )

!if ( (poollt%poolpft(lpc13) .gt. dzero) .or. &
!     (poollt%poolpft(wpc13) .gt. dzero) .or. &
!     (pooldt%poollu(cdbpc13) .gt. dzero) .or. &
!     (pooldt%poollu(metlpc13) .gt. dzero) .or. &
!     (pooldt%poollu(strlpc13) .gt. dzero)) then

if (tmpval .gt. dzero) then
      poollt%rcpoolfire = & 
         ( (poollt%poolpft(lpc13) &
            + poollt%poolpft(wpc13) &
            + pooldt%poollu(cdbpc13) &
            + pooldt%poollu(metlpc13) &
            + pooldt%poollu(strlpc13)) &
         + (poollt%poolpft_dgain(lpc13,1) &
            + poollt%poolpft_dgain(wpc13,1) &
            + pooldt%poollu_dgain(cdbpc13,1) &
            + pooldt%poollu_dgain(metlpc13,1) &
            + pooldt%poollu_dgain(strlpc13,1)) &
         - (poollt%poolpft_dloss(lpc13,1) &
            + poollt%poolpft_dloss(wpc13,1) &
            + pooldt%poollu_dloss(cdbpc13,1) &                
            + pooldt%poollu_dloss(metlpc13,1) &
            + pooldt%poollu_dloss(strlpc13,1)) )/ &
         ( (poollt%poolpft(lp) &
            + poollt%poolpft(wp) &
            + pooldt%poollu(cdbp) &
            + pooldt%poollu(metlp) &
            + pooldt%poollu(strlp)) &
         + (poollt%poolpft_dgain(lp,1) &
            + poollt%poolpft_dgain(wp,1) &
            + pooldt%poollu_dgain(cdbp,1) &
            + pooldt%poollu_dgain(metlp,1) &
            + pooldt%poollu_dgain(strlp,1)) &
         - (poollt%poolpft_dloss(lp,1) &
            + poollt%poolpft_dloss(wp,1) &
            + pooldt%poollu_dloss(cdbp,1) &              
            + pooldt%poollu_dloss(metlp,1) &
            + pooldt%poollu_dloss(strlp,1)) )
endif

end subroutine c13_iso_calc
