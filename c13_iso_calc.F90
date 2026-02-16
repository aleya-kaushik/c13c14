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
    dtsib, wt_daily, year
use module_poolinfo
use module_sibconst, only: &
    npoolpft, nsoil, npoollu, &
    varc13_switch, varco2_switch
!    npoollu, &
!    npoolsfc, npoolsoil, &
!    npoolsfcc13, npoolsoilc13
use module_phosib, only: &
    c4, gah2o, co2cap, &
    resp_cas!, &
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
real(r8) :: tmpval,tmpvdr
real(r8) :: resp_casc13, resp_casc12
real(r8) :: nzero=1.E-14

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



!! print statements to check poolpft_loss and poolpft_lay
!print*,' '
!print*,'code: c13_iso_calc'
!print*,'poolpft_dloss(7,1/2/3):',&
!    poollt%poolpft_dloss(7,1),poollt%poolpft_dloss(7,2),poollt%poolpft_dloss(7,3)
!print*,'poolpft_lay(7,1/2/3) :',&
!    poollt%poolpft_lay(7,1),poollt%poolpft_lay(7,2),poollt%poolpft_lay(7,3)
!print*,' '
 

!...Initialize values
!fract%d13cpool_leafc13 = nzero
!fract%d13cpool_frootc13 = nzero
!fract%d13cpool_crootc13 = nzero
!fract%d13cpool_stwdc13 = nzero
!fract%d13cpool_prodc13 = nzero
!fract%d13cpool_cdbc13 = nzero
!fract%d13cpool_metlc13 = nzero
!fract%d13cpool_strlc13 = nzero
!fract%d13cpool_slitc13 = nzero
!fract%d13cpool_slowc13 = nzero
!fract%d13cpool_armc13 = nzero

!-------Calculate isotope signatures of pools and respiration fluxes--------

!...Isotope signatures of pools
!if (poollt%poolpft(lpc13) .GT. nzero) then
!  fract%d13cpool_leafc13 = dble(( ( (poollt%poolpft(lpc13)/ &
!                              (poollt%poolpft(lp)-poollt%poolpft(lpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!
!if (poollt%poolpft(frpc13) .GT. nzero) then
!  fract%d13cpool_frootc13 = dble(( ( (poollt%poolpft(frpc13)/ &
!                              (poollt%poolpft(frp)-poollt%poolpft(frpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!
!!print*,' '
!!print*,'d13cpool_leafc13: ',fract%d13cpool_leafc13
!!print*,' '
!!print*,' '
!!print*,'d13cpool_frootc13: ',fract%d13cpool_frootc13
!!print*,' '
!
!if (poollt%poolpft(crpc13) .GT. nzero) then
!  fract%d13cpool_crootc13 = dble(( ( (poollt%poolpft(crpc13)/ &
!                              (poollt%poolpft(crp)-poollt%poolpft(crpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!if (poollt%poolpft(wpc13) .GT. nzero) then
!  fract%d13cpool_stwdc13 = dble(( ( (poollt%poolpft(wpc13)/ &
!                              (poollt%poolpft(wp)-poollt%poolpft(wpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!if (poollt%poolpft(ppc13) .GT. nzero) then
!  fract%d13cpool_prodc13 = dble(( ( (poollt%poolpft(ppc13)/ &
!                              (poollt%poolpft(pp)-poollt%poolpft(ppc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!
!if (pooldt%poollu(cdbpc13) .GT. nzero) then
!  fract%d13cpool_cdbc13 = dble(( ( (pooldt%poollu(cdbpc13)/ &
!                              (pooldt%poollu(cdbp)-pooldt%poollu(cdbpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!if (pooldt%poollu(metlpc13) .GT. nzero) then
!  fract%d13cpool_metlc13 = dble(( ( (pooldt%poollu(metlpc13)/ &
!                              (pooldt%poollu(metlp)-pooldt%poollu(metlpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif 
!if (pooldt%poollu(strlpc13) .GT. nzero) then
!  fract%d13cpool_strlc13 = dble(( ( (pooldt%poollu(strlpc13)/ &
!                              (pooldt%poollu(strlp)-pooldt%poollu(strlpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!if (pooldt%poollu(slitpc13) .GT. nzero) then
!  fract%d13cpool_slitc13 = dble(( ( (pooldt%poollu(slitpc13)/ &
!                              (pooldt%poollu(slitp)-pooldt%poollu(slitpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!if (pooldt%poollu(slowpc13) .GT. nzero) then
!  fract%d13cpool_slowc13 = dble(( ( (pooldt%poollu(slowpc13)/ &
!                              (pooldt%poollu(slowp)-pooldt%poollu(slowpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif
!if (pooldt%poollu(armpc13) .GT. nzero) then
!  fract%d13cpool_armc13 = dble(( ( (pooldt%poollu(armpc13)/ &
!                              (pooldt%poollu(armp)-pooldt%poollu(armpc13)))/pdb) &
!                               - 1.0D0)*1000.0D0)
!endif

!...Isotope signatures of respiration fluxes
!...also calculate weights for diagnostic output at finest timestep

if (poollt%resp_auto .GT. nzero) then
  fract%d13cresp_autoc13 = &
      (((poollt%resp_autoc13/(poollt%resp_auto-poollt%resp_autoc13))/pdb)-1.0D0)*1000.0D0  
!  fract%d13crautoxrauto = fract%d13cresp_autoc13 * poollt%resp_auto
endif

if (poollt%resp_grow .GT. nzero) then
  fract%d13cresp_growc13 = &
      (((poollt%resp_growc13/(poollt%resp_grow-poollt%resp_growc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crgrowxrgrow = fract%d13cresp_growc13 * poollt%resp_grow
endif

if (poollt%resp_leaf .GT. nzero) then
  fract%d13cresp_leafc13 = &
      (((poollt%resp_leafc13/(poollt%resp_leaf-poollt%resp_leafc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crleafxrleaf = fract%d13cresp_leafc13 * poollt%resp_leaf
endif

if (poollt%resp_mntn .GT. nzero) then
  fract%d13cresp_mntnc13 = &
      (((poollt%resp_mntnc13/(poollt%resp_mntn-poollt%resp_mntnc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crmntnxrmntn = fract%d13cresp_mntnc13 * poollt%resp_mntn
endif

if (poollt%resp_root .GT. nzero) then
  fract%d13cresp_rootc13 = &
      (((poollt%resp_rootc13/(poollt%resp_root-poollt%resp_rootc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crrootxrroot = fract%d13cresp_rootc13 * poollt%resp_root
endif

if (poollt%resp_fire .GT. nzero) then
  fract%d13cresp_firec13 = &
      (((poollt%resp_firec13/(poollt%resp_fire-poollt%resp_firec13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crfirexrfire = fract%d13cresp_firec13 * poollt%resp_fire
endif

if (poollt%resp_grz .GT. nzero) then
  fract%d13cresp_grzc13 = &
      (((poollt%resp_grzc13/(poollt%resp_grz-poollt%resp_grzc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crgrzxrgrz = fract%d13cresp_grzc13 * poollt%resp_grz
endif

if (poollt%resp_hrvst .GT. nzero) then
  fract%d13cresp_hrvstc13 = &
      (((poollt%resp_hrvstc13/(poollt%resp_hrvst-poollt%resp_hrvstc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crhrvstxrhrvst = fract%d13cresp_hrvstc13 * poollt%resp_hrvst
endif

if (poollt%resp_nveg .GT. nzero) then
  fract%d13cresp_nvegc13 = &
      (((poollt%resp_nvegc13/(poollt%resp_nveg-poollt%resp_nvegc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crnvegxrnveg = fract%d13cresp_nvegc13 * poollt%resp_nveg
endif


if (pooldt%resp_het .GT. nzero) then
  fract%d13cresp_hetc13 = &
      (((pooldt%resp_hetc13/(pooldt%resp_het-pooldt%resp_hetc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crhetxrhet = fract%d13cresp_hetc13 * pooldt%resp_het
endif

if (pooldt%resp_soil .GT. nzero) then
  fract%d13cresp_soilc13 = &
      (((pooldt%resp_soilc13/(pooldt%resp_soil-pooldt%resp_soilc13))/pdb)-1.0D0)*1000.0D0
!  fract%d13crsoilxrsoil = fract%d13cresp_soilc13 * pooldt%resp_soil
endif

tottemp = (poollt%resp_auto + pooldt%resp_het + poollt%resp_nveg &
           + poollt%resp_grz + poollt%resp_hrvst)

totc13temp = (poollt%resp_autoc13 + pooldt%resp_hetc13 + poollt%resp_nvegc13 &
           + poollt%resp_grzc13 + poollt%resp_hrvstc13)

fract%c13resptot=totc13temp
fract%c12resptot=(tottemp-totc13temp)

if (tottemp .GT. nzero) then
  fract%d13cresp_totc13 = &
            (((totc13temp/(tottemp-totc13temp))/pdb)-1.0D0)*1000.0D0
!  fract%d13crtotxrtot = fract%d13cresp_totc13 * tottemp
endif

!!.. Update d13cca
!!.. following similar calculation in phosib for CO2

!if (varc13_switch .or. varco2_switch) then

resp_casc13 = poollt%resp_autoc13 + pooldt%resp_hetc13
resp_casc12 = resp_cas - resp_casc13

poollt%resp_leafc12 = poollt%resp_leaf-poollt%resp_leafc13

!if ((resp_casc13 .gt. nzero) .and. (resp_casc12 .gt. nzero)) then
!!  fract%c13cca = dble(fract%c13ca + dble(dtsib / co2cap) * dble(resp_casc13 - &
!!                    (fract%c13assim - poollt%resp_leafc13) + (fract%c13cm*gah2o ))) &
!!                    / (1.0D0 + dble(dtsib*gah2o / co2cap))
!!  fract%c12cca = dble(fract%c12ca + dble(dtsib / co2cap) * dble(resp_casc12 - &
!!                     (fract%c12assim - poollt%resp_leafc12) + (fract%c12cm*gah2o ))) &
!!                     / (1.0D0 +dble(dtsib*gah2o / co2cap))
!
!  fract%c13cca = dble(fract%c13ca + dble(dtsib / co2cap) * dble(resp_casc13 - &
!                    (fract%c13assim) + (fract%c13cm*gah2o))) &
!                    / (1.0D0 + dble(dtsib*gah2o / co2cap))
!  fract%c12cca = dble(fract%c12ca + dble(dtsib / co2cap) * dble(resp_casc12 - &
!                     (fract%c12assim) + (fract%c12cm*gah2o))) &
!                     / (1.0D0 +dble(dtsib*gah2o / co2cap))
!
!  fract%d13cca_updated = ((fract%c13cca/fract%c12cca - pdb) / pdb) *1000.D0
! if (fract%d13cca .lt. -1000. .or. fract%d13cassim .lt. -400.) then ! &
!       ! fract%d13cca .gt. 10. .or. fract%d13cassim .lt. -400.) then
! !if (fract%d13cca .lt. -1000. .or.  &
! !    fract%d13cca .gt. 10. .or. fract%d13cassim .lt. -400.) then
!  print*,'FROM C13_ISO_CALC: '
!  print*,'d13cca_updated: ',fract%d13cca_updated
!  print*,'d13cassim: ',fract%d13cassim
!  print*,'kiecps: ' ,fract%kiecps
!  print*,'c13assim: ',fract%c13assim
!  print*,'c12assim: ',fract%c12assim
!  print*,'c13ca: ',fract%c13cca
!  print*,'c12ca: ',fract%c12cca
!  print*,'resp_casc13:',resp_casc13
!  print*,'resp_casc13:',resp_casc12
!  print*,'c13cm: ',fract%c13cm
! endif
!endif
!...Calculate rcpool for each pool/layer
!...uses tmpval to check for mass in total pools first

do n=npoolpft/2+1,npoolpft !(6,10) C13 pools
   tcref=n-5 !(1,5) totC live pools
   !tsref=n-5 !ntpool for totC live pools
   sref=n+npoolpft/2+1 !12,16 ntpool

   !if (poollt%poolpft(tcref) .gt. nzero) then
   if (pool_indx_lay(sref) .eq. 1) then
      tmpvdr=( poollt%poolpft_lay(tcref,1) &
              + poollt%poolpft_dgain(tcref,1) &
              - poollt%poolpft_dloss(tcref,1) )
      if (tmpvdr .ne. 0.) then
        tmpval=( poollt%poolpft_lay(n,1) &
              + poollt%poolpft_dgain(n,1) &
              - poollt%poolpft_dloss(n,1) )/ &
            ( poollt%poolpft_lay(tcref,1) &
              + poollt%poolpft_dgain(tcref,1) &
              - poollt%poolpft_dloss(tcref,1) )
         if ((tmpval .lt. 0.011) .and. (tmpval .gt. 0.)) then
            poollt%rcpoolpft_lay(n,1) = tmpval

            poollt%rcpoolpft(n) = poollt%rcpoolpft_lay(n,1)
            poollt%curpoolpft(tcref) = &
                ( poollt%poolpft_lay(tcref,1) &
                  + poollt%poolpft_dgain(tcref,1) &
                  - poollt%poolpft_dloss(tcref,1) )
            poollt%curpoolpft(n) = &
                ( poollt%poolpft_lay(n,1) &
                  + poollt%poolpft_dgain(n,1) &
                  - poollt%poolpft_dloss(n,1) )
     !else
     !   poollt%rcpoolpft(n) = (fract%rcassim/(fract%rcassim+1.0D0))
         endif !tmpval <0.011 and >0
       endif ! tmpvdr ne 0
   else
       do s=1,pool_indx_lay(sref) 
         tmpvdr=( poollt%poolpft_lay(tcref,s) &
                + poollt%poolpft_dgain(tcref,s) &
                - poollt%poolpft_dloss(tcref,s) )
         if (tmpvdr .ne. 0.) then
         tmpval=( poollt%poolpft_lay(n,s) &
                + poollt%poolpft_dgain(n,s) &
                - poollt%poolpft_dloss(n,s) )/ &
              ( poollt%poolpft_lay(tcref,s) &
                + poollt%poolpft_dgain(tcref,s) &
                - poollt%poolpft_dloss(tcref,s) )
           if ((tmpval .lt. 0.011) .and. (tmpval .gt. 0.)) then      
           !pool_indx_lay the same sref for live pools in sequence
             poollt%rcpoolpft_lay(n,s) = tmpval
 
              if ( (poollt%rcpoolpft_lay(n,s) .gt. 1.0) .or. &
                    (poollt%rcpoolpft_lay(n,s) .lt. 0.0) ) then
                 
                print*,' '
                print*,'code: c13_iso_calc'
                print*,'poollt%rcpoolpft_lay(n,s):',poollt%rcpoolpft_lay(n,s)
                if (s .gt. 1) then
                  print*,'poollt%rcpoolpft_lay(n,s-1):',poollt%rcpoolpft_lay(n,s-1)
                endif
                print*,'poollt%rcpoolpft_lay(n,s+1):',poollt%rcpoolpft_lay(n,s+1)
                print*,'poollt%poolpft_lay(n,s):',poollt%poolpft_lay(n,s)
                if (s .gt. 1) then
                  print*,'poollt%poolpft_lay(n,s-1):',poollt%poolpft_lay(n,s-1)
                endif
                print*,'poollt%poolpft_lay(n,s+1):',poollt%poolpft_lay(n,s+1)
                print*,'poollt%poolpft_dgain(n,s):',poollt%poolpft_dgain(n,s)
                if (s .gt. 1) then
                  print*,'poollt%poolpft_dgain(n,s-1):',poollt%poolpft_dgain(n,s-1)
                endif
                print*,'poollt%poolpft_dgain(n,s+1):',poollt%poolpft_dgain(n,s+1)
                print*,'poollt%poolpft_dloss(n,s):',poollt%poolpft_dloss(n,s)
                if (s .gt. 1) then
                  print*,'poollt%poolpft_dloss(n,s-1):',poollt%poolpft_dloss(n,s-1)
                endif
                print*,'poollt%poolpft_dloss(n,s+1):',poollt%poolpft_dloss(n,s+1)
                print*,'poollt%poolpft_lay(tcref,s):',poollt%poolpft_lay(tcref,s)
                if (s .gt. 1) then
                  print*,'poollt%poolpft_lay(tcref,s-1):',poollt%poolpft_lay(tcref,s-1)
                endif
                print*,'poollt%poolpft_lay(tcref,s+1):',poollt%poolpft_lay(tcref,s+1)
                print*,'poollt%poolpft_dgain(tcref,s):',poollt%poolpft_dgain(tcref,s)
                if (s .gt. 1) then
                  print*,'poollt%poolpft_dgain(tcref,s-1):',poollt%poolpft_dgain(tcref,s-1)
                endif
                print*,'poollt%poolpft_dgain(tcref,s+1):',poollt%poolpft_dgain(tcref,s+1)
                print*,'poollt%poolpft_dloss(tcref,s):',poollt%poolpft_dloss(tcref,s)
                if (s .gt. 1) then
                  print*,'poollt%poolpft_dloss(tcref,s-1):',poollt%poolpft_dloss(tcref,s-1)
                endif
                print*,'poollt%poolpft_dloss(tcref,s+1):',poollt%poolpft_dloss(tcref,s+1)   
                print*,' '
              endif

        !else 
         !  poollt%rcpoolpft_lay(n,s) = (fract%rcassim/(fract%rcassim+1.0D0))
           endif !tmpval <0.011 and >0
         endif ! tmpvdr ne 0
       enddo !soil layers

       tmpvdr=( sum(poollt%poolpft_lay(tcref,:)) &
              + sum(poollt%poolpft_dgain(tcref,:)) &
              - sum(poollt%poolpft_dloss(tcref,:)) )
       if (tmpvdr .ne. 0.) then
         tmpval=( sum(poollt%poolpft_lay(n,:)) &
                + sum(poollt%poolpft_dgain(n,:)) &
                - sum(poollt%poolpft_dloss(n,:)) )/ &
              ( sum(poollt%poolpft_lay(tcref,:)) &
                + sum(poollt%poolpft_dgain(tcref,:)) &
                - sum(poollt%poolpft_dloss(tcref,:)) )
         if ((tmpval .lt. 0.011) .and. (tmpval .gt. 0.)) then
             poollt%rcpoolpft(n) = tmpval
             poollt%curpoolpft(tcref) = &
                 ( sum(poollt%poolpft_lay(tcref,:)) &
                  + sum(poollt%poolpft_dgain(tcref,:)) &
                  - sum(poollt%poolpft_dloss(tcref,:)) )
             poollt%curpoolpft(n) = &
                ( sum(poollt%poolpft_lay(n,:)) &
                  + sum(poollt%poolpft_dgain(n,:)) &
                  - sum(poollt%poolpft_dloss(n,:)) )
           !else
           !  poollt%rcpoolpft(n) = (fract%rcassim/(fract%rcassim+1.0D0))
         endif !tmpval <0.011 and >0
       endif ! tmpvdr ne 0
   endif ! pool_indx_lay
enddo !poolpft loop

!... Calculate live pool isotope ratios
if ((poollt%rcpoolpft(lpc13) .gt. nzero) .and. & 
    (poollt%rcpoolpft(lpc13) .lt. done)) then
  fract%d13cpool_leafc13 = &
   (((1.0D0/((1.0D0/poollt%rcpoolpft(lpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((poollt%rcpoolpft(frpc13) .gt. nzero) .and. &
    (poollt%rcpoolpft(frpc13) .lt. done)) then
  !if (year .eq. 2013) then
  ! print*,'poollt%rcpoolpft(frpc13): ',poollt%rcpoolpft(frpc13)
  !endif
  fract%d13cpool_frootc13 = &
   (((1.0D0/((1.0D0/poollt%rcpoolpft(frpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((poollt%rcpoolpft(crpc13) .gt. nzero) .and. & 
    (poollt%rcpoolpft(crpc13) .lt. done)) then
  fract%d13cpool_crootc13 = &
   (((1.0D0/((1.0D0/poollt%rcpoolpft(crpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((poollt%rcpoolpft(wpc13) .gt. nzero) .and. & 
    (poollt%rcpoolpft(wpc13) .lt. done)) then
  fract%d13cpool_stwdc13 = &
   (((1.0D0/((1.0D0/poollt%rcpoolpft(wpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((poollt%rcpoolpft(ppc13) .gt. nzero) .and. & 
    (poollt%rcpoolpft(ppc13) .lt. done)) then
  fract%d13cpool_prodc13 = &
   (((1.0D0/((1.0D0/poollt%rcpoolpft(ppc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
!... Calculate weights for diagnostic output at finest timestep
!fract%d13cpleafxpleaf = fract%d13cpool_leafc13 * poollt%curpoolpft(lp)
!fract%d13cpfrootxpfroot = fract%d13cpool_frootc13 * poollt%curpoolpft(frp)
!fract%d13cpcrootxpcroot = fract%d13cpool_crootc13 * poollt%curpoolpft(crp)
!fract%d13cpstwdxpstwd = fract%d13cpool_stwdc13 * poollt%curpoolpft(wp)
!fract%d13cpprodxpprod = fract%d13cpool_prodc13 * poollt%curpoolpft(pp)

do n=npoollu/2+1,npoollu !(7,12) C13 pools
   tcref=n-6 !(1,6) totC dead pools
   !tsref=n !ntpool for totC dead pools
   sref=n+npoolpft !17,22 ntpool
   !if (pooldt%poollu(tcref) .gt. nzero) then
   if (pool_indx_lay(sref) .eq. 1) then
     tmpvdr = ( pooldt%poollu_lay(tcref,1) &
              + pooldt%poollu_dgain(tcref,1) &
              - pooldt%poollu_dloss(tcref,1) )
     if (tmpvdr .ne. 0.) then
     tmpval= ( pooldt%poollu_lay(n,1) &
            + pooldt%poollu_dgain(n,1) &
            - pooldt%poollu_dloss(n,1) )/ &
            ( pooldt%poollu_lay(tcref,1) &
              + pooldt%poollu_dgain(tcref,1) &
              - pooldt%poollu_dloss(tcref,1) )
         if ((tmpval .lt. 0.011) .and. (tmpval .gt. 0.)) then
           pooldt%rcpoollu_lay(n,1) = tmpval
           pooldt%rcpoollu(n) = pooldt%rcpoollu_lay(n,1)
           pooldt%curpoollu(tcref) = ( pooldt%poollu_lay(tcref,1) &
                + pooldt%poollu_dgain(tcref,1) &
                - pooldt%poollu_dloss(tcref,1) )
           pooldt%curpoollu(n) = &
              ( pooldt%poollu_lay(n,1) &
                + pooldt%poollu_dgain(n,1) &
                - pooldt%poollu_dloss(n,1) )
         !else
         !  pooldt%rcpoollu(n) = (fract%rcassim/(fract%rcassim+1.0D0))
         endif !tmpval <0.011 and >0
     endif !tmpvdr ne 0
   else
       do s=1,pool_indx_lay(sref)
       !pool_indx_lay the same sref for dead pools in sequence
         tmpvdr = ( pooldt%poollu_lay(tcref,s) &
                  + pooldt%poollu_dgain(tcref,s) &
                  - pooldt%poollu_dloss(tcref,s) )
         if (tmpvdr .ne. 0.) then
         tmpval= ( pooldt%poollu_lay(n,s) &
            + pooldt%poollu_dgain(n,s) &
            - pooldt%poollu_dloss(n,s) )/ &
            ( pooldt%poollu_lay(tcref,s) &
              + pooldt%poollu_dgain(tcref,s) &
              - pooldt%poollu_dloss(tcref,s) )
           if ((tmpval .lt. 0.011) .and. (tmpval .gt. 0.)) then
             pooldt%rcpoollu_lay(n,s) = tmpval
           endif !tmpval < 0.011 and >0
         !else
         !  pooldt%rcpoollu_lay(n,s) = (fract%rcassim/(fract%rcassim+1.0D0))
         endif !tmpvdr ne 0
       enddo !soil layers
       tmpvdr=( sum(pooldt%poollu_lay(tcref,:)) &
              + sum(pooldt%poollu_dgain(tcref,:)) &
              - sum(pooldt%poollu_dloss(tcref,:)) )
       if (tmpvdr .ne. 0.) then
       tmpval = &
            ( sum(pooldt%poollu_lay(n,:)) &
              + sum(pooldt%poollu_dgain(n,:)) &
              - sum(pooldt%poollu_dloss(n,:)) )/ &
            ( sum(pooldt%poollu_lay(tcref,:)) &
              + sum(pooldt%poollu_dgain(tcref,:)) &
              - sum(pooldt%poollu_dloss(tcref,:)) )
         if ((tmpval .lt. 0.011) .and. (tmpval .gt. 0.)) then
           pooldt%rcpoollu(n) = tmpval
           pooldt%curpoollu(tcref) = ( sum(pooldt%poollu_lay(tcref,:)) &
                + sum(pooldt%poollu_dgain(tcref,:)) &
                - sum(pooldt%poollu_dloss(tcref,:)) )
           pooldt%curpoollu(n) = &
              ( sum(pooldt%poollu_lay(n,:)) &
                + sum(pooldt%poollu_dgain(n,:)) &
                - sum(pooldt%poollu_dloss(n,:)) )           
       !else
         endif !tmpval <0.011 and >0
       !   pooldt%rcpoollu(n) = (fract%rcassim/(fract%rcassim+1.0D0))  
       endif !tmpvdr ne 0
   endif !pool_indx_lay
enddo !poollu loop

!... Calculate dead pool isotope ratios
if ((pooldt%rcpoollu(cdbpc13) .gt. nzero) .and. &
    (pooldt%rcpoollu(cdbpc13) .lt. done)) then
  fract%d13cpool_cdbc13 = &
   (((1.0D0/((1.0D0/pooldt%rcpoollu(cdbpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((pooldt%rcpoollu(metlpc13) .gt. nzero) .and. &
    (pooldt%rcpoollu(metlpc13) .lt. done)) then
  fract%d13cpool_metlc13 = &
   (((1.0D0/((1.0D0/pooldt%rcpoollu(metlpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((pooldt%rcpoollu(strlpc13) .gt. nzero) .and. &
    (pooldt%rcpoollu(strlpc13) .lt. done)) then
  fract%d13cpool_strlc13 = &
   (((1.0D0/((1.0D0/pooldt%rcpoollu(strlpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((pooldt%rcpoollu(slitpc13) .gt. nzero) .and. &
    (pooldt%rcpoollu(slitpc13) .lt. done)) then
  fract%d13cpool_slitc13 = &
   (((1.0D0/((1.0D0/pooldt%rcpoollu(slitpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((pooldt%rcpoollu(slowpc13) .gt. nzero) .and. &
    (pooldt%rcpoollu(slowpc13) .lt. done)) then
  fract%d13cpool_slowc13 = &
   (((1.0D0/((1.0D0/pooldt%rcpoollu(slowpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
if ((pooldt%rcpoollu(armpc13) .gt. nzero) .and. &
    (pooldt%rcpoollu(armpc13) .lt. done)) then
  fract%d13cpool_armc13 = &
   (((1.0D0/((1.0D0/pooldt%rcpoollu(armpc13))-1.0D0))/pdb)-1.0D0)*1000.0D0
endif
!... Calculate weights for diagnostic output at finest timestep
!fract%d13cpcdbxpcdb = fract%d13cpool_cdbc13 * pooldt%curpoollu(cdbp)
!fract%d13cpmetlxpmetl = fract%d13cpool_metlc13 * pooldt%curpoollu(metlp)
!fract%d13cpstrlxpstrl = fract%d13cpool_strlc13 * pooldt%curpoollu(strlp)
!fract%d13cpslitxpslit = fract%d13cpool_slitc13 * pooldt%curpoollu(slitp)
!fract%d13cpslowxpslow = fract%d13cpool_slowc13 * pooldt%curpoollu(slowp)
!fract%d13cparmxparm = fract%d13cpool_armc13 * pooldt%curpoollu(armp)

!...Update poolpft_min for the isotope pools
tmpmin = nzero
!...Update poolpft_min for isotope pools based on current pool sizes
do n=npoolpft/2+1,npoolpft !(6,10)
   tcref=n-5
   if (poollt%rcpoolpft(n) .gt. nzero) then
     tmpmin = (poollt%rcpoolpft(n)) * poolcont%poolpft_min(tcref)
   else
     tmpmin = fract%rcassimfac * poolcont%poolpft_min(tcref)    
   endif
     poolcont%poolpft_min(n) = tmpmin

     !..update new diagnostics
     poollt%poolpftmin(tcref) = poolcont%poolpft_min(tcref)
     poollt%poolpftmin(n) = tmpmin
enddo


!.. Get grid cell pool ratio for (lp+wp+cdb+metl+strl) for fire

tmpval = ( (poollt%poolpft_lay(lp,1) &
            + poollt%poolpft_lay(wp,1) &
            + pooldt%poollu_lay(cdbp,1) &
            + pooldt%poollu_lay(metlp,1) &
            + pooldt%poollu_lay(strlp,1)) &
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

!if ( (poollt%poolpft(lpc13) .gt. nzero) .or. &
!     (poollt%poolpft(wpc13) .gt. nzero) .or. &
!     (pooldt%poollu(cdbpc13) .gt. nzero) .or. &
!     (pooldt%poollu(metlpc13) .gt. nzero) .or. &
!     (pooldt%poollu(strlpc13) .gt. nzero)) then

if (tmpval .gt. nzero) then
  fract%rcpoolfire = & 
     ( (poollt%poolpft_lay(lpc13,1) &
        + poollt%poolpft_lay(wpc13,1) &
        + pooldt%poollu_lay(cdbpc13,1) &
        + pooldt%poollu_lay(metlpc13,1) &
        + pooldt%poollu_lay(strlpc13,1)) &
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
     ( (poollt%poolpft_lay(lp,1) &
        + poollt%poolpft_lay(wp,1) &
        + pooldt%poollu_lay(cdbp,1) &
        + pooldt%poollu_lay(metlp,1) &
        + pooldt%poollu_lay(strlp,1)) &
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

   fract%poolemistotC = &
     ( (poollt%poolpft_lay(lp,1) &
        + poollt%poolpft_lay(wp,1) &
        + pooldt%poollu_lay(cdbp,1) &
        + pooldt%poollu_lay(metlp,1) &
        + pooldt%poollu_lay(strlp,1)) &
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

   fract%poolemisc13 = &
     ( (poollt%poolpft_lay(lpc13,1) &
        + poollt%poolpft_lay(wpc13,1) &
        + pooldt%poollu_lay(cdbpc13,1) &
        + pooldt%poollu_lay(metlpc13,1) &
        + pooldt%poollu_lay(strlpc13,1)) &
     + (poollt%poolpft_dgain(lpc13,1) &
        + poollt%poolpft_dgain(wpc13,1) &
        + pooldt%poollu_dgain(cdbpc13,1) &
        + pooldt%poollu_dgain(metlpc13,1) &
        + pooldt%poollu_dgain(strlpc13,1)) &
     - (poollt%poolpft_dloss(lpc13,1) &
        + poollt%poolpft_dloss(wpc13,1) &
        + pooldt%poollu_dloss(cdbpc13,1) &
        + pooldt%poollu_dloss(metlpc13,1) &
        + pooldt%poollu_dloss(strlpc13,1)) )
!else
!  poollt%rcpoolfire = (fract%rcassim/(fract%rcassim+1.0D0))  
endif

end subroutine c13_iso_calc
