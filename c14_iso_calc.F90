!  c14_iso_calc calculates 13C and 12C 
!  for respiration fluxes and C13 pools
!
subroutine c14_iso_calc( &
           poolcont, &
           poollt, pooldt, &
           fract)

use kinds
use module_fractsib
use module_sib, only: &
    pooll_type, poold_type, &
    fract_type
use module_pparams,only: &
    pdb, c14taumean
use module_time, only: &
    dtsib, wt_daily, year
use module_poolinfo
use module_sibconst, only: &
    npoolpft, nsoil, npoollu, &
    varciso_switch, varco2_switch, &
    spinup
use module_phosib, only: &
    c4, gah2o, co2cap, &
    resp_cas
use module_param, only: &
   pool_param

implicit none

!Input Variables
type(fract_type), intent(inout) :: fract
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(pool_param), intent(inout) :: poolcont

integer(i4) :: n,s,sref,tcref,nref

integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc14,frpc14,crpc14,wpc14,ppc14
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc14, metlpc14, strlpc14, slitpc14, slowpc14, armpc14

real(r8) :: tottemp, totc14temp
real(r8) :: tmpmin
real(r8) :: tmpval,tmpvdr
real(r8) :: resp_casc14, resp_casc12
real(r8) :: nzero=1.E-14

!...reset variables
poollt%loss_raddecay_lay(:,:) = dzero
pooldt%loss_raddecay_lay(:,:) = dzero

!...with C14: ntpool=33, npoolpft=15, npoollu=18 
lp =  pool_indx_leaf !ntpool index 1
frp = pool_indx_froot !ntpool index 2
crp = pool_indx_croot !ntpool index 3
wp =  pool_indx_stwd !ntpool index 4
pp =  pool_indx_prod !ntpool index 5

lpc14 =  pool_indx_leaf_c14-(2*npoollu/3) !ntpool index 23, npoolpft index 11
frpc14 = pool_indx_froot_c14-(2*npoollu/3) !ntpool index 24, npoolpft index 12
crpc14 = pool_indx_croot_c14-(2*npoollu/3) !ntpool index 25, npoolpft index 13
wpc14 =  pool_indx_stwd_c14-(2*npoollu/3) !ntpool index 26, npoolpft index 14
ppc14 =  pool_indx_prod_c14-(2*npoollu/3) !ntpool index 27, npoolpft index 15

cdbp  = pool_indx_cdb-npoolpft/3 !ntpool index 6, npoollu index 1
metlp = pool_indx_metl-npoolpft/3 !ntpool index 7, npoollu index 2
strlp = pool_indx_strl-npoolpft/3 !ntpool index 8, npoollu index 3
slitp = pool_indx_slit-npoolpft/3 !ntpool index 9, npoollu index 4
slowp = pool_indx_slow-npoolpft/3 !ntpool index 10, npoollu index 5
armp  = pool_indx_arm-npoolpft/3 !ntpool index 11, npoollu index 6

cdbpc14  = pool_indx_cdb_c14 - npoolpft !ntpool index 28, npoollu index 13
metlpc14  = pool_indx_metl_c14 - npoolpft !ntpool index 29, npoollu index 14
strlpc14  = pool_indx_strl_c14 - npoolpft !ntpool index 30, npoollu index 15
slitpc14  = pool_indx_slit_c14 - npoolpft !ntpool index 31, npoollu index 16
slowpc14 = pool_indx_slow_c14 - npoolpft !ntpool index 32, npoollu index 17
armpc14  = pool_indx_arm_c14 - npoolpft !ntpool index 33, npoollu index 18


!! print statements to check poolpft_loss and poolpft_lay
!print*,' '
!print*,'code: c14_iso_calc'
!print*,'poolpft_dloss(7,1/2/3):',&
!    poollt%poolpft_dloss(7,1),poollt%poolpft_dloss(7,2),poollt%poolpft_dloss(7,3)
!print*,'poolpft_lay(7,1/2/3) :',&
!    poollt%poolpft_lay(7,1),poollt%poolpft_lay(7,2),poollt%poolpft_lay(7,3)
!print*,' '

!... decay only from all c14 pools... and decay only if not spinning up...
!... decay from live pools
IF (.not. spinup) THEN

   do n=1,npoolpft/3 !1,5 live totC 
      nref=n !1,5 ntpool
      do s=1,pool_indx_lay(nref) !ntpool
           poollt%loss_raddecay_lay(n,s) = dzero !npoolpft,nsoil
      enddo
   enddo
   
   do n=npoolpft/3+1,2*npoolpft/3 !6,10 live C13
      nref=n+npoolpft/3+1 !12,16 ntpool
      do s=1,pool_indx_lay(nref) !ntpool
           poollt%loss_raddecay_lay(n,s) = dzero !npoolpft,nsoil
      enddo
   enddo
   
   do n=2*npoolpft/3+1,npoolpft !11,15 live C14
      nref=n+2*npoolpft/3+2 !23,27 ntpool
    !only decay if pool mass is non-zero
    if (poollt%poolpft_lay(n,s) .gt. dzero) then 
      do s=1,pool_indx_lay(nref) !ntpool
           poollt%loss_raddecay_lay(n,s) = &
                 (1.0D0/c14taumean)*poollt%poolpft_lay(n,s) !npoolpft,nsoil
      enddo
    endif !only decay if pool mass is non-zero
   enddo
   
   !... decay from dead pools
   do n=1,npoollu/3 !1,6 dead totC 
      nref=n+npoolpft/3 !6,10 ntpool
      do s=1,pool_indx_lay(nref) !ntpool
           pooldt%loss_raddecay_lay(n,s) = dzero !npoollu,nsoil
      enddo
   enddo
   
   do n=npoollu/3+1,2*npoollu/3 !7,12 dead totC 
      nref=n+2*npoolpft/3 !17,22 ntpool
      do s=1,pool_indx_lay(nref) !ntpool
           pooldt%loss_raddecay_lay(n,s) = dzero !npoollu,nsoil
      enddo
   enddo
   
   do n=2*npoollu/3+1,npoollu !13,18 dead totC 
      nref=n+npoolpft !28,33 ntpool
     !only decay if pool mass is non-zero
     if (pooldt%poollu_lay(n,s) .gt. dzero) then
      do s=1,pool_indx_lay(nref) !ntpool
           pooldt%loss_raddecay_lay(n,s) = & 
                (1.0D0/c14taumean)*pooldt%poollu_lay(n,s) !npoollu,nsoil
      enddo
    endif !only decay if pool mass is non-zero
   enddo


ENDIF !spinup

!...Accumulate losses for pool updates
poollt%poolpft_dloss = poollt%poolpft_dloss & !(npoolpft,nsoil)
                      + poollt%loss_raddecay_lay*dtsib
pooldt%poollu_dloss = pooldt%poollu_dloss & !(npoollu,nsoil)
                      + pooldt%loss_raddecay_lay*dtsib

!...Initialize values
!fract%d13cpool_leafc14 = nzero
!fract%d13cpool_frootc14 = nzero
!fract%d13cpool_crootc14 = nzero
!fract%d13cpool_stwdc14 = nzero
!fract%d13cpool_prodc14 = nzero
!fract%d13cpool_cdbc14 = nzero
!fract%d13cpool_metlc14 = nzero
!fract%d13cpool_strlc14 = nzero
!fract%d13cpool_slitc14 = nzero
!fract%d13cpool_slowc14 = nzero
!fract%d13cpool_armc14 = nzero

!-------Calculate isotope signatures of pools and respiration fluxes--------

!...Isotope signatures of respiration fluxes
!... implement this after taking C14 std into account

!if (poollt%resp_auto .GT. nzero) then
!  fract%d13cresp_autoc14 = &
!      (((poollt%resp_autoc14/(poollt%resp_auto-poollt%resp_autoc14))/pdb)-1.0D0)*1000.0D0  
!endif

!if (poollt%resp_grow .GT. nzero) then
!  fract%d13cresp_growc14 = &
!      (((poollt%resp_growc14/(poollt%resp_grow-poollt%resp_growc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_leaf .GT. nzero) then
!  fract%d13cresp_leafc14 = &
!      (((poollt%resp_leafc14/(poollt%resp_leaf-poollt%resp_leafc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_mntn .GT. nzero) then
!  fract%d13cresp_mntnc14 = &
!      (((poollt%resp_mntnc14/(poollt%resp_mntn-poollt%resp_mntnc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_root .GT. nzero) then
!  fract%d13cresp_rootc14 = &
!      (((poollt%resp_rootc14/(poollt%resp_root-poollt%resp_rootc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_fire .GT. nzero) then
!  fract%d13cresp_firec14 = &
!      (((poollt%resp_firec14/(poollt%resp_fire-poollt%resp_firec14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_grz .GT. nzero) then
!  fract%d13cresp_grzc14 = &
!      (((poollt%resp_grzc14/(poollt%resp_grz-poollt%resp_grzc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_hrvst .GT. nzero) then
!  fract%d13cresp_hrvstc14 = &
!      (((poollt%resp_hrvstc14/(poollt%resp_hrvst-poollt%resp_hrvstc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (poollt%resp_nveg .GT. nzero) then
!  fract%d13cresp_nvegc14 = &
!      (((poollt%resp_nvegc14/(poollt%resp_nveg-poollt%resp_nvegc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (pooldt%resp_het .GT. nzero) then
!  fract%d13cresp_hetc14 = &
!      (((pooldt%resp_hetc14/(pooldt%resp_het-pooldt%resp_hetc14))/pdb)-1.0D0)*1000.0D0
!endif
!
!if (pooldt%resp_soil .GT. nzero) then
!  fract%d13cresp_soilc14 = &
!      (((pooldt%resp_soilc14/(pooldt%resp_soil-pooldt%resp_soilc14))/pdb)-1.0D0)*1000.0D0
!endif


tottemp = (poollt%resp_auto + pooldt%resp_het + poollt%resp_nveg &
           + poollt%resp_grz + poollt%resp_hrvst)

totc14temp = (poollt%resp_autoc14 + pooldt%resp_hetc14 + poollt%resp_nvegc14 &
           + poollt%resp_grzc14 + poollt%resp_hrvstc14)

fract%c14resptot=totc14temp
fract%c12resptot=(tottemp-totc14temp-fract%c13resptot)

!if (tottemp .GT. nzero) then
!  fract%d13cresp_totc14 = &
!            (((totc14temp/(tottemp-totc14temp))/pdb)-1.0D0)*1000.0D0
!endif

!!.. Update d13cca
!!.. following similar calculation in phosib for CO2

!if (varciso_switch .or. varco2_switch) then

resp_casc14 = poollt%resp_autoc14 + pooldt%resp_hetc14
resp_casc12 = resp_cas - resp_casc14 - (poollt%resp_autoc13 + pooldt%resp_hetc13)

poollt%resp_leafc12 = poollt%resp_leaf-poollt%resp_leafc14-poollt%resp_leafc13

!  print*,'FROM C14_ISO_CALC: '
!  print*,'d14cca_updated: ',fract%d14cca_updated
!  print*,'d14cassim: ',fract%d14cassim
!  print*,'c14alpha: ' ,fract%c14alpha
!  print*,'c14assim: ',fract%c14assim
!  print*,'c12assim: ',fract%c12assim
!  print*,'c14ca: ',fract%c14cca
!  print*,'c12ca: ',fract%c12cca
!  print*,'resp_casc14:',resp_casc14
!  print*,'resp_casc14:',resp_casc12
!  print*,'c14cm: ',fract%c14cm
! endif
!endif


!...Calculate rcpool for each pool/layer
!...uses tmpval to check for mass in total pools first

do n=2*npoolpft/3+1,npoolpft !(11,15) C14 live pools
   tcref=n-10 !(1,5) totC live pools
   !tsref=n-5 !ntpool for totC live pools
   sref=n+12 !23,27 ntpool

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
                print*,'code: c14_iso_calc'
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
!if ((poollt%rcpoolpft(lpc14) .gt. nzero) .and. & 
!    (poollt%rcpoolpft(lpc14) .lt. done)) then
!  fract%d13cpool_leafc14 = &
!   (((1.0D0/((1.0D0/poollt%rcpoolpft(lpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((poollt%rcpoolpft(frpc14) .gt. nzero) .and. &
!    (poollt%rcpoolpft(frpc14) .lt. done)) then
  !if (year .eq. 2013) then
  ! print*,'poollt%rcpoolpft(frpc14): ',poollt%rcpoolpft(frpc14)
  !endif
!  fract%d13cpool_frootc14 = &
!   (((1.0D0/((1.0D0/poollt%rcpoolpft(frpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((poollt%rcpoolpft(crpc14) .gt. nzero) .and. & 
!    (poollt%rcpoolpft(crpc14) .lt. done)) then
!  fract%d13cpool_crootc14 = &
!   (((1.0D0/((1.0D0/poollt%rcpoolpft(crpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((poollt%rcpoolpft(wpc14) .gt. nzero) .and. & 
!    (poollt%rcpoolpft(wpc14) .lt. done)) then
!  fract%d13cpool_stwdc14 = &
!   (((1.0D0/((1.0D0/poollt%rcpoolpft(wpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((poollt%rcpoolpft(ppc14) .gt. nzero) .and. & 
!    (poollt%rcpoolpft(ppc14) .lt. done)) then
!  fract%d13cpool_prodc14 = &
!   (((1.0D0/((1.0D0/poollt%rcpoolpft(ppc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!... Calculate weights for diagnostic output at finest timestep
!fract%d13cpleafxpleaf = fract%d13cpool_leafc14 * poollt%curpoolpft(lp)
!fract%d13cpfrootxpfroot = fract%d13cpool_frootc14 * poollt%curpoolpft(frp)
!fract%d13cpcrootxpcroot = fract%d13cpool_crootc14 * poollt%curpoolpft(crp)
!fract%d13cpstwdxpstwd = fract%d13cpool_stwdc14 * poollt%curpoolpft(wp)
!fract%d13cpprodxpprod = fract%d13cpool_prodc14 * poollt%curpoolpft(pp)

do n=2*npoollu/3+1,npoollu !(13,18) C14 dead pools
   tcref=n-12 !(1,6) totC dead pools
   !tsref=n !ntpool for totC dead pools
   sref=n+npoolpft !28,33 ntpool
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
!if ((pooldt%rcpoollu(cdbpc14) .gt. nzero) .and. &
!    (pooldt%rcpoollu(cdbpc14) .lt. done)) then
!  fract%d13cpool_cdbc14 = &
!   (((1.0D0/((1.0D0/pooldt%rcpoollu(cdbpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((pooldt%rcpoollu(metlpc14) .gt. nzero) .and. &
!    (pooldt%rcpoollu(metlpc14) .lt. done)) then
!  fract%d13cpool_metlc14 = &
!   (((1.0D0/((1.0D0/pooldt%rcpoollu(metlpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((pooldt%rcpoollu(strlpc14) .gt. nzero) .and. &
!    (pooldt%rcpoollu(strlpc14) .lt. done)) then
!  fract%d13cpool_strlc14 = &
!   (((1.0D0/((1.0D0/pooldt%rcpoollu(strlpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((pooldt%rcpoollu(slitpc14) .gt. nzero) .and. &
!    (pooldt%rcpoollu(slitpc14) .lt. done)) then
!  fract%d13cpool_slitc14 = &
!   (((1.0D0/((1.0D0/pooldt%rcpoollu(slitpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((pooldt%rcpoollu(slowpc14) .gt. nzero) .and. &
!    (pooldt%rcpoollu(slowpc14) .lt. done)) then
!  fract%d13cpool_slowc14 = &
!   (((1.0D0/((1.0D0/pooldt%rcpoollu(slowpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!if ((pooldt%rcpoollu(armpc14) .gt. nzero) .and. &
!    (pooldt%rcpoollu(armpc14) .lt. done)) then
!  fract%d13cpool_armc14 = &
!   (((1.0D0/((1.0D0/pooldt%rcpoollu(armpc14))-1.0D0))/pdb)-1.0D0)*1000.0D0
!endif
!... Calculate weights for diagnostic output at finest timestep
!fract%d13cpcdbxpcdb = fract%d13cpool_cdbc14 * pooldt%curpoollu(cdbp)
!fract%d13cpmetlxpmetl = fract%d13cpool_metlc14 * pooldt%curpoollu(metlp)
!fract%d13cpstrlxpstrl = fract%d13cpool_strlc14 * pooldt%curpoollu(strlp)
!fract%d13cpslitxpslit = fract%d13cpool_slitc14 * pooldt%curpoollu(slitp)
!fract%d13cpslowxpslow = fract%d13cpool_slowc14 * pooldt%curpoollu(slowp)
!fract%d13cparmxparm = fract%d13cpool_armc14 * pooldt%curpoollu(armp)

!...Update poolpft_min for the isotope pools
tmpmin = nzero
!...Update poolpft_min for isotope pools based on current pool sizes
do n=2*npoolpft/3+1,npoolpft !(11,15)
   tcref=n-10 !(1,5) totC live pools
   if (poollt%rcpoolpft(n) .gt. nzero) then
     tmpmin = (poollt%rcpoolpft(n)) * poolcont%poolpft_min(tcref)
   else
     tmpmin = fract%rcassimfacc14 * poolcont%poolpft_min(tcref)    
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

!if ( (poollt%poolpft(lpc14) .gt. nzero) .or. &
!     (poollt%poolpft(wpc14) .gt. nzero) .or. &
!     (pooldt%poollu(cdbpc14) .gt. nzero) .or. &
!     (pooldt%poollu(metlpc14) .gt. nzero) .or. &
!     (pooldt%poollu(strlpc14) .gt. nzero)) then

if (tmpval .gt. nzero) then
  fract%rcpoolfirec14 = & 
     ( (poollt%poolpft_lay(lpc14,1) &
        + poollt%poolpft_lay(wpc14,1) &
        + pooldt%poollu_lay(cdbpc14,1) &
        + pooldt%poollu_lay(metlpc14,1) &
        + pooldt%poollu_lay(strlpc14,1)) &
     + (poollt%poolpft_dgain(lpc14,1) &
        + poollt%poolpft_dgain(wpc14,1) &
        + pooldt%poollu_dgain(cdbpc14,1) &
        + pooldt%poollu_dgain(metlpc14,1) &
        + pooldt%poollu_dgain(strlpc14,1)) &
     - (poollt%poolpft_dloss(lpc14,1) &
        + poollt%poolpft_dloss(wpc14,1) &
        + pooldt%poollu_dloss(cdbpc14,1) &                
        + pooldt%poollu_dloss(metlpc14,1) &
        + pooldt%poollu_dloss(strlpc14,1)) )/ &
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

   fract%poolemisc14 = &
     ( (poollt%poolpft_lay(lpc14,1) &
        + poollt%poolpft_lay(wpc14,1) &
        + pooldt%poollu_lay(cdbpc14,1) &
        + pooldt%poollu_lay(metlpc14,1) &
        + pooldt%poollu_lay(strlpc14,1)) &
     + (poollt%poolpft_dgain(lpc14,1) &
        + poollt%poolpft_dgain(wpc14,1) &
        + pooldt%poollu_dgain(cdbpc14,1) &
        + pooldt%poollu_dgain(metlpc14,1) &
        + pooldt%poollu_dgain(strlpc14,1)) &
     - (poollt%poolpft_dloss(lpc14,1) &
        + poollt%poolpft_dloss(wpc14,1) &
        + pooldt%poollu_dloss(cdbpc14,1) &
        + pooldt%poollu_dloss(metlpc14,1) &
        + pooldt%poollu_dloss(strlpc14,1)) )
endif

end subroutine c14_iso_calc
