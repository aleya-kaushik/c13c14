!--------------------------------------------------------------
subroutine fire_interp(indx, lon, lat, &
           sibg)

!--------------------------------------------------------------
!
! This subroutine interpolates the sibdrv fire emissions
!     between their read times
!
!--------------------------------------------------------------

use kinds
use module_io, only: &
   fire_step, fire_seccur, fire_secnext
use module_oparams, only: &
    fire_leaff, fire_stwdf, &
    fire_cdbf, fire_metlf, fire_strlf
use module_param, only: poolcon
use module_pparams, only: &
     mwc, mol_to_umol, month_names
use module_pftinfo, only: pft_num
use module_poolinfo
use module_sib, only: gridcell_type
use module_sibconst, only: &
   npoolpft, npoollu, nsoil, &
   fireb_print, fireb_stop, fireb_thresh
use module_time, only: &
   month, day, year, &
   dtisib, dtsib, sec_tot, wt_daily

implicit none

!...input variables
integer(i4), intent(in) :: indx
real(r4), intent(in) :: lon, lat
type (gridcell_type), intent(inout) :: sibg

!...parameters
integer(i4), parameter :: isave=3
real(r8) :: dnzero=1.E-10

!...interpolation variables
real(r8) :: facsibdrv  ! scaling factor between driver data points
real(r8) :: totemis, curemis, pcemis
real(r8) :: curemisc13

!...distribution variables
integer(i4) :: ntpft
integer(i4), dimension(:), allocatable :: tpref, tpnum
real(r4), dimension(:), allocatable :: tparea
real(r8), dimension(:), allocatable :: tpagb, tpagbtemp
real(r8), dimension(:), allocatable :: tpbiomass

real(r8) :: tfarea
integer(i4), dimension(:), allocatable :: tsortref
real(r8), dimension(:), allocatable :: flossb
real(r8), dimension(:,:), allocatable :: flosspft, flosslu
real(r8), dimension(:), allocatable :: flossbc13

!...net balance variables
logical :: fb_err

!...misc variables
integer :: ido, idoc, l, p, s, myl
real(r8) :: myemis, tempemis, tempemisc13
integer(i4), dimension(1) :: tempstore
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13
real(r8) :: rcpool, rcpoolf

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot
wp =  pool_indx_stwd
pp =  pool_indx_prod

lpc13 =  pool_indx_leaf_c13-6
frpc13 = pool_indx_froot_c13-6
crpc13 = pool_indx_croot_c13-6
wpc13 =  pool_indx_stwd_c13-6
ppc13 =  pool_indx_prod_c13-6

cdbp  = pool_indx_cdb-npoolpft/2
metlp = pool_indx_metl-npoolpft/2
strlp = pool_indx_strl-npoolpft/2
slitp = pool_indx_slit-npoolpft/2
slowp = pool_indx_slow-npoolpft/2
armp  = pool_indx_arm-npoolpft/2

cdbpc13  = pool_indx_cdb_c13 - npoolpft
metlpc13  = pool_indx_metl_c13 - npoolpft
strlpc13  = pool_indx_strl_c13 - npoolpft
slitpc13  = pool_indx_slit_c13 - npoolpft
slowpc13 = pool_indx_slow_c13 - npoolpft
armpc13  = pool_indx_arm_c13 - npoolpft

!-----------------------------------------------
! reset values
do l=1, sibg%g_nlu
   sibg%l(l)%poollt%resp_fire = dzero
   sibg%l(l)%poollt%rmmd_fire = dzero
   sibg%l(l)%poollt%resp_firec13 = dzero
   sibg%l(l)%poollt%rmmd_firec13 = dzero
   !sibg%l(l)%poollt%loss_fire_lay(1:5,:) = dzero
   !sibg%l(l)%pooldt%loss_fire_lay(1:6,:) = dzero
   sibg%l(l)%poollt%loss_fire_lay(:,:) = dzero
   sibg%l(l)%pooldt%loss_fire_lay(:,:) = dzero
enddo
   
! only continue if fire emissions are being used
IF (fire_step .le. izero) RETURN

! get scaling factors
facsibdrv = dble(fire_seccur-sec_tot) / dble(fire_step)

! only continue if fire emissions are valid at this time
IF ((facsibdrv .GT. 1) .OR. (facsibdrv .LT. 0)) THEN
   sibg%gprogt%firec = dzero
   sibg%gprogt%fireco2 = dzero
   RETURN
ENDIF

! interpolate fire C
sibg%gprogt%firec = facsibdrv*sibg%gprogt%firec1 + &
                    (1.-facsibdrv) * sibg%gprogt%firec2

! interpolate fire CO2
sibg%gprogt%fireco2 = facsibdrv*sibg%gprogt%fireco21 &
         + (1.-facsibdrv) * sibg%gprogt%fireco22

! distribute fire emissions per PFT/land unit
if (sibg%gprogt%firec .gt. dzero) then
   totemis = sibg%gprogt%firec * dtsib
   ntpft = sibg%g_nlu 
   ! number of land units/PFTs per cell -> =1 for Hyy when run as 1.0ENF
   allocate(tpref(ntpft),tpnum(ntpft),tparea(ntpft))
   allocate(tpagb(ntpft),tpbiomass(ntpft))
   tpref(:) = sibg%l(1:ntpft)%ipft   
   tpnum(:) = pft_num(tpref)
   tparea(:) = sibg%l(1:ntpft)%larea
   tpagb(:) = dzero
   tpbiomass(:) = dzero

   do l=1, ntpft
     !...calculate total above-ground biomass (m-2)
      do p=1, npoolpft/2 ! 1,5 npoolpft
         if (pool_indx_lay(p) .eq. 1) then !1,5 ntpool
            tpagb(l) = tpagb(l) + tparea(l) &
                 * (sibg%l(l)%poollt%poolpft(p) &
                 - sum(sibg%l(l)%poollt%poolpft_dloss(p,:)) &
                 - poolcon(tpnum(l))%poolpft_min(p))
         endif
      enddo
      do p=1, npoollu/2 !1,6 npoollu
         if (pool_indx_lay(p+npoolpft/2) .eq. 1) then ! 6,11 ntpool
            tpagb(l) = tpagb(l) + tparea(l) &
                 * (sibg%l(l)%pooldt%poollu(p) &
                    - sum(sibg%l(l)%pooldt%poollu_dloss(p,:)))
         endif
      enddo
      
      !...calculate total biomass (m-2)
      tpbiomass(l) = tpbiomass(l) + sum(sibg%l(l)%poollt%poolpft(1:5)) &
           - sum(sibg%l(l)%poollt%poolpft_dloss(1:5,:)) &
           - sum(poolcon(tpnum(l))%poolpft_min(1:5))
      tpbiomass(l) = tpbiomass(l) + sum(sibg%l(l)%pooldt%poollu(1:6)) &
           - sum(sibg%l(l)%pooldt%poollu_dloss(1:6,:))
   enddo

   !...rank PFTs by above ground biomass
   allocate(tpagbtemp(ntpft),tsortref(ntpft))
   tpagbtemp = tpagb
   do l=1,ntpft
       tempstore = maxloc(tpagbtemp)
       tsortref(l) = tempstore(1)
       tpagbtemp(tsortref(l)) = -1
   enddo
   
   !...check area of top PFTs
   ido=MIN(isave,ntpft)
   idoc=ido
   do l=1,ido
      myl = tsortref(l)
      if (tpagb(myl) .lt. totemis) then
         idoc=MAX(ione, idoc-1)
      endif
   enddo
   ido=idoc
   tfarea = SUM(tparea(tsortref(1:ido)))

   !...remove carbon from top PFTs
   allocate(flosspft(ntpft,npoolpft))
   allocate(flosslu(ntpft,npoollu))
   flosspft(:,:) = dzero
   flosslu(:,:) = dzero
   allocate(flossb(ntpft))
   flossb(:) = dzero
  
   allocate(flossbc13(ntpft))
   flossbc13(:) = dzero
 
   do l=ido,1,-1
      myl = tsortref(l)
      sibg%l(myl)%poollt%nd_fire = sibg%l(myl)%poollt%nd_fire + wt_daily
      
      tempemis = totemis*(tparea(myl)/tfarea)
      curemis = tempemis

      rcpoolf =( (sibg%l(myl)%poollt%poolpft(lpc13) &
                  + sibg%l(myl)%poollt%poolpft(wpc13) &
                  + sibg%l(myl)%pooldt%poollu(cdbpc13) &
                  + sibg%l(myl)%pooldt%poollu(metlpc13) &
                  + sibg%l(myl)%pooldt%poollu(strlpc13)) &
               + (sibg%l(myl)%poollt%poolpft_dgain(lpc13,1) &
                  + sibg%l(myl)%poollt%poolpft_dgain(wpc13,1) &
                  + sibg%l(myl)%pooldt%poollu_dgain(cdbpc13,1) &
                  + sibg%l(myl)%pooldt%poollu_dgain(metlpc13,1) &
                  + sibg%l(myl)%pooldt%poollu_dgain(strlpc13,1)) &
               - (sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) &
                  + sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) &
                  + sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) &                
                  + sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) &
                  + sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1)) )/ &
               ( (sibg%l(myl)%poollt%poolpft(lp) &
                  + sibg%l(myl)%poollt%poolpft(wp) &
                  + sibg%l(myl)%pooldt%poollu(cdbp) &
                  + sibg%l(myl)%pooldt%poollu(metlp) &
                  + sibg%l(myl)%pooldt%poollu(strlp)) &
               + (sibg%l(myl)%poollt%poolpft_dgain(lp,1) &
                  + sibg%l(myl)%poollt%poolpft_dgain(wp,1) &
                  + sibg%l(myl)%pooldt%poollu_dgain(cdbp,1) &
                  + sibg%l(myl)%pooldt%poollu_dgain(metlp,1) &
                  + sibg%l(myl)%pooldt%poollu_dgain(strlp,1)) &
               - (sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
                  + sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
                  + sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) &              
                  + sibg%l(myl)%pooldt%poollu_dloss(metlp,1) &
                  + sibg%l(myl)%pooldt%poollu_dloss(strlp,1)) )

      tempemisc13 = rcpoolf*tempemis
      curemisc13 = tempemisc13
      !rcpool = sibg%l(myl)%fract%rcassim/(sibg%l(myl)%fract%rcassim+1.0D0)

      !....remove C from leaf pool
      rcpool=( sibg%l(myl)%poollt%poolpft(lpc13) &
               + sibg%l(myl)%poollt%poolpft_dgain(lpc13,1) &
               - sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) )/ &
             ( sibg%l(myl)%poollt%poolpft(lp) &
               + sibg%l(myl)%poollt%poolpft_dgain(lp,1) &
               - sibg%l(myl)%poollt%poolpft_dloss(lp,1) )

      myemis = MIN(fire_leaff*tempemis, MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(lp) &
             - sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
             - poolcon(tpnum(myl))%poolpft_min(lp)))
      flosspft(myl,lp) = myemis
      sibg%l(myl)%poollt%loss_fire_lay(lp,1) = myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lp,1) &
           + myemis
      curemis = curemis - myemis

      flosspft(myl,lpc13) = rcpool*myemis
      sibg%l(myl)%poollt%loss_fire_lay(lpc13,1) = rcpool*myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) &
           + rcpool*myemis
      curemisc13 = curemisc13 - rcpool*myemis

      !....remove C from wood pool
      rcpool=( sibg%l(myl)%poollt%poolpft(wpc13) &
               + sibg%l(myl)%poollt%poolpft_dgain(wpc13,1) &
               - sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) )/ &
             ( sibg%l(myl)%poollt%poolpft(wp) &
               + sibg%l(myl)%poollt%poolpft_dgain(wp,1) &
               - sibg%l(myl)%poollt%poolpft_dloss(wp,1) )

      myemis = MIN(fire_stwdf*tempemis, MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(wp) &
             - sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
             - poolcon(tpnum(myl))%poolpft_min(wp)))
      sibg%l(myl)%poollt%loss_fire_lay(wp,1) = myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wp,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wp,1) &
           + myemis
      flosspft(myl,wp) = myemis
      curemis = curemis - myemis

      flosspft(myl,wpc13) = rcpool*myemis
      sibg%l(myl)%poollt%loss_fire_lay(wpc13,1) = rcpool*myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) &
           + rcpool*myemis
      curemisc13 = curemisc13 - rcpool*myemis

      !....remove C from metabolic litter
      rcpool=( sibg%l(myl)%pooldt%poollu(metlpc13) &
               + sibg%l(myl)%pooldt%poollu_dgain(metlpc13,1) &
               - sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) )/ &
             ( sibg%l(myl)%pooldt%poollu(metlp) &
               + sibg%l(myl)%pooldt%poollu_dgain(metlp,1) &
               - sibg%l(myl)%pooldt%poollu_dloss(metlp,1) )

      myemis = MIN(fire_metlf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(metlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(metlp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(metlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlp,1) &
           + myemis
      flosslu(myl,metlp) = myemis
      curemis = curemis - myemis

      flosslu(myl,metlpc13) = rcpool*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(metlpc13,1) = rcpool*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) &
           + rcpool*myemis
      curemisc13 = curemisc13 - rcpool*myemis

      !....remove C from structural litter
      rcpool=( sibg%l(myl)%pooldt%poollu(strlpc13) &
               + sibg%l(myl)%pooldt%poollu_dgain(strlpc13,1) &
               - sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) )/ &
             ( sibg%l(myl)%pooldt%poollu(strlp) &
               + sibg%l(myl)%pooldt%poollu_dgain(strlp,1) &
               - sibg%l(myl)%pooldt%poollu_dloss(strlp,1) )

      myemis = MIN(fire_strlf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(strlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(strlp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(strlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlp,1) &
           + myemis
      flosslu(myl,strlp) = myemis
      curemis = curemis - myemis

      flosslu(myl,strlpc13) = rcpool*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(strlpc13,1) = rcpool*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) &
           + rcpool*myemis
      curemisc13 = curemisc13 - rcpool*myemis

      !....remove C from coarse dead biomass
      rcpool=( sibg%l(myl)%pooldt%poollu(cdbpc13) &
               + sibg%l(myl)%pooldt%poollu_dgain(cdbpc13,1) &
               - sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) )/ &
             ( sibg%l(myl)%pooldt%poollu(cdbp) &
               + sibg%l(myl)%pooldt%poollu_dgain(cdbp,1) &
               - sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) )

      myemis = MIN(fire_cdbf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(cdbp) &
           - sibg%l(myl)%pooldt%poollu_dloss(cdbp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(cdbp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) &
           + myemis
      flosslu(myl,cdbp) = myemis
      curemis = curemis - myemis

      flosslu(myl,cdbpc13) = rcpool*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(cdbpc13,1) = rcpool*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) &
           + rcpool*myemis
      curemisc13 = curemisc13 - rcpool*myemis

      !...remove C from product
      IF (curemis .gt. dnzero) THEN
         rcpool=( sibg%l(myl)%poollt%poolpft(ppc13) &
                  + sibg%l(myl)%poollt%poolpft_dgain(ppc13,1) &
                  - sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) )/ &
                ( sibg%l(myl)%poollt%poolpft(pp) &
                  + sibg%l(myl)%poollt%poolpft_dgain(pp,1) &
                  - sibg%l(myl)%poollt%poolpft_dloss(pp,1) )

         myemis = MIN(curemis, MAX(dzero, &
              sibg%l(myl)%poollt%poolpft(pp) &
              - sibg%l(myl)%poollt%poolpft_dloss(pp,1) &
              - poolcon(tpnum(myl))%poolpft_min(pp)))
         sibg%l(myl)%poollt%loss_fire_lay(pp,1) = myemis * dtisib
         sibg%l(myl)%poollt%poolpft_dloss(pp,1) = &
              sibg%l(myl)%poollt%poolpft_dloss(pp,1) &
              + myemis
         flosspft(myl,pp) = myemis
         curemis = curemis - myemis

         flosspft(myl,ppc13) = rcpool*myemis
         sibg%l(myl)%poollt%loss_fire_lay(ppc13,1) = rcpool*myemis * dtisib
         sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) = &
             sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) &
             + rcpool*myemis
         curemisc13 = curemisc13 - rcpool*myemis
      ENDIF
      
      !...remove C from roots
      IF (curemis .gt. dnzero) THEN
         rcpool=( sibg%l(myl)%poollt%poolpft(frpc13) &
                  + sum(sibg%l(myl)%poollt%poolpft_dgain(frpc13,:)) &
                  - sum(sibg%l(myl)%poollt%poolpft_dloss(frpc13,:)) )/ &
                ( sibg%l(myl)%poollt%poolpft(frp) &
                  + sum(sibg%l(myl)%poollt%poolpft_dgain(frp,:)) &
                  - sum(sibg%l(myl)%poollt%poolpft_dloss(frp,:)) )

          myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%poollt%poolpft(frp) &
               - sum(sibg%l(myl)%poollt%poolpft_dloss(frp,:)) &
               - poolcon(tpnum(myl))%poolpft_min(frp)))

          DO s=1,nsoil 
             sibg%l(myl)%poollt%loss_fire_lay(frp,s) = &
                  myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(frp,s)
             sibg%l(myl)%poollt%poolpft_dloss(frp,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(frp,s) &
                  + myemis * sibg%l(myl)%poollt%poolpft_flay(frp,s)

             sibg%l(myl)%poollt%loss_fire_lay(frpc13,s) = &
                  rcpool*myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
             sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) &
                  + rcpool*myemis * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
          ENDDO
          flosspft(myl,frp) = myemis
          flosspft(myl,frpc13) = rcpool*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - rcpool*myemis
       ENDIF

       IF (curemis .gt. dnzero) THEN        
           rcpool=( sibg%l(myl)%poollt%poolpft(crpc13) &
                    + sum(sibg%l(myl)%poollt%poolpft_dgain(crpc13,:)) &
                    - sum(sibg%l(myl)%poollt%poolpft_dloss(crpc13,:)) )/ &
                  ( sibg%l(myl)%poollt%poolpft(crp) &
                    + sum(sibg%l(myl)%poollt%poolpft_dgain(crp,:)) &
                    - sum(sibg%l(myl)%poollt%poolpft_dloss(crp,:)) )

           myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%poollt%poolpft(crp) &
               - sum(sibg%l(myl)%poollt%poolpft_dloss(crp,:)) & 
               - poolcon(tpnum(myl))%poolpft_min(crp)))

          DO s=1,nsoil 
             sibg%l(myl)%poollt%loss_fire_lay(crp,s) = &
                  myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(crp,s)
             sibg%l(myl)%poollt%poolpft_dloss(crp,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(crp,s) &
                  + myemis * sibg%l(myl)%poollt%poolpft_flay(crp,s)

             sibg%l(myl)%poollt%loss_fire_lay(crpc13,s) = &
                  rcpool*myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
             sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) &
                  + rcpool*myemis * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
          ENDDO
          flosspft(myl,crp) = myemis
          flosspft(myl,crpc13) = rcpool*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - rcpool*myemis
       ENDIF

       !...remove C from soil litter
       IF (curemis .gt. dnzero) THEN
           rcpool=( sibg%l(myl)%pooldt%poollu(slitpc13) &
                 + sum(sibg%l(myl)%pooldt%poollu_dgain(slitpc13,:)) &
                 - sum(sibg%l(myl)%pooldt%poollu_dloss(slitpc13,:)) )/ &
                 ( sibg%l(myl)%pooldt%poollu(slitp) &
                 + sum(sibg%l(myl)%pooldt%poollu_dgain(slitp,:)) &
                 - sum(sibg%l(myl)%pooldt%poollu_dloss(slitp,:)) )

           myemis = MIN(curemis, MAX(dzero, &
                sibg%l(myl)%pooldt%poollu(slitp) &
                - sum(sibg%l(myl)%pooldt%poollu_dloss(slitp,:))))

           DO s=1,nsoil 
             sibg%l(myl)%pooldt%loss_fire_lay(slitp,s) = &
                  myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitp,s)
             sibg%l(myl)%pooldt%poollu_dloss(slitp,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slitp,s) &
                  + myemis * sibg%l(myl)%pooldt%poollu_flay(slitp,s)

             sibg%l(myl)%pooldt%loss_fire_lay(slitpc13,s) = &
                  rcpool*myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) &
                  + rcpool*myemis * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
          ENDDO
          flosslu(myl,slitp) = myemis
          flosslu(myl,slitpc13) = rcpool*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - rcpool*myemis
      ENDIF
      
      !...remove C from soil slow
      IF (curemis .gt. dnzero) THEN
           rcpool=( sibg%l(myl)%pooldt%poollu(slowpc13) &
                 + sum(sibg%l(myl)%pooldt%poollu_dgain(slowpc13,:)) &
                 - sum(sibg%l(myl)%pooldt%poollu_dloss(slowpc13,:)) )/ &
                 ( sibg%l(myl)%pooldt%poollu(slowp) &
                 + sum(sibg%l(myl)%pooldt%poollu_dgain(slowp,:)) &
                 - sum(sibg%l(myl)%pooldt%poollu_dloss(slowp,:)) )
          
           myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%pooldt%poollu(slowp) &
               - sum(sibg%l(myl)%pooldt%poollu_dloss(slowp,:))))

          DO s=1,nsoil 
             sibg%l(myl)%pooldt%loss_fire_lay(slowp,s) = &
                  myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowp,s)
             sibg%l(myl)%pooldt%poollu_dloss(slowp,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slowp,s) &
                  + myemis * sibg%l(myl)%pooldt%poollu_flay(slowp,s)

             sibg%l(myl)%pooldt%loss_fire_lay(slowpc13,s) = &
                  rcpool*myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) &
                  + rcpool*myemis * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
          ENDDO
          flosslu(myl,slowp) = myemis
          flosslu(myl,slowpc13) = rcpool*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - rcpool*myemis
      ENDIF

      !...remove C from soil passive
      IF (curemis .gt. dnzero) THEN
          rcpool=( sibg%l(myl)%pooldt%poollu(armpc13) &
                 + sum(sibg%l(myl)%pooldt%poollu_dgain(armpc13,:)) &
                 - sum(sibg%l(myl)%pooldt%poollu_dloss(armpc13,:)) )/ &
                 ( sibg%l(myl)%pooldt%poollu(armp) &
                 + sum(sibg%l(myl)%pooldt%poollu_dgain(armp,:)) &
                 - sum(sibg%l(myl)%pooldt%poollu_dloss(armp,:)) )
          myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%pooldt%poollu(armp) &
               - sum(sibg%l(myl)%pooldt%poollu_dloss(armp,:))))

         DO s=1,nsoil 
             sibg%l(myl)%pooldt%loss_fire_lay(armp,s) = &
                  myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(armp,s)
             sibg%l(myl)%pooldt%poollu_dloss(armp,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(armp,s) &
                  + myemis * sibg%l(myl)%pooldt%poollu_flay(armp,s)

             sibg%l(myl)%pooldt%loss_fire_lay(armpc13,s) = &
                  rcpool*myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) &
                  + rcpool*myemis * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
          ENDDO
          flosslu(myl,armp) = myemis
          flosslu(myl,armpc13) = rcpool*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - rcpool*myemis
      ENDIF

      !...save what C was emitted but not removed
      flossb(myl) = curemis
      sibg%l(myl)%poollt%rmmd_fire = flossb(myl)*dtisib
      sibg%l(myl)%poollt%resp_fire = tempemis * dtisib

      flossbc13(myl) = curemisc13
      sibg%l(myl)%poollt%rmmd_firec13 = flossbc13(myl)*dtisib
      sibg%l(myl)%poollt%resp_firec13 = tempemisc13 * dtisib

   enddo !cycling through top PFTs

   !...check to make sure all emissions have been taken from somewhere
   curemis = totemis - sum(flosspft(:,1:5)) - sum(flosslu(:,1:6)) - sum(flossb)

   IF (curemis .gt. fireb_thresh) THEN
      fb_err = .true.
   else
      fb_err = .false.
   endif
   
   !...Print Results
   IF ((fb_err) .OR. (fireb_print)) THEN
      print*,''
      print*,'---FIRE CARBON---'
      IF (fb_err) THEN
         print('(a)'),'!!Fire Carbon Imbalance!!'
         print*,'Fire Emissions Mismatch (mol C/m2): ', curemis
      ENDIF

      print('(a,a,i3,a,i4)'), &
          '      Date: ', trim(month_names(month)), day, ', ', year
      print('(a,i6,2f8.2)'),  '      Point/Lon/Lat: ', indx, lon, lat
      print('(a,i14)'),       '      Current Model Second: ', sec_tot
      print('(a,2i12)'),      '      Current/Next Fire Second: ', &
              fire_seccur, fire_secnext

      print*,''
      print('(a,i4)'),     '                                ntpft: ', &
           ntpft
      print*,''
      print('(a,f18.8)'),     '                            tfarea: ', &
           tfarea

      print*,''
      print('(a,f18.8)'),     '             sumtotal biomass(m-2): ', &
           sum(tpbiomass)
!      print*,''
!      print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!           'biomass(m-2): ', &
!           tpbiomass(1),' ',tpbiomass(2),' ',tpbiomass(3),' ',&
!           tpbiomass(4),' ',tpbiomass(5)
      print*,''
      print('(a,f18.8)'),     'sumtotal above-ground biomass(m-2): ', &
           sum(tpagb)
!      print*,''
!      print('(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a,f6.2)'),&
!           'above-ground biomass(m-2): ', &
!           tpagb(1),' ',tpagb(2),' ',tpagb(3),' ',&
!           tpagb(4),' ',tpagb(5)

      print*,''
      print('(a,f18.8)'),     '      Fire C Emissions (umol/m2/s): ', &
           sibg%gprogt%firec*mol_to_umol
      print('(a,f18.8)'),     '      Time-Step C Losses (mol/m2):', totemis

      tempemis = sum(flosspft(:,1:5)) + sum(flosslu(:,1:6))
      print('(a,f18.8)'),     '      SiB4 C Removal (mol/m2):    ', tempemis
      print('(a)'),           '         PFT    Loss          %-BioBurned   Fraction'
         do l=1,ido
            myl = tsortref(l)
            curemis = sum(flosslu(myl,1:6)) + sum(flosspft(myl,1:5))
            if (tpbiomass(myl) .gt. dzero) then
               pcemis = curemis/tpbiomass(myl)*100.
            else
               pcemis = dzero
            endif
            print('(a,i2,2f14.8,a,f6.2)'), '          ', tpref(myl),  &
               curemis, pcemis, '  ',tparea(myl)/tfarea
         enddo
      print('(a,f12.8)'),     '      Non-Matched C Respired: ', sum(flossb)

      IF (fb_err .AND. (fireb_stop)) STOP
   ENDIF  !print

   deallocate(tpref,tpnum,tparea)
   deallocate(tpagb,tpagbtemp)
   deallocate(flosspft,flosslu,flossb,flossbc13)

endif  !firec > 0

end subroutine fire_interp

