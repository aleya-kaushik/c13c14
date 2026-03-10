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
     mwc, mol_to_umol, &
     month_names, pdb
!use module_pftinfo, only: pft_num
use module_pftinfo
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
integer(i4), parameter :: isave=5
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
character(len=clen), dimension(:), allocatable :: tpname

real(r8) :: tfarea
integer(i4), dimension(:), allocatable :: tsortref
real(r8), dimension(:), allocatable :: flossb
real(r8), dimension(:,:), allocatable :: flosspft, flosslu
real(r8), dimension(:), allocatable :: flossbc13, flossbc14

!...net balance variables
logical :: fb_err

!...misc variables
integer :: ido, idoc, l, p, s, myl
integer :: loc_c3g, loc_c4g
real(r8) :: myemis, tempemis, tempemisc13, tempemisc14
real(r8) :: tmppooltot, tmppoolc13, tempemispoolc13, rcemispoolc13
real(r8) :: tempemispoolc14, rcemispoolc14
integer(i4), dimension(1) :: tempstore
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: lpc14,frpc14,crpc14,wpc14,ppc14
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13
integer(i4) :: cdbpc14, metlpc14, strlpc14, slitpc14, slowpc14, armpc14
!real(r8) :: rcpooltest

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot
wp =  pool_indx_stwd
pp =  pool_indx_prod

lpc13 =  pool_indx_leaf_c13-npoollu/3
frpc13 = pool_indx_froot_c13-npoollu/3
crpc13 = pool_indx_croot_c13-npoollu/3
wpc13 =  pool_indx_stwd_c13-npoollu/3
ppc13 =  pool_indx_prod_c13-npoollu/3

lpc14 =  pool_indx_leaf_c14-(2*(npoollu/3)) !ntpool index 23, npoolpft index 11
frpc14 = pool_indx_froot_c14-(2*(npoollu/3)) !ntpool index 24, npoolpft index 12
crpc14 = pool_indx_croot_c14-(2*(npoollu/3)) !ntpool index 25, npoolpft index 13
wpc14 =  pool_indx_stwd_c14-(2*(npoollu/3)) !ntpool index 26, npoolpft index 14
ppc14 =  pool_indx_prod_c14-(2*(npoollu/3)) !ntpool index 27, npoolpft index 15

cdbp  = pool_indx_cdb-npoolpft/3
metlp = pool_indx_metl-npoolpft/3
strlp = pool_indx_strl-npoolpft/3
slitp = pool_indx_slit-npoolpft/3
slowp = pool_indx_slow-npoolpft/3
armp  = pool_indx_arm-npoolpft/3

cdbpc13  = pool_indx_cdb_c13 - (2*(npoolpft/3))
metlpc13  = pool_indx_metl_c13 - (2*(npoolpft/3))
strlpc13  = pool_indx_strl_c13 - (2*(npoolpft/3))
slitpc13  = pool_indx_slit_c13 - (2*(npoolpft/3))
slowpc13 = pool_indx_slow_c13 - (2*(npoolpft/3))
armpc13  = pool_indx_arm_c13 - (2*(npoolpft/3))

cdbpc14  = pool_indx_cdb_c14 - npoolpft !ntpool index 28, npoollu index 13
metlpc14  = pool_indx_metl_c14 - npoolpft !ntpool index 29, npoollu index 14
strlpc14  = pool_indx_strl_c14 - npoolpft !ntpool index 30, npoollu index 15
slitpc14  = pool_indx_slit_c14 - npoolpft !ntpool index 31, npoollu index 16
slowpc14 = pool_indx_slow_c14 - npoolpft !ntpool index 32, npoollu index 17
armpc14  = pool_indx_arm_c14 - npoolpft !ntpool index 33, npoollu index 18


!! print statements to check poolpft_loss and poolpft_lay
!print*,' '
!print*,'code: fire_interp'
!print*,'poolpft_dloss(1,1/2/3):',&
!    sibg%l(l)%poollt%poolpft_dloss(1,1),sibg%l(l)%poollt%poolpft_dloss(1,2),sibg%l(l)%poollt%poolpft_dloss(1,3)
!print*,'poolpft_dloss(2,1/2/3):',&
!    sibg%l(l)%poollt%poolpft_dloss(2,1),sibg%l(l)%poollt%poolpft_dloss(2,2),sibg%l(l)%poollt%poolpft_dloss(2,3)
!print*,'poolpft(1) :', sibg%l(l)%poollt%poolpft(1)
!print*,'poolpft(2) :', sibg%l(l)%poollt%poolpft(2)
!print*,'poolpft(3) :', sibg%l(l)%poollt%poolpft(3)
!print*,' '
 


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
   allocate(tpname(ntpft))
   tpref(:) = sibg%l(1:ntpft)%ipft   
   tpnum(:) = pft_num(tpref)
   !print*,'tpnum: ',tpnum
   !print*,'pft_name: ',pft_name(tpnum)
   do l=1,ntpft
     tpname(l) = pft_name(tpnum(l))
   enddo
!   print*,'tpname: ',tpname
   tparea(:) = sibg%l(1:ntpft)%larea
   tpagb(:) = dzero
   tpbiomass(:) = dzero

   do l=1, ntpft
     !...calculate total above-ground biomass (m-2)
      do p=1, npoolpft/3 ! 1,5 npoolpft
         if (pool_indx_lay(p) .eq. 1) then !1,5 ntpool
            tpagb(l) = tpagb(l) + tparea(l) &
                 * (sibg%l(l)%poollt%poolpft(p) &
                 - sum(sibg%l(l)%poollt%poolpft_dloss(p,:)) &
                 - poolcon(tpnum(l))%poolpft_min(p))
         endif
      enddo
      do p=1, npoollu/3 !1,6 npoollu
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
   !print*,'tpagbtemp: ',tpagbtemp

   ! case: both c3g and c4g present
   if ( (ANY(tpname .eq. "c3g")) .and. &
        (ANY(tpname .eq. "c4g")) ) then
     !print*,'found both grasses'
     ! first find the array locs for grasses
     do l=1,ntpft
       if (tpname(l) .eq. "c3g") loc_c3g = l
       !print*,'assigned c3g: ',loc_c3g
       if (tpname(l) .eq. "c4g") loc_c4g = l
       !print*,'assigned c4g: ',loc_c4g
     enddo     
     ! sort grasses to the top, c4g first
     tsortref(1) = loc_c4g
     tsortref(2) = loc_c3g
     !print*,'tsortref after grass locs: ',tsortref(1),tsortref(2)
     ! set tpagbtemp (loc=grasses) to -1
     tpagbtemp(loc_c4g) = -1
     tpagbtemp(loc_c3g) = -1
     !print*,'tpagbtemp after grass loc: ',tpagbtemp 
     do l=3,ntpft
       ! add a statement to check for grass pft, 
       ! if grass, then continue
       !if (tpname(l) .eq. "c3g") CYCLE
       !if (tpname(l) .eq. "c4g") CYCLE  
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif ! end case selection: c3g & c4g

   !print*,'tsortref after grass loop: ',tsortref

   ! case: c3g only
   if ( (ANY(tpname .eq. "c3g")) .and. &
        (ALL(tpname .ne. "c4g")) ) then
    !print*,'found c3g'
     do l=1,ntpft
       if (tpname(l) .eq. "c3g") loc_c3g = l
     enddo
     tsortref(1) = loc_c3g
     tpagbtemp(loc_c3g) = -1
     do l=2,ntpft
       !if (tpname(l) .eq. "c3g") CYCLE
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif! end case selection: c3g

   ! case: c4g only
   if ( (ALL(tpname .ne. "c3g")) .and. &
        (ANY(tpname .eq. "c4g")) ) then
    !print*,'found c4g'
     do l=1,ntpft
       if (tpname(l) .eq. "c4g") loc_c4g = l
     enddo
     tsortref(1) = loc_c4g
     tpagbtemp(loc_c4g) = -1
     do l=2,ntpft
       !if (tpname(l) .eq. "c4g") CYCLE
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif! end case selection: c4g

   ! case: no grasses, default to original scheme
   if ( (ALL(tpname .ne. "c3g")) .and. &
        (ALL(tpname .ne. "c4g")) ) then 
     do l=1,ntpft
         tempstore = maxloc(tpagbtemp)
         tsortref(l) = tempstore(1)
         tpagbtemp(tsortref(l)) = -1
     enddo
   endif! end case selection: no grasses

   !print*,'tsortref: ',tsortref

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

   allocate(flossbc14(ntpft))
   flossbc14(:) = dzero
 
   do l=ido,1,-1
      myl = tsortref(l)
      sibg%l(myl)%poollt%nd_fire = sibg%l(myl)%poollt%nd_fire + wt_daily
      
      tempemis = totemis*(tparea(myl)/tfarea)
      curemis = tempemis

      tempemisc13 = sibg%l(myl)%fract%rcpoolfirec13*tempemis
      curemisc13 = tempemisc13
      
      tempemisc14 = sibg%l(myl)%fract%rcpoolfirec14*tempemis
      curemisc14 = tempemisc14

      !... Iniital pool totals come from c13_iso_calc 
      !... for sum of (lp+wp+cbdp+metlp+strlp)
      tmppooltot = sibg%l(myl)%fract%poolemistotC
      tmppoolc13 = sibg%l(myl)%fract%poolemisc13
      tmppoolc14 = sibg%l(myl)%fract%poolemisc14

      !....remove C from leaf pool
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

      flosspft(myl,lpc13) = sibg%l(myl)%poollt%rcpoolpft(lpc13)*myemis
      sibg%l(myl)%poollt%loss_fire_lay(lpc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(lpc13,1)*myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) &
           + sibg%l(myl)%poollt%rcpoolpft_lay(lpc13,1)*myemis
      curemisc13 = curemisc13 - sibg%l(myl)%poollt%rcpoolpft(lpc13)*myemis

      flosspft(myl,lpc14) = sibg%l(myl)%poollt%rcpoolpft(lpc14)*myemis
      sibg%l(myl)%poollt%loss_fire_lay(lpc14,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(lpc14,1)*myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lpc14,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lpc14,1) &
           + sibg%l(myl)%poollt%rcpoolpft_lay(lpc14,1)*myemis
      curemisc14 = curemisc14 - sibg%l(myl)%poollt%rcpoolpft(lpc14)*myemis

      !....remove C from wood pool
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

      flosspft(myl,wpc13) = sibg%l(myl)%poollt%rcpoolpft(wpc13)*myemis
      sibg%l(myl)%poollt%loss_fire_lay(wpc13,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(wpc13,1)*myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) &
           + sibg%l(myl)%poollt%rcpoolpft_lay(wpc13,1)*myemis
      curemisc13 = curemisc13 - sibg%l(myl)%poollt%rcpoolpft(wpc13)*myemis

      flosspft(myl,wpc14) = sibg%l(myl)%poollt%rcpoolpft(wpc14)*myemis
      sibg%l(myl)%poollt%loss_fire_lay(wpc14,1) = &
           sibg%l(myl)%poollt%rcpoolpft_lay(wpc14,1)*myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wpc14,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wpc14,1) &
           + sibg%l(myl)%poollt%rcpoolpft_lay(wpc14,1)*myemis
      curemisc14 = curemisc14 - sibg%l(myl)%poollt%rcpoolpft(wpc14)*myemis

      !....remove C from metabolic litter
      myemis = MIN(fire_metlf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(metlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(metlp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(metlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlp,1) &
           + myemis
      flosslu(myl,metlp) = myemis
      curemis = curemis - myemis

      flosslu(myl,metlpc13) = sibg%l(myl)%pooldt%rcpoollu(metlpc13)*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(metlpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(metlpc13,1)*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(metlpc13,1)*myemis
      curemisc13 = curemisc13 - sibg%l(myl)%pooldt%rcpoollu(metlpc13)*myemis

      flosslu(myl,metlpc14) = sibg%l(myl)%pooldt%rcpoollu(metlpc14)*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(metlpc14,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(metlpc14,1)*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlpc14,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlpc14,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(metlpc14,1)*myemis
      curemisc14 = curemisc14 - sibg%l(myl)%pooldt%rcpoollu(metlpc14)*myemis

      !....remove C from structural litter
      myemis = MIN(fire_strlf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(strlp) &
           - sibg%l(myl)%pooldt%poollu_dloss(strlp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(strlp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlp,1) &
           + myemis
      flosslu(myl,strlp) = myemis
      curemis = curemis - myemis

      flosslu(myl,strlpc13) = sibg%l(myl)%pooldt%rcpoollu(strlpc13)*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(strlpc13,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(strlpc13,1)*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(strlpc13,1)*myemis
      curemisc13 = curemisc13 - sibg%l(myl)%pooldt%rcpoollu(strlpc13)*myemis

      flosslu(myl,strlpc14) = sibg%l(myl)%pooldt%rcpoollu(strlpc14)*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(strlpc14,1) = &
           sibg%l(myl)%pooldt%rcpoollu_lay(strlpc14,1)*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlpc14,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlpc14,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(strlpc14,1)*myemis
      curemisc14 = curemisc14 - sibg%l(myl)%pooldt%rcpoollu(strlpc14)*myemis

      !....remove C from coarse dead biomass
      myemis = MIN(fire_cdbf*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(cdbp) &
           - sibg%l(myl)%pooldt%poollu_dloss(cdbp,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(cdbp,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbp,1) &
           + myemis
      flosslu(myl,cdbp) = myemis
      curemis = curemis - myemis

      flosslu(myl,cdbpc13) = sibg%l(myl)%pooldt%rcpoollu(cdbpc13)*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(cdbpc13,1) = &
             sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc13,1)*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc13,1)*myemis
      curemisc13 = curemisc13 - sibg%l(myl)%pooldt%rcpoollu(cdbpc13)*myemis

      flosslu(myl,cdbpc14) = sibg%l(myl)%pooldt%rcpoollu(cdbpc14)*myemis
      sibg%l(myl)%pooldt%loss_fire_lay(cdbpc14,1) = &
             sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc14,1)*myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbpc14,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbpc14,1) &
           + sibg%l(myl)%pooldt%rcpoollu_lay(cdbpc14,1)*myemis
      curemisc14 = curemisc14 - sibg%l(myl)%pooldt%rcpoollu(cdbpc14)*myemis

      !...remove C from product
      IF (curemis .gt. dnzero) THEN

         tmppooltot = tmppooltot + sibg%l(myl)%poollt%curpoolpft(pp)
         tmppoolc13 = tmppoolc13 + sibg%l(myl)%poollt%curpoolpft(ppc13)
         tmppoolc14 = tmppoolc14 + sibg%l(myl)%poollt%curpoolpft(ppc14)

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

         flosspft(myl,ppc13) = sibg%l(myl)%poollt%rcpoolpft(ppc13)*myemis
         sibg%l(myl)%poollt%loss_fire_lay(ppc13,1) = &
               sibg%l(myl)%poollt%rcpoolpft_lay(ppc13,1)*myemis * dtisib
         sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) = &
             sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) &
             + sibg%l(myl)%poollt%rcpoolpft_lay(ppc13,1)*myemis
         curemisc13 = curemisc13 - sibg%l(myl)%poollt%rcpoolpft(ppc13)*myemis

         flosspft(myl,ppc14) = sibg%l(myl)%poollt%rcpoolpft(ppc14)*myemis
         sibg%l(myl)%poollt%loss_fire_lay(ppc14,1) = &
               sibg%l(myl)%poollt%rcpoolpft_lay(ppc14,1)*myemis * dtisib
         sibg%l(myl)%poollt%poolpft_dloss(ppc14,1) = &
             sibg%l(myl)%poollt%poolpft_dloss(ppc14,1) &
             + sibg%l(myl)%poollt%rcpoolpft_lay(ppc14,1)*myemis
         curemisc14 = curemisc14 - sibg%l(myl)%poollt%rcpoolpft(ppc14)*myemis

      ENDIF
      
      !...remove C from roots
      IF (curemis .gt. dnzero) THEN
         tmppooltot = tmppooltot + sibg%l(myl)%poollt%curpoolpft(frp)
         tmppoolc13 = tmppoolc13 + sibg%l(myl)%poollt%curpoolpft(frpc13)
         tmppoolc14 = tmppoolc14 + sibg%l(myl)%poollt%curpoolpft(frpc14)

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
                  sibg%l(myl)%poollt%rcpoolpft_lay(frpc13,s)*myemis &
                  * dtisib * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
             sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) &
                  + sibg%l(myl)%poollt%rcpoolpft_lay(frpc13,s)*myemis &
                  * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)

             sibg%l(myl)%poollt%loss_fire_lay(frpc14,s) = &
                  sibg%l(myl)%poollt%rcpoolpft_lay(frpc14,s)*myemis &
                  * dtisib * sibg%l(myl)%poollt%poolpft_flay(frpc14,s)
             sibg%l(myl)%poollt%poolpft_dloss(frpc14,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(frpc14,s) &
                  + sibg%l(myl)%poollt%rcpoolpft_lay(frpc14,s)*myemis &
                  * sibg%l(myl)%poollt%poolpft_flay(frpc14,s)
          ENDDO
          flosspft(myl,frp) = myemis
          flosspft(myl,frpc13) = sibg%l(myl)%poollt%rcpoolpft(frpc13)*myemis
          flosspft(myl,frpc14) = sibg%l(myl)%poollt%rcpoolpft(frpc14)*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - sibg%l(myl)%poollt%rcpoolpft(frpc13)*myemis
          curemisc14 = curemisc14 - sibg%l(myl)%poollt%rcpoolpft(frpc14)*myemis
       ENDIF

       IF (curemis .gt. dnzero) THEN        
         tmppooltot = tmppooltot + sibg%l(myl)%poollt%curpoolpft(crp)
         tmppoolc13 = tmppoolc13 + sibg%l(myl)%poollt%curpoolpft(crpc13)
         tmppoolc14 = tmppoolc14 + sibg%l(myl)%poollt%curpoolpft(crpc14)

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
                  sibg%l(myl)%poollt%rcpoolpft_lay(crpc13,s)*myemis &
                  * dtisib * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
             sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) &
                  + sibg%l(myl)%poollt%rcpoolpft_lay(crpc13,s)*myemis &
                  * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)

             sibg%l(myl)%poollt%loss_fire_lay(crpc14,s) = &
                  sibg%l(myl)%poollt%rcpoolpft_lay(crpc14,s)*myemis &
                  * dtisib * sibg%l(myl)%poollt%poolpft_flay(crpc14,s)
             sibg%l(myl)%poollt%poolpft_dloss(crpc14,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(crpc14,s) &
                  + sibg%l(myl)%poollt%rcpoolpft_lay(crpc14,s)*myemis &
                  * sibg%l(myl)%poollt%poolpft_flay(crpc14,s)
          ENDDO
          flosspft(myl,crp) = myemis
          flosspft(myl,crpc13) = sibg%l(myl)%poollt%rcpoolpft(crpc13)*myemis
          flosspft(myl,crpc14) = sibg%l(myl)%poollt%rcpoolpft(crpc14)*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - sibg%l(myl)%poollt%rcpoolpft(crpc13)*myemis
          curemisc14 = curemisc14 - sibg%l(myl)%poollt%rcpoolpft(crpc14)*myemis
       ENDIF

       !...remove C from soil litter
       IF (curemis .gt. dnzero) THEN
         tmppooltot = tmppooltot + sibg%l(myl)%pooldt%curpoollu(slitp)
         tmppoolc13 = tmppoolc13 + sibg%l(myl)%pooldt%curpoollu(slitpc13)
         tmppoolc14 = tmppoolc14 + sibg%l(myl)%pooldt%curpoollu(slitpc14)

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
                  sibg%l(myl)%pooldt%rcpoollu_lay(slitpc13,s)*myemis &
                  * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) &
                  + sibg%l(myl)%pooldt%rcpoollu_lay(slitpc13,s)*myemis &
                  * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)

             sibg%l(myl)%pooldt%loss_fire_lay(slitpc14,s) = &
                  sibg%l(myl)%pooldt%rcpoollu_lay(slitpc14,s)*myemis &
                  * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitpc14,s)
             sibg%l(myl)%pooldt%poollu_dloss(slitpc14,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slitpc14,s) &
                  + sibg%l(myl)%pooldt%rcpoollu_lay(slitpc14,s)*myemis &
                  * sibg%l(myl)%pooldt%poollu_flay(slitpc14,s)
          ENDDO
          flosslu(myl,slitp) = myemis
          flosslu(myl,slitpc13) = sibg%l(myl)%pooldt%rcpoollu(slitpc13)*myemis
          flosslu(myl,slitpc14) = sibg%l(myl)%pooldt%rcpoollu(slitpc14)*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - sibg%l(myl)%pooldt%rcpoollu(slitpc13)*myemis
          curemisc14 = curemisc14 - sibg%l(myl)%pooldt%rcpoollu(slitpc14)*myemis
      ENDIF
      
      !...remove C from soil slow
      IF (curemis .gt. dnzero) THEN
         tmppooltot = tmppooltot + sibg%l(myl)%pooldt%curpoollu(slowp)
         tmppoolc13 = tmppoolc13 + sibg%l(myl)%pooldt%curpoollu(slowpc13)
         tmppoolc14 = tmppoolc14 + sibg%l(myl)%pooldt%curpoollu(slowpc14)

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
                  sibg%l(myl)%pooldt%rcpoollu_lay(slowpc13,s)*myemis &
                  * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) &
                  + sibg%l(myl)%pooldt%rcpoollu_lay(slowpc13,s)*myemis &
                  * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)

             sibg%l(myl)%pooldt%loss_fire_lay(slowpc14,s) = &
                  sibg%l(myl)%pooldt%rcpoollu_lay(slowpc14,s)*myemis &
                  * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowpc14,s)
             sibg%l(myl)%pooldt%poollu_dloss(slowpc14,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slowpc14,s) &
                  + sibg%l(myl)%pooldt%rcpoollu_lay(slowpc14,s)*myemis &
                  * sibg%l(myl)%pooldt%poollu_flay(slowpc14,s)
          ENDDO
          flosslu(myl,slowp) = myemis
          flosslu(myl,slowpc13) = sibg%l(myl)%pooldt%rcpoollu(slowpc13)*myemis
          flosslu(myl,slowpc14) = sibg%l(myl)%pooldt%rcpoollu(slowpc14)*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - sibg%l(myl)%pooldt%rcpoollu(slowpc13)*myemis
          curemisc14 = curemisc14 - sibg%l(myl)%pooldt%rcpoollu(slowpc14)*myemis
      ENDIF

      !...remove C from soil passive
      IF (curemis .gt. dnzero) THEN
         tmppooltot = tmppooltot + sibg%l(myl)%pooldt%curpoollu(armp)
         tmppoolc13 = tmppoolc13 + sibg%l(myl)%pooldt%curpoollu(armpc13)
         tmppoolc14 = tmppoolc14 + sibg%l(myl)%pooldt%curpoollu(armpc14)

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
                  sibg%l(myl)%pooldt%rcpoollu_lay(armpc13,s)*myemis &
                  * dtisib * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) &
                  + sibg%l(myl)%pooldt%rcpoollu_lay(armpc13,s)*myemis &
                  * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)

             sibg%l(myl)%pooldt%loss_fire_lay(armpc14,s) = &
                  sibg%l(myl)%pooldt%rcpoollu_lay(armpc14,s)*myemis &
                  * dtisib * sibg%l(myl)%pooldt%poollu_flay(armpc14,s)
             sibg%l(myl)%pooldt%poollu_dloss(armpc14,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(armpc14,s) &
                  + sibg%l(myl)%pooldt%rcpoollu_lay(armpc14,s)*myemis &
                  * sibg%l(myl)%pooldt%poollu_flay(armpc14,s)
          ENDDO
          flosslu(myl,armp) = myemis
          flosslu(myl,armpc13) = sibg%l(myl)%pooldt%rcpoollu(armpc13)*myemis
          flosslu(myl,armpc14) = sibg%l(myl)%pooldt%rcpoollu(armpc14)*myemis
          curemis = curemis - myemis
          curemisc13 = curemisc13 - sibg%l(myl)%pooldt%rcpoollu(armpc13)*myemis
          curemisc14 = curemisc14 - sibg%l(myl)%pooldt%rcpoollu(armpc14)*myemis
      ENDIF

      !...save what C was emitted but not removed
      flossb(myl) = curemis
      sibg%l(myl)%poollt%rmmd_fire = flossb(myl)*dtisib
      sibg%l(myl)%poollt%resp_fire = tempemis * dtisib

      flossbc13(myl) = curemisc13
      sibg%l(myl)%poollt%rmmd_firec13 = flossbc13(myl)*dtisib

      flossbc14(myl) = curemisc14
      sibg%l(myl)%poollt%rmmd_firec14 = flossbc14(myl)*dtisib

      !... Orig. calculation based on GFED4 emissions & rcpoolfire only
      !sibg%l(myl)%poollt%resp_firec13 = tempemisc13 * dtisib
      !... Now updated with actual fire isotope ratio based on pools burned

      if (tmppooltot .GT. dnzero) then
          rcemispoolc13 = (tmppoolc13/tmppooltot)
          tempemispoolc13 = rcemispoolc13*tempemis
          sibg%l(myl)%poollt%resp_firec13 = tempemispoolc13 * dtisib
   
          sibg%l(myl)%fract%d13cemisfire_pool = &
           (((1.0D0/((1.0D0/rcemispoolc13)-1.0D0))/pdb)-1.0D0)*1000.0D0

          rcemispoolc14 = (tmppoolc14/tmppooltot)
          tempemispoolc14 = rcemispoolc14*tempemis
          sibg%l(myl)%poollt%resp_firec14 = tempemispoolc14 * dtisib

      !..notes: 1. post-processing weighting would be done with sum(fire_loss_lay)
      !........ 2. also in post-processing, assign NaN to this isotope value if
      !........    resp_fire < 0 for any timestep 
      endif
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
   deallocate(flosspft,flosslu,flossb,flossbc13,flossbc14)

endif  !firec > 0

end subroutine fire_interp

