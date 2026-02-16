!--------------------------------------------------------------
subroutine fire_interp_c13(indx, lon, lat, &
           sibg, physcont)

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
    fire_leaffc13, fire_stwdfc13, &
    fire_cdbfc13, fire_metlfc13, fire_strlfc13
use module_param, only: poolcon
use module_pparams, only: &
     mwc, mol_to_umol, &
     month_names, rpoolinitc3, &
     rpoolinitc4
use module_pftinfo, only: pft_num
use module_poolinfo
use module_sib, only: &
   gridcell_type
use module_sibconst, only: &
   npoolpft, npoollu, nsoil, &
   fireb_print, fireb_stop, fireb_thresh
use module_time, only: &
   month, day, year, &
   dtisib, dtsib, sec_tot, wt_daily
!use module_phosib, only: c4
use module_param, only: &
    phys_param, physcon

implicit none

!...input variables
integer(i4), intent(in) :: indx
real(r4), intent(in) :: lon, lat
type (gridcell_type), intent(inout) :: sibg
type(phys_param), intent(in) :: physcont

!...parameters
integer(i4), parameter :: isave=3
real(r8) :: dnzero=1.E-10

!...interpolation variables
real(r8) :: facsibdrv  ! scaling factor between driver data points
real(r8) :: totemisc13, curemis, pcemis

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

!...net balance variables
logical :: fb_err

!...misc variables
integer :: ido, idoc, l, p, s, myl
real(r8) :: myemis, tempemis
integer(i4), dimension(1) :: tempstore
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13
real(r8) :: c4_flag, isofactor

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
   sibg%l(l)%poollt%resp_firec13 = dzero
   sibg%l(l)%poollt%rmmd_firec13 = dzero
   sibg%l(l)%poollt%loss_fire_lay(6:10,:) = dzero
   sibg%l(l)%pooldt%loss_fire_lay(7:12,:) = dzero
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

   !if (c4_flag .EQ. dzero) then
   !  totemisc13 = rpoolinitc3*sibg%gprogt%firec * dtsib
   !else
   !  totemisc13 = rpoolinitc4*sibg%gprogt%firec * dtsib
   !endif

   !.... in order to use this below, need a grid-cell rcassim??
   !isofactor=fract%rcassim/(fract%rcassim+1.0D0)
   !totemisc13 = isofactor * sibg%gprogt%firec * dtsib

   ntpft = sibg%g_nlu !number of PFTs/land units in a given grid cell
   allocate(tpref(ntpft),tpnum(ntpft),tparea(ntpft))
   allocate(tpagb(ntpft),tpbiomass(ntpft))
   tpref(:) = sibg%l(1:ntpft)%ipft   
   tpnum(:) = pft_num(tpref)
   tparea(:) = sibg%l(1:ntpft)%larea
   tpagb(:) = dzero
   tpbiomass(:) = dzero

   do l=1, ntpft
     !...calculate total above-ground biomass (m-2) for C13 pools
      do p=npoolpft/2+1,npoolpft !6,10 npoolpft
         if (pool_indx_lay(p+npoolpft/2+1) .eq. 1) then !12,16 ntpool

           c4_flag = dble(physcon(p)%c4flag) 
            if (c4_flag .EQ. dzero) then
              totemisc13 = rpoolinitc3*sibg%gprogt%firec * dtsib
            else
              totemisc13 = rpoolinitc4*sibg%gprogt%firec * dtsib
            endif

           tpagb(l) = tpagb(l) + tparea(l) &
                 * (sibg%l(l)%poollt%poolpft(p) &
                 - sum(sibg%l(l)%poollt%poolpft_dloss(p,:)) &
                 !- poolcon(tpnum(l))%poolpft_min(p))
                 - sibg%l(l)%poollt%poolpftmin_updated(p)) 
         endif
      enddo
      do p=npoollu/2+1,npoollu !7,12 npoollu
         if (pool_indx_lay(p+npoolpft) .eq. 1) then !17,22 ntpool
            tpagb(l) = tpagb(l) + tparea(l) &
                 * (sibg%l(l)%pooldt%poollu(p) &
                    - sum(sibg%l(l)%pooldt%poollu_dloss(p,:)))
         endif
      enddo
      
      !...calculate total biomass (m-2)
      tpbiomass(l) = tpbiomass(l) + sum(sibg%l(l)%poollt%poolpft(6:10)) &
           - sum(sibg%l(l)%poollt%poolpft_dloss(6:10,:)) &
           !- sum(poolcon(tpnum(l))%poolpft_min(6:10))
           - sum(sibg%l(l)%poollt%poolpftmin_updated(6:10))
      tpbiomass(l) = tpbiomass(l) + sum(sibg%l(l)%pooldt%poollu(7:12)) &
           - sum(sibg%l(l)%pooldt%poollu_dloss(7:12,:))
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
      if (tpagb(myl) .lt. totemisc13) then
         idoc=MAX(ione, idoc-1)
      endif
   enddo
   ido=idoc
   tfarea = SUM(tparea(tsortref(1:ido)))

   !...remove carbon from top PFTs
   allocate(flosspft(ntpft,npoolpft/2))
   allocate(flosslu(ntpft,npoollu/2))
   flosspft(:,:) = dzero
   flosslu(:,:) = dzero
   allocate(flossb(ntpft))
   flossb(:) = dzero
   
   do l=ido,1,-1
      myl = tsortref(l)
      sibg%l(myl)%poollt%nd_fire = &
             sibg%l(myl)%poollt%nd_fire + wt_daily
      
      tempemis = totemisc13*(tparea(myl)/tfarea)
      curemis = tempemis

      !....remove C from leaf pool
      myemis = MIN(fire_leaffc13*tempemis, MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(lpc13) &
             - sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) &
            ! - poolcon(tpnum(myl))%poolpft_min(lpc13)))
             - sibg%l(myl)%poollt%poolpftmin_updated(lpc13)))
      flosspft(myl,lp) = myemis
      sibg%l(myl)%poollt%loss_fire_lay(lpc13,1) = myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(lpc13,1) &
           + myemis
      curemis = curemis - myemis

      !....remove C from wood pool
      myemis = MIN(fire_stwdfc13*tempemis, MAX(dzero, &
           sibg%l(myl)%poollt%poolpft(wpc13) &
             - sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) &
            ! - poolcon(tpnum(myl))%poolpft_min(wpc13)))
             - sibg%l(myl)%poollt%poolpftmin_updated(wpc13)))
      sibg%l(myl)%poollt%loss_fire_lay(wpc13,1) = myemis * dtisib
      sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) = &
           sibg%l(myl)%poollt%poolpft_dloss(wpc13,1) &
           + myemis
      flosspft(myl,wp) = myemis
      curemis = curemis - myemis

      !....remove C from metabolic litter
      myemis = MIN(fire_metlfc13*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(metlpc13) &
           - sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(metlpc13,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(metlpc13,1) &
           + myemis
      flosslu(myl,metlp) = myemis
      curemis = curemis - myemis

      !....remove C from structural litter
      myemis = MIN(fire_strlfc13*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(strlpc13) &
           - sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(strlpc13,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(strlpc13,1) &
           + myemis
      flosslu(myl,strlp) = myemis
      curemis = curemis - myemis

      !....remove C from coarse dead biomass
      myemis = MIN(fire_cdbfc13*tempemis, MAX(dzero, &
           sibg%l(myl)%pooldt%poollu(cdbpc13) &
           - sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1)))
      sibg%l(myl)%pooldt%loss_fire_lay(cdbpc13,1) = myemis * dtisib
      sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) = &
           sibg%l(myl)%pooldt%poollu_dloss(cdbpc13,1) &
           + myemis
      flosslu(myl,cdbp) = myemis
      curemis = curemis - myemis

      !...remove C from product
      IF (curemis .gt. dnzero) THEN
         myemis = MIN(curemis, MAX(dzero, &
              sibg%l(myl)%poollt%poolpft(ppc13) &
              - sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) &
             ! - poolcon(tpnum(myl))%poolpft_min(ppc13)))
              - sibg%l(myl)%poollt%poolpftmin_updated(ppc13)))
         sibg%l(myl)%poollt%loss_fire_lay(ppc13,1) = myemis * dtisib
         sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) = &
              sibg%l(myl)%poollt%poolpft_dloss(ppc13,1) &
              + myemis
         flosspft(myl,pp) = myemis
         curemis = curemis - myemis
      ENDIF
      
      !...remove C from roots
      IF (curemis .gt. dnzero) THEN
          myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%poollt%poolpft(frpc13) &
               - sum(sibg%l(myl)%poollt%poolpft_dloss(frpc13,:)) &
              ! - poolcon(tpnum(myl))%poolpft_min(frpc13)))
               - sibg%l(myl)%poollt%poolpftmin_updated(frpc13)))
          DO s=1,nsoil 
             sibg%l(myl)%poollt%loss_fire_lay(frpc13,s) = &
                  myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
             sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(frpc13,s) &
                  + myemis * sibg%l(myl)%poollt%poolpft_flay(frpc13,s)
          ENDDO
          flosspft(myl,frp) = myemis
          curemis = curemis - myemis
       ENDIF

       IF (curemis .gt. dnzero) THEN        
           myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%poollt%poolpft(crpc13) &
               - sum(sibg%l(myl)%poollt%poolpft_dloss(crpc13,:)) & 
              ! - poolcon(tpnum(myl))%poolpft_min(crpc13)))
               - sibg%l(myl)%poollt%poolpftmin_updated(crpc13)))
          DO s=1,nsoil 
             sibg%l(myl)%poollt%loss_fire_lay(crpc13,s) = &
                  myemis * dtisib * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
             sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) = &
                  sibg%l(myl)%poollt%poolpft_dloss(crpc13,s) &
                  + myemis * sibg%l(myl)%poollt%poolpft_flay(crpc13,s)
          ENDDO
          flosspft(myl,crp) = myemis
          curemis = curemis - myemis
       ENDIF

       !...remove C from soil litter
       IF (curemis .gt. dnzero) THEN
           myemis = MIN(curemis, MAX(dzero, &
                sibg%l(myl)%pooldt%poollu(slitpc13) &
                - sum(sibg%l(myl)%pooldt%poollu_dloss(slitpc13,:))))
           DO s=1,nsoil 
             sibg%l(myl)%pooldt%loss_fire_lay(slitpc13,s) = &
                  myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slitpc13,s) &
                  + myemis * sibg%l(myl)%pooldt%poollu_flay(slitpc13,s)
          ENDDO
          flosslu(myl,slitp) = myemis
         curemis = curemis - myemis
      ENDIF
      
      !...remove C from soil slow
      IF (curemis .gt. dnzero) THEN
          myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%pooldt%poollu(slowpc13) &
               - sum(sibg%l(myl)%pooldt%poollu_dloss(slowpc13,:))))
          DO s=1,nsoil 
             sibg%l(myl)%pooldt%loss_fire_lay(slowpc13,s) = &
                  myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(slowpc13,s) &
                  + myemis * sibg%l(myl)%pooldt%poollu_flay(slowpc13,s)
          ENDDO
          flosslu(myl,slowp) = myemis
          curemis = curemis - myemis
      ENDIF

      !...remove C from soil passive
      IF (curemis .gt. dnzero) THEN
          myemis = MIN(curemis, MAX(dzero, &
               sibg%l(myl)%pooldt%poollu(armpc13) &
               - sum(sibg%l(myl)%pooldt%poollu_dloss(armpc13,:))))
         DO s=1,nsoil 
             sibg%l(myl)%pooldt%loss_fire_lay(armpc13,s) = &
                  myemis * dtisib * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
             sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) = &
                  sibg%l(myl)%pooldt%poollu_dloss(armpc13,s) &
                  + myemis * sibg%l(myl)%pooldt%poollu_flay(armpc13,s)
          ENDDO
          flosslu(myl,armp) = myemis
          curemis = curemis - myemis
      ENDIF

      !...save what C was emitted but not removed
      flossb(myl) = curemis
      sibg%l(myl)%poollt%rmmd_firec13 = flossb(myl)*dtisib
      sibg%l(myl)%poollt%resp_firec13 = tempemis * dtisib
   enddo

   !...check to make sure all emissions have been taken from somewhere
   curemis = totemisc13 - sum(flosspft) - sum(flosslu) - sum(flossb)

   IF (curemis .gt. fireb_thresh) THEN
      fb_err = .true.
   else
      fb_err = .false.
   endif
   
   !...Print Results
   IF ((fb_err) .OR. (fireb_print)) THEN
      print*,''
      print*,'---FIRE C13 CARBON---'
      IF (fb_err) THEN
         print('(a)'),'!!Fire C13 Carbon Imbalance!!'
         print*,'Fire Emissions Mismatch (mol C/m2): ', curemis
      ENDIF

      print('(a,a,i3,a,i4)'), &
          '      Date: ', trim(month_names(month)), day, ', ', year
      print('(a,i6,2f8.2)'),  '      Point/Lon/Lat: ', indx, lon, lat
      print('(a,i14)'),       '      Current Model Second: ', sec_tot
      print('(a,2i12)'),      '      Current/Next Fire Second: ', &
              fire_seccur, fire_secnext

      print*,''
      if (c4_flag .EQ. dzero) then
          print('(a,f18.8)'),     '      Fire C13 Emissions (umol/m2/s): ', &
              rpoolinitc3*sibg%gprogt%firec*mol_to_umol
      else
          print('(a,f18.8)'),     '      Fire C13 Emissions (umol/m2/s): ', &
              rpoolinitc4*sibg%gprogt%firec*mol_to_umol
      endif
      print('(a,f18.8)'),     '      Time-Step C13 Losses (mol/m2):', totemisc13

      tempemis = sum(flosspft) + sum(flosslu)
      print('(a,f18.8)'),     '      SiB4 C Removal (mol/m2):    ', tempemis
      print('(a)'),           '         PFT    Loss          %-BioBurned   Fraction'
         do l=1,ido
            myl = tsortref(l)
            curemis = sum(flosslu(myl,:)) + sum(flosspft(myl,:))
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
   deallocate(flosspft,flosslu,flossb)

endif  !firec > 0

end subroutine fire_interp_c13

