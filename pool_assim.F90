!==========================================================================
subroutine pool_assim( &
    sibpt, lonsib, latsib, pref, &
    hhti, hlti, shti, slti, gr_frac, &
    assim, c13assim, rstfac2, tm, &
    poollt, co2t, fract)
!==========================================================================
! -----------
! Updates live carbon pools from photosynthetic gains
! -----------

use kinds
use module_poolinfo, only: pool_indx_lay
use module_sib, only: &
    pooll_type, co2_type, &
    fract_type
use module_sibconst, only: npoolpft, npoollu
use module_time, only: dtsib

implicit none

!...input variables
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: lonsib, latsib
real(r4), intent(in) :: hhti, hlti, shti, slti
real(r4), dimension(npoolpft), intent(in) :: gr_frac
real(r8), intent(in) :: assim, c13assim, rstfac2, tm
type(pooll_type), intent(inout) :: poollt
type(co2_type), intent(in) :: co2t
type(fract_type), intent(in) :: fract

!...local variables
integer(i4) :: k, p, kref
integer(i4) :: ialloc, iallocc13
real(r8) :: deltac

!--------------------------------------------------------------------
!Reset variables
poollt%resp_grow = dzero
poollt%resp_growc13 = dzero
poollt%resp_nveg = dzero
poollt%resp_nvegc13 = dzero
poollt%gain_assim = dzero
poollt%loss_gresp = dzero
poollt%resp_nveg_tmp = dzero
poollt%resp_nvegc13_tmp = dzero

!Check for assimilation to allocate
IF (assim .le. dzero) RETURN

!Check phenology stage allocation fractions,
!...if not one, then in dormancy and respire
!...back out the 'assimilated' carbon
ialloc = nint(sum(poollt%alloc_phen(1:5)))
!iallocc13 = nint(sum(poollt%alloc_phen(6:10)))
IF (ialloc .ne. ione) THEN ! first 5 are total C, next 5 are C-13 live pools
  !do p=1, npoolpft/2
  !  poollt%resp_nveg_tmp(p) = assim
  !enddo
  !poollt%resp_nveg = sum(poollt%resp_nveg_tmp)
  poollt%resp_nveg = assim
  poollt%resp_nvegc13 = c13assim
  RETURN
ENDIF
!IF (iallocc13 .ne. ione) THEN
!  do p=npoolpft/2+1,npoolpft
!    poollt%resp_nvegc13_tmp(p) = assim 
!    poollt%resp_nvegc13_tmp(p) = c13assim    
!  enddo
!  poollt%resp_nvegc13 = sum(poollt%resp_nvegc13_tmp)
!  RETURN
!ENDIF

!Calculate allocation fraction adjustments
IF ((poollt%aadj_moist) .or. (poollt%aadj_temp)) THEN
    call pool_alloc( &
         sibpt, lonsib, latsib, pref, &
         hhti, hlti, shti, slti, &
         rstfac2, tm, &
         poollt%aadj_moist, poollt%aadj_temp, &
         poollt%alloc_phen, poollt%alloc_moist, &
         poollt%alloc_temp, poollt%alloc)
ELSE
    poollt%alloc_moist(:) = dzero
    poollt%alloc_temp(:) = dzero
    poollt%alloc(:) = poollt%alloc_phen
ENDIF

!Assign photosynthate to live pools
do p=1, npoolpft/2 !npoolpft is 10, first 5 are total C, next 5 are C-13
!  if (p <= 5) then ! total C pools
    deltac = co2t%assim*poollt%alloc(p)
    poollt%gain_assim(p) = deltac
    poollt%loss_gresp(p) = deltac*gr_frac(p)
    !pool_indx_lay(ntpool)
    do k=1,pool_indx_lay(p) !1,pool_indx_lay(1-5) is either 1 or 10
       poollt%poolpft_dgain(p,k) = poollt%poolpft_dgain(p,k) &
             + poollt%gain_assim(p)*poollt%poolpft_flay(p,k)*dtsib
       poollt%poolpft_dloss(p,k) = poollt%poolpft_dloss(p,k) &
             + poollt%loss_gresp(p)*poollt%poolpft_flay(p,k)*dtsib
    enddo
enddo
!  else !do the same calculations but for C-13 pools
do p=npoolpft/2+1,npoolpft !6,10
    kref=p+npoolpft/2+1 !12,16 ntpool
    deltac = fract%c13assim*poollt%alloc(p)
!    deltac = fract%rcpoolfac*co2t%assim*poollt%alloc(p)
!    deltac = (poollt%rcpoolpft(p))*assim*poollt%alloc(p)
    poollt%gain_assim(p) = deltac
    poollt%loss_gresp(p) = deltac*gr_frac(p)
    do k=1,pool_indx_lay(kref)
       !poolpft_dgain (npoolpft,nsoil)
       poollt%poolpft_dgain(p,k) = poollt%poolpft_dgain(p,k) &
             + poollt%gain_assim(p)*poollt%poolpft_flay(p,k)*dtsib
       poollt%poolpft_dloss(p,k) = poollt%poolpft_dloss(p,k) &
             + poollt%loss_gresp(p)*poollt%poolpft_flay(p,k)*dtsib

       if ( (poollt%poolpft_dloss(p,k) .gt. 10.) .or. &
            (poollt%poolpft_dloss(p,k) .lt. -10.)) then
          print*,' '
          print*,'code: pool_assim'
          print*,'p,k:',p,k
          print*,'poollt%poolpft_dloss(p,k): ',poollt%poolpft_dloss(p,k)
          print*,'poollt%loss_gresp(p): ',poollt%loss_gresp(p)
          print*,'poollt%poolpft_dloss(p-5,k): ',poollt%poolpft_dloss(p-5,k)
       endif

    enddo
!  endif
enddo
poollt%resp_grow = sum(poollt%loss_gresp(1:5))
poollt%resp_growc13 = sum(poollt%loss_gresp(6:10))


!...Additional diagnostics
poollt%assim_dgain = poollt%poolpft_dgain
poollt%assim_dloss = poollt%poolpft_dloss

end subroutine pool_assim


!==========================================================================
subroutine pool_alloc( &
    sibpt, lonsib, latsib, pref, &
    hhti, hlti, shti, slti,  &
    rstfac2, tm, &
    adj_moist, adj_temp, &
    alloc_phen, alloc_moist, alloc_temp, &
    alloc)
!==========================================================================
! -----------
!Calculates allocation fractions
!  -Phenology-based allocation from phenology stage
!  -Climate-based adjustments from Friedlingstein et al. (1999)
! -----------

use kinds
use module_oparams, only: &
    aadjustmin, &
    lftit, lftif, lgrw_min, &
    moistmult, wftit, wftif
use module_poolinfo
use module_sibconst, only: &
    npoolpft

implicit none

!...input variables
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: lonsib, latsib
real(r4), intent(in) :: hhti, hlti, shti, slti
real(r8), intent(in) :: rstfac2, tm
logical, intent(in) :: adj_moist, adj_temp
real(r8), dimension(npoolpft), intent(in) :: alloc_phen
real(r8), dimension(npoolpft), intent(inout) :: &
    alloc_moist, alloc_temp, alloc

!...local variables
integer(i4) :: lp,frp,crp,wp,pp,lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: npallow, npallowc13
real(r8) :: atot, atotc13, aadjust, aadjustc13, aadd, aaddc13
real(r8), dimension(npoolpft) :: alloc_now
logical, dimension(npoolpft) :: alloc_allow

real(r8) :: leaf_grw_frz  !scaling factor for leaf growth from temperature
                          !..leaf growth decreases when leaves cold/frozen
real(r8) :: wood_grw_frz  !scaling factor for wood growth from temperature
                          !..no wood growth for freezing temperatures
real(r8) :: wood_grw_temp !scaling factor for wood growth from temperature
                          !..wood growth decreases for low temperatures
real(r8) :: wood_grw_tot  !total wood growth temp scaling factor
real(r8) :: wood_grw_moist !scaling factor for wood growth from moisture
                           !..wood growth decreases with moisture stress

!...misc variables
integer(i4) :: p

!--------------------------------------------------------------------

!Set local variables
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

!Reset allocations
alloc_moist(:) = dzero
alloc_temp(:) = dzero
alloc(:) = dzero
alloc_now(:) = alloc_phen(:)
alloc_allow(:) = .false.

!--------------------------------
!-----ALLOCATION ADJUSTMENTS-----
!--------------------------------

!...scale the factors for temperature stress
IF ((adj_temp) .AND. &
    ((alloc_now(lp) .GT. 0.) .OR. (alloc_now(wp) .GT. 0.))) THEN
   !.....leaf growth decline for cold temperatures
   leaf_grw_frz = 1./(1.+EXP(lftif*(lftit-tm)))
   leaf_grw_frz = MAX(lgrw_min,leaf_grw_frz)

   !.....wood growth decline for cool temperatures
   wood_grw_temp = 1./(1.+EXP(slti*(hlti-tm))) / &
                   (1.+EXP(shti*(tm-hhti)))
   wood_grw_frz = 1./(1.+EXP(wftif*(wftit-tm)))
   wood_grw_tot = wood_grw_temp * wood_grw_frz

   aadjust = MAX(0.,alloc_now(wp)*(1.-wood_grw_tot))  + &
                  MAX(0.,alloc_now(crp)*(1.-wood_grw_tot)) + &
                  MAX(0.,alloc_now(lp)*(1.-leaf_grw_frz))

   aadjustc13 = MAX(0.,alloc_now(wpc13)*(1.-wood_grw_tot))  + &
                  MAX(0.,alloc_now(crpc13)*(1.-wood_grw_tot)) + &
                  MAX(0.,alloc_now(lpc13)*(1.-leaf_grw_frz))

   if (aadjust >= aadjustmin) then
      alloc_temp(lp) = -1. * MAX(0.,alloc_now(lp)*(1.-leaf_grw_frz))
      alloc_now(lp) = alloc_now(lp) + alloc_temp(lp)
      alloc_temp(wp) = -1. * MAX(0.,alloc_now(wp)*(1-wood_grw_tot))
      alloc_now(wp) = alloc_now(wp) + alloc_temp(wp)
      alloc_temp(crp)= -1. * MAX(0.,alloc_now(crp)*(1-wood_grw_tot))
      alloc_now(crp) = alloc_now(crp) + alloc_temp(crp)

   !if (aadjustc13 >= aadjustmin) then
      ! same calculations for C-13 pools
      alloc_temp(lpc13) = -1. * MAX(0.,alloc_now(lpc13)*(1.-leaf_grw_frz))
      alloc_now(lpc13) = alloc_now(lpc13) + alloc_temp(lpc13)
      alloc_temp(wpc13) = -1. * MAX(0.,alloc_now(wpc13)*(1-wood_grw_tot))
      alloc_now(wpc13) = alloc_now(wpc13) + alloc_temp(wpc13)
      alloc_temp(crpc13)= -1. * MAX(0.,alloc_now(crpc13)*(1-wood_grw_tot))
      alloc_now(crpc13) = alloc_now(crpc13) + alloc_temp(crpc13)

      npallow=0
      do p=1,npoolpft/2
         IF ((p .ne. lp) .and. (p .ne. wp) .and. &
            (alloc_now(p) .gt. 0.)) THEN
             alloc_allow(p) = .true.
             npallow=npallow+1
         ELSE
             alloc_allow(p) = .false.
         ENDIF
      enddo
      !.. same as above but for c13
      npallowc13=0
      do p=npoolpft/2+1,npoolpft 
         IF ((p .ne. lpc13) .and. (p .ne. wpc13) .and. &
             (alloc_now(p) .gt. 0.)) THEN
             alloc_allow(p) = .true.
             npallowc13=npallowc13+1
         ELSE
             alloc_allow(p) = .false.
         ENDIF
      enddo

      if (npallow .EQ. 0) then
         !adjustment fine roots
         alloc_temp(frp) = aadjust
         alloc_now(frp) = aadjust
      else
         aadd = aadjust/dble(npallow)
         do p=1,npoolpft/2
             IF (alloc_allow(p)) THEN
               alloc_temp(p) = &
                  alloc_temp(p) + aadd
               alloc_now(p) = alloc_now(p) + aadd
             ENDIF
         enddo
      endif
      !.. same as above but for c13
      if (npallowc13 .EQ. 0) then
         !adjustment fine roots
         alloc_temp(frpc13) = aadjustc13
         alloc_now(frpc13) = aadjustc13
      else
         aaddc13 = aadjustc13/dble(npallowc13)
         do p=npoolpft/2+1,npoolpft
             IF (alloc_allow(p)) THEN
               alloc_temp(p) = &
                  alloc_temp(p) + aaddc13
               alloc_now(p) = alloc_now(p) + aaddc13
             ENDIF
         enddo
      endif



      !check totals
      atot = sum(alloc_temp(1:5))
      if ((atot .lt. -0.01) .or. (atot .gt. 0.01)) then
           print*,'---Error with temperature adjustment allocation factors---'
           print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
           print*,'   Temperature Adjustment Total Allocation Factor: ',atot
           print*,'   Temperature Adjustment Contributions: '
           print*, alloc_temp(1:5)
           stop
       endif

       atot = sum(alloc_now(1:5))
       if ((atot .lt. 0.99) .or. (atot .gt. 1.01)) then
           print*,'---Error with total allocation following temp adjustments---'
           print*,'   SiB Point/Lon/Lat/PFT: ',sibpt, lonsib, latsib, pref
           print*,'   Total Allocation Factor: ',atot
           print*,'   Total Allocation Contributions: '
           print*, alloc_now(1:5)
           stop
       endif

      atotc13 = sum(alloc_temp(6:10))
      if ((atotc13 .lt. -0.01) .or. (atotc13 .gt. 0.01)) then
           print*,'---Error with C13 temperature adjustment allocation factors---'
           print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
           print*,'   Temperature Adjustment Total Allocation Factor: ',atotc13
           print*,'   Temperature Adjustment Contributions: '
           print*, alloc_temp(6:10)
           stop
       endif

       atotc13 = sum(alloc_now(6:10))
       if ((atotc13 .lt. 0.99) .or. (atotc13 .gt. 1.01)) then
           print*,'---Error with C13 total allocation following temp adjustments---'
           print*,'   SiB Point/Lon/Lat/PFT: ',sibpt, lonsib, latsib, pref
           print*,'   Total Allocation Factor: ',atotc13
           print*,'   Total Allocation Contributions: '
           print*, alloc_now(6:10)
           stop
       endif


   endif !(aadjust > adjustmin)
ENDIF !use_temp

!----------------
IF ((adj_moist) .AND. &
    ((alloc_now(lp) .GT. 0.) .OR. &
     (alloc_now(wp) .GT. 0.) .OR. & 
     (alloc_now(pp) .GT. 0.))) THEN

   !...scale the factors for moisture stress
   wood_grw_moist = MIN(1.0,MAX(0.,1. - rstfac2*moistmult))
   aadjust = MIN(1.0,MAX(0., &
       alloc_now(lp)*wood_grw_moist + &
       alloc_now(wp)*wood_grw_moist + &
       alloc_now(pp)*wood_grw_moist))

   aadjustc13 = MIN(1.0,MAX(0., &
       alloc_now(lpc13)*wood_grw_moist + &
       alloc_now(wpc13)*wood_grw_moist + &
       alloc_now(ppc13)*wood_grw_moist))


   if (aadjust >= aadjustmin) then
      alloc_moist(lp) = -1. * alloc_now(lp)*wood_grw_moist
      alloc_now(lp) = alloc_now(lp) + alloc_moist(lp)

      alloc_moist(wp) = -1. * alloc_now(wp)*wood_grw_moist
      alloc_now(wp) = alloc_now(wp) + alloc_moist(wp)

      alloc_moist(pp) = -1. * alloc_now(pp)*wood_grw_moist
      alloc_now(pp) = alloc_now(pp) + alloc_moist(pp)
   !same calculations for C-13 pools
   !if (aadjustc13 >= aadjustmin) then
      alloc_moist(lpc13) = -1. * alloc_now(lpc13)*wood_grw_moist
      alloc_now(lpc13) = alloc_now(lpc13) + alloc_moist(lp)

      alloc_moist(wpc13) = -1. * alloc_now(wpc13)*wood_grw_moist
      alloc_now(wpc13) = alloc_now(wpc13) + alloc_moist(wpc13)

      alloc_moist(ppc13) = -1. * alloc_now(ppc13)*wood_grw_moist
      alloc_now(ppc13) = alloc_now(ppc13) + alloc_moist(ppc13)

      !put in fine roots
      alloc_moist(frp) = aadjust
      alloc_now(frp) = alloc_now(frp) + alloc_moist(frp)
      alloc_moist(frpc13) = aadjustc13
      alloc_now(frpc13) = alloc_now(frpc13) + alloc_moist(frpc13)

       !check totals
       atot = sum(alloc_moist(1:5))
       if ((atot .lt. -0.01) .or. (atot .gt. 0.01)) then
            print*,'---Error with moisutre adjustment allocation factors---'
            print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
            print*,'   Moisture Adjustment Total Allocation Factor: ',atot
            print*,'   Moisture Adjustment Contributions: '
            print*, alloc_moist(1:5)
            stop
        endif

        atot = sum(alloc_now(1:5))
        if ((atot .lt. 0.99) .or. (atot .gt. 1.01)) then
            print*,'---Error with total allocation following moisture adjustments---'
            print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
            print*,'   Total Allocation Factor: ',atot
            print*,'   Total Allocation Contributions: '
            print*, alloc_now(1:5)
            stop
        endif

       atotc13 = sum(alloc_moist(6:10))
       if ((atotc13 .lt. -0.01) .or. (atotc13 .gt. 0.01)) then
            print*,'---Error with C13 moisutre adjustment allocation factors---'
            print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
            print*,'   Moisture Adjustment Total Allocation Factor: ',atotc13
            print*,'   Moisture Adjustment Contributions: '
            print*, alloc_moist(6:10)
            stop
        endif

        atotc13 = sum(alloc_now(6:10))
        if ((atotc13 .lt. 0.99) .or. (atotc13 .gt. 1.01)) then
            print*,'---Error with C13 total allocation following moisture adjustments---'
            print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
            print*,'   Total Allocation Factor: ',atotc13
            print*,'   Total Allocation Contributions: '
            print*, alloc_now(6:10)
            stop
        endif

   endif  !aadjust > aadjustmin

ENDIF !use_moist

!-----------------------------
!---SET ALLOCATION FACTORS---
alloc(:) = alloc_now(:)

!...check total allocation
atot = sum(alloc(1:5))
if ((atot .lt. 0.99) .or. (atot .gt. 1.01)) then
   print*,'---Error with final allocation factors---'
   print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
   print*,'   Total Allocation Factor: ',atot
   print*,'   Allocation Contributions: '
   print*, alloc(1:5)
   stop
endif

atotc13 = sum(alloc(6:10))
if ((atotc13 .lt. 0.99) .or. (atotc13 .gt. 1.01)) then
   print*,'---Error with C13 final allocation factors---'
   print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
   print*,'   Total Allocation Factor: ',atotc13
   print*,'   Allocation Contributions: '
   print*, alloc(6:10)
   stop
endif

!...check for negative allocations
do p=1,npoolpft/2
   if (alloc(p) < 0.) then
      if (alloc(p) > -.001) then
         alloc(p) = 0.
      else
          print*,'---Negative final allocation factors---'
          print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
          print*,'   Total Allocation Factor: ',atot
          print*,'   Allocation Contributions: '
          print*, alloc(1:5)
          stop
      endif
   endif
enddo

do p=npoolpft/2+1,npoolpft
   if (alloc(p) < 0.) then
      if (alloc(p) > -.001) then
         alloc(p) = 0.
      else
          print*,'---Negative C13 final allocation factors---'
          print*,'   SiB Point/Lon/Lat/PFT: ', sibpt, lonsib, latsib, pref
          print*,'   Total Allocation Factor: ',atotc13
          print*,'   Allocation Contributions: '
          print*, alloc(6:10)
          stop
      endif
   endif
enddo

end subroutine pool_alloc




