!==========================================================================
subroutine pool_graze( &
   poolcont, &
   clim_lai, lai, &
   poollt, pooldt, fract)
!==========================================================================
! ----------------------------------------------------------------
! Calculates grazing, removing carbon from live above-ground pools 
!    and transferring it to dead pools
! ------------------------------------------------------------

use kinds
use module_oparams, only: &
   graze_cfracp, graze_cfracd, &
   graze_minlai, graze_climlai
use module_pparams, only: mol_to_umol
use module_poolinfo, only: &
   pool_indx_can, pool_indx_lay, &
   pool_indx_canc13, pool_indx_canc14
use module_param, only: pool_param
use module_sib, only: &
   poold_type, pooll_type, &
   fract_type
use module_sibconst, only: &
   npoolcan, npoolpft, npoollu, &
   ntpool, npoolcanc13, npoolcanc14
use module_time, only: &
   dtsib, dtisib, wt_daily

implicit none

! Input variables
type(pool_param), intent(in) :: poolcont
real(r8), intent(in) :: clim_lai, lai
!real(r8), dimension(npoolpft), intent(in) :: poolpft_min
!real(r4), dimension(npoollu+2), intent(in) :: graze_trans
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(fract_type), intent(in) :: fract

! Local variables
integer(i4) :: cp, p, cpc13, cpc14!, canpc13
integer(i4) :: m, mref, s, tcref
integer(i4) :: isorefc13, isorefc14
real(r8) :: cpooltot, pool_avail
real(r8) :: cpooltotc13, pool_availc13
real(r8) :: cpooltotc14, pool_availc14
real(r8) :: grzf, grzc, grzc13, grzc14


!canpc13 = pool_indx_canc13-6

!! print statements to check poolpft_loss and poolpft_lay
!print*,' '
!print*,'code: pool_graze'
!print*,'poolpft_dloss(7,1/2/3):',&
!    poollt%poolpft_dloss(7,1),poollt%poolpft_dloss(7,2),poollt%poolpft_dloss(7,3)
!print*,'poolpft_lay(7,1/2/3) :',&
!    poollt%poolpft_lay(7,1),poollt%poolpft_lay(7,2),poollt%poolpft_lay(7,3)
!print*,' '
 


!----------------------------------------------
!Reset variables
poollt%resp_grz = dzero
poollt%resp_grzc13 = dzero
poollt%resp_grzc14 = dzero
poollt%loss_grz = dzero
poollt%loss_grzc13 = dzero
poollt%loss_grzc14 = dzero

pooldt%gain_grz_lay = dzero


!Check to make sure canopy/leaf pools are large
!enough to support grazing
cpooltot = sum(poollt%poolpft(pool_indx_can)) !1,4,5 live pool & ntpool
pool_avail = sum(poolcont%poolpft_min(pool_indx_can))
!cpooltotc13 = sum(poollt%poolpft(pool_indx_canc13-6)) !6,9,10 live pool, 12,15,16 ntpool 
!pool_availc13 = sum(poolpft_min(pool_indx_canc13-6))
IF ((lai .le. graze_minlai) .or. &
    (cpooltot .lt. pool_avail)) RETURN

!Set grazing rate
IF (clim_lai .gt. graze_climlai) THEN
   grzf = graze_cfracp * wt_daily * dtisib
ELSE
   grzf = graze_cfracd * wt_daily * dtisib
ENDIF

!Increment number of grazing days
poollt%nd_grz = poollt%nd_grz + wt_daily

!----------------------------------------------
!Calculate grazing loss and live removal
DO p=1, npoolcan
   cp = pool_indx_can(p) !live totC pool indexed by canopy (1,4,5)
   cpc13 = pool_indx_can(p)+5 !(6,9,10) live C13 pool indexed by canopy
   cpc14 = pool_indx_can(p)+10 !(11,14,15) live C14 pool indexed by canopy

   pool_avail = poollt%poolpft(cp) - poolcont%poolpft_min(cp)
   !..for total C live pools
   poollt%loss_grz(p) = MIN(pool_avail, &
        poollt%poolpft(cp) * grzf)
   poollt%poolpft_dloss(cp,1) = &
       poollt%poolpft_dloss(cp,1) + poollt%loss_grz(p)*dtsib
 
   !...same as above for C13
   poollt%loss_grzc13(p) = poollt%rcpoolpft(cpc13) * poollt%loss_grz(p)
   poollt%poolpft_dloss(cpc13,1) = &
       poollt%poolpft_dloss(cpc13,1) + poollt%loss_grzc13(p)*dtsib
   !...same as above for C14
   poollt%loss_grzc14(p) = poollt%rcpoolpft(cpc14) * poollt%loss_grz(p)
   poollt%poolpft_dloss(cpc14,1) = &
       poollt%poolpft_dloss(cpc14,1) + poollt%loss_grzc14(p)*dtsib
ENDDO


!Determine where grazed carbon loss goes
!   (respired, removed, transferred to dead)
grzc = sum(poollt%loss_grz(1:3))
poollt%resp_grz = grzc * poolcont%graze_trans(1)
grzc13 = sum(poollt%loss_grzc13(1:3))
poollt%resp_grzc13 = grzc13 * poolcont%graze_trans(1)
grzc14 = sum(poollt%loss_grzc14(1:3))
poollt%resp_grzc14 = grzc14 * poolcont%graze_trans(1)
!do m=npoolpft+1,ntpool !goes from (6,11)
!   mref=m-npoolpft !(m-5, i.e. 1-6 with npoolpft=5)
! with C-13 addition, npoolpft numbering is different
! npoolpft is now 10 [used to be 5], ntpool is 22 [used to be 11]
! live pool indexing: 1-5, 12-16
! dead pool indexing: 6-11, 17-22 (1-6,7-12)
! with C-14 addition, npoolpft numbering is different
! npoolpft is now 15 [used to be 10], ntpool is 33 [used to be 10]
! live pool indexing: 1-5, 12-16, 23-27 ntpool
! dead pool indexing: 6-11, 17-22, 28-33 ntpool (1-6, 7-12, 13-18 npoollu)

do m=npoollu/3,ntpool/3 ! goes from (6,11) ntpool dead pools
   mref=m-npoolpft/3 ! (m-5, i.e. 1-6 npoollu dead pools, with npoolpft=15)
   if (poolcont%graze_trans(mref+2) .gt. dzero) then !3-8
       do s=1,pool_indx_lay(m) !1,pool_indx_lay(6,11) is either 1 or 10
           !pool_indx_lay(ntpool)
           pooldt%gain_grz_lay(mref,s) = grzc & !(npoollu,nsoil)
                 * poolcont%graze_trans(mref+2) &
                 * pooldt%poollu_flay(mref,s)  !(npoollu,nsoil)
       enddo
       pooldt%poollu_dgain(mref,:) = pooldt%poollu_dgain(mref,:) &
                 + pooldt%gain_grz_lay(mref,:)*dtsib
    endif
enddo

!...same as above for C13
do m=npoolpft+2,2*ntpool/3 ! goes from (17,22) ntpool dead pools
   mref=m-2*npoolpft/3 ! (m-10, 7-12 npoollu dead pools, with npoolpft=105
   if (poolcont%graze_trans(mref+2) .gt. dzero) then !9-14
       do s=1,pool_indx_lay(m) !pool_indx_lay(17,22) is ether 1 or 10
           pooldt%gain_grz_lay(mref,s) = grzc13 &
                 * poolcont%graze_trans(mref+2) &
                 * pooldt%poollu_flay(mref,s)
       enddo
       pooldt%poollu_dgain(mref,:) = pooldt%poollu_dgain(mref,:) &
                 + pooldt%gain_grz_lay(mref,:)*dtsib
    endif
enddo

!...same as above for C14
do m=ntpool-npoolpft/3,ntpool ! goes from (28,33) ntpool dead pools
   mref=m-npoolpft ! (m-15, 13-18 npoollu dead pools, with npoolpft=15)
   if (poolcont%graze_trans(mref+2) .gt. dzero) then !15-20
       do s=1,pool_indx_lay(m) !pool_indx_lay(17,22) is ether 1 or 10
           pooldt%gain_grz_lay(mref,s) = grzc14 &
                 * poolcont%graze_trans(mref+2) &
                 * pooldt%poollu_flay(mref,s)
       enddo
       pooldt%poollu_dgain(mref,:) = pooldt%poollu_dgain(mref,:) &
                 + pooldt%gain_grz_lay(mref,:)*dtsib
    endif
enddo

end subroutine pool_graze
