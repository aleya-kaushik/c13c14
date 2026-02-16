!==========================================================================
subroutine pool_graze( &
   poolcont, &
   clim_lai, lai, &
   poollt, pooldt)
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
   pool_indx_canc13
use module_param, only: pool_param
use module_sib, only: &
   poold_type, pooll_type
use module_sibconst, only: &
   npoolcan, npoolpft, npoollu, &
   ntpool, npoolcanc13
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

! Local variables
integer(i4) :: cp, p, cpc13!, canpc13
integer(i4) :: m, mref, s
real(r8) :: cpooltot, pool_avail
real(r8) :: cpooltotc13, pool_availc13
real(r8) :: grzf, grzc, grzc13


!canpc13 = pool_indx_canc13-6

!----------------------------------------------
!Reset variables
poollt%resp_grz = dzero
poollt%resp_grzc13 = dzero
poollt%loss_grz = dzero
poollt%loss_grzc13 = dzero
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
   cp = pool_indx_can(p) !live pool indexed by canopy
   pool_avail = poollt%poolpft(cp) - poolcont%poolpft_min(cp)
   poollt%loss_grz(p) = MIN(pool_avail, &
        poollt%poolpft(cp) * grzf)
   poollt%poolpft_dloss(cp,1) = &
       poollt%poolpft_dloss(cp,1) + poollt%loss_grz(p)*dtsib
ENDDO

DO p=1, npoolcanc13
   cpc13 = pool_indx_canc13(p)-6
   pool_availc13 = poollt%poolpft(cpc13) - poolcont%poolpft_min(cpc13)
!   pool_availc13 = poollt%poolpft(cpc13) - poollt%poolpftmin_updated(cpc13)
   poollt%loss_grzc13(p) = MIN(pool_availc13, &
        poollt%poolpft(cpc13) * grzf)
   poollt%poolpft_dloss(cpc13,1) = &
       poollt%poolpft_dloss(cpc13,1) + poollt%loss_grzc13(p)*dtsib
ENDDO

!Determine where grazed carbon loss goes
!   (respired, removed, transferred to dead)
grzc = sum(poollt%loss_grz(1:3))
poollt%resp_grz = grzc * poolcont%graze_trans(1)
grzc13 = sum(poollt%loss_grzc13(1:3))
poollt%resp_grzc13 = grzc13 * poolcont%graze_trans(1)
!do m=npoolpft+1,ntpool !goes from (6,11)
!   mref=m-npoolpft !(m-5, i.e. 1-6 with npoolpft=5)
! with C-13 addition, npoolpft numbering is different
! npoolpft is now 10 [used to be 5], ntpool is 22 [used to be 11]
! live pool indexing: 1-5, 12-16
! dead pool indexing: 6-11, 17-22 (1-6,7-12)
do m=npoolpft/2+1,ntpool/2 ! goes from (6,11) ntpool dead pools
   mref=m-npoolpft/2 ! (m-5, i.e. 1-6 npoollu dead pools, with npoolpft=10)
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

do m=npoolpft+2+npoolpft/2,ntpool ! goes from (17,22) ntpool dead pools
   mref=m-npoolpft ! (m-10, 7-12 npoollu dead pools, with npoolpft=10)
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

end subroutine pool_graze
