!======================================================================
subroutine pool_auto_tran( poolcont, &
     daylen, daylendt, daylenmax, &
     tc, pawfrw, rootf_lay, poollt, &
     gain_transl_lay, poollu_dgain)
!======================================================================

! Description
! ------------
! Calculates the autotrophic transfer from the live pools. 
!
use kinds
use module_param, only: pool_param
use module_poolinfo, only: & 
   pool_indx_lay, pool_indx_leaf, &
   pool_indx_leaf_c13
use module_sib, only: pooll_type
use module_sibconst, only: &
   npoolpft, npoollu, &
   ntpool, nsoil
use module_time, only: &
   dtisib, dtsib, steps_per_day

implicit none

!...input variables
type(pool_param), intent(in) :: poolcont
real(r4), intent(in) :: daylenmax
real(r8), intent(in) :: daylen, daylendt
real(r8), intent(in) :: tc, pawfrw
real(r8), dimension(nsoil), intent(in) :: rootf_lay
type(pooll_type), intent(inout) :: poollt
real(r8), dimension(npoollu,nsoil), intent(inout) :: &
    gain_transl_lay, poollu_dgain

!...local variables
integer(byte) :: lp, lpc13
real(r8) :: qt, tfraclp, tfraclpc13
real(r8), dimension(npoolpft) :: tfrac !transfer fraction
real(r8), dimension(npoolpft,nsoil) :: &
   poolpft_avail_lay, & !available pool carbon (mol C/m2)
   tloss_lay   !pool transfer loss (mol C/m2/timestep)
real(r8), dimension(npoollu,nsoil) :: lutemp_gainl_lay

!...misc values
integer(i4) :: n,s,nref,tcref
integer(i4) :: m,mref

!-----------------------------------------
!...Reset transfer variables
lp = pool_indx_leaf
lpc13 = pool_indx_leaf_c13-6
tfrac(:) = dzero
poolpft_avail_lay(:,:) = dzero
tloss_lay(:,:) = dzero
lutemp_gainl_lay(:,:) = dzero
poollt%tfl_daylen = dzero
poollt%tfl_freeze = dzero
poollt%tfl_dry = dzero
poollt%tfl_total = dzero
poollt%tfl_totalc13 = dzero
poollt%tf_turnover(:) = dzero
poollt%loss_trans_lay(:,:) = dzero
gain_transl_lay(:,:) = dzero

!...Only transfer if there is pool carbon available
IF (sum(poollt%poolpft(1:5)) .gt. sum(poolcont%poolpft_min(1:5))) THEN

do n=1,npoolpft/2 !(1,5)
   do s=1,pool_indx_lay(n) !1,pool_indx_lay(1,5) either 1 or 10
       poolpft_avail_lay(n,s) = poollt%poolpft_lay(n,s) &
         - poolcont%poolpft_min(n)  &
         + poollt%poolpft_dgain(n,s) &
         - poollt%poolpft_dloss(n,s)
       poollt%pftavail_lay(n,s) = poolpft_avail_lay(n,s)
   enddo
enddo

do n=npoolpft/2+1,npoolpft !(6,10)
   tcref=n-5
   nref=n+npoolpft/2+1 !(12,16)
   do s=1,pool_indx_lay(nref) !1,pool_indx_lay(12-16) either 1 or 10
        poolpft_avail_lay(n,s) = &
             poollt%rcpoolpft_lay(n,s) * poolpft_avail_lay(tcref,s)
!       poolpft_avail_lay(n,s) = poollt%poolpft_lay(n,s) & !npoolpft,nsoil
!        - poolcont%poolpft_min(n)  & !npoolpft
!!         - poollt%poolpftmin_updated(n) &
!         + poollt%poolpft_dgain(n,s) & !npoolpft,nsoil
!         - poollt%poolpft_dloss(n,s) !npoolpft,nsoil
       poollt%pftavail_lay(n,s) = poolpft_avail_lay(n,s)
   enddo
enddo


!-----------------------------------------------------------------
!-----Turnover Fractions------------------------------------------
!For leaves, this is supplemented by transfer from shortening days,
! freezing temperatures, and lack of water.
!For all other live pools, this is the only source of transfer.
do n=1, npoolpft/2 !1,5
   do s=1,pool_indx_lay(n) !1,pool_indx_lay(1,5) either 1 or 10
      tloss_lay(n,s) = (1 - poolcont%lresp_eff(n))  & !lresp_eff:npoolpft
           * poollt%krater_lay(n,s) * poolpft_avail_lay(n,s) * dtsib
            !krater_lay:npoolpft,nsoil

      if (poolpft_avail_lay(n,s) .gt. dzero) then
           tfrac(n) = tfrac(n) &
               + poollt%poolpft_flay(n,s) * tloss_lay(n,s) &
               / poolpft_avail_lay(n,s)
      endif
   enddo
   poollt%tf_turnover(n) = tfrac(n)*steps_per_day
enddo

do n=npoolpft/2+1,npoolpft !6,10
   tcref=n-5
   nref=n+npoolpft/2+1 !(12,16)
   do s=1,pool_indx_lay(nref) !1,pool_indx_lay(12,16) either 1 or 10
      !tloss_lay(n,s) = (1 - poolcont%lresp_eff(n))  &
      !     * poollt%krater_lay(n,s) * poolpft_avail_lay(n,s) * dtsib
      !tloss_lay(npoolpft,nsoil)
      tloss_lay(n,s) = poollt%rcpoolpft_lay(n,s) * tloss_lay(tcref,s) 
      if (poolpft_avail_lay(n,s) .gt. dzero) then !(npoolpft,nsoil)
      !     tfrac(n) = tfrac(n) &
      !         + poollt%poolpft_flay(n,s) * tloss_lay(n,s) &
      !         / poolpft_avail_lay(n,s)
          tfrac(n) = tfrac(tcref)
      endif
   enddo
   poollt%tf_turnover(n) = tfrac(n)*steps_per_day
enddo


!----Additional Leaf Pool Loss Fractions-----
!Daylength
if (daylendt .lt. dzero) then
    tfraclp = MIN(poolcont%lt_dmax, MAX(dzero, &
             poolcont%lt_dcoef * &
             (daylenmax - daylen)* &
             (daylenmax-poolcont%lt_dref)))
    tfrac(lp) = tfrac(lp) + tfraclp / steps_per_day
!    tfrac(lpc13) = tfrac(lpc13) + tfraclp / steps_per_day
    poollt%tfl_daylen = tfraclp
else
    poollt%tfl_daylen = dzero    
endif

!Freezing
if (tc .lt. poolcont%lt_fref) then
    qt = 0.01 * (poolcont%lt_fref - tc)
    tfraclp = MIN(poolcont%lt_fmax, MAX(dzero, &
               poolcont%lt_fq10**qt - done))
    tfrac(lp) = tfrac(lp) + tfraclp / steps_per_day
!    tfrac(lpc13) = tfrac(lpc13) + tfraclp / steps_per_day
    poollt%tfl_freeze = tfraclp
else
    poollt%tfl_freeze = dzero
endif

!Phenology stage
tfrac(lp) = tfrac(lp) + poollt%tfl_pstage / steps_per_day
!tfrac(lpc13) = tfrac(lpc13) + poollt%tfl_pstage / steps_per_day

!Water Deficiency
if (pawfrw .lt. poolcont%lt_wref) then
   tfraclp = MAX(dzero, MIN(poolcont%lt_wmax, &
        poolcont%lt_wcoef*(poolcont%lt_wbase &
        ** (10.*(pawfrw - poolcont%lt_wref)) - 1.)))
   tfrac(lp) = tfrac(lp) + tfraclp / steps_per_day
!   tfrac(lpc13) = tfrac(lpc13) + tfraclp / steps_per_day
   poollt%tfl_dry = tfraclp
else
   poollt%tfl_dry = dzero
endif

!Combined factors
tfraclp = MIN(done, MAX(dzero, tfrac(lp)))
tfraclpc13 = tfraclp
!tfraclpc13 = MIN(done, MAX(dzero, tfrac(lpc13)))
poollt%tfl_total = tfraclp * steps_per_day
poollt%tfl_totalc13 = tfraclpc13 * steps_per_day
tloss_lay(lp,1) = poolpft_avail_lay(lp,1) * tfraclp
tloss_lay(lpc13,1) = poolpft_avail_lay(lpc13,1) * tfraclpc13

poollt%tfrac_lp = tfraclp
poollt%tfrac_lpc13 = tfraclpc13

!-----Live Pool Transfers------
do n=1,npoolpft/2 !1,5
    do s=1,pool_indx_lay(n) !1,pool_indx_lay(1,5) either 1 or 10
       !poollt%loss_trans_lay(npoolpft,nsoil)
       poollt%loss_trans_lay(n,s) = tloss_lay(n,s) * dtisib
       poollt%poolpft_dloss(n,s) = &
              poollt%poolpft_dloss(n,s) + tloss_lay(n,s)
    enddo !s=1,pool_indx_lay
enddo !n=1,npoolpft

do n=npoolpft/2+1,npoolpft !6,10
   nref=n+npoolpft/2+1 !12,16
   do s=1,pool_indx_lay(nref) !1,pool_indx_lay(12,16) either 1 or 10
       poollt%loss_trans_lay(n,s) = tloss_lay(n,s) * dtisib
       poollt%poolpft_dloss(n,s) = &
              poollt%poolpft_dloss(n,s) + tloss_lay(n,s)
   enddo !s=1,pool_index_lay
enddo !n=1,npoolpft


!----Transfer To Dead Pools---
!...n is the sending/from pool
!...m is the receiving/to pool
tloss_lay = tloss_lay * dtisib !convert to mol C/m2/s 
!tloss_lay(npoolpft,nsoil)??
!mref=1
!... transfer from live total C pools to dead total C pools
do n=1,npoolpft/2 !npoolpft was 5, now 10 with C-13 pools, (1,5) live total C
   !do m=npoolpft+1,ntpool !was (6,11) for dead pools 
   do m=npoolpft/2+1,ntpool/2 !now (6,11) 
      mref=m-npoolpft/2 !(m-5), i.e. 1-6 for npoollu dead pools
      if (poolcont%pool_trans_frac(n,m) > dzero) then
          !pool_trans_frac(ntpool,ntpool) (n=1-5,m=6,11)
          !.....transfer from single-layer canopy/soil to surface
          if ((pool_indx_lay(n) .eq. 1) .and. &
              (pool_indx_lay(m) .eq. 1)) then
                 lutemp_gainl_lay(mref,1) = & !lutemp_gainl_lay(npoollu,nsoil)
                   lutemp_gainl_lay(mref,1) + &
                   tloss_lay(n,1) * poolcont%pool_trans_frac(n,m) 
                   !pool_trans_frac(ntpool,ntpool) n=1,5 livepool m=1,6 deadpool
                   !tloss_lay(npoolpft,nsoil)

          !.....transfer from single-lay canopy/soil to soil
          elseif ((pool_indx_lay(n) .eq. 1) .and. &
                   (pool_indx_lay(m) .eq. nsoil)) then
                   do s=1,nsoil
                     lutemp_gainl_lay(mref,s) = &
                         lutemp_gainl_lay(mref,s) + &
                         tloss_lay(n,1) * rootf_lay(s) * &
                         poolcont%pool_trans_frac(n,m)
                  enddo

          !.....transfer from soil to soil
           elseif ((pool_indx_lay(n) .eq. nsoil) .and. &
                    (pool_indx_lay(m) .eq. nsoil)) then
                    do s=1,nsoil
                        lutemp_gainl_lay(mref,s) = &
                           lutemp_gainl_lay(mref,s) + &
                           tloss_lay(n,s) * &
                           poolcont%pool_trans_frac(n,m)
                    enddo
           else
                print*, 'Mismatching levels between pool transfers.'
                print*, 'Stopping in pool_trans_auto.'
                stop
            endif  !dead pool transfer
        endif !trans_frac > 0.
   enddo  !m=npoolpft/2+1,ntpool/2
enddo

!... transfer from live C-13 pools to dead C-13 pools
!do n=ntpool/2+1,ntpool/2+5 ! (12,16) live C-13 pools
do n=npoolpft/2+1,npoolpft !6,10 poolpft
   nref=n+npoolpft/2+1 !12,16 ntpool
   do m=npoolpft/2+1+ntpool/2,ntpool !(17,22) ntpool for C-13 dead pools 
      mref=m-npoolpft !(m-10, 7-12 npoollu dead pools with npoolpft=10)
      !pool_trans_frac(ntpool,ntpool)
      if (poolcont%pool_trans_frac(nref,m) > dzero) then
          !.....transfer from single-layer canopy/soil to surface
          !pool_indx_lay(ntpool)
          if ((pool_indx_lay(nref) .eq. 1) .and. &
              (pool_indx_lay(m) .eq. 1)) then
                 lutemp_gainl_lay(mref,1) = &
                   lutemp_gainl_lay(mref,1) + &
                   tloss_lay(n,1) * poolcont%pool_trans_frac(nref,m)
                  !tloss_lay(npoolpft,nsoil)
          !.....transfer from single-lay canopy/soil to soil
          elseif ((pool_indx_lay(nref) .eq. 1) .and. &
                   (pool_indx_lay(m) .eq. nsoil)) then
                   do s=1,nsoil
                     lutemp_gainl_lay(mref,s) = &
                         lutemp_gainl_lay(mref,s) + &
                         tloss_lay(n,1) * rootf_lay(s) * &
                         poolcont%pool_trans_frac(nref,m)
                  enddo

          !.....transfer from soil to soil
           elseif ((pool_indx_lay(nref) .eq. nsoil) .and. &
                    (pool_indx_lay(m) .eq. nsoil)) then
                    do s=1,nsoil
                        lutemp_gainl_lay(mref,s) = &
                           lutemp_gainl_lay(mref,s) + &
                           tloss_lay(n,s) * &
                           poolcont%pool_trans_frac(nref,m)
                    enddo
           else
                print*, 'Mismatching levels between pool transfers.'
                print*, 'Stopping in pool_trans_auto.'
                stop
            endif  !dead pool transfer
        endif !trans_frac > 0.
   enddo  !m=npoolpft+2+npoolpft/2,ntpool
enddo !n=1,npoolpft
ENDIF !pools > 0

!...Calculate dead pool gains
gain_transl_lay = lutemp_gainl_lay
poollu_dgain = poollu_dgain + lutemp_gainl_lay*dtsib

!...Additional diagnostics
poollt%autotran_dloss = poollt%poolpft_dloss

end subroutine pool_auto_tran
