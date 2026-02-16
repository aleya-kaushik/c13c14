!======================================================================
subroutine pool_auto_resp( poolcont, &
     clim_assim, clim_lai, &
     assimd, c13assimd, lai, tc, td_lay, &
     rootf_lay, poollt, resp_soil, resp_soil_lay, &
     resp_soilc13, resp_soilc13_lay, fract)
!======================================================================

! Description
! ------------
! Calculates the autotrophic respiration from the live pools. 
!
! Notes
! -------
! Studies suggest auto resp may scale predominantly 
! with assimilation (e.g. Flexas et al., 2006; Kirschbaum 1988;
!    Meir, 2008; Molchanov, 2009)
!
use kinds
use module_param, only: &
   pool_param
use module_poolinfo, only: & 
   pool_indx_froot, pool_indx_croot, &
   pool_indx_leaf, pool_indx_lay, &
   pool_indx_froot_c13, pool_indx_croot_c13, &
   pool_indx_leaf_c13
use module_sib, only: &
   soil_type, pooll_type, &
   fract_type
use module_sibconst, only: &
   npoolpft, nsoil
use module_time, only: dtsib

implicit none


!...input variables
type(pool_param), intent(inout) :: poolcont
real(r8), intent(in) :: clim_assim, clim_lai
real(r8), intent(in) :: assimd, c13assimd, lai, tc
real(r8), dimension(nsoil), intent(in) :: &
      rootf_lay, td_lay
type(pooll_type), intent(inout) :: poollt
real(r8), intent(inout) :: resp_soil, resp_soilc13
real(r8), dimension(nsoil), intent(inout) :: &
      resp_soil_lay, resp_soilc13_lay
type(fract_type), intent(in) :: fract

!...local variables
real(r8) :: slope, yint
real(r8) :: qt, mtemp
real(r8) :: pool_valid, pool_updated
real(r8) :: temp_mrespr, temp_mrespl
real(r8) :: temp_mresprc13, temp_mresplc13
real(r8) :: tmpv,tmpv1,tmpv2
!real(r8) :: isofactor

!...misc values
integer(i4) :: n,s,nref,tcref,isoref
integer(i4) :: lp,frp,crp
integer(i4) :: lpc13,frpc13,crpc13

!...alias the pool indices
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot

lpc13 =  pool_indx_leaf_c13-6
frpc13 = pool_indx_froot_c13-6
crpc13 = pool_indx_croot_c13-6

!-----------------------------------------
!...Reset respiration variables
poollt%mcr_assim = done
poollt%mcr_freeze = done
poollt%mcr_hot = done
poollt%mcr_scale = done
poollt%mrr_freeze_lay(:) = done !nsoil
poollt%mrr_hot_lay(:) = done !nsoil
poollt%mrr_scale_lay(:) = done !nsoil
poollt%mrr_assim = done
poollt%mrr_freeze = dzero
poollt%mrr_hot = dzero
poollt%mrr_lai = done
poollt%mrr_scale = dzero
poollt%krater_lay(:,:) = dzero !(npoolpft,nsoil)
poollt%loss_mresp_lay(:,:) = dzero !(npoolpft,nsoil)

!...Only respire if pools are greater than 
!...required minimum
IF (sum(poollt%poolpft(1:5)) .GT. sum(poolcont%poolpft_min(1:5))) THEN !npoolpft

!-----------------------------------------------------------------
!-----Canopy Autotrophic Respiration Scaling Factors-----
!Assimilation Rate Scalar
if (assimd .lt. clim_assim*poolcont%cr_aml) then
   poollt%mcr_assim = poolcont%cr_amin
elseif (assimd .lt. clim_assim*poolcont%cr_amh) then
   slope = (poolcont%cr_amax - poolcont%cr_amin) &
       /(clim_assim*poolcont%cr_amh - clim_assim*poolcont%cr_aml)
   yint = poolcont%cr_amin - slope*clim_assim*poolcont%cr_aml
   poollt%mcr_assim = assimd*slope + yint
else
   poollt%mcr_assim = poolcont%cr_amax
endif

!Freeze Inhibition
poollt%mcr_freeze = MAX(poolcont%cr_fmin, MIN(1.0, &
           EXP(poolcont%cr_fmul*(tc - poolcont%cr_fref))))

!High Temperature Exponential
qt = 0.1 * (tc - poolcont%cr_href)
mtemp = poolcont%cr_hq10**qt
poollt%mcr_hot = MAX(done, &
    MIN(poolcont%cr_hmax, mtemp))

!Combined factors
poollt%mcr_scale = poollt%mcr_assim * poollt%mcr_freeze &
                 * poollt%mcr_hot

!----Soil Autotrophic Respiration Scaling Factors----
!Assimilation Rate Scalar
if (assimd .lt. clim_assim*poolcont%rrt_aml) then
   poollt%mrr_assim = poolcont%rrt_amin
elseif (assimd .lt. clim_assim*poolcont%rrt_amh) then
   slope = (poolcont%rrt_amax - poolcont%rrt_amin) &
       /(clim_assim*poolcont%rrt_amh - clim_assim*poolcont%rrt_aml)
   yint = poolcont%rrt_amin - slope*clim_assim*poolcont%rrt_aml
   poollt%mrr_assim = assimd*slope + yint
else
   poollt%mrr_assim = poolcont%rrt_amax
endif

!LAI
poollt%mrr_lai = MIN(poolcont%rrt_laimax, &
     MAX(poolcont%rrt_laimin, lai / clim_lai))

do s=1,nsoil
   !...Soil Freeze Inhibition
   poollt%mrr_freeze_lay = MAX(poolcont%rrt_fmin, MIN(1.0, &
        EXP(poolcont%rrt_fmul*(td_lay(s) - poolcont%rrt_fref))))
   poollt%mrr_freeze = poollt%mrr_freeze &
       + poollt%mrr_freeze_lay(s) * rootf_lay(s)

   !...Soil High Temp Exponential
    qt = 0.1 * (td_lay(s) - poolcont%rrt_href)
    mtemp = poolcont%rrt_hq10**qt
    poollt%mrr_hot_lay(s) = MAX(done, &
        MIN(poolcont%rrt_hmax, mtemp))
    poollt%mrr_hot = poollt%mrr_hot &
       + poollt%mrr_hot_lay(s) * rootf_lay(s)
enddo

!Combined Factors
poollt%mrr_scale_lay(:) = poollt%mrr_freeze_lay(:) &
   * poollt%mrr_hot_lay(:) * poollt%mrr_lai * poollt%mrr_assim

DO s=1,nsoil
   poollt%mrr_scale = poollt%mrr_scale  &
       + poollt%mrr_scale_lay(s) * rootf_lay(s)  
ENDDO


!-----------------------------------------------------------------
!-----Autotrophic Respiration Rate----
do n=1,npoolpft/2 !1,5
    isoref=n+5
    !...only resp if above min pool value
    pool_valid = poolcont%poolpft_min(n)
!    poollt%poolpftmin(n) = pool_valid
    pool_updated = poollt%poolpft(n) &
        + sum(poollt%poolpft_dgain(n,:)) &
        - sum(poollt%poolpft_dloss(n,:))
    IF (pool_updated .LE. pool_valid) CYCLE

    !...calculate loss rate
    if (pool_indx_lay(n) .eq. 1) then !1,5 ntpool
        !krater_lay(npoolpft,nsoil), k_rate(ntpool)
        poollt%krater_lay(n,1) = poollt%mcr_scale &
              * poolcont%k_rate(n)
        poollt%krater_lay(isoref,1) = poollt%mcr_scale &
              * poolcont%k_rate(isoref)
    else
        do s=1,pool_indx_lay(n)
           poollt%krater_lay(n,s) = poollt%mrr_scale_lay(s) &
               * poolcont%k_rate(n)
           poollt%krater_lay(isoref,s) = poollt%mrr_scale_lay(s) &
               * poolcont%k_rate(isoref)
        enddo
    endif

    !.....calculate/check maintenance resp
    do s=1,pool_indx_lay(n) !1,5 ntpool, indx_lay same for isoref
       temp_mrespr = poollt%poolpft_lay(n,s) * &
                     poollt%krater_lay(n,s) * &
                     poolcont%lresp_eff(n)
       temp_mrespl = temp_mrespr*dtsib
       pool_updated = pool_updated - temp_mrespl
       IF (pool_updated .LT. pool_valid) CYCLE

       poollt%loss_mresp_lay(n,s) = temp_mrespr
       poollt%poolpft_dloss(n,s) = temp_mrespl &
           + poollt%poolpft_dloss(n,s)

       !poollt%loss_mresp_lay(isoref,s) = fract%rcpoolfac*temp_mrespr
       !poollt%poolpft_dloss(isoref,s) = fract%rcpoolfac*temp_mrespl & !(npoolpft,nsoil)
       !    + poollt%poolpft_dloss(isoref,s)
!       temp_mresprc13 = poollt%poolpft_lay(isoref,s) * &
!                     poollt%krater_lay(isoref,s) * &
!                     poolcont%lresp_eff(isoref)
!       temp_mresplc13 = temp_mrespr*dtsib

       poollt%loss_mresp_lay(isoref,s) = poollt%rcpoolpft_lay(isoref,s)*temp_mrespr
       poollt%poolpft_dloss(isoref,s) = poollt%rcpoolpft_lay(isoref,s)*temp_mrespl &
           + poollt%poolpft_dloss(isoref,s)

       !if ( (poollt%loss_mresp_lay(isoref,s) .gt. 1) .or. &
       !     (poollt%loss_mresp_lay(isoref,s) .lt. -1) ) then
       if ( (poollt%poolpft_dloss(isoref,s) .gt. 10.) .or. &
            (poollt%poolpft_dloss(isoref,s) .lt. -10.)) then
         print*,'code: pool_auto_resp'
         print*,'isoref: ',isoref
         print*,'s: ',s
         print*,'temp_mresprc13: ',temp_mresprc13
         print*,'temp_mrespr: ',temp_mrespr
         print*,'poollt%loss_mresp_lay(isoref-5,s):',poollt%loss_mresp_lay(isoref-5,s)
         !print*,'poollt%loss_mresp_lay(isoref-5,s-1):',poollt%loss_mresp_lay(isoref-5,s-1)
         !print*,'poollt%loss_mresp_lay(isoref-5,s+1):',poollt%loss_mresp_lay(isoref-5,s+1)
         print*,'poollt%loss_mresp_lay(isoref,s): ',poollt%loss_mresp_lay(isoref,s)
         !print*,'poollt%loss_mresp_lay(isoref,s-1):',poollt%loss_mresp_lay(isoref,s-1)
         !print*,'poollt%loss_mresp_lay(isoref,s+1):',poollt%loss_mresp_lay(isoref,s+1)
         print*,'poollt%rcpoolpft_lay(isoref,s): ',poollt%rcpoolpft_lay(isoref,s)
         !print*,'poollt%rcpoolpft_lay(isoref,s-1):',poollt%rcpoolpft_lay(isoref,s-1)
         !print*,'poollt%rcpoolpft_lay(isoref,s+1):',poollt%rcpoolpft_lay(isoref,s+1)
         print*,'temp_mrespr: ',temp_mrespr

         tmpv=( poollt%poolpft_lay(isoref,s) &
                + poollt%poolpft_dgain(isoref,s) &
                - poollt%poolpft_dloss(isoref,s) )
         if (s .gt. 1) then
           tmpv1=( poollt%poolpft_lay(isoref,s-1) &
                + poollt%poolpft_dgain(isoref,s-1) &
                - poollt%poolpft_dloss(isoref,s-1) )
         endif
         tmpv2=( poollt%poolpft_lay(isoref,s+1) &
                + poollt%poolpft_dgain(isoref,s+1) &
                - poollt%poolpft_dloss(isoref,s+1) )
         print*,'tmpv(isoref,s):',tmpv
         if (s .gt. 1) then
            print*,'tmpv(isoref,s-1):',tmpv1
         endif
         print*,'tmpv(isoref,s+1):',tmpv2


         print*,'poollt%poolpft_lay(isoref,s):',poollt%poolpft_lay(isoref,s)
         if (s .gt. 1) then
           print*,'poollt%poolpft_lay(isoref,s-1):',poollt%poolpft_lay(isoref,s-1)
         endif
         print*,'poollt%poolpft_lay(isoref,s+1):',poollt%poolpft_lay(isoref,s+1)
         print*,'poollt%poolpft_dgain(isoref,s):',poollt%poolpft_dgain(isoref,s)
         if (s .gt. 1) then
           print*,'poollt%poolpft_dgain(isoref,s-1):',poollt%poolpft_dgain(isoref,s-1)
         endif
         print*,'poollt%poolpft_dgain(isoref,s+1):',poollt%poolpft_dgain(isoref,s+1)
         print*,'poollt%poolpft_dloss(isoref,s):',poollt%poolpft_dloss(isoref,s)
         if (s .gt. 1) then
           print*,'poollt%poolpft_dloss(isoref,s-1):',poollt%poolpft_dloss(isoref,s-1)
         endif
         print*,'poollt%poolpft_dloss(isoref,s+1):',poollt%poolpft_dloss(isoref,s+1)
       endif

    enddo !s=1,pool_indx_lay
    
enddo !n=1,npoolpft
!ENDIF !pools are greater than required minimum

!..same as above but for C13 pools
!do n=npoolpft/2+1,npoolpft !6,10 live pools
!    !...only resp if above min pool value
!    !...for C13, still use totalC as the trigger to cycle or not
!    tcref=n-5
!    pool_valid = poolcont%poolpft_min(tcref)
!!    pool_valid = poollt%poolpftmin_updated(n)
! !   poollt%poolpftmin(n) = poolcont%poolpft_min(n)
!    pool_updated = poollt%poolpft(tcref) &
!        + sum(poollt%poolpft_dgain(tcref,:)) &
!        - sum(poollt%poolpft_dloss(tcref,:))
!!    pool_valid = poolcont%poolpft_min(n)
!!    poollt%poolpftmin(n) = poolcont%poolpft_min(n)
!!    pool_updated = poollt%poolpft(n) &
!!        + sum(poollt%poolpft_dgain(n,:)) &
!!        - sum(poollt%poolpft_dloss(n,:))
!
!    IF (pool_updated .LE. pool_valid) CYCLE
!
!    !...calculate loss rate
!    nref=n+npoolpft/2+1 !12,16 ntpool
!    if (pool_indx_lay(nref) .eq. 1) then
!       !krater_lay(npoolpft,nsoil), k_rate(ntpool)
!        poollt%krater_lay(n,1) = poollt%mcr_scale &
!              * poolcont%k_rate(nref)
!    else
!        do s=1,pool_indx_lay(nref)
!           poollt%krater_lay(n,s) = poollt%mrr_scale_lay(s) &
!               * poolcont%k_rate(nref)
!        enddo
!    endif
!
!    !.....calculate/check maintenance resp
!    do s=1,pool_indx_lay(nref) !12,16 ntpool, s(nref) same as for 1,5 ntpool
!       temp_mrespr = poollt%poolpft_lay(tcref,s) * & !(npoolpft,nsoil)
!                     poollt%krater_lay(tcref,s) * & !(npoolpft,nsoil)
!                     poolcont%lresp_eff(tcref) !(npoolpft)
!       temp_mrespl = temp_mrespr*dtsib
!       pool_updated = pool_updated - temp_mrespl
!       IF (pool_updated .LT. pool_valid) CYCLE
!
!       !if (poollt%rcpoolpft_lay(n,s) .gt. dzero) then
!       !temp_mresprc13 = fract%rcpoolfac*temp_mrespr
!       !else
!       !temp_mresprc13 = poollt%poolpft_lay(n,s) * & !(npoolpft,nsoil)
!       !              poollt%krater_lay(n,s) * & !(npoolpft,nsoil)
!       !              poolcont%lresp_eff(n) !(npoolpft)
!       !endif
!       !temp_mresplc13 = temp_mresprc13*dtsib
!       !pool_updated = pool_updated - temp_mresplc13
!
!       !IF (pool_updated .LT. pool_valid) CYCLE
!       poollt%loss_mresp_lay(n,s) = fract%rcpoolfac*temp_mrespr
!       !poollt%loss_mresp_lay(n,s) = temp_mresprc13 !(npoolpft,nsoil)
!       poollt%poolpft_dloss(n,s) = fract%rcpoolfac*temp_mrespl & !(npoolpft,nsoil)
!           + poollt%poolpft_dloss(n,s)
!    enddo !s=1,pool_indx_lay
!
!enddo !n=1,npoolpft
ENDIF !pools are greater than required minimum

!print*,'pool_valid lpc13: ',poollt%poolpftmin(6)

!...Additional diagnostics
poollt%autoresp_dgain = poollt%poolpft_dgain
poollt%autoresp_dloss = poollt%poolpft_dloss

!----Save respirations
poollt%resp_auto = sum(poollt%loss_gresp(1:5)) &
                   + sum(poollt%loss_mresp_lay(1:5,:))
poollt%resp_leaf = poollt%loss_gresp(lp) &
                   + sum(poollt%loss_mresp_lay(lp,:))
poollt%resp_mntn = sum(poollt%loss_mresp_lay(1:5,:))
poollt%resp_root = &
    sum(poollt%loss_mresp_lay(frp,:)) + &
    sum(poollt%loss_mresp_lay(crp,:)) + &
    poollt%loss_gresp(frp) + &
    poollt%loss_gresp(crp) 
resp_soil = poollt%resp_root
resp_soil_lay(:) = &
    poollt%loss_mresp_lay(frp,:) + &
    poollt%loss_mresp_lay(crp,:) + &
    poollt%loss_gresp(frp)*rootf_lay(:) + &
    poollt%loss_gresp(crp)*rootf_lay(:) 

!----Similar calculations for carbon-13 respiration
poollt%resp_autoc13 = sum(poollt%loss_gresp(6:10)) &
                     + sum(poollt%loss_mresp_lay(6:10,:))
poollt%resp_leafc13 = poollt%loss_gresp(lpc13) & 
                   + sum(poollt%loss_mresp_lay(lpc13,:))
poollt%resp_mntnc13 = sum(poollt%loss_mresp_lay(6:10,:))
poollt%resp_rootc13 = &
    sum(poollt%loss_mresp_lay(frpc13,:)) + &
    sum(poollt%loss_mresp_lay(crpc13,:)) + &
    poollt%loss_gresp(frpc13) + &
    poollt%loss_gresp(crpc13)
resp_soilc13 = poollt%resp_rootc13
resp_soilc13_lay(:) = &
    poollt%loss_mresp_lay(frpc13,:) + &
    poollt%loss_mresp_lay(crpc13,:) + &
    poollt%loss_gresp(frpc13)*rootf_lay(:) + &
    poollt%loss_gresp(crpc13)*rootf_lay(:)

end subroutine pool_auto_resp
