!======================================================================
subroutine pool_het_retr(poolcont, &
    zm, woptzm, wsat, seas_precip, clim_precip, &
    clim_assim, assimd, rootf_lay, pawfrac_lay, &
    td_lay, satfrac_lay, pooldt, fract)
!======================================================================

! Description
! ------------
! Calculates the heterotrophic respiration from the dead
!   surface and soil pools.
! 
! Currently the vertical distribution of carbon follows
!   the rooting profile (lrootf)
!
use kinds
use module_oparams, only: &
     rt_moist_exp, rt_moist_range
use module_param, only: pool_param
use module_poolinfo, only: &
     pool_indx_lay, pool_indx_sfc, &
     pool_indx_soil, &
     pool_indx_sfcc13, pool_indx_soilc13
use module_sib, only: &
     poold_type, fract_type
use module_sibconst, only: &
     nsoil, npoolpft, npoollu, &
     npoolsfc, npoolsoil, &
     npoolsfcc13, npoolsoilc13
use module_time, only: dtsib

implicit none


!...input variables
type(pool_param), intent(in) :: poolcont
real(r8), intent(in) :: zm, woptzm, wsat
real(r8), intent(in) :: seas_precip, clim_precip
real(r8), intent(in) :: clim_assim, assimd
real(r8), dimension(nsoil), intent(in) :: rootf_lay, pawfrac_lay
real(r8), dimension(nsoil), intent(in) :: td_lay
real(r8), dimension(nsoil), intent(in) :: satfrac_lay
type(poold_type), intent(inout) :: pooldt
type(fract_type), intent(in) :: fract

!...local variables
real(r8) :: wet_exp  !soil wetness exponent (-)
                     !(b in eqn 8 Denning et al [1996])
real(r8) :: &
    temp_loss_lay(npoollu,nsoil)   !pool loss from het resp (mole/m2/s)
real(r8) :: slope, yint

!...misc variables
integer(i4) :: j,jref,k,kref,n,nref,s,tcref
real(r8) :: qt, mtemp

!-----------------------------------------
!...Reset Variables
pooldt%mhrt_sfc_assim = done
pooldt%mhrt_sfc_freeze = done
pooldt%mhrt_sfc_hot = done
pooldt%mhrt_sfc_precip = done
pooldt%mhrt_sfc_scale = done
pooldt%mhrt_soil_freeze_lay = done
pooldt%mhrt_soil_hot_lay = done
pooldt%mhrt_soil_moist_lay = done
pooldt%mhrt_soil_pawf_lay = done
pooldt%mhrt_soil_scale_lay = done
pooldt%mhrt_soil_assim = dzero
pooldt%mhrt_soil_freeze = dzero
pooldt%mhrt_soil_hot = dzero
pooldt%mhrt_soil_moist = dzero
pooldt%mhrt_soil_pawfrw = dzero
pooldt%mhrt_soil_scale = dzero
pooldt%kratert_lay = dzero
pooldt%gain_transd_lay = dzero
pooldt%loss_resp_lay = dzero
pooldt%loss_trans_lay = dzero
pooldt%resp_het = dzero
pooldt%resp_soilnr = dzero
pooldt%resp_soilnr_lay = dzero
pooldt%resp_hetc13 = dzero
pooldt%resp_soilnrc13 = dzero
pooldt%resp_soilnrc13_lay = dzero

!...Set Local Variables
temp_loss_lay(:,:) = dzero

!...Only respire/transfer if pools exist
IF (sum(pooldt%poollu(1:6)) .gt. dzero) THEN


!-----Heterotrophic Respiration Scaling Factors-----
!...Surface scaling factors

!Assimilation Rate Scalar
if (assimd .lt. clim_assim*poolcont%hrt_sfc_aml) then
   pooldt%mhrt_sfc_assim = poolcont%hrt_sfc_amin
elseif (assimd .lt. clim_assim*poolcont%hrt_sfc_amh) then
   slope = (poolcont%hrt_sfc_amax - poolcont%hrt_sfc_amin) &
       /(clim_assim*poolcont%hrt_sfc_amh - clim_assim*poolcont%hrt_sfc_aml)
   yint = poolcont%hrt_sfc_amin - slope*clim_assim*poolcont%hrt_sfc_aml
   pooldt%mhrt_sfc_assim = assimd*slope + yint
else
   pooldt%mhrt_sfc_assim = poolcont%hrt_sfc_amax
endif
   
!...Freeze Inhibition
pooldt%mhrt_sfc_freeze = MAX(poolcont%hrt_sfc_fmin, MIN(1.0, &
           EXP(poolcont%hrt_sfc_fmul*(td_lay(1) - poolcont%hrt_sfc_fref))))

!...High Temperature Exponential
!...Q10 with increasing temperature
pooldt%mhrt_sfc_hot = MIN(poolcont%hrt_sfc_hmax, &
    MAX(done, poolcont%hrt_sfc_hq10** &
        ((td_lay(1)-poolcont%hrt_sfc_href)/10.0)))

!...Precipitation Inhbition
pooldt%mhrt_sfc_precip = MIN(done, MAX(poolcont%hrt_sfc_pmin, &
    seas_precip / (clim_precip*poolcont%hrt_sfc_pml)))

!surface combined scaling factor
pooldt%mhrt_sfc_scale = pooldt%mhrt_sfc_freeze &
     * pooldt%mhrt_sfc_hot * pooldt%mhrt_sfc_precip &
     * pooldt%mhrt_sfc_assim

!...Soil scaling factors non-depth related
!Assimilation Rate Scalar
if (assimd .lt. clim_assim*poolcont%hrt_soil_aml) then
   pooldt%mhrt_soil_assim = poolcont%hrt_soil_amin
elseif (assimd .lt. clim_assim*poolcont%hrt_soil_amh) then
   slope = (poolcont%hrt_soil_amax - poolcont%hrt_soil_amin) &
       /(clim_assim*poolcont%hrt_soil_amh - clim_assim*poolcont%hrt_soil_aml)
   yint = poolcont%hrt_soil_amin - slope*clim_assim*poolcont%hrt_soil_aml
   pooldt%mhrt_soil_assim = assimd*slope + yint
else
   pooldt%mhrt_soil_assim = poolcont%hrt_soil_amax
endif

!...Soil scaling factors per soil layer
do s=1,nsoil
   !...Soil Freeze Inhibition
   pooldt%mhrt_soil_freeze_lay(s) = &
      MAX(poolcont%hrt_soil_fmin, MIN(1.0, &
          EXP(poolcont%hrt_soil_fmul &
             *(td_lay(s) - poolcont%hrt_soil_fref))))
   pooldt%mhrt_soil_freeze = pooldt%mhrt_soil_freeze &
       + pooldt%mhrt_soil_freeze_lay(s) * rootf_lay(s)

   !...Soil High Temp Exponential
    qt = 0.1 * (td_lay(s) - poolcont%hrt_soil_href)
    mtemp = poolcont%hrt_soil_hq10**qt
    pooldt%mhrt_soil_hot_lay(s) = MIN(poolcont%hrt_soil_hmax, &
       MAX(done, mtemp))
    pooldt%mhrt_soil_hot = pooldt%mhrt_soil_hot &
       + pooldt%mhrt_soil_hot_lay(s) * rootf_lay(s)

   !...Soil Moisture Factor
   !.....Same formulation as het resp moisture factor
   !.....(Li et al., 2013; Flexas et al., 2005/2006)
   wet_exp = rt_moist_exp
   if (satfrac_lay(s) > dzero) then
       wet_exp = ((satfrac_lay(s)**zm-woptzm) / &
                  (1.-woptzm))**2
   endif
   wet_exp = min(wet_exp,rt_moist_exp)
   pooldt%mhrt_soil_moist_lay(s) = MIN(done, &
        MAX(poolcont%hrt_soil_mmin, &
            rt_moist_range * wsat**wet_exp))
   pooldt%mhrt_soil_moist = pooldt%mhrt_soil_moist  &
        + pooldt%mhrt_soil_moist_lay(s) * rootf_lay(s)

   !.....Plant Available Water (PAW)
   pooldt%mhrt_soil_pawf_lay(s) = MIN(done, &
         MAX(poolcont%hrt_soil_pawmin, pawfrac_lay(s)))
   pooldt%mhrt_soil_pawfrw = pooldt%mhrt_soil_pawfrw &
         + pooldt%mhrt_soil_pawf_lay(s) * rootf_lay(s)
enddo  !s=1,nsoil

!Combined
pooldt%mhrt_soil_scale_lay(:) = pooldt%mhrt_soil_freeze_lay(:) &
   * pooldt%mhrt_soil_hot_lay(:) * pooldt%mhrt_soil_moist_lay(:) &
   * pooldt%mhrt_soil_pawf_lay(:)
DO s=1,nsoil
   pooldt%mhrt_soil_scale_lay(s) = pooldt%mhrt_soil_scale_lay(s) &
        * pooldt%mhrt_soil_assim
   pooldt%mhrt_soil_scale = pooldt%mhrt_soil_scale  &
       + pooldt%mhrt_soil_scale_lay(s) * rootf_lay(s)
ENDDO

!-----Scaled Surface/Soil Pool Decay Rate Constants-----
!...surface pools
do n=1,npoolsfc
    nref = pool_indx_sfc(n) !6-8
    do s=1,pool_indx_lay(nref)
        !pooldt%kratert_lay(npoollu,nsoil)
        pooldt%kratert_lay(nref-npoolpft/2,s) = & ! -5 (npoolpft=10), indexes 1-3
            pooldt%mhrt_sfc_scale*poolcont%k_rate(nref)
        !k_rate is indexed as (ntpool)
     enddo
enddo

!...soil pools
do n=1,npoolsoil
    nref=pool_indx_soil(n) ! 2,3 and 9-11
    if (nref .gt. npoolpft/2) then ! gt 5
        do s=1,pool_indx_lay(nref)
            pooldt%kratert_lay(nref-npoolpft/2,s) = & !4-6 dead pools
                pooldt%mhrt_soil_scale_lay(s)*poolcont%k_rate(nref)
        enddo
    endif
enddo

!...surface C13 pools
do n=1,npoolsfcc13
    nref = pool_indx_sfcc13(n) !17-19
    do s=1,pool_indx_lay(nref)
        pooldt%kratert_lay(nref-npoolpft,s) = & ! -10, indexes 7-9 dead pools
            pooldt%mhrt_sfc_scale*poolcont%k_rate(nref)
     enddo
enddo

!...soil C-13 pools
do n=1,npoolsoilc13
    nref=pool_indx_soilc13(n) !13,14 and 20-22 
    if (nref .gt. (npoolpft+6)) then !gt 16
        do s=1,pool_indx_lay(nref)
            pooldt%kratert_lay(nref-npoolpft,s) = & !10-12 dead pools
                pooldt%mhrt_soil_scale_lay(s)*poolcont%k_rate(nref)
        enddo
    endif
enddo

!----Surface/Soil Losses from Heterotrophic Respiration-----
do n=1,npoollu/2 !npoollu=12 dead pools, 1-6 total C dead pools
    do s=1,pool_indx_lay(npoolpft/2+n) !npoolpft=10, index 5+n: 6-11 ntpool
        temp_loss_lay(n,s) = pooldt%poollu_lay(n,s) &
                           * pooldt%kratert_lay(n,s)
    enddo
enddo

!...loop through C-13 pools
do n=npoollu/2+1,npoollu !12 dead pools, 7-12 C-13 dead pools
    tcref=n-6 !totC pools (1,6)
    do s=1,pool_indx_lay(npoolpft+n) !index 10+n, pool_indx_lay(17-22) ntpool
!        temp_loss_lay(n,s) = pooldt%poollu_lay(n,s) & !(npoollu,nsoil)
!                           * pooldt%kratert_lay(n,s) !(npoollu,nsoil)
        !temp_loss_lay(n,s) = fract%rcpoolfac * temp_loss_lay(tcref,s)
        temp_loss_lay(n,s) = pooldt%rcpoollu_lay(n,s) * temp_loss_lay(tcref,s)
    enddo
enddo

!-----Respiration Fluxes and Carbon Transfers-----
!...j is the sending/from pool
!...k is the recieving/to pool
do j=1,npoollu/2 !npoollu=12, (1,6) dead pools
   jref=j+npoolpft/2 ! npoolpft=10, (6,11) from ntpool 

    do k=1,npoollu/2 ! (1,6) dead pools
       kref=k+npoolpft/2 ! (6,11) ntpool
       if (poolcont%pool_trans_frac(jref,kref) > dzero) then
          !.....transfer/respiration losses
         do s=1,pool_indx_lay(jref)
             pooldt%loss_resp_lay(j,s) = &
                 pooldt%loss_resp_lay(j,s) &
                 + temp_loss_lay(j,s) &
                 * poolcont%dresp_eff(j,k) &
                 * poolcont%pool_trans_frac(jref,kref)

             pooldt%loss_trans_lay(j,s) = &
                 pooldt%loss_trans_lay(j,s) &
                  + temp_loss_lay(j,s) &
                  * (1. - poolcont%dresp_eff(j,k)) &
                  * poolcont%pool_trans_frac(jref,kref)
          enddo

         !.....transfer gains
         if ((pool_indx_lay(jref) .eq. 1) .and. &
             (pool_indx_lay(kref) .eq. 1)) then
             pooldt%gain_transd_lay(k,1) = & 
                pooldt%gain_transd_lay(k,1) &
                + temp_loss_lay(j,1)  &
                 * (1. - poolcont%dresp_eff(j,k)) &
                 * poolcont%pool_trans_frac(jref,kref)
         elseif ((pool_indx_lay(jref) .eq. 1) .and. &
                 (pool_indx_lay(kref) .eq. nsoil)) then
                do s=1,nsoil
                   pooldt%gain_transd_lay(k,s) =  &
                      pooldt%gain_transd_lay(k,s) &
                      + temp_loss_lay(j,1) * rootf_lay(s)  &
                      * (1. - poolcont%dresp_eff(j,k)) &
                      * poolcont%pool_trans_frac(jref,kref)
                enddo
         elseif ((pool_indx_lay(jref) .eq. nsoil) .and. &
                 (pool_indx_lay(kref) .eq. nsoil)) then
                 pooldt%gain_transd_lay(k,:) = &
                      pooldt%gain_transd_lay(k,:) &
                      + temp_loss_lay(j,:) &
                      * (1. - poolcont%dresp_eff(j,k)) &
                      * poolcont%pool_trans_frac(jref,kref)
          else
                print*,'Mismatching levels between pool transfers.'
                print*,'Stopping in pool_resp_het.'
                stop
          endif !transfer gains
        endif !trans frac > 0
    enddo  !k=1,npoollu
enddo  !j=1,npoollu

!...same calculations as above but for C-13 pools
do j=npoollu/2+1,npoollu !(7,12) npoollu dead pools
   jref=j+npoolpft !(17,22) ntpool

    do k=npoollu/2+1,npoollu !7,12 npoollu dead pools
       kref=k+npoolpft !(17,22) ntpool
       if (poolcont%pool_trans_frac(jref,kref) > dzero) then
          !.....transfer/respiration losses
         do s=1,pool_indx_lay(jref)
             pooldt%loss_resp_lay(j,s) = & !loss_resp_lay(npoollu,nsoil)
                 pooldt%loss_resp_lay(j,s) &
                 + temp_loss_lay(j,s) & !temp_loss_lay(npoollu,nsoil)
                 * poolcont%dresp_eff(j,k) & !dresp_eff(npoollu,npoollu)
                 * poolcont%pool_trans_frac(jref,kref) !*frac(ntpool,ntpool)

             pooldt%loss_trans_lay(j,s) = & !npoollu,nsoil
                 pooldt%loss_trans_lay(j,s) &
                  + temp_loss_lay(j,s) &
                  * (1. - poolcont%dresp_eff(j,k)) &
                  * poolcont%pool_trans_frac(jref,kref)
          enddo

         !.....transfer gains
         if ((pool_indx_lay(jref) .eq. 1) .and. &
             (pool_indx_lay(kref) .eq. 1)) then
             pooldt%gain_transd_lay(k,1) = & !npoollu,nsoil
                pooldt%gain_transd_lay(k,1) &
                + temp_loss_lay(j,1)  &
                 * (1. - poolcont%dresp_eff(j,k)) &
                 * poolcont%pool_trans_frac(jref,kref)
         elseif ((pool_indx_lay(jref) .eq. 1) .and. &
                 (pool_indx_lay(kref) .eq. nsoil)) then
                do s=1,nsoil
                   pooldt%gain_transd_lay(k,s) =  &
                      pooldt%gain_transd_lay(k,s) &
                      + temp_loss_lay(j,1) * rootf_lay(s)  &
                      * (1. - poolcont%dresp_eff(j,k)) &
                      * poolcont%pool_trans_frac(jref,kref)
                enddo
         elseif ((pool_indx_lay(jref) .eq. nsoil) .and. &
                 (pool_indx_lay(kref) .eq. nsoil)) then
                 pooldt%gain_transd_lay(k,:) = &
                      pooldt%gain_transd_lay(k,:) &
                      + temp_loss_lay(j,:) &
                      * (1. - poolcont%dresp_eff(j,k)) &
                      * poolcont%pool_trans_frac(jref,kref)
          else
                print*,'Mismatching levels between C13 pool transfers.'
                print*,'Stopping in c13 pool_resp_het.'
                stop
          endif !transfer gains
        endif !trans frac > 0
    enddo  !k=1,npoollu
enddo  !j=1,npoollu
ENDIF !respire/transfer only if pools > dzero

!-----Accumulate/Combine Fluxes-----
!...Total heterotrophic respiration
do j=1,npoollu/2 !(1,6)
   pooldt%resp_het = pooldt%resp_het  &
          + sum(pooldt%loss_resp_lay(j,:))
   IF (pool_indx_lay(j+npoolpft/2) .gt. 1) THEN !6,11
      pooldt%resp_soil = pooldt%resp_soil &
          + sum(pooldt%loss_resp_lay(j,:))
      pooldt%resp_soil_lay(:) = pooldt%resp_soil_lay(:) &
          + pooldt%loss_resp_lay(j,:)
      pooldt%resp_soilnr = pooldt%resp_soilnr &
          + sum(pooldt%loss_resp_lay(j,:))
      pooldt%resp_soilnr_lay(:) = pooldt%resp_soilnr_lay(:) &
          + pooldt%loss_resp_lay(j,:)
   ENDIF
enddo

!...same calculations but for C-13 pools
do j=npoollu/2+1,npoollu !(7,12)
   pooldt%resp_hetc13 = pooldt%resp_hetc13  &
          + sum(pooldt%loss_resp_lay(j,:))
   IF (pool_indx_lay(j+npoolpft) .gt. 1) THEN !17,22
      pooldt%resp_soilc13 = pooldt%resp_soilc13 &
          + sum(pooldt%loss_resp_lay(j,:))
      pooldt%resp_soilc13_lay(:) = pooldt%resp_soilc13_lay(:) &
          + pooldt%loss_resp_lay(j,:)
      pooldt%resp_soilnrc13 = pooldt%resp_soilnrc13 &
          + sum(pooldt%loss_resp_lay(j,:))
      pooldt%resp_soilnrc13_lay(:) = pooldt%resp_soilnrc13_lay(:) &
          + pooldt%loss_resp_lay(j,:)
   ENDIF
enddo


!...Accumulate losses for pool updates
pooldt%poollu_dgain = pooldt%poollu_dgain & !(npoollu,nsoil)
   + pooldt%gain_transd_lay*dtsib
pooldt%poollu_dloss = pooldt%poollu_dloss & !(npoollu,nsoil)
   + pooldt%loss_resp_lay*dtsib + pooldt%loss_trans_lay*dtsib


end subroutine pool_het_retr
