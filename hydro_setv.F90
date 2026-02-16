!=====================SUBROUTINE HYDRO_SETV=================================

subroutine hydro_setv( fc_eff, wp_eff, &
                      rootf, sscolt, hydrovt)
!=======================================================================
!     UPDATING OF VEGETATION (PFT) HYDROLOGICAL VARIABLES.  
!=======================================================================

use kinds
use module_sibconst, only: &
    nsoil, nsoiltop
use module_sib, only: &
    sscol_type, hydrov_type
use module_time, only: wt_clim

implicit none

!----------------------------------------------------------------------
!...input variables
real(r8), dimension(nsoil), intent(in) :: fc_eff, wp_eff, rootf
type(sscol_type), intent(in) :: sscolt
type(hydrov_type), intent(inout) :: hydrovt

!...misc variables
integer(i4) :: i

!----------------------------------------------------------------------

     !---------------------------------
     !...Available Water (PAW and TAW)
     hydrovt%pawfrw  = dzero
     hydrovt%pawftop = dzero
     hydrovt%pawfzw = dzero
     hydrovt%tawfrw  = dzero
     hydrovt%tawftop = dzero
     hydrovt%tawfzw = dzero

     do i=1,nsoil
        hydrovt%paw_lay(i) = max(dzero, &
              sscolt%vol_liq(i) - wp_eff(i))
        hydrovt%pawmax_lay(i) = fc_eff(i) - wp_eff(i)
        if (hydrovt%pawmax_lay(i) .le. dzero) then 
           print*,'Error with PAW in hydro_pft.F90'
           print*,'Stopping.'
           STOP
        endif

        hydrovt%pawfrac_lay(i) = MIN(done, &
              hydrovt%paw_lay(i) / hydrovt%pawmax_lay(i))
        hydrovt%pawfrw = hydrovt%pawfrw + &
               hydrovt%pawfrac_lay(i) * rootf(i)
        hydrovt%pawfzw = hydrovt%pawfzw + &
             hydrovt%pawfrac_lay(i) * &
             (sscolt%dz(i)/sscolt%layer_z(nsoil))
        
        hydrovt%taw_lay(i) = max(dzero, &
             (sscolt%vol_liq(i) + sscolt%vol_ice(i)) - wp_eff(i))
        hydrovt%tawfrac_lay(i) = MIN(done, &
               hydrovt%taw_lay(i) / hydrovt%pawmax_lay(i))
        hydrovt%tawfrw = hydrovt%tawfrw + &
             hydrovt%tawfrac_lay(i) * rootf(i)
        hydrovt%tawfzw = hydrovt%tawfzw + &
             hydrovt%tawfrac_lay(i) * &
             (sscolt%dz(i)/sscolt%layer_z(nsoil))

      enddo  !i=1,nsoil
  
      hydrovt%pawftop = sum(hydrovt%pawfrac_lay(1:nsoiltop))/real(nsoiltop)
      hydrovt%tawftop = sum(hydrovt%tawfrac_lay(1:nsoiltop))/real(nsoiltop)

      !----------------------------------------
      !...Long Term / Climatological PAW Updates
      hydrovt%clim_pawfrw = (1.-wt_clim)*hydrovt%clim_pawfrw &
                          + wt_clim*hydrovt%pawfrw
      hydrovt%clim_tawfrw = (1.-wt_clim)*hydrovt%clim_tawfrw &
                          + wt_clim*hydrovt%tawfrw

end subroutine hydro_setv
