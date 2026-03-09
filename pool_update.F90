!!Subroutine to update the carbon pools
!!
!!===========================================
subroutine pool_update( &
    sibpt, siblon, siblat, pref, &
    assimin, c13assimin, c14assimin, & 
    laiin, fparin, &
    equibdt, pooldt, equiblt, poollt, grz_transfer, &
    co2t, fract)
!!==========================================

use kinds
use module_sib, only: &
   equibd_type, poold_type, &
   equibl_type, pooll_type, &
   co2_type, fract_type
use module_sibconst, only: &
   npoollu, npoolpft

implicit none

!...input values
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: siblon, siblat
real(r8), intent(in) :: assimin, laiin, fparin, c13assimin, c14assimin
real(r4), dimension(npoollu+2), intent(in) :: grz_transfer

type(equibd_type), intent(inout) :: equibdt
type(equibl_type), intent(inout) :: equiblt
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(co2_type), intent(in) :: co2t
type(fract_type), intent(in) :: fract

!...local variables
integer(i4) :: p,tcref
real(r8) :: tmptotc,tmpc13,tmpdiff,tmpc14
real(r8) :: dnrzero=1.D-10
!------------------------------------------------
!Update pools and equilibrium total gains/losses
!Check carbon balance
!Reset variables


!...dead pools

do p=1,npoollu
  pooldt%poollu_lay(p,:) = pooldt%poollu_lay(p,:) &
     + pooldt%poollu_dgain(p,:) - pooldt%poollu_dloss(p,:)
enddo

do p=1,npoollu
   equibdt%poollu_totgain(p) = equibdt%poollu_totgain(p) &
       + sum(pooldt%poollu_dgain(p,:))
   equibdt%poollu_totloss(p) = equibdt%poollu_totloss(p) &
       + sum(pooldt%poollu_dloss(p,:))

   pooldt%poollu(p) = sum(pooldt%poollu_lay(p,:))
   equibdt%poollu_max(p) = MAX(equibdt%poollu_max(p), &
       pooldt%poollu(p))
   equibdt%poollu_min(p) = MIN(equibdt%poollu_min(p), &
       pooldt%poollu(p))
enddo

!...live pools
do p=1,npoolpft
  poollt%poolpft_lay(p,:) = poollt%poolpft_lay(p,:) &
   + poollt%poolpft_dgain(p,:) - poollt%poolpft_dloss(p,:)
enddo

do p=1,npoolpft
   equiblt%poolpft_totgain(p) = equiblt%poolpft_totgain(p) &
       + sum(poollt%poolpft_dgain(p,:))
   equiblt%poolpft_totloss(p) = equiblt%poolpft_totloss(p) &
       + sum(poollt%poolpft_dloss(p,:))

   poollt%poolpft(p) = sum(poollt%poolpft_lay(p,:))
   equiblt%poolpft_max(p) = MAX(equiblt%poolpft_max(p), &
       poollt%poolpft(p))
   equiblt%poolpft_min(p) = MIN(equiblt%poolpft_min(p), &
       poollt%poolpft(p))

!if (sum(poollt%rcpoolpft_lay(p,:)) .gt. 1.) then
!if ( (sum(poollt%poolpft_dloss(p,:)) .gt. 10) .or. &
!     (sum(poollt%poolpft_dloss(p,:)) .lt. -10)) then
!    print*,' '
!    print*,'code: pool_update'
!    print*,'p: ',p
!    print*,'poolpft_dloss(p,p+1):',&
!        sum(poollt%poolpft_dloss(p,:)),sum(poollt%poolpft_dloss(p+1,:))
!    print*,'poolpft(p,p+1) :',&
!        sum(poollt%poolpft_lay(p,:)),sum(poollt%poolpft_lay(p+1,:))
!    print*,' '
!endif

enddo

!...check carbon balance
call balan_carbon(sibpt, siblon, siblat, pref, &
     assimin, c13assimin, c14assimin, laiin, fparin, &
     pooldt, poollt, grz_transfer, co2t, fract)

!...reset
pooldt%poollu_dgain = dzero
pooldt%poollu_dloss = dzero
poollt%poolpft_dgain = dzero
poollt%poolpft_dloss = dzero


end subroutine pool_update
