
!===================SUBROUTINE SORTIN===================================

subroutine phosort( numic, ic, gamma, range, eyy, pco2y)

!=======================================================================
!
!     ARRANGES SUCCESSIVE PCO2/ERROR PAIRS IN ORDER OF INCREASING PCO2.
!       ESTIMATES NEXT GUESS FOR PCO2 USING COMBINATION OF LINEAR AND
!       QUADRATIC FITS.
!
!=======================================================================


    use kinds

    implicit none

    !...INPUT VARIABLES
    integer(i4), intent(in) :: numic ! number of total iterations
    integer(i4), intent(in) :: ic    ! iteration count
    real(r8),intent(in) :: range     ! inititial guess of pco2y (Pa)
    real(r8),intent(in) :: gamma     ! CO2 photocompensation point (Pa)

    real(r8), dimension(numic), intent(inout) ::  &
           eyy, &  !difference between pco2y estimates (Pa)
           pco2y   !pco2y estimates (Pa)

    !..LOCAL VARIABLES
    integer ::  i,j,n,i1,i2,i3,isp,is,ix
    logical :: bitx

    real(r8) :: pmin        !
    real(r8) :: emin        !
    real(r8) :: a           !
    real(r8) :: b           !
    real(r8) :: pco2yl      !
    real(r8) :: pco2yq      !
    real(r8) :: ac1         !
    real(r8) :: ac2         !
    real(r8) :: bc1         !
    real(r8) :: bc2         !
    real(r8) :: cc1         !
    real(r8) :: cc2         !
    real(r8) :: aterm       !
    real(r8) :: bterm       !
    real(r8) :: cterm       !
    real(r8) :: pco2b       !
    real(r8) :: eyyisp      !
    real(r8) :: eyyis       !
    real(r8) :: eyyi1       !
    real(r8) :: eyyi2       !
    real(r8) :: eyyi3       !
    real(r8) :: pco2yisp    !
    real(r8) :: pco2yis     !
    real(r8) :: pco2yi1     !
    real(r8) :: pco2yi2     !
    real(r8) :: pco2yi3     !


    if( ic < 4 ) then
        pco2y(1) = gamma + 0.5_r8*range
        pco2y(2) = gamma                                             &
            + range*( 0.5_r8 - 0.3_r8*sign(done,eyy(1)) )
        pco2y(3) = pco2y(1)- ((pco2y(1)-pco2y(2))                      &
            /(eyy(1)-eyy(2)+1.e-10_r8))*eyy(1)
        pmin = min( pco2y(1), pco2y(2) )
        emin = min(   eyy(1),   eyy(2) )
        if ( emin > 0. .and. pco2y(3) > pmin )                        &
            pco2y(3) = gamma
    else

        n = ic - 1
        bitx = abs(eyy(n)) > 0.1
        if(.not. bitx) pco2y(ic) = pco2y(n)
        if(bitx) then
            do j = 2, n
                a = eyy(j)
                b = pco2y(j)
                do i = j-1,1,-1
                    if(eyy(i) <= a ) go to 100
                     ! if triggered, assigning eyy(2) as eyy(1)
                    eyy(i+1) = eyy(i)
                    pco2y(i+1) = pco2y(i)
                enddo ! i loop
                i = 0
                100        continue
                ! assigning eyy(1) as eyy(2)
                eyy(i+1) = a
                pco2y(i+1) = b
            enddo  ! j loop
        endif
        !.. should result in ordering these pairs in ascending order
!-----------------------------------------------------------------------

        if(bitx) then
            pco2b = 0.
            is    = 1
        endif

        do ix = 1, n
            if(bitx) then
                if( eyy(ix) < 0. )  then
                    pco2b = pco2y(ix)
                    is = ix
                endif
            endif
        enddo

        if(bitx) then
            i1 = is-1
            i1 = MAX(1, i1)
            i1 = min(n-2, i1)
            i2 = i1 + 1
            i3 = i1 + 2
            isp   = is + 1
            isp = min0( isp, n )
            is = isp - 1
            eyyisp = eyy(isp)
            eyyis = eyy(is)
            eyyi1 = eyy(i1)
            eyyi2 = eyy(i2)
            eyyi3 = eyy(i3)
            pco2yisp = pco2y(isp)
            pco2yis = pco2y(is)
            pco2yi1 = pco2y(i1)
            pco2yi2 = pco2y(i2)
            pco2yi3 = pco2y(i3)
        endif

        if(bitx) then

            !...Patch to check for zero in the denominator...
            if(eyyis /= eyyisp)then
                pco2yl=pco2yis - (pco2yis-pco2yisp) / (eyyis-eyyisp)*eyyis
            else
                pco2yl = pco2yis * 1.01
            endif

            !   METHOD USING A QUADRATIC FIT

            ac1 = eyyi1*eyyi1 - eyyi2*eyyi2
            ac2 = eyyi2*eyyi2 - eyyi3*eyyi3
            bc1 = eyyi1 - eyyi2
            bc2 = eyyi2 - eyyi3
            cc1 = pco2yi1 - pco2yi2
            cc2 = pco2yi2 - pco2yi3

            !...Patch to prevent zero in denominator...
            if(bc1*ac2-ac1*bc2 /= 0.0 .and. ac1 /= 0.0_r8)then
                bterm = (cc1*ac2-cc2*ac1)/(bc1*ac2-ac1*bc2)
                aterm = (cc1-bc1*bterm)/ac1
                cterm = pco2yi2-aterm*eyyi2*eyyi2-bterm*eyyi2
                pco2yq= cterm
                pco2yq= MAX( pco2yq, pco2b )
                pco2y(ic) = ( pco2yl+pco2yq)/2.0_r8
            else
                pco2y(ic) = pco2y(ic) * 1.01_r8
            endif

        endif

    endif
!
! make sure pco2 does not fall below compensation point
    pco2y(ic) = MAX(pco2y(ic),gamma+0.01_r8)

end
