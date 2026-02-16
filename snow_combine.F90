!----------------------------------------------------------------------
subroutine snow_combine(gref, glon, glat, pref, lai, &
      poros, hydrost, sscolt)

!----------------------------------------------------------------------
!
!   Based on CLM subroutine CLM_COMBIN
!
!   CLM web info: http://clm.gsfc.nasa.gov
!
!   Description:
!   This subroutine checks for elements which are below the prescribed 
!   minimum for thickness or mass.  If the snow element thickness or 
!   mass is less than a prescribed minimum, then it is combined with a 
!   neighboring element.  The subroutine subdivide_snow then executes 
!   the combination of mass and energy.
!
!   Revision History:
!   15 September 1999: Yongjiu Dai; initial code
!   15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!   30 January   2002: Ian Baker; SiB integration
!
!----------------------------------------------------------------------

use kinds

use module_pparams, only: &
   denice, denh2o, month_names
use module_sibconst, only: &
    nsnow, &
    snocmbn_print, snocmbn_stop, &
    snocmbn_thresh
use module_sib, only: &
    hydros_type, &
    sscol_type
use module_time, only: &
    day, month, year, sec_tot

implicit none

!----------------------------------------------------------------------
!...input variables
integer(i4), intent(in) :: gref, pref
real(r4), intent(in) :: glon, glat
real(r8), intent(in) :: lai, poros
type(hydros_type), intent(inout) :: hydrost
type(sscol_type), intent(inout) :: sscolt

!...local variables
integer(byte) :: nsnow_old
real(r8) :: depth_old, mass_old
real(r8) :: capacg_old, roffo_old
real(r8), dimension(-nsnow+1:1) :: &
    td_old, ice_old, liq_old
real(r8) :: zwice  ! snow-total ice (kg m^-2)
real(r8) :: zwliq  ! snow-total liquid (kg m^-2)
real(r8) :: dzmin(5)  ! minimum depth for snow layers (m)

real(r8) :: waterin, waterout, wbal
logical :: wbal_err

!...misc variables
integer(byte) :: ii, jj, j, k, m
integer(byte) :: jold, jnew
integer(byte) :: mssi, neibor

data dzmin/0.010,0.015,0.025,0.055,0.115/
!------------------------------------------------------------------

    !...set local variables
    nsnow_old = sscolt%nsl
    depth_old = hydrost%snow_gdepth
    mass_old  = hydrost%snow_gmass
    capacg_old = hydrost%capacg
    roffo_old = hydrost%roffo

    td_old = sscolt%td(-nsnow+1:1)
    ice_old = sscolt%www_ice(-nsnow+1:1)
    liq_old = sscolt%www_liq(-nsnow+1:1)
 
    !-------------------------
    do jold = nsnow_old+bone, bzero

         !...Threshold amount to maintain a 
         !.......snow layer has been 0.1 (kg/m^2 water)
        if (sscolt%www_ice(jold) <= 0.05 ) then

            sscolt%nsl = sscolt%nsl + bone
            jnew = jold + bone

            if (jnew .eq. bone) then
                !...puddle on surface if snow gone
                hydrost%capacg = hydrost%capacg + &
                      sscolt%www_liq(jold) + sscolt%www_ice(jold)
            else
                !...combine layers
                if ((sscolt%www_ice(jnew)/(sscolt%dz(jnew)*denice)) + &
                    (sscolt%www_liq(jnew)/(sscolt%dz(jnew)*denh2o)) > poros) then
                    !.....prevent supersaturating the soil column
                    hydrost%roffo = hydrost%roffo + &
                         (sscolt%www_ice(jold) + sscolt%www_liq(jold))
                else
                    sscolt%www_liq(jnew) = sscolt%www_liq(jnew) &
                                            + sscolt%www_liq(jold)
                    sscolt%www_ice(jnew) = sscolt%www_ice(jnew) &
                         + sscolt%www_ice(jold)
                 endif

            endif
                
            !...shift all layers above this down one
            if ((jold > nsnow_old+bone) .and. (sscolt%nsl < bzero)) then
                do jj = jold, sscolt%nsl+bone, bnegone
                    sscolt%dz(jj) = sscolt%dz(jj-bone)
                    sscolt%td(jj) = sscolt%td(jj-bone)
                    sscolt%www_liq(jj) = sscolt%www_liq(jj-bone)
                    sscolt%www_ice(jj) = sscolt%www_ice(jj-bone)
                enddo
            endif

        endif  !removed snow layer due to www_ice
   enddo  !jold

   !...calculate new snow depth
   hydrost%snow_gdepth = dzero
   hydrost%snow_gmass  = dzero
   zwice      = dzero
   zwliq      = dzero

   do j=sscolt%nsl+bone,bzero
      hydrost%snow_gmass  = hydrost%snow_gmass + sscolt%www_ice(j) + sscolt%www_liq(j)
      hydrost%snow_gdepth = hydrost%snow_gdepth + sscolt%dz(j)
      zwice      = zwice + sscolt%www_ice(j)
      zwliq      = zwliq + sscolt%www_liq(j)
   enddo

   !...check snow depth,
   !.....removing layers below minimal depth
   if (hydrost%snow_gdepth < 1.d-6) then !snow gone!
      sscolt%nsl = bzero
      hydrost%snow_gmass = dzero
      hydrost%snow_gdepth = dzero
      hydrost%capacg = hydrost%capacg + zwliq + zwice
   else  
      !...check layer depths
      !.....for two or more layers
      if (sscolt%nsl < -bone) then
          jj=sscolt%nsl
          mssi = 1
          do k = jj+bone, bzero
             if (sscolt%dz(k) < dzmin(mssi)) then

                !.....top node removed, combine with bottom neighbor
                if (k == sscolt%nsl+bone) then
                    neibor = k + bone

                !.....bottom node removed, combine with top neighbor
                elseif (k == 0) then
                    neibor = k - bone

                !.....mid node removed, combine with thinnest neighbor
                else
                    neibor = k + bone
                    if ((sscolt%dz(k-bone) + sscolt%dz(k)) < &
                        (sscolt%dz(k+bone) + sscolt%dz(k))) neibor = k-bone
                endif

                !node jold and jnew are combined, stored as jnew
                if (neibor > k) then
                    jnew = neibor
                    jold = k
                else
                    jnew = k
                    jold = neibor
                endif

                !.....combine thickness, water, and temperature
                call clm_combo(sscolt%dz(jnew), sscolt%www_liq(jnew), &
                               sscolt%www_ice(jnew), sscolt%td(jnew), &
                               sscolt%dz(jold), sscolt%www_liq(jold), &
                               sscolt%www_ice(jold), sscolt%td(jold) )

                !.....shift all layers above this down one
                if (jnew-1 > sscolt%nsl+1) then
                    do m = jnew-bone, sscolt%nsl+btwo, bnegone
                       sscolt%dz(m) = sscolt%dz(m-1)
                       sscolt%td(m) = sscolt%td(m-1)
                       sscolt%www_ice(m) = sscolt%www_ice(m-1)
                       sscolt%www_liq(m) = sscolt%www_liq(m-1)
                    enddo
                endif

                sscolt%nsl = sscolt%nsl + bone
                if (sscolt%nsl >= bone) cycle

             else   
                mssi = mssi + bone
             endif  !dz > dzmin
          enddo   !snow layer loop

      endif !two or more layers condition
   endif !snow depth check


    !...Update snow and clear old layers
    if (nsnow_old .lt. sscolt%nsl) then

        !...set old layers to zero
        ii=-nsnow+1
        jj=sscolt%nsl 
        sscolt%dz(ii:jj) = dzero
        sscolt%node_z(ii:jj)  = dzero
        sscolt%layer_z(ii:jj-1) = dzero

        sscolt%td(ii:jj) = dzero
        sscolt%www_liq(ii:jj) = dzero
        sscolt%vol_liq(ii:jj) = dzero
        sscolt%www_ice(ii:jj) = dzero
        sscolt%vol_ice(ii:jj) = dzero

        !...reset the node depth and 
         !.....the depth of layer interface
         do j=bzero,sscolt%nsl+bone,bnegone
           sscolt%node_z(j) = sscolt%layer_z(j) - 0.5 * sscolt%dz(j)
           sscolt%layer_z(j-1)  = sscolt%layer_z(j) - sscolt%dz(j)
        enddo
   endif

   !...Calculate balance info
   waterin = sum(ice_old) + sum(liq_old) &
             + capacg_old + roffo_old
   waterout = sum(sscolt%www_ice(-nsnow+1:1)) + &
              sum(sscolt%www_liq(-nsnow+1:1)) + &
              hydrost%capacg + hydrost%roffo
   wbal = waterout - waterin
 
   wbal_err = .false.
   if (abs(wbal) .gt. snocmbn_thresh) wbal_err = .true.

   !...Print info if requested
   if (((snocmbn_print) .and. (nsnow_old .lt. sscolt%nsl)) &
        .or. (wbal_err)) then
        print*,''
        print('(a)'), '-----SNOW COMBINE INFO-----'
        print('(a,g12.4)'), &
              'H2O Balance (kg/m2)= ',wbal
        print('(2(a,g12.4))'), &
              '  H2O In= ',waterin,' H2O Out= ',waterout
        print('(a,a,i3,a,i4)'), &
              'Date: ', trim(month_names(month)), day, ', ', year
        print'(a,i12)','Second=',sec_tot
        print'(a)',''
        print'(a,2f12.4,1i7,2i4)','Point (Lon/Lat/Ref/PFT): ', &
                 glon, glat, gref, pref
        print'(a,f8.4)',    'Current LAI: ',lai
        
        print('(a,2i6)'),    &
             'Snow Layers (in/out): ',nsnow_old,sscolt%nsl
        print('(2(a,g12.4))'), &
              'Snow Depth (in/out): ',depth_old, '   ',hydrost%snow_gdepth
        print('(2(a,g12.4))'), &
              'Snow Mass (in/out) : ',mass_old, '   ',hydrost%snow_gmass
        print('(2(a,g12.4))'), &
              'Capacg (in/out)    : ',capacg_old, '   ',hydrost%capacg
        print('(2(a,g12.4))'), &
              'Runoff (in/out)    ; ',roffo_old, '   ',hydrost%roffo

        print*,''
        print('(a)'), &
              '   Lay    ICE_ORIG     LIQ_ORIG       ICE_NEW     LIQ_NEW'
        do k=-nsnow+1, 1
           print('(i6,4g14.4)'), &
              k,ice_old(k),liq_old(k),sscolt%www_ice(k),sscolt%www_liq(k)
        enddo

        print*,''

        if ((wbal_err) .and. (snocmbn_stop)) then
           print*,'==>Water Balance Error, Stopping.'
           print*,''
           stop
        endif

     endif  !print info

end subroutine snow_combine



!=========================================================================
!=========================================================================
!----------------------------------------------------------------------
subroutine CLM_COMBO(dz1,liq1,ice1,temp1,dz2,liq2,ice2,temp2)

    use kinds
    use module_pparams, only: &
        tice, lfus, cpice, cpliq

    implicit none
    real(r8),intent(inout) :: dz1
    real(r8),intent(inout) :: liq1
    real(r8),intent(inout) :: ice1
    real(r8),intent(inout) :: temp1

    real(r8),intent(in) :: dz2
    real(r8),intent(in) :: liq2
    real(r8),intent(in) :: ice2
    real(r8),intent(in) :: temp2


    !...local variables...
    real(r8) :: dzc
    real(r8) :: wicec
    real(r8) :: wliqc
    real(r8) :: tc
    real(r8) :: h1
    real(r8) :: h2
    real(r8) :: hc

    !
    !   Code based on CLM subroutine CLM_COMBO, modified for use
    !   with SiB
    !
    !   CLM Web Info:  http://clm.gsfc.nasa.gov
    !
    !   Description
    !   Combines two elements, and returns dz (thickness) temperature
    !   www_liq and www_ice.
    !
    !   Revision History:
    !   15 September 1999; Yongjiu Dai, original code
    !   15 December 1999;  Paul Houser and Jon Radakovich, F90 revision
    !   01 March 2002;     Ian Baker, SiB integration

    dzc   = dz1 + dz2
    wicec = ice1 + ice2
    wliqc = liq1 + liq2

    h1    = (cpice*ice1+cpliq*liq1) &
        *(temp1-tice)+lfus*liq1
    h2    = (cpice*ice2+cpliq*liq2) &
        *(temp2-tice)+lfus*liq2
    hc    = h1 + h2

    if(hc < 0.0) then
        tc = tice + hc/(cpice*wicec+cpliq*wliqc)
    elseif(hc <= lfus*wliqc) then
        tc = tice
    else
        tc = tice + (hc - lfus*wliqc)/(cpice*wicec &
            +cpliq*wliqc)
    endif

    dz1   = dzc
    ice1  = wicec
    liq1  = wliqc
    temp1 = tc

    return
end subroutine clm_combo
