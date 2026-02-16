!----------------------------------------------------------------------
subroutine sib_control ( gnum, gref,  &
       lonsib, latsib, daynew, daylmax, &
       doy, sibg )

    use kinds
    use module_param, only: &
       phencon, physcon, poolcon
    use module_pftinfo
    use module_sib, only: &
        gridcell_type
    use module_sibconst, only: nsoil

    implicit none

    !---INPUT VARIABLES
    integer(i4), intent(in) :: gnum, gref, doy
    real(r4), intent(in) :: lonsib, latsib, daylmax
    logical, intent(in) :: daynew
    type (gridcell_type), intent(inout) :: sibg

    !---LOCAL VARIABLES
    integer(i4) :: l,pref,pnum
    integer(i4) :: vegtype
    logical :: isbare, iscrop, isgrass, isveg
    
!-------------------------------------------------------------------------

! Loop over land units
do l=1,sibg%g_nlu

   if ((.not. sibg%l(l)%equibdt%lupft_spunup)) then

       pref = sibg%l(l)%ipft
       pnum = pft_num(pref)
       vegtype = pft_type(pnum)

       isbare = (vegtype .eq. type_bare)
       iscrop = (vegtype .eq. type_crop)
       isgrass = (vegtype .eq. type_grass)
       isveg = (.not. isbare)


       if (isveg) then
          !...Update phenology
          !...(potentials every time-step, stage once daily)
          call phen_update( &
               gnum, l, pnum, sibg%l(l)%ipft, &
               phencon(pnum), poolcon(pnum), &
               daynew, daylmax, sibg%gdiagt%daylen, &
               sibg%gdiagt%daylendt, &
               sibg%gprogt%clim_cupr, sibg%gprogt%clim_precip, &
               sibg%gprogt%cupr, sibg%gprogt%lspr, sibg%gdiagt%tmdf, &
               sibg%gprogt%tm, sibg%l(l)%co2t%assim, &
               sibg%l(l)%co2t%rstfac(2), sibg%l(l)%co2t%rstfac(4), &
               sibg%l(l)%vegt%lai, sibg%l(l)%hydrovt, sibg%l(l)%phent, &
               sibg%l(l)%pooldt, sibg%l(l)%poollt, sibg%l(l)%vegt%vmax, &
               sibg%l(l)%fract%c13assim, physcon(pnum), sibg%l(l)%fract)

          if (daynew) then

             !...Update pools (once daily)
             call pool_update(gref, lonsib, latsib, pref, &
                  sibg%l(l)%co2t%assim, sibg%l(l)%fract%c13assim, &
                  sibg%l(l)%vegt%lai, sibg%l(l)%vegt%fpar, &
                  sibg%l(l)%equibdt, sibg%l(l)%pooldt, &
                  sibg%l(l)%equiblt, sibg%l(l)%poollt, &
                  poolcon(pnum)%graze_trans, sibg%l(l)%co2t, &
                  sibg%l(l)%fract)

             !...Update vegetation (once daily)
             call veg_update(doy, &
                  gref, lonsib, latsib, &
                  pnum, pref, iscrop, isgrass, &
                  physcon(pnum), sibg%l(l)%pooldt%poollu, &
                  sibg%l(l)%hydrost%snow_cvfc, &
                  sibg%l(l)%sscolt%node_z(1:nsoil), &
                  sibg%l(l)%poollt, sibg%l(l)%vegt)
          endif
       endif

       !...Call SiB Main
       call sib_main(gref, lonsib, latsib, &
            pref, pnum, physcon(pnum), &
            sibg%gprogt, sibg%gdiagt, sibg%l(l))

       !...Calculate Pool Assimilation Gains
       call pool_assim( &
            gref, lonsib, latsib, pref, &
            physcon(pnum)%hhti, physcon(pnum)%hlti, &
            physcon(pnum)%shti, physcon(pnum)%slti, &
            poolcon(pnum)%gr_frac, sibg%l(l)%co2t%assim, &
            sibg%l(l)%fract%c13assim, sibg%l(l)%co2t%rstfac(2), &
            sibg%gprogt%tm, sibg%l(l)%poollt, sibg%l(l)%co2t, &
            sibg%l(l)%fract)

       !...Autotrophic Respiration
       call pool_auto_resp(poolcon(pnum), &
            sibg%l(l)%co2t%clim_assim, sibg%l(l)%vegt%clim_lai, &
            sibg%l(l)%co2t%assimd, sibg%l(l)%fract%c13assimd, & 
            sibg%l(l)%vegt%lai, &
            sibg%l(l)%cast%tc, sibg%l(l)%sscolt%td(1:nsoil), &
            sibg%l(l)%vegt%rootf, sibg%l(l)%poollt, &
            sibg%l(l)%pooldt%resp_soil, sibg%l(l)%pooldt%resp_soil_lay, &
            sibg%l(l)%pooldt%resp_soilc13, sibg%l(l)%pooldt%resp_soilc13_lay, &
            sibg%l(l)%fract)

       !...Autotrophic Transfer
       call pool_auto_tran(poolcon(pnum), &
            sibg%gdiagt%daylen, sibg%gdiagt%daylendt, &
            daylmax, sibg%l(l)%cast%tc, sibg%l(l)%hydrovt%pawfrw, &
            sibg%l(l)%vegt%rootf, sibg%l(l)%poollt, &
            sibg%l(l)%pooldt%gain_transl_lay, &
            sibg%l(l)%pooldt%poollu_dgain, sibg%l(l)%fract) 

       !...Grazing
       IF (physcon(pnum)%pftgraze) THEN
          !call pool_graze(poolcon(pnum)%poolpft_min, &
          !     poolcon(pnum)%graze_trans, &
          !     sibg%l(l)%vegt%clim_lai, sibg%l(l)%vegt%lai, &
          !     sibg%l(l)%poollt, sibg%l(l)%pooldt)
          call pool_graze(poolcon(pnum), &
               sibg%l(l)%vegt%clim_lai, sibg%l(l)%vegt%lai, &
               sibg%l(l)%poollt, sibg%l(l)%pooldt, sibg%l(l)%fract)
       ENDIF

       !...Heterotrophic Respiration/Transfer
       call pool_het_retr(poolcon(pnum), &
            sibg%l(l)%soilt%zm, sibg%l(l)%soilt%woptzm, &
            sibg%l(l)%soilt%wsat, &
            sibg%gprogt%seas_precip, sibg%gprogt%clim_precip, &
            sibg%l(l)%co2t%clim_assim, sibg%l(l)%co2t%assimd, &
            sibg%l(l)%vegt%rootf, &
            sibg%l(l)%hydrovt%pawfrac_lay(1:nsoil), &
            sibg%l(l)%sscolt%td(1:nsoil), &
            sibg%l(l)%sscolt%satfrac_lay, sibg%l(l)%pooldt, &
            sibg%l(l)%fract)   


       !...C13 isotope signature calculations
       call c13_iso_calc(poolcon(pnum),sibg%l(l)%poollt, &
            sibg%l(l)%pooldt, sibg%l(l)%fract)

       !...Save diagnostics
       call diagnostic_save_lall(gnum, l, pnum, &
            sibg%l(l)%soilt, &
            sibg%l(l)%cast, sibg%l(l)%co2t,  &
            sibg%l(l)%cost, sibg%l(l)%fract, sibg%l(l)%fluxt, &
            sibg%l(l)%hydrost, sibg%l(l)%hydrovt, &
            sibg%l(l)%phent, sibg%l(l)%pooldt, &
            sibg%l(l)%poollt, sibg%l(l)%radt, &
            sibg%l(l)%sift, sibg%l(l)%sscolt, &
            sibg%l(l)%vegt)


   endif !.not. lu_spunup

enddo !l=1,g_nlu
 
end subroutine sib_control
