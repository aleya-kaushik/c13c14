subroutine flux_vrbrd( &
    pnum, z2, nsl,     &
    rst, td, cast, vegt, &
    gprogt, gdiagt,  &
    fluxt, hydrost)

    use kinds
    use module_local
    use module_pftinfo, only: &
        pft_group, group_ndlfor
    use module_pparams, only: &
        snofac, vkrmn
    use module_sibconst, only:  &
        nsnow, nsoil
    use module_sib, only: &
        gprog_type, gdiag_type, &
        cas_type, veg_type, &
        flux_type, hydros_type
    use module_time, only: dtsib

    implicit none

    !----------------------------------------------------------------------

    !...INPUT VARIABLES
    integer(i4), intent(in) :: pnum
    real(r4), intent(in) :: z2
    integer(byte), intent(in) :: nsl
    real(r8), intent(in) :: rst, td
    type(gprog_type), intent(in) :: gprogt
    type(gdiag_type), intent(in) :: gdiagt
    type(cas_type), intent(in) :: cast
    type(veg_type), intent(in) :: vegt
    type(flux_type), intent(inout) :: fluxt
    type(hydros_type), intent(inout) :: hydrost

    !...LOCAL VARIABLES
    integer(byte) :: dirg
    real(r8) :: u2 
    real(r8) :: cuni 
    real(r8) :: temv
    real(r8) :: epct
    real(r8) :: eps, epsc, epsg
    logical  :: isnforest
    
    !----------------------------------------------------------------------  

    !...set local variables
    eps    = 1. / snofac                                   
    isnforest = (pft_group(pnum) .eq. group_ndlfor)
    
    !...ventilation mass flux
    call flux_vmf(vegt%zzwind, vegt%zztemp, vegt%z0, &
       gdiagt%ros, gprogt%spdm, gprogt%sh, gdiagt%thm, &
       cast%shcas, cast%thcas, &
       fluxt%ustar, cuni, fluxt%cu, fluxt%ct, fluxt%ventmf)

    !...aerodynamic resistance
    fluxt%ra = gdiagt%ros / fluxt%ventmf 
    fluxt%drag = gdiagt%ros * fluxt%cu * fluxt%ustar
    temv = (z2 - vegt%zp_dispd) / vegt%z0
    temv = max(temv,1.00_r8)
    temv = log(temv) 
    u2     = gprogt%spdm / (cuni * vkrmn)
    u2 = MAX(1.0_r8, u2 * temv)

    !...calculate canopy-CAS and ground-CAS resistance
    call flux_rbrd(z2, u2, &
               cast%tc, cast%tcas, td, &
               hydrost%snow_cvfc, vegt%lai, &
               vegt%cc1, vegt%cc2,   &
               fluxt%rbc, fluxt%rdc, &
               fluxt%rb, fluxt%rd )

    epsc = 1.
    epsg = 1. 
    !.....this only makes sense for canopy leaves, 
    !.....since there can only be water OR snow, not both. 
    if (hydrost%capacc_snow > 0.0) epsc = eps
    if (nsl < 0)   epsg = eps

    !...compute resistance terms
    if (cast%eacas > etg) then 
        dirg = 0 
    else 
        dirg = 1
    endif
    fluxt%rc  = rst + fluxt%rb + fluxt%rb
    fluxt%rds = hydrost%rsoil * real(dirg) + fluxt%rd

    !...compute evaporation/transpiration terms
    if (isnforest) then
       !epct = MIN(1.1, MAX(0.5, &
       !         hydrost%satfrac*2. - 0.4))
       epct = MIN(0.8, MAX(0.2, &
                 hydrost%satfrac*2. - 0.4))
    else
       epct = done
    endif

    gect = epct * (1. - hydrost%wetfracc) / fluxt%rc
    geci = epsc * hydrost%wetfracc / (fluxt%rb + fluxt%rb)
    gegs =  (1. - hydrost%wetfracg) / fluxt%rds
    gegi = epsg * hydrost%wetfracg / fluxt%rd

    coc = gect + geci

    !...calculate ecmass -- canopy evapotranspiration
    !...temporary value to be passed into phostress
    hydrost%ecmass = (etc - cast%eacas) * coc *  &
        gdiagt%ros  * 0.622e0_r8 /gprogt%ps * dtsib
                        
end subroutine flux_vrbrd
