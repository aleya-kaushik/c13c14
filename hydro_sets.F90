!
!================SUBROUTINE HYDRO_SETS=================================
!At the beginning of every time-step, update hydrological diagnostics
!  set from soil moisture.
!====================================================================== 
!
!References
!      Community Land Model (CLM)
!      Farouki (1981), Jordan (1991), de Vires (1963)
!      Sellers (1995), Sellers (1986)
!
!======================================================================

subroutine hydro_sets( &
    soilt, hydrost, sscolt)

    use kinds
    use module_pparams, only: &
        cpice, cpliq, grav, &
        tkice, tkair, tice, & 
        h2ohc, denice, denh2o
    use module_sibconst, only: &
        nsoil, nsnow
    use module_sib, only: &
        soil_type, hydros_type, &
        sscol_type
 
    implicit none

    !...Input Variables
    type(soil_type), intent(in) :: soilt
    type(hydros_type), intent(inout) :: hydrost
    type(sscol_type), intent(inout) :: sscolt

    !...Local Variables
    real(r8) :: satw   !fraction of saturation
    real(r8) :: fl     !fraction of liquid  
    real(r8) :: dke    !Kersten number
    real(r8) :: bw     !partial density of water (ice+liquid)
    real(r8) :: facorig  !original factor used for soil resistance
    real(r8) :: facnew   !fraction of saturation of top soil layer
    real(r8) :: fac      !combined factor for soil resistance
    real(r8) :: psit   !moisture potential of top soil layer
    real(r8) :: argg   !RH of air at soil surface
    real(r8) :: wbar   !column total fraction of saturation
    real(r8) :: w1     !layer one fraction of saturation
    real(r8), dimension(nsoil) :: &
        dksat  ! soil thermal conductivity, saturated soil (W/m/K)
    real(r8), dimension(-nsnow+1:nsoil) :: &
        thk    ! soil/snow thermal conductivity of layer (W/m/K)
    
    !...Misc Variables
    integer(i4) :: j

    !-------------------------------------------
    !...reset runoff
    hydrost%roff = dzero
    hydrost%roffo = dzero

    !...reset soil column saturation fraction
    hydrost%satfrac = dzero

    !-----------------------------
    !...SURFACE VARIABLES
    !...soil surface layer resistance
    facorig = min(done, sscolt%www_liq(1) / &
             (soilt%poros*sscolt%dz(1)*denh2o) + &
             sscolt%www_ice(1) / &
             (soilt%poros*sscolt%dz(1)*denice))
    facorig = max( facorig, 0.02*done)

    wbar = 0.
    do j=1,nsoil
       wbar = wbar + (sscolt%www_liq(j) &
              / (soilt%poros * sscolt%dz(j) * denh2o)) &
              * (sscolt%dz(j) / sscolt%layer_z(nsoil))
    enddo
    w1 = (sscolt%www_liq(1) &
          / (soilt%poros * sscolt%dz(1) * denh2o)) 
    facnew = min(done, max(dzero, &
             w1 / max(0.001, wbar) ))
    fac = max(0.02*done, min(facnew, facorig))
    hydrost%rsoil = exp(8.206 - 4.255 * fac)

    !...soil surface layer relative humidity
    psit = soilt%phsat * fac ** (-soilt%bee)
    argg = max(-10.*done,(psit*grav/(461.5*sscolt%td(sscolt%nsl+1))))
    hydrost%rhsoil = exp(argg)

    !-------------------
    !...SNOW VARIABLES
    if (sscolt%nsl < 0) then
        do j = sscolt%nsl+1,0
            !...thermal conductivity of snow
            bw = (sscolt%www_liq(j) + sscolt%www_ice(j))/sscolt%dz(j)
            thk(j) = tkair + (7.75E-5*bw + 1.105E-6*bw*bw)*(tkice-tkair) 

            !.....heat capacity, snow
            sscolt%shcap(j) = cpliq*sscolt%www_liq(j) + &
                               cpice*sscolt%www_ice(j)
        enddo
    endif

    !------------------
    !...SOIL VARIABLES
    do j=1,nsoil

        !...fractional expression of saturation
        satw = ((sscolt%www_liq(j)/denh2o) + (sscolt%www_ice(j)/denice))/ &
                (sscolt%dz(j) * soilt%poros)
        satw = min(done,satw)

        !...thermal conductivity of soil
        if (satw > 1.0E-6) then   ! water present in soil
            fl = sscolt%www_liq(j) / &
                  (sscolt%www_liq(j) + sscolt%www_ice(j)) ! frac liq

            if (sscolt%td(j) >= tice) then ! unfrozen soil
                ! kersten number
                dke = max(dzero, log10(satw) + done)
                dksat(j) = soilt%tksat
            else ! frozen soil
                dke = satw
                dksat(j) = soilt%tkmg * 0.249**(fl*soilt%poros) &
                    *2.29**soilt%poros
            endif

            thk(j) = dke*dksat(j) + (1.0-dke)*soilt%tkdry

        else ! soil is very dry
            thk(j) = soilt%tkdry
        endif

        !...heat capacity, soil
        sscolt%shcap(j) = soilt%csolid*(1.0-soilt%poros)*sscolt%dz(j) &
            + sscolt%www_ice(j)*cpice + sscolt%www_liq(j)*cpliq

        !...volume liq/ice and effective porosity
        sscolt%vol_ice(j) = min(soilt%poros, &
                sscolt%www_ice(j)/(sscolt%dz(j)*denice))
        sscolt%eff_poros(j) = soilt%poros - sscolt%vol_ice(j)
        sscolt%vol_liq(j) = min(sscolt%eff_poros(j), &
                sscolt%www_liq(j)/(sscolt%dz(j)*denh2o))

       !...saturation capacity
       sscolt%satfrac_lay(j) = &
            ((sscolt%www_liq(j)/denh2o) + (sscolt%www_ice(j)/denice)) / &
             (sscolt%dz(j) * soilt%poros)
       hydrost%satfrac = hydrost%satfrac + sscolt%satfrac_lay(j)

    enddo  !nsoil loop
    hydrost%satfrac = hydrost%satfrac/nsoil
    sscolt%shcap(1) = sscolt%shcap(1) + &
                      (hydrost%capacg/denh2o) * h2ohc

    !------------------------------
    !...SNOW/SOIL COLUMN VARIABLES
    do j = sscolt%nsl+1,nsoil-1
        !...snow/soil thermal conductivity
        sscolt%tksoil(j) = thk(j)*thk(j+1)* &
            (sscolt%node_z(j+1)-sscolt%node_z(j))  &
            /(thk(j)*(sscolt%node_z(j+1)-sscolt%layer_z(j)) + &
            thk(j+1)*(sscolt%layer_z(j)-sscolt%node_z(j)))

         !...snow/soil heat flux term
         sscolt%slamda(j) = sscolt%tksoil(j) / &
             (sscolt%node_z(j+1) - sscolt%node_z(j))
    enddo
    sscolt%tksoil(nsoil) = sscolt%tksoil(nsoil-1)
    sscolt%slamda(nsoil) = sscolt%slamda(nsoil-1)

end subroutine hydro_sets
