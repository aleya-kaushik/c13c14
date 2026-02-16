module fcos_solver_Ogee

!!!!!!!!!!!!!!!!!!!
!!! Calculation of COS soil flux after Ogee et al., 2016.
!!! Added by Linda Kooijmans, July 2020.
!!!
!!!  Ogée, J., J. Sauze, J. Kesselmeier, B. Genty, H. Van Diest, T. Launois, and
!!!        L.Wingate (2016). A new mechanistic framework to predict OCS fluxes
!!!        from soils. Biogeosciences, 13(8):2221-2240. doi: 10.5194/bg-13-2221-2016.
!!!!!!!!!!!!!!!!!!!


implicit none

contains

!! Calculation COS soil Flux after Ogee et al., 2015
subroutine cos_flux_solver_Ogee( &
    coscas, T_s, poros, theta_w, pressure, P, &
    kuncat, fca, &
    Fcos)

    real(8), intent(in) :: coscas
    real(8), intent(in) :: T_s, poros, theta_w, pressure, P, kuncat, fca
    real(8), intent(out) :: Fcos

    ! intermediate variables
    real(8) :: &
         !pH = 4.5, &
         D0l25 = 1.94d-9, &
         !pKw = 14.0, &
         R = 8.3145, & ! J/(mol K)
         DHa = 40.0d3, & ! J/mol
         DHd = 200.0d3, & ! J/mol
         DSd = 660.0, & ! J/(mol K)
         zp = 1.0

    real(8) :: phi, ea, D0a, D0lT0, D0l, bt, ta, tl, Deffa, &
    Deffl, xCA, xCAref, Kh, B, D, mV, Ca, &
    k, z1

    phi = poros ! 0.4
    ea = MAX(0.01, phi - theta_w)

    ! Diffusivity, page 2224 Ogee
    D0a = 1.27d-5*(T_s/298.15)**1.5
    D0lT0 = D0l25/((298.15/216. - 1.0)**2.)
    D0l = D0lT0*((T_s/216. - 1.)**2.)

    ! Tortuosities, torttype = "Mol03r"i
    !bt = 3.
    !Moldrup et al., 2003
    !ta = (ea**(1.+3./bt))/(phi**(3./bt))  
    !tl = (theta_w**(bt/3.))/(phi**(bt/3.-1.)) 
    
    !Deepagoda et al., 2011, undisturbed soil
    ta = (0.2*(ea/phi)**2+0.004)/phi
    
    !Millington and Quirk, 1961.
    tl = (theta_w**(7/3))/(phi**2)

    Deffa = D0a*ta*ea
    Deffl = D0l*tl*theta_w

    ! Kuncat
    !kuncat = 2.15d-5 * exp(-14050.0*(1./T_s - 1./298.)) + &
    !    12.7*10**(-pKw+pH)*exp(-6040.0*(1./T_s - 1./298.))

    ! xCA
    xCA = exp(-DHa/(R*T_s))/(1. + exp(-DHd/(R*T_s) + DSd/R))

    xCAref = exp(-DHa/(R*293.))/(1. + exp(-DHd/(R*293.) + DSd/R))

    ! Kh, Henry's law constant, p 2223 Ogee
    Kh = 0.021d-2*exp(24900.*(1./T_s - 1./298.15)/R)

    ! B, related to Henry's law constant, p 2223 Ogee
    B = Kh*R*T_s

    D = Deffa + Deffl*B
    mV = R*T_s/pressure ! m3/mol
    Ca = coscas/mV !mol/m3

    k = fCA*kuncat*xCA/xCAref

    ! Production
    ! P = V_max_cos_source*Q10**((T_s-298.15)/10.)

    z1 = sqrt(D/(k*B*theta_w))
    Fcos = sqrt(k*B*theta_w*D)*(Ca - (z1*z1*P/D)*(1.-exp(-zp/z1)))

    !print*, 'Ogee solver: '
    !print*, 'D: ', D
    !print*, 'k: ', k
    !print*, 'B: ', B
    !print*, 'theta_w: ', theta_w
    !print*, 'Ca: ', Ca
    !print*, 'z1: ', z1
    !print*, 'P: ', P
    !print*, 'zp: ', zp

    !print *, "Fcos",Fcos

end subroutine


end module fcos_solver_Ogee
