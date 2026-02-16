!
!=====================SUBROUTINE CYCALC=================================
subroutine phocycalc( toa_par, par, aparkk,    &
        vmax, vmaxss, atheta, btheta, gamma, &
        rrkk, omss, c3, c4, pco2i,   &
        omc, ome, oms,               &
        assim)
!=======================================================================
! calculates smoothed minimum photosynthetic rate equivelent to steps
! in figure 4 of SE-92A.  C4 calculation based on CO-92
!
! OUTPUT:   
!       OMC    RUBISCO LIMITED ASSIMILATION (MOL M-2 S-1): EQUATION (11), SE-92A
!       OME    LIGHT LIMITED ASSIMILATION (MOL M-2 S-1): EQUATION (12), SE-92A
!       OMS    SINK LIMITED ASSIMILATION (MOL M-2 S-1): EQUATION (13), SE-92A
!       ASSIM  GROSS ASSIMILATION (MOL M-2 S-1): EQUATION (14,15), SE-92A
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use kinds

    implicit none

    !...INPUT VARIABLES
    real(r8),intent(in) :: &
        toa_par, & !Top-of-Atmosphere PAR (mol/m2/s)
        par,    & ! Photosynthetically active radiation (mol/m2/s)
        aparkk, & ! PAR-Use Parameter (-)
        vmax,   & ! Catalytic capacity of Rubisco (mol/m2/s) 
        vmaxss, & ! Temperature-Stressed catalytic capacity of Rubisco (mol/m2/s)
        atheta, & ! Coupling parameter (-)
        btheta, & ! Coupling parameter (-)
        gamma,  & ! CO2 photocompensation point (Pa)
        rrkk,   & ! Portion of sink-limited assimilation (mol/m2/s)
        omss,   & ! Portion of sink-limited assimilation (mol/m2/s)
        c3,     & ! Fraction of C3 vegetation (0-1)
        c4,     & ! Fraction of C4 vegetation (0-1)
        pco2i     ! Chloroplast CO2 conentration (Pa)

    !...OUTPUT VARIABLES
    real(r8),intent(inout) :: &
        omc,    & ! Rubisco-limited assimilation (mol/m2/s)
        ome,    & ! Light-limited assimilation (mol/m2/s)
        oms,    & ! Sink-limited assimilation (mol/m2/s)
        assim     ! Gross assimilation (mol/m2/s)

    !...LOCAL VARIABLES
    real(r8) :: &
        omp,   & ! Smoothed combo of rubisco and light limitations (mol/m2/s)
        sqrtin   ! Factor for calculating smoothed assimilation (mol/m2/s)

    real(r8) :: assimpot, lightpot

    !-----------------------------------------------------------------------
    !     CALCULATE ASSIMILATION RATE
    !-----------------------------------------------------------------------

    !...vmax assimilation
    assimpot = vmax *(pco2i-gamma)/(pco2i + rrkk)*c3    &
               + vmax * c4
    omc = vmaxss *(pco2i-gamma)/(pco2i + rrkk)*c3    &
        + vmaxss * c4

    !...light assimilation
    lightpot = toa_par*(pco2i-gamma)/(pco2i+2.*gamma)*c3 &
        + toa_par*c4
    ome = par*(pco2i-gamma)/(pco2i+2.*gamma)*c3 &
        + par * c4

    !...transfer assimilation
    oms  = omss * c3 + omss*pco2i * c4

    !...combined assimilation
    sqrtin= MAX( 0.0_r8, ( (ome+omc)**2 - 4.*atheta*ome*omc ) )
    omp  = ( ( ome+omc ) - SQRT( sqrtin ) ) / ( 2.*atheta )
    sqrtin= MAX( 0.0_r8, ( (omp+oms)**2 - 4.*btheta*omp*oms ) ) 
    assim = (( ( oms+omp ) - SQRT( sqrtin ) ) / ( 2.*btheta ))*aparkk

    return

end subroutine phocycalc   
