!==================SUBROUTINE RNLOAD====================================
subroutine radabs( dlwbot, gdiagt, radt )

    !=======================================================================
    !
    !    Calculation of absorption of radiation by surface.  
    !       Note that the output from this calculation, radc3, only accounts
    !       for the absorption of incident longwave and shortwave fluxes.
    !       
    !       The total net radiation calculation is in netrad.
    !
    !    Output: radc3c, radc3g - sum of absorbed radiative fluxes (W/m2)
    !=======================================================================

    use kinds
    use module_sib, only: &
       gdiag_type, &
       rad_type
 
    implicit none

    ! input variables
    real(r8), intent(in) :: dlwbot
    type(gdiag_type), intent(in) :: gdiagt
    type(rad_type), intent(inout) :: radt    

    ! local variables
    integer(i4) :: iwave, irad
    real(r8)    :: radc3c, radc3g, radn(2,2)

    !...reset variables
    radc3c = 0.
    radc3g = 0.
    radn(1,1) = gdiagt%radvbc
    radn(1,2) = gdiagt%radvdc
    radn(2,1) = gdiagt%radnbc
    radn(2,2) = gdiagt%radndc

    do iwave=1,2
       do irad=1,2
          radc3c = radc3c +  &
                   radt%radfacc(iwave,irad) * radn(iwave,irad)
          radc3g = radc3g + &
                   radt%radfacg(iwave,irad) * radn(iwave,irad)
       enddo
    enddo

    !...absorb downwelling radiation 
    radt%radc3c = radc3c + dlwbot * radt%effgc
    radt%radc3g = radc3g + dlwbot * (1.-radt%effgc)

end subroutine radabs
