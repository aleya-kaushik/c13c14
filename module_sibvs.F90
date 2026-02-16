module module_sibvs

!----------------------------------------------------------------------
!
!   SiB4 Vegetation Structure Module
!
!----------------------------------------------------------------------

use kinds
implicit none

!******************************************************************

!------------------------------------------------------------------
! SiB4 Vegetation Information
!------------------------------------------------------------------
type sib_vs_vars
  
  integer(byte) gnlu !Number of land units (lu)
                      ! per grid cell
                       
                      !Land Units consist of consistent
                      !  soil column (Soil), vegetation (PFT),
                      !  canopy (Can), and canopy air space (CAS)
 
  real(r4), dimension(:), allocatable :: & !(nlu)
          larea  !Fractional coverage of each land unit
  integer(i4), dimension(:), allocatable :: & !(nlvc)
          pftref !PFT reference for each lvc unit

  real(r8), dimension(:), allocatable :: & !(nlu)
        sandfrac, & 
        clayfrac, &
        soref_vis, &
        soref_nir

end type sib_vs_vars


!***********************************************************************
!-----------------------------------------------------------------------
! SiB4 Vegetation Structure Variables
  type(sib_vs_vars), allocatable :: sibvs(:)

!-----------------------------------------------------------------------

end module module_sibvs
