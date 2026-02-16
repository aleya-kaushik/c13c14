module module_isodata

!----------------------------------------------------------------------
!   SiB4 C-13 Variables    
!----------------------------------------------------------------------

    use kinds

    implicit none
    
    !...Variables from iso data file in /params/c_iso_time_series.dat
    real(r8), dimension(:), allocatable :: isoyr
    real(r8), dimension(:), allocatable :: nhc14
    real(r8), dimension(:), allocatable :: tropc14
    real(r8), dimension(:), allocatable :: shc14
    real(r8), dimension(:), allocatable :: globc13

end module module_isodata
