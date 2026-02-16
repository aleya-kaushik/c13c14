module module_poolinfo

!----------------------------------------------------------------------
!
!   SiB4 Pool Informational Module
!
!----------------------------------------------------------------------

use kinds
implicit none


  character(len=23), dimension(:), allocatable :: pool_name_long ! pool names
  character(len=8), dimension(:), allocatable ::  pool_name      ! short pool names
  character(len=4), dimension(:), allocatable ::  pool_type      ! pool types (live, dead)
  character(len=10), dimension(:), allocatable ::  pool_loc       ! pool locations (soil, surface, canopy)
  integer(byte), dimension(:), allocatable ::     pool_indx_lay  ! index for lowest soil layer with carbon

  !...Pool indices
  integer(byte), dimension(:), allocatable :: &
       pool_indx_can,  &  ! pool index numbers for all canopy pools
       pool_indx_sfc,  &  ! pool index numbers for all surface pools
       pool_indx_soil,  &   ! pool index numbers for all soil pools
       !.. same as above but for C13 pools
       pool_indx_canc13,  &  ! pool index numbers for all canopy pools
       pool_indx_sfcc13,  &  ! pool index numbers for all surface pools
       pool_indx_soilc13     ! pool index numbers for all soil pools

       !.. same as above but for C14 pools
       pool_indx_canc14,  &  ! pool index numbers for all canopy pools
       pool_indx_sfcc14,  &  ! pool index numbers for all surface pools
       pool_indx_soilc14     ! pool index numbers for all soil pools

  integer(byte) pool_indx_leaf    ! pool index number for leaf pool
  integer(byte) pool_indx_froot   ! pool index number for fine root pool
  integer(byte) pool_indx_croot   ! pool index number for coarse root pool
  integer(byte) pool_indx_stwd    ! pool index number for stem/wood pool
  integer(byte) pool_indx_prod    ! pool index number for product pool
  integer(byte) pool_indx_cdb     ! pool index number for coarse dead biomass pool
  integer(byte) pool_indx_metl    ! pool index number for metabolic litter pool
  integer(byte) pool_indx_strl    ! pool index number for structural litter pool
  integer(byte) pool_indx_slit    ! pool index number for soil litter pool
  integer(byte) pool_indx_slow    ! pool index number for soil slow pool
  integer(byte) pool_indx_arm     ! pool index number for soil armored pool

  integer(byte) pool_indx_leaf_c13    ! pool index number for leaf c13 pool
  integer(byte) pool_indx_froot_c13   ! pool index number for fine root c13 pool
  integer(byte) pool_indx_croot_c13   ! pool index number for coarse root c13 pool
  integer(byte) pool_indx_stwd_c13    ! pool index number for stem/wood c13 pool
  integer(byte) pool_indx_prod_c13    ! pool index number for product c13 pool
  integer(byte) pool_indx_cdb_c13     ! pool index number for coarse dead biomass c13 pool
  integer(byte) pool_indx_metl_c13    ! pool index number for metabolic litter c13 pool
  integer(byte) pool_indx_strl_c13    ! pool index number for structural litter c13 pool
  integer(byte) pool_indx_slit_c13    ! pool index number for soil litter c13 pool
  integer(byte) pool_indx_slow_c13    ! pool index number for soil slow c13 pool
  integer(byte) pool_indx_arm_c13     ! pool index number for soil armored c13 pool

  integer(byte) pool_indx_leaf_c14    ! pool index number for leaf c14 pool
  integer(byte) pool_indx_froot_c14   ! pool index number for fine root c14 pool
  integer(byte) pool_indx_croot_c14   ! pool index number for coarse root c14 pool
  integer(byte) pool_indx_stwd_c14    ! pool index number for stem/wood c14 pool
  integer(byte) pool_indx_prod_c14    ! pool index number for product c14 pool
  integer(byte) pool_indx_cdb_c14     ! pool index number for coarse dead biomass c14 pool
  integer(byte) pool_indx_metl_c14    ! pool index number for metabolic litter c14 pool
  integer(byte) pool_indx_strl_c14    ! pool index number for structural litter c14 pool
  integer(byte) pool_indx_slit_c14    ! pool index number for soil litter c14 pool
  integer(byte) pool_indx_slow_c14    ! pool index number for soil slow c14 pool
  integer(byte) pool_indx_arm_c14     ! pool index number for soil armored c14 pool

end module module_poolinfo

