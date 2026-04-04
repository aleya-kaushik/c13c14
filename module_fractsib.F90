module module_fractsib

!----------------------------------------------------------------------
!   SiB4 C-13 Variables    
!----------------------------------------------------------------------

    use kinds

    implicit none

    !...cfrax variables
    real(r8) :: d13cca   ! isotope effects
    real(r8) :: d13cm     !
    real(r8) :: c13cca
    real(r8) :: c12cca 
    real(r8) :: c13ca     !
    real(r8) :: c12ca     !
    real(r8) :: c13cm     !
    real(r8) :: c12cm     !
    real(r8) :: kiecps     ! fractionation during photosynthesis
    real(r8) :: kiecps_nog ! fractionation during photosynthesis, no gamma term
    real(r8) :: rcassim    ! isotope ratio value of assimilation 
    real(r8) :: rcassimfac
    real(r8) :: d13cassim   ! delta value of assimilation
    real(r8) :: d13cassim_nog   ! delta value of assimilation, no gamma term
    real(r8) :: d14cassim   ! delta value of assimilation
    real(r8) :: c13assim    ! recently assimilated carbon-13 (mol C/m2/s)
    real(r8) :: c13assimd   ! daily (24-hour running mean) assimilated carbon-13 (mol C/m2/s)
    real(r8) :: c12assim    ! recently assimilated carbon-12 (mol C/m2/s)
    real(r8) :: c13resptot    ! recently respired carbon-13 (mol C/m2/s)
    real(r8) :: c12resptot    ! recently respired carbon-12 (mol C/m2/s)

    real(r8) :: kiecps_k1   ! fractionation term 1
    real(r8) :: kiecps_k2   ! fractionation term 2
    real(r8) :: kiecps_k3   ! fractionation term 3
    real(r8) :: kiecps_k4   ! fractionation term 4
    real(r8) :: kiecps_k5   ! fractionation term 5

    !c14 variables as above
    real(r8) :: d14cca   ! isotope effects
    real(r8) :: d14cm     ! 
    real(r8) :: c14ca     !
    real(r8) :: c14cca
    real(r8) :: c14cm     !
    real(r8) :: rcassimc14    ! isotope ratio value of assimilation 
    real(r8) :: rcassimfacc14
    real(r8) :: d14cassim   ! delta value of assimilation
    real(r8) :: c14assim    ! recently assimilated carbon-13 (mol C/m2/s)
    real(r8) :: c14assimn   ! recently net assimilated carbon-14 (mol C/m2/s)
    real(r8) :: c14assimd   ! daily (24-hour running mean) assimilated carbon-14 (mol C/m2/s)
    real(r8) :: c14resptot    ! recently respired carbon-13 (mol C/m2/s)

    real(r8) :: c14alpha  ! 14alpha_ph = (13alpha_ph)^2, alpha_ab = 1 + wtkiecps/1000.

    !...pool variables
    real(r8) :: d13cpool_leafc13 
    real(r8) :: d13cpool_frootc13
    real(r8) :: d13cpool_crootc13
    real(r8) :: d13cpool_stwdc13
    real(r8) :: d13cpool_prodc13
    real(r8) :: d13cpool_cdbc13
    real(r8) :: d13cpool_metlc13
    real(r8) :: d13cpool_strlc13
    real(r8) :: d13cpool_slitc13
    real(r8) :: d13cpool_slowc13
    real(r8) :: d13cpool_armc13

    real(r8) :: d14cpool_frootc14
    real(r8) :: d14cpool_crootc14
    real(r8) :: d14cpool_stwdc14
    real(r8) :: d14cpool_prodc14
    real(r8) :: d14cpool_cdbc14
    real(r8) :: d14cpool_metlc14
    real(r8) :: d14cpool_strlc14
    real(r8) :: d14cpool_slitc14
    real(r8) :: d14cpool_slowc14
    real(r8) :: d14cpool_armc14

    !...resp variables
    real(r8) :: d13cresp_autoc13
    real(r8) :: d13cresp_growc13
    real(r8) :: d13cresp_leafc13
    real(r8) :: d13cresp_mntnc13
    real(r8) :: d13cresp_rootc13
    real(r8) :: d13cresp_hetc13
    real(r8) :: d13cresp_firec13
    real(r8) :: d13cresp_grzc13
    real(r8) :: d13cresp_hrvstc13
    real(r8) :: d13cresp_nvegc13
    real(r8) :: d13cresp_soilc13
    real(r8) :: d13cresp_totc13

    real(r8) :: d14cresp_autoc14
    real(r8) :: d14cresp_growc14
    real(r8) :: d14cresp_leafc14
    real(r8) :: d14cresp_mntnc14
    real(r8) :: d14cresp_rootc14
    real(r8) :: d14cresp_hetc14
    real(r8) :: d14cresp_firec14
    real(r8) :: d14cemisfire_pool
    real(r8) :: d14cresp_grzc14
    real(r8) :: d14cresp_hrvstc14
    real(r8) :: d14cresp_nvegc14
    real(r8) :: d14cresp_soilc14
    real(r8) :: d14cresp_totc14

    !...met varaibles

    real(r8) :: co2m_cfrax
    real(r8) :: press_cfrax
    real(r8) :: press_cfraxps

    !...fire variables

    real(r8) :: rcpoolfirec13
    real(r8) :: rcpoolfirec14
    real(r8) :: poolemistotC
    real(r8) :: poolemisc13
    real(r8) :: poolemisc14

end module module_fractsib
