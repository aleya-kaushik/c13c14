module cos_soil_production

implicit none

! #########################
!    COS soil emissions
! #########################
real(8) :: &

! COS soil emission parameters after Meredith et al., 2018, SS
! Meredith, L.K.; Boye, K.; Youngerman, C.; Whelan, M.; Ogée, J.; Sauze, J.;
! Wingate, L. Coupled Biological and Abiotic Mechanisms Driving Carbonyl Sulfide
! Production in Soils. Soil Syst. 2018, 2, 37.

!temperature ENF
alfa_enf = 4.86, & !pmol m-3 s-1
beta_enf = 0.1011, &

!temperate CRO
alfa_cro = 9.59, &   !pmol m-3 s-1
beta_cro = 0.1039, &

!temperate DBF
alfa_dbf = 4.94, & !pmol m-3 s-1
beta_dbf = 0.1074, &

!mediteranean GRA
alfa_gra = 2.20, & !pmol m-3 s-1
beta_gra = 0.0960, &

!arid desert
alfa_des = 5.60, & !pmol m-3 s-1
beta_des = 0.0499


contains

function COS_enf_production(soil_temp) result(forest_prod)
  real(8), intent(in) :: soil_temp
  real(8) :: forest_prod
  forest_prod = alfa_enf * exp(beta_enf * soil_temp)
end function

function COS_cro_production(soil_temp) result(ag_prod)
  real(8), intent(in) :: soil_temp
  real(8) :: ag_prod
  ag_prod = alfa_cro * exp(beta_cro * soil_temp)
end function

function COS_gra_production(soil_temp) result(grass_prod)
  real(8), intent(in) :: soil_temp
  real(8) :: grass_prod
  grass_prod = alfa_gra * exp(beta_gra * soil_temp)
end function

function COS_dbf_production(soil_temp) result(dbf_prod)
  real(8), intent(in) :: soil_temp
  real(8) :: dbf_prod
  dbf_prod = alfa_dbf * exp(beta_dbf * soil_temp)
end function

function COS_des_production(soil_temp) result(des_prod)
  real(8), intent(in) :: soil_temp
  real(8) :: des_prod
  des_prod = alfa_des * exp(beta_des * soil_temp)
end function


end module cos_soil_production
