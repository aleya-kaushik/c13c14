module module_pftinfo

!----------------------------------------------------------------------
!
!   SiB4 PFT Informational Module
!
!----------------------------------------------------------------------

use kinds
implicit none

   !...Parameters
   integer(i4), parameter :: alen=8  !Abreviated type/group name lengths
   integer(i4), parameter :: clen=3  !PFT name length
   integer(i4), parameter :: llen=40 !Long name lengths
   integer(i4), parameter :: plen=20 !Phenology method/stage name lengths
   
   !...Map information
   character(len=llen) :: pft_source    !PFT Reference
   character(len=llen) :: crop_source   !Crop Reference
   character(len=llen) :: soil_source   !Soil Type Reference
   character(len=llen) :: soref_source  !Soil Reflectance Reference

   !...Type Information (bare, evg, decid, grass, crop)
   character(len=llen), dimension(:), allocatable :: &
        type_name_long !Full type name
   character(len=alen), dimension(:), allocatable :: &
        type_name !Abbreviated type name

   !......type references
   integer(byte) :: type_bare    !Type number for bare ground
   integer(byte) :: type_evg     !Type number for evergreen
   integer(byte) :: type_decid   !Type number for deciduous
   integer(byte) :: type_grass   !Type number for grass
   integer(byte) :: type_crop    !Type number for crops
!   integer(byte) :: type_shrub   !Type number for shrubs

   !...Group Information (bare, ndlfor, bdlfor, shrub, grass, crop)
   character(len=llen), dimension(:), allocatable :: &
        group_name_long !Full type name
   character(len=alen), dimension(:), allocatable  :: &
        group_name      !Abbreviated type name 

   !......group references
   integer(byte) :: group_bare    !Group number for bare ground
   integer(byte) :: group_ndlfor  !Group number for needle forests
   integer(byte) :: group_bdlfor  !Group number for broadleaf forests
   integer(byte) :: group_shrub   !Group number for shrub
   integer(byte) :: group_grass   !Group number for grass
   integer(byte) :: group_crop    !Group number for crop

   !...Phenology Method (non-veg, stage-based, gdd-based)
   character(len=plen), dimension(:), allocatable :: &
                    pmeth_name  !Method name

   !.....method references
   integer(byte) :: pmeth_nvg !Method # for non-vegetation
   integer(byte) :: pmeth_stg !Method # for stage-based phenology
   integer(byte) :: pmeth_gdd !Method # for gdd-based phenology

   !...PFT Information
   character(len=clen), dimension(:), allocatable :: pft_name ! PFT names
   character(len=llen), dimension(:), allocatable :: pft_name_long
   integer(i4), dimension(:), allocatable :: pft_ref     ! PFT reference numbers
   integer(i4), dimension(:), allocatable :: pft_num     ! PFT sequential numbers (no skips)
   integer(byte), dimension(:), allocatable :: pft_type  ! PFT types (1-5) 
   integer(byte), dimension(:), allocatable :: pft_group ! PFT group
   integer(byte), dimension(:), allocatable :: pft_pmeth ! PFT phenology method

   !......PFT type counts
   integer(byte) :: type_nbare    !Number of bare ground PFTs
   integer(byte) :: type_nevg     !Number of evergreen PFTs
   integer(byte) :: type_ndecid   !Number of deciduous PFTs
   integer(byte) :: type_ngrass   !Number of grass PFTs
   integer(byte) :: type_ncrop    !Number of crop PFTs

   !......PFT group counts
   integer(byte) :: group_nbare   !Number of bare PFTs
   integer(byte) :: group_nndlfor !Number of needle forest PFTs
   integer(byte) :: group_nbdlfor !Number of broadleaf forest PFTs
   integer(byte) :: group_nshb    !Number of shrub PFTs
   integer(byte) :: group_ngrass  !Number of grass PFTs
   integer(byte) :: group_ncrop   !Number of crop PFTs

   !.....PFT phenology method counts and indices
   integer(byte) :: npft_nvg !Number of non-veg PFTs
   integer(byte) :: npft_stg !Number of stage-based PFTs
   integer(byte) :: npft_gdd !Number of gdd-based PFTs

   integer(byte), dimension(:), allocatable :: &  !(max PFT ref)
         nvgindx_pftref,  &  !nvg index number for PFT references
         gddindx_pftref,  &  !gdd index number for PFT references
         stgindx_pftref      !stg index number for PFT references

    integer(byte), dimension(:), allocatable :: & !(npft_nvg)
         pftnum_nvgindx  !PFT numbers for non-veg PFTs
    integer(byte), dimension(:), allocatable :: & !(npft_gdd)
         pftnum_gddindx  !PFT numbers for gdd-based PFTs
    integer(byte), dimension(:), allocatable :: & !(npft_stg)
         pftnum_stgindx !PFT numbers for stage-based PFTs

   !.....PFT Numbers
   integer(byte) :: pft_dbg  !PFT number for desert and bare ground
   integer(byte) :: pft_en2  !PFT number for temperate evergreen needleleaf forest
   integer(byte) :: pft_en3  !PFT number for boreal evergreen needleleaf forest
   integer(byte) :: pft_dnf  !PFT number for boreal deciduous needleleaf forest
   integer(byte) :: pft_eb1  !PFT number for tropical evergreen broadleaf forest
   integer(byte) :: pft_eb2  !PFT number for temperate evergreen broadleaf forest
   integer(byte) :: pft_db1  !PFT number for tropical deciduous broadleaf forest
   integer(byte) :: pft_db2  !PFT number for temperate deciduous broadleaf forest
   integer(byte) :: pft_db3  !PFT number for boreal deciduous broadleaf forest

   integer(byte) :: pft_shb  !PFT number for shrubs (non-tundra)
   integer(byte) :: pft_sha  !PFT number for tundra shrubs
   integer(byte) :: pft_c3a  !PFT number for tundra grasslands
   integer(byte) :: pft_c3g  !PFT number for C3 grasslands (non-tundra)
   integer(byte) :: pft_c4g  !PFT number for C4 grasslands

   integer(byte) :: pft_c3c  !PFT number for C3 crops (generic)
   integer(byte) :: pft_c4c  !PFT number for C4 crops (generic)
   integer(byte) :: pft_mze  !PFT number for maize
   integer(byte) :: pft_soy  !PFT number for soybeans
   integer(byte) :: pft_wwt  !PFT number for winter wheat


end module module_pftinfo
