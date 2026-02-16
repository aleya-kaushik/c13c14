
! Opens and reads in PFT information.
subroutine read_pftinfo()

    use kinds
    use module_sibconst, only: &
        npft, ntype, ngroup, &
        npmeth, npstgmax,    &
        print_pftinfo, print_stop
    use module_io, only: &
        pft_info, pftid
    use module_pftinfo

    implicit none

    !...file variables
    integer(i4) :: ref, prefmax
    character(len=30) :: trash

    !...misc variables
    integer(byte) :: ibyte
    integer(i4) :: i

!---------------------------------
!...Open file
print*,''
print*,'Reading PFT Informational File: '
print*,'  ',trim(pft_info)
open(unit=pftid,file=trim(pft_info),form='formatted')

read(pftid,*) npft
read(pftid,*) ntype 
read(pftid,*) ngroup
read(pftid,*) npmeth
read(pftid,*) npstgmax

!...Read in map information
read(pftid,*) trash
read(pftid,*) trash
read(pftid,*) trash

read(pftid,*) trash, pft_source
read(pftid,*) trash, crop_source
read(pftid,*) trash, soil_source
read(pftid,*) trash, soref_source

!...Read in type information
read(pftid,*) trash
read(pftid,*) trash
read(pftid,*) trash
allocate(type_name_long(ntype))
allocate(type_name(ntype))
do i=1,ntype
    read(pftid,'(i6,a20,a8)') ref, type_name_long(i), &
                type_name(i)
enddo

!...Read in group information
read(pftid,*) trash
read(pftid,*) trash
read(pftid,*) trash
allocate(group_name_long(ngroup))
allocate(group_name(ngroup))
do i=1,ngroup
   read(pftid,'(i6,a20,a8)') ref, group_name_long(i), &
               group_name(i)
enddo

!...Read in phenology method information
read(pftid,*) trash
read(pftid,*) trash
read(pftid,*) trash
allocate(pmeth_name(npmeth))
do i=1,npmeth
   read(pftid,'(i6,a20)') ref, pmeth_name(i)
enddo

!...Read in PFT information
read(pftid,*) trash
read(pftid,*) trash
read(pftid,*) trash
allocate(pft_ref(npft))
allocate(pft_name(npft))
allocate(pft_name_long(npft))
allocate(pft_type(npft))
allocate(pft_group(npft))
allocate(pft_pmeth(npft))
do i=1,npft
   read(pftid,'(i6,a3,3i6,a40)') &
        pft_ref(i), pft_name(i), &
        pft_type(i), pft_group(i), &
        pft_pmeth(i), pft_name_long(i)
enddo

close(pftid)

!------------------------------------
!...Set PFT cross-reference/numbers
prefmax = maxval(pft_ref)
allocate(pft_num(prefmax))
do i=1,npft
   pft_num(pft_ref(i)) = i
enddo

!...Set Type References
type_nbare = izero
type_nevg = izero
type_ndecid = izero
type_ncrop = izero

do i=1,npft
   select case (pft_type(i))
   case (1)
       type_nbare = type_nbare + bone       
   case (2) 
       type_nevg = type_nevg + bone 
   case (3) 
       type_ndecid = type_ndecid + bone
   case (4)
       type_ngrass = type_ngrass + bone
   case (5) 
       type_ncrop = type_ncrop + bone
   case default
      print*,'Unknown Type Number: ',pft_type(i)
      print*,'Stopping'
      stop
      
   end select 
enddo

do i=1,ntype
   ibyte=int(real(i),KIND=byte)
   select case (trim(type_name(i)))
   case ('bare')
       type_bare = ibyte
   case ('evg')
       type_evg = ibyte
   case ('decid')
       type_decid = ibyte
   case ('grass')
       type_grass = ibyte
   case ('crop')
      type_crop = ibyte
   case default
      print*,'Unknown Type Name: ',trim(type_name(i))
      print*,'Stopping'
      stop
   end select
enddo

!------------------------------------
!...Set Group References
group_nbare = izero
group_nndlfor = izero
group_nbdlfor = izero
group_nshb = izero
group_ngrass = izero
group_ncrop = izero

do i=1,npft
   select case (pft_group(i))
   case (1)
       group_nbare = group_nbare + bone       
   case (2) 
       group_nndlfor = group_nndlfor + bone 
   case (3) 
       group_nbdlfor = group_nbdlfor + bone
   case (4)
       group_nshb = group_nshb + bone
   case (5) 
       group_ngrass = group_ngrass + bone
   case (6) 
       group_ncrop = group_ncrop + bone
    case default
       print*,'Unknown Group Number: ',pft_group(i)
       print*,'Stopping'
       stop
   end select 
enddo

do i=1,ngroup
   ibyte=int(i,kind=byte)
   select case (trim(group_name(i)))
   case ('bare')
       group_bare = ibyte
   case ('ndlfor')
       group_ndlfor = ibyte
   case ('bdlfor')
       group_bdlfor = ibyte
   case ('shrub')
       group_shrub = ibyte
   case ('grass')
       group_grass = ibyte
   case ('crop')
      group_crop = ibyte
   case default
      print*,'Unknown Group Name: ',trim(group_name(i))
      print*,'Stopping'
      stop      
   end select
enddo


!-----------------------------------
!...Set Phenology Method References
pmeth_nvg = bzero
pmeth_stg = bzero
pmeth_gdd = bzero

do i=1,npmeth
   ibyte=int(i,kind=byte)
   select case (trim(pmeth_name(i)))
   case ('Non-Veg')
      pmeth_nvg = ibyte
   case ('Stage-Based')
      pmeth_stg = ibyte
   case ('GDD-Based')
      pmeth_gdd = ibyte
   case default
      print*,'Unknown Phenology Method: ',trim(pmeth_name(i))
      print*,'Stopping'
      stop
   end select
enddo

npft_nvg = bzero
npft_stg = bzero
npft_gdd = bzero
allocate(nvgindx_pftref(prefmax))
nvgindx_pftref(:) = bzero
allocate(stgindx_pftref(prefmax))
stgindx_pftref(:) = bzero
allocate(gddindx_pftref(prefmax))
gddindx_pftref(:) = bzero

do i=1,npft
   ibyte= pft_pmeth(i)
   ref  = pft_ref(i)
   if (ibyte .eq. pmeth_nvg) then
       npft_nvg = npft_nvg + bone
       nvgindx_pftref(ref) = npft_nvg 
   elseif (ibyte .eq. pmeth_stg) then 
       npft_stg = npft_stg + bone
       stgindx_pftref(ref) = npft_stg
   elseif (ibyte .eq. pmeth_gdd) then 
       npft_gdd = npft_gdd + bone
       gddindx_pftref(ref) = npft_gdd
   else
      print*,'Unknown Phenology Method: ',ibyte
      print*,'Stopping'
      stop
   endif
enddo

allocate(pftnum_nvgindx(npft_nvg))
allocate(pftnum_stgindx(npft_stg))
allocate(pftnum_gddindx(npft_gdd))
npft_nvg = bzero
npft_stg = bzero
npft_gdd = bzero
do i=1, npft
   ibyte = pft_pmeth(i)
   if (ibyte .eq. pmeth_nvg) then
       npft_nvg = npft_nvg + bone
       pftnum_nvgindx(npft_nvg) = int(i,kind=byte)
   elseif (ibyte .eq. pmeth_gdd) then
      npft_gdd = npft_gdd + bone
      pftnum_gddindx(npft_gdd) = int(i,kind=byte)
   elseif (ibyte .eq. pmeth_stg) then
      npft_stg = npft_stg + bone
      pftnum_stgindx(npft_stg) = int(i,kind=byte)
   else
      print*,'Unknown Phenology Method: ',ibyte
      print*,'Stopping'
      stop
   endif
enddo

!------------------------------------
!...Set PFT references
do i=1,npft
   ibyte=int(i,kind=byte)

   select case (pft_name(i))
   case ('dbg')
        pft_dbg=ibyte
   case ('DBG')
        pft_dbg=ibyte

   case ('en2')
        pft_en2=ibyte
   case ('EN2')
        pft_en2=ibyte
   case ('en3')
        pft_en3=ibyte
   case ('EN3')
        pft_en3=ibyte

   case ('dnf')
        pft_dnf=ibyte
   case ('DNF')
        pft_dnf=ibyte

   case ('eb1')
        pft_eb1=ibyte
   case ('EB1')
        pft_eb1=ibyte
   case ('eb2')
        pft_eb2=ibyte
   case ('EB2')
        pft_eb2=ibyte

   case ('db1')
        pft_db1=ibyte
   case ('DB1')
        pft_db1=ibyte
   case ('db2')
        pft_db2=ibyte
   case ('DB2')
        pft_db2=ibyte
   case ('db3')
        pft_db3=ibyte
   case ('DB3')
        pft_db3=ibyte

   case ('shb')
        pft_shb=ibyte
   case ('SHB')
        pft_shb=ibyte

   case ('sha')
        pft_sha=ibyte
   case ('SHA')
        pft_sha=ibyte

   case ('c3a')
        pft_c3a=ibyte
   case ('C3A')
        pft_c3a=ibyte

   case ('c3g')
        pft_c3g=ibyte
   case ('C3G')
        pft_c3g=ibyte

   case ('c4g')
        pft_c4g=ibyte
   case ('C4G')
        pft_c4g=ibyte

   case ('c3c')
        pft_c3c=ibyte
   case ('C3C') 
        pft_c3c=ibyte

   case ('c4c')
        pft_c4c=ibyte
   case ('C4C')
        pft_c4c=ibyte

   case ('mze')
        pft_mze=ibyte
   case ('MZE')
        pft_mze=ibyte

   case ('soy')
        pft_soy=ibyte
   case ('SOY')
        pft_soy=ibyte

   case ('wwt')
        pft_wwt=ibyte
   case ('WWT')
        pft_wwt=ibyte

   case default
      print*,'Unknown PFT Name: ', pft_name(i)
      print*,'Stopping.'
      stop
   end select

enddo


!----------------------
!...Print Information
if (print_pftinfo) then
   print('(a,i4)'),'   Number of PFTs:', npft
   print('(a,i3)'),'   Number of Types:', ntype
   print('(2a,5i4)'),'      Type Disbtribution ', &
         '(evg/decid/grass/crop): ', &
         type_nevg,type_ndecid,type_ngrass,type_ncrop
   print('(a,i3)'),'   Number of Groups:', ngroup
   print('(2a,5i4)'),'      Group Distribution ', &
        '(ndlfor/bdlfor/shb/grass/crop):', &
         group_nndlfor, group_nbdlfor, group_nshb, &
         group_ngrass, group_ncrop
   print('(a,i3)'),'  Number of Phenology Methods: ', npmeth
   print('(2a,3i4)'),'      Method Distribution ', &
        '(non-veg/stage/gdd):', &
        npft_nvg, npft_stg, npft_gdd

   if (print_stop) stop
endif

end subroutine read_pftinfo
