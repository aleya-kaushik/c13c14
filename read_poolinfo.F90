
! Opens and reads in pool information.
subroutine read_poolinfo()

    use kinds
    use module_sibconst, only: &
         npoolpft, npoollu, ntpool, &
         npoolcan, npoolsfc, npoolsoil, &
         npoolcanc13, npoolsfcc13, &
         npoolsoilc13, &
         npoolcanc14, npoolsfcc14, &
         npoolsoilc14
    use module_io, only: &
         pool_info, piid
    use module_poolinfo

    implicit none

    !...file variables
    integer(i4) :: num
    character(len=30) :: trash

    !...misc variables
    integer(byte) :: iref
    integer(i4) :: i
!---------------------------------
!...Open file
print*,''
print*,'Reading Pool Informational File: '
print*,'  ',trim(pool_info)
open(unit=piid,file=trim(pool_info),form='formatted')

read(piid,*) npoolpft
read(piid,*) npoollu
ntpool = npoolpft + npoollu

read(piid,*) trash
read(piid,*) trash
read(piid,*) trash
read(piid,*) trash
allocate(pool_name_long(ntpool))
allocate(pool_name(ntpool))
allocate(pool_type(ntpool))
allocate(pool_loc(ntpool))
allocate(pool_indx_lay(ntpool))
do i=1,ntpool
   read(piid,'(i6,a23,a2,a8,a2,a4,a4,a10,a2,i2)') &
        num, pool_name_long(i), trash, &
        pool_name(i), trash, &
        pool_type(i), trash, &
        pool_loc(i), trash, &
        pool_indx_lay(i)
enddo

close(piid)

!-----Categorize pools-----------------------------------------
!...Clear out all index variables
npoolcan=0
npoolsfc=0
npoolsoil=0

npoolcanc13=0
npoolsfcc13=0
npoolsoilc13=0

npoolcanc14=0
npoolsfcc14=0
npoolsoilc14=0

pool_indx_leaf=0
pool_indx_froot=0
pool_indx_croot=0
pool_indx_stwd=0
pool_indx_prod=0
pool_indx_cdb=0
pool_indx_metl=0
pool_indx_strl=0
pool_indx_slit=0
pool_indx_slow=0
pool_indx_arm=0

pool_indx_leaf_c13=0
pool_indx_froot_c13=0
pool_indx_croot_c13=0
pool_indx_stwd_c13=0
pool_indx_prod_c13=0
pool_indx_cdb_c13=0
pool_indx_metl_c13=0
pool_indx_strl_c13=0
pool_indx_slit_c13=0
pool_indx_slow_c13=0
pool_indx_arm_c13=0

pool_indx_leaf_c14=0
pool_indx_froot_c14=0
pool_indx_croot_c14=0
pool_indx_stwd_c14=0
pool_indx_prod_c14=0
pool_indx_cdb_c14=0
pool_indx_metl_c14=0
pool_indx_strl_c14=0
pool_indx_slit_c14=0
pool_indx_slow_c14=0
pool_indx_arm_c14=0

!...scan through pool information and count pools of each type
pool_type(:) = adjustl(pool_type(:))
pool_loc(:) = adjustl(pool_loc(:))
pool_name(:) = adjustl(pool_name(:))

do i=1,ntpool
   if(trim(pool_loc(i))=='canopy') &
       npoolcan = npoolcan + bone
   if(trim(pool_loc(i))=='surface') &
       npoolsfc = npoolsfc + bone
   if(trim(pool_loc(i))=='soil') &
        npoolsoil = npoolsoil + bone

   if(trim(pool_loc(i))=='canopyc13') &
       npoolcanc13 = npoolcanc13 + bone
   if(trim(pool_loc(i))=='surfacec13') &
       npoolsfcc13 = npoolsfcc13 + bone
   if(trim(pool_loc(i))=='soilc13') &
        npoolsoilc13 = npoolsoilc13 + bone

   if(trim(pool_loc(i))=='canopyc14') &
       npoolcanc14 = npoolcanc14 + bone
   if(trim(pool_loc(i))=='surfacec14') &
       npoolsfcc14 = npoolsfcc14 + bone
   if(trim(pool_loc(i))=='soilc14') &
        npoolsoilc14 = npoolsoilc14 + bone

    iref=int(i,kind=byte)
    if(trim(pool_name(i))=='leaf') pool_indx_leaf=iref
    if(trim(pool_name(i))=='froot') pool_indx_froot=iref
    if(trim(pool_name(i))=='croot') pool_indx_croot=iref
    if(trim(pool_name(i))=='stwd') pool_indx_stwd=iref
    if(trim(pool_name(i))=='prod') pool_indx_prod=iref
    if(trim(pool_name(i))=='cdb')     pool_indx_cdb=iref
    if(trim(pool_name(i))=='metl')  pool_indx_metl=iref
    if(trim(pool_name(i))=='strl')  pool_indx_strl=iref
    if(trim(pool_name(i))=='slit')    pool_indx_slit=iref
    if(trim(pool_name(i))=='slow')    pool_indx_slow=iref
    if(trim(pool_name(i))=='arm')     pool_indx_arm=iref

    if(trim(pool_name(i))=='leafc13') pool_indx_leaf_c13=iref
    if(trim(pool_name(i))=='frootc13') pool_indx_froot_c13=iref
    if(trim(pool_name(i))=='crootc13') pool_indx_croot_c13=iref
    if(trim(pool_name(i))=='stwdc13') pool_indx_stwd_c13=iref
    if(trim(pool_name(i))=='prodc13') pool_indx_prod_c13=iref
    if(trim(pool_name(i))=='cdbc13')     pool_indx_cdb_c13=iref
    if(trim(pool_name(i))=='metlc13')  pool_indx_metl_c13=iref
    if(trim(pool_name(i))=='strlc13')  pool_indx_strl_c13=iref
    if(trim(pool_name(i))=='slitc13')    pool_indx_slit_c13=iref
    if(trim(pool_name(i))=='slowc13')    pool_indx_slow_c13=iref
    if(trim(pool_name(i))=='armc13')     pool_indx_arm_c13=iref

    if(trim(pool_name(i))=='leafc14') pool_indx_leaf_c14=iref
    if(trim(pool_name(i))=='frootc14') pool_indx_froot_c14=iref
    if(trim(pool_name(i))=='crootc14') pool_indx_croot_c14=iref
    if(trim(pool_name(i))=='stwdc14') pool_indx_stwd_c14=iref
    if(trim(pool_name(i))=='prodc14') pool_indx_prod_c14=iref
    if(trim(pool_name(i))=='cdbc14')     pool_indx_cdb_c14=iref
    if(trim(pool_name(i))=='metlc14')  pool_indx_metl_c14=iref
    if(trim(pool_name(i))=='strlc14')  pool_indx_strl_c14=iref
    if(trim(pool_name(i))=='slitc14')    pool_indx_slit_c14=iref
    if(trim(pool_name(i))=='slowc14')    pool_indx_slow_c14=iref
    if(trim(pool_name(i))=='armc14')     pool_indx_arm_c14=iref

enddo


!...set the pool type/location indices
allocate(pool_indx_can(npoolcan))
allocate(pool_indx_canc13(npoolcanc13))
allocate(pool_indx_canc14(npoolcanc14))
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'canopy') then 
       pool_indx_can(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'canopyc13') then
       pool_indx_canc13(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'canopyc14') then
       pool_indx_canc14(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo


allocate(pool_indx_sfc(npoolsfc))
allocate(pool_indx_sfcc13(npoolsfcc13))
allocate(pool_indx_sfcc14(npoolsfcc14))
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'surface') then
       pool_indx_sfc(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'surfacec13') then
       pool_indx_sfcc13(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'surfacec14') then
       pool_indx_sfcc14(iref)=int(i,kind=byte)
       iref = iref + bone
   endif
enddo

!debugging checks
!print*,'pool_indx_sfc :',pool_indx_sfc(:)
!print*,'pool_indx_sfcc13 :',pool_indx_sfcc13(:)
!print*,'pool_indx_sfcc14 :',pool_indx_sfcc14(:)


allocate(pool_indx_soil(npoolsoil))
allocate(pool_indx_soilc13(npoolsoilc13))
allocate(pool_indx_soilc14(npoolsoilc14))
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'soil') then
      pool_indx_soil(iref)=int(i,kind=byte)
      iref=iref+bone
   endif
enddo
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'soilc13') then
      pool_indx_soilc13(iref)=int(i,kind=byte)
      iref=iref+bone
   endif
enddo
iref=bone
do i=1, ntpool
   if (trim(pool_loc(i)) == 'soilc14') then
      pool_indx_soilc14(iref)=int(i,kind=byte)
      iref=iref+bone
   endif
enddo

!debugging checks
!print*,'pool_indx_soil :',pool_indx_soil(:)
!print*,'pool_indx_soilc13 :',pool_indx_soilc13(:)
!print*,'pool_indx_soilc14 :',pool_indx_soilc14(:)

!----------------------
!...Print Information
print('(a,3i4)'),'   Pool Number (Tot/PFT/LU):         ', ntpool, npoolpft, npoollu
print('(a,3i4)'),'   Pool Location (Can/Sfc/Soil):     ', npoolcan, npoolsfc, npoolsoil
print('(a,3i4)'),'   C13 Pool Location (Can/Sfc/Soil): ', &
                          npoolcanc13, npoolsfcc13,npoolsoilc13
print('(a,3i4)'),'   C14 Pool Location (Can/Sfc/Soil): ', &
                          npoolcanc14, npoolsfcc14,npoolsoilc14
print*,''

end subroutine read_poolinfo
