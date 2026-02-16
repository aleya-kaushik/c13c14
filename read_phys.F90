! Opens and reads in SiB physiological parameters.
subroutine read_phys ()

    use kinds
    use module_sibconst, only: &
        grazing_switch, npft, nsoil
    use module_io, only: &
        phys_file, physid
    use module_param, only: &
        physcon
    use module_pftinfo, only: &
        clen

    !...file variables
    integer(i4) :: finpft
    character(len=clen) :: pftname
    logical :: iscomment1, iscomment2

    !...local variables
    integer(i4) :: i, num
    character(len=16) :: trash
    integer(i4) :: tempgrzswitch

!----------------------------
!...Allocate the physiological parameters
allocate(physcon(npft))

!...Open file
print*,'Reading Physiological File: '
print*,'  ',trim(phys_file)
open(unit=physid,file=trim(phys_file),form='formatted')

read(physid,*) finpft
if (finpft /= npft) then
    print*,''
    print('(a)'),'!!!Physiological file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Phys file npft: ',finpft,' Sim npft: ',npft
    print*,''
    stop
endif

!...Read in comments
iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(physid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else 
           iscomment2=.false.
        endif
    endif
enddo

!...Read in vegetation information
read(physid,*) trash
read(physid,*) trash
read(physid,*) trash

do i=1, npft
   read(physid,*) num, pftname, physcon(i)%sla,  &
              physcon(i)%laimin, physcon(i)%laisat, &
              physcon(i)%fparsat, tempgrzswitch

  physcon(i)%pftgraze = .false.
  if ((grazing_switch) .and. (tempgrzswitch .gt. 0)) &
        physcon(i)%pftgraze=.true.
     
enddo

read(physid,*) trash
read(physid,*) trash

do i=1,npft
   read(physid,*) num, pftname, &
               physcon(i)%c4flag, physcon(i)%chil, &
               physcon(i)%z1,physcon(i)%z2,        &
               physcon(i)%kroot, physcon(i)%rootd
enddo

!...Read in temperature stress parameters
do i=1,3
   read(physid,*) trash
enddo

do i=1,npft
   read(physid,*) num,pftname, &
              physcon(i)%slti, physcon(i)%shti, &
              physcon(i)%hlti, physcon(i)%hhti, physcon(i)%hfti,  &
              physcon(i)%sfti
enddo

!...Read in soil moisture stress parameters
do i=1,3
   read(physid,*) trash
enddo

do i=1,npft
   allocate(physcon(i)%wp_min(nsoil/2))
   read(physid,*) num, pftname, &
             physcon(i)%wssp, physcon(i)%wp_min(:)
enddo

do i=1,3
   read(physid,*) trash
enddo

do i=1,npft
   allocate(physcon(i)%fc_min(nsoil/2))
   read(physid,*) num, pftname, &
        physcon(i)%fc_min(:)
enddo

!...Read in photosynthesis parameters
do i=1,3
   read(physid,*) trash
enddo

do i=1,npft
   read(physid,*) num, pftname, &
              physcon(i)%effcon, physcon(i)%gmeso, &
              physcon(i)%binter, physcon(i)%gradm, &
              physcon(i)%atheta, physcon(i)%btheta, physcon(i)%gmin
enddo

do i=1,3
   read(physid,*) trash
enddo

!...Read in radiation parameters
do i=1,npft
   read(physid,*) num, pftname,  &
              physcon(i)%tran(1,1), physcon(i)%tran(2,1), &
              physcon(i)%tran(1,2), physcon(i)%tran(2,2),  &
              physcon(i)%ref(1,1), physcon(i)%ref(2,1),    &
              physcon(i)%ref(1,2), physcon(i)%ref(2,2)
enddo

close(physid)


!------Perform conversions and calculate dependent parameters----

!...check fparsat
if (maxval(physcon%fparsat) .gt. 1.0) then
   print*,'Error In Phenology Parameters: fparsat > 1.0'
   stop
endif

!...change units on SLA from cm2/g to m2/g
physcon%sla = physcon%sla * .0001

end subroutine read_phys

