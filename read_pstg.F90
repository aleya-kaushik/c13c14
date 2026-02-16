
! Opens and reads in SiB phenological parameters
!   for the stage-based method.
subroutine read_pstg()

use kinds
use module_sibconst, only: &
    npft, npoolpft, npoollu
use module_io, only: &
    pstg_file, pstgid
use module_param, only: &
    phencon
use module_pftinfo, only: &
    clen, pft_dbg, pft_ref, pft_name, &
    npft_stg, pftnum_stgindx
use module_time, only: steps_per_day

implicit none

!...file variables
integer(i4) :: finpft, finpftstg, finpoolpft
integer(i4) :: npstg
character(len=clen) :: pftname

character(len=16) :: trash
logical :: iscomment1, iscomment2

!...allocation check variables
real(r4) :: asum, asumc13
logical :: aok, aokc13

!...misc variables
integer(i4) :: i, iref, ref, s

!---------------------

!...Open file
print*,'Reading Stage-Based Phenology File: '
print*,'  ',trim(pstg_file)
open(unit=pstgid,file=trim(pstg_file),form='formatted')

read(pstgid,*) finpft
read(pstgid,*) finpftstg
read(pstgid,*) finpoolpft
read(pstgid,*) npstg

do i=1, npft_stg
   iref = pftnum_stgindx(i)
   phencon(iref)%npstg = npstg
enddo

!.....check file vs simulation constants
if (finpft /= npft) then
    print*,''
    print('(a)'),'!!!Stage-Phen File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Phen file npft: ', &
          finpft,' Sim npft: ',npft
    print*,''
    stop
endif

if (finpftstg /= npft_stg) then
    print*,''
    print('(a)'),'!!!Stage-Phen File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Phen file npft_stg: ', &
          finpftstg,' Sim npft_stg: ',npft_stg
    print*,''
    stop
endif

if (finpoolpft /= npoolpft) then
    print*,''
    print('(a)'),'!!!Stage-Phen File Mismatch!!!'
    print('(a,i4,a,i4)'),'  Phen file npoolpft: ', &
          finpoolpft, ' Sim npoolpft: ',npoolpft
    print*,''
    stop
endif

!...Read in comments
iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(pstgid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else 
           iscomment2=.false.
        endif
    endif
enddo

!...Read in the growing season start parameters
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i) 
   read(pstgid,*) ref, pftname, &
            phencon(iref)%daylen_mini, phencon(iref)%daylen_offd, &
            phencon(iref)%tawftop_min, phencon(iref)%tawftop_len, &
            phencon(iref)%tm_min, phencon(iref)%tm_max, phencon(iref)%tm_len
enddo

!...Read in the growing season reset parameters
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i) 
   read(pstgid,*) ref, pftname, &
            phencon(iref)%assim_resetl, phencon(iref)%assim_resetv
enddo

!...Read in the phenology stage daylength potential parameters
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, &
        phencon(iref)%psdayl_ref, phencon(iref)%psdayl_mul, &
        phencon(iref)%psdayl_min
enddo

!...Read in the phenology stage growth potential parameters
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, &
        phencon(iref)%climp_a, &
        phencon(iref)%climp_b, phencon(iref)%climp_c, &
        phencon(iref)%climp_d, &
        phencon(iref)%climp_min, phencon(iref)%climp_max, &
        phencon(iref)%cwa_type
enddo

read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, &
        phencon(iref)%psg_min, phencon(iref)%clai_coef, &
        phencon(iref)%clai_offl, phencon(iref)%clai_offg
enddo

!...Read in the phenology stage weather potential parameters
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, &
        phencon(iref)%pswx_type, phencon(iref)%pswx_rml
enddo

!...Read in phenology stage index thresholds
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, phencon(iref)%threshp(1:npstg-1)
enddo

!...Read in allocation parameters
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash

do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, phencon(iref)%allocp(1:5,1:3)
   phencon(iref)%allocp(6:10,1:3) = phencon(iref)%allocp(1:5,1:3)
enddo

read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, phencon(iref)%allocp(1:5,4:5), &
        phencon(iref)%adj_moist, phencon(iref)%adj_temp
   phencon(iref)%allocp(6:10,4:5) = phencon(iref)%allocp(1:5,4:5)
enddo

!...Read in the leaf pool transfer fractions
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, &
        phencon(iref)%lptransfer(1:npstg)
enddo

!...Read in the vmax
read(pstgid,*) trash
read(pstgid,*) trash
read(pstgid,*) trash
do i=1, npft_stg
   iref = pftnum_stgindx(i)
   read(pstgid,*) ref, pftname, &
        phencon(iref)%vmax0(1:npstg)
enddo

close(pstgid)


!------Perform conversions and calculate dependent parameters----
do i=1, npft_stg
   iref = pftnum_stgindx(i)

   !...set averaging weights
   phencon(iref)%wt_assim = 1./MAX(1.0,real(phencon(iref)%assim_resetl)*real(steps_per_day))
   phencon(iref)%wt_tawftop = 1./MAX(1.0,real(phencon(iref)%tawftop_len)*real(steps_per_day))
   phencon(iref)%wt_tm = 1./MAX(1.0,real(phencon(iref)%tm_len)*real(steps_per_day))
   phencon(iref)%wt_pswx = 1./MAX(1.0,real(phencon(iref)%pswx_rml)*real(steps_per_day))
enddo

!------Check allocation fractions----
do i=1, npft_stg
   IF (i .ne. pft_dbg) THEN
       do s=1, phencon(i)%npstg
          aok = .true.
          asum = sum(phencon(i)%allocp(1:5,s))
          if (asum .eq. rzero) then
              if ((s .gt. 1) .and. &
                  (s .lt. phencon(i)%npstg)) then
                  aok = .false.
              endif
          elseif ((asum .lt. 0.99) .or. (asum .gt. 1.01) .or. &
                  (minval(phencon(i)%allocp(1:5,s)) .lt. rzero)) then
                  aok = .false.
          endif

          if (.not. aok) then
              print*,'---Error with Allocation Factors---'
              print('(a,2i5,a3,a3)'),'   PFT Num/Ref/Name: ',i,pft_ref(i),'   ',pft_name(i)
              print('(a,i3)'),'   Phenology Stage: ', s
              print('(a,f8.3)'),'   Phenology Factor Total: ', asum
              print('(a)'), &
                  '   Phenology Factor Contribuations: (leaf/froot/croot/wood/prod):'
              print('(a,5f8.3)'),'  ', phencon(i)%allocp(1:5,s)

              print('(a)'),'Stopping.'
              stop
         endif
      enddo !s=1,npstg
   ENDIF !not pft_dbg
enddo !i=1,npft
!same check for c13
do i=1, npft_stg
   IF (i .ne. pft_dbg) THEN
       do s=1, phencon(i)%npstg
          aokc13 = .true.
          asumc13 = sum(phencon(i)%allocp(6:10,s))
          if (asumc13 .eq. rzero) then
              if ((s .gt. 1) .and. &
                  (s .lt. phencon(i)%npstg)) then
                  aokc13 = .false.
              endif
          elseif ((asumc13 .lt. 0.99) .or. (asumc13 .gt. 1.01) .or. &
                  (minval(phencon(i)%allocp(6:10,s)) .lt. rzero)) then
                  aokc13 = .false.
          endif

          if (.not. aokc13) then
              print*,'---Error with C13 Allocation Factors---'
              print('(a,2i5,a3,a3)'),'   PFT Num/Ref/Name: ',i,pft_ref(i),'   ',pft_name(i)
              print('(a,i3)'),'   Phenology Stage: ', s
              print('(a,f8.3)'),'   Phenology Factor Total: ', asumc13
              print('(a)'), &
                  '   Phenology Factor Contribuations: (leaf/froot/croot/wood/prod):'
              print('(a,5f8.3)'),'  ', phencon(i)%allocp(6:10,s)

              print('(a)'),'Stopping.'
              stop
         endif
      enddo !s=1,npstg
   ENDIF !not pft_dbg
enddo !i=1,npft

end subroutine read_pstg

