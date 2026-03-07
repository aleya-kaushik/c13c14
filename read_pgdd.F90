
! Opens and reads in SiB phenological parameters
!   for the gpp-based method.
subroutine read_pgdd()

use kinds
use module_io, only: &
    pgdd_file, pgddid
use module_param, only: &
    phencon
use module_pftinfo, only:  &
    clen, npft_gdd, pftnum_gddindx
use module_poolinfo, only: &
    pool_indx_leaf, pool_indx_froot, &
    pool_indx_stwd, pool_indx_prod, &
    pool_indx_leaf_c13, pool_indx_froot_c13, &
    pool_indx_stwd_c13, pool_indx_prod_c13, &
    pool_indx_leaf_c14, pool_indx_froot_c14, &
    pool_indx_stwd_c14, pool_indx_prod_c14
use module_pparams, only: mwc
use module_sibconst, only: &
    npft, npoolpft, npstgmax, npoollu
use module_time, only: steps_per_day

implicit none

!...file variables
integer(i4) :: finpft, finpftgdd, finpoolpft
character(len=clen) :: pftname

character(len=16) :: trash
logical :: iscomment1, iscomment2

!...misc variables
integer(i4) :: i, iref, ref
integer(i4) :: lp,frp,wp,pp
integer(i4) :: lpc13,frpc13,wpc13,ppc13
integer(i4) :: lpc14,frpc14,wpc14,ppc14

!...alias the pool indices
lp =  pool_indx_leaf
frp = pool_indx_froot
wp =  pool_indx_stwd
pp =  pool_indx_prod

lpc13 =  pool_indx_leaf_c13-6
frpc13 = pool_indx_froot_c13-6
wpc13 =  pool_indx_stwd_c13-6
ppc13 =  pool_indx_prod_c13-6

lpc14 =  pool_indx_leaf_c14-12
frpc14 = pool_indx_froot_c14-12
wpc14 =  pool_indx_stwd_c14-12
ppc14 =  pool_indx_prod_c14-12

!...Initialize parameters to zero.
allocate(phencon(npft))
call init_phencon(npft, npoolpft, npstgmax, phencon)

!...Open file
print*,'Reading GDD-Based Phenology File: '
print*,'  ',trim(pgdd_file)
open(unit=pgddid,file=trim(pgdd_file),form='formatted')

read(pgddid,*) finpft
read(pgddid,*) finpftgdd
read(pgddid,*) finpoolpft

!.....check file vs simulation constants
if (finpft /= npft) then
    print*,''
    print('(a)'),'!!!GDD phenological file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Phen file npft: ',finpft,' Sim npft: ',npft
    print*,''
    stop
endif

if (finpftgdd /= npft_gdd) then
    print*,''
    print('(a)'),'!!!GDD phenological file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Phen file npft_gdd: ',finpftgdd,  &
          ' Sim npft_gdd: ',npft_gdd
    print*,''
    stop
endif

if (finpoolpft /= npoolpft) then
    print*,''
    print('(a)'),'!!!GDD phenological file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Phen file npoolpft: ',finpoolpft,  &
          ' Sim npoolpft: ',npoolpft
    print*,''
    stop
endif

!...Read in comments
iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(pgddid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else 
           iscomment2=.false.
        endif
    endif
enddo


!...Read in growing season start parameters
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, &
        phencon(iref)%daylen_mini, phencon(iref)%daylen_offd, &
        phencon(iref)%precip_len, phencon(iref)%precip_bef, &
        phencon(iref)%precip_aft, &
        phencon(iref)%tawftop_min, phencon(iref)%tawftop_len, &
        phencon(iref)%tm_min, phencon(iref)%tm_max, phencon(iref)%tm_len
enddo

!...Read in growing season reset parameters
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, &
        phencon(iref)%assim_resetl, phencon(iref)%assim_resetv
enddo

!...Read in stage and seed information
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%gdd_or_pd, &
         phencon(iref)%gdd_tbase, phencon(iref)%gdd_tmax, &
         phencon(iref)%npstg, phencon(iref)%seed_carbon, &
         phencon(iref)%seed_release, phencon(iref)%gslmax
enddo

!...Read in stage thresholds
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%threshp
enddo

!...Read in allocation fractions
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%allocp(lp,:), &
       phencon(iref)%adj_moist, phencon(iref)%adj_temp
!   backspace (pgddid)
!   read(pgddid,*) ref, pftname, phencon(iref)%allocp(lpc13,:), &
!       phencon(iref)%adj_moist, phencon(iref)%adj_temp
   phencon(iref)%allocp(lpc13,:) = &
               phencon(iref)%allocp(lp,:)
   phencon(iref)%allocp(lpc14,:) = &
               phencon(iref)%allocp(lp,:)
enddo

read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%allocp(frp,:)
   phencon(iref)%allocp(frpc13,:) = &
               phencon(iref)%allocp(frp,:)
   phencon(iref)%allocp(frpc14,:) = &
               phencon(iref)%allocp(frp,:)
enddo

read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%allocp(wp,:)
   phencon(iref)%allocp(wpc13,:) = &
               phencon(iref)%allocp(wp,:)
   phencon(iref)%allocp(wpc14,:) = &
               phencon(iref)%allocp(wp,:)
enddo

read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%allocp(pp,:)
   phencon(iref)%allocp(ppc13,:) = &
               phencon(iref)%allocp(pp,:)
   phencon(iref)%allocp(ppc14,:) = &
               phencon(iref)%allocp(pp,:)
enddo

!...Read in leaf transfer fractions
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, phencon(iref)%lptransfer(:)
enddo

!...Read in the vmax
read(pgddid,*) trash
read(pgddid,*) trash
read(pgddid,*) trash
do i=1, npft_gdd
   iref = pftnum_gddindx(i)
   read(pgddid,*) ref, pftname, &
        phencon(iref)%vmax0(:)
enddo

close(pgddid)

!------Perform conversions and calculate dependent parameters----
do i=1, npft_gdd
    iref = pftnum_gddindx(i)

   !...convert seed information from grams to moles
   phencon(iref)%seed_carbon = phencon(iref)%seed_carbon / real(mwc)
   phencon(iref)%seed_release = phencon(iref)%seed_release / real(mwc)

   !...set averaging weights
   phencon(iref)%wt_assim = 1./MAX(1.0,real(phencon(iref)%assim_resetl)*real(steps_per_day))
   phencon(iref)%wt_precip = 1./MAX(1.0, real(phencon(iref)%precip_len)*real(steps_per_day))
   phencon(iref)%wt_tawftop = 1./MAX(1.0,real(phencon(iref)%tawftop_len)*real(steps_per_day))
   phencon(iref)%wt_tm = 1./MAX(1.0,real(phencon(iref)%tm_len)*real(steps_per_day))

enddo


end subroutine read_pgdd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_phencon(npft, npoolpft, npstgmax, &
      phencon)

  use kinds
  use module_param, only: &
      phen_param
  use module_time, only: wt_daily

  implicit none

  !...input variables
  integer(i4) :: npft, npoolpft, npstgmax
  type(phen_param), dimension(npft), intent(inout) :: phencon

  !...local variables
  integer(i4) :: iref

  phencon%tawftop_len = izero
  phencon%tm_len = izero
  phencon%daylen_mini = rzero
  phencon%daylen_offd = rzero
  phencon%precip_len = izero
  phencon%precip_bef = izero
  phencon%precip_aft = izero
  phencon%tawftop_min = rzero
  phencon%tm_min = rzero
  phencon%tm_max = rzero

  phencon%assim_resetl = izero
  phencon%assim_resetv = rzero

  phencon%npstg = izero

  do iref=1,npft
     allocate(phencon(iref)%allocp(npoolpft,npstgmax))
     phencon(iref)%allocp = rzero

     allocate(phencon(iref)%threshp(npstgmax-1))
     phencon(iref)%threshp = rzero

     allocate(phencon(iref)%lptransfer(npstgmax))
     phencon(iref)%lptransfer = rzero

     allocate(phencon(iref)%vmax0(npstgmax))
     phencon(iref)%vmax0 = 0.1E-4
  enddo
  phencon%adj_moist = .false.
  phencon%adj_temp = .false.
  
  phencon%psdayl_ref = rzero
  phencon%psdayl_mul = rzero
  phencon%psdayl_min = rzero

  phencon%clai_coef = rzero
  phencon%clai_offl = rzero
  phencon%clai_offg = rzero
  phencon%climp_a = rzero
  phencon%climp_b = rzero
  phencon%climp_c = rzero
  phencon%climp_d = rzero
  phencon%climp_min = rzero
  phencon%climp_max = rzero
  phencon%cwa_type = izero
  phencon%psg_min = rzero

  phencon%pswx_rml = ione
  phencon%pswx_type = ione

  phencon%gdd_or_pd = bzero
  phencon%gdd_tbase = rzero
  phencon%gdd_tmax = rzero

  phencon%gslmax = izero
  phencon%seed_carbon = rzero
  phencon%seed_release = rzero

  phencon%wt_assim = wt_daily
  phencon%wt_precip = wt_daily
  phencon%wt_tawftop = wt_daily
  phencon%wt_tm = wt_daily
  phencon%wt_pswx = wt_daily
  
end subroutine init_phencon
