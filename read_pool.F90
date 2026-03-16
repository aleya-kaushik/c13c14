
! Opens and reads in SiB pool parameters.
subroutine read_pool()

use kinds
use module_sibconst, only: &
    cornsoy_switch, &
    npft, ngroup, &
    npoolpft, npoollu, ntpool, &
    nisodatayr, varcisom_switch, &
    varciso_switch, varco2_switch
use module_io, only: &
    pool_file, poolid
use module_pparams, only: &
    secs_per_day, days_per_year, &
    drytoc, mwc, pdb, stdC14
use module_pftinfo, only: &
    clen,npft_gdd, &
    pft_mze, pft_soy, &
    pft_group, group_grass, group_crop
use module_poolinfo, only: &
    pool_indx_leaf, &
    pool_indx_leaf_c13, &
    pool_indx_leaf_c14
use module_param, only: &
    physcon, poolcon
use module_time, only: &
    year, startyear
use module_isodata, only: &
    isoyr, globc13, globc14

!use module_phosib, only: c4

implicit none

!...file variables
integer(i4) :: finpft, fingddpft
integer(i4) :: finpool, fingroup
integer(i4) :: num
character(len=clen), dimension(npft) :: pftname
character(len=10), dimension(ntpool) :: poolname
character(len=100) :: trash
logical :: iscomment1, iscomment2

real(r4), dimension(npoollu+2) :: graze_trans, harvest_trans

!...misc variables
integer(i4) :: i,j,iiso
integer(byte) :: groupref
real(r8) :: poolval
integer(i4) :: lp, lpc13, lpc14
integer(i4) :: yrnow, loc
real(r8) :: d_13cm, d_14cm, d_13cca, d_14cca, N
real(r8) :: r_c13a, r_c13assim, r_c13poolinitc3, r_c13poolinitc4
real(r8) :: r_c14a, r_c14assim, r_c14poolinitc3, r_c14poolinitc4

!...alias the pool indices
lp =  pool_indx_leaf
lpc13 =  pool_indx_leaf_c13-6
lpc14 =  pool_indx_leaf_c14-12

!-------------------
!...Initialize the pool variables
allocate(poolcon(npft))
call init_poolcon(npft, npoolpft, npoollu, poolcon)

!-------------------

!...Open file
print*,'Reading Pool Parameter File: '
print*,'  ',trim(pool_file)
open(unit=poolid,file=trim(pool_file),form='formatted')

read(poolid,*) finpft
if (finpft /= npft) then
    print*,''
    print('(a)'),'!!!Pool File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Pool file npft: ', &
          finpft,' Sim npft: ',npft
    print*,''
    stop
endif

read(poolid,*) fingddpft
if (fingddpft /= npft_gdd) then
    print*,''
    print('(a)'),'!!!Pool File Mistmatch!!!'
    print('(a,i4,a,i4)'),'  Pool file npft_gdd: ', &
          fingddpft,' Sim npft_gdd: ',npft_gdd
    print*,''
    stop
endif
read(poolid,*) finpool
read(poolid,*) fingroup

if (finpool /= ntpool) then
    print*,''
    print('(a)'),'!!!Pool param file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Pool file npool: ',finpool,' Sim npool: ',ntpool
    print*,''
    stop
endif

if (fingroup /= ngroup) then
    print*,''
    print('(a)'),'!!!Pool param file does not match simulation!!!'
    print('(a,i4,a,i4)'),'  Pool file ngroup: ',fingroup,' Sim ngroup: ',ngroup
    print*,''
    stop
endif

iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(poolid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else 
           iscomment2=.false.
        endif
    endif
enddo

!...Read in the variables
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%gr_frac(1:5)
   poolcon(i)%gr_frac(6:10)=poolcon(i)%gr_frac(1:5)
   poolcon(i)%gr_frac(11:15)=poolcon(i)%gr_frac(1:5)
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%lresp_eff(1:5)
   poolcon(i)%lresp_eff(6:10)=poolcon(i)%lresp_eff(1:5)
   poolcon(i)%lresp_eff(11:15)=poolcon(i)%lresp_eff(1:5)
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%cr_aml,  poolcon(i)%cr_amh,  &
        poolcon(i)%cr_amin, poolcon(i)%cr_amax, &
        poolcon(i)%cr_fmul, poolcon(i)%cr_fref, &
        poolcon(i)%cr_fmin, &
        poolcon(i)%cr_hq10, poolcon(i)%cr_href, &
        poolcon(i)%cr_hmax
enddo


read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%lt_fref, poolcon(i)%lt_fq10, &
        poolcon(i)%lt_fmax, &
        poolcon(i)%lt_dref, poolcon(i)%lt_dcoef, &
        poolcon(i)%lt_dmax, &
        poolcon(i)%lt_wref, poolcon(i)%lt_wbase, &
        poolcon(i)%lt_wcoef, poolcon(i)%lt_wmax
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%rrt_aml, poolcon(i)%rrt_amh, &
        poolcon(i)%rrt_amin, poolcon(i)%rrt_amax, &
        poolcon(i)%rrt_fmul, poolcon(i)%rrt_fref, &
        poolcon(i)%rrt_fmin, &
        poolcon(i)%rrt_hq10, poolcon(i)%rrt_href, &
        poolcon(i)%rrt_hmax, poolcon(i)%rrt_laimin, &
        poolcon(i)%rrt_laimax
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%hrt_sfc_aml, poolcon(i)%hrt_sfc_amh, &
        poolcon(i)%hrt_sfc_amin, poolcon(i)%hrt_sfc_amax, &
        poolcon(i)%hrt_sfc_fmul, poolcon(i)%hrt_sfc_fref, &
        poolcon(i)%hrt_sfc_fmin, &
        poolcon(i)%hrt_sfc_hq10, poolcon(i)%hrt_sfc_href, &
        poolcon(i)%hrt_sfc_hmax, &
        poolcon(i)%hrt_sfc_pml, poolcon(i)%hrt_sfc_pmin
enddo

read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%hrt_soil_aml, poolcon(i)%hrt_soil_amh, &
        poolcon(i)%hrt_soil_amin, poolcon(i)%hrt_soil_amax, &
        poolcon(i)%hrt_soil_fmul, poolcon(i)%hrt_soil_fref, &
        poolcon(i)%hrt_soil_fmin, &
        poolcon(i)%hrt_soil_hq10, poolcon(i)%hrt_soil_href, &
        poolcon(i)%hrt_soil_hmax, &
        poolcon(i)%hrt_soil_mmin, poolcon(i)%hrt_soil_pawmin
enddo


read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1, npft
   read(poolid,*) num, pftname(i), &
        poolcon(i)%turnover(1:11)
   poolcon(i)%turnover(12:22) = poolcon(i)%turnover(1:11)
   poolcon(i)%turnover(23:33) = poolcon(i)%turnover(1:11)
enddo
        
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash

!.....Biological efficiencies and transfers
do i=1,npft
   read(poolid,*) trash
   read(poolid,*) trash
   read(poolid,*) trash
   read(poolid,*) trash
   do j=1, npoollu/3 !npoollu is 18, first 6 are C, next 6 are C-13, last 6 are C-14
     read(poolid,*) num, poolname(j), poolcon(i)%dresp_eff(1:6,j)
   enddo
   poolcon(i)%dresp_eff(7:12,7:12)=poolcon(i)%dresp_eff(1:6,1:6)
   poolcon(i)%dresp_eff(1:6,7:12)=dzero
   poolcon(i)%dresp_eff(7:12,1:6)=dzero

   poolcon(i)%dresp_eff(13:18,13:18)=poolcon(i)%dresp_eff(1:6,1:6)
   poolcon(i)%dresp_eff(1:6,13:18)=dzero
   poolcon(i)%dresp_eff(13:18,1:6)=dzero

   poolcon(i)%dresp_eff(7:12,13:18)=dzero
   poolcon(i)%dresp_eff(13:18,7:12)=dzero

   read(poolid,*) trash
   read(poolid,*) trash
   read(poolid,*) trash
   do j=1, ntpool/3 !ntpool is 33, first 11 are C, next 11 are C-13, last 11 are C-14            
     read(poolid,*) num, poolname(j), poolcon(i)%pool_trans_frac(1:11,j)
   enddo
   poolcon(i)%pool_trans_frac(12:22,12:22)=poolcon(i)%pool_trans_frac(1:11,1:11)
   poolcon(i)%pool_trans_frac(1:11,12:22)=dzero
   poolcon(i)%pool_trans_frac(12:22,1:11)=dzero

   poolcon(i)%pool_trans_frac(23:33,23:33)=poolcon(i)%pool_trans_frac(1:11,1:11)
   poolcon(i)%pool_trans_frac(1:11,23:33)=dzero
   poolcon(i)%pool_trans_frac(23:33,1:11)=dzero

   poolcon(i)%pool_trans_frac(12:22,23:33)=dzero
   poolcon(i)%pool_trans_frac(23:33,12:22)=dzero

enddo

!print*,'trans_frac (1,1:11): ',poolcon(2)%pool_trans_frac(1,1:11)
!print*,'trans_frac (12,12:22): ',poolcon(2)%pool_trans_frac(12,12:22)
!print*,'trans_frac (1,12:22): ',poolcon(2)%pool_trans_frac(1,12:22)
!print*,'trans_frac (12,1:11): ',poolcon(2)%pool_trans_frac(12,1:11)

!.....Specialty transfers
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
read(poolid,*) trash
do i=1,npoollu/3+2 ! 1,8 first 2 listings are respire and remove, next 6 are dead pools
    read(poolid,*) num, poolname(i), graze_trans(i), harvest_trans(i)
enddo
graze_trans(9:14) = graze_trans(3:8)
harvest_trans(9:14) = harvest_trans(3:8)

graze_trans(15:20) = graze_trans(3:8)
harvest_trans(15:20) = harvest_trans(3:8)

close(poolid)

!...Save the parameters
do i=1,npft
    groupref = pft_group(i)

   !...Set specialty transfers
   if (groupref .eq. group_grass) then
      poolcon(i)%graze_trans(:) = graze_trans(:)
   endif

   if (groupref .eq. group_crop) then
      if ( (sum(harvest_trans(1:8)) .gt. 0.999) .and. &
           (sum(harvest_trans(1:8)) .lt. 1.001) ) then
          poolcon(i)%harvest_trans(:) = harvest_trans(:)
      else
         print*,''
         print*,'--Incorrect Harvest Transfer Fractions--'
         print*,'Must Sum To 1, Current Sum: ',sum(harvest_trans(1:8))
         print*,'Stopping.'
         print*,''
         stop
     endif
   endif

   !...Set calculated parameters
   do j=1,ntpool
      if (poolcon(i)%turnover(j) .gt. 1.e-12) then
         poolcon(i)%k_rate(j) = 1./ &
           (poolcon(i)%turnover(j)*real(secs_per_day)*real(days_per_year))
     endif
   enddo
   
  !...set pool minimum values
  if (physcon(i)%sla .gt. 1.E-12) then
     poolval = dble(physcon(i)%laimin / physcon(i)%sla &
                / real(drytoc) / real(mwc))
     poolcon(i)%poolpft_min(lp) = dble(poolval)

       !...initialize the rpoolinit for c13 and c14
       
       if (varcisom_switch) then ! update the c13,c14 values from input file
          !..Update d13cm,d14cm from isodata input
          yrnow=year
          !idx=findloc(isoyrtmp,yrnow+0.5)
          do iiso=1,nisodatayr
            if (floor(isoyr(iiso)) .eq. yrnow) then
             loc=iiso
             exit
            endif
          enddo
         d_13cm = dble(globc13(loc))
         d_14cm = dble(globc14(loc))
       else ! use the startyear to find loc
          do iiso=1,nisodatayr
            if (floor(isoyr(iiso)) .eq. startyear) then
             loc=iiso
             exit
            endif
          enddo
         d_13cm = dble(globc13(loc))
         d_14cm = dble(globc14(loc))
       endif
       
       if (varciso_switch .or. varco2_switch) then
          d_13cca = d_13cm
          d_14cca = d_14cm
       else !use the startyear to find values
          do iiso=1,nisodatayr
            if (floor(isoyr(iiso)) .eq. startyear) then
             loc=iiso
             exit
            endif
          enddo
         d_13cca = dble(globc13(loc))
         d_14cca = dble(globc14(loc))
       endif
       
       r_c13a = ((d_13cca/1000.0D0) + 1.0D0)*pdb

       N = dble((1+dble(-0.025))**2.0) / (dble(1.0D0+d_13cca/1000.0D0)**2.0)
       r_c14a = (dble(d_14cca+1.0D0)*stdC14)/N

       if (physcon(i)%c4flag .EQ. dzero) then !c3 plants
          r_c13assim = r_c13a*((-18.0D0/1000.0D0) + 1.0D0)
          r_c13poolinitc3 = (r_c13assim/(r_c13assim+1.0D0))

          r_c14assim = r_c14a*((1.0D0+(-18.0D0/1000.0D0))**2)
          r_c14poolinitc3 = r_c14assim
          poolcon(i)%poolpft_min(lpc13) = dble(r_c13poolinitc3*poolval) ! based on rcassim equiv to -26
          poolcon(i)%poolpft_min(lpc14) = dble(r_c14poolinitc3*poolval)
       else !c4 plants
          r_c13assim = r_c13a*((-4.4D0/1000.0D0) + 1.0D0)
          r_c13poolinitc4 = (r_c13assim/(r_c13assim+1.0D0))

          r_c14assim = r_c14a*((1.0D0+(-4.4D0/1000.0D0))**2)
          r_c14poolinitc4 = r_c14assim
          poolcon(i)%poolpft_min(lpc13) = dble(r_c13poolinitc4*poolval) ! based on rcassim equiv to -12.4          
          poolcon(i)%poolpft_min(lpc14) = dble(r_c14poolinitc4*poolval)
       endif

   !print*,'lp min from readpool :',poolcon(i)%poolpft_min(lp)
   !print*,'lpc13 min from readpool :',poolcon(i)%poolpft_min(lpc13)

   endif !sla > 0

enddo

if (cornsoy_switch) then
   poolval = dble(max(poolcon(pft_mze)%poolpft_min(1), &
                 poolcon(pft_soy)%poolpft_min(1)))
   poolcon(pft_mze)%poolpft_min(lp) = dble(poolval)
   poolcon(pft_soy)%poolpft_min(lp) = dble(poolval)
   poolcon(pft_mze)%poolpft_min(lpc13) = dble(r_c13poolinitc4*poolval) ! same as above
   poolcon(pft_soy)%poolpft_min(lpc13) = dble(r_c13poolinitc3*poolval) ! same as above
   poolcon(pft_mze)%poolpft_min(lpc14) = dble(r_c14poolinitc4*poolval) ! same as above
   poolcon(pft_soy)%poolpft_min(lpc14) = dble(r_c14poolinitc3*poolval) ! same as above
endif


end subroutine read_pool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initializes the pool parameters.
subroutine init_poolcon(npft, npoolpft, npoollu, &
     poolcon)

  use kinds
  use module_param, only: pool_param

  !...input variables
  integer(i4), intent(in) :: npft, npoolpft, npoollu
  type(pool_param), dimension(npft), intent(inout) :: poolcon

  !...local variables
  integer(i4) :: i, ntpool

  !...set local variables
  ntpool = npoolpft + npoollu

  !...initialize pool parameters
  do i=1,npft
     allocate(poolcon(i)%lresp_eff(npoolpft))
     poolcon(i)%lresp_eff(:) = dzero
     allocate(poolcon(i)%gr_frac(npoolpft))
     poolcon(i)%gr_frac(:) = dzero

     poolcon(i)%cr_aml = dzero
     poolcon(i)%cr_amh = dzero
     poolcon(i)%cr_amin = done
     poolcon(i)%cr_amax = done
     poolcon(i)%cr_fmul = dzero
     poolcon(i)%cr_fref = dzero
     poolcon(i)%cr_fmin = done
     poolcon(i)%cr_hq10 = dzero
     poolcon(i)%cr_href = dzero
     poolcon(i)%cr_hmax = done

     poolcon(i)%lt_fq10 = dzero
     poolcon(i)%lt_fref = dzero
     poolcon(i)%lt_fmax = dzero
     poolcon(i)%lt_dcoef = dzero
     poolcon(i)%lt_dref = dzero
     poolcon(i)%lt_dmax = dzero
     poolcon(i)%lt_wref = dzero
     poolcon(i)%lt_wbase = dzero
     poolcon(i)%lt_wcoef = dzero
     poolcon(i)%lt_wmax = dzero

     poolcon(i)%rrt_aml = dzero
     poolcon(i)%rrt_amh = dzero
     poolcon(i)%rrt_amin = done
     poolcon(i)%rrt_amax = done
     poolcon(i)%rrt_fmul = dzero
     poolcon(i)%rrt_fref = dzero
     poolcon(i)%rrt_fmin = done
     poolcon(i)%rrt_hq10 = dzero
     poolcon(i)%rrt_href = dzero
     poolcon(i)%rrt_hmax = done
     poolcon(i)%rrt_laimin = done
     poolcon(i)%rrt_laimax = done
     
     poolcon(i)%hrt_sfc_aml = dzero
     poolcon(i)%hrt_sfc_amh = dzero
     poolcon(i)%hrt_sfc_amin = done
     poolcon(i)%hrt_sfc_amax = done
     poolcon(i)%hrt_sfc_fmul = done
     poolcon(i)%hrt_sfc_fref = done
     poolcon(i)%hrt_sfc_fmin = done
     poolcon(i)%hrt_sfc_hq10 = done
     poolcon(i)%hrt_sfc_href = done
     poolcon(i)%hrt_sfc_hmax = done
     poolcon(i)%hrt_sfc_pml = dzero
     poolcon(i)%hrt_sfc_pmin = done

     poolcon(i)%hrt_soil_aml = dzero
     poolcon(i)%hrt_soil_amh = dzero
     poolcon(i)%hrt_soil_amin = done
     poolcon(i)%hrt_soil_amax = done
     poolcon(i)%hrt_soil_fmul = done
     poolcon(i)%hrt_soil_fref = done
     poolcon(i)%hrt_soil_fmin = done
     poolcon(i)%hrt_soil_hq10 = done
     poolcon(i)%hrt_soil_href = done
     poolcon(i)%hrt_soil_hmax = done
     poolcon(i)%hrt_soil_mmin = done
     poolcon(i)%hrt_soil_pawmin = done

     allocate(poolcon(i)%dresp_eff(npoollu,npoollu))
     poolcon(i)%dresp_eff(:,:) = dzero
     allocate(poolcon(i)%graze_trans(npoollu+2))
     allocate(poolcon(i)%harvest_trans(npoollu+2))
     poolcon(i)%graze_trans(:) = dzero
     poolcon(i)%harvest_trans(:) = dzero

     allocate(poolcon(i)%turnover(ntpool))
     poolcon(i)%turnover(:) = dzero
     allocate(poolcon(i)%pool_trans_frac(ntpool,ntpool))
     poolcon(i)%pool_trans_frac(:,:) = dzero
     allocate(poolcon(i)%k_rate(ntpool))
     poolcon(i)%k_rate(:) = dzero

     allocate(poolcon(i)%poolpft_min(npoolpft))
     poolcon(i)%poolpft_min(:) = dzero

  enddo !i=1,npft

end subroutine init_poolcon
