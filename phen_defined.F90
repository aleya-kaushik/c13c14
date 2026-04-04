!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Defined Phenology Based on User-Defined Stage Lengths
!   Number Of Stages = npstg 
!   1:            Planting (prior to emergence)
!   2:            Emergence (transfer of carbon from seed)
!   3-??:         Growth and Development
!   ??-(npstg-1): Browning and Drying
!   npstg:        Harvest/Dormant
!
!   Stages are determined based on either:
!    - Growing Degree Days (GDD)
!    - Days After Planting Date (DAPD)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine phen_defined( &
    subpt, subl, ipft, pnum, &
    phencont, poolcont, &
    lai, tmdf, phent, pooldt, &
    poollt, physcont)

use kinds
use module_pftinfo, only: &
    pft_c4c, pft_mze, pft_soy, &
    pft_wwt, pft_ref
use module_poolinfo, only: &
    pool_indx_lay, pool_name
use module_pparams, only: &
    mwc, month_names, tffrz, &
    pdb, stdC14
use module_io, only: &
    npbp, pbp_gref, pbp_pref
use module_param, only: &
    phen_param, pool_param, &
    phys_param
use module_sib, only: &
    phen_type, poold_type, &
    pooll_type
use module_sibconst, only: &
    cornsoy_switch, &
    print_harvest, print_stop, &
    npoolpft, ntpool, &
    subset, sublonsib, sublatsib, subpref, &
    nisodatayr, varcisom_switch, &
    varciso_switch, varco2_switch
use module_time, only: &
    doy, day, month, year, &
    steps_per_day, dtisib, dtsib, &
    startyear
use module_isodata, only: &
    isoyr, globc13, globc14
!use module_phosib, only: c4

implicit none

!...parameters
integer(i4), parameter :: igrowstg=2

!...input variables
integer(i4), intent(in) :: subpt, subl, pnum
integer(i4), intent(inout) :: ipft
type(phen_param), intent(in) :: phencont
type(pool_param), intent(in) :: poolcont
real(r8), intent(in) :: lai, tmdf
type(phen_type), intent(inout) :: phent
type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(phys_param), intent(in) :: physcont

!...local variables
integer(i4) :: i, k, m, mref, p, kref, tcref
integer(i4) :: ips, icount
integer(i4) :: yrnow, loc
real(r8) :: deltac, hrvstc, hrvstc13, tempc, tempc13
real(r8) :: hrvstt, hrvsttc13, ratio, ratioc13
real(r8) :: tempc14, ratioc14, hrvstc14, hrvsttc14
real(r8) :: c4_flag

real(r8) :: d_13cm, d_14cm, d_13cca, d_14cca, N
real(r8) :: r_c13a, r_c13assim, r_c13poolinitc3, r_c13poolinitc4 
real(r8) :: r_c14a, r_c14assim, r_c14poolinitc3, r_c14poolinitc4

!------set C4 flag------------
c4_flag = dble(physcont%c4flag)

!------initialize poolinit for c13,c14--------

if (varcisom_switch) then ! update the c13,c14 values from input file
   !..Update d13cm,d14cm from isodata input
   yrnow=year
   !idx=findloc(isoyrtmp,yrnow+0.5)
   do i=1,nisodatayr
     if (floor(isoyr(i)) .eq. yrnow) then
      loc=i
      exit
     endif
   enddo
  d_13cm = dble(globc13(loc))
  d_14cm = dble(globc14(loc))
else ! use the startyear to find loc
   do i=1,nisodatayr
     if (floor(isoyr(i)) .eq. startyear) then
      loc=i
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
   do i=1,nisodatayr
     if (floor(isoyr(i)) .eq. startyear) then
      loc=i
      exit
     endif
   enddo
  d_13cca = dble(globc13(loc))
  d_14cca = dble(globc14(loc))
endif

r_c13a = ((d_13cca/1000.0D0) + 1.0D0)*pdb
N = dble((1+dble(-0.025))**2.0) / (dble(1.0D0+d_13cca/1000.0D0)**2.0)
r_c14a = (dble(1.0D0 + d_14cca/1000.0D0)*stdC14)/N

if (c4_flag .EQ. dzero) then !c3 plants
   r_c13assim = r_c13a*((-18.0D0/1000.0D0) + 1.0D0)
   r_c13poolinitc3 = (r_c13assim/(r_c13assim+1.0D0))
   r_c14assim = r_c14a*(1.0D0+(-18.0D0/1000.0D0))**2
   r_c14poolinitc3 = r_c14assim
else !c4 plants
   r_c13assim = r_c13a*((-4.4D0/1000.0D0) + 1.0D0)
   r_c13poolinitc4 = (r_c13assim/(r_c13assim+1.0D0))
   r_c14assim = r_c14a*(1.0D0+(-4.4D0/1000.0D0))**2
   r_c14poolinitc4 = r_c14assim
endif

   !print*,'r_c13poolinitc3 from phen_defined :',r_c13poolinitc3
   !print*,'r_c13poolinitc4 from phen_defined :',r_c13poolinitc4
   !print*,'r_c14poolinitc3 from phen_defined :',r_c14poolinitc3
   !print*,'r_c14poolinitc4 from phen_defined :',r_c14poolinitc4

!------------------------------
!Reset variables
poollt%gain_seed = dzero
pooldt%gain_hrvst_lay = dzero
poollt%loss_hrvst_lay = dzero
poollt%resp_hrvst = dzero
poollt%resp_hrvstc13 = dzero
poollt%resp_hrvstc14 = dzero
poollt%rmvd_hrvst = dzero
poollt%rmvd_hrvstc13 = dzero
poollt%rmvd_hrvstc14 = dzero
phent%phen_istage = phencont%npstg

if ((phent%phenflag_gsspass) .and. &
    (phent%phenflag_precip) .and. &
    (phent%ipd .eq. izero)) then
    phent%ipd = doy
    phent%seed_pool = phencont%seed_carbon
endif

!Only set phenology stage index when planted
IF (phent%ipd .le. izero) RETURN

!Increment days after planting date
phent%dapd = phent%dapd + ione
IF (tmdf .gt. phencont%gdd_tbase) &
!IF (tmdf .gt. tffrz) &
   phent%dapdaf = phent%dapdaf + ione

!Increment growing degree days
IF (tmdf .ge. phencont%gdd_tmax) THEN
   phent%gdd = phent%gdd + (phencont%gdd_tmax - phencont%gdd_tbase)
ELSEIF (tmdf .ge. phencont%gdd_tbase) THEN
    phent%gdd = phent%gdd + (tmdf - phencont%gdd_tbase)
ENDIF

!Choose between using GDD or DAPD for phenology
if (phencont%gdd_or_pd .eq. 1) then
    phent%phen_pi = phent%gdd
elseif (phencont%gdd_or_pd .eq. 2) then
    phent%phen_pi = phent%dapd
elseif (phencont%gdd_or_pd .eq. 3) then
    phent%phen_pi = phent%dapdaf
else
    print*,'Unexpected Metric For GDD Phenology: ',phencont%gdd_or_pd
    print*,'Expecting 1 (GDD) or 2 (DAPD).'
    print*,'Stopping.'
    stop
endif

!Set Phenology Stage
ips = izero
icount = 1
do while (ips .eq. izero)
   if (phent%phen_pi .lt. phencont%threshp(icount)) then
       ips = icount
   else
       icount = icount + ione
   endif

   if (icount .eq. phencont%npstg) then
       ips = icount
   endif
enddo
phent%phen_istage = ips

!Ensure Ample Emergence and Initial Growth Stages
IF ((pnum .eq. pft_mze) .or. (pnum .eq. pft_c4c)) THEN
   IF ((phent%phen_istage .GT. 2) .AND. &
        (phent%gdd .LT. phencont%threshp(3)) .AND. &
        (lai .LT. 0.8)) THEN
         phent%gdd = phencont%threshp(2)
         phent%phen_pi = phent%gdd
         phent%phen_istage = 2
    ELSEIF ((phent%phen_istage .GT. 3) .AND. &
         (phent%gdd .LT. phencont%threshp(4)) .AND. &
         (lai .LT. 2.2)) THEN
         phent%gdd = phencont%threshp(3)
         phent%phen_pi = phent%gdd
         phent%phen_istage = 3
    ENDIF
ENDIF

IF (pnum .eq. pft_wwt) THEN
   IF ((phent%phen_istage .GT. 3) .AND. &
       (phent%dapdaf .LT. phencont%threshp(4)) .AND. &
       (lai .LT. 1.8)) THEN
      phent%dapdaf = phent%dapdaf - ione
      phent%phen_pi = phent%dapdaf
      phent%phen_istage = 3
   ENDIF
ENDIF
     

!--------Specialty Pool Transfers---------!
!----Harvest----
if ((phent%phen_pi .ge. phencont%threshp(phencont%npstg-1)) .or. &
    (phent%dapd .ge. phencont%gslmax)) then
    !...calculate harvested carbon and
    !.....remove from pools
    hrvstc = dzero
    do p=1,npoolpft/3 !1,5 npoolpft (npoolpft=15 with C14)
       do k=1,pool_indx_lay(p) !pool_indx_lay(1,5) either 1 or 10
          tempc = poollt%poolpft_lay(p,k) &
               - poolcont%poolpft_min(p) &
               + poollt%poolpft_dgain(p,k) &
               - poollt%poolpft_dloss(p,k)
          hrvstc = hrvstc + tempc
          poollt%loss_hrvst_lay(p,k) = tempc
        enddo
     enddo
    !... same as above but for C-13 pools
    hrvstc13 = dzero
    do p=npoolpft/3+1,2*(npoolpft/3) !6,10 npoolpft
       tcref=p-npoolpft/3 !(p-5) index the original totC pools and apply rc factor to convert to c13
       kref=p+npoolpft/3+1 !12,16 ntpool pool_indx_lay(12,16) 1 or 10
       do k=1,pool_indx_lay(kref)
          tempc = poollt%poolpft_lay(tcref,k) &
               - poolcont%poolpft_min(tcref) &
               + poollt%poolpft_dgain(tcref,k) &
               - poollt%poolpft_dloss(tcref,k)
          !tempc13 = fract%rcpoolfac*tempc
          tempc13 = poollt%rcpoolpft_lay(p,k)*tempc
          hrvstc13 = hrvstc13 + tempc13
          poollt%loss_hrvst_lay(p,k) = tempc13

          !if ((poollt%rcpoolpft_lay(p,k) .gt. 1.)) then
           if ( (poollt%poolpft_dloss(p,k) .gt. 10.) .or. &
                (poollt%poolpft_dloss(p,k) .lt. -10.) ) then
               print*,' '
               print*,'code: phen_defined'
               print*,'p,k: ',p,k
               print*,'poolpft_dloss(p,k-1/k/k+1):',&
                   poollt%poolpft_dloss(p,k-1),poollt%poolpft_dloss(p,k),poollt%poolpft_dloss(p,k+1)
               print*,'poolpft_lay(p,k-1/k/k+1) :',&
                   poollt%poolpft_lay(p,k-1),poollt%poolpft_lay(p,k),poollt%poolpft_lay(p,k+1)
               print*,' '
           endif

       enddo
    enddo
    !... same as above but for C-14 pools
    hrvstc14 = dzero
    do p=2*(npoolpft/3)+1,npoolpft !11,15 npoolpft
       tcref=p-(2*npoolpft/3) !(p-10)
       kref=p+2*npoolpft/3+2 !23,27 ntpool pool_indx_lay(23,27) 1 or 10
       do k=1,pool_indx_lay(kref)
          tempc = poollt%poolpft_lay(tcref,k) &
               - poolcont%poolpft_min(tcref) &
               + poollt%poolpft_dgain(tcref,k) &
               - poollt%poolpft_dloss(tcref,k)
          !tempc13 = fract%rcpoolfac*tempc
          tempc14 = poollt%rcpoolpft_lay(p,k)*tempc
          hrvstc14 = hrvstc14 + tempc14
          poollt%loss_hrvst_lay(p,k) = tempc14
       enddo
    enddo

     poollt%poolpft_dloss = poollt%poolpft_dloss + &
          poollt%loss_hrvst_lay
     poollt%loss_hrvst_lay = poollt%loss_hrvst_lay*dtisib

     !...set carbon respired from harvest
     poollt%resp_hrvst = hrvstc * poolcont%harvest_trans(1) * dtisib
     poollt%resp_hrvstc13 = hrvstc13 * poolcont%harvest_trans(1) * dtisib
     poollt%resp_hrvstc14 = hrvstc14 * poolcont%harvest_trans(1) * dtisib

     !...set carbon removed from harvest
     poollt%rmvd_hrvst = hrvstc * poolcont%harvest_trans(2)
     poollt%rmvd_hrvstc13 = hrvstc13 * poolcont%harvest_trans(2)
     poollt%rmvd_hrvstc14 = hrvstc14 * poolcont%harvest_trans(2)

!print*,'phen_defined poollt%rcpoolpft_lay 6: ',poollt%rcpoolpft_lay(6,:)
!print*,'phen_defined poollt%rcpoolpft_lay 11: ',poollt%rcpoolpft_lay(11,:)
!
!print*,'phen_defined poollt%rcpoolpft_lay 7: ',poollt%rcpoolpft_lay(7,:)
!print*,'phen_defined poollt%rcpoolpft_lay 12: ',poollt%rcpoolpft_lay(12,:)
!
!print*,'phen_defined poollt%rcpoolpft_lay 10: ',poollt%rcpoolpft_lay(10,:)
!print*,'phen_defined poollt%rcpoolpft_lay 15: ',poollt%rcpoolpft_lay(15,:)
!
!print*,'phen_defined ratio c13 loss_hrvst_lay 6:',poollt%loss_hrvst_lay(6,:)/poollt%loss_hrvst_lay(1,:)
!print*,'phen_defined ratio c14 loss_hrvst_lay 11:',poollt%loss_hrvst_lay(11,:)/poollt%loss_hrvst_lay(1,:)
!
!print*,'phen_defined ratio c13 loss_hrvst_lay 7:',poollt%loss_hrvst_lay(7,:)/poollt%loss_hrvst_lay(2,:)
!print*,'phen_defined ratio c14 loss_hrvst_lay 12:',poollt%loss_hrvst_lay(12,:)/poollt%loss_hrvst_lay(2,:)
!
!print*,'phen_defined ratio c13 loss_hrvst_lay 8:',poollt%loss_hrvst_lay(8,:)/poollt%loss_hrvst_lay(3,:)
!print*,'phen_defined ratio c14 loss_hrvst_lay 13:',poollt%loss_hrvst_lay(13,:)/poollt%loss_hrvst_lay(3,:)
!
!print*,'phen_defined ratio c13 loss_hrvst_lay 9:',poollt%loss_hrvst_lay(9,:)/poollt%loss_hrvst_lay(4,:)
!print*,'phen_defined ratio c14 loss_hrvst_lay 14:',poollt%loss_hrvst_lay(14,:)/poollt%loss_hrvst_lay(4,:)
!
!print*,'phen_defined ratio c13 loss_hrvst_lay 10:',poollt%loss_hrvst_lay(10,:)/poollt%loss_hrvst_lay(5,:)
!print*,'phen_defined ratio c14 loss_hrvst_lay 15:',poollt%loss_hrvst_lay(15,:)/poollt%loss_hrvst_lay(5,:)



     !...transfer harvested carbon to dead pools
     do m=npoolpft/3+1,ntpool/3 ! goes from (6,11) ntpool to 1-6 dead pools
        mref = m - npoolpft/3 ! (m-5, i.e. 1-6 dead pools, with npoolpft=15)
        if (poolcont%harvest_trans(mref+2) .gt. rzero) then !3,8
           do k=1,pool_indx_lay(m) !6,11 ntpool
              pooldt%gain_hrvst_lay(mref,k) = hrvstc &
                        * poolcont%harvest_trans(mref+2) &
                        * pooldt%poollu_flay(mref,k)
              pooldt%poollu_dgain(mref,k) = pooldt%poollu_dgain(mref,k) &
                        + pooldt%gain_hrvst_lay(mref,k)
           enddo
        endif
     enddo
     ! same as above but for C13
     do m=npoolpft+2,2*ntpool/3 ! goes from (17,22) ntpool to 7-12 dead pools
        mref = m - 2*npoolpft/3 ! (m-10, 7-12 dead pools, with npoolpft=15)
        if (poolcont%harvest_trans(mref+2) .gt. rzero) then !9,14
           do k=1,pool_indx_lay(m) !17,22 ntpool
              pooldt%gain_hrvst_lay(mref,k) = hrvstc13 &
                        * poolcont%harvest_trans(mref+2) &
                        * pooldt%poollu_flay(mref,k)
              pooldt%poollu_dgain(mref,k) = pooldt%poollu_dgain(mref,k) &
                        + pooldt%gain_hrvst_lay(mref,k)
           enddo
        endif
     enddo
     ! same as above but for C14
     do m=ntpool-npoolpft/3,ntpool ! goes from (28,33) ntpool to 13-18 dead pools
        mref = m - npoolpft ! (m-15, 13-18 dead pools, with npoolpft=15)
        if (poolcont%harvest_trans(mref+2) .gt. rzero) then !15,20
           do k=1,pool_indx_lay(m) !28,33 ntpool
              pooldt%gain_hrvst_lay(mref,k) = hrvstc14 &
                        * poolcont%harvest_trans(mref+2) &
                        * pooldt%poollu_flay(mref,k)
              pooldt%poollu_dgain(mref,k) = pooldt%poollu_dgain(mref,k) &
                        + pooldt%gain_hrvst_lay(mref,k)
           enddo
        endif
     enddo


     hrvstt = hrvstc * sum(poolcont%harvest_trans(3:8))
     if (sum(pooldt%gain_hrvst_lay(1:6,:)) .gt. dzero) then
       ratio = hrvstt/sum(pooldt%gain_hrvst_lay(1:6,:))
     endif
     hrvsttc13 = hrvstc13 * sum(poolcont%harvest_trans(9:14))
     if (sum(pooldt%gain_hrvst_lay(7:12,:)) .gt. dzero) then
       ratioc13 = hrvsttc13/sum(pooldt%gain_hrvst_lay(7:12,:))
     endif
     hrvsttc14 = hrvstc14 * sum(poolcont%harvest_trans(15:20))
     if (sum(pooldt%gain_hrvst_lay(13:18,:)) .gt. dzero) then
       ratioc14 = hrvsttc14/sum(pooldt%gain_hrvst_lay(13:18,:))
     endif


     pooldt%gain_hrvst_lay(1:6,:) = pooldt%gain_hrvst_lay(1:6,:)*ratio      
     pooldt%gain_hrvst_lay(1:6,:) = pooldt%gain_hrvst_lay(1:6,:)*dtisib
     pooldt%gain_hrvst_lay(7:12,:) = pooldt%gain_hrvst_lay(7:12,:)*ratioc13
     pooldt%gain_hrvst_lay(7:12,:) = pooldt%gain_hrvst_lay(7:12,:)*dtisib
     pooldt%gain_hrvst_lay(13:18,:) = pooldt%gain_hrvst_lay(13:18,:)*ratioc14
     pooldt%gain_hrvst_lay(13:18,:) = pooldt%gain_hrvst_lay(13:18,:)*dtisib


!print*,'phen_defined ratio: ',ratio
!print*,'phen_defined ratioc13: ',ratioc13
!print*,'phen_defined ratioc14: ',ratioc14
!
!print*,'phen_defined ratio c13 gain_hrvst_lay 7:',pooldt%gain_hrvst_lay(7,:)/pooldt%gain_hrvst_lay(1,:)
!print*,'phen_defined ratio c14 gain_hrvst_lay 13:',pooldt%gain_hrvst_lay(13,:)/pooldt%gain_hrvst_lay(1,:)
!
!print*,'phen_defined ratio c13 gain_hrvst_lay 8:',pooldt%gain_hrvst_lay(8,:)/pooldt%gain_hrvst_lay(2,:)
!print*,'phen_defined ratio c14 gain_hrvst_lay 14:',pooldt%gain_hrvst_lay(14,:)/pooldt%gain_hrvst_lay(2,:)
!
!print*,'phen_defined ratio c13 gain_hrvst_lay 9:',pooldt%gain_hrvst_lay(9,:)/pooldt%gain_hrvst_lay(3,:)
!print*,'phen_defined ratio c14 gain_hrvst_lay 15:',pooldt%gain_hrvst_lay(15,:)/pooldt%gain_hrvst_lay(3,:)
!
!print*,'phen_defined ratio c13 gain_hrvst_lay 10:',pooldt%gain_hrvst_lay(10,:)/pooldt%gain_hrvst_lay(4,:)
!print*,'phen_defined ratio c14 gain_hrvst_lay 16:',pooldt%gain_hrvst_lay(16,:)/pooldt%gain_hrvst_lay(4,:)
!
!print*,'phen_defined ratio c13 gain_hrvst_lay 11:',pooldt%gain_hrvst_lay(11,:)/pooldt%gain_hrvst_lay(5,:)
!print*,'phen_defined ratio c14 gain_hrvst_lay 17:',pooldt%gain_hrvst_lay(17,:)/pooldt%gain_hrvst_lay(5,:)
!
!print*,'phen_defined ratio c13 gain_hrvst_lay 12:',pooldt%gain_hrvst_lay(12,:)/pooldt%gain_hrvst_lay(6,:)
!print*,'phen_defined ratio c14 gain_hrvst_lay 18:',pooldt%gain_hrvst_lay(18,:)/pooldt%gain_hrvst_lay(6,:)

     !...reset growing season variables
     phent%nd_gs = dzero
     phent%nd_stg(:) = dzero

     phent%ipd = izero
     phent%dapd = izero
     phent%dapdaf = izero
     phent%gdd = dzero
     phent%phen_pi = dzero

     !...switch crops if requested
     if (cornsoy_switch) then
        IF (pnum == pft_mze) THEN
           ipft = pft_ref(pft_soy)
           subpref(subpt,subl) = pft_ref(pft_soy)

           do i=1,npbp
              if (pbp_gref(i) .eq. subpt) then
                 pbp_pref(i,subl) = pft_ref(pft_soy)
             endif
           enddo
        ELSEIF (pnum == pft_soy) THEN
           ipft = pft_ref(pft_mze)
           subpref(subpt,subl) = pft_ref(pft_mze)
           do i=1,npbp
              if (pbp_gref(i) .eq. subpt) then
                 pbp_pref(i,subl) = pft_ref(pft_mze)
              endif
           enddo
        ENDIF !mze/soy pfts
    endif !switch corn and soybeans

    !--------Print Info---------!
    if (print_harvest) then
       print('(a)'), '---Harvesting---'
       print('(a,a,i3,a,i4)'), &
             'DATE: ', trim(month_names(month)), day, ', ', year
       print('(a,i6,2f8.3,i3)'), &
           '   SiB Point/Lon/Lat/PFT: ', &
            subset(subpt), sublonsib(subpt), &
            sublatsib(subpt), pft_ref(pnum)
       print('(a,f14.4)'), '  Carbon Harvested (g C/m2): ', &
          hrvstc*mwc
       do p=1,npoolpft/3
          print('(3a,f7.3)'),'      ',pool_name(p),': ', &
             sum(poollt%loss_hrvst_lay(p,:))*mwc*dtsib
       enddo

       print('(a,f14.4)'), '  Carbon Moved (g C/m2): ', &
          (poollt%rmvd_hrvst + poollt%resp_hrvst*dtsib &
            + sum(pooldt%gain_hrvst_lay(1:6,:)))*mwc*dtsib
       print('(a,f14.4)'), '      Respired   : ',poollt%resp_hrvst*dtsib*mwc
       print('(a,f14.4)'), '      Removed    : ',poollt%rmvd_hrvst*mwc
       print('(a,f14.4)'), '      Transferred: ', &
            sum(pooldt%gain_hrvst_lay(1:6,:))*mwc*dtsib
   
   if (print_stop) stop
endif

endif !harvest


!----Seed Growth----
!....beginning growth after emergence
deltac = dzero
IF ((phent%seed_pool .gt. dzero) .and. &
    (ips .ge. igrowstg)) then

    !...calculate carbon to add
    deltac = phencont%seed_release 
    if (phent%seed_pool .le. deltac) then
       deltac = phent%seed_pool
       phent%seed_pool = dzero
    else
       phent%seed_pool = phent%seed_pool - deltac
    endif

    !...add seed carbon to pools
    do p=1, npoolpft/3 !1,5 npoolpft
       poollt%gain_seed(p) = deltac*phencont%allocp(p,ips)
       do k=1,pool_indx_lay(p) !1,5 ntpool
          poollt%poolpft_dgain(p,k) = poollt%poolpft_dgain(p,k) &
              + poollt%gain_seed(p) * poollt%poolpft_flay(p,k)
       enddo
    enddo

    !...add seed carbon to C13 pools
    do p=npoolpft/3+1, 2*npoolpft/3 !6,10 npoolpft with npoolpft=15
       kref=p+npoolpft/3+1 !12,16 ntpool
       if (c4_flag .EQ. dzero) then !c3 plants
          if (poollt%rcpoolpft(p) .gt. dzero) then
             poollt%gain_seed(p) = poollt%rcpoolpft(p)*deltac*phencont%allocp(p,ips)
          else
             poollt%gain_seed(p) = r_c13poolinitc3*deltac*phencont%allocp(p,ips) !based on d13c=-26
          endif
       else !c4 plants
          if (poollt%rcpoolpft(p) .gt. dzero) then
             poollt%gain_seed(p) = poollt%rcpoolpft(p)*deltac*phencont%allocp(p,ips)
          else
             poollt%gain_seed(p) = r_c13poolinitc4*deltac*phencont%allocp(p,ips) !based on d13c=-12.4              
          endif
       endif

       do k=1,pool_indx_lay(kref) !12,16 ntpool
          poollt%poolpft_dgain(p,k) = poollt%poolpft_dgain(p,k) &
              + poollt%gain_seed(p) * poollt%poolpft_flay(p,k)
       enddo
    enddo    

    !...add seed carbon to C14 pools
    do p=2*npoolpft/3+1, npoolpft !11,15 npoolpft with npoolpft=15
       kref=p+2*npoolpft/3+2 !23,27 ntpool
       if (c4_flag .EQ. dzero) then !c3 plants
          if (poollt%rcpoolpft(p) .gt. dzero) then
             poollt%gain_seed(p) = poollt%rcpoolpft(p)*deltac*phencont%allocp(p,ips)
          else
             poollt%gain_seed(p) = r_c14poolinitc3*deltac*phencont%allocp(p,ips) !based on d14c=equiv to d13c above
          endif
       else !c4 plants
          if (poollt%rcpoolpft(p) .gt. dzero) then
             poollt%gain_seed(p) = poollt%rcpoolpft(p)*deltac*phencont%allocp(p,ips)
          else
             poollt%gain_seed(p) = r_c14poolinitc4*deltac*phencont%allocp(p,ips) !based on d14c=equiv to d13c above              
          endif
       endif

       do k=1,pool_indx_lay(kref) !23,27 ntpool
          poollt%poolpft_dgain(p,k) = poollt%poolpft_dgain(p,k) &
              + poollt%gain_seed(p) * poollt%poolpft_flay(p,k)
       enddo
    enddo

    poollt%gain_seed = poollt%gain_seed / dble(steps_per_day) * dtisib
ENDIF

!print*,'phen_defined ratio c13 poollt%gain_seed 6: ', poollt%gain_seed(6)/poollt%gain_seed(1)
!print*,'phen_defined ratio c14 poollt%gain_seed 11: ', poollt%gain_seed(11)/poollt%gain_seed(1)
!
!print*,'phen_defined ratio c13 poollt%gain_seed 7: ', poollt%gain_seed(7)/poollt%gain_seed(2)
!print*,'phen_defined ratio c14 poollt%gain_seed 12: ', poollt%gain_seed(12)/poollt%gain_seed(2)
!
!print*,'phen_defined ratio c13 poollt%gain_seed 10: ', poollt%gain_seed(10)/poollt%gain_seed(5)
!print*,'phen_defined ratio c14 poollt%gain_seed 15: ', poollt%gain_seed(15)/poollt%gain_seed(5)


end subroutine phen_defined
