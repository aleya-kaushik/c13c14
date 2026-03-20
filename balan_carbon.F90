!...calculation to determine if carbon
!...balance is maintained
subroutine balan_carbon( &
    sibpt, siblon, siblat, pref, &
    assimin, c13assimin, c14assimin, & 
    laiin, fparin, &
    pooldt, poollt, grz_transfer, co2t, fract)

use kinds

use module_poolinfo
use module_pparams, only: &
    mol_to_mg, mol_to_umol, &
    month_names, molc13_to_mg, &
    molc14_to_mg, pdb
use module_sib, only: &
    poold_type, pooll_type, &
    co2_type, fract_type
use module_sibconst, only: &
    npoolcan, npoollu, npoolpft, &
    carbonb_print, carbonb_stop, &
    carbonb_thresh, carbonb_threshc13, &
    npoolcanc13, npoolcanc14
use module_time, only: &
    dtsib, day, month, year

implicit none

!...input variables
integer(i4), intent(in) :: sibpt, pref
real(r4), intent(in) :: siblon, siblat
real(r8), intent(in) :: assimin, c13assimin, c14assimin
real(r8), intent(in) :: fparin, laiin
real(r4), dimension(npoollu+2), intent(in) :: grz_transfer

type(poold_type), intent(inout) :: pooldt
type(pooll_type), intent(inout) :: poollt
type(co2_type), intent(in) :: co2t
type(fract_type), intent(in) :: fract

!...switches
logical, parameter :: checkmol=.true.
real(r8) :: myconvert

!...assimilation variables
real(r8) :: assimtot
real(r8) :: assimtotc13
real(r8) :: assimtotc14

!...dead pool variables
real(r8), dimension(npoollu) :: dgain, dloss
real(r8), dimension(npoolpft) :: lgain, lloss

real(r8) :: dcarbonb, dcarbonbc13, dcarbonbc14
real(r8) :: dcarbonin, dcarbonout, dcarbons
real(r8) :: dcarboninc13, dcarbonoutc13, dcarbonsc13
real(r8) :: dcarboninc14, dcarbonoutc14, dcarbonsc14
real(r8) :: dgaingrz, dgainhrv, &
            dgaintd, dgaintl
real(r8) :: dgaingrzc13, dgainhrvc13, &
            dgaintdc13, dgaintlc13
real(r8) :: dgaingrzc14, dgainhrvc14, &
            dgaintdc14, dgaintlc14
real(r8) :: dlossr, dlosstd, dlossf, dlossrdd
real(r8) :: dlossrc13, dlosstdc13, dlossfc13, dlossrddc13
real(r8) :: dlossrc14, dlosstdc14, dlossfc14, dlossrddc14
real(r8) :: dpoolinit, dpoolend, dpoolendc, dpoolchange
real(r8) :: dpoolendc13,dpoolinitc13, &
            dpoolchangec13
real(r8) :: dpoolendc14,dpoolinitc14, &
            dpoolchangec14

!...flux variables
real(r8) :: fassim, fresp
real(r8) :: respleaf, resproot
real(r8) :: respauto, respgrow, respmntn
real(r8) :: resphet, respsoil
real(r8) :: fassimc13, frespc13
real(r8) :: respleafc13, resprootc13
real(r8) :: respautoc13, respgrowc13, respmntnc13
real(r8) :: resphetc13, respsoilc13
real(r8) :: fassimc14, frespc14
real(r8) :: respleafc14, resprootc14
real(r8) :: respautoc14, respgrowc14, respmntnc14
real(r8) :: resphetc14, respsoilc14

!...grazing variables
real(r8) :: grzcarbonb
real(r8) :: grzcarbonbc13
real(r8) :: grzcarbonbc14

!...harvest variables
real(r8) :: hrvcarbonb
real(r8) :: hrvresp, hrvrmvd
real(r8) :: hrvcarbonbc13
real(r8) :: hrvcarbonbc14
real(r8) :: hrvrespc13, hrvrmvdc13
real(r8) :: hrvrespc14, hrvrmvdc14

!...live carbon pool variables
real(r8) :: lcarbonb, lcarbonbc13, lcarbonbc14
real(r8) :: lcarbonin, lcarbonout, lcarbons
real(r8) :: lcarboninc13, lcarbonoutc13, lcarbonsc13
real(r8) :: lcarboninc14, lcarbonoutc14, lcarbonsc14
real(r8) :: lgaina, lgains
real(r8) :: lgainac13, lgainsc13
real(r8) :: lgainac14, lgainsc14
real(r8) :: llossfire, llossfirec13, llossfirec14
real(r8) :: llossrdd, llossrddc13, llossrddc14
real(r8) :: llossgrz, llossrgrz, llosstgrz
real(r8) :: llossgrzc13, llossrgrzc13, llosstgrzc13
real(r8) :: llossgrzc14, llossrgrzc14, llosstgrzc14
real(r8) :: llosshrv, llosshrvc13, llosshrvc14
real(r8) :: llossrg, llossrm, llossrnveg
real(r8) :: llossrgc13, llossrmc13, llossrnvegc13
real(r8) :: llossrgc14, llossrmc14, llossrnvegc14
real(r8) :: llosstd, llosstdc13, llosstdc14
real(r8) :: lpoolinit, lpoolend, lpoolendc, lpoolchange
real(r8) :: lpoolendc13, lpoolinitc13, &
            lpoolchangec13
real(r8) :: lpoolendc14, lpoolinitc14, &
            lpoolchangec14

real(r8) :: netcarbonb
real(r8) :: netcarbonbc13
real(r8) :: netcarbonbc14

real(r8) :: tcarbonb
real(r8) :: tlivedead, gdeadlive
real(r8) :: tcarbonbc13
real(r8) :: tlivedeadc13, gdeadlivec13
real(r8) :: tcarbonbc14
real(r8) :: tlivedeadc14, gdeadlivec14

!...net balance variables
logical :: cb_err, loc_printout
logical :: cb_errc13, loc_printoutc13
logical :: cb_errc14, loc_printoutc14

!...reference indices
integer(i4) :: lp, wp, crp, frp, pp
integer(i4) :: cdb, strl, metl, slit, slow, arm 
integer(i4) :: lpc13, wpc13, crpc13, frpc13, ppc13
integer(i4) :: cdbc13, strlc13, metlc13, slitc13, slowc13, armc13
integer(i4) :: lpc14, wpc14, crpc14, frpc14, ppc14
integer(i4) :: cdbc14, strlc14, metlc14, slitc14, slowc14, armc14

!...misc variables
integer(i4) :: cref, p, crefc13, crefc14
real(r8) :: tempb, tempbc13, tempbc14
logical :: seed_printout
!real(r8) :: carbonb_threshc13

!---------------------------------------------
!...set local variables
!carbonb_threshc13=carbonb_thresh/1000.
cb_err = .false.
cb_errc13=.false.

lp =  pool_indx_leaf !ntpool index 1
frp = pool_indx_froot !ntpool index 2
crp = pool_indx_croot !ntpool index 3
wp =  pool_indx_stwd !ntpool index 4
pp =  pool_indx_prod !ntpool index 5

lpc13 =  pool_indx_leaf_c13-npoollu/3 !ntpool index 12, npoolpft index 6
frpc13 = pool_indx_froot_c13-npoollu/3 !ntpool index 13, npoolpft index 7
crpc13 = pool_indx_croot_c13-npoollu/3 !ntpool index 14, npoolpft index 8
wpc13 =  pool_indx_stwd_c13-npoollu/3 !ntpool index 15, npoolpft index 9
ppc13 =  pool_indx_prod_c13-npoollu/3 !ntpool index 16, npoolpft index 10

lpc14 =  pool_indx_leaf_c14-2*npoollu/3 !ntpool index 23, npoolpft index 11
frpc14 = pool_indx_froot_c14-2*npoollu/3 !ntpool index 24, npoolpft index 12
crpc14 = pool_indx_croot_c14-2*npoollu/3 !ntpool index 25, npoolpft index 13
wpc14 =  pool_indx_stwd_c14-2*npoollu/3 !ntpool index 26, npoolpft index 14
ppc14 =  pool_indx_prod_c14-2*npoollu/3 !ntpool index 27, npoolpft index 15

cdb  = pool_indx_cdb-npoolpft/3 !ntpool index 6, npoollu index 1
metl = pool_indx_metl-npoolpft/3 !ntpool index 7, npoollu index 2
strl = pool_indx_strl-npoolpft/3 !ntpool index 8, npoollu index 3
slit = pool_indx_slit-npoolpft/3 !ntpool index 9, npoollu index 4
slow = pool_indx_slow-npoolpft/3 !ntpool index 10, npoollu index 5
arm  = pool_indx_arm-npoolpft/3 !ntpool index 11, npoollu index 6

cdbc13  = pool_indx_cdb_c13 - 2*npoolpft/3 !ntpool index 17, npoollu index 7
metlc13  = pool_indx_metl_c13 - 2*npoolpft/3 !ntpool index 18, npoollu index 8
strlc13  = pool_indx_strl_c13 - 2*npoolpft/3 !ntpool index 19, npoollu index 9
slitc13  = pool_indx_slit_c13 - 2*npoolpft/3 !ntpool index 20, npoollu index 10
slowc13 = pool_indx_slow_c13 - 2*npoolpft/3 !ntpool index 21, npoollu index 11
armc13  = pool_indx_arm_c13 - 2*npoolpft/3 !ntpool index 22, npoollu index 12

cdbc14  = pool_indx_cdb_c14 - npoolpft !ntpool index 28, npoollu index 13
metlc14  = pool_indx_metl_c14 - npoolpft !ntpool index 29, npoollu index 14
strlc14  = pool_indx_strl_c14 - npoolpft !ntpool index 30, npoollu index 15
slitc14  = pool_indx_slit_c14 - npoolpft !ntpool index 31, npoollu index 16
slowc14 = pool_indx_slow_c14 - npoolpft !ntpool index 32, npoollu index 17
armc14 = pool_indx_arm_c14 - npoolpft !ntpool index 33, npoollu index 18

if (checkmol) then
   myconvert = done/mol_to_umol
else
   myconvert = done
endif
seed_printout = .false.


!----- DAILY CARBON BALANCE-----
!...dead pools
dpoolendc = (sum(pooldt%poollup(1:6)) + sum(pooldt%poollu_dgain(1:6,:)) &
     - sum(pooldt%poollu_dloss(1:6,:)))*mol_to_mg
tempb = sum(pooldt%poollu(1:6))*mol_to_mg - dpoolendc
dpoolendc13 = (sum(pooldt%poollup(7:12)) + sum(pooldt%poollu_dgain(7:12,:)) &
     - sum(pooldt%poollu_dloss(7:12,:)))*molc13_to_mg
tempbc13 = sum(pooldt%poollu(7:12))*molc13_to_mg - dpoolendc13
dpoolendc14 = (sum(pooldt%poollup(13:18)) + sum(pooldt%poollu_dgain(13:18,:)) &
     - sum(pooldt%poollu_dloss(13:18,:)))*molc14_to_mg
tempbc14 = sum(pooldt%poollu(13:18))*molc14_to_mg - dpoolendc14
IF (tempb .GT. carbonb_thresh) THEN
   print*,''
   print('(a)'),'!!Dead Pool Daily Error!!'
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,e14.6)'),'Mismatch (Mg C/ha): ',tempb
   print('(a,e14.6)'),'Previous Day Pool: ',sum(pooldt%poollup(1:6))*mol_to_mg
   print('(a,2e14.6)'),'Previous Day Loss/Gain: ', &
        sum(pooldt%poollu_dloss(1:6,:))*mol_to_mg,sum(pooldt%poollu_dgain(1:6,:))*mol_to_mg
   print('(a,2e14.6)'),'New Day Calculated/Saved: ',dpoolendc,sum(pooldt%poollu(1:6))*mol_to_mg
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,''
   print('(a)'),'!!Dead C13 Pool Daily Error!!'
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,e14.6)'),'Mismatch (Mg C/ha): ',tempbc13
   print('(a,e14.6)'),'Previous Day Pool: ',sum(pooldt%poollup(7:12))*molc13_to_mg
   print('(a,2e14.6)'),'Previous Day Loss/Gain: ', &
        sum(pooldt%poollu_dloss(7:12,:))*molc13_to_mg,sum(pooldt%poollu_dgain(7:12,:))*molc13_to_mg
   print('(a,2e14.6)'),'New Day Calculated/Saved: ',dpoolendc,sum(pooldt%poollu(7:12))*molc13_to_mg
   if (carbonb_stop) stop
ENDIF

!...live pools
lpoolendc = (sum(poollt%poolpftp(1:5)) + sum(poollt%poolpft_dgain(1:5,:)) &
     - sum(poollt%poolpft_dloss(1:5,:)))*mol_to_mg
tempb = sum(poollt%poolpft(1:5))*mol_to_mg - lpoolendc
lpoolendc13 = (sum(poollt%poolpftp(6:10)) + sum(poollt%poolpft_dgain(6:10,:)) &
     - sum(poollt%poolpft_dloss(6:10,:)))*molc13_to_mg
tempbc13 = sum(poollt%poolpft(6:10))*molc13_to_mg - lpoolendc13
lpoolendc14 = (sum(poollt%poolpftp(11:15)) + sum(poollt%poolpft_dgain(11:15,:)) &
     - sum(poollt%poolpft_dloss(11:15,:)))*molc14_to_mg
tempbc14 = sum(poollt%poolpft(11:15))*molc14_to_mg - lpoolendc14
IF (tempb .GT. carbonb_thresh) THEN
   print*,''
   print('(a)'),'!!Live Pool Daily Error!!'
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,e14.6)'),'Mismatch (Mg C/ha): ',tempb
   print('(a,e14.6)'),'Previous Day Pool: ',sum(poollt%poolpftp(1:5))*mol_to_mg
   print('(a,2e14.6)'),'Previous Day Loss/Gain: ',&
      sum(poollt%poolpft_dloss(1:5,:))*mol_to_mg, sum(poollt%poolpft_dgain(1:5,:))*mol_to_mg
   print('(a,2e14.6)'),'New Day Calculated/Saved: ',lpoolendc,sum(poollt%poolpft(1:5))*mol_to_mg
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,''
   print('(a)'),'!!Live C13 Pool Daily Error!!'
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,e14.6)'),'C13 Mismatch (Mg C/ha): ',tempbc13
   print('(a,e14.6)'),'C13 Previous Day Pool: ',sum(poollt%poolpftp(6:10))*molc13_to_mg
   print*,'poollt%poolpftp(6): ',poollt%poolpftp(6)*molc13_to_mg
   print*,'poollt%poolpftp(7): ',poollt%poolpftp(7)*molc13_to_mg
   print*,'poollt%poolpftp(8): ',poollt%poolpftp(8)*molc13_to_mg
   print*,'poollt%poolpftp(9): ',poollt%poolpftp(9)*molc13_to_mg
   print*,'poollt%poolpftp(10): ',poollt%poolpftp(10)*molc13_to_mg

   print('(a,2e14.6)'),'C13 Previous Day Loss/Gain: ',&
      sum(poollt%poolpft_dloss(6:10,:))*molc13_to_mg, sum(poollt%poolpft_dgain(6:10,:))*molc13_to_mg
   print*,'poollt%poolpft_dloss(6): ',poollt%poolpft_dloss(6,:)*molc13_to_mg
   print*,'poollt%poolpft_dloss(7): ',poollt%poolpft_dloss(7,:)*molc13_to_mg
   print*,'poollt%poolpft_dloss(8): ',poollt%poolpft_dloss(8,:)*molc13_to_mg
   print*,'poollt%poolpft_dloss(9): ',poollt%poolpft_dloss(9,:)*molc13_to_mg
   print*,'poollt%poolpft_dloss(10): ',poollt%poolpft_dloss(10,:)*molc13_to_mg
   print*,'poollt%poolpft_dgain(6): ',poollt%poolpft_dgain(6,:)*molc13_to_mg
   print*,'poollt%poolpft_dgain(7): ',poollt%poolpft_dgain(7,:)*molc13_to_mg
   print*,'poollt%poolpft_dgain(8): ',poollt%poolpft_dgain(8,:)*molc13_to_mg
   print*,'poollt%poolpft_dgain(9): ',poollt%poolpft_dgain(9,:)*molc13_to_mg
   print*,'poollt%poolpft_dgain(10): ',poollt%poolpft_dgain(10,:)*molc13_to_mg
   print*,' '
   print('(a,2e14.6)'),'New Day Calculated/Saved: ',lpoolendc13,sum(poollt%poolpft(6:10))*molc13_to_mg
   print*,' '
   print*,'rcpoolpft: ',poollt%rcpoolpft(6:10)
   !print*,'d13c: ',(((1.0D0/(1.0D0/poollt%rcpoolpft(6)-1.0D0))-pdb)/pdb)*1000, &
   !                (((1.0D0/(1.0D0/poollt%rcpoolpft(7)-1.0D0))-pdb)/pdb)*1000, &
   !                (((1.0D0/(1.0D0/poollt%rcpoolpft(8)-1.0D0))-pdb)/pdb)*1000, &
   !                (((1.0D0/(1.0D0/poollt%rcpoolpft(9)-1.0D0))-pdb)/pdb)*1000
   print*,' '

   print('(a,e14.6)'),'TotC Mismatch (Mg C/ha): ',tempb
   print('(a,e14.6)'),'TotC Previous Day Pool: ',sum(poollt%poolpftp(1:5))*mol_to_mg
   print*,'poollt%poolpftp(1): ',poollt%poolpftp(1)*mol_to_mg
   print*,'poollt%poolpftp(2): ',poollt%poolpftp(2)*mol_to_mg
   print*,'poollt%poolpftp(3): ',poollt%poolpftp(3)*mol_to_mg
   print*,'poollt%poolpftp(4): ',poollt%poolpftp(4)*mol_to_mg
   print*,'poollt%poolpftp(5): ',poollt%poolpftp(5)*mol_to_mg
   print('(a,2e14.6)'),'TotC Previous Day Loss/Gain: ',&
      sum(poollt%poolpft_dloss(1:5,:))*mol_to_mg, sum(poollt%poolpft_dgain(1:5,:))*mol_to_mg
   print('(a,2e14.6)'),'TotC New Day Calculated/Saved:',lpoolendc,sum(poollt%poolpft(1:5))*mol_to_mg

  if (carbonb_stop) stop
ENDIF


!----- TIME-STEP CARBON BALANCE-----
!...set balance variables
!.....assimilation and fluxes
fassim = assimin*mol_to_umol
assimtot = sum(poollt%gain_assim(1:5))*dtsib &
             + poollt%resp_nveg*dtsib
tempb = abs(assimin*dtsib - assimtot)
fassimc13 = c13assimin*mol_to_umol
assimtotc13 = sum(poollt%gain_assim(6:10))*dtsib &
             + poollt%resp_nvegc13*dtsib
tempbc13 = abs(c13assimin*dtsib - assimtotc13)
fassimc14 = c14assimin*mol_to_umol
assimtotc14 = sum(poollt%gain_assim(11:15))*dtsib &
             + poollt%resp_nvegc14*dtsib
tempbc14 = abs(c14assimin*dtsib - assimtotc14)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Assimilation Error!!'
   print('(a,i6,2f10.2,i4)'),  &
        ' Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print*,'Mismatch (mol C/m2): ',tempb, &
        assimin*dtsib - poollt%resp_nveg*dtsib
   print*,'Assim from phosib: ',co2t%assim*dtsib
   print*,'Assimilation In: ',assimin*dtsib
   print*,'Assimilation C13 In: ',c13assimin*dtsib
   print*,'Live Pool Assim Gain: ',sum(poollt%gain_assim(1:5))*dtsib
   print*,'Non-Veg Resp Loss:',poollt%resp_nveg*dtsib
   print*,'poollt%alloc_phen(1:5)',poollt%alloc_phen(1),poollt%alloc_phen(2),&
        poollt%alloc_phen(3),poollt%alloc_phen(4),poollt%alloc_phen(5)
   print*,'ialloc: ',nint(sum(poollt%alloc_phen(1:5)))
   print*,'poollt%gain_assim(1:5)',poollt%gain_assim(1),poollt%gain_assim(2),&
        poollt%gain_assim(3),poollt%gain_assim(4),poollt%gain_assim(5)
   print*,'poollt%gain_assim(6:10)',poollt%gain_assim(6),poollt%gain_assim(7),&
        poollt%gain_assim(8),poollt%gain_assim(9),poollt%gain_assim(10)
   print*,'Live Pool Assim Gain: ',sum(poollt%gain_assim(1:5))*dtsib
   print*,'Non-Veg Resp Loss:',poollt%resp_nveg*dtsib
   print*,'Live Pool C13 Assim Gain: ',sum(poollt%gain_assim(6:10))*dtsib
   print*,'Non-Veg Resp C13 Loss:',poollt%resp_nvegc13*dtsib
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Assimilation C13 Error!!'
   print('(a,i6,2f10.2,i4)'),  &
        ' Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
   print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
   print*,'Mismatch (mol C/m2): ',tempbc13, &
        c13assimin*dtsib - poollt%resp_nvegc13*dtsib
   print*,'Assimilation C13 In: ',c13assimin*dtsib
   print*,'Assimilation In: ',assimin*dtsib
   print*,' '
   print*,'Assim from phosib: ',co2t%assim
   print*,'Assim C13 from cfrax: ',fract%c13assim   
   print*,' '
!   print*,'poollt%alloc(6:10)',poollt%alloc(6),poollt%alloc(7),&
!        poollt%alloc(8),poollt%alloc(9),poollt%alloc(10)
   print*,'poollt%gain_assim(1:5)',poollt%gain_assim(1),poollt%gain_assim(2),&
        poollt%gain_assim(3),poollt%gain_assim(4),poollt%gain_assim(5)
   print*,'poollt%gain_assim(6:10)',poollt%gain_assim(6),poollt%gain_assim(7),&
        poollt%gain_assim(8),poollt%gain_assim(9),poollt%gain_assim(10)
   print*,'Live Pool Assim Gain: ',sum(poollt%gain_assim(1:5))*dtsib
   print*,'Non-Veg Resp Loss:',poollt%resp_nveg*dtsib
   print*,'Live Pool C13 Assim Gain: ',sum(poollt%gain_assim(6:10))*dtsib
   print*,'Non-Veg Resp C13 Loss:',poollt%resp_nvegc13*dtsib
   if (carbonb_stop) stop
ENDIF


respgrow = sum(poollt%loss_gresp(1:5))*mol_to_umol
tempb = abs(poollt%resp_grow*mol_to_umol - respgrow)
respgrowc13 = sum(poollt%loss_gresp(6:10))*mol_to_umol
tempbc13 = abs(poollt%resp_growc13*mol_to_umol - respgrowc13)
respgrowc14 = sum(poollt%loss_gresp(11:15))*mol_to_umol
tempbc14 = abs(poollt%resp_growc14*mol_to_umol - respgrowc14)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Growth Respiration Error!!'
   print*,'  Resp_Grow: ',poollt%resp_grow*mol_to_umol
   print*,'  Loss_Gresp:',respgrow
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Growth C13 Respiration Error!!'
   print*,'  Resp_Grow C13: ',poollt%resp_growc13*mol_to_umol
   print*,'  Loss_Gresp C13:',respgrowc13
   if (carbonb_stop) stop
ENDIF

respmntn = sum(poollt%loss_mresp_lay(1:5,:))*mol_to_umol
tempb = abs(poollt%resp_mntn*mol_to_umol - respmntn)
respmntnc13 = sum(poollt%loss_mresp_lay(6:10,:))*mol_to_umol
tempbc13 = abs(poollt%resp_mntnc13*mol_to_umol - respmntnc13)
respmntnc14 = sum(poollt%loss_mresp_lay(11:15,:))*mol_to_umol
tempbc14 = abs(poollt%resp_mntnc14*mol_to_umol - respmntnc14)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Maintenance Respiration Error!!'
   print*,'  Resp_Mntn: ',poollt%resp_mntn*mol_to_umol
   print*,'  Loss_Mresp:',respmntn
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Maintenance C13 Respiration Error!!'
   print*,'tempb:',tempb
   print*,'tempbc13: ',tempbc13
   print*,'  Resp_Mntn C13: ',poollt%resp_mntnc13*mol_to_umol
   print*,'  Loss_Mresp C13:',respmntnc13
   print*,'  Resp_Mntn: ',poollt%resp_mntn*mol_to_umol
   print*,'  Loss_Mresp:',respmntn
   print*,'poollt%loss_mresp_lay(1,:)',poollt%loss_mresp_lay(1,:)
   print*,'poollt%loss_mresp_lay(2,:)',poollt%loss_mresp_lay(2,:)
   print*,'poollt%loss_mresp_lay(3,:)',poollt%loss_mresp_lay(3,:)
   print*,'poollt%loss_mresp_lay(4,:)',poollt%loss_mresp_lay(4,:)
   print*,'poollt%loss_mresp_lay(5,:)',poollt%loss_mresp_lay(5,:)
   print*,'poollt%loss_mresp_lay(6,:)',poollt%loss_mresp_lay(6,:)
   print*,'poollt%loss_mresp_lay(7,:)',poollt%loss_mresp_lay(7,:)
   print*,'poollt%loss_mresp_lay(8,:)',poollt%loss_mresp_lay(8,:)
   print*,'poollt%loss_mresp_lay(9,:)',poollt%loss_mresp_lay(9,:)
   print*,'poollt%loss_mresp_lay(10,:)',poollt%loss_mresp_lay(10,:)

   if (carbonb_stop) stop
ENDIF

respleaf = (poollt%loss_gresp(lp) &
    + poollt%loss_mresp_lay(lp,1)) * mol_to_umol
tempb = abs(poollt%resp_leaf*mol_to_umol - respleaf)
respleafc13 = (poollt%loss_gresp(lpc13) &
    + poollt%loss_mresp_lay(lpc13,1)) * mol_to_umol
tempbc13 = abs(poollt%resp_leafc13*mol_to_umol - respleafc13)
respleafc14 = (poollt%loss_gresp(lpc14) &
    + poollt%loss_mresp_lay(lpc14,1)) * mol_to_umol
tempbc14 = abs(poollt%resp_leafc14*mol_to_umol - respleafc14)
IF (tempb .GT. carbonb_thresh) THEN
    print*,'!!Leaf Respiration Error!!'
    print*,'  Resp_Leaf: ',poollt%resp_leaf*mol_to_umol
    print*,'  Loss_Leaf:',respleaf
    if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
    print*,'!!Leaf C13 Respiration Error!!'
    print*,'  Resp_Leaf C13: ',poollt%resp_leafc13*mol_to_umol
    print*,'  Loss_Leaf C13:',respleafc13
    if (carbonb_stop) stop
ENDIF

resproot = (poollt%loss_gresp(frp) + poollt%loss_gresp(crp) &
    + sum(poollt%loss_mresp_lay(frp,:)) &
    + sum(poollt%loss_mresp_lay(crp,:))) * mol_to_umol
tempb = abs(poollt%resp_root*mol_to_umol - resproot)
resprootc13 = (poollt%loss_gresp(frpc13) + poollt%loss_gresp(crpc13) &
    + sum(poollt%loss_mresp_lay(frpc13,:)) &
    + sum(poollt%loss_mresp_lay(crpc13,:))) * mol_to_umol
tempbc13 = abs(poollt%resp_rootc13*mol_to_umol - resprootc13)
resprootc14 = (poollt%loss_gresp(frpc14) + poollt%loss_gresp(crpc14) &
    + sum(poollt%loss_mresp_lay(frpc14,:)) &
    + sum(poollt%loss_mresp_lay(crpc14,:))) * mol_to_umol
tempbc14 = abs(poollt%resp_rootc14*mol_to_umol - resprootc14)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Root Respiration Error!!'
   print*,'  Resp_Root: ',poollt%resp_root*mol_to_umol
   print*,'  Loss_Root: ',resproot
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,' '
   print*,'tempb: ',tempb
   print*,'tempbc13: ',tempbc13
   print*,'!!Root C13 Respiration Error!!'
   print*,'  Resp_Root: ',poollt%resp_root*mol_to_umol
   print*,'  Loss_Root: ',resproot
   print*,'  Resp_Root C13: ',poollt%resp_rootc13*mol_to_umol
   print*,'  Loss_Root C13: ',resprootc13
   print*,'poollt%loss_gresp(frpc13): ',poollt%loss_gresp(frpc13)
   print*,'poollt%loss_gresp(crpc13): ',poollt%loss_gresp(crpc13)
   print*,'sum(poollt%loss_mresp_lay(frpc13,:)) ',sum(poollt%loss_mresp_lay(frpc13,:))
   print*,'sum(poollt%loss_mresp_lay(crpc13,:))',sum(poollt%loss_mresp_lay(crpc13,:))
   print*, ' '
   if (carbonb_stop) stop
ENDIF

respauto = sum(poollt%loss_gresp(1:5))*mol_to_umol &
     + sum(poollt%loss_mresp_lay(1:5,:))*mol_to_umol
tempb = abs(poollt%resp_auto*mol_to_umol - respauto)
respautoc13 = sum(poollt%loss_gresp(6:10))*mol_to_umol &
     + sum(poollt%loss_mresp_lay(6:10,:))*mol_to_umol
tempbc13 = abs(poollt%resp_autoc13*mol_to_umol - respautoc13)
respautoc14 = sum(poollt%loss_gresp(11:15))*mol_to_umol &
     + sum(poollt%loss_mresp_lay(11:15,:))*mol_to_umol
tempbc14 = abs(poollt%resp_autoc14*mol_to_umol - respautoc14)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Autotrophic Respiration Error!!'
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Autotrophic C13 Respiration Error!!'
   if (carbonb_stop) stop
ENDIF

resphet = sum(pooldt%loss_resp_lay(1:6,:))*mol_to_umol
tempb = abs(pooldt%resp_het*mol_to_umol - resphet)
resphetc13 = sum(pooldt%loss_resp_lay(7:12,:))*mol_to_umol
tempbc13 = abs(pooldt%resp_hetc13*mol_to_umol - resphetc13)
resphetc14 = sum(pooldt%loss_resp_lay(13:18,:))*mol_to_umol
tempbc14 = abs(pooldt%resp_hetc14*mol_to_umol - resphetc14)
IF (tempb .GT. carbonb_thresh) THEN
    print*,'!!Heterotrophic Respiration Error!!'
    print*,'Loss/Resp/Diff: ', &
        resphet, pooldt%resp_het*mol_to_umol, tempb
    if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
    print*,'!!Heterotrophic C13 Respiration Error!!'
    print*,'Loss/Resp/Diff: ', &
        resphet, pooldt%resp_het*mol_to_umol, tempb
    print*,'Loss/Resp/Diff C13: ', &
        resphetc13, pooldt%resp_hetc13*mol_to_umol, tempbc13
    if (carbonb_stop) stop
ENDIF

respsoil = dzero
do p=1,npoolpft/3 !1-5
   if (pool_indx_lay(p) .GT. 1) then !1-5, pool_indx_lay(ntpool=22)
       respsoil = respsoil &
          + poollt%loss_gresp(p)*mol_to_umol &
          + sum(poollt%loss_mresp_lay(p,:))*mol_to_umol
   endif
enddo
do p=1,npoollu/3 !1-6
   if (pool_indx_lay(p+npoolpft/3) .GT. 1) then !6-11
       respsoil = respsoil &
         + sum(pooldt%loss_resp_lay(p,:))*mol_to_umol
    endif
enddo
tempb = abs(pooldt%resp_soil*mol_to_umol - respsoil)

!...same as above but for C-13 pools
respsoilc13 = dzero
do p=npoolpft/3+1,2*npoolpft/3 !6-10
   if (pool_indx_lay(p+npoolpft/3+1) .GT. 1) then !ntpool 12,16
       respsoilc13 = respsoilc13 &
          + poollt%loss_gresp(p)*mol_to_umol &
          + sum(poollt%loss_mresp_lay(p,:))*mol_to_umol
   endif
enddo
do p=npoollu/3+1,2*npoollu/3 !7-12
   if (pool_indx_lay(p+2*npoolpft/3) .GT. 1) then !ntpool 17,22
       respsoilc13 = respsoilc13 &
         + sum(pooldt%loss_resp_lay(p,:))*mol_to_umol
    endif
enddo
tempbc13 = abs(pooldt%resp_soilc13*mol_to_umol - respsoilc13)
IF (tempb .GT. carbonb_thresh) THEN
   print*,'!!Soil Respiration Error!!'
   print*,'   Soil_Resp: ',pooldt%resp_soil*mol_to_umol
   print*,'   Soil_RLoss:',respsoil
   if (carbonb_stop) stop
ENDIF
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Soil C13 Respiration Error!!'
   print*,'   Soil_Resp C13: ',pooldt%resp_soilc13*mol_to_umol
   print*,'   Soil_RLoss C13:',respsoilc13
   if (carbonb_stop) stop
ENDIF

!...same as above but for C-14 pools
respsoilc14 = dzero
do p=2*npoolpft/3+1,npoolpft !11-15
   if (pool_indx_lay(p+2*npoolpft/3+2) .GT. 1) then !23,27
       respsoilc14 = respsoilc14 &
          + poollt%loss_gresp(p)*mol_to_umol &
          + sum(poollt%loss_mresp_lay(p,:))*mol_to_umol
   endif
enddo
do p=2*npoollu/3+1,npoollu !13-18
   if (pool_indx_lay(p+npoolpft) .GT. 1) then !28,33
       respsoilc14 = respsoilc14 &
         + sum(pooldt%loss_resp_lay(p,:))*mol_to_umol
    endif
enddo
tempbc14 = abs(pooldt%resp_soilc14*mol_to_umol - respsoilc14)


fresp = (pooldt%resp_het + poollt%resp_auto &
        + poollt%resp_nveg + poollt%resp_fire  &
        + poollt%resp_grz + poollt%resp_hrvst) &
        * mol_to_umol

frespc13 = (pooldt%resp_hetc13 + poollt%resp_autoc13 &
        + poollt%resp_nvegc13 + poollt%resp_firec13  &
        + poollt%resp_grzc13 + poollt%resp_hrvstc13) &
        * mol_to_umol

frespc14 = (pooldt%resp_hetc14 + poollt%resp_autoc14 &
        + poollt%resp_nvegc14 + poollt%resp_firec14  &
        + poollt%resp_grzc14 + poollt%resp_hrvstc14) &
        * mol_to_umol

!.....dead pools
do p=1,npoollu/3 !1,6
   dgain(p) = sum(pooldt%gain_grz_lay(p,:))*mol_to_mg &
            + sum(pooldt%gain_hrvst_lay(p,:))*mol_to_mg &
            + sum(pooldt%gain_transd_lay(p,:))*mol_to_mg &
            + sum(pooldt%gain_transl_lay(p,:))*mol_to_mg
   dloss(p) = & 
            + sum(pooldt%loss_resp_lay(p,:))*mol_to_mg &
            + sum(pooldt%loss_trans_lay(p,:))*mol_to_mg &
            + sum(pooldt%loss_fire_lay(p,:))*mol_to_mg &
            + sum(pooldt%loss_raddecay_lay(p,:))*mol_to_mg
enddo
do p=npoollu/3+1,2*npoollu/3 !7,12
   dgain(p) = sum(pooldt%gain_grz_lay(p,:))*molc13_to_mg &
            + sum(pooldt%gain_hrvst_lay(p,:))*molc13_to_mg &
            + sum(pooldt%gain_transd_lay(p,:))*molc13_to_mg &
            + sum(pooldt%gain_transl_lay(p,:))*molc13_to_mg
   dloss(p) = &
            + sum(pooldt%loss_resp_lay(p,:))*molc13_to_mg &
            + sum(pooldt%loss_trans_lay(p,:))*molc13_to_mg &
            + sum(pooldt%loss_fire_lay(p,:))*molc13_to_mg &
            + sum(pooldt%loss_raddecay_lay(p,:))*molc13_to_mg
enddo
do p=2*npoollu/3+1,npoollu !13,18
   dgain(p) = sum(pooldt%gain_grz_lay(p,:))*molc14_to_mg &
            + sum(pooldt%gain_hrvst_lay(p,:))*molc14_to_mg &
            + sum(pooldt%gain_transd_lay(p,:))*molc14_to_mg &
            + sum(pooldt%gain_transl_lay(p,:))*molc14_to_mg
   dloss(p) = &
            + sum(pooldt%loss_resp_lay(p,:))*molc14_to_mg &
            + sum(pooldt%loss_trans_lay(p,:))*molc14_to_mg &
            + sum(pooldt%loss_fire_lay(p,:))*molc14_to_mg &
            + sum(pooldt%loss_raddecay_lay(p,:))*molc14_to_mg
enddo

dgain = dgain*dtsib
dloss = dloss*dtsib
dpoolinit = (sum(pooldt%poollu(1:6))*mol_to_mg &
            - sum(dgain(1:6)) + sum(dloss(1:6)))
dpoolend = sum(pooldt%poollu(1:6))*mol_to_mg
dpoolchange = dpoolend - dpoolinit

dpoolinitc13 = (sum(pooldt%poollu(7:12))*molc13_to_mg &
            - sum(dgain(7:12)) + sum(dloss(7:12)))
dpoolendc13 = sum(pooldt%poollu(7:12))*molc13_to_mg
dpoolchangec13 = dpoolendc13 - dpoolinitc13

dpoolinitc14 = (sum(pooldt%poollu(13:18))*molc14_to_mg &
            - sum(dgain(13:18)) + sum(dloss(13:18)))
dpoolendc14 = sum(pooldt%poollu(13:18))*molc14_to_mg
dpoolchangec14 = dpoolendc14 - dpoolinitc14

dgaingrz = sum(pooldt%gain_grz_lay(1:6,:))*mol_to_umol*dtsib
dgainhrv = sum(pooldt%gain_hrvst_lay(1:6,:))*mol_to_umol*dtsib
dgaintd = sum(pooldt%gain_transd_lay(1:6,:))*mol_to_umol*dtsib
dgaintl = sum(pooldt%gain_transl_lay(1:6,:))*mol_to_umol*dtsib
dcarbonin = dgaingrz + dgainhrv &
              + dgaintd + dgaintl

dlossr = sum(pooldt%loss_resp_lay(1:6,:))*mol_to_umol*dtsib
dlosstd = sum(pooldt%loss_trans_lay(1:6,:))*mol_to_umol*dtsib
dlossf = sum(pooldt%loss_fire_lay(1:6,:))*mol_to_umol*dtsib
dlossrdd = sum(pooldt%loss_raddecay_lay(1:6,:))*mol_to_umol*dtsib
dcarbonout = dlossr + dlosstd + dlossf + dlossrdd

dcarbons = dpoolchange/mol_to_mg*mol_to_umol
dcarbonb = dcarbonin - dcarbonout - dcarbons
tempb = abs(dcarbonb)/mol_to_umol
!if (tempb .GT. carbonb_thresh) cb_err = .true.
IF (tempb .GT. carbonb_thresh) THEN
   cb_err = .true.
   print*,'!!Dead Pool Store/In/Out Balance Error!!'
   print*,'   DPool_store: ',dcarbons
   print*,'   DPool_gain:  ',dcarbonin
   print*,'   DPool_loss:  ',dcarbonout
   if (carbonb_stop) stop
ENDIF

!...same as above but for C13 pools
dgaingrzc13 = sum(pooldt%gain_grz_lay(7:12,:))*mol_to_umol*dtsib
dgainhrvc13 = sum(pooldt%gain_hrvst_lay(7:12,:))*mol_to_umol*dtsib
dgaintdc13 = sum(pooldt%gain_transd_lay(7:12,:))*mol_to_umol*dtsib
dgaintlc13 = sum(pooldt%gain_transl_lay(7:12,:))*mol_to_umol*dtsib
dcarboninc13 = dgaingrzc13 + dgainhrvc13 &
              + dgaintdc13 + dgaintlc13

dlossrc13 = sum(pooldt%loss_resp_lay(7:12,:))*mol_to_umol*dtsib
dlosstdc13 = sum(pooldt%loss_trans_lay(7:12,:))*mol_to_umol*dtsib
dlossfc13 = sum(pooldt%loss_fire_lay(7:12,:))*mol_to_umol*dtsib
dlossrddc13 = sum(pooldt%loss_raddecay_lay(7:12,:))*mol_to_umol*dtsib
dcarbonoutc13 = dlossrc13 + dlosstdc13 + dlossfc13 + dlossrddc13

dcarbonsc13 = dpoolchangec13/molc13_to_mg*mol_to_umol
dcarbonbc13 = dcarboninc13 - dcarbonoutc13 - dcarbonsc13
tempbc13 = abs(dcarbonbc13)/mol_to_umol
!if (tempbc13 .GT. carbonb_threshc13) cb_errc13 = .true.
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Dead C13 Pool Store/In/Out Balance Error!!'
   print*,'   '
   print*,'   DPool_store: ',dcarbons
   print*,'   DPool_gain:  ',dcarbonin
   print*,'   DPool_loss:  ',dcarbonout
   print*,'   '
   print*,'   tempbc13: ',tempbc13
   print*,'   DPool_store C13: ',dcarbonsc13
   print*,'   DPool_gain C13:  ',dcarboninc13
   print*,'   DPool_loss C13:  ',dcarbonoutc13
   print*,'   '
   print*,'   dgaingrzc13: ',dgaingrzc13
   print*,'   dgainhrvc13: ',dgaingrzc13
   print*,'   dgaintdc13: ',dgaingrzc13
   print*,'   dgaintlc13: ',dgaingrzc13
   print*,'   '
   print*,'   dlossrc13: ',dlossrc13
   print*,'   dlosstdc13: ',dlosstdc13
   print*,'   dlossfc13: ',dlossfc13
   print*,'   dlossrddc13: ',dlossrddc13

   if (carbonb_stop) stop
ENDIF

!...same as above but for C14 pools
dgaingrzc14 = sum(pooldt%gain_grz_lay(13:18,:))*mol_to_umol*dtsib
dgainhrvc14 = sum(pooldt%gain_hrvst_lay(13:18,:))*mol_to_umol*dtsib
dgaintdc14 = sum(pooldt%gain_transd_lay(13:18,:))*mol_to_umol*dtsib
dgaintlc14 = sum(pooldt%gain_transl_lay(13:18,:))*mol_to_umol*dtsib
dcarboninc14 = dgaingrzc14 + dgainhrvc14 &
              + dgaintdc14 + dgaintlc14

dlossrc14 = sum(pooldt%loss_resp_lay(13:18,:))*mol_to_umol*dtsib
dlosstdc14 = sum(pooldt%loss_trans_lay(13:18,:))*mol_to_umol*dtsib
dlossfc14 = sum(pooldt%loss_fire_lay(13:18,:))*mol_to_umol*dtsib
dlossrddc14 = sum(pooldt%loss_raddecay_lay(13:18,:))*mol_to_umol*dtsib

dcarbonoutc14 = dlossrc14 + dlosstdc14 + dlossfc14 + dlossrddc14

dcarbonsc14 = dpoolchangec14/molc14_to_mg*mol_to_umol
dcarbonbc14 = dcarboninc14 - dcarbonoutc14 - dcarbonsc14
tempbc14 = abs(dcarbonbc14)/mol_to_umol

!.....live pools
do p=1,npoolpft/3 !1,5
   lgain(p) = poollt%gain_assim(p)*mol_to_mg*dtsib &
            + poollt%gain_seed(p)*mol_to_mg
   lloss(p) = &
            sum(poollt%loss_fire_lay(p,:))*mol_to_mg*dtsib &
            + sum(poollt%loss_hrvst_lay(p,:))*mol_to_mg*dtsib &
            + poollt%loss_gresp(p)*mol_to_mg*dtsib &
            + sum(poollt%loss_mresp_lay(p,:))*mol_to_mg*dtsib &
            !+ poollt%resp_nveg*mol_to_mg*dtsib &
            + sum(poollt%loss_trans_lay(p,:))*mol_to_mg*dtsib &
            + sum(poollt%loss_raddecay_lay(p,:))*mol_to_mg*dtsib
enddo
do p=npoolpft/3+1,2*npoolpft/3 !6,10
   lgain(p) = poollt%gain_assim(p)*molc13_to_mg*dtsib &
            + poollt%gain_seed(p)*molc13_to_mg
   lloss(p) = &
            sum(poollt%loss_fire_lay(p,:))*molc13_to_mg*dtsib &
            + sum(poollt%loss_hrvst_lay(p,:))*molc13_to_mg*dtsib &
            + poollt%loss_gresp(p)*molc13_to_mg*dtsib &
            + sum(poollt%loss_mresp_lay(p,:))*molc13_to_mg*dtsib &
            !+ poollt%resp_nveg*molc13_to_mg*dtsib &
            + sum(poollt%loss_trans_lay(p,:))*molc13_to_mg*dtsib &
            + sum(poollt%loss_raddecay_lay(p,:))*molc13_to_mg*dtsib
enddo
do p=2*npoolpft/3+1,npoolpft !11,15
   lgain(p) = poollt%gain_assim(p)*molc14_to_mg*dtsib &
            + poollt%gain_seed(p)*molc14_to_mg
   lloss(p) = &
            sum(poollt%loss_fire_lay(p,:))*molc14_to_mg*dtsib &
            + sum(poollt%loss_hrvst_lay(p,:))*molc14_to_mg*dtsib &
            + poollt%loss_gresp(p)*molc14_to_mg*dtsib &
            + sum(poollt%loss_mresp_lay(p,:))*molc14_to_mg*dtsib &
            !+ poollt%resp_nveg*molc14_to_mg*dtsib &
            + sum(poollt%loss_trans_lay(p,:))*molc14_to_mg*dtsib &
            + sum(poollt%loss_raddecay_lay(p,:))*molc14_to_mg*dtsib
enddo

do p=1,npoolcan
   cref = pool_indx_can(p) !(1,4,5)
   lloss(cref) = lloss(cref) + poollt%loss_grz(p)*mol_to_mg*dtsib
   crefc13 = pool_indx_can(p)+5 !(6,9,10) live C13 pool indexed by canopy
   lloss(crefc13) = lloss(crefc13) + poollt%loss_grzc13(p)*molc13_to_mg*dtsib
   crefc14 = pool_indx_can(p)+10 !(11,14,15) live C14 pool indexed by canopy
   lloss(crefc14) = lloss(crefc14) + poollt%loss_grzc14(p)*molc14_to_mg*dtsib
enddo

lpoolinit = (sum(poollt%poolpft(1:5))*mol_to_mg &
             - sum(lgain(1:5)) + sum(lloss(1:5)))
lpoolend = sum(poollt%poolpft(1:5))*mol_to_mg
lpoolchange = lpoolend - lpoolinit

lpoolinitc13 = (sum(poollt%poolpft(6:10))*molc13_to_mg &
             - sum(lgain(6:10)) + sum(lloss(6:10)))
lpoolendc13 = sum(poollt%poolpft(6:10))*molc13_to_mg
lpoolchangec13 = lpoolendc13 - lpoolinitc13

lpoolinitc14 = (sum(poollt%poolpft(11:15))*molc14_to_mg &
             - sum(lgain(11:15)) + sum(lloss(11:15)))
lpoolendc14 = sum(poollt%poolpft(11:15))*molc14_to_mg
lpoolchangec14 = lpoolendc14 - lpoolinitc14

!lgaina = sum(poollt%gain_assim)*mol_to_umol*dtsib
lgaina = assimin*mol_to_umol*dtsib
lgains = sum(poollt%gain_seed(1:5))*mol_to_umol
lcarbonin = lgaina + lgains

llossfire = sum(poollt%loss_fire_lay(1:5,:))*mol_to_umol*dtsib
llossrdd = sum(poollt%loss_raddecay_lay(1:5,:))*mol_to_umol*dtsib
llossgrz = sum(poollt%loss_grz(1:3))*mol_to_umol*dtsib !npoolcan
llossrgrz = poollt%resp_grz*mol_to_umol*dtsib
llosstgrz = sum(pooldt%gain_grz_lay(1:6,:))*mol_to_umol*dtsib

llosshrv = sum(poollt%loss_hrvst_lay(1:5,:))*mol_to_umol*dtsib
llossrg = sum(poollt%loss_gresp(1:5))*mol_to_umol*dtsib
llossrm = sum(poollt%loss_mresp_lay(1:5,:))*mol_to_umol*dtsib
llossrnveg = poollt%resp_nveg*mol_to_umol*dtsib
llosstd = sum(poollt%loss_trans_lay(1:5,:))*mol_to_umol*dtsib

lcarbonout = llossgrz + llosshrv + llossfire + llossrdd &
     + llossrg + llossrm + llossrnveg + llosstd
lcarbons = lpoolchange/mol_to_mg*mol_to_umol

lcarbonb = lcarbonin - lcarbonout - lcarbons
tempb = abs(lcarbonb)*myconvert

if (tempb .GT. carbonb_thresh) then
    cb_err = .true.
   print*,'!!Live Pool Store/In/Out Balance Error!!'
   print*,'   Live Pool_store: ',lcarbons
   print*,'   Live Pool_gain:  ',lcarbonin
   print*,'   Live Pool_loss:  ',lcarbonout
!   print*,'   Live C13 Pool_store: ',lcarbonsc13
!   print*,'   Live C13 Pool_gain:  ',lcarboninc13
!   print*,'   Live C13 Pool_loss:  ',lcarbonoutc13
   if (carbonb_stop) stop
endif

!...same as above but for C13 pools
!lgaina = sum(poollt%gain_assim)*mol_to_umol*dtsib
lgainac13 = c13assimin*mol_to_umol*dtsib
lgainsc13 = sum(poollt%gain_seed(6:10))*mol_to_umol
lcarboninc13 = lgainac13 + lgainsc13

llossfirec13 = sum(poollt%loss_fire_lay(6:10,:))*mol_to_umol*dtsib
llossrddc13 = sum(poollt%loss_raddecay_lay(6:10,:))*mol_to_umol*dtsib
llossgrzc13 = sum(poollt%loss_grzc13(1:3))*mol_to_umol*dtsib !npoolcanc13
llossrgrzc13 = poollt%resp_grzc13*mol_to_umol*dtsib
llosstgrzc13 = sum(pooldt%gain_grz_lay(7:12,:))*mol_to_umol*dtsib

llosshrvc13 = sum(poollt%loss_hrvst_lay(6:10,:))*mol_to_umol*dtsib
llossrgc13 = sum(poollt%loss_gresp(6:10))*mol_to_umol*dtsib
llossrmc13 = sum(poollt%loss_mresp_lay(6:10,:))*mol_to_umol*dtsib
llossrnvegc13 = poollt%resp_nvegc13*mol_to_umol*dtsib
llosstdc13 = sum(poollt%loss_trans_lay(6:10,:))*mol_to_umol*dtsib

lcarbonoutc13 = llossgrzc13 + llosshrvc13 + llossfirec13 + llossrddc13 &
     + llossrgc13 + llossrmc13 + llossrnvegc13 + llosstdc13
lcarbonsc13 = lpoolchangec13/molc13_to_mg*mol_to_umol

lcarbonbc13 = lcarboninc13 - lcarbonoutc13 - lcarbonsc13
tempbc13 = abs(lcarbonbc13)*myconvert
!if (tempbc13 .GT. carbonb_threshc13) cb_errc13 = .true.

!if (tempb .GT. carbonb_thresh) then
!    cb_err = .true.
!   print*,'!!Live Pool Store/In/Out Balance Error!!'
!   print*,'   Live Pool_store: ',lcarbons
!   print*,'   Live Pool_gain:  ',lcarbonin
!   print*,'   Live Pool_loss:  ',lcarbonout
!   print*,'   Live C13 Pool_store: ',lcarbonsc13
!   print*,'   Live C13 Pool_gain:  ',lcarboninc13
!   print*,'   Live C13 Pool_loss:  ',lcarbonoutc13
!   if (carbonb_stop) stop
!endif


IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Live C13 Pool Store/In/Out Balance Error!!'
   print*,'   Live C13 Pool_store: ',lcarbonsc13
   print*,'   Live C13 Pool_gain:  ',lcarboninc13
   print*,'   Live C13 Pool_loss:  ',lcarbonoutc13
!   print*,'   lgain seed c13:      ',lgainsc13
!   print*,'   lgain assim c13:     ',lgainac13
   print*,' '
!   print*,'   Live Pool_store: ',lcarbons
!   print*,'   Live Pool_gain:  ',lcarbonin
!   print*,'   Live Pool_loss:  ',lcarbonout
!   print*,'   lgain seed:      ',lgains
!   print*,'   lgain assim:     ',lgaina

   if (carbonb_stop) stop
ENDIF

!...same as above but for C14 pools
!lgaina = sum(poollt%gain_assim)*mol_to_umol*dtsib
lgainac14 = c14assimin*mol_to_umol*dtsib
lgainsc14 = sum(poollt%gain_seed(11:15))*mol_to_umol
lcarboninc14 = lgainac14 + lgainsc14

llossfirec14 = sum(poollt%loss_fire_lay(11:15,:))*mol_to_umol*dtsib
llossrddc14 = sum(poollt%loss_raddecay_lay(11:15,:))*mol_to_umol*dtsib
llossgrzc14 = sum(poollt%loss_grzc14(1:3))*mol_to_umol*dtsib !npoolcanc14
llossrgrzc14 = poollt%resp_grzc14*mol_to_umol*dtsib
llosstgrzc14 = sum(pooldt%gain_grz_lay(13:18,:))*mol_to_umol*dtsib

llosshrvc14 = sum(poollt%loss_hrvst_lay(11:15,:))*mol_to_umol*dtsib
llossrgc14 = sum(poollt%loss_gresp(11:15))*mol_to_umol*dtsib
llossrmc14 = sum(poollt%loss_mresp_lay(11:15,:))*mol_to_umol*dtsib
llossrnvegc14 = poollt%resp_nvegc14*mol_to_umol*dtsib
llosstdc14 = sum(poollt%loss_trans_lay(11:15,:))*mol_to_umol*dtsib

lcarbonoutc14 = llossgrzc14 + llosshrvc14 + llossfirec14 + llossrddc14 &
     + llossrgc14 + llossrmc14 + llossrnvegc14 + llosstdc14
lcarbonsc14 = lpoolchangec14/molc14_to_mg*mol_to_umol

lcarbonbc14 = lcarboninc14 - lcarbonoutc14 - lcarbonsc14
tempbc14 = abs(lcarbonbc14)*myconvert

!.....grazing carbon balance
grzcarbonb = llossgrz - llossrgrz - llosstgrz
tempb = abs(grzcarbonb)*myconvert

grzcarbonbc13 = llossgrzc13 - llossrgrzc13 - llosstgrzc13
tempbc13 = abs(grzcarbonbc13)*myconvert

grzcarbonbc14 = llossgrzc14 - llossrgrzc14 - llosstgrzc14
tempbc14 = abs(grzcarbonbc14)*myconvert

if (tempb .GT. carbonb_thresh) then 
   cb_err = .true.
   print*,'!!Grazing Balance Error!!'
   print*,'   llosgrz:    ',llossgrz
   print*,'   llossrgrz:  ',llossrgrz
   print*,'   llosstgrz:  ',llosstgrz
   print*,' '
   print*,'   llosgrzc13:    ',llossgrzc13
   print*,'   llossrgrzc13:  ',llossrgrzc13
   print*,'   llosstgrzc13:  ',llosstgrzc13
   if (carbonb_stop) stop
endif

grzcarbonbc13 = llossgrzc13 - llossrgrzc13 - llosstgrzc13
tempbc13 = abs(grzcarbonbc13)*myconvert
!if (tempbc13 .GT. carbonb_threshc13) cb_errc13 = .true.
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Grazing C13 Balance Error!!'
   print*,'   llosgrzc13:    ',llossgrzc13
   print*,'   llossrgrzc13:  ',llossrgrzc13
   print*,'   llosstgrzc13:  ',llosstgrzc13
!   print*,'   grztransfer3-8: ',grz_transfer(3:8)
!   print*,'   grztransfer9-14: ',grz_transfer(9:14)
   print*,' '
!   print*,'   llosgrz:       ',llossgrz
!   print*,'   llossgrzr:     ',llossrgrz
!   print*,'   llossgrzt:     ',llosstgrz
   if (carbonb_stop) stop
ENDIF

grzcarbonbc14 = llossgrzc14 - llossrgrzc14 - llosstgrzc14
tempbc14 = abs(grzcarbonbc14)*myconvert

!.....harvest carbon balance
hrvresp = poollt%resp_hrvst*dtsib*mol_to_umol
hrvrmvd = poollt%rmvd_hrvst*mol_to_umol
hrvcarbonb = llosshrv &
                 - dgainhrv - hrvresp - hrvrmvd
tempb = abs(hrvcarbonb)*myconvert
if (tempb .GT. carbonb_thresh) then
   cb_err = .true.
   print*,'!!Harvest Balance Error!!'
   print*,'   Live Pool_store: ',lcarbons
   print*,'   Live Pool_gain:  ',lcarbonin
   print*,'   Live Pool_loss:  ',lcarbonout
   print*,'   llosshrv:    ',llosshrv
   print*,'   dgainhrv:    ',dgainhrv
   print*,'   hrvresp:     ',hrvresp
   print*,'   hrvrmvd:     ',hrvrmvd
   print*,'   hrvcarbonb:  ',hrvcarbonb
   print*,'   DPool_store: ',dcarbons
   print*,'   DPool_gain:  ',dcarbonin
   print*,'   DPool_loss:  ',dcarbonout
   if (carbonb_stop) stop
endif

hrvrespc13 = poollt%resp_hrvstc13*dtsib*mol_to_umol
hrvrmvdc13 = poollt%rmvd_hrvstc13*mol_to_umol
hrvcarbonbc13 = llosshrvc13 &
                 - dgainhrvc13 - hrvrespc13 - hrvrmvdc13
tempbc13 = abs(hrvcarbonbc13)*myconvert
!if (tempbc13 .GT. carbonb_threshc13) cb_errc13 = .true.
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Harvest C13 Balance Error!!'
   print*,'   hrvrespc13:    ',hrvrespc13
   print*,'   hrvrmvdc13:  ',hrvrmvdc13
   print*,'   hrvcarbonbc13:  ',hrvcarbonbc13
   if (carbonb_stop) stop
ENDIF

hrvrespc14 = poollt%resp_hrvstc14*dtsib*mol_to_umol
hrvrmvdc14 = poollt%rmvd_hrvstc14*mol_to_umol
hrvcarbonbc14 = llosshrvc14 &
                 - dgainhrvc14 - hrvrespc14 - hrvrmvdc14
tempbc14 = abs(hrvcarbonbc14)*myconvert

!.....live-to-dead transfers
tlivedead = llosstd + llosstgrz + (llosshrv - hrvresp - hrvrmvd)
gdeadlive = dgaintl + dgaingrz + dgainhrv
tcarbonb = (tlivedead - gdeadlive) &
         + (llosstd - dgaintl) + (dlosstd - dgaintd)  
tempb = abs(tcarbonb)*myconvert
if (tempb .GT. carbonb_thresh) then 
   cb_err = .true.
   print*,'!!Live-to-Dead Balance Error!!'
   if (carbonb_stop) stop
endif

tlivedeadc13 = llosstdc13 + llosstgrzc13 + &
                (llosshrvc13 - hrvrespc13 - hrvrmvdc13)
gdeadlivec13 = dgaintlc13 + dgaingrzc13 + dgainhrvc13
tcarbonbc13 = (tlivedeadc13 - gdeadlivec13) &
         + (llosstdc13 - dgaintlc13) + (dlosstdc13 - dgaintdc13)
tempbc13 = abs(tcarbonbc13)*myconvert
!if (tempbc13 .GT. carbonb_threshc13) cb_errc13 = .true.
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Live-to-Dead C13 Balance Error!!'
   if (carbonb_stop) stop
ENDIF

tlivedeadc14 = llosstdc14 + llosstgrzc14 + &
                (llosshrvc14 - hrvrespc14 - hrvrmvdc14)
gdeadlivec14 = dgaintlc14 + dgaingrzc14 + dgainhrvc14
tcarbonbc14 = (tlivedeadc14 - gdeadlivec14) &
         + (llosstdc14 - dgaintlc14) + (dlosstdc14 - dgaintdc14)
tempbc14 = abs(tcarbonbc14)*myconvert

!.....net carbon balance
netcarbonb = dcarbonb + lcarbonb + tcarbonb &
           + grzcarbonb + hrvcarbonb
tempb = abs(netcarbonb)*myconvert
IF ((tempb .GT. carbonb_thresh) .and. &
     (.not. seed_printout)) cb_err = .true.

netcarbonbc13 = dcarbonbc13 + lcarbonbc13 + tcarbonbc13 &
           + grzcarbonbc13 + hrvcarbonbc13
tempbc13 = abs(netcarbonbc13)*myconvert
!IF ((tempbc13 .GT. carbonb_threshc13) .and. &
!     (.not. seed_printout)) cb_errc13 = .true.
IF (tempbc13 .GT. carbonb_threshc13) THEN
   print*,'!!Net Carbon C13 Balance Error!!'
   if (carbonb_stop) stop
ENDIF

netcarbonbc14 = dcarbonbc14 + lcarbonbc14 + tcarbonbc14 &
           + grzcarbonbc14 + hrvcarbonbc14
tempbc14 = abs(netcarbonbc14)*myconvert

!.....print results
loc_printout = cb_err .or. carbonb_print
loc_printoutc13 = cb_errc13 .or. carbonb_print

if (loc_printout) then
    print('(a)'),'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    print('(a)'), 'TIME-STEP CARBON CYCLE BALANCE'
    print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
    print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
    print('(a,f12.6,a,f12.6)'),'LAI/FPAR: ', laiin, '   ', fparin

    print('(a)'),''
    print('(a)'),'--CARBON FLUXES (micromoles C/m2/s)--'
    print('(a,f12.6)'),'  In:  ',fassim
    print('(a,f12.6)'),'  Out: ',fresp
    print('(a)'),      '       Resp        Saved        Calculated'
    print('(a,f12.6,a,f12.6)'), &
                       '       Auto    ', &
          poollt%resp_auto*mol_to_umol, '  ', respauto
    print('(a,f12.6,a,f12.6)'), &
                       '       Het     ', &
          pooldt%resp_het*mol_to_umol, '  ', resphet
    print('(a,f12.6,a,f12.6)'), &
                       '       NVeg    ', &
          poollt%resp_nveg*mol_to_umol, '  ', poollt%resp_nveg*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       Graze   ', &
          poollt%resp_grz*mol_to_umol, '  ', poollt%resp_grz*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       Harvest ', &
          poollt%resp_hrvst*mol_to_umol, '  ', poollt%resp_hrvst*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       Growth  ', &
          poollt%resp_grow*mol_to_umol, '  ', respgrow
    print('(a,f12.6,a,f12.6)'), &
                       '       Leaf    ', &
          poollt%resp_leaf*mol_to_umol, '  ', respleaf
    print('(a,f12.6,a,f12.6)'), &
                       '       Mntn    ', &
          poollt%resp_mntn*mol_to_umol, '  ', respmntn
    print('(a,f12.6,a,f12.6)'), &
                       '       Roots   ', &
          poollt%resp_root*mol_to_umol, '  ', resproot
    print('(a,f12.6,a,f12.6)'), &
                       '       Soil    ', &
          pooldt%resp_soil*mol_to_umol, '  ', respsoil
   print('(a,f12.6,a)'),'  NEE: ',fresp-fassim, '  (Out-In)'


    print('(a)'),''
    print('(a)'),'--CARBON POOLS (Mg C/ha)--'
    print('(a,f14.8,a,e14.6,f14.8)'),'  Live Carbon Init/Change/End: ', &
            lpoolinit, '   ', lpoolchange, lpoolend
    print('(a,f14.8,a,e14.6,f14.8)'),'  Dead Carbon Init/Change/End: ', &
            dpoolinit, '   ', dpoolchange, dpoolend

    print('(a)'),''
    print('(a)'),'  Pool     Orig            Gain            Loss              New'
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Leaf  ', &
            (poollt%poolpft(lp) - sum(poollt%poolpft_dgain(lp,:)) &
               + sum(poollt%poolpft_dloss(lp,:)))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dgain(lp,:))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dloss(lp,:))*mol_to_mg, &
               '  ', poollt%poolpft(lp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  FRoot ', &
            (poollt%poolpft(frp)-sum(poollt%poolpft_dgain(frp,:)) &
              + sum(poollt%poolpft_dloss(frp,:)))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dgain(frp,:))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dloss(frp,:))*mol_to_mg, &
              '  ', poollt%poolpft(frp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  CRoot ', &
            (poollt%poolpft(crp)-sum(poollt%poolpft_dgain(crp,:)) &
               + sum(poollt%poolpft_dloss(crp,:)))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dgain(crp,:))*mol_to_mg, &
               '  ', sum(poollt%poolpft_dloss(crp,:))*mol_to_mg, &
               '  ', poollt%poolpft(crp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  StWd  ', &
            (poollt%poolpft(wp)-sum(poollt%poolpft_dgain(wp,:)) &
              + sum(poollt%poolpft_dloss(wp,:)))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dgain(wp,:))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dloss(wp,:))*mol_to_mg, &
              '  ', poollt%poolpft(wp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Prod  ', &
            (poollt%poolpft(pp)-sum(poollt%poolpft_dgain(pp,:)) &
              + sum(poollt%poolpft_dloss(pp,:)))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dgain(pp,:))*mol_to_mg, &
              '  ', sum(poollt%poolpft_dloss(pp,:))*mol_to_mg, &
              '  ', poollt%poolpft(pp)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  CDB   ', &
            (pooldt%poollu(cdb)-sum(pooldt%poollu_dgain(cdb,:)) &
               + sum(pooldt%poollu_dloss(cdb,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(cdb,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(cdb,:))*mol_to_mg, &
              '  ', pooldt%poollu(cdb)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  MetL  ', &
            (pooldt%poollu(metl)-sum(pooldt%poollu_dgain(metl,:)) &
              + sum(pooldt%poollu_dloss(metl,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(metl,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(metl,:))*mol_to_mg, &
              '  ', pooldt%poollu(metl)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  StrL  ', &
            (pooldt%poollu(strl)-sum(pooldt%poollu_dgain(strl,:)) &
               + sum(pooldt%poollu_dloss(strl,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(strl,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(strl,:))*mol_to_mg, &
              '  ', pooldt%poollu(strl)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  SLit  ', &
            (pooldt%poollu(slit)-sum(pooldt%poollu_dgain(slit,:)) &
               + sum(pooldt%poollu_dloss(slit,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(slit,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(slit,:))*mol_to_mg, &
              '  ', pooldt%poollu(slit)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Slow  ', &
            (pooldt%poollu(slow)-sum(pooldt%poollu_dgain(slow,:)) &
               + sum(pooldt%poollu_dloss(slow,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(slow,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(slow,:))*mol_to_mg, &
              '  ', pooldt%poollu(slow)*mol_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Arm   ', &
            (pooldt%poollu(arm)-sum(pooldt%poollu_dgain(arm,:)) &
               + sum(pooldt%poollu_dloss(arm,:)))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dgain(arm,:))*mol_to_mg, &
              '  ', sum(pooldt%poollu_dloss(arm,:))*mol_to_mg, &
              '  ', pooldt%poollu(arm)*mol_to_mg

    print('(a)'),''
    print('(a)'),'--LIVE CARBON BALANCE (micromoles C/m2)--'
    print('(a,e14.6)'),'  Balance Error:  ',lcarbonb
    if (seed_printout) print('(a)'),'  Non-Minimal Error Due To Crop Seed Release'
    print('(a)'),       '  Gain:      Assim          Seed_Release'
    print('(a,e14.6,a,e14.6)'),         '           ', &
         lgaina, ' ', lgains
    print('(a)'),       '  Loss:      Fire           Harvest        Tran_Dead'
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
         llossfire, ' ', llosshrv, ' ', llosstd
    print('(a)'),       '  Resp:      Growth         Maintenance    Non-Veg'
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
         llossrg, ' ', llossrm, ' ', llossrnveg
    IF (abs(grzcarbonb) .GT. dzero) THEN
         print('(a)'), ' Graze:      Loss           Resp           Tran_Dead'
         print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
               llossgrz, '  ', llossrgrz, '  ', llosstgrz
    ENDIF
    print('(a)'),       '  Net:       In             Out            Stored'
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
          lcarbonin, '   ', lcarbonout, '  ', lcarbons
    
    print*,''
    print('(a)'),'--DEAD CARBON BALANCE (micromoles C/m2)--'
    print('(a,e14.6)'),'  Balance Error:  ',dcarbonb
    print('(a)'),       '  In:   Grazing      Harvest      Trans_Dead   Trans_Live'
    print('(a,e14.6,a,e14.6,a,e14.6,a,e14.6)'),'      ', &
          dgaingrz, '  ', dgainhrv, ' ', dgaintd, ' ', dgaintl
    print('(a)'),       '  Out:   Loss_Fire   Loss_Resp    Tran_Dead'
    print('(a,3(e14.6,a))'),        '      ', &
          dlossf, ' ', dlossr, ' ', dlosstd
    print('(a)'),       '  Check:   Init      End         Change'
    print('(a,3(e14.6,a))'),        '      ', &
          dpoolinit, ' ', dpoolend, ' ', dpoolchange
    print('(a)'),       '  Net:   In          Out          Stored'
    print('(a,e14.6,a,e14.6,a,e14.6)'),        '      ', &
          dcarbonin, '  ', dcarbonout, '  ', dcarbons

    print*,''
    print('(a)'),'--TRANSFER CARBON BALANCE (micromoles C/m2)--'
    tcarbonb = (llosstd - dgaintl) + (dlosstd - dgaintd)  

    print('(a,e14.6)'),'  Balance Error:  ', tcarbonb
    print('(a)'),'  Live->Dead:  Diff          LTran_Dead     Gain_DeadL'
    print('(a,e14.6,2(a,e14.6))'),    '              ', &
              llosstd-dgaintl, '  ', llosstd, '  ', dgaintl
    print('(a)'),'  TLive->Dead: Tot Diff      Tot_LLoss      Tot_DGain'
    print('(a,e14.6,2(a,e14.6))'),     '              ', &
              tlivedead - gdeadlive, '', tlivedead, '', gdeadlive
    print('(a)'),'  Dead->Dead:  Diff          DTran_Dead     Gain_DeadD'
    print('(a,e14.6,2(a,e14.6))'),    '              ', &
              dlosstd-dgaintd, '  ', dlosstd, '  ', dgaintd


    IF (abs(hrvcarbonb) .GT. dzero) THEN
       print*,''
       print('(a)'),'--HARVEST CARBON BALANCE (micromoles C/m2)--'
       print('(a,e14.6)'),'  Balance Error:  ', hrvcarbonb
       print('(a)'),'  Live_Loss     Dead_Gain   Resp_Hrvst   Rmvd_Hrvst'
       print('(a,4(e14.6,a))'),'       ', &
             llosshrv, dgainhrv, hrvresp, hrvrmvd
    ENDIF

    print('(a)'),''
    print('(a,e14.6)'),'==>Maximum Carbon Imbalance: ', netcarbonb
             

    if (cb_err) then
       print('(a)'),''
       print('(a)'),'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print('(a)'),'  !!!Error in Carbon Balance!!!'
       print('(a)'),'  !!!       Stopping.       !!!'
       print('(a)'),'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print('(a)'),''
    endif

    print('(a)'),'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    if (cb_err .and. carbonb_stop) stop
endif



if (loc_printoutc13) then
    print('(a)'),'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    print('(a)'), 'TIME-STEP CARBON C13 CYCLE BALANCE'
    print('(a,i6,2f10.2,i4)'),  &
        'Point/Lon/Lat/PFT: ',sibpt, siblon,siblat,pref
    print('(a,a,i3,a,i4)'), &
        'Date: ', trim(month_names(month)), day, ', ', year
    print('(a,f12.6,a,f12.6)'),'LAI/FPAR: ', laiin, '   ', fparin

    print('(a)'),''
    print('(a)'),'--CARBON FLUXES (micromoles C/m2/s)--'
    print('(a,f12.6)'),'  In:  ',fassimc13
    print('(a,f12.6)'),'  Out: ',frespc13
    print('(a)'),      '       RespC13        Saved        Calculated'
    print('(a,f12.6,a,f12.6)'), &
                       '       AutoC13 ', &
          poollt%resp_autoc13*mol_to_umol, '  ', respautoc13
    print('(a,f12.6,a,f12.6)'), &
                       '       HetC13  ', &
          pooldt%resp_hetc13*mol_to_umol, '  ', resphetc13
    print('(a,f12.6,a,f12.6)'), &
                       '       NVegC13 ', &
          poollt%resp_nvegc13*mol_to_umol, '  ', poollt%resp_nvegc13*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       GrazeC13', &
          poollt%resp_grzc13*mol_to_umol, '  ', poollt%resp_grzc13*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       HrvstC13', &
          poollt%resp_hrvstc13*mol_to_umol, '  ',poollt%resp_hrvstc13*mol_to_umol
    print('(a,f12.6,a,f12.6)'), &
                       '       GrwthC13', &
          poollt%resp_growc13*mol_to_umol, '  ', respgrowc13
    print('(a,f12.6,a,f12.6)'), &
                       '       LeafC13 ', &
          poollt%resp_leafc13*mol_to_umol, '  ', respleafc13
    print('(a,f12.6,a,f12.6)'), &
                       '       MntnC13 ', &
          poollt%resp_mntnc13*mol_to_umol, '  ', respmntnc13
    print('(a,f12.6,a,f12.6)'), &
                       '       RootsC13', &
          poollt%resp_rootc13*mol_to_umol, '  ', resprootc13
    print('(a,f12.6,a,f12.6)'), &
                       '       SoilC13 ', &
          pooldt%resp_soilc13*mol_to_umol, '  ', respsoilc13
   print('(a,f12.6,a)'),'  NEEC13: ',frespc13-fassimc13, '  (Out-In)'


    print('(a)'),''
    print('(a)'),'--CARBON C13 POOLS (Mg C/ha)--'
    print('(a,f14.8,a,e14.6,f14.8)'),'  Live C13 Carbon Init/Change/End: ', &
            lpoolinitc13, '   ', lpoolchangec13, lpoolendc13
    print('(a,f14.8,a,e14.6,f14.8)'),'  Dead C13 Carbon Init/Change/End: ', &
            dpoolinitc13, '   ', dpoolchangec13, dpoolendc13
    print('(a,f14.8,a,e14.6,f14.8)'),'  Live Carbon Init/Change/End: ', &
            lpoolinit, '   ', lpoolchange, lpoolend
    print('(a,f14.8,a,e14.6,f14.8)'),'  Dead Carbon Init/Change/End: ', &
            dpoolinit, '   ', dpoolchange, dpoolend

    print('(a)'),''
    print('(a)'),'  Pool     Orig            Gain            Loss              New'   
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Leaf  ', &
            (poollt%poolpft(lpc13) - sum(poollt%poolpft_dgain(lpc13,:)) &
               + sum(poollt%poolpft_dloss(lpc13,:)))*molc13_to_mg, &
               '  ', sum(poollt%poolpft_dgain(lpc13,:))*molc13_to_mg, &
               '  ', sum(poollt%poolpft_dloss(lpc13,:))*molc13_to_mg, &
               '  ', poollt%poolpft(lpc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  FRoot ', &
            (poollt%poolpft(frpc13)-sum(poollt%poolpft_dgain(frpc13,:)) &
              + sum(poollt%poolpft_dloss(frpc13,:)))*molc13_to_mg, &
              '  ', sum(poollt%poolpft_dgain(frpc13,:))*molc13_to_mg, &
              '  ', sum(poollt%poolpft_dloss(frpc13,:))*molc13_to_mg, &
              '  ', poollt%poolpft(frpc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  CRoot ', &
            (poollt%poolpft(crpc13)-sum(poollt%poolpft_dgain(crpc13,:)) &
               + sum(poollt%poolpft_dloss(crpc13,:)))*molc13_to_mg, &
               '  ', sum(poollt%poolpft_dgain(crpc13,:))*molc13_to_mg, &
               '  ', sum(poollt%poolpft_dloss(crpc13,:))*molc13_to_mg, &
               '  ', poollt%poolpft(crpc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  StWd  ', &
            (poollt%poolpft(wpc13)-sum(poollt%poolpft_dgain(wpc13,:)) &
              + sum(poollt%poolpft_dloss(wpc13,:)))*molc13_to_mg, &
              '  ', sum(poollt%poolpft_dgain(wpc13,:))*molc13_to_mg, &
              '  ', sum(poollt%poolpft_dloss(wpc13,:))*molc13_to_mg, &
              '  ', poollt%poolpft(wpc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Prod  ', &
            (poollt%poolpft(ppc13)-sum(poollt%poolpft_dgain(ppc13,:)) &
              + sum(poollt%poolpft_dloss(ppc13,:)))*molc13_to_mg, &
              '  ', sum(poollt%poolpft_dgain(ppc13,:))*molc13_to_mg, &
              '  ', sum(poollt%poolpft_dloss(ppc13,:))*molc13_to_mg, &
              '  ', poollt%poolpft(ppc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  CDB   ', &
            (pooldt%poollu(cdbc13)-sum(pooldt%poollu_dgain(cdbc13,:)) &
               + sum(pooldt%poollu_dloss(cdbc13,:)))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dgain(cdbc13,:))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dloss(cdbc13,:))*molc13_to_mg, &
              '  ', pooldt%poollu(cdbc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  MetL  ', &
            (pooldt%poollu(metlc13)-sum(pooldt%poollu_dgain(metlc13,:)) &
              + sum(pooldt%poollu_dloss(metlc13,:)))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dgain(metlc13,:))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dloss(metlc13,:))*molc13_to_mg, &
              '  ', pooldt%poollu(metlc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  StrL  ', &
            (pooldt%poollu(strlc13)-sum(pooldt%poollu_dgain(strlc13,:)) &
               + sum(pooldt%poollu_dloss(strlc13,:)))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dgain(strlc13,:))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dloss(strlc13,:))*molc13_to_mg, &
              '  ', pooldt%poollu(strlc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  SLit  ', &
            (pooldt%poollu(slitc13)-sum(pooldt%poollu_dgain(slitc13,:)) &
               + sum(pooldt%poollu_dloss(slitc13,:)))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dgain(slitc13,:))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dloss(slitc13,:))*molc13_to_mg, &
              '  ', pooldt%poollu(slitc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Slow  ', &
            (pooldt%poollu(slowc13)-sum(pooldt%poollu_dgain(slowc13,:)) &
               + sum(pooldt%poollu_dloss(slowc13,:)))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dgain(slowc13,:))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dloss(slowc13,:))*molc13_to_mg, &
              '  ', pooldt%poollu(slowc13)*molc13_to_mg
    print('(a,f14.9,a,e14.6,a,e14.6,a,f14.9)'), &
            '  Arm   ', &
            (pooldt%poollu(armc13)-sum(pooldt%poollu_dgain(armc13,:)) &
               + sum(pooldt%poollu_dloss(armc13,:)))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dgain(armc13,:))*molc13_to_mg, &
              '  ', sum(pooldt%poollu_dloss(armc13,:))*molc13_to_mg, &
              '  ', pooldt%poollu(armc13)*molc13_to_mg

    print('(a)'),''
    print('(a)'),'--LIVE C13 CARBON BALANCE (micromoles C/m2)--'
    print('(a,e14.6)'),'  Balance Error:  ',lcarbonbc13
    if (seed_printout) print('(a)'),'  Non-Minimal Error Due To Crop Seed Release'
    print('(a)'),       '  Gain:      Assim          Seed_Release'
    print('(a,e14.6,a,e14.6)'),         '           ', &
         lgainac13, ' ', lgainsc13
    print('(a)'),       '  Loss:      Fire           Harvest        Tran_Dead'     
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
         llossfirec13, ' ', llosshrvc13, ' ', llosstdc13
    print('(a)'),       '  Resp:      Growth         Maintenance    Non-Veg' 
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
         llossrgc13, ' ', llossrmc13, ' ', llossrnvegc13
    IF (abs(grzcarbonbc13) .GT. dzero) THEN
         print('(a)'), ' Graze:      Loss           Resp           Tran_Dead'
         print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
               llossgrzc13, '  ', llossrgrzc13, '  ', llosstgrzc13
    ENDIF
    print('(a)'),       '  Net:       In             Out            Stored'            
    print('(a,e14.6,a,e14.6,a,e14.6)'), '           ', &
          lcarboninc13, '   ', lcarbonoutc13, '  ', lcarbonsc13
    
    print*,''
    print('(a)'),'--DEAD C13 CARBON BALANCE (micromoles C/m2)--'
    print('(a,e14.6)'),'  Balance Error:  ',dcarbonbc13
    print('(a)'),       '  In:   Grazing      Harvest      Trans_Dead   Trans_Live'
    print('(a,e14.6,a,e14.6,a,e14.6,a,e14.6)'),'      ', &
          dgaingrzc13, '  ', dgainhrvc13, ' ', dgaintdc13, ' ', dgaintlc13
    print('(a)'),       '  Out:   Loss_Fire   Loss_Resp    Tran_Dead'
    print('(a,3(e14.6,a))'),        '      ', &
          dlossfc13, ' ', dlossrc13, ' ', dlosstdc13
    print('(a)'),       '  Check:   Init      End         Change'
    print('(a,3(e14.6,a))'),        '      ', &
          dpoolinitc13, ' ', dpoolendc13, ' ', dpoolchangec13
    print('(a)'),       '  Net:   In          Out          Stored'
    print('(a,e14.6,a,e14.6,a,e14.6)'),        '      ', &
          dcarboninc13, '  ', dcarbonoutc13, '  ', dcarbonsc13

    print*,''
    print('(a)'),'--TRANSFER C13 CARBON BALANCE (micromoles C/m2)--'
    tcarbonbc13 = (llosstdc13 - dgaintlc13) + (dlosstdc13 - dgaintdc13)

    print('(a,e14.6)'),'  Balance Error:  ', tcarbonbc13
    print('(a)'),'  Live->Dead:  Diff          LTran_Dead     Gain_DeadL'
    print('(a,e14.6,2(a,e14.6))'),    '              ', &
              llosstdc13-dgaintlc13, '  ', llosstdc13, '  ', dgaintlc13
    print('(a)'),'  TLive->Dead: Tot Diff      Tot_LLoss      Tot_DGain'
    print('(a,e14.6,2(a,e14.6))'),     '              ', &
              tlivedeadc13 - gdeadlivec13, '', tlivedeadc13, '', gdeadlivec13
    print('(a)'),'  Dead->Dead:  Diff          DTran_Dead     Gain_DeadD'
    print('(a,e14.6,2(a,e14.6))'),    '              ', &
              dlosstdc13-dgaintdc13, '  ', dlosstdc13, '  ', dgaintdc13


    IF (abs(hrvcarbonbc13) .GT. dzero) THEN
       print*,''
       print('(a)'),'--HARVEST C13 CARBON BALANCE (micromoles C/m2)--'
       print('(a,e14.6)'),'  Balance Error:  ', hrvcarbonbc13
       print('(a)'),'  Live_Loss     Dead_Gain   Resp_Hrvst   Rmvd_Hrvst'
       print('(a,4(e14.6,a))'),'       ', &
             llosshrvc13, dgainhrvc13, hrvrespc13, hrvrmvdc13
    ENDIF

    print('(a)'),''
    print('(a,e14.6)'),'==>Maximum Carbon Imbalance: ', netcarbonbc13


    if (cb_err) then
       print('(a)'),''
       print('(a)'),'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print('(a)'),'  !!!Error in Carbon C13 Balance!!!'
       print('(a)'),'  !!!!!       Stopping.       !!!!!'
       print('(a)'),'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print('(a)'),''
    endif

    print('(a)'),'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    if (cb_errc13 .and. carbonb_stop) stop
endif



!...update variables
pooldt%poollup(:) = pooldt%poollu(:) 
poollt%poolpftp(:) = poollt%poolpft(:)                              
        
    
end subroutine balan_carbon    
