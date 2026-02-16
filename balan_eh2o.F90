!...calculation to determine if water and energy
!...balance is maintained
subroutine balan_eh2o(gref, glon, glat, pref, &
                 lspr, cupr, lai, &
                 fluxt, sscolt, hydrost, &
                 ebalnum, wbalnum)

use kinds
use module_local
use module_pftinfo, only: &
     group_crop, pft_group, pft_num
use module_pparams, only: &
    lvap, month_names, secs_per_day, snofac
use module_sib, only: &
    flux_type, hydros_type, &
    sscol_type
use module_sibconst, only: &
    nsoil, nsnow, &
    bnum_allow,   &
    energyb_print, energyb_stop, &
    energyb_thresh, &
    waterb_print, waterb_stop, &
    waterb_thresh
use module_time, only: &
    dtsib, dtisib, &
    day, month, year, &
    sec_tot, init_second

implicit none

!-----------------------------------------------------------------------
!...WATER BALANCE...
!
!     WATER IN + CHANGE IN STORAGE = WATER OUT
!
!     inputs:
!       precipitation (cupr,lspr)    (kg/m2)
!     outputs:
!       CAS-mixed layer latent heat flux (fws)   (W/m2)
!       runoff
!          surface runoff       (roffo)       (mm)
!          subsurface runoff    (roff)        (mm)
!
!     storage:
!          liquid interception reservoirs
!             canopy interception   (capacc_liq)  (kg/m2)
!             ground interception   (capacg)      (kg/m2)
!          soil
!             soil water            (www_liq)     (kg/m2)
!             soil ice              (www_ice)     (kg/m2)
!          snow
!             snow water            (www_liq)     (kg/m2)
!             snow ice              (www_ice)     (kg/m2)
!
!          ** note: indices for soil and ice water
!             will be 1 to 10 for soil (1=top, 10=bottom)
!                    -4 to 0 for snow (-4=top, 0=surface)
!
!----------------------------------------------------------------------
!...ENERGY BALANCE
!
!    ENERGY IN + CHANGE IN STORAGE = ENERGY OUT
!
!    inputs:
!      radiation absorbed by canopy    (radtc)   (W/m2)
!      radiation absorbed by ground    (radtg)   (W/m2)
!
!    outputs:   
!      latent heat 
!       transpiration                       (ect)      (J/m2)
!       evaporation
!          from canopy interception storage (eci)      (J/m2)
!          from ground interception storage (egi)      (J/m2)
!          from surface soil layer          (egs)      (J/m2)
!          from snow                        (es)       (J/m2)  
!      sensible heat
!          canopy                           (hc)       (J/m2)
!          ground                           (hg)       (J/m2)
!          snow                             (hs)       (J/m2)
!
!    storage
!      canopy storage flux                  (storhc)      (W/m2)
!      ground storage flux                  (storhg)      (W/m2)
!
!
!----------------------------------------------------------------------
!...input variables
integer(i4), intent(in) :: gref, pref
real(r4), intent(in) :: glon, glat
real(r8), intent(in) :: lspr, cupr, lai
type(flux_type), intent(in) :: fluxt
type(sscol_type), intent(in) :: sscolt
type(hydros_type), intent(inout) :: hydrost
integer(i4), intent(inout) :: ebalnum, wbalnum

!...water balance variables
real(r8) :: wbal     ! water balance (kg/m2)
real(r8) :: wbalalt  ! alternate water balance (kg/m2)

!.....storage
real(r8) :: dstor    ! net change in storage (kg/m2)
                     !   --canopy and ground interception
                     !   --soil column
                     !   --snow columh
real(r8) :: dqint    ! change in storage, 
                     !     canopy/ground interception (kg/m2)
real(r8) :: dqsoil   ! change in storage, soil (kg/m2)
real(r8) :: dqsnow   ! change in storage, snow (kg/m2)

!.....cas
real(r8) :: dcas       ! change in CAS from latent heat flux (kg/m2)
real(r8) :: dcas_eci   ! canopy interception (kg/m2)
real(r8) :: dcas_ect   ! canopy transpiration (kg/m2)
real(r8) :: dcas_egi   ! ground interception (kg/m2)
real(r8) :: dcas_egs   ! ground evaporation (kg/m2)
real(r8) :: dcas_es    ! snow evaporation (kg/m2)
real(r8) :: dcas_stor  ! storage change (kg/m2)
real(r8) :: dcas_sum   ! sum of components (kg/m2)

!.....runoff
real(r8) :: runoff   ! net runoff (kg/m2)
                     !   --surface + subsurface

!.....input
real(r8) :: precip   ! precip (kg/m2)

!.....holders
real(r8) :: sbegi,sendi ! beginning and ending ice (kg/m2)
real(r8) :: sbegl,sendl ! beginning and ending liquid (kg/m2)

!...energy balance variables
real(r8) :: casbal   ! canopy-air-space (CAS) balance (W/m2)
real(r8) :: cbal     ! canopy balance (W/m2)
real(r8) :: gbal     ! ground balance (W/m2)
real(r8) :: ebal     ! energy balance (W/m2)
real(r8) :: ein      ! energy in (W/m2)
real(r8) :: edstor   ! energy stored (W/m2)
real(r8) :: eout     ! energy out (W/m2)
real(r8) :: etot, htot

!...local variables
integer(i4) :: i        ! loop index
logical :: loc_printout
logical :: eb_err, wb_err

!---------------------------------------------
!----- ENERGY BALANCE....units of W/m2-----

    ein = (fluxt%hg + fluxt%hc + fluxt%hs + &
          fluxt%ect + fluxt%eci + fluxt%egi + &
          fluxt%egs  + fluxt%es) * dtisib
    edstor = cas_w_storage + &  ! change in CAS water storage
          cas_e_storage         ! change in CAS heat storage

    eout = fluxt%fss       +      &  ! sensible heat flux
           fluxt%fws                 ! latent heat flux

    casbal = ein - edstor - eout
    cbal = radttc - fluxt%storhc - &
         ( fluxt%ect + fluxt%eci + fluxt%hc ) * dtisib
    gbal = radttg +  radtts - fluxt%storhg &
         - (fluxt%hg + fluxt%hs + fluxt%egi + fluxt%egs + fluxt%es ) * dtisib

    ein = radttc + radttg + radtts
    edstor = cas_e_storage + cas_w_storage + &
             fluxt%storhc + fluxt%storhg 
    eout  = fluxt%fss + fluxt%fws

    ebal  = ein - edstor - eout

    etot = (fluxt%ect+fluxt%eci+fluxt%egi+fluxt%egs+fluxt%es)*dtisib
    htot = (fluxt%hc + fluxt%hg + fluxt%hs)*dtisib

    !...print if requested
     eb_err = .false.
    if ((abs(ebal) > energyb_thresh) .and. &
        (sec_tot > init_second)) then
         if (ebalnum .lt. bnum_allow) then
             ebalnum = ebalnum + 1
         else
             eb_err = .true.
         endif
    else
         ebalnum = 0
    endif
    loc_printout = eb_err .or. energyb_print

     if (loc_printout) then

           print('(a)'),'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
           print'(2(a,g16.6))','ENERGY IMBALANCE    : ',ebal,' kg/m2'
           print('(a,a,i3,a,i4)'), &
                 'Date: ', trim(month_names(month)), day, ', ', year
           print('(a,i12)'),    'Second =',sec_tot
           print('(a,i8,2f10.2,i4)'),  &
                   'Point/lon/lat/PFT: ',gref, glon, glat, pref
           print('(a,i6,g14.6)'), &
                   'Snow Layers and Ground Cover: ',-sscolt%nsl, hydrost%snow_gvfc
           print('(a)'),''
           print('(a,g14.6,a)'),'canopy air space balance=',casbal,' (W/m2)'
           print('(a,g14.6,a)'),'canopy balance=',cbal,' (W/m2)'
           print('(a,g14.6,a)'),'ground balance=',gbal,' (W/m2)'
           print('(a)'),''
           print('(a)'),'TERMS:  ENERGY IN = STORAGE CHANGE + ENERGY OUT'
           print'(3(a,g14.6))', ' in=',ein,' dstor=',edstor,' out=',eout
           print('(a)'),''
           print'(a,g14.6)','IN: ',radttc+radttg+radtts
           print'(3(a,g14.6))',' radttc=',radttc,' radttg=', &
               radttg, ' radtts=',radtts
           print'(a,g14.6)','DSTOR: ', cas_w_storage + cas_e_storage + &
                             fluxt%storhc + fluxt%storhg
           print'(2(a,g14.6))',' w_storage=',cas_w_storage,' e_storage=',cas_e_storage
           print'(2(a,g14.6))',' storhc=',fluxt%storhc,' storhg=',fluxt%storhg
           print'(1(a,g14.6))','NET OUT: ',eout
           print'(2(a,g14.6))',' fss=',fluxt%fss,' fws=',fluxt%fws
           print'(2(a,g14.6))',' ect=',fluxt%ect*dtisib,' eci=',fluxt%eci*dtisib
           print'(2(a,g14.6))',' egi=',fluxt%egi*dtisib,' egs=',fluxt%egs*dtisib
           print'(2(a,g14.6))',' es=',fluxt%es*dtisib, ' etot=',etot
           print'(2(a,g14.6))',' hc=',fluxt%hc*dtisib, ' hg=', fluxt%hg*dtisib
           print'(2(a,g14.6))',' hs=',fluxt%hs*dtisib, ' htot=',htot
           print'(a)','EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'

           if ((energyb_stop) .and. (eb_err)) stop

    endif  !energy_balance_printout


!-----------------------------------------------------------------------
!-----WATER BALANCE....units of kg/m2 water-----

    !Storage Terms
    !...change in canopy and surface interception storage 
    dqint = (hydrost%capacc_liq - capaccl_old) +  & ! canopy liquid interception
            (hydrost%capacc_snow - capaccs_old) + & ! canopy snow interception
            (hydrost%capacg - capacg_old)           ! surface interception

    !...change in soil water
    dqsoil = 0.0
    do i=1,nsoil
        dqsoil = dqsoil + (sscolt%www_liq(i) - wwwliq_old(i)) + &
            (sscolt%www_ice(i) - wwwice_old(i))
    enddo

    !...change in snow mass
    dqsnow = hydrost%snow_gmass - snow_gmass_old

    !...net change in storage
    dstor = dqint + dqsoil + dqsnow

    !CAS Terms
    dcas = fluxt%fws * dtsib / lvap  ! net CAS water vapor flux 
    dcas_eci = fluxt%eci / lvap  ! canopy interception
    dcas_ect = fluxt%ect / lvap  ! canopy transpiration
    dcas_egi = fluxt%egi / lvap  ! ground interception
    dcas_egs = fluxt%egs / lvap  ! ground evaporation 
    dcas_es  = fluxt%es*snofac / lvap  ! snow evaporation
    dcas_stor = cas_w_storage * dtsib / lvap  
                                  ! cas storage change
    dcas_sum = dcas_eci + dcas_ect + dcas_egi + dcas_egs + &
               dcas_es

    !...runoff, surface + subsurface
    runoff = hydrost%roffo + hydrost%roff

    !...precip
    precip = (lspr + cupr) * dtsib

    !...WATER BALANCE...
    wbal =  precip - dcas - runoff - dstor - dcas_stor
    wbalalt = precip - dcas_sum - runoff - dstor - dcas_stor
    
    !!...stash a quick fix to take from runoff....
    if ((wbal .lt. dzero) .and. &
        (hydrost%roff .gt. dzero)) then
       hydrost%roff = MAX(dzero, hydrost%roff + wbal)
 
       runoff = hydrost%roffo + hydrost%roff
       wbal = precip - dcas - runoff - dstor - dcas_stor
    endif

    !!...determine total water balance
    wb_err = .false.
    if ((abs(wbal)    > waterb_thresh) .and. &
        (abs(wbalalt) > waterb_thresh) .and. &
        (sec_tot > init_second)) then
         if (wbalnum .lt. bnum_allow) then
             wbalnum = wbalnum + 1
         else
             wb_err = .true.
         endif
    else
         wbalnum = 0
    endif
    loc_printout = wb_err .or. waterb_print

    !!...allow unbalance for 'irrigated' crops
    if (wb_err) THEN
       IF (pft_group(pft_num(pref)) .EQ. group_crop) THEN
          wbalnum = wbalnum + 1
          !IF (wbalnum .lt. 3*bnum_allow) THEN
          wb_err = .false.
          loc_printout = .false.
          !ENDIF
       ENDIF
    endif
    
     !!...print if requested
     if (loc_printout) then
           print'(a)','WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW'
           print'(2(a,g16.6))','WATER IMBALANCE    : ',wbal,' kg/m2'
           print'(2(a,g16.6))','ALT WATER IMBALANCE: ',wbalalt, ' kg/m2'
           print('(a,a,i3,a,i4)'), &
                 'Date: ', trim(month_names(month)), day, ', ', year
           print'(a,i12)','Second=',sec_tot
           print'(a)',''
           print'(a,2f12.4,1i7,2i4)','Point (Lon/Lat/Ref/PFT): ', &
                 glon, glat, gref, pref
           print'(a,f8.4)',    'Current LAI: ',lai
           print'(a,i4,f8.4)', 'Snow Layers and Ground Cover: ', &
                    -sscolt%nsl, hydrost%snow_gvfc
           print'(a,2f8.4)', 'Wetness Fractions (Can/Grnd): ', &
                    hydrost%wetfracc, hydrost%wetfracg

           print'(a)',''
           print'(a)','TERMS:  precip = flux + runoff + storage (all units kg/m2)'
           print'(a,g14.6)','input: ',precip
           print'(a,g14.6)','   lspr =',lspr*dtsib
           print'(a,g14.6)','   cupr =',cupr*dtsib
           print'(a)',''
           print'(a,g14.6)','output: ',dcas+dstor+runoff+dcas_stor
           print'(a,g14.6)','   flux    =',dcas
           print'(a,g14.6)','   runoff  =',runoff
           print'(a,g14.6)','   storage =',dstor+dcas_stor   

           print'(a)','' 
           print'(a)','Flux Components:'
           print'(a,g14.6)', '  Net CAS Flux              = ',dcas
           print'(a,2g14.6)','  Canopy Sfc_Evap/Transpire = ',dcas_eci,dcas_ect
           print'(a,2g14.6)','  Ground Sfc_Evap/Soil_Evap = ',dcas_egi,dcas_egs
           print'(a,g14.6)', '  Snow Evaporate            = ',dcas_es
           print'(a,2g14.6)', '  Sum, Difference           = ',dcas_sum, dcas - dcas_sum

           print'(a)',''
           print'(a,g14.6)','Runoff Amounts: ',hydrost%roffo+hydrost%roff
           print'(a,g14.6)','   overland = ',hydrost%roffo
           print'(a,g14.6)','   subsfc   = ',hydrost%roff

           print'(a)',''
           print'(a,g14.6)','Storage Changes: ', dstor+dcas_stor
           print'(a,g14.6)', '  cas          = ', dcas_stor
           print'(a,g14.6)', '  interception = ', dqint
           print'(a,g14.6)', '  soil         = ', dqsoil
           print'(a,g14.6)', '  snow         = ', dqsnow

           print('(a)'),''
           print('(a)'),'   ->Interception'
           print'(2a,3G14.6)','     ', &
                    'capacc_liq  (old/new/diff): ',capaccl_old, &
                    hydrost%capacc_liq, hydrost%capacc_liq - capaccl_old
           print'(2a,3G14.6)','     ', &
                    'capacc_snow (old/new/diff): ',capaccs_old, &
                    hydrost%capacc_snow, hydrost%capacc_snow - capaccs_old
           print'(2a,3G14.6)','     ', &
                    'capacg      (old/new/diff): ',capacg_old, &
                    hydrost%capacg, hydrost%capacg - capacg_old

           print'(a)',''
           print'(a)','  ->Soil Moisture'
           sbegl = 0.0
           sendl = 0.0
           sbegi = 0.0
           sendi = 0.0
           print'(2a)','      Lev  Beg_Ice       End_Ice  ', &
                       '    Beg_Liq       End_Liq       TD'
           do i=1,nsoil
              print'(a,i4,5g14.6)','    ', &
                   i, wwwice_old(i), sscolt%www_ice(i), &
                   wwwliq_old(i),sscolt%www_liq(i), &
                   sscolt%td(i)

               sbegl = sbegl + wwwliq_old(i)
               sendl = sendl + sscolt%www_liq(i)
               sbegi = sbegi + wwwice_old(i)
               sendi = sendi + sscolt%www_ice(i)
           enddo
           print'(a,4g14.6)','  Total:', &
                sbegi,sendi,sbegl,sendl
           print'(2a,g14.6)','  ', &
                'Total Ice Change=',sendi-sbegi
           print'(2a,g14.6)','  ', &
                'Total Liq Change=',sendl-sbegl
           print'(2a,g14.6)','  ', &
                'Total H2O Change=',(sendi+sendl)-(sbegi+sbegl)

           if ((nsl_old < 0) .or. &
               (sscolt%nsl < 0)) then
               print'(a)',''
               print'(a)','  ->Snow'
               print'(a,2(a,i3))','     ','  nsl_old=',nsl_old, ' nsl=',sscolt%nsl
               print'(a,2(a,g12.4))','     ', &
                     ' gmass_old=',snow_gmass_old, &
                     ' gmass=',hydrost%snow_gmass

               sbegi = 0.0
               sendi = 0.0
               sbegl = 0.0
               sendl = 0.0
               do i=-nsnow+1,0
                   print'(a,i4,5(a,g14.6))','     ', &
                         i,' beg ice=',wwwice_old(i), &
                         ' end ice=' ,sscolt%www_ice(i), &
                         ' beg liq=',wwwliq_old(i),' end liq=' ,sscolt%www_liq(i), &
                         ' td=',sscolt%td(i)
                   sbegi = sbegi+ wwwice_old(i)
                   sendi = sendi + sscolt%www_ice(i)
                   sbegl = sbegl + wwwliq_old(i)
                   sendl = sendl + sscolt%www_liq(i)
               enddo

               print'(a,3(a,g14.6))','     ', &
                     '   Total Beg Ice=',sbegi,' Total End Ice=',sendi, &
                     ' Total Ice Change=',sendi-sbegi
               print'(a,3(a,g14.6))','     ', &
                    '   Total Beg Liq=',sbegl,' Total End Liq=',sendl, &
                    ' Total Liq Change=',sendl-sbegl
               print'(a,3(a,g14.6))','     ', &
                    '   Total Beg H2O=',sbegl+sbegi, &
                    ' Total End H2O=',sendl+sendi, &
                    ' Total H2O Change=',sendl-sbegl+sendi-sbegi
            else
               print'(a)',''
               print'(a)','   ->No Snow.'
            endif

            print'(a)',''
            print'(a)','WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW'
            print'(a)',''

            if ((waterb_stop) .and. (wb_err)) stop

    endif  !water_balance_printout


end subroutine balan_eh2o
