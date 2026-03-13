!==========================================================================
subroutine set_continue_spinup()
!==========================================================================
! Resets the pools for continuing spinup runs.

use kinds
use module_sibconst, only: &
    subcount, spinup_done, &
    npoolpft, npoollu, &
    spinup_continue
use module_io, only: requib_writef
use module_pftinfo, only: pft_num
use module_poolinfo, only: pool_indx_lay
use module_sib, only: sib

implicit none

!...local variables
integer :: i,l,n,s
integer(i4) :: pref, pnum

!... call and run continue version of equipools_calc()
call equipools_calc_continue()

!...Reset equilibrium switch
requib_writef = .false.

!...Set the pools for the next spin-up iteration
if (.not. spinup_done) then

   do i=1,subcount
      if (.not. sib%g(i)%gdiagt%gridcell_spunup) then
          do l=1,sib%g(i)%g_nlu
             if (.not. sib%g(i)%l(l)%equibdt%lupft_spunup) then
                pref = sib%g(i)%l(l)%ipft
                pnum = pft_num(pref)

                call equipools_restart_continue(pnum, &
                    sib%g(i)%l(l)%equiblt%poolpft_end,   &
                    sib%g(i)%l(l)%equiblt%poolpft_equib, &
                    sib%g(i)%l(l)%poollt%poolpft_flay,   &
                    sib%g(i)%l(l)%poollt%poolpft_lay)

                !...reset live pools
                do n=1,npoolpft
                   sib%g(i)%l(l)%equiblt%poolpft_init(n) = &
                       sum(sib%g(i)%l(l)%poollt%poolpft_lay(n,:))
                   sib%g(i)%l(l)%poollt%poolpft(n) = &
                       sib%g(i)%l(l)%equiblt%poolpft_init(n)
                   sib%g(i)%l(l)%equiblt%poolpft_min(n) = &
                       sib%g(i)%l(l)%equiblt%poolpft_init(n)
                   sib%g(i)%l(l)%equiblt%poolpft_max(n) = &
                        sib%g(i)%l(l)%equiblt%poolpft_init(n)

                   sib%g(i)%l(l)%equiblt%poolpft_totgain(n) = dzero
                   sib%g(i)%l(l)%equiblt%poolpft_totloss(n) = dzero

               enddo !n=1,npoolpft
               sib%g(i)%l(l)%poollt%poolpftp(:) = sib%g(i)%l(l)%poollt%poolpft(:)

               !...reset dead pools
               do n=1,npoollu/3 !1,6 (npoollu=18 with C14)
                  sib%g(i)%l(l)%equibdt%poollu_init(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%pooldt%poollu(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_min(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_max(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)

                  do s=1,pool_indx_lay(n+npoolpft/3) !(6,11)->pool_indx_lay either 1 or 10
                  !npoolpft is 15 with C14
                  !pool_indx_lay(ntpool), poollu_lay(npoollu,nsoil)
                      sib%g(i)%l(l)%pooldt%poollu_lay(n,s) = &
                              sib%g(i)%l(l)%equibdt%poollu_equib(n) * &
                              sib%g(i)%l(l)%pooldt%poollu_flay(n,s)
                   enddo
                   sib%g(i)%l(l)%equibdt%poollu_totgain(n) = dzero
                   sib%g(i)%l(l)%equibdt%poollu_totloss(n) = dzero
               enddo !n=1,npoollu/2

               !..same as above but for C13 pools
               do n=npoollu/3+1,2*(npoollu/3) !7,12 (npoollu=18 with C14)
                  sib%g(i)%l(l)%equibdt%poollu_init(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%pooldt%poollu(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_min(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_max(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)

                  do s=1,pool_indx_lay(n+(2*(npoolpft/3))) !(17,22)->pool_indx_lay either 1 or 10
                      sib%g(i)%l(l)%pooldt%poollu_lay(n,s) = &
                              sib%g(i)%l(l)%equibdt%poollu_equib(n) * &
                              sib%g(i)%l(l)%pooldt%poollu_flay(n,s)
                   enddo
                   sib%g(i)%l(l)%equibdt%poollu_totgain(n) = dzero
                   sib%g(i)%l(l)%equibdt%poollu_totloss(n) = dzero
               enddo !n=1,npoollu

               !..same as above but for C14 pools
               do n=2*(npoollu/3)+1,npoollu !13,18
                  sib%g(i)%l(l)%equibdt%poollu_init(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%pooldt%poollu(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_min(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_max(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)

                  do s=1,pool_indx_lay(n+npoolpft) !(28,33)->pool_indx_lay either 1 or 10
                      sib%g(i)%l(l)%pooldt%poollu_lay(n,s) = &
                              sib%g(i)%l(l)%equibdt%poollu_equib(n) * &
                              sib%g(i)%l(l)%pooldt%poollu_flay(n,s)
                   enddo
                   sib%g(i)%l(l)%equibdt%poollu_totgain(n) = dzero
                   sib%g(i)%l(l)%equibdt%poollu_totloss(n) = dzero
               enddo !n=1,npoollu

               sib%g(i)%l(l)%pooldt%poollup(:) = sib%g(i)%l(l)%pooldt%poollu(:)

             endif  !.not. lu_spunup
          enddo  !l=1,g_nlu

      endif  !.not. gridcell_spunup
   enddo  !subcount
endif  !spinup_done

end subroutine set_continue_spinup



!==========================================================================
subroutine equipools_calc_continue()
!==========================================================================
! Calculates the quasi-equilibrium carbon pools.

use kinds
use module_param, only: poolcon
use module_pftinfo, only: &
    pft_num, pft_type, type_bare
use module_poolinfo
use module_sib, only: sib
use module_sibconst, only: &
    single_pt, subcount,   &
    npoolpft, npoollu, &
    spinup_threshold, spinup_done, &
    spinup_continue

implicit none

!...local variables
real(r8) :: ave_gain    !(mol/m2/s) average external inputs per pool
real(r8) :: ave_loss    !(mol/m2/s) average external outputs per pool
real(r8) :: ave_k_rate  !(1/s) average scaled decay rate constant
real(r8) :: pool_init, pool_end, init_ratio, end_ratio
real(r8) :: pdiffr, pdiffi, pdiffe
real(r8) :: pdiffrc13, pdiffrc14

!...misc variables
integer(byte) :: ptype
integer(i4) :: i,l,n
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: lpc13,frpc13,crpc13,wpc13,ppc13
integer(i4) :: lpc14,frpc14,crpc14,wpc14,ppc14
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: cdbpc13, metlpc13, strlpc13, slitpc13, slowpc13, armpc13
integer(i4) :: cdbpc14, metlpc14, strlpc14, slitpc14, slowpc14, armpc14
integer(i4) :: pref, pnum

!-------Calculatae quasi-equlibrium pools------------
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

cdbp  = pool_indx_cdb-npoolpft/3 !ntpool index 6, npoollu index 1
metlp = pool_indx_metl-npoolpft/3 !ntpool index 7, npoollu index 2
strlp = pool_indx_strl-npoolpft/3 !ntpool index 8, npoollu index 3
slitp = pool_indx_slit-npoolpft/3 !ntpool index 9, npoollu index 4
slowp = pool_indx_slow-npoolpft/3 !ntpool index 10, npoollu index 5
armp  = pool_indx_arm-npoolpft/3 !ntpool index 11, npoollu index 6

cdbpc13  = pool_indx_cdb_c13 - 2*npoolpft/3 !ntpool index 17, npoollu index 7
metlpc13  = pool_indx_metl_c13 - 2*npoolpft/3 !ntpool index 18, npoollu index 8
strlpc13  = pool_indx_strl_c13 - 2*npoolpft/3 !ntpool index 19, npoollu index 9
slitpc13  = pool_indx_slit_c13 - 2*npoolpft/3 !ntpool index 20, npoollu index 10
slowpc13 = pool_indx_slow_c13 - 2*npoolpft/3 !ntpool index 21, npoollu index 11
armpc13  = pool_indx_arm_c13 - 2*npoolpft/3 !ntpool index 22, npoollu index 12

cdbpc14  = pool_indx_cdb_c14 - npoolpft !ntpool index 28, npoollu index 13
metlpc14  = pool_indx_metl_c14 - npoolpft !ntpool index 29, npoollu index 14
strlpc14  = pool_indx_strl_c14 - npoolpft !ntpool index 30, npoollu index 15
slitpc14  = pool_indx_slit_c14 - npoolpft !ntpool index 31, npoollu index 16
slowpc14 = pool_indx_slow_c14 - npoolpft !ntpool index 32, npoollu index 17
armpc14 = pool_indx_arm_c14 - npoolpft !ntpool index 33, npoollu index 18

spinup_done = .true.

!...loop through gridcell points
do i=1,subcount

   if (.not. sib%g(i)%gdiagt%gridcell_spunup) then
      sib%g(i)%gdiagt%gridcell_spunup = .true.

      !...loop through landunits
      do l=1,sib%g(i)%g_nlu
         if (.not. sib%g(i)%l(l)%equibdt%lupft_spunup) then

             sib%g(i)%l(l)%equibdt%lupft_spunup = .true.
             pref=sib%g(i)%l(l)%ipft
             pnum=pft_num(pref)
             ptype=pft_type(pnum)

             !--------------
             !--DEAD POOLS--
             do n=1,npoollu
                pool_init = sib%g(i)%l(l)%equibdt%poollu_init(n)
                pool_end = sib%g(i)%l(l)%pooldt%poollu(n)

                !...calculate time average decay rates

                !...for continuing spinup, this avoids pool_end
                !...being set as pool_equib below?? (doesn't work)
                !if (spinup_continue) then
                !   sib%g(i)%l(l)%equibdt%poollu_totgain(n) = dzero
                !   sib%g(i)%l(l)%equibdt%poollu_totloss(n) = dzero
                !endif

                ave_gain = sib%g(i)%l(l)%equibdt%poollu_totgain(n)
                ave_loss = sib%g(i)%l(l)%equibdt%poollu_totloss(n)

                if (pool_end .gt. dzero) then
                     ave_k_rate = ave_loss / pool_end
                else
                     ave_k_rate = dzero
                endif

                !...save the net gains/losses/end
                sib%g(i)%l(l)%equibdt%poollu_gain(n) = ave_gain
                sib%g(i)%l(l)%equibdt%poollu_loss(n) = ave_loss
                sib%g(i)%l(l)%equibdt%poollu_end(n) = pool_end

                !...solve for ratio of input/output
                if (ave_loss .gt. dzero) then
                    sib%g(i)%l(l)%equibdt%poollu_ratio(n) = ave_gain / ave_loss
                endif
                pdiffr=abs(sib%g(i)%l(l)%equibdt%poollu_ratio(n) - done)

                !...solve for equilibrium pool sizes
                if (ave_k_rate .gt. dzero) then
                     sib%g(i)%l(l)%equibdt%poollu_equib(n) = &
                          ave_gain / ave_k_rate
                elseif ((ave_gain > dzero) .or. (ave_loss > dzero)) then
                     sib%g(i)%l(l)%equibdt%poollu_equib(n) = pool_end
                else
                     sib%g(i)%l(l)%equibdt%poollu_equib(n) = dzero
                     sib%g(i)%l(l)%equibdt%poollu_ratio(n) = done
                endif !zero input/output

                !...calculate starting/ending ratios for spinup constraints
                if (sib%g(i)%l(l)%equibdt%poollu_equib(n) .gt. dzero) then
                    init_ratio = pool_init/sib%g(i)%l(l)%equibdt%poollu_equib(n)
                    end_ratio = pool_end/sib%g(i)%l(l)%equibdt%poollu_equib(n)
                else
                    init_ratio = done
                    end_ratio = done
                endif
                pdiffi = abs(init_ratio - done)
                pdiffe = abs(end_ratio - done)

                !...determine if the dead pools are spunup
                if ((pdiffr <= spinup_threshold) .and. &
                    ((pdiffi <= spinup_threshold) .or. &
                     (pdiffe <= spinup_threshold))) then
                     sib%g(i)%l(l)%equibdt%poollu_notdone(n) = .false.
                else
                     sib%g(i)%l(l)%equibdt%poollu_notdone(n) = .true.
                endif
             enddo  !n=1,npoollu

             !...Determine if the surface pools are spun-up
             sib%g(i)%l(l)%equibdt%deadsfc_init  = &
                    sib%g(i)%l(l)%equibdt%poollu_init(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(strlp)
             sib%g(i)%l(l)%equibdt%deadsfc_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(strlp)
             sib%g(i)%l(l)%equibdt%deadsfc_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(strlp)
             sib%g(i)%l(l)%equibdt%deadsfc_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(cdbp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(metlp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(strlp)

             if (sib%g(i)%l(l)%equibdt%deadsfc_loss > 0.) then
                    sib%g(i)%l(l)%equibdt%deadsfc_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsfc_gain / &
                        sib%g(i)%l(l)%equibdt%deadsfc_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsfc_ratio = 1.0
             endif

             pdiffr = abs(sib%g(i)%l(l)%equibdt%deadsfc_ratio - 1.0)
             if (pdiffr <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsfc_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsfc_notdone = .true.
             endif

             !...same spinup check calculations for C-13 sfc pools
             sib%g(i)%l(l)%equibdt%deadsfc_c13_init  = &
                    sib%g(i)%l(l)%equibdt%poollu_init(cdbpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(metlpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(strlpc13)
             sib%g(i)%l(l)%equibdt%deadsfc_c13_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(cdbpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(metlpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(strlpc13)
             sib%g(i)%l(l)%equibdt%deadsfc_c13_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(cdbpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(metlpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(strlpc13)
             sib%g(i)%l(l)%equibdt%deadsfc_c13_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(cdbpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(metlpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(strlpc13)

             if (sib%g(i)%l(l)%equibdt%deadsfc_c13_loss > 0.) then
                    sib%g(i)%l(l)%equibdt%deadsfc_c13_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsfc_c13_gain / &
                        sib%g(i)%l(l)%equibdt%deadsfc_c13_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsfc_c13_ratio = 1.0
             endif

             pdiffrc13 = abs(sib%g(i)%l(l)%equibdt%deadsfc_c13_ratio - 1.0)
             if (pdiffrc13 <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsfc_c13_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsfc_c13_notdone = .true.
             endif

             !...same spinup check calculations for C-14 sfc pools
             sib%g(i)%l(l)%equibdt%deadsfc_c14_init  = &
                    sib%g(i)%l(l)%equibdt%poollu_init(cdbpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(metlpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(strlpc14)
             sib%g(i)%l(l)%equibdt%deadsfc_c14_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(cdbpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(metlpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(strlpc14)
             sib%g(i)%l(l)%equibdt%deadsfc_c14_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(cdbpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(metlpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(strlpc14)
             sib%g(i)%l(l)%equibdt%deadsfc_c14_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(cdbpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(metlpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(strlpc14)

             if (sib%g(i)%l(l)%equibdt%deadsfc_c14_loss > 0.) then
                    sib%g(i)%l(l)%equibdt%deadsfc_c14_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsfc_c14_gain / &
                        sib%g(i)%l(l)%equibdt%deadsfc_c14_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsfc_c14_ratio = 1.0
             endif

             pdiffrc14 = abs(sib%g(i)%l(l)%equibdt%deadsfc_c14_ratio - 1.0)
             if (pdiffrc14 <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsfc_c14_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsfc_c14_notdone = .true.
             endif


             !...Determine if the soil pools are spunup
             sib%g(i)%l(l)%equibdt%deadsoil_init = &
                    sib%g(i)%l(l)%equibdt%poollu_init(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(armp)
             sib%g(i)%l(l)%equibdt%deadsoil_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(armp)
             sib%g(i)%l(l)%equibdt%deadsoil_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(armp)
             sib%g(i)%l(l)%equibdt%deadsoil_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slitp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slowp) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(armp)

             if (sib%g(i)%l(l)%equibdt%deadsoil_loss > dzero) then
                    sib%g(i)%l(l)%equibdt%deadsoil_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsoil_gain / &
                        sib%g(i)%l(l)%equibdt%deadsoil_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsoil_ratio = 1.0
             endif

             pdiffr = abs(sib%g(i)%l(l)%equibdt%deadsoil_ratio - 1.0)
             if (pdiffr <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsoil_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsoil_notdone = .true.
             endif

             !...same spinup check calculations for C-13 soil pools
             sib%g(i)%l(l)%equibdt%deadsoil_c13_init = &
                    sib%g(i)%l(l)%equibdt%poollu_init(slitpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(slowpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(armpc13)
             sib%g(i)%l(l)%equibdt%deadsoil_c13_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(slitpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(slowpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(armpc13)
             sib%g(i)%l(l)%equibdt%deadsoil_c13_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slitpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slowpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(armpc13)
             sib%g(i)%l(l)%equibdt%deadsoil_c13_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slitpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slowpc13) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(armpc13)
             if (sib%g(i)%l(l)%equibdt%deadsoil_c13_loss > dzero) then
                    sib%g(i)%l(l)%equibdt%deadsoil_c13_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsoil_c13_gain / &
                        sib%g(i)%l(l)%equibdt%deadsoil_c13_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsoil_c13_ratio = 1.0
             endif

             pdiffrc13 = abs(sib%g(i)%l(l)%equibdt%deadsoil_c13_ratio - 1.0)
             if (pdiffrc13 <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsoil_c13_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsoil_c13_notdone = .true.
             endif

             !...same spinup check calculations for C-14 soil pools
             sib%g(i)%l(l)%equibdt%deadsoil_c14_init = &
                    sib%g(i)%l(l)%equibdt%poollu_init(slitpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(slowpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_init(armpc14)
             sib%g(i)%l(l)%equibdt%deadsoil_c14_end = &
                    sib%g(i)%l(l)%equibdt%poollu_end(slitpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(slowpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_end(armpc14)
             sib%g(i)%l(l)%equibdt%deadsoil_c14_gain = &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slitpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(slowpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_gain(armpc14)
             sib%g(i)%l(l)%equibdt%deadsoil_c14_loss = &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slitpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(slowpc14) + &
                    sib%g(i)%l(l)%equibdt%poollu_loss(armpc14)
             if (sib%g(i)%l(l)%equibdt%deadsoil_c14_loss > dzero) then
                    sib%g(i)%l(l)%equibdt%deadsoil_c14_ratio = &
                        sib%g(i)%l(l)%equibdt%deadsoil_c14_gain / &
                        sib%g(i)%l(l)%equibdt%deadsoil_c14_loss
             else
                   sib%g(i)%l(l)%equibdt%deadsoil_c14_ratio = 1.0
             endif

             pdiffrc14 = abs(sib%g(i)%l(l)%equibdt%deadsoil_c14_ratio - 1.0)
             if (pdiffrc14 <= spinup_threshold) then
                 sib%g(i)%l(l)%equibdt%deadsoil_c14_notdone = .false.
             else
                 sib%g(i)%l(l)%equibdt%deadsoil_c14_notdone = .true.
             endif


             !--------------
             !--LIVE POOLS--
             sib%g(i)%l(l)%equiblt%poolpft_end(:) = &
                    sib%g(i)%l(l)%poollt%poolpft(:)

                !...For bare ground do not calculate equilibrium pools,
                !...instead set them to end pools
                if (ptype == type_bare) then
                    sib%g(i)%l(l)%equibdt%poollu_notdone(:) = .false.
                    sib%g(i)%l(l)%equibdt%poollu_equib(:) = &
                          sib%g(i)%l(l)%equibdt%poollu_end(:)

                    sib%g(i)%l(l)%equibdt%deadsfc_notdone = .false.
                    sib%g(i)%l(l)%equibdt%deadsoil_notdone = .false.

                    sib%g(i)%l(l)%equibdt%deadsfc_c13_notdone = .false.
                    sib%g(i)%l(l)%equibdt%deadsoil_c13_notdone = .false.
                    sib%g(i)%l(l)%equibdt%deadsfc_c14_notdone = .false.
                    sib%g(i)%l(l)%equibdt%deadsoil_c14_notdone = .false.

                    sib%g(i)%l(l)%equiblt%poolpft_notdone(:) = .false.
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = &
                           sib%g(i)%l(l)%equiblt%poolpft_end(:)
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = 0.0

                    sib%g(i)%l(l)%equiblt%live_notdone = .false.
                    sib%g(i)%l(l)%equiblt%live_c13_notdone = .false.
                    sib%g(i)%l(l)%equiblt%live_c14_notdone = .false.
                else
                    !...Calculate equilibrium PFT pools
                    do n=1,npoolpft

                        !...calculate time average decay rates
                        ave_gain = sib%g(i)%l(l)%equiblt%poolpft_totgain(n)
                        ave_loss = sib%g(i)%l(l)%equiblt%poolpft_totloss(n)
                        IF (sib%g(i)%l(l)%poollt%poolpft(n) .gt. dzero) THEN
                            ave_k_rate = ave_loss /sib%g(i)%l(l)%poollt%poolpft(n)
                        ELSE
                            ave_k_rate = dzero
                        ENDIF

                        !...save net gains/losses
                        sib%g(i)%l(l)%equiblt%poolpft_gain(n) = ave_gain
                        sib%g(i)%l(l)%equiblt%poolpft_loss(n) = ave_loss
                        if (ave_loss > dzero) then
                           sib%g(i)%l(l)%equiblt%poolpft_ratio(n) = ave_gain/ave_loss
                        else
                           sib%g(i)%l(l)%equiblt%poolpft_ratio(n) = done
                        endif

                        !...solve for equilibrium pool sizes
                        if (ave_k_rate > dzero) then
                            sib%g(i)%l(l)%equiblt%poolpft_equib(n) = &
                                 MAX( poolcon(pnum)%poolpft_min(n), &
                                 ave_gain / ave_k_rate )
                        elseif ((ave_gain > dzero) .or. (ave_loss > dzero)) then
                            sib%g(i)%l(l)%equiblt%poolpft_equib(n) = &
                                MAX( poolcon(pnum)%poolpft_min(n), &
                                sib%g(i)%l(l)%equiblt%poolpft_end(n))
                        else
                           if (((single_pt) .or. (subcount == 1)) .and. &
                                (n .eq. lp)) then
                               print*,'   !!!Error in Equilibrium - Leaf Pool Zero In/Out!!!'
                               print*,'   !!!Setting to Minimum Value!!!'
                               print('(a,2F10.4,a,a)'),'    Input/Output/Pool:', &
                                    ave_gain,ave_k_rate,'  ', pool_name(n)
                                    !stop
                            endif
                            sib%g(i)%l(l)%equiblt%poolpft_equib(n) = poolcon(pnum)%poolpft_min(n)
                            sib%g(i)%l(l)%equiblt%poolpft_ratio(n) = 1.0
                        endif !zero input/output

                         !...determine if the pools are spun-up
                         pdiffr=abs(sib%g(i)%l(l)%equiblt%poolpft_ratio(n) - 1.0)
                         if (pdiffr <= spinup_threshold) then
                             sib%g(i)%l(l)%equiblt%poolpft_notdone(n) = .false.
                         else
                             sib%g(i)%l(l)%equiblt%poolpft_notdone(n) = .true.
                         endif

                    enddo !n=1,npoolpft

                    !...determine if using total live carbon is spun-up
                    sib%g(i)%l(l)%equiblt%live_init = &
                         sib%g(i)%l(l)%equiblt%poolpft_init(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(pp)
                    sib%g(i)%l(l)%equiblt%live_end  = &
                         sib%g(i)%l(l)%equiblt%poolpft_end(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(pp)
                    sib%g(i)%l(l)%equiblt%live_gain = &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(pp)
                    sib%g(i)%l(l)%equiblt%live_loss = &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(lp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(wp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(frp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(crp) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(pp)
                    if (sib%g(i)%l(l)%equiblt%live_loss > dzero) then
                          sib%g(i)%l(l)%equiblt%live_ratio = &
                                 sib%g(i)%l(l)%equiblt%live_gain / &
                                 sib%g(i)%l(l)%equiblt%live_loss
                    else
                          sib%g(i)%l(l)%equiblt%live_ratio = 1.
                    endif

                    !...test for spinup determination
                    pdiffr = abs(sib%g(i)%l(l)%equiblt%live_ratio - 1.0)
                    if (pdiffr < spinup_threshold) then
                          sib%g(i)%l(l)%equiblt%live_notdone = .false.
                    else
                          sib%g(i)%l(l)%equiblt%live_notdone = .true.
                    endif

                    !...same as above but for C-13 pools
                    !...determine if using total live carbon is spun-up
                    sib%g(i)%l(l)%equiblt%live_c13_init = &
                         sib%g(i)%l(l)%equiblt%poolpft_init(lpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(wpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(frpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(crpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(ppc13)
                    sib%g(i)%l(l)%equiblt%live_c13_end  = &
                         sib%g(i)%l(l)%equiblt%poolpft_end(lpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(wpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(frpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(crpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(ppc13)
                    sib%g(i)%l(l)%equiblt%live_c13_gain = &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(lpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(wpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(frpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(crpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(ppc13)
                    sib%g(i)%l(l)%equiblt%live_c13_loss = &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(lpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(wpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(frpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(crpc13) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(ppc13)
                    if (sib%g(i)%l(l)%equiblt%live_c13_loss > dzero) then
                          sib%g(i)%l(l)%equiblt%live_c13_ratio = &
                                 sib%g(i)%l(l)%equiblt%live_c13_gain / &
                                 sib%g(i)%l(l)%equiblt%live_c13_loss
                    else
                          sib%g(i)%l(l)%equiblt%live_c13_ratio = 1.
                    endif

                    !...test for spinup determination
                    pdiffrc13 = abs(sib%g(i)%l(l)%equiblt%live_c13_ratio - 1.0)
                    if (pdiffrc13 < spinup_threshold) then
                          sib%g(i)%l(l)%equiblt%live_c13_notdone = .false.
                    else
                          sib%g(i)%l(l)%equiblt%live_c13_notdone = .true.
                    endif

                    !...same as above but for C-14 pools
                    !...determine if using total live carbon is spun-up
                    sib%g(i)%l(l)%equiblt%live_c14_init = &
                         sib%g(i)%l(l)%equiblt%poolpft_init(lpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(wpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(frpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(crpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_init(ppc14)
                    sib%g(i)%l(l)%equiblt%live_c14_end  = &
                         sib%g(i)%l(l)%equiblt%poolpft_end(lpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(wpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(frpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(crpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_end(ppc14)
                    sib%g(i)%l(l)%equiblt%live_c14_gain = &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(lpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(wpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(frpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(crpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_gain(ppc14)
                    sib%g(i)%l(l)%equiblt%live_c14_loss = &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(lpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(wpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(frpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(crpc14) + &
                         sib%g(i)%l(l)%equiblt%poolpft_loss(ppc14)
                    if (sib%g(i)%l(l)%equiblt%live_c14_loss > dzero) then
                          sib%g(i)%l(l)%equiblt%live_c14_ratio = &
                                 sib%g(i)%l(l)%equiblt%live_c14_gain / &
                                 sib%g(i)%l(l)%equiblt%live_c14_loss
                    else
                          sib%g(i)%l(l)%equiblt%live_c14_ratio = 1.
                    endif

                    !...test for spinup determination
                    pdiffrc14 = abs(sib%g(i)%l(l)%equiblt%live_c14_ratio - 1.0)
                    if (pdiffrc14 < spinup_threshold) then
                          sib%g(i)%l(l)%equiblt%live_c14_notdone = .false.
                    else
                          sib%g(i)%l(l)%equiblt%live_c14_notdone = .true.
                    endif

                endif !pft not bare

                !!!set spinup_done using individual pools!!!
                if  (sib%g(i)%l(l)%equiblt%poolpft_notdone(lp)  .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(frp) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(crp) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(wp)  .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(cdbp) .or.  &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(metlp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(strlp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slitp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slowp) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(armp) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(lpc13)  .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(frpc13) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(crpc13) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(wpc13)  .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(cdbpc13) .or.  &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(metlpc13) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(strlpc13) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slitpc13) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slowpc13) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(armpc13)  .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(lpc14)  .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(frpc14) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(crpc14) .or. &
                     sib%g(i)%l(l)%equiblt%poolpft_notdone(wpc14)  .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(cdbpc14) .or.  &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(metlpc14) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(strlpc14) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slitpc14) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(slowpc14) .or. &
                     sib%g(i)%l(l)%equibdt%poollu_notdone(armpc14)) then
                         sib%g(i)%l(l)%equibdt%lupft_spunup = .false.
                         sib%g(i)%gdiagt%gridcell_spunup = .false.
                         spinup_done = .false.
                    endif

                !!!set spinup_done for PFT pools using live total!!!
                !if (sib%g(i)%l(l)%equibdt%deadsfc_notdone .or. &
                !     sib%g(i)%l(l)%equibdt%deadsoil_notdone .or. &
                !     sib%g(i)%l(l)%equiblt%live_notdone) then
                !        sib%g(i)%l(l)%equibdt%lupft_spunup = .false.
                !        sib%g(i)%gdiagt%gridcell_spunup = .false.
                !        spinup_done = .false.
                !endif

         endif  !.not. lupft_spunup
      enddo !landunit

   endif  !.not. gridcell_spunup
enddo  !subcount

end subroutine equipools_calc_continue


!============================================
subroutine equipools_restart_continue(pnum, &
     poolpft_end, poolpft_equib, poolpft_flay, poolpft_out)
!============================================
! Sets the equilibrium pools to values that can
!   be used to restart a simulation.

use kinds
use module_sibconst, only: &
   npoolpft, nsoil, &
   spinup
use module_pftinfo, only: &
   pft_type, pft_group, &
   type_decid, type_grass, type_crop
use module_poolinfo, only: &
   pool_indx_leaf, pool_indx_stwd, &
   pool_indx_prod, pool_indx_lay, &
   pool_indx_leaf_c13, pool_indx_stwd_c13, &
   pool_indx_prod_c13, &
   pool_indx_leaf_c14, pool_indx_stwd_c14, &
   pool_indx_prod_c14

implicit none

!...input variables
integer(i4), intent(in) :: pnum
real(r8), dimension(npoolpft), intent(in) :: poolpft_end, poolpft_equib
real(r8), dimension(npoolpft,nsoil), intent(in) :: poolpft_flay
real(r8), dimension(npoolpft,nsoil), intent(inout) :: poolpft_out

!...local variables
integer(i4) :: n,s
integer(byte) :: ptype, pgroup
real(r8) :: poolrestart
integer(i4) :: lp,wp,pp
integer(i4) :: lpc13,wpc13,ppc13
integer(i4) :: lpc14,wpc14,ppc14

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
wp =  pool_indx_stwd
pp =  pool_indx_prod

lpc13 =  pool_indx_leaf_c13-6
wpc13 =  pool_indx_stwd_c13-6
ppc13 =  pool_indx_prod_c13-6

lpc14 =  pool_indx_leaf_c14-12
wpc14 =  pool_indx_stwd_c14-12
ppc14 =  pool_indx_prod_c14-12

!---------------------------
!...Set local variables
ptype = pft_type(pnum)
pgroup = pft_group(pnum)
poolpft_out(:,:) = 0.


!...Set output PFT pools
do n=1,npoolpft/3 !1,5
    poolrestart = poolpft_equib(n)

    !if (.not. spinup) then
      !For deciduous PFTs, set restart to end value
      !for leaf and product pools
      if (ptype == type_decid) then
          if ((n == lp) .or. &
              (n == pp)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For grass PFTs, set restart to end value
      ! for leaf, stem, and product pools
      if (ptype == type_grass) then
          if ((n == lp) .or. &
              (n == wp) .or. &
              (n == pp)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For crop PFTs, set restart to end value
      ! for all live pools
      if (ptype == type_crop) then
         poolrestart = MAX(dzero, poolpft_end(n))
      endif
    !endif

    do s=1,pool_indx_lay(n) !for 1,5 ntpool
       poolpft_out(n,s) = poolrestart * &
            poolpft_flay(n,s)
    enddo
enddo

do n=npoolpft/3+1,2*npoolpft/3 !6,10
    poolrestart = poolpft_equib(n)

    !if (.not. spinup) then
      !For deciduous PFTs, set restart to end value
      !for leaf and product pools
      if (ptype == type_decid) then
          if ((n == lpc13) .or. &
              (n == ppc13)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For grass PFTs, set restart to end value
      ! for leaf, stem, and product pools
      if (ptype == type_grass) then
          if ((n == lpc13) .or. &
              (n == wpc13) .or. &
              (n == ppc13)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For crop PFTs, set restart to end value
      ! for all live pools
      if (ptype == type_crop) then
         poolrestart = MAX(dzero, poolpft_end(n))
      endif
    !endif

    do s=1,pool_indx_lay(n+npoolpft/3+1) !for 12,16 ntpool
       poolpft_out(n,s) = poolrestart * &
            poolpft_flay(n,s)
    enddo
enddo

do n=2*npoolpft/3+1,npoolpft !11,15
    poolrestart = poolpft_equib(n)

    !if (.not. spinup) then
      !For deciduous PFTs, set restart to end value
      !for leaf and product pools
      if (ptype == type_decid) then
          if ((n == lpc14) .or. &
              (n == ppc14)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For grass PFTs, set restart to end value
      ! for leaf, stem, and product pools
      if (ptype == type_grass) then
          if ((n == lpc14) .or. &
              (n == wpc14) .or. &
              (n == ppc14)) then
               poolrestart = MAX(dzero, poolpft_end(n))
          endif
      endif

      !For crop PFTs, set restart to end value
      ! for all live pools
      if (ptype == type_crop) then
         poolrestart = MAX(dzero, poolpft_end(n))
      endif
    !endif

    do s=1,pool_indx_lay(n+2*npoolpft/3+2) !for 23,27 ntpool
       poolpft_out(n,s) = poolrestart * &
            poolpft_flay(n,s)
    enddo
enddo


end subroutine equipools_restart_continue
