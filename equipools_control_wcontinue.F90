!==========================================================================
subroutine equipools_control()
!==========================================================================
! Controls pool spin-up
!    - Calculates equilibrium pools
!    - Prints pools to a file and/or to screen
!    - Writes an equilibrium restart file
!    - Updates carbon pools

use kinds
use module_io, only: requib_writef
use module_sibconst, only: &
    spinup, spinup_done, spinup_lnum, &
    spinup_continue

implicit none

if (requib_writef) then
   call equipools_calc()
   call equipools_print()
   call restart_write(.true.)
   call equipools_reset()

   IF ((spinup) .and. (.not. spinup_done) .and. &
       (.not. spinup_continue)) THEN
        print*,''
        print'(a,i4,a)','***Spinup Repitition Number ',spinup_lnum,'***'
   ELSEIF ((spinup) .and. (.not. spinup_done) .and. &
           (spinup_continue)) THEN
           print*,''
           print'(a,i4,a)','***Spinup Repitition Number ', &
                   (spinup_lnum),'***'
   ENDIF
endif

end subroutine equipools_control


!==========================================================================
subroutine equipools_calc()
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
    spinup_threshold, spinup_done

implicit none

!...local variables
real(r8) :: ave_gain    !(mol/m2/s) average external inputs per pool
real(r8) :: ave_loss    !(mol/m2/s) average external outputs per pool
real(r8) :: ave_k_rate  !(1/s) average scaled decay rate constant
real(r8) :: pool_init, pool_end, init_ratio, end_ratio
real(r8) :: pdiffr, pdiffi, pdiffe

!...misc variables
integer(byte) :: ptype
integer(i4) :: i,l,n
integer(i4) :: lp,frp,crp,wp,pp
integer(i4) :: cdbp, metlp, strlp, slitp, slowp, armp
integer(i4) :: pref, pnum

!-------Calculatae quasi-equlibrium pools------------
lp =  pool_indx_leaf
frp = pool_indx_froot
crp = pool_indx_croot
wp =  pool_indx_stwd
pp =  pool_indx_prod

cdbp  = pool_indx_cdb-npoolpft
metlp = pool_indx_metl-npoolpft
strlp = pool_indx_strl-npoolpft
slitp = pool_indx_slit-npoolpft
slowp = pool_indx_slow-npoolpft
armp  = pool_indx_arm-npoolpft

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

                    sib%g(i)%l(l)%equiblt%poolpft_notdone(:) = .false.
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = &
                           sib%g(i)%l(l)%equiblt%poolpft_end(:)
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = 0.0
  
                    sib%g(i)%l(l)%equiblt%live_notdone = .false.
                else
                    !...Calculate equilibrium PFT pools
                    do n=1,npoolpft
  
                        !...calculate time average decay rates
                        ave_gain = sib%g(i)%l(l)%equiblt%poolpft_totgain(n)
                        ave_loss = sib%g(i)%l(l)%equiblt%poolpft_totloss(n)
                        IF (sib%g(i)%l(l)%poollt%poolpft(n) .gt. dzero) THEN
                            ave_k_rate = ave_loss / sib%g(i)%l(l)%poollt%poolpft(n)
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
                               print('(a,2F10.4,a,a)'),'    Input/Output/Pool: ', &
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
                         sum(sib%g(i)%l(l)%equiblt%poolpft_init)
                    sib%g(i)%l(l)%equiblt%live_end  = & 
                         sum(sib%g(i)%l(l)%equiblt%poolpft_end)
                    sib%g(i)%l(l)%equiblt%live_gain = &
                         sum(sib%g(i)%l(l)%equiblt%poolpft_gain)
                    sib%g(i)%l(l)%equiblt%live_loss = &
                         sum(sib%g(i)%l(l)%equiblt%poolpft_loss)
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
                     sib%g(i)%l(l)%equibdt%poollu_notdone(armp)) then
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

end subroutine equipools_calc

!============================================
subroutine equipools_restart(pnum, &
     poolpft_end, poolpft_equib, poolpft_flay, poolpft_out)
!============================================
! Sets the equilibrium pools to values that can
!   be used to restart a simulation.

use kinds
use module_sibconst, only: &
   npoolpft, nsoil
use module_pftinfo, only: &
   pft_type, pft_group, &
   type_decid, type_grass, type_crop
use module_poolinfo, only: &
   pool_indx_leaf, pool_indx_stwd, &
   pool_indx_prod, pool_indx_lay

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

!---------------------------
!...Set local variables
ptype = pft_type(pnum)
pgroup = pft_group(pnum)
poolpft_out(:,:) = 0.


!...Set output PFT pools
do n=1,npoolpft
    poolrestart = poolpft_equib(n)

    !For deciduous PFTs, set restart to end value
    !for leaf and product pools
    if (ptype == type_decid) then
        if ((n == pool_indx_leaf) .or. &
            (n == pool_indx_prod)) then
             poolrestart = MAX(dzero, poolpft_end(n))
        endif
    endif

    !For grass PFTs, set restart to end value
    ! for leaf, stem, and product pools
    if (ptype == type_grass) then
        if ((n == pool_indx_leaf) .or. &
            (n == pool_indx_stwd) .or. &
            (n == pool_indx_prod)) then
             poolrestart = MAX(dzero, poolpft_end(n))
        endif
    endif

    !For crop PFTs, set restart to end value
    ! for all live pools
    if (ptype == type_crop) then
       poolrestart = MAX(dzero, poolpft_end(n))
    endif
    
    do s=1,pool_indx_lay(n)
       poolpft_out(n,s) = poolrestart * &
            poolpft_flay(n,s)
    enddo
enddo


end subroutine equipools_restart

!==========================================================================
subroutine equipools_print()
!==========================================================================
! Control routine for printing the equilibrium pools 
!   to standard output and/or to a file

use kinds
use module_io, only: &
    rank, nchunks,   &
    out_path, requib_filename
use module_sibconst, only: &
    subcount,    &
    spinup, spinup_done, &
    spinup_continue, &
    spinup_lnum, spinup_maxiter, &
    spinup_writetxtf
use module_time, only: end_year

implicit none

!...parameters
integer(i4), parameter :: spinup_nprints=10

!...local variables
character(len=3) :: sname
character(len=4) :: curyear
character(len=6) :: curpref

!...Save filename for output
write(curyear, '(i4.4,i2.2)') end_year
if (rank .gt. 1) then
    write(curpref, '(a,i5.5)') 'p', rank
elseif (rank .ne. nchunks) then
    write(curpref, '(a,i5.5)') 'p', rank
else
    write(curpref, '(a)') ''
endif

if ((spinup) .and. (.not. spinup_continue))  then
   write(sname, '(a1,i2.2)') 's', spinup_lnum-1
   requib_filename = trim(out_path)//"sib_requib"//trim(curpref)//trim(sname)//".nc"
elseif (spinup .and. spinup_continue) then
   write(sname, '(a1,i2.2)') 's', spinup_lnum-1
   requib_filename = trim(out_path)//"sib_requib"//trim(curpref)//trim(sname)//".nc"
elseif (.not. spinup) then
   requib_filename = trim(out_path)//"sib_requib"//trim(curyear)//trim(curpref)//".nc"
endif

!...Print out message to standard output or file if requested
if (.not. spinup) spinup_done = .true.
if (spinup) then
    if (spinup_done) then
        print*,''
        print'(a,i4,a)','Spinup Complete After ',spinup_lnum-1, ' iterations.'
        print*,''

        if (subcount .lt. spinup_nprints) then
            call equipools_prints()
        endif
     elseif (spinup_lnum > spinup_maxiter) then
        print*,''
        print'(a,i4,a)','Spinup Force Stopped After ',spinup_lnum-1,' iterations.'
        print*,''
        spinup_done = .true.

        if (subcount .lt. spinup_nprints) then
            call equipools_prints()
         endif

        if (spinup_writetxtf) call equipools_printf()
     else
        if (subcount .lt. spinup_nprints) then
            call equipools_prints()
        endif
     endif
else
     if (subcount .lt. spinup_nprints) then
          call equipools_prints()
     endif
endif !spinup runs

end subroutine equipools_print


!============================================================
!============================================================
subroutine equipools_printf()
!============================================================

!------------------------
!Prints the equilibrium carbon pools to a file.
!------------------------

use kinds
use module_io, only: out_path, &
    rank, nchunks
use module_pparams, only: &
    mol_to_mg
use module_pftinfo, only: &
    pft_name, pft_num
use module_poolinfo, only: pool_name
use module_sib, only: sib
use module_sibconst, only: &
    subcount, &
    subset, sublatsib, sublonsib,  &
    npoolpft, npoollu

implicit none

! Parameters
!...print pools in moles C? (vs Mg C)
logical, parameter :: eq_pfmol=.false.

! Misc Variables
integer :: i,l,n
integer :: pref, pnum
real(r8) :: pool_convert
character*120 :: filename

! Check for printed units
IF (eq_pfmol) THEN
   pool_convert = 1.0
ELSE
   pool_convert = mol_to_mg
ENDIF

! Print Results
if (rank .ne. nchunks) then
   write(filename,'(a,a,a,i5.5,a)') trim(out_path), 'spinup_info','p',rank,'.txt'
else
   write(filename,'(a,a)') trim(out_path), 'spinup_info.txt'
endif
open(unit=12,file=filename,form='formatted')

write(12,*) ''
write(12,'(a)') '------------------------------------------------'
write(12,'(a)') 'Equilibrium Pool Estimates For Non-Spunup Grid Cells'

do i=1, subcount
  if (.not. sib%g(i)%gdiagt%gridcell_spunup) then
     write(12,'(a)') ''
     write(12,'(a)') ''
     write(12,'(a,i8,2F10.3,a)') '****SiB Point/Lat/Lon=', &
           subset(i),sublatsib(i),sublonsib(i),'****'

     do l=1, sib%g(i)%g_nlu
        if (.not. sib%g(i)%l(l)%equibdt%lupft_spunup) then
           pref=sib%g(i)%l(l)%ipft
           pnum=pft_num(pref)

           write(12,'(a)') ''
           write(12,'(a,i2,a3,a)') '  PFT: ',pref,'   ',pft_name(pnum)
           write(12,'(2a)') '       Pool     Init      End       Equib   Set', &
                            '   Min       Max      Input     Output     In/Out'

           !...print pft test messages
          do n=1,npoolpft
             write(12,'(i4,a8,3f10.4,a,l,a,2f10.4,2f12.4,1f7.4)') &
                 n, trim(pool_name(n)), &
                 (sib%g(i)%l(l)%equiblt%poolpft_init(n)*pool_convert), &
                 (sib%g(i)%l(l)%equiblt%poolpft_end(n)*pool_convert), &
                 (sib%g(i)%l(l)%equiblt%poolpft_equib(n)*pool_convert), &
                 '  ',(.not. sib%g(i)%l(l)%equiblt%poolpft_notdone(n)), &
                 '  ', &
                 sib%g(i)%l(l)%equiblt%poolpft_min(n)*pool_convert,   &
                 sib%g(i)%l(l)%equiblt%poolpft_max(n)*pool_convert,   &
                 sib%g(i)%l(l)%equiblt%poolpft_gain(n)*pool_convert,  &
                 sib%g(i)%l(l)%equiblt%poolpft_loss(n)*pool_convert,  &
                 sib%g(i)%l(l)%equiblt%poolpft_ratio(n)
          enddo
          write(12,'(a,2f10.4,a,l,a,2f12.4,1f7.4)') '-->Live:    ', &
                 sib%g(i)%l(l)%equiblt%live_init*pool_convert, &
                 sib%g(i)%l(l)%equiblt%live_end*pool_convert,   &
                 '               ', & 
                 (.not. sib%g(i)%l(l)%equiblt%live_notdone), &
                 '   ', &
                 sib%g(i)%l(l)%equiblt%live_gain*pool_convert, &
                 sib%g(i)%l(l)%equiblt%live_loss*pool_convert, &
                 sib%g(i)%l(l)%equiblt%live_ratio

          !...print land unit test messages
          write(12,'(a)') ''
          do n=1+npoolpft,npoolpft+npoollu
              write(12,'(i4,a8,3f10.4,a,l,a,2f10.4,2f12.4,1f7.4)') &
                  n, trim(pool_name(n)), &
                  sib%g(i)%l(l)%equibdt%poollu_init(n-npoolpft)*pool_convert, &
                  sib%g(i)%l(l)%equibdt%poollu_end(n-npoolpft)*pool_convert, &
                  sib%g(i)%l(l)%equibdt%poollu_equib(n-npoolpft)*pool_convert, &
                  ' ', &
                  (.not. sib%g(i)%l(l)%equibdt%poollu_notdone(n-npoolpft)), &
                  ' ', &
                  sib%g(i)%l(l)%equibdt%poollu_min(n-npoolpft)*pool_convert, &
                  sib%g(i)%l(l)%equibdt%poollu_max(n-npoolpft)*pool_convert, &
                  sib%g(i)%l(l)%equibdt%poollu_gain(n-npoolpft)*pool_convert, &
                  sib%g(i)%l(l)%equibdt%poollu_loss(n-npoolpft)*pool_convert, &
                  sib%g(i)%l(l)%equibdt%poollu_ratio(n-npoolpft)
           enddo

          write(12,'(a,2f10.4,a,l,a,2f12.4,1f7.4)') '--> Sfc CWD:', &
                 sib%g(i)%l(l)%equibdt%deadsfc_init*pool_convert, &
                 sib%g(i)%l(l)%equibdt%deadsfc_end*pool_convert,   &
                 '               ',  &
                 (.not. sib%g(i)%l(l)%equibdt%deadsfc_notdone), &
                 '   ', &
                 sib%g(i)%l(l)%equibdt%deadsfc_gain*pool_convert, &
                 sib%g(i)%l(l)%equibdt%deadsfc_loss*pool_convert, &
                 sib%g(i)%l(l)%equibdt%deadsfc_ratio
          write(12,'(a,2f10.4,a,l,a,2f12.4,1f7.4)') '-->Soil Tot:', &
                 sib%g(i)%l(l)%equibdt%deadsoil_init*pool_convert, &
                 sib%g(i)%l(l)%equibdt%deadsoil_end*pool_convert,   &
                 '               ', &
                 (.not. sib%g(i)%l(l)%equibdt%deadsoil_notdone), &
                 '   ', &
                 sib%g(i)%l(l)%equibdt%deadsoil_gain*pool_convert, &
                 sib%g(i)%l(l)%equibdt%deadsoil_loss*pool_convert, &
                 sib%g(i)%l(l)%equibdt%deadsoil_ratio
      endif !LU not spunup
   enddo !l=1,g_nlu
 endif  !Grid Cell not spunup
enddo !i=1,nsib

write(12,'(a)') '------------------------------------------------'
write(12,'(a)') ''

end subroutine equipools_printf



!============================================================
!============================================================
subroutine equipools_prints()
!============================================================

!------------------------
!Prints the equilibrium carbon pools to the 
!  screen (standard output).
!------------------------

use kinds
use module_pftinfo, only: &
    pft_name, pft_num, pft_pmeth, pmeth_stg
use module_poolinfo
use module_pparams, only: &
    mol_to_mg
use module_sib, only: sib
use module_sibconst, only: & 
    subcount, &
    subset, sublatsib, sublonsib,  &
    npoolpft, npoollu

implicit none

! Parameters
!...print pools in moles C? (vs Mg C)
logical, parameter :: eq_psmol=.false.

! Misc Variables
integer(i4) :: i,l,n,nref
integer(i4) :: lp,wp,pp
integer(i4) :: pref, pnum
real(r8) :: pool_convert
character(len=12) :: pool_outputt

!==============================
! Set local variables
lp=pool_indx_leaf
wp=pool_indx_stwd
pp=pool_indx_prod

! Check for printed units
IF (eq_psmol) THEN
   pool_outputt='mol C/m2'
   pool_convert = 1.0
ELSE
   pool_outputt='Mg C/ha'
   pool_convert = mol_to_mg
ENDIF

! Print Results
print*,''
print('(a)'),'------------------------------------------------'
print('(a,a)'),'Equilibrium Pool Estimates in ',trim(pool_outputt)

do i=1, subcount
  print*,''
  print*,''
  print('(a,i8,2F10.3,a)'),'****SiB Point/Lat/Lon=', &
        subset(i),sublatsib(i),sublonsib(i),'****'

  do l=1, sib%g(i)%g_nlu
     pref=sib%g(i)%l(l)%ipft
     pnum=pft_num(pref)

     print*,''
     print'(a,i2,a3,a)','  PFT: ',pref,'   ',pft_name(pnum)
     print*,'       Pool     Init      End      Equib   Set  ', &
            '   Min       Max         Input       Output   In/Out'

     !...Print live pool test messages
     do n=1,npoolpft
        print'(i4,a8,3f10.4,a,l,a,2f10.4,2f12.4,1f9.4)', &
            n, trim(pool_name(n)), &
            (sib%g(i)%l(l)%equiblt%poolpft_init(n)*pool_convert), &
            (sib%g(i)%l(l)%equiblt%poolpft_end(n)*pool_convert), &
            (sib%g(i)%l(l)%equiblt%poolpft_equib(n)*pool_convert), &
            '   ',(.not. sib%g(i)%l(l)%equiblt%poolpft_notdone(n)),'  ', &
            sib%g(i)%l(l)%equiblt%poolpft_min(n)*pool_convert, &
            sib%g(i)%l(l)%equiblt%poolpft_max(n)*pool_convert, &
            sib%g(i)%l(l)%equiblt%poolpft_gain(n)*pool_convert, &
            sib%g(i)%l(l)%equiblt%poolpft_loss(n)*pool_convert,  &
            sib%g(i)%l(l)%equiblt%poolpft_ratio(n)
     enddo
     print'(a,2f10.4,a,l,a,2f12.4,1f9.4)','-->Live:    ', &
            sib%g(i)%l(l)%equiblt%live_init*pool_convert, &
            sib%g(i)%l(l)%equiblt%live_end*pool_convert,   &
            '             ', &
            (.not. sib%g(i)%l(l)%equiblt%live_notdone), &
            '                      ', &
            sib%g(i)%l(l)%equiblt%live_gain*pool_convert, &
            sib%g(i)%l(l)%equiblt%live_loss*pool_convert, &
            sib%g(i)%l(l)%equiblt%live_ratio

      !...Print dead pool test messages
      print*,''
      do n=1+npoolpft,npoolpft+npoollu
          nref=n-npoolpft
          print'(i4,a8,3f10.4,a,l,a,2f10.4,2f12.4,1f9.4)', n, trim(pool_name(n)), &
              sib%g(i)%l(l)%equibdt%poollu_init(nref)*pool_convert, &
              sib%g(i)%l(l)%equibdt%poollu_end(nref)*pool_convert, &
              sib%g(i)%l(l)%equibdt%poollu_equib(nref)*pool_convert, &
              '   ',(.not. sib%g(i)%l(l)%equibdt%poollu_notdone(nref)),'  ', &
               sib%g(i)%l(l)%equibdt%poollu_min(nref)*pool_convert, &
               sib%g(i)%l(l)%equibdt%poollu_max(nref)*pool_convert, &
               sib%g(i)%l(l)%equibdt%poollu_gain(nref)*pool_convert, &
               sib%g(i)%l(l)%equibdt%poollu_loss(nref)*pool_convert,  &
               sib%g(i)%l(l)%equibdt%poollu_ratio(nref)
       enddo

      print'(a,2f10.4,a,l,a,2f12.4,1f9.4)','--> Sfc Tot:', &
             sib%g(i)%l(l)%equibdt%deadsfc_init*pool_convert, &
             sib%g(i)%l(l)%equibdt%deadsfc_end*pool_convert,   &
             '             ', &
             (.not. sib%g(i)%l(l)%equibdt%deadsfc_notdone), &
             '                      ', &
             sib%g(i)%l(l)%equibdt%deadsfc_gain*pool_convert, &
             sib%g(i)%l(l)%equibdt%deadsfc_loss*pool_convert, &
             sib%g(i)%l(l)%equibdt%deadsfc_ratio
      print'(a,2f10.4,a,l,a,2f12.4,1f9.4)','-->Soil Tot:', &
             sib%g(i)%l(l)%equibdt%deadsoil_init*pool_convert, &
             sib%g(i)%l(l)%equibdt%deadsoil_end*pool_convert,   &
             '             ', &
             (.not. sib%g(i)%l(l)%equibdt%deadsoil_notdone), &
             '                      ', &
             sib%g(i)%l(l)%equibdt%deadsoil_gain*pool_convert, &
             sib%g(i)%l(l)%equibdt%deadsoil_loss*pool_convert, &
             sib%g(i)%l(l)%equibdt%deadsoil_ratio

      !...Print climatology information
      IF (pft_pmeth(pnum) .eq. pmeth_stg) THEN
         print*,''
         print*,' Climatology Water Availability: '
         print'(a,f8.4)','    Mean CUPR (mm/day): ', &
             sib%g(i)%gprogt%clim_cupr
         print'(a,f8.4)','    Mean Precip (mm/day): ', &
             sib%g(i)%gprogt%clim_precip
         print'(a,f8.4)','    Mean PAW Fraction (-): ', &
            sib%g(i)%l(l)%hydrovt%clim_pawfrw
         print'(a,f8.4)','    Mean TAW Fraction (-): ', &
            sib%g(i)%l(l)%hydrovt%clim_tawfrw
         print'(a,f8.4)','    ClimP: ', &
            sib%g(i)%l(l)%phent%phenc_climp
      ENDIF
      
  enddo !l=1,g_nlu
enddo !i=1,subcount

print*,'------------------------------------------------'
print*,''

end subroutine equipools_prints


!==========================================================================
subroutine equipools_reset()
!==========================================================================
! Resets the pools for spinup runs.

use kinds
use module_sibconst, only: &
    subcount, spinup_done, &
    npoolpft, npoollu
use module_io, only: requib_writef
use module_pftinfo, only: pft_num
use module_poolinfo, only: pool_indx_lay
use module_sib, only: sib

implicit none

!...local variables
integer :: i,l,n,s
integer(i4) :: pref, pnum

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
                   
                call equipools_restart(pnum, &
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
               do n=1,npoollu
                  sib%g(i)%l(l)%equibdt%poollu_init(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%pooldt%poollu(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_min(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)
                  sib%g(i)%l(l)%equibdt%poollu_max(n) = &
                      sib%g(i)%l(l)%equibdt%poollu_equib(n)

                  do s=1,pool_indx_lay(n+npoolpft)
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

end subroutine equipools_reset

