
!======================SUBROUTINE SIBSLV=================================

subroutine sibslv(cast, fluxt, radt, sscolt)

!========================================================================
!
!     Calculation of time increments in Tc, Tgs, Theta-m and Qm using an
!        implicit backwards method with explicit coefficients.  
!     Similar to equations (10-15), SA-92B. 
!
!     Longwave feedbacks are now really included
!
!======================================================================== 

!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       DEA   CAS VAPOR PRESSURE INCREMENT (hPa or mb)
!       DTA   CAS TEMPERATURE INCREMENT (K)
!       DTC   CANOPY TEMPERATURE INCREMENT (K)
!       DTD   DEEP SOIL TEMPERATURE INCREMENT (K)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use kinds
use module_local
use module_sibconst, only:  &
    nsoil, nsnow, ntot, &
    print_avec, print_stop

use module_sib, only: &
    cas_type, &
    flux_type, &
    rad_type,  &
    sscol_type
use module_time, only: &
    dtisib, sec_tot

IMPLICIT none

!----------------------------------------------------------------------
!...this routine sets up the coupled system of partial 
!...differential equations described in Sato et al., 
!...with the exception that now Ta and ea are prognostic
!...variables, and so we have two more equations, reflected
!...in two more lines and two more columns as compared 
!...to the old sibslv.F (used for no prognistic CAS calculations)

!...this matrix is expandable, size determined by the number of 
!...snow layers (0 to 5). All energy is passed through the snow layer
!...if snow exists. Partial snowcover is not dealt with. We know
!...that this is physically incorrect, it is an issue that we hope
!...to deal with in the future...

!...VARIABLES
!...1 : TREF (lowest model layer/driver data temperature)
!...2 : EREF (lowest model layer/driver data water vapor)
!...3 : TA   (CAS temperature) 
!...4 : EA   (CAS water vapor)
!...5 : TC   (vegetation temperature)
!...6 : TG   (Ground or snow surface temperature)
!...7-14 to 19 : TD (interior soil layers temperature)
!...last member: TD (bottom soil level)
!
!********************************************************************

!...input variables
type(cas_type), intent(in)   :: cast
type(flux_type), intent(in)  :: fluxt
type(rad_type), intent(in)   :: radt
type(sscol_type), intent(in) :: sscolt

!...local variables
integer :: error_flag,i,k

!...matrix arrays - variable, due to potential snow layers
real(r8),dimension(ntot-sscolt%nsl, ntot-sscolt%nsl) :: avec
real(r8),dimension(ntot-sscolt%nsl) :: bvec
integer(i4),dimension(ntot-sscolt%nsl) :: cvec   
!-------------------------------------------------------

!.....zero all of 'em out first, then fill 'em in...
avec(:,:) = dzero

!1     TREF EQUATION       
avec(1,1) =  done                               
bvec(1)   =  dzero

!2     EREF EQUATION  
avec(2,2) =  done
bvec(2)   =  dzero        

!3     TA EQUATION 
if (sscolt%nsl == 0) then
    avec(3,1)  =  hadth                                    
    avec(3,3)  =  cast%hcapcas   * dtisib  & 
                    + hadta  - hcdta - hgdta 
    avec(3,5)  =  - hcdtc     
    avec(3,6)  = - hgdtg
    bvec(3)    =  fluxt%hc * dtisib - fluxt%fss * dtisib &
                    + fluxt%hg * dtisib 
else
    avec(3,1)  =  hadth                                    
    avec(3,3)  =  cast%hcapcas   * dtisib  & 
                    + hadta  - hcdta - hsdta
    avec(3,5)  =  - hcdtc     
    avec(3,6)  = - hsdts
    bvec(3)    =  fluxt%hc * dtisib - fluxt%fss * dtisib &
                    + fluxt%hs * dtisib    
endif

!4    EA EQUATION  
if (sscolt%nsl == 0) then
    avec(4,2)  =  eadem                                  
    avec(4,4)  =  cast%vcapcas   * dtisib  &
                    + eadea  - ecdea - egdea
    avec(4,5)  =  - ecdtc               
    avec(4,6)  =  - egdtg         
    bvec(4)    =  fluxt%ec * dtisib - fluxt%fws * dtisib &
                    + fluxt%eg * dtisib 
else
    avec(4,2)  =  eadem                                  
    avec(4,4)  =  cast%vcapcas   * dtisib  &
                    + eadea  - ecdea - esdea  
    avec(4,5)  =  - ecdtc               
    avec(4,6)  =  - esdts 
    bvec(4)    =  fluxt%ec * dtisib - fluxt%fws * dtisib &
                    + fluxt%es * dtisib   
endif

!5    TC EQUATION 
avec(5,3)  =  hcdta                                    
avec(5,4)  =  ecdea          
avec(5,5)  =  cast%hcapc * dtisib + hcdtc    &
               + ecdtc + lcdtc
avec(5,6)  =  lcdtg
bvec(5)    =  radt%radtc - fluxt%hc * dtisib - fluxt%ec * dtisib


!6    TOP SOIL LAYER (TD1 OR TG) EQUATION
if (sscolt%nsl == 0) then !NO SNOW CASE
     avec(6,3)  =  hgdta                                  
     avec(6,4)  =  egdea          
     avec(6,5)  =  lgdtc               
     avec(6,6)  =  sscolt%slamda(1) &
                    + sscolt%shcap(1) * dtisib    &
                    + hgdtg + egdtg + lgdtg          
     avec(6,7)  =  -sscolt%slamda(1) 
     bvec(6)    =  radt%radtg - fluxt%hg * dtisib &
                    - fluxt%eg * dtisib  &
                    - sscolt%slamda(1) * (sscolt%td(1) &
                    - sscolt%td(2))
else   ! SNOW CASE
     avec(6,3)  =  hsdta                            
     avec(6,4)  =  esdea       
     avec(6,5)  =  lsdtc   
     avec(6,6)  =  sscolt%slamda(sscolt%nsl+1) &
                    + sscolt%shcap(sscolt%nsl+1) * dtisib &
                    + hsdts + esdts + lsdts 
     avec(6,7)  =  -sscolt%slamda(sscolt%nsl+1)
     bvec(6)    =  radt%radts &
                    - fluxt%hs * dtisib &
                    - fluxt%es * dtisib &
                    - sscolt%slamda(sscolt%nsl+1) * (sscolt%td(sscolt%nsl+1) &
                    - sscolt%td(sscolt%nsl+2))
endif

!7-(NSOIL+SSCOLT%NSL)    INTERIOR SOIL LAYERS
do i = 7, (ntot-1) - sscolt%nsl   ! matrix indices
    k = i - nsnow + sscolt%nsl    ! soil layer indices
     avec(i,i-1) = -sscolt%slamda(k-1)
     avec(i,i)   = sscolt%shcap(k)*dtisib &
                    + sscolt%slamda(k) + sscolt%slamda(k-1)
     avec(i,i+1) = -sscolt%slamda(k)
     bvec(i)     = sscolt%slamda(k)*(sscolt%td(k+1) - sscolt%td(k))  &
                    - sscolt%slamda(k-1) * (sscolt%td(k) - sscolt%td(k-1))
enddo

!    BOTTOM SOIL LAYER
i = ntot - sscolt%nsl
avec(i,i-1) =  sscolt%slamda(nsoil - 1) &
                    * (sscolt%td(nsoil) - sscolt%td(nsoil - 1))             
avec(i,i) =  sscolt%shcap(nsoil) * dtisib &
                    + sscolt%slamda(nsoil - 1 ) 
bvec(i) =  - sscolt%slamda(nsoil - 1) &
                    * (sscolt%td(nsoil) - sscolt%td(nsoil - 1))


!     PRINT INFO IF REQUESTED
if (print_avec) then
    print*,'Second: ',sec_tot
    do k=1,8
        print*,' avec,bvec,k: ',avec(k,k),bvec(k),k
    enddo
    if (print_stop) stop
endif

!     SOLVE MATRIX EQUATION   
call dgesv( ntot - sscolt%nsl, 1, avec,      &
            ntot - sscolt%nsl, cvec, bvec, &
            ntot - sscolt%nsl, error_flag )
 
    dth = bvec(1)           
    dqm = bvec(2)
    dta = bvec(3)
    dea = bvec(4)
    dtc = bvec(5)        

    do i=6,ntot-sscolt%nsl
        dtd(i-5+sscolt%nsl) = bvec(i)
    enddo

end subroutine sibslv
