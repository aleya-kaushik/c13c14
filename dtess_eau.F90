!------------------------------------------------------------------------------
subroutine dtess_eau (len, pl, tl, ess, dtess)
!------------------------------------------------------------------------------
!eau_sat computes the saturation mixing ratio, vapor pressure, and saturation,
!and their derivatives with respect to temperature over water, ice and mixed-
!phase clouds. the parameterization is that used in the Community Climate Com-
!munity Climate Model at NCAR.
!Laura D. Fowler /slikrock (07-01-01).

!----------------

use kinds
use module_pparams
      
implicit none

!input arguments:
!-----------------
integer, intent(in) :: len                !length of vector.
real(r8), intent(in), dimension(len):: &
   pl,               &!pressure                                           (Pa).
   tl                 !temperature                                         (K).

!output arguments:
!-----------------
real (r8), intent(out), dimension(len) :: &
   ess,              &!saturation vapor pressure                          (Pa).
   dtess              !derivative of es with respect to temperature     (Pa/K).

!local variables:
!-----------------
integer:: i
real (r8) :: tstl , t0tl
real (r8), dimension(len):: &
      esw ,     dtesw ,      esi ,     dtesi , &
      esm ,     dtesm ,      tl0 ,     wghtm 
real (r8):: &
       e1 ,   e2 ,     f ,    f1 , &
       f2 ,   f3 ,    f4 ,    f5 , &
       lphase , term1 , term2 , term3     

!------------------------------------------------------------------------------
!initialization of different arrays:

tl0    = tl
esw    = 0.0_r8
esi    = 0.0_r8
esm    = 0.0_r8
dtesw  = 0.0_r8
dtesi  = 0.0_r8
dtesm  = 0.0_r8

ess    = 0.0_r8
dtess  = 0.0_r8

!saturation over water:

do i = 1, len

   tl0(i)    = max(twmin,tl0(i))
   tl0(i)    = min(twmax,tl0(i))
   tstl      = tsref / tl0(i)
   e1        = 11.344*(1.0 - tl0(i)/tsref)
   e2        = -3.49149*(tstl - 1.0)
   f1        = -7.90298*(tstl - 1.0)
   f2        = 5.02808*log10(tstl)
   f3        = -1.3816*(10.0**e1-1.0)/10000000.0
   f4        = 8.1328*(10.0**e2-1.0)/1000.0
   f5        = log10(psref)
   f         = f1 + f2 + f3 + f4 + f5

   esw(i)    = (10.0**f)*1.e+02
   esw(i)    = min(esw(i),pl(i)*0.9)
   dtesw(i)  = lvap*esw(i)/(rv*tl0(i)*tl0(i))

   ess(i)    = esw(i)
   dtess(i)  = dtesw(i)


!saturation over ice:

   if(tl0(i)<timax) then

      tl0(i)    = max(tl0(i),timin)
      t0tl      = tice / tl0(i)
      term1     = 2.01889049/(t0tl)
      term2     = 3.56654*log(t0tl)
      term3     = 20.947031*(t0tl)

      esi(i)    = 575.185606e10*exp(-(term1 + term2 + term3))
      esi(i)    = min(esi(i),pl(i)*0.9)
      dtesi(i)  = lsub*esi(i)/(rv*tl0(i)*tl0(i))

      ess(i)    = esi(i)
      dtess(i)  = dtesi(i)

   endif

!interpolated saturation variables:

   if(tl0(i)<tbgmax .and. tl0(i)>=tbgmin) then

      wghtm(i)  = (tl0(i)-tbgmin)/(tbgmax-tbgmin)
      lphase    = lvap*wghtm(i)+lsub*(1.-wghtm(i))
      esm(i)    = wghtm(i)*esw(i) + (1.-wghtm(i))*esi(i)
      esm(i)    = min(esm(i),pl(i)*0.9)
      dtesm(i)  = lphase*esm(i)/(rv*tl0(i)*tl0(i))

      ess(i)    = esm(i)
      dtess(i)  = dtesm(i)

   endif
enddo

end subroutine dtess_eau
