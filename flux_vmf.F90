subroutine flux_vmf(zzwind,zztemp,z0, &
    ros, spdm, sh, thm, sha, tha, &
    ustar, cuni, cu, ct, ventmf)

    use kinds
    use module_oparams, only:  &
        bunstablM, bunstablT, &
        cunstablM, cunstablT, &
        bstabl, cstabl
    use module_pparams, only:  &
        delta, grav, vkrmn

    implicit none

    !----------------------------------------------------------------------

    !...INPUT VARIABLES
    real(r8), intent(in) ::  &
        zzwind, & ! adjusted wind measurement height (m)       
        zztemp    ! adjusted temp measurement height (m)

    real(r8), intent(in) :: z0
    real(r8), intent(in) :: ros, spdm, sh, thm
    real(r8), intent(in) :: sha, tha

    !...OUTPUT VARIABLES
    real(r8), intent(inout) :: ustar
    real(r8), intent(inout) :: cuni, cu, ct, ventmf

    !----------------------------------------------------------------------  

!*****************************************************************************                                                                       
!     VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum 
!     and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
!     and Boville(1993), and by Beljaars and Holtslag (1991).              
!  
!     References:
!       Beljars and Holtslag (1991): Flux parameterization over land surfaces
!              for atmospheric models. J. Appl. Meteo., 30, 327-341.
!       Holtslag and Boville (1993): Local versus nonlocal boundary-layer 
!              diffusion in a global climate model. J. of Climate, 6, 1825-
!              1842.
!       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
!              atmosphere. Boundary-Layer Meteo., 17, 187-202.
!       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
!              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
!              Layer parameterization, ECMWF, 59-79.
!
!     General formulation:
!        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.) 
!     Transfer coefficients for mommentum and heat fluxes are:
!        CU = CUN*Fm, and
!        CT = CTN*Fh
!        where  CUN and CTN are nutral values of momentum and heat transfers,
!           and Fm and Fh are stability functions derived from surface
!           similarity relationships.     
!*****************************************************************************

    !...PATCH WARNING: an unjustified patch has been put in the code, 
    !...whereupon when cuni=1/cun is calculated, the square root is taken.
    !...this is a patch that makes the results better, but it is 
    !...unjustified scientifically.

    real(r8) ::   &
        wgm,     & ! moisture mixing ratio deficit, 
                   !  CAS to reference layer (kg/kg)
        thgm,    & ! temperature difference (theta) CAS-ref level (K)
        thvgm,   & ! combined temperature and moisture difference
                   !  CAS to reference layer (-)
        z1z0u,   & ! ratio of reference height to roughness length
        z1z0urt, & ! square root of z1z0u
        z1z0t,   & ! ratio of reference height to roughness length
        z1z0trt, & ! square root of z1z0t
        !...currently, z1z0u and z1z0t are identical. theoretically, they 
        !...can be changed for different wind/temp measurement heights. 
        !...they both use zzwind right now.
        cun,     & ! momentum transfer coefficient (?)
        ctn,     & ! thermal transfer coefficient (?)
        temv,    & ! part of Richardson No. calculation
        zrib,    & ! part of Richardson No. calculation
        rib,     & ! Richardson Number
        fmomn,   & !
        fheat      !

    real(r8) ::  &
        ribtemp, & !
        dm,      & !
        dh         !

    zrib = zzwind **2 / zztemp                                                    
                                   
    !                                                                       
    ! SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE         
    ! WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.    
    !  
    wgm   = sha - sh                                                                     
    thgm  = tha  - thm     
    thvgm = thgm + tha * delta * wgm       

    !   Ratio of reference height (zwind/ztemp) and roughness length:
    z1z0u = zzwind/z0
    z1z0urt = sqrt( z1z0U )
    z1z0u = log( z1z0U )
    z1z0t = zzwind/z0
    z1z0trt = sqrt( z1z0t )
    z1z0t = log( z1z0t )

    !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:
    cun = vkrmn*vkrmn / (z1z0u*z1z0u )   !neutral Cm & Ct
    ctn = vkrmn*vkrmn / (z1z0t*z1z0t )

    !...PATCH-when 1/cun is calculated, the square root is taken.
    cuni = z1z0u / vkrmn

    !                                                                       
    !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.            
    !   RIB IS THE BULK RICHARDSON NUMBER, between reference
    !   height and surface.

    temv = tha * spdm * spdm   
    temv = max(0.000001_r8,temv)
    rib = -thvgm * grav * zrib / temv 

    !   The stability functions for momentum and heat/moisture fluxes as
    !   derived from the surface-similarity theory by Luis (1079, 1982), and
    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag 
    !   (1991).
    if (rib >= 0.0) then                                           

        !  THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT              

        rib   = min( rib, 0.5_r8)                   
        fmomn = (1. + cstabl * rib * (1.+ bstabl * rib))
        fmomn = 1. / fmomn
        fmomn = max(0.0001_r8,fmomn)
        fheat = fmomn

    else                                  

        !  THE UNSTABLE CASE.    

        ribtemp = abs(rib)
        ribtemp = sqrt( ribtemp )
        dm      = 1. + cunstablM * cun * z1z0Urt * ribtemp
        dh      = 1. + cunstablT * ctn * z1z0Trt * ribtemp
        fmomn   = 1. - (bunstablM * rib ) / dm
        fheat   = 1. - (bunstablT * rib ) / dh

    endif    

    !   surface-air transfer coefficients for momentum CU, for heat and 
    !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

    cu = cun * fmomn 
    ct = ctn * fheat

    !   Ustar and ventlation mass flux: note that the ustar and ventlation 
    !   are calculated differently from the Deardoff's methods due to their
    !   differences in define the CU and CT.

    ustar  = spdm * spdm * cu 
    ustar  = sqrt( ustar ) 
    ventmf = ros * ct * spdm  
    !                                                                       
    ! Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION  
    ! VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
    ! free convection conditions.                                            

end subroutine flux_vmf
