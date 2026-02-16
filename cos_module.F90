!================SUBROUTINE COS=========================
!    
!    Module to set pcosm using either bcosm (500 ppt) or 
!    the TM5 mixing ratio obtained through inversion that
!    covers seasonal and spatial variability. 
! 
!    Choice for either one is based on tm5mr_switch in namel_sibdrv
!    tm5mr_switch = .false.  = 500 ppt
!    tm5mr_switch = .true. = uses tm5 mixing ratio. 
!   
!
!-------------------------------------------------------

subroutine set_cos(gprogt)

    use module_pparams, only: p0_sfc, bcosm
    use module_sib, only: sib, gprog_type
    use module_sibconst, only: tm5mr_switch, subcount
    use module_phosib, only: &
        pressure

    implicit none

!    type(gprog_type), dimension(subcount), intent(inout) :: gprogt
    integer :: n
    type(gprog_type), intent(inout) :: gprogt
!    type (gprog_type) :: gprogt
    pressure = dble(gprogt%ps) * 100.0

    if (tm5mr_switch) then
      do n=1,subcount
         ! pcosm in pa, cosm_tm5 in ppt
         ! gprogt%pcosm = (gprogt%cosm_tm5*p0_sfc)/1.E12
         sib%g(n)%gprogt%pcosm = dble(sib%g(n)%gprogt%cosm_tm5*pressure)/1.E12
      enddo
    else
      do n=1,subcount
         ! pcosm in pa, bcosm in ppt
         ! gprogt%pcosm = (bcosm*p0_sfc)/1.E12
         sib%g(n)%gprogt%pcosm = dble(bcosm*pressure)/1.E12
      enddo
    endif

    !enddo

end subroutine set_cos
