subroutine restart_read()

    !============================================================
    !  This subroutine reads in the initial conditions file and
    !    stores the variables into the sib structure.
    !=============================================================

    use netcdf 
    use kinds
    use module_io
    use module_pparams, only:  &
       bco2m, bcosm, p0_sfc, secs_per_day
    use module_sibconst, only: &
       eqclear_switch,   &
       nsib, npft, nlu, &
       npoolpft, npoollu, npoolcan, &
       nsoil, nsnow, ntot, &
       subcount, subset, &
       spinup, spinup_default, &
       npoolcanc13, spinup_continue
    use module_sib, only: &
       sib, gprog_type
    use module_phosib, only: &
        pressure

    implicit none

    !...netcdf variables
    integer(i4) :: ncid, status
    integer(i4) :: dimid, varid
    integer(i4) :: numdims
 
    !...local variables
    integer(i4) :: nsibf, ntest
    logical :: nsibr_flag
    logical :: file_exist, var_exist
    
    integer(i4) :: i,l,p,v
    integer(i4) :: gref

    integer(byte) :: btvar
    integer(i4), dimension(:,:), allocatable :: ivar2d

    !real(r4), dimension(:), allocatable :: tvar
    !real(r4), dimension(:,:), allocatable :: tvar2d
    !real(r4), dimension(:,:,:), allocatable :: tvar3d
    !real(r4), dimension(:,:,:,:), allocatable :: tvar4d

    real(r8), dimension(:), allocatable :: dvar
    real(r8), dimension(:,:), allocatable :: dvar2d
    real(r8), dimension(:,:,:), allocatable :: dvar3d
    real(r8), dimension(:,:,:,:), allocatable :: dvar4d

!    type(gprog_type), intent(in) :: gprogt
!    pressure = gprogt%ps * 100.0    
!--------------------------------------------
print*,'Setting SiB4 Restart Variables'

! Enquire if initial conditions file exists
file_exist = .true.
inquire(file=trim(ic_file),exist=file_exist)
if (spinup .and. (spinup_default) .and. &
    (.not. (spinup_continue))) file_exist=.false.
if (spinup .and. (spinup_default) .and. &
    (spinup_continue)) file_exist=.true.

! Open and read initial conditions
if (.not. file_exist) then
   print*,'   Using Default Initial Conditions'
   ncid = -1
else
   print*,'   Reading IC File: '
   print*,'     ',trim(ic_file)
   status = nf90_open( trim(ic_file), nf90_nowrite, ncid )
   if (status /= nf90_noerr) then
      print*,'Error Opening IC File.'
      stop
   endif
   
   ! Read and check dimensions
   status = nf90_inq_dimid( ncid, 'nsib', dimid )
   if (status == nf90_noerr) then 
       status = nf90_inquire_dimension( ncid, dimid, len=nsibf )
   endif
   !....currently require restart file to have nsib or subcount elements...
   if (nsibf .eq. nsib) then
      nsibr_flag=.true.
   elseif (nsibf .eq. subcount) then
      nsibr_flag=.false.
   else
       print*, 'NSIB INCORRECT IN RESTART FILE (IC_PATH)'
       print*,' Sim nsib=',nsib, ' File nsib=',ntest, 'Subcount=',subcount
       stop
   endif

   status = nf90_inq_dimid( ncid, 'npft', dimid )
   if (status == nf90_noerr) then
      status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (npft /= ntest) then
       print*, 'NPFT IN RESTART FILE (IC_PATH) INCORRECT'
       print*,'  Sim npft=',npft, ' File npft=',ntest
       !stop
   endif

   status = nf90_inq_dimid( ncid, 'nlu', dimid )
   if (status == nf90_noerr) then
       status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (nlu /= ntest) then
       print*, 'NLU IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim nlu=',nlu, ' File nlu=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'npoolpft', dimid )
   if (status == nf90_noerr) then
       status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (npoolpft /= ntest) then
       print*, 'NPOOLPFT IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim npoolpft=',npoolpft, ' File npoolpft=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'npoollu', dimid )
   if (status == nf90_noerr) then
      status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (npoollu /= ntest) then
       print*, 'NPOOLLU IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim npoollu=',npoollu, ' File npoollu=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'nsoil', dimid )
   if (status == nf90_noerr) then
      status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (nsoil /= ntest) then
       print*, 'NSOIL IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim nsoil=',nsoil, ' File nsoil=',ntest
       stop
   endif

   status = nf90_inq_dimid( ncid, 'nsnow', dimid )
   if (status == nf90_noerr) then
       status = nf90_inquire_dimension( ncid, dimid, len=ntest )
   endif
   if (nsnow /= ntest) then
       print*, 'NSNOW IN RESTART FILE (IC_PATH) INCORRECT'
       print*,' Sim nsnow=',nsnow, ' File nsnow=',ntest
       stop
   endif
endif !ic file exists test

!---------------------------------------------------------------
! Loop over SiB4 variables, setting where appropriate
do v = 1, sibr_nvar
   status = nf90_inq_varid(ncid, trim(sibr_vname(v)), varid)
   var_exist = (status .eq. nf90_noerr)

   select case (sibr_vref(v))

   case (1)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%tmd = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%tmd =  sib%g(i)%gprogt%tm
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%tmd = sibr_vd(v)
              enddo
           endif
       endif

       
   case (2)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%seas_precip = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%seas_precip = MAX(0.01, &
                    (sib%g(i)%gprogt%lspr + sib%g(i)%gprogt%cupr)*secs_per_day)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%seas_precip = MAX(0.01, sibr_vd(v))
              enddo
           endif
       endif

   case (3)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%seas_tm = dvar(gref)
           enddo
           deallocate(dvar)
       else
           do i=1, subcount
              if (sibr_vd(v) .lt. rzero) then
                  sib%g(i)%gprogt%seas_tm = sibr_vd(v)
              else
                  sib%g(i)%gprogt%seas_tm = sib%g(i)%gprogt%tm
              endif
           enddo
       endif
       
   case (4)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              print*,'  npts/varid: ',nsibf,varid
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%clim_cupr = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%clim_cupr = MAX(0.01, &
                    sib%g(i)%gprogt%cupr*secs_per_day)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%clim_cupr = MAX(0.01, sibr_vd(v))
              enddo
           endif
       endif

   case (5)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%clim_precip = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%clim_precip = MAX(0.01, &
                    (sib%g(i)%gprogt%lspr + sib%g(i)%gprogt%cupr)*secs_per_day)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%clim_precip = MAX(0.01, sibr_vd(v))
              enddo
           endif
       endif
       
   case (6)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%clim_tm = dvar(gref)
           enddo
           deallocate(dvar)
       else
           do i=1, subcount
              if (sibr_vd(v) .lt. rzero) then
                  sib%g(i)%gprogt%clim_tm = sibr_vd(v)
              else
                  sib%g(i)%gprogt%clim_tm = sib%g(i)%gprogt%tm
              endif
           enddo
       endif

   case (10)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%eacas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%eacas = sibr_vd(v)
              enddo
           enddo
       endif

   case (11)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%shcas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%shcas = sib%g(i)%gprogt%sh
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%eacas = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (12)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tcas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%tcas = sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%tcas = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (13)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tkecas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%cast%tkecas = sibr_vd(v)
              enddo
           enddo
       endif

   case (14)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tc = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%tc = sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%tc = sibr_vd(v)
                 enddo
              enddo
           endif
       endif
       
   case (15)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cast%tcmin = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%cast%tcmin = sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%cast%tcmin = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (20)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2m = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%co2t%co2m = dble(bco2m)
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%co2t%co2m = sibr_vd(v)
                 enddo
              enddo
           endif
       endif
       
   !case (20)
   !    if (var_exist) then
   !        allocate(dvar(nsibf))
   !        status = nf90_get_var(ncid, varid, dvar)
   !        if (status .ne. nf90_noerr) then
   !           print*,'Error reading restart variable: ', &
   !               trim(sibr_vname(v))
   !           stop
   !        endif

   !        do i=1, subcount
   !           if (nsibr_flag) then
   !               gref = subset(i)
   !           else
   !               gref = i
   !           endif

   !           sib%g(i)%gprogt%co2m = dvar(gref)
   !        enddo
   !        deallocate(dvar)
   !    else
   !        if (sibr_vd(v) .lt. rzero) then
   !           do i=1,subcount
   !              !sib%g(i)%gprogt%pco2m = &
   !              !    (bco2m * p0_sfc/1.E6)
   !              !sib%g(i)%gprogt%pco2m = &
   !              !    dble(bco2m * pressure/1.E6)
   !              sib%g(i)%gprogt%co2m = dble(bco2m) 
   !              !print*,'sibr_vd(v) :',sibr_vd(v)
   !              !print*,'bco2m :',bco2m
   !              !print*,'p0_sfc :',p0_sfc
   !              !print*,'pressure :',pressure
   !              !print*,'pco2m rr 1 :',sib%g(i)%gprogt%pco2m
   !           enddo
   !        else
   !           do i=1, subcount
   !              sib%g(i)%gprogt%co2m = sibr_vd(v)
   !              !print*,'sibr_vd(v) :',sibr_vd(v)
   !              !print*,'pco2m rr 2 :',sib%g(i)%gprogt%pco2m
   !           enddo
   !        endif
   !    endif
       
   case (21)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2cas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2cas = sibr_vd(v)
              enddo
            enddo
       endif

   case (22)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2c = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2c = sibr_vd(v)
              enddo
            enddo
         endif

    case (23)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2i = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2i = sibr_vd(v)
              enddo
            enddo
         endif

   case (24)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%pco2s = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%pco2s = sibr_vd(v)
              enddo
            enddo
       endif
       
   case (25)
       if (var_exist) then
           allocate(dvar(nsibf))
           status = nf90_get_var(ncid, varid, dvar)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              sib%g(i)%gprogt%pcosm = dvar(gref)
           enddo
           deallocate(dvar)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 sib%g(i)%gprogt%pcosm = &
                      dble(bcosm * p0_sfc/1.E12)
              enddo
           else
              do i=1, subcount
                 sib%g(i)%gprogt%pcosm = sibr_vd(v)
              enddo
           endif
       endif
       
   case (26)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cost%coscasp = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%cost%coscasp = sibr_vd(v)
              enddo
            enddo
       endif

   case (27)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%cost%coscas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%cost%coscas = sibr_vd(v)
              enddo
            enddo
       endif

   case (28)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2cas = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2cas = sibr_vd(v)
              enddo
            enddo
       endif
       
   case (29)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif  
              
           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif
                  
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2c = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2c = sibr_vd(v)
              enddo
            enddo
       endif

   case (30)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif  
              
           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif
                  
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2i = dvar2d(gref,l)
              enddo
           enddo  
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2i = sibr_vd(v)
              enddo
            enddo
       endif

   case (31)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif  
              
           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif
                  
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%co2s = dvar2d(gref,l)
              enddo
           enddo  
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%co2s = sibr_vd(v)
              enddo
            enddo
       endif

   case (35)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fluxt%ra = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fluxt%ra = sibr_vd(v)
              enddo
           enddo
       endif

   case (36)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fluxt%rb = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fluxt%rb = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (37)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%rst = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%rst = sibr_vd(v)
              enddo
           enddo
       endif

   case (40)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%assim = sibr_vd(v)
              enddo
           enddo
       endif

       
   case (41)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%assimd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%assimd = sibr_vd(v)
              enddo
           enddo
       endif

   case (42)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_nveg = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_nveg = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (43)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_leaf = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_leaf = sibr_vd(v)
              enddo
           enddo
       endif

   case (44)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_root = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_root = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (45)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grow = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grow = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (46)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_mntn = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_mntn = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (47)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_auto = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_auto = sibr_vd(v)
              enddo
           enddo
       endif

   case (48)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_soil = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_soil = sibr_vd(v)
              enddo
           enddo
       endif

   case (49)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_het = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_het = sibr_vd(v)
              enddo
           enddo
       endif
      
   case (50)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grz = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grz = sibr_vd(v)
              enddo
           enddo
       endif
 
   case (60)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrost%capacc_liq = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrost%capacc_liq = sibr_vd(v)
              enddo
           enddo
       endif

   case (61)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrost%capacc_snow = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrost%capacc_snow = sibr_vd(v)
              enddo
           enddo
       endif

   case (62)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrost%capacg = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrost%capacg = sibr_vd(v)
              enddo
           enddo
       endif

   case (63)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrovt%clim_pawfrw = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrovt%clim_pawfrw = sibr_vd(v)
              enddo
           enddo
       endif

   case (64)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%hydrovt%clim_tawfrw = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%hydrovt%clim_tawfrw = sibr_vd(v)
              enddo
           enddo
       endif

   case (70)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  btvar=bnegone*abs(int(ivar2d(gref,l),kind=byte))
                  sib%g(i)%l(l)%sscolt%nsl = btvar
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%nsl = bnegone*abs(int(sibr_vd(v),kind=byte))
              enddo
           enddo
       endif

   case (71)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,nsnow))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%dz(-nsnow+1:0) = dvar3d(gref,l,1:nsnow)
              enddo
           enddo
           deallocate(dvar3d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%dz(-nsnow+1:0) = sibr_vd(v)
              enddo
           enddo
       endif

   case (72)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,ntot))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%td(-nsnow+1:nsoil) = dvar3d(gref,l,1:ntot)
              enddo
           enddo
           deallocate(dvar3d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    btvar = sib%g(i)%l(l)%sscolt%nsl
                    sib%g(i)%l(l)%sscolt%td(btvar+1:nsoil) = &
                         sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    btvar = sib%g(i)%l(l)%sscolt%nsl
                    sib%g(i)%l(l)%sscolt%td(btvar+1:nsoil) = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (73)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,ntot))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%www_liq(-nsnow+1:nsoil) = dvar3d(gref,l,1:ntot)
              enddo
           enddo
           deallocate(dvar3d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%www_liq(-nsnow+1:nsoil) = sibr_vd(v)
              enddo
           enddo
       endif

   case (74)
       if (var_exist) then
           allocate(dvar3d(nsibf,nlu,ntot))
           status = nf90_get_var(ncid, varid, dvar3d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%sscolt%www_ice(-nsnow+1:nsoil) = dvar3d(gref,l,1:ntot)
              enddo
           enddo
           deallocate(dvar3d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%sscolt%www_ice(-nsnow+1:nsoil) = sibr_vd(v)
              enddo
           enddo
       endif

   case (80)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%co2t%clim_assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%co2t%clim_assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (81)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%vegt%clim_lai = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%vegt%clim_lai = sibr_vd(v)
              enddo
           enddo
       endif

   case (82)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (83)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_assimsm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_assimsm = sibr_vd(v)
              enddo
           enddo
       endif

   case (84)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_pr = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_pr = sibr_vd(v)
              enddo
           enddo
       endif

   case (85)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_prsm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_prsm = sibr_vd(v)
              enddo
           enddo
       endif

   case (86)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_prsdoy = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_prsdoy = sibr_vd(v)
              enddo
           enddo
       endif

   case (87)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_prcdoy = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_prcdoy = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (88)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_tawftop = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_tawftop = sibr_vd(v)
              enddo
           enddo
       endif

   case (89)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_tm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           if (sibr_vd(v) .lt. rzero) then
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%phent%phenave_tm = &
                        sib%g(i)%gprogt%tm
                 enddo
              enddo
           else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%phent%phenave_tm = sibr_vd(v)
                 enddo
              enddo
           endif
       endif

   case (90)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenc_climp = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenc_climp = sibr_vd(v)
              enddo
           enddo
       endif

   case (91)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenc_laimax = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenc_laimax = sibr_vd(v)
              enddo
           enddo
       endif

   case (92)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenc_laimin = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenc_laimin = sibr_vd(v)
              enddo
           enddo
       endif

   case (93)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_env = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_env = sibr_vd(v)
              enddo
           enddo
       endif

   case (94)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_wa = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_wa = sibr_vd(v)
              enddo
           enddo
       endif

   case (95)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phenave_wacsm = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phenave_wacsm = sibr_vd(v)
              enddo
           enddo
       endif

   case (96)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phen_pi = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phen_pi = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (97)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%phen_istage = ivar2d(gref,l)
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%phen_istage = int(sibr_vd(v))
              enddo
           enddo
       endif
       
   case (100)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%dapd = ivar2d(gref,l)
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%dapd = int(sibr_vd(v))
              enddo
           enddo
       endif

   case (101)
       if (var_exist) then
           allocate(ivar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, ivar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%dapdaf = ivar2d(gref,l)
              enddo
           enddo
           deallocate(ivar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%dapdaf = int(sibr_vd(v))
              enddo
           enddo
       endif
       
   case (102)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%gdd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%gdd = sibr_vd(v)
              enddo
           enddo
       endif

   case (103)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%phent%seed_pool = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%phent%seed_pool = sibr_vd(v)
              enddo
           enddo
       endif

    case (110)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'poolpft_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          status = nf90_inquire_variable(ncid, varid, ndims=numdims)

          if (numdims .eq. 3) then
             allocate(dvar3d(nsibf,nlu,npoolpft))
             status = nf90_get_var(ncid, varid, dvar3d)
             if (status .ne. nf90_noerr) then
                print*,'Error reading restart variable: ', &
                    trim(sibr_vname(v))
                stop
             endif

             do i=1, subcount
                if (nsibr_flag) then
                    gref = subset(i)
                else
                    gref = i
                endif

                do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%poollt%poolpft(:) = dvar3d(gref,l,:)
                enddo
             enddo
             deallocate(dvar3d)
          endif

          if (numdims .eq. 4) then
              allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
              status = nf90_get_var(ncid, varid, dvar4d)
              if (status .ne. nf90_noerr) then
                   print*,'Error reading poolpft_lay.'
                   stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                  do l=1,sib%g(i)%g_nlu
                     do p=1,npoolpft
                        sib%g(i)%l(l)%poollt%poolpft_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                         sib%g(i)%l(l)%poollt%poolpft(p) = &
                              sum(dvar4d(gref,l,p,:))
                       enddo
                    enddo
                enddo
                deallocate(dvar4d)
             endif
          else  !var does not exist
           do i=1, subcount
               do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (111)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_dgain.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%poolpft_dgain(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft_dgain(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (112)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_dloss.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%poolpft_dloss(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft_dloss(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (113)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_assim', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading gain_assim.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%gain_assim = dvar3d(gref,l,:)
              enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%gain_assim(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (114)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_seed', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading gain_seed.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%gain_seed = dvar3d(gref,l,:)
              enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%gain_seed(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (115)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_gresp', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolpft))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_gresp.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%loss_gresp = dvar3d(gref,l,:)
              enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_gresp(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (116)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_mresp_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_mresp_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_mresp_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_mresp_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (117)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_trans_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_trans_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_trans_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_trans_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (118)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_fire_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_fire_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_fire_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_fire_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (119)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_grz', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolcan))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_grz.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%loss_grz(:) = dvar3d(gref,l,:)
             enddo
           enddo
           deallocate(dvar3d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_grz(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (120)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_hrvst_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_hrvst_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   sib%g(i)%l(l)%poollt%loss_hrvst_lay(p,:) = dvar4d(gref,l,p,:)
                enddo
             enddo
           enddo
           deallocate(dvar4d)
             
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_hrvst_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif
       
   case (121)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_totgain(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_totgain(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

   case (122)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_totloss(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_totloss(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

    case (123)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_lay.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%poolpft_lay(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%poolpft_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (124)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_grzc13', varid)
          var_exist = (status .eq. nf90_noerr)
       endif

       if (var_exist) then
          allocate(dvar3d(nsibf,nlu,npoolcanc13))
          status = nf90_get_var(ncid, varid, dvar3d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading loss_grzc13.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

             do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%poollt%loss_grzc13(:) = dvar3d(gref,l,:)
             enddo
           enddo
           deallocate(dvar3d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%loss_grzc13(:) = sibr_vd(v)
               enddo
           enddo
       endif

   case (125)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_equib(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_equib(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (126)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%poollt%poolpftmin(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%poollt%poolpftmin(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (127)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_init(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_init(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (128)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equiblt%poolpft_end(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equiblt%poolpft_end(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (130)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'poollu_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
          status = nf90_inquire_variable(ncid, varid, ndims=numdims)

         if (numdims .eq. 3) then
             allocate(dvar3d(nsibf,nlu,npoollu))
             status = nf90_get_var(ncid, varid, dvar3d)
             if (status .ne. nf90_noerr) then
                print*,'Error reading restart variable: ', &
                    trim(sibr_vname(v))
                stop
             endif

             do i=1, subcount
                if (nsibr_flag) then
                    gref = subset(i)
                else
                    gref = i
                endif

                do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%pooldt%poollu(:) = dvar3d(gref,l,:)
                enddo
             enddo
             deallocate(dvar3d)
          endif

          if (numdims .eq. 4) then
              allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
              status = nf90_get_var(ncid, varid, dvar4d)
              if (status .ne. nf90_noerr) then
                   print*,'Error reading poollu_lay.'
                   stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                  do l=1,sib%g(i)%g_nlu
                     do p=1,npoollu
                        sib%g(i)%l(l)%pooldt%poollu_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                         sib%g(i)%l(l)%pooldt%poollu(p) = &
                              sum(dvar4d(gref,l,p,:))
                       enddo
                    enddo
                enddo
                deallocate(dvar4d)
             endif
            
       else  !var does not exist
           do i=1, subcount
               do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%poollu(:) = sibr_vd(v)
               enddo
           enddo
       endif

    case (131)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_dgain.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%poollu_dgain(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollu_dgain(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (132)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_dloss.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%poollu_dloss(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollu_dloss(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (133)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_grz_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_grz_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_grz_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_grz_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (134)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_hrvst_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_hrvst_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_hrvst_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_hrvst_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (135)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_transl_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_transl_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_transl_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_transl_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (136)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'gain_transd_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading gain_transd_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%gain_transd_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%gain_transd_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (137)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_fire_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_fire_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_fire_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_fire_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (138)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_resp_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_resp_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_resp_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_resp_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

    case (139)
       if (.not. var_exist) then
          status = nf90_inq_varid(ncid, 'loss_trans_lay', varid)
          var_exist = (status .eq. nf90_noerr)
       endif
       
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading loss_trans_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%loss_trans_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)
            
       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%loss_trans_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif
       
   case (140)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_totgain(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_totgain(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

   case (141)
       !if (.not. eqclear_switch) then
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_totloss(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_totloss(:) = sibr_vd(v)
                 enddo
              enddo
          endif
      !endif !not eqclear_switch

    case (142)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%poollu_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%poollu_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

   case (143)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_equib(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_equib(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (144)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_init(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_init(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (145)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%equibdt%poollu_end(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%equibdt%poollu_end(:) = sibr_vd(v)
                 enddo
              enddo
          endif

   case (147)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%kiecps = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%kiecps = sibr_vd(v)
              enddo
           enddo
       endif

   case (148)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%d13cca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%d13cca = sibr_vd(v)
              enddo
           enddo
       endif

   case (149)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcassim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcassim = sibr_vd(v)
              enddo
           enddo
       endif

   case (150)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (151)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13assimd = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13assimd = sibr_vd(v)
              enddo
           enddo
       endif

   case (152)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_nvegc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_nvegc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (153)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_leafc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_leafc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (154)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_rootc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_rootc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (155)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_growc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_growc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (156)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_mntnc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_mntnc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (157)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_autoc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_autoc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (158)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_soilc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_soilc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (159)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%pooldt%resp_hetc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%pooldt%resp_hetc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (160)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%resp_grzc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%poollt%resp_grzc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (161)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoolpft))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%poollt%rcpoolpft(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%poollt%rcpoolpft(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (162)
       if (var_exist) then
          allocate(dvar4d(nsibf,nlu,npoolpft,nsoil))
          status = nf90_get_var(ncid, varid, dvar4d)
          if (status .ne. nf90_noerr) then
               print*,'Error reading poolpft_dgain.'
               stop
          endif

          do i=1, subcount
             if (nsibr_flag) then
                 gref = subset(i)
             else
                 gref = i
             endif

              do l=1,sib%g(i)%g_nlu
                 do p=1,npoolpft
                    sib%g(i)%l(l)%poollt%rcpoolpft_lay(p,:) = &
                          dvar4d(gref,l,p,:)
                  enddo
               enddo
           enddo
           deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%poollt%rcpoolpft_lay(:,:) = sibr_vd(v)
               enddo
           enddo
       endif

   case (163)
          if (var_exist) then
              allocate(dvar3d(nsibf,nlu,npoollu))
              status = nf90_get_var(ncid, varid, dvar3d)
              if (status .ne. nf90_noerr) then
                 print*,'Error reading restart variable: ', &
                     trim(sibr_vname(v))
                 stop
              endif

              do i=1, subcount
                 if (nsibr_flag) then
                     gref = subset(i)
                 else
                     gref = i
                 endif

                 do l=1,sib%g(i)%g_nlu
                     sib%g(i)%l(l)%pooldt%rcpoollu(:) = dvar3d(gref,l,:)
                 enddo
              enddo
              deallocate(dvar3d)
          else
              do i=1, subcount
                 do l=1,sib%g(i)%g_nlu
                    sib%g(i)%l(l)%pooldt%rcpoollu(:) = sibr_vd(v)
                 enddo
              enddo
          endif

    case (164)
       if (var_exist) then
           allocate(dvar4d(nsibf,nlu,npoollu,nsoil))
           status = nf90_get_var(ncid, varid, dvar4d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading poollu_lay.'
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

               do l=1,sib%g(i)%g_nlu
                  do p=1,npoollu
                     sib%g(i)%l(l)%pooldt%rcpoollu_lay(p,:) = &
                              dvar4d(gref,l,p,:)
                   enddo
               enddo
            enddo
            deallocate(dvar4d)

       else  !var does not exist
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                sib%g(i)%l(l)%pooldt%rcpoollu_lay(:,:) = sibr_vd(v)
              enddo
           enddo
       endif

   case (165)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcpoolfire = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcpoolfire = sibr_vd(v)
              enddo
           enddo
       endif

   case (166)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%rcassimfac = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%rcassimfac = sibr_vd(v)
              enddo
           enddo
       endif

   case (167)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%poolemistotC = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%poolemistotC = sibr_vd(v)
              enddo
           enddo
       endif

   case (168)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%poolemisc13 = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%poolemisc13 = sibr_vd(v)
              enddo
           enddo
       endif

   case (170)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c12ca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c12ca = sibr_vd(v)
              enddo
           enddo
       endif

   case (171)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13ca = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13ca = sibr_vd(v)
              enddo
           enddo
       endif

   case (172)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c12assim = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c12assim = sibr_vd(v)
              enddo
           enddo
       endif

   case (173)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c13resptot = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c13resptot = sibr_vd(v)
              enddo
           enddo
       endif

   case (174)
       if (var_exist) then
           allocate(dvar2d(nsibf,nlu))
           status = nf90_get_var(ncid, varid, dvar2d)
           if (status .ne. nf90_noerr) then
              print*,'Error reading restart variable: ', &
                  trim(sibr_vname(v))
              stop
           endif

           do i=1, subcount
              if (nsibr_flag) then
                  gref = subset(i)
              else
                  gref = i
              endif

              do l=1,sib%g(i)%g_nlu
                  sib%g(i)%l(l)%fract%c12resptot = dvar2d(gref,l)
              enddo
           enddo
           deallocate(dvar2d)
       else
           do i=1, subcount
              do l=1,sib%g(i)%g_nlu
                 sib%g(i)%l(l)%fract%c12resptot = sibr_vd(v)
              enddo
           enddo
       endif

   end select
enddo  !v=1,sibr_nvar

status = nf90_close(ncid)

end subroutine restart_read


