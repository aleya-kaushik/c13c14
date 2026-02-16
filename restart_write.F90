#include "nc_util.h"

subroutine restart_write( equib )
    !---------------------------------------------
    ! This subroutine creates a netcdf restart file
    ! for all prognostic variables.
    !---------------------------------------------

    use kinds
    use netcdf

    use module_io
    use module_pftinfo, only: &
        clen, pft_name, pft_ref
    use module_sibconst, only: &
        slen, npft, nlu, &
        npoolpft, npoollu, npoolcan, &
        nsoil, nsnow, ntot, &
        subcount, sublatsib, sublonsib, &
        subsitename, subpref, sublarea, &
        npoolcanc13
    use module_sib, only: sib

    implicit none
    
    !---input variables
    logical :: equib  !equilibrium pool IC file

    !---netcdf variables
    integer(i4) :: ncid, status
    integer(i4) :: varid
    
    integer(i4) :: nsibdid, npftdid, nludid, &
                   nppdid, nlpdid, ncpdid, &
                   nsdid, nsnowdid, ntotdid, &
                   nclendid, nslendid, &
                   ncpdidc13
    integer(i4) :: lonid, latid, &
                   pnameid, prefid, snameid, &
                   areaid, pftrefid

    integer(i4), dimension(:,:), allocatable :: ivar2d
    !real(r4), dimension(:), allocatable :: tvar
    !real(r4), dimension(:,:), allocatable :: tvar2d
    !real(r4), dimension(:,:,:), allocatable :: tvar3d
    !real(r4), dimension(:,:,:,:), allocatable :: tvar4d

    real(r8), dimension(:), allocatable :: dvar
    real(r8), dimension(:,:), allocatable :: dvar2d
    real(r8), dimension(:,:,:), allocatable :: dvar3d
    real(r8), dimension(:,:,:,:), allocatable :: dvar4d
    
    !---local variables
    integer(i4) :: i,l,p,v
    logical :: savesite
    logical, dimension(:), allocatable :: dovar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Determine if saving site names
if (subsitename(1) .ne. '') then
    savesite = .true.
else
    savesite = .false.
endif
    
! Write out the data to a file
allocate(dovar(sibr_nvar))
if (equib) then
    status = nf90_create( requib_filename, ncid=ncid, &
             cmode=or(nf90_clobber,nf90_64bit_offset))
    if (status .ne. nf90_noerr) then
        print('(a)'),'Error creating requib file: '
        print('(2a)'),'  ',requib_filename
        stop
    else
        print('(a)'),'Writing Equilibrium Restart File: '
        print('(2a)'),'  ', trim(requib_filename)
    endif
    dovar = sibreq_doref
else
    status = nf90_create( restart_filename, ncid=ncid, &
             cmode=or(nf90_clobber,nf90_64bit_offset))
    if (status .ne. nf90_noerr) then
        print('(a)'),'Error creating file: '
        print('(2a)'),'  ',restart_filename
        stop
    !else
    !    print('(a)'),'Writing Restart File: '
    !    print('(2a)'),'  ', trim(restart_filename)
    endif
    dovar = sibr_doref
endif


!...dimensions...
CHECK( nf90_def_dim( ncid, 'nsib', subcount, nsibdid ) )
CHECK( nf90_def_dim( ncid, 'nlu', nlu, nludid ) )
CHECK( nf90_def_dim( ncid, 'npft', npft, npftdid ) )
CHECK( nf90_def_dim( ncid, 'npoolcan', npoolcan, ncpdid ) )
CHECK( nf90_def_dim( ncid, 'npoolcanc13', npoolcanc13, ncpdidc13 ) )
CHECK( nf90_def_dim( ncid, 'npoolpft', npoolpft, nppdid ) )
CHECK( nf90_def_dim( ncid, 'npoollu', npoollu, nlpdid ) )
CHECK( nf90_def_dim( ncid, 'nsnow', nsnow, nsnowdid ) )
CHECK( nf90_def_dim( ncid, 'nsoil', nsoil, nsdid ) )
CHECK( nf90_def_dim( ncid, 'ntot', ntot, ntotdid ) )
CHECK( nf90_def_dim( ncid, 'clen', clen, nclendid ) )
if (savesite) CHECK( nf90_def_dim( ncid, 'slen', slen, nslendid ) )

!...define and write info variables...
CHECK( nf90_def_var( ncid, 'lonsib', nf90_float, nsibdid, lonid ) )
CHECK( nf90_def_var( ncid, 'latsib', nf90_float, nsibdid, latid ) )
CHECK( nf90_def_var( ncid, 'pftnames', nf90_char, (/nclendid,npftdid/), pnameid))
CHECK( nf90_def_var( ncid, 'pftrefs', nf90_int, (/npftdid/), prefid ) )
if (savesite) then
    CHECK( nf90_def_var( ncid, 'site_names', nf90_char, (/nslendid,nsibdid/), snameid))
endif
CHECK( nf90_def_var( ncid, 'lu_area', nf90_float, (/nsibdid,nludid/), areaid ) )
CHECK( nf90_def_var( ncid, 'lu_pftref', nf90_int, (/nsibdid,nludid/), pftrefid ) )

CHECK( nf90_enddef(ncid) )
CHECK( nf90_put_var( ncid, lonid, sublonsib ) )
CHECK( nf90_put_var( ncid, latid, sublatsib ) )
CHECK( nf90_put_var( ncid, pnameid, pft_name ) ) 
CHECK( nf90_put_var( ncid, prefid, pft_ref ) )
if (savesite) CHECK( nf90_put_var( ncid, snameid, subsitename ) )
CHECK( nf90_put_var( ncid, areaid, sublarea ) )
CHECK( nf90_put_var( ncid, pftrefid, subpref ) )


!...define variables...
do v=1, sibr_nvar

   if (dovar(v)) then
      select case (sibr_vref(v))

      case (1)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%tmd
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (2)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%seas_precip
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

      case (3)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%seas_tm
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (4)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%clim_cupr
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

      case (5)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%clim_precip
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (6)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%clim_tm
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)

     case (10)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%eacas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (11)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%shcas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (12)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tcas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (13)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tkecas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (14)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tc
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (15)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cast%tcmin
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (20)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2m
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
         
      case (21)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2cas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (22)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2c
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (23)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2i
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (24)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%pco2s
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
    case (25)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar(subcount))
          do i=1,subcount
             dvar(i) = sib%g(i)%gprogt%pcosm
          enddo
          status = nf90_put_var(ncid,varid,dvar)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar)
          
      case (26)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cost%coscasp
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

       case (27)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%cost%coscas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)             

      case (28)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2cas
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (29)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2c
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (30)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
              
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2i
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (31)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF
              
          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%co2s
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (35)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fluxt%ra
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (36)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fluxt%rb
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (37)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%rst
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (40)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (41)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%assimd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (42)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_nveg
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (43)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_leaf
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (44)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_root
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (45)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grow
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (46)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_mntn
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (47)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_auto
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (48)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_soil
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (49)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_het
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (50)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grz
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (60)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%capacc_liq
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (61)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%capacc_snow
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (62)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrost%capacg
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (63)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%clim_pawfrw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (64)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%hydrovt%clim_tawfrw
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (70)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = int(sib%g(i)%l(l)%sscolt%nsl)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (71)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nsnowdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,nsnow))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%dz(-nsnow+1:0)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (72)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ntotdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,ntot))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%td(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (73)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ntotdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,ntot))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%www_liq(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (74)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ntotdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,ntot))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%sscolt%www_ice(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (80)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%co2t%clim_assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (81)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%vegt%clim_lai
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)


      case (82)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (83)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_assimsm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (84)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_pr
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (85)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prsm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (86)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prsdoy
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (87)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_prcdoy
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)
          
      case (88)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_tawftop
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (89)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_tm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (90)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenc_climp
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (91)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenc_laimax
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (92)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenc_laimin
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (93)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_env
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (94)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_wa
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (95)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phenave_wacsm
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (96)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%phen_pi
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)          

      case (97)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%phen_istage
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (100)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%dapd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)

      case (101)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_int, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(ivar2d(subcount,nlu))
          ivar2d(:,:) = izero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                ivar2d(i,l) = sib%g(i)%l(l)%phent%dapdaf
             enddo
          enddo
          status = nf90_put_var(ncid,varid,ivar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(ivar2d)
          
      case (102)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%gdd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (103)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%phent%seed_pool
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (110)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)
          
      case (111)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_dgain(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (112)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_dloss(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (113)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%gain_assim(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (114)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%gain_seed(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)
          
      case (115)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_gresp(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (116)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_mresp_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (117)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_trans_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (118)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_fire_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (119)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ncpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolcan))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_grz(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)
          
      case (120)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%loss_hrvst_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (121)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_totgain(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (122)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_totloss(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (123)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%poolpft_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (124)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,ncpdidc13/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolcanc13))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%loss_grzc13(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (125)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_equib(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (126)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%poolpftmin(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (127)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_init(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (128)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equiblt%poolpft_end(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (130)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                    dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_lay(p,:)
                 enddo
              enddo   
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (131)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_dgain(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (132)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_dloss(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (133)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_grz_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (134)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_hrvst_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (135)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_transl_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (136)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%gain_transd_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (137)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_fire_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (138)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_resp_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (139)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%loss_trans_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)
          
      case (140)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_totgain(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (141)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_totloss(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (142)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%poollu_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (143)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_equib(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (144)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_init(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (145)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%equibdt%poollu_end(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (147)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%kiecps
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (148)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%d13cca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (149)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcassim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (150)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (151)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13assimd
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (152)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_nvegc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (153)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_leafc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (154)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_rootc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (155)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_growc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (156)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_mntnc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (157)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_autoc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (158)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_soilc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (159)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%pooldt%resp_hetc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (160)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%poollt%resp_grzc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (161) 
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoolpft))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%poollt%rcpoolpft(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (162) 
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nppdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoolpft,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoolpft
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%poollt%rcpoolpft_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (163)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar3d(subcount,nlu,npoollu))
          dvar3d(:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar3d(i,l,:) = sib%g(i)%l(l)%pooldt%rcpoollu(:)
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar3d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar3d)

      case (164) 
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid,nlpdid,nsdid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar4d(subcount,nlu,npoollu,nsoil))
          dvar4d(:,:,:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                do p=1,npoollu
                   dvar4d(i,l,p,:) = sib%g(i)%l(l)%pooldt%rcpoollu_lay(p,:)
                enddo
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar4d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar4d)

      case (165)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcpoolfire
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (166)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%rcassimfac
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (167)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%poolemistotC
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (168)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%poolemisc13
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (170)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c12ca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (171)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13ca
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (172)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c12assim
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (173)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c13resptot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      case (174)
          CHECK( nf90_redef(ncid) )
          status = nf90_def_var(ncid, trim(sibr_vname(v)), nf90_double, &
                        (/nsibdid,nludid/), varid)
          IF (status .ne. nf90_noerr) THEN
              print*,'Error Defining Variable: ',sibr_vname(v)
              stop
          ENDIF

          CHECK( nf90_enddef(ncid) )
          allocate(dvar2d(subcount,nlu))
          dvar2d(:,:) = dzero
          do i=1,subcount 
             do l=1,sib%g(i)%g_nlu
                dvar2d(i,l) = sib%g(i)%l(l)%fract%c12resptot
             enddo
          enddo
          status = nf90_put_var(ncid,varid,dvar2d)
          IF (status .ne. nf90_noerr) THEN
             print*,'Error Putting Variable: ',sibr_vname(v)
             stop
          ENDIF
          deallocate(dvar2d)

      !case default
      !   print*,'Invalid Restart Variable: ', sibr_vname(v)
      !   stop

      end select
   endif !dovar
enddo


!...close the file
CHECK( nf90_close( ncid ) )


end subroutine restart_write
