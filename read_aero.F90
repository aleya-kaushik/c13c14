#include "nc_util.h"

! Opens and reads in SiB aerodynamic tables.
subroutine read_aero ()

    use kinds
    use netcdf

    use module_sibconst, only: npft
    use module_io, only: aero_file
    use module_param

    !...local variables
    integer(i4) ::  ncid         !  file id#
    integer(i4) ::  dimid        !  dimension id#
    integer(i4) ::  varid        !  variable id#
    character(len=12) ::  name   !  variable name

    integer(i4) :: temp
    real(r4), dimension(:,:,:), allocatable :: &
         aero_zo, aero_zp, aero_rbc, aero_rdc

    !...open file
    print*,''
    print *, 'Reading Aerodynamic File: '
    print *, '  ',trim(aero_file)
    CHECK( nf90_open ( trim(aero_file), nf90_nowrite, ncid ) )

    !...check npft
    CHECK( nf90_inq_dimid ( ncid, 'npft', dimid ) )
    CHECK( nf90_inquire_dimension ( ncid, dimid, name, temp ) )

    if (temp /= npft) then
        print*,''
        print('(a)'),'!!!Aerodynamic file does not match simulation!!!'
        print('(a,i4,a,i4)'),'  Aero file npft: ',temp,' Sim npft: ',npft
        print*,''
        stop
    endif

    CHECK( nf90_inq_dimid ( ncid, 'ngrid', dimid ) )
    CHECK( nf90_inquire_dimension ( ncid, dimid, name, ngrid ) )

    allocate(LAIgrid(ngrid), fVCovergrid(ngrid))
    allocate(aerovar(npft,ngrid,ngrid))
    allocate(aero_zo(npft,ngrid,ngrid),aero_zp(npft,ngrid,ngrid), &
        aero_rbc(npft,ngrid,ngrid),aero_rdc(npft,ngrid,ngrid))

    !...read in values
    ENSURE_VAR ( ncid, 'laigrid', varid )
    CHECK( nf90_get_var ( ncid, varid, laigrid ) )

    ENSURE_VAR ( ncid, 'fvcovergrid', varid )
    CHECK( nf90_get_var( ncid, varid, fvcovergrid ) )

    ENSURE_VAR ( ncid, 'aero_zo', varid )
    CHECK( nf90_get_var ( ncid, varid, aero_zo ) )

    ENSURE_VAR ( ncid, 'aero_zp', varid )
    CHECK( nf90_get_var ( ncid, varid, aero_zp ) )

    ENSURE_VAR ( ncid, 'aero_rbc', varid )
    CHECK( nf90_get_var ( ncid, varid, aero_rbc ) )

    ENSURE_VAR ( ncid, 'aero_rdc', varid )
    CHECK( nf90_get_var ( ncid, varid, aero_rdc ) )
    CHECK( nf90_close(ncid) )

    !------------Assign values from tables into data structures--------
    aerovar%zo = aero_zo
    aerovar%zp_disp = aero_zp
    aerovar%rbc = aero_rbc
    aerovar%rdc = aero_rdc

end subroutine read_aero

