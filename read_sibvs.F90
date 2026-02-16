#include "nc_util.h"

! Reads in vegetation structure file.
subroutine read_sibvs()

use kinds
use netcdf
use module_sibconst, only: &
   nsib, npft, nlu, &
   single_pt, &
   latsib, lonsib, &
   slen, sitenamesib
use module_io, only: vs_file
use module_sibvs

implicit none

!...variables for regional/global runs
integer(i4) :: ncid, status
integer(i4) :: dimid, varid, dimlen

integer(byte), dimension(:), allocatable :: gnlu
real(r4), dimension(:,:), allocatable :: larea
integer(byte), dimension(:,:), allocatable :: pftref
real(r4), dimension(:,:), allocatable :: sandfrac, clayfrac
real(r4), dimension(:,:), allocatable :: soref_vis, soref_nir

!...local variables
integer :: i,l
integer :: finsib,finpft,finlu
integer :: sibpt
character(len=100) :: trash
logical :: exist

!--------------------------------
! Allocate data
allocate( sibvs(nsib) )
allocate( latsib(nsib) )
allocate( lonsib(nsib) )
allocate( sitenamesib(nsib) )

! Initialize sibvs data
call init_sibvs(nlu, nsib, sibvs)

!Read in sibvs file
inquire(file=trim(vs_file),exist=exist)
IF (.not. exist) THEN
   print*,'Missing/Non-Existent SiB VS File!'
   print*,'Please check vs_file in namel_sibdrv.'
   print*,''
   STOP
ELSE
   print*,'Reading SiB VS File: '
   print*,'  ',trim(vs_file)
ENDIF

!...Expecting text file for single point runs
if (single_pt) then
    ! Open the file
    open(unit=32,file=trim(vs_file),form='formatted')

    read(32,*) finsib
    if (finsib /= nsib) then
        print*,''
        print('(a)'),'!!!SiB VS file does not match simulation!!!'
        print('(a,i6,a,i6)'),' VS file nsib: ',finsib, ' Sim nsib: ',nsib
        print*,''
        stop
    endif

    read(32,*) finpft
    if (finpft /= npft) then
        print*,''
        print*,'!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!'
        print*,'!!!SiB VS file does not match simulation!!!'
        print('(a,i4,a,i4)'),'   VS file npft: ',finpft, ' Sim npft: ',npft
        print*,''
        stop
    endif

    read(32,*) finlu
    if (finlu /= nlu) then
        print*,''
        print*,'!!!SiB VS file does not match simulation!!!'
        print('(a,i4,a,i4)'),'   VS file nlu: ',finlu, ' Sim nlu: ',nlu
        print*,''
        stop
    endif
    
    !...read in grid point location and soil information
    read(32,*) latsib(1)
    read(32,*) lonsib(1)
    read(32,*) sibpt
    read(32,*) sitenamesib

    !...read in structure information
    read(32,*) sibvs(1)%gnlu

    do i=1,sibvs(1)%gnlu
       read(32,*) trash
       read(32,*) sibvs(1)%pftref(i)
       read(32,*) sibvs(1)%larea(i)
       read(32,*) sibvs(1)%sandfrac(i)
       read(32,*) sibvs(1)%clayfrac(i)
       read(32,*) sibvs(1)%soref_vis(i)
       read(32,*) sibvs(1)%soref_nir(i)       
    enddo

    close(32)

else
    !Read in sibvs file
    !...Expecting netcdf file for regional/global runs
    CHECK( nf90_open(trim(vs_file),nf90_nowrite,ncid) )

    CHECK( nf90_inq_dimid( ncid, 'nsib', dimid ) )
    CHECK( nf90_inquire_dimension( ncid, dimid, len=dimlen ) )
    if ( dimlen /= nsib ) then         
        print*,''
        print('(a)'),'!!!SiB VS file does not match simulation!!!'
        print*,'  VS file nsib: ',dimlen,'Sim nsib: ', nsib
        stop
    endif

    CHECK( nf90_inq_dimid( ncid, 'npft', dimid ) )
    CHECK( nf90_inquire_dimension( ncid, dimid, len=dimlen ) )
    if ( dimlen /= npft) then
        print*,''
        print('(a)'),'   !!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!'
        print('(a)'),'   !!!SiB VS file does not match simulation!!!'
        print('(a,i4,a,i4)'),  &
                     '     File npft: ',dimlen, ' Sim npft: ',npft
        print*,''
        stop
    endif

    CHECK( nf90_inq_dimid( ncid, 'nlu', dimid ) )
    CHECK( nf90_inquire_dimension( ncid, dimid, len=dimlen ) )
    if ( dimlen /= nlu) then
        print*,''
        print('(a)'),'!!!SiB VS file does not match simulation!!!'
        print*,'   File nlu: ',dimlen, ' Sim nlu: ',nlu
        stop
    endif

    allocate(gnlu(nsib))
    allocate(larea(nsib,nlu),pftref(nsib,nlu))
    allocate(sandfrac(nsib,nlu),clayfrac(nsib,nlu))
    allocate(soref_vis(nsib,nlu),soref_nir(nsib,nlu))

    ENSURE_VAR( ncid, 'latsib', varid )
    CHECK( nf90_get_var( ncid, varid, latsib ) )
    ENSURE_VAR( ncid, 'lonsib', varid )
    CHECK( nf90_get_var( ncid, varid, lonsib ) )

    status = nf90_inq_varid(ncid, 'site_names', varid )
    IF (status == nf90_noerr) THEN
        CHECK( nf90_inq_dimid( ncid, 'slen', dimid ) )
        CHECK( nf90_inquire_dimension( ncid, dimid, len=dimlen ) )
        IF (dimlen .ne. slen) THEN
            print*,''
            print('(a)'),'!!!SiB VS file does not match simulation!!!'
            print*,'   File slen: ',dimlen, ' Sim slen: ',slen
            stop
        ENDIF
        CHECK( nf90_get_var( ncid, varid, sitenamesib ) )
    ELSE
       sitenamesib(:) = ''
    ENDIF

    ENSURE_VAR( ncid, 'g_nlu', varid )
    CHECK( nf90_get_var( ncid, varid, gnlu ) )
    ENSURE_VAR( ncid, 'larea', varid )
    CHECK( nf90_get_var ( ncid, varid, larea ) )
    ENSURE_VAR( ncid, 'pftref', varid )
    CHECK( nf90_get_var( ncid, varid, pftref ) )
    ENSURE_VAR( ncid, 'sand_frac', varid )
    CHECK( nf90_get_var( ncid, varid, sandfrac ) )
    ENSURE_VAR( ncid, 'clay_frac', varid )
    CHECK( nf90_get_var( ncid, varid, clayfrac ) )
    ENSURE_VAR( ncid, 'soref_vis', varid )
    CHECK( nf90_get_var( ncid, varid, soref_vis ) )
    ENSURE_VAR( ncid, 'soref_nir', varid )
    CHECK( nf90_get_var( ncid, varid, soref_nir ) )

    !...close file
    CHECK( nf90_close(ncid) )

    !...put the variables into the sibvs stucture
    do i=1, nsib
       sibvs(i)%gnlu = gnlu(i)
       do l=1,sibvs(i)%gnlu
           sibvs(i)%larea(l) = larea(i,l)
           sibvs(i)%pftref(l) = pftref(i,l)
           sibvs(i)%sandfrac(l) = sandfrac(i,l)
           sibvs(i)%clayfrac(l) = clayfrac(i,l)
           sibvs(i)%soref_vis(l) = soref_vis(i,l)
           sibvs(i)%soref_nir(l) = soref_nir(i,l)           
       enddo !l=1,gnlu
    enddo !i=1,nsib

endif

end subroutine read_sibvs

!-------------------------------------
! Initializes the vegetation structure file.
subroutine init_sibvs(nlu, nsib, sibvs)

use kinds
use module_sibvs, only: sib_vs_vars
implicit none

!...input variables
integer(i4), intent(in) :: nlu, nsib
type(sib_vs_vars), intent(inout), dimension(nsib) :: sibvs

!...local variables
integer(i4) :: i

sibvs(:)%gnlu = bzero
do i=1,nsib
   allocate(sibvs(i)%larea(nlu))
   sibvs(i)%larea(:)=rzero
   allocate(sibvs(i)%pftref(nlu))
   sibvs(i)%pftref(:) = izero

   allocate(sibvs(i)%sandfrac(nlu))
   sibvs(i)%sandfrac(:) = dzero
   allocate(sibvs(i)%clayfrac(nlu))
   sibvs(i)%clayfrac(:) = dzero
   allocate(sibvs(i)%soref_vis(nlu))
   sibvs(i)%soref_vis(:) = dzero
   allocate(sibvs(i)%soref_nir(nlu))
   sibvs(i)%soref_nir(:) = dzero
enddo

end subroutine init_sibvs
