
! Opens and reads in C isotope delta-c13 from c_iso_time_series.dat
subroutine read_ciso()

use kinds
use module_isodata
use module_sibconst, only: &
    nisodatayr
!use module_sib, only: &
!    fract_type
use module_io, only: &
    isodata_file, cisoid

integer(i4) :: i
!type(fract_type), intent(inout) :: fract
character(len=100) :: trash
logical :: iscomment1, iscomment2

!integer(i4), intent(in) :: nisodatayr
real(r8), dimension(171) :: &
   isoyrtmp, nhc14tmp, tropc14tmp, &
   shc14tmp, globc13tmp
!real(r8), intent(in) :: &
!   isoyr, nhc14, tropc14, shc14, globc13
!CHARACTER(LEN=30), PARAMETER :: FMT1 = "(F5.2,F5,2,F5.2,F5.2,F5.2)"

!-------------------
!...Initialize the iso variables (done in sibtype_init)
!allocate(isoyr(nisodatayr))
!isoyr(:)=dzero
!allocate(nhc14(nisodatayr))
!nhc14(:)=dzero
!allocate(tropc14(nisodatayr))
!tropc14(:)=dzero
!allocate(shc14(nisodatayr))
!shc14(:)=dzero
!allocate(globc13(nisodatayr))
!globc13(:)=dzero
!fract%isoyr(:)=dzero
!fract%nhc14(:)=dzero
!fract%tropc14(:)=dzero
!fract%shc14(:)=dzero
!fract%globc13(:)=dzero

!...Switch this to namelist file later
!iso_file='/home/ccg/kaushik/sib4/sib4v2_corral-master-lk/input-c13/params/c_iso_time_series.dat'

!...Open file
print*,'Reading Carbon Isotope records File: '
print*,'  ',trim(isodata_file)
open(unit=cisoid,file=trim(isodata_file),form='formatted')

!...top line of file = # yrs data
read(cisoid,*) nisodatayr
!print('(a,i4)'),'nisodatayr: ', nisodatayr

!...ignore text in between lines of ***
iscomment1=.true.
iscomment2=.true.
do while ((iscomment1) .or. (iscomment2))
    read(cisoid,*) trash
    if (index(trash,'****') .gt. 0) then
        if (iscomment1) then
           iscomment1=.false.
        else
           iscomment2=.false.
        endif
    endif
enddo

!...read in data values
read(cisoid,*) trash
read(cisoid,*) trash
do i=1,nisodatayr
   read(cisoid,*) isoyrtmp(i), nhc14tmp(i), &
                  tropc14tmp(i), shc14tmp(i), &
                  globc13tmp(i)
   !print*,''
   !print('(a,f4.2,f4.2)'),'ciso_file isoyr: ', isoyrtmp(i),globc13tmp(i)
   !write(*,*) isoyrtmp(i),globc13tmp(i)
   !print*,''
enddo

close(cisoid)

isoyr=isoyrtmp
nhc14=nhc14tmp
tropc14=tropc14tmp
shc14=shc14tmp 
globc13=globc13tmp

end subroutine read_ciso
