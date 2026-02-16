!===============================================================================
subroutine read_routopts()
!===============================================================================

use module_sibconst, only: &
    nlu, nsoil, nsnow, ntot
use module_io

implicit none

integer(i4) :: count
integer(i4) :: ireftem
real(r4) :: default
character(len=1) :: dortem, doreqtem
character (len=20) :: nametem

    !---------------------------------------------------------------------------
    ! read sib_routopts and count number of restart variables to be output
    !---------------------------------------------------------------------------
 
    sibr_nvar = 0

    print*,'Reading Restart Output File: '
    print*,'  ',trim(out_rinfo)
    open(unit=3,file=trim(out_rinfo),form='formatted')

    !...read in comments
    call read_comments(3)

    do 
       read(3,*,end=1200) dortem, doreqtem, nametem, &
                         ireftem, default
       sibr_nvar = sibr_nvar + 1
    enddo

    1200 continue
    rewind 3
    call read_comments(3)

    allocate(sibr_doref(sibr_nvar))
    sibr_doref(:) = .false.
    allocate(sibreq_doref(sibr_nvar))
    sibreq_doref(:) = .false.

    allocate(sibr_vname(sibr_nvar))
    allocate(sibr_vref(sibr_nvar))
    allocate(sibr_vd(sibr_nvar))

    count=0

    do
       read(3,*,end=944) dortem, doreqtem, nametem, &
                         ireftem, default
       count = count + 1
       sibr_vname(count) = nametem
       sibr_vref(count) = ireftem
       sibr_vd(count) = default
       if (dortem == 't') then
          sibr_doref(count) = .true.
       endif
       if (doreqtem == 't') then
          sibreq_doref(count) = .true.
       endif   
    enddo

    944 continue
    close(3)

    !print('(a,i4)'),  &
    !   '  Number of restart variables: ',  sibr_nvar

end subroutine read_routopts


!==================
subroutine read_comments( fid )

    use kinds

    implicit none

    integer(i4), intent(in) :: fid
    logical :: iscomment1, iscomment2
    character(len=100) :: trash

    iscomment1=.true.
    iscomment2=.true.
    do while ((iscomment1) .or. (iscomment2))
       read(fid,*) trash
       if (index(trash,'*****') .gt. 0) then
          if (iscomment1) then
             iscomment1=.false.
          else
             iscomment2=.false.
          endif
       endif
     enddo

end subroutine read_comments



