#include "nc_util.h"

!--------------------------------------------
!Includes routines to write data to 
!   both global and land unit output files.
!---------------------------------------------

subroutine diagnostic_writeg(  &
    npts, ntpersave, nsaveperout, &
    numvars, fileid, write_step,  &
    varrefs, varids, sif_satcount, &
    outsave)

    use netcdf
    use kinds
    use module_io, only: missing
    
    implicit none
    
    ! input variables
    integer(i4), intent(in) :: npts
    integer(i4), intent(in) :: ntpersave    ! # of timesteps per saved steps
    integer(i4), intent(in) :: nsaveperout  ! # of saved steps per output
    integer(i4), intent(in) :: numvars      ! # of variables written to file
    integer(i4), intent(in) :: fileid       ! file id #
    integer(i4), intent(in) :: write_step   ! step to write
    
    integer(i4), dimension(numvars), intent(in) :: varrefs   !variable ref #s
    integer(i4), dimension(numvars), intent(in) :: varids   ! variable id #s
    integer(i4), dimension(npts, nsaveperout, 2), intent(in) :: sif_satcount
    real(r8), dimension(npts, numvars, nsaveperout), intent(in) :: &
          outsave ! variable values

    ! netcdf variables
    integer(i4) :: status

    ! local variables
    integer(i4) :: nv  !index variable
    real(r4), dimension(npts,nsaveperout) :: outtemp

    !--------------------------------
    ! write out data variables
    do nv = 1, numvars
       IF (varrefs(nv) .eq. 909) THEN
          where(sif_satcount(:,:,1) .gt. 0) 
             outtemp(:,:) = real(outsave(:,nv,:)) * ntpersave / sif_satcount(:,:,1)
          elsewhere
             outtemp(:,:) = missing
          endwhere

       ELSEIF (varrefs(nv) .eq. 910) THEN
          where(sif_satcount(:,:,2) .gt. 0)
             outtemp(:,:) = real(outsave(:,nv,:)) * ntpersave / sif_satcount(:,:,2)
          elsewhere
             outtemp(:,:) = missing
          endwhere

       ELSE
           outtemp = real(outsave(:,nv,:))
       ENDIF
       
       status = nf90_put_var(fileid,varids(nv),outtemp, &
                start=(/1,write_step/), &
                count=(/npts,nsaveperout/))
        if (status .ne. nf90_noerr) then
              print*,'Error Writing Grid Cell Diagnostic Output'
              print*,'Variable ID: ',varids(nv)
              print*,'Variable Min/Max: ', &
                      minval(outtemp),maxval(outtemp)
              print*,'Time-Steps: ', &
                     write_step,write_step+nsaveperout-1
              print*,'Stopping.'
              stop
        endif
    enddo

end subroutine diagnostic_writeg

!-----------------------------------------------------------------------

subroutine diagnostic_writelu( &
    npts, nsaveperout, fileid, &
    numvars, nlevs, varids, &
    write_step, outsave)

    use netcdf
    use kinds

    ! input variables
    integer(i4), intent(in) :: npts        ! # of points to save
    integer(i4), intent(in) :: nsaveperout ! # of timesteps per output
    integer(i4), intent(in) :: fileid      ! file id #
    integer(i4), intent(in) :: numvars     ! # of variables written to file
    integer(i4), intent(in) :: nlevs       ! # levels
    integer(i4), dimension(numvars), intent(in) :: varids  ! variable id #s
    integer(i4), intent(in) :: write_step  ! starting timestep #
    real(r8), dimension(npts, nlevs, numvars, nsaveperout), intent(in) :: &
           outsave ! variable values

    ! netcdf variables
    integer(i4) :: status

    ! local variables
    integer(i4) :: nv
    real(r4), dimension(npts,nlevs,nsaveperout) :: outtemp

    !--------------------------------
    ! write out data variables
     do nv = 1, numvars
        outtemp = real(outsave(:,:,nv,:))
        status = nf90_put_var(fileid,varids(nv),outtemp, &
                 start=(/1,1,write_step/), &
                 count=(/npts,nlevs,nsaveperout/))

        if (status .ne. nf90_noerr) then
            print*,'Error Writing Land Unit Diagnostic Output'
            print*,'Variable ID: ',varids(nv)
            print*,'Variable Min/Max: ', &
                    minval(outtemp),maxval(outtemp)
            print*,'Time-Steps: ', &
                   write_step,write_step+nsaveperout-1
            print*,'Stopping.'
            stop
        endif
     enddo


end subroutine diagnostic_writelu

