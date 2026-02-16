subroutine output_closer()

!----------------------------
!Close any open files
!----------------------------

    use kinds
    use module_io, only:  &
        hr_idg, hr_idlu,  &
        qp_idg, qp_idlu,  &
        pbp_idg, pbp_idlu, &
        driverid

    use netcdf

    implicit none

    integer(i4) :: status

    ! close all output files 
    status = nf90_noerr
    if (qp_idg > 0) status = nf90_close( qp_idg )
    IF (status /= nf90_noerr) print*,'Error Closing QP Gridcell File'
    if (qp_idlu > 0) status = nf90_close( qp_idlu )    
    IF (status /= nf90_noerr) print*,'Error Closing QP Land Unit File'

    if (pbp_idg > 0) status = nf90_close( pbp_idg )
    IF (status /= nf90_noerr) print*,'Error Closing PBP Gridcell File'
    if (pbp_idlu > 0) status = nf90_close( pbp_idlu )
    IF (status /= nf90_noerr) print*,'Error Closing PBP Land Unit File'

    if (hr_idg > 0) status = nf90_close( hr_idg ) 
    IF (status /= nf90_noerr) print*,'Error Closing HR Gridcell File'
    if (hr_idlu > 0) status = nf90_close( hr_idlu ) 
    IF (status /= nf90_noerr) print*,'Error Closing HR Land Unit File'

    ! close driver data
    status = nf90_close( driverid )

end subroutine output_closer
