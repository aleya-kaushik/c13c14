
! Sets up grid.
subroutine grid_init( )

use netcdf
use kinds

use module_pftinfo, only: &
   pft_name, pft_num
use module_sibconst
use module_io
use module_sibvs

! local grid variables
integer(i4) :: i,j,l  ! index variables
integer(i4) :: gref

integer(i4) :: pbp_index, sibpt_index   ! index value for pbp's
real(r8)    :: pbp_diff    ! differences between pbp and input
integer(i4), allocatable, dimension(:) :: &
   temp_pbp,  &  ! temporary array of pbp coordinates on subgrid
   temp_sibpt    ! temporary array of pbp coordinates on nsib

integer(i4) :: subgridcount   ! total number of points in subgrid
integer(i4) :: start_count    ! starting count for parallel
integer(i4) :: stop_count     ! stoping count for parallel
integer(i4) :: olength
integer(i4) :: nrresidual ! LK

integer(i4) :: start_count_LK    ! starting count for parallel
integer(i4) :: stop_count_LK     ! stoping count for parallel
integer(i4) :: olength_LK

integer(i4) :: start_count_ES    ! starting count for parallel
integer(i4) :: stop_count_ES     ! stoping count for parallel

integer(i4), dimension(:), allocatable :: &
   subgridref ! index reference for subgrid
integer(i4), dimension(:), allocatable :: &
   index_ref  ! index reference for parallelization

! Variables to do PFT-based splitting
!   integer(i4), dimension(:), allocatable :: &
!       subgridpft ! index reference for number of pfts per subgrid
!   integer(i4), dimension(:), allocatable :: &
!       split_points ! index reference for start/stop split points
!   integer(i4) :: subgridchunk  ! PFT chunk size for each subgrid
!   integer(i4) :: split_index
!   integer(i4) :: pft_count

!-----------------------------------------------------------------------
print*,'Setting Grid'

if (single_pt) then
    !...Single Point Simulations
    if (nsib .ne. 1) then
        print*,'Expecting nsib value of 1.  Stopping.'
        stop
    endif
    print('(a,i5)'), '     SiB points simulated =', nsib

    subcount = 1
    allocate(subset(1))
    subset(1) = 1

    allocate(sublonsib(1))
    allocate(sublatsib(1))
    sublatsib(1) = latsib(1)
    sublonsib(1) = lonsib(1)

    allocate(subsitename(1))
    subsitename(:) = sitenamesib(1)

    allocate(subpref(1,nlu))
    allocate(sublarea(1,nlu))
    subpref(1,:) = sibvs(1)%pftref(:)
    sublarea(1,:) = sibvs(1)%larea(:)
 
    allocate(pbp_lat(1))
    allocate(pbp_lon(1))
    allocate(pbp_sibpt(1))

    pbp_lat(1) = sublatsib(1)
    pbp_lon(1) = sublonsib(1)
    pbp_sibpt(1) = 1

    npbp=1       
    allocate(pbp_sitename(npbp))
    pbp_sitename(1) = subsitename(1)
    allocate(pbp_gref(1))
    pbp_gref(1) = subset(1)
    allocate(pbp_outref(1))
    pbp_outref(1) = 1

    allocate(pbp_pref(npbp,nlu))
    allocate(pbp_larea(npbp,nlu))
    pbp_pref(1,:) = sibvs(1)%pftref(:)
    pbp_larea(1,:) = sibvs(1)%larea(:)

    print('(2a)'),'     Site Name: ',subsitename(1)
    print('(a,2F8.2)'),'     Site Lon/Lat: ',sublonsib(1),sublatsib(1)

    j=0
    do i=1,nlu
       if (subpref(1,i) .gt. 0) then
           j=j+1
           if (pft_num(subpref(1,i)) .eq. 0) then
              print*,''
              print('(a,i4)'),'Unknown PFT: ',subpref(1,i)
              print('(a)'),'Stopping.'
              print*,''
              stop
           endif
       endif
    enddo
    print('(a,10a5)'),'     Site PFTs: ',pft_name(pft_num(subpref(1,1:j)))
else
    !...Regional/Global Simulations

    !...calculate subcount for parallelization
    subgridcount=0
    do i=1,nsib
       if ((lonsib(i) >= minlon) .and. (lonsib(i) <= maxlon) .AND. &
           (latsib(i) >= minlat) .and. (latsib(i) <= maxlat)) then
           subgridcount = subgridcount + 1
       endif
    enddo

    allocate(subgridref(subgridcount))
    subcount=1
    do i=1,nsib
       if ((lonsib(i) >= minlon) .and. (lonsib(i) <= maxlon) .AND. &
           (latsib(i) >= minlat) .and. (latsib(i) <= maxlat)) then
           subgridref(subcount) = i
           subcount=subcount+1 
       endif
    enddo

    olength = nint(subgridcount / dble(nchunks))
    start_count=(rank-1)*olength+1
    stop_count=rank*olength
    if (rank .eq. nchunks) stop_count=subgridcount

    !...count number of landpoints in subdomain
    subcount=stop_count-start_count+1
    if (subcount .lt. 1) then
        print*,'No SiB points to simulate.  Stopping.'
        stop
    endif

    allocate(index_ref(subcount))
    index_ref(:) = subgridref(start_count:stop_count)

    print('(a,i5)'),'   SiB points simulated=', subcount

    !...allocate subset and assign values
    allocate( subset(subcount) )
    subset(:) = index_ref( 1 : subcount )

    allocate(sublatsib(subcount), sublonsib(subcount))
    allocate(sublarea(subcount,nlu), subpref(subcount,nlu))

    do i=1, subcount
       sublatsib(i) = latsib(subset(i))
       sublonsib(i) = lonsib(subset(i))
       
       sublarea(i,:) = sibvs(subset(i))%larea(:)
       subpref(i,:) = sibvs(subset(i))%pftref(:)

    enddo

    allocate(subsitename(subcount))
    do i=1,subcount
        subsitename(i) = sitenamesib(subset(i))
    enddo

    if (subcount < 2) then
        print*,'      SiB Point Number=',subset(1)
        print*,'      SiB Lon/Lat=',sublonsib(1),sublatsib(1)    
    endif

    !...find pbp indices and remove duplicates
    if (npbp .GT. 0) then
       allocate(temp_pbp(npbp),temp_sibpt(npbp))
       temp_pbp = 0
       temp_sibpt = 0
       pbp_index = 0
       sibpt_index = 0
       do i = 1, npbp
          ! find closest point
          pbp_diff = minval( abs(pbp_lonlat(2,i) - latsib) + &
                             abs(pbp_lonlat(1,i) - lonsib) )

          do j=1,subcount
             if (pbp_diff == (abs (pbp_lonlat(2,i) - sublatsib(j)) + &
                              abs (pbp_lonlat(1,i) - sublonsib(j)))) then
                 pbp_index = j
                !print*,'PBP Point lon/lat: ',pbp_lonlat(1,i),pbp_lonlat(2,i),i
                !print*,' Sim pt/lon/lat: ',j,lonsib(j),latsib(j)
                exit
             endif
           enddo

           do j=1,nsib
              if (pbp_diff == (abs(pbp_lonlat(2,i) - latsib(j)) + &
                               abs(pbp_lonlat(1,i) - lonsib(j)))) then
                  sibpt_index = j
                  exit
              endif
           enddo

           !if (pbp_index == 0) then
           !    print*, 'Point ', pbp_lonlat(1,i), pbp_lonlat(2,i)
           !    print*, '  is not inside the grid, ignoring.'
           !endif
           temp_pbp(i) = pbp_index
           temp_sibpt(i) = sibpt_index
       enddo  !npbp

       !...remove duplicates and values outside of domain
       do i = 1, npbp-1
          if (temp_pbp(i) /= 0) then
             do j = i+1, npbp
                 if ( temp_pbp(i) == temp_pbp(j) ) then
                     ! duplicate, set second instance to zero
                     temp_pbp(j) = 0
                     temp_sibpt(j) = 0
                     !print *, 'Point', pbp_lonlat(1,j), pbp_lonlat(2,j),  &
                     !    ' is a duplicate of point', pbp_lonlat(1,i), pbp_lonlat(2,i)
                     !print *, 'This point has been removed'
                  endif
              enddo
          endif
       enddo !npbp

       !...count remaining number of pbp's, allocate new array and copy subset index
       j = 0
       do i = 1, npbp
          if ( temp_pbp(i) /= 0 ) j = j + 1
       enddo
    elseif (npbp == -1) then
       !...include all points in pbp's
       j = subcount
       npbp = subcount
       allocate(temp_pbp(subcount),temp_sibpt(subcount))
       temp_sibpt(:) = subset(:)
       do i = 1, npbp
          temp_pbp(i) = i
       enddo

       pbp_lonlat(1,1:subcount) = sublonsib(:)
       pbp_lonlat(2,1:subcount) = sublatsib(:)
    else
        j=0
    endif  !specific or all pbps

    if ((pbp_dtsib .ne. 0) .and. &
         (j > 0)) then
       print('(a,i5)'),'   SiB PBPs simulated= ', j

       !...there is at least one pbp in the subdomain
       allocate( pbp_outref(subcount) )
       allocate( pbp_sibpt(j) )
       allocate( pbp_gref(j) )
       allocate( pbp_lat(j) )
       allocate( pbp_lon(j) )
       allocate( pbp_larea(j,nlu) )
       allocate( pbp_pref(j,nlu) )
       allocate( pbp_sitename(j))

       pbp_outref(:) = izero
       pbp_sibpt(:) = izero
       pbp_gref(:) = izero
       pbp_lat(:) = rzero
       pbp_lon(:) = rzero
       pbp_pref(:,:) = izero
       pbp_larea(:,:) = rzero
       pbp_sitename(:) = ''

       j = 0
       do i = 1, npbp
          if ( temp_pbp(i) /= 0 ) then
              j = j + 1
              pbp_lat(j) = pbp_lonlat(2,i)
              pbp_lon(j) = pbp_lonlat(1,i)

              gref = temp_sibpt(i)
              pbp_sibpt(j) = gref
              pbp_outref(temp_pbp(i)) = j
              pbp_gref(j) = temp_pbp(i)
              pbp_sitename(j) = sitenamesib(gref)

              if ((pbp_dtsib .ne. 0) .and. (npbp .lt. 100)) then
                  print'(a,2f10.2,a,a6)','       PBP Lon/Lat: ', &
                        pbp_lon(j),pbp_lat(j),'   ',trim(pbp_sitename(j))
              endif
                  
              do l=1,sibvs(gref)%gnlu
                 pbp_larea(j,l) = sibvs(gref)%larea(l)
                 pbp_pref(j,l) = sibvs(gref)%pftref(l)
              enddo  !l=1,g_nlu
          endif !temp_pbp(i) /= 0 
       enddo !i=1,npbp
       npbp = j

       deallocate(temp_pbp)
       deallocate(index_ref)
       deallocate(subgridref)

    else !dtsibpbp == 0 and/or j == 0
          pbp_savegf = .false.
          pbp_saveluf = .false.
          npbp = 0
    endif

endif  !if not a single point run

if (allocated(pbp_lonlat)) deallocate(pbp_lonlat)

end subroutine grid_init

