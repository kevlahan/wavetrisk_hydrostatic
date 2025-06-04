program vertices2wavetrisk
  ! Converts generic icosahedral refined triangle vertex data in .xyz format to WAVETRISK format
  ! gen is number of subdivisions of icosahedron (icosahedron is gen = 0)
  use geom_mod
  implicit none
  ! Command line inputs
  integer            :: lev_beg  ! first level to convert
  integer            :: lev_end  ! last  level to convert
  character(len=999) :: data_pre ! prefix for dataset including path: path/data_pre

  integer                                :: i, count, iostat, jmin, l, num_vertices, stat
  real(dp)                               :: dmax, dmin, dx
  type(coord)                            :: vertHR
  type(coord), dimension(:), allocatable :: distance, vert_data
  character(len=999)                     :: arg, datafile, grid_file
  
  radius = 1.0_dp

  ! Read input data
  count = command_argument_count ()
  if (count == 0) then
     write (6,'(/,a,/)') "Usage:"
     write (6,'(a,/)') "./vertices2wavetrisk data_pre lev_beg lev_end"
     write (6,'(a,/)') "Example to process levels 5 to 7 grids with prefix HR95HK in directory ~/icos_grids:"
     write (6,'(a)') "./vertices2wavetrisk ~/icos_grids/HR95JT 5 7"
     stop
  end if
  call get_command_argument (1, arg);  data_pre = arg
  call get_command_argument (2, arg) ; read (arg, *, iostat=stat) lev_beg
  call get_command_argument (3, arg) ; read (arg, *, iostat=stat) lev_end
  write (6,'(a,a,1x,2(i1,1x))') "Input data = ", trim (data_pre), lev_beg, lev_end

  do l = lev_beg, lev_end
     dx = dx_avg (l)
     
     ! Read vertex data to be converted
     write (datafile, '(a,a,i3.3,a)')  trim(data_pre), "_", l, ".xyz"
     open (unit = 10, file = trim (datafile), status = "old", form = "formatted", action = "read")

     read (10,*) num_vertices
     allocate (vert_data(1:num_vertices), distance(1:num_vertices))

     do i = 1, num_vertices
        read (10,*) vert_data(i)
     end do
     close (10)

     ! File for converted data
     write (datafile, '(a,a,i3.3)')  trim (data_pre), "_WT_", l
     open (unit = 20, file = trim(datafile), status = "replace", form = "formatted", action = "write")

     ! Heikes & Randall data in WAVETRISK form (does not include poles)
     write (grid_file, '(a,i3.3)')  "grids/HR95JT_WT_", l
     open (unit = 10, file = trim (grid_file), status = "old", form = "formatted", action = "read")

     i = 0
     dmax = 0d0
     do
        read (10, *, iostat=iostat) vertHR ! read each Heikes & Randall vertex coordinate
        if (iostat /= 0) exit              ! end of file
        i = i + 1

        call min_dist (vertHR, vert_data, dmin, jmin) ! find closest vertex to Heikes & Randall vertex
        
        if (dmin > dmax) dmax = dmin

        write (20,'(3(es24.17e2,1x))') vert_data(jmin)
     end do
     write (6,'(a,i1,a,es8.2)') "Max diff in vertices at level ", l, " is ", dmax
     if (dmax > dx) write (6,'(a,es10.4)') "WARNING: minimum distance > dx_avg = ", dx
     
     close (10); close (20); deallocate (distance, vert_data)
  end do
end program vertices2wavetrisk
