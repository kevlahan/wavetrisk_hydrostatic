!to compile: execute
!                   gfortran trisk2vtksinglelayer.f90 -o trisk2vtksinglelayer
!to run:
!                   ./trisk2vtksinglelayer fort.1 0 0 5 7 vtkoutput

program trisk2vtk
  ! Converts standard ASCII trisk data files from wavetrisk code to BINARY .vtk format for paraview
  !
  ! Usage: see below

  implicit none
  integer :: icell, icoord, ivert, itime, n_cells, n_vertices, tstart, tend
  integer :: stat
  integer :: j, jmin, jmax, n_cells_old
  integer, parameter :: iunit=10
  integer, parameter :: vtk_type = 7 ! vtk polygonal data
  integer, dimension(:), allocatable :: tmp_mask, mask, tmp_level, level

  character(2)   :: j_lev
  character(3)   :: s_time
  character(12)  :: str1, str2
  character(5+6+2) :: filename_in
  character(32)  :: file_vtk, arg, command, filename_out, file_base
  character(1), parameter :: lf=char(10) ! line feed character
  
  real(4), dimension (:,:),   allocatable :: tmp_outv, outv
  real(4), dimension (:,:,:), allocatable :: tmp_vertices, vertices

  ! Get input parameters
  call get_command_argument(1, arg); file_base = trim(arg)
  call get_command_argument(2, arg); read (arg,'(I12)') tstart
  call get_command_argument(3, arg); read (arg,'(I12)') tend
  call get_command_argument(4, arg); read (arg,'(I12)') jmin
  call get_command_argument(5, arg); read (arg,'(I12)') jmax
  call get_command_argument(6, arg); file_vtk =  trim(arg)

  if (file_base .eq. " ") then
     write(*,*) " "
     write(*,*) "Usage: trisk2vtk file_base tbegin tend jmin jmax file_vtk"
     write(*,*) " "
     write(*,*) "file_base = base name for files"
     write(*,*) "tstart    = number of first file to read"
     write(*,*) "tend      = number of last file to read"
     write(*,*) "jmin      = minimum scale to save"
     write(*,*) "jmax      = maximum scale to save"
     write(*,*) "file_vtk  = file for vtk data (without extension)"
     stop
  else
     write(*,*) "Were read in:"
     write(*,*) " "
     write(*,*) "file_base = base name for files  ", file_base
     write(*,*) "tstart    = number of first file to read  ", tstart
     write(*,*) "tend      = number of last file to read  ", tend
     write(*,*) "jmin      = minimum scale to save  ", jmin
     write(*,*) "jmax      = maximum scale to save  ", jmax
     write(*,*) "file_vtk  = file for vtk data (without extension)  ", file_vtk
  end if

  do itime = tstart, tend
     PRINT *, '-------'
     write (s_time,'(I3.3)') itime

     ! Loop through scales
     n_cells = 0
     do j = jmin, jmax
        write (j_lev,'(I2.2)') j
        !PRINT *, 's_time j_lev', s_time, j_lev
        filename_in = trim(file_base) // s_time // j_lev 
        write(*,'("Reading from file ",A)') filename_in

        ! First pass to determine number of cells
        open (unit=iunit, file=filename_in, form="formatted")
        
        n_cells_old = n_cells
        stat = 0
        do while (stat>=0)
           read (iunit, IOSTAT=stat,FMT=*)
           if (stat==0) then
              n_cells = n_cells+1
           elseif (stat>0) then
              write(*,*) "Error reading file"
              stop
           end if
        end do
        close (iunit)
        
        write (*,'("Number of cells found at level ", i2, " = ",i12)') j, n_cells-n_cells_old

        n_vertices=6

        if (j>jmin) then ! Need to increase size of arrays and keep previous values
           allocate(tmp_vertices(1:n_cells,1:n_vertices,1:3))
           tmp_vertices(1:size(vertices,1),:,:) = vertices
           deallocate(vertices)
           allocate(vertices(1:n_cells,1:n_vertices,1:3))
           vertices = tmp_vertices
           deallocate(tmp_vertices)
           
           allocate(tmp_outv(1:n_cells,1:4))
           tmp_outv(1:size(outv,1),:) = outv
           deallocate(outv)
           allocate(outv(1:n_cells,1:4))
           outv = tmp_outv
           deallocate(tmp_outv)
           
           allocate(tmp_mask(1:n_cells))
           tmp_mask(1:size(mask)) = mask
           deallocate(mask)
           allocate(mask(1:n_cells))
           mask = tmp_mask
           deallocate(tmp_mask)
           
           allocate(tmp_level(1:n_cells))
           tmp_level(1:size(level)) = level
           deallocate(level)
           allocate(level(1:n_cells))
           level = tmp_level
           deallocate(tmp_level)  
        else
           allocate (vertices(1:n_cells,1:n_vertices,1:3), mask(1:n_cells), level(1:n_cells), outv(1:n_cells,1:4))
        end if
        
        ! Second pass to actually read in data
        open (unit=iunit, file=filename_in, form="formatted")
           do icell = n_cells_old+1, n_cells
              read (iunit, fmt='(18(E14.5E2, 1X), 4(E14.5E2, 1X), I3, 1X, I3)') &
                   ((vertices(icell,ivert,icoord),icoord=1,3),ivert=1,n_vertices), &
                   outv(icell,1:4), mask(icell), level(icell)
           end do
        close(iunit)   
     end do
          
     ! Write out data in old vtk format
     write (s_time, '(I3.3)')  itime
     filename_out = trim(file_vtk) // s_time // ".vtk"

     write (*,'("File ", A," opened to write vtk output field")') filename_out
     
     open (unit=iunit, file=filename_out, form="unformatted", access='stream',convert='BIG_ENDIAN')
     write(iunit) '# vtk DataFile Version 3.0'//lf
     write(iunit) 'vtk output'//lf                
     write(iunit) 'BINARY'//lf                    
     write(iunit) 'DATASET UNSTRUCTURED_GRID'//lf  

     write(str1(1:12),'(i12)') n_cells*n_vertices
     write(iunit) 'POINTS '//str1//' float'//lf
     
     ! Write out vertices
     do icell = 1, n_cells
        do ivert = 1, n_vertices
           write(iunit) (vertices(icell,ivert,icoord),icoord=1,3)
        end do
     end do
     write(iunit) lf
     
     ! Write out cells
     write(str1(1:12),'(i12)') n_cells
     write(str2(1:12),'(i12)') n_cells*(1+n_vertices)
     write(iunit) 'CELLS '//str1//str2//lf
     do icell = 1, n_cells
        write (iunit) n_vertices, ((icell-1)*n_vertices+ivert-1, ivert=1,n_vertices)
     end do
     write(iunit) lf
     
     ! Write out type of each cell (all polygons!)
     write(iunit) 'CELL_TYPES '//str1//lf
     do icell = 1, n_cells
        write (iunit) vtk_type
     end do
     write(iunit) lf

     ! Write out cell data for each cell
     write(iunit) 'CELL_DATA '//str1//lf
        write(iunit) 'SCALARS temp float'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) outv(icell,1)
        end do
        write(iunit) lf
        
        write(iunit) 'SCALARS column float'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) outv(icell,2)
        end do
        write(iunit) lf
        
        write(iunit) 'SCALARS mass float'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) outv(icell,3)
        end do
        write(iunit) lf

        write(iunit) 'SCALARS mass float'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) outv(icell,4)
        end do
        write(iunit) lf
        
        write(iunit) 'SCALARS mask int'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) mask(icell)
        end do
        write(iunit) lf
        
        write(iunit) 'SCALARS level int'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) level(icell)
        end do
     
     close(iunit)

    deallocate (vertices, mask, level, outv)
  end do
end program trisk2vtk
 
