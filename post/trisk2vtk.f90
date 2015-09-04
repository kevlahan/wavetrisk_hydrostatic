program trisk2vtk
  ! Converts standard ASCII data files from wavetrisk code to BINARY .vtk format for paraview
  !
  ! Usage: trisk2vtk filetype iwrite jmin jmax file_vtk
  !
  ! filetype = primal (hexagons) or dual (triangles)
  ! iwrite   = file to read
  ! jmin     = minimum scale to save
  ! jmax     = maximum scale to save
  ! file_vtk = file for vtk data

  implicit none
  integer :: fid, i, icell, icoord, ipts, ivert, iwrite, n_cells, n_vertices, u
  integer :: stat
  integer :: j, jmin, jmax, n_cells_old
  integer, parameter :: iunit=10
  integer, parameter :: vtk_type = 7 ! vtk polygonal data
  integer, dimension(:), allocatable :: tmp_mask, mask, tmp_level, level
  
  character(12)  :: str1, str2
  character(6)   :: filetype
  character(5+6) :: filename_in
  character(32)  :: file_vtk, arg
  character(1), parameter :: lf=char(10) ! line feed character
  
  real(4), dimension (:,:),   allocatable :: tmp_outv, outv
  real(4), dimension (:,:,:), allocatable :: tmp_vertices, vertices

  ! Get input parameters
  call get_command_argument(1, arg); filetype = trim(arg)
  call get_command_argument(2, arg); read (arg,'(I12)') iwrite
  call get_command_argument(3, arg); read (arg,'(I12)') jmin
  call get_command_argument(4, arg); read (arg,'(I12)') jmax
  call get_command_argument(5, arg); file_vtk =  trim(arg)//'.vtk'
  
  n_cells = 0
  do j = jmin, jmax
     if (filetype .eq. "primal") then
        n_vertices = 6 ! Hexagonal cells (primal grid)
        u = 100000+100*iwrite
     elseif (filetype .eq. "dual") then
        n_vertices = 3 ! Triangular cells (dual grid)
        u = 200000+100*iwrite
     elseif (filetype .eq." ") then
        write(*,*) " "
        write(*,*) "Usage: trisk2vtk filetype iwrite jmin jmax file_vtk"
        write(*,*) " "
        write(*,*) "filetype = primal (hexagons) or dual (triangles)"
        write(*,*) "iwrite   = file to read"
        write(*,*) "jmin     = minimum scale to save"
        write(*,*) "jmax     = maximum scale to save"
        write(*,*) "file_vtk = file for vtk data"
        stop
     end if
     fid = u+j
     write (filename_in, '(A,I6)')  "fort.", fid

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
     if (filetype.eq."primal") then
        ! outv(:,1) is topography, outv(:,2) is penalization mask, outv(:,3) is height, outv(:,4) is kinetic energy
        do icell = n_cells_old+1, n_cells
           read (iunit, fmt='(18(E14.5E2, 1X), 4(E14.5E2, 1X), I3, 1X, I3)') &
                ((vertices(icell,ivert,icoord),icoord=1,3),ivert=1,n_vertices), &
                outv(icell,1:4), mask(icell), level(icell)
        end do
     elseif (filetype.eq."dual") then
        ! outv(:,1) is relative vorticity and outv(:,2) is the dual chi mask
        do icell = n_cells_old+1, n_cells
           read (iunit, fmt='(9(E14.5E2, 1X), 2(E14.5E2, 1X), I3)') &
                ((vertices(icell,ivert,icoord),icoord=1,3),ivert=1,n_vertices), &
                outv(icell,1), outv(icell,2), level(icell)
        end do
     end if
     close(iunit)   
  end do
     
  ! Write out data in old vtk format
  write (*,'("File ", A," opened to write vtk output field")') trim(file_vtk)

  open (unit=iunit, file=file_vtk, form="unformatted", access='stream',convert='BIG_ENDIAN')
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
  if (filetype.eq."primal") then
     write(iunit) 'SCALARS topography float'//lf
     write(iunit) 'LOOKUP_TABLE default'//lf
     do icell = 1, n_cells
        write (iunit) outv(icell,1)
     end do
     write(iunit) lf

     write(iunit) 'SCALARS penalization float'//lf
     write(iunit) 'LOOKUP_TABLE default'//lf
     do icell = 1, n_cells
        write (iunit) outv(icell,2)
     end do
     write(iunit) lf

     write(iunit) 'SCALARS height float'//lf
     write(iunit) 'LOOKUP_TABLE default'//lf
     do icell = 1, n_cells
        write (iunit) outv(icell,3)
     end do
     write(iunit) lf

     write(iunit) 'SCALARS KE float'//lf
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
  elseif (filetype.eq."dual") then
     write(iunit) 'SCALARS relvort float'//lf
     write(iunit) 'LOOKUP_TABLE default'//lf
     do icell = 1, n_cells
        write (iunit) outv(icell,1)
     end do
     write(iunit) lf

     write(iunit) 'SCALARS penal float'//lf
     write(iunit) 'LOOKUP_TABLE default'//lf
     do icell = 1, n_cells
        write (iunit) outv(icell,2)
     end do
     write(iunit) lf

     write(iunit) 'SCALARS level int'//lf
     write(iunit) 'LOOKUP_TABLE default'//lf
     do icell = 1, n_cells
        write (iunit) level(icell)
     end do
  end if
  
  close(iunit)
end program trisk2vtk
 
