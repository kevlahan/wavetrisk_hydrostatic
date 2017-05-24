program trisk2vtk
  ! Converts standard ASCII trisk data files from wavetrisk code to BINARY .vtk format for paraview
  !
  ! Usage: trisk2vtk filetype tstart tend jmin jmax file_vtk
  !
  ! filetype = primal (hexagons) or dual (triangles)
  ! tstart   = first file to read
  ! tend     = last file to read
  ! jmin     = minimum scale to save
  ! jmax     = maximum scale to save
  ! file_vtk = file for vtk data

  implicit none
  integer :: fid, i, icell, icoord, ipts, ivert, itime, n_cells, n_vertices, tstart, tend, u
  integer :: stat
  integer :: j, jmin, jmax, n_cells_old
  integer, parameter :: iunit=10
  integer, parameter :: vtk_type = 7 ! vtk polygonal data
  integer, dimension(:), allocatable :: tmp_mask, mask, tmp_level, level

  character(2)   :: j_lev
  character(3)   :: s_time
  character(12)  :: str1, str2
  character(6)   :: filetype
  character(32)  :: file_vtk, arg, command, filename_in, filename_out, file_base
  character(1), parameter :: lf=char(10) ! line feed character
  
  real(4), dimension (:,:),   allocatable :: tmp_outv, outv
  real(4), dimension (:,:,:), allocatable :: tmp_vertices, vertices

  real(4), parameter :: Hdim = 3.344175893265152e+03, g = 9.80665
  real(4), parameter :: Ldim = 6371e3
  real(4), parameter :: Udim = sqrt(Hdim*g)
  real(4), parameter :: Tdim = Ldim/Udim

  logical :: file_exists

  ! Get input parameters
  call get_command_argument(1, arg); file_base = arg
  call get_command_argument(2, arg); filetype = arg
  call get_command_argument(3, arg); read (arg,'(I12)') tstart
  call get_command_argument(4, arg); read (arg,'(I12)') tend
  call get_command_argument(5, arg); file_vtk =  arg

  if (filetype .eq. " ") then
     write(*,*) " "
     write(*,*) "Usage: trisk2vtk file_base filetype tbegin tend jmin jmax file_vtk"
     write(*,*) " "
     write(*,*) "file_base = base name for files"
     write(*,*) "filetype  = primal (hexagons) or dual (triangles)"
     write(*,*) "tstart    = number of first file to read"
     write(*,*) "tend      = number of last file to read"
     write(*,*) "file_vtk  = file for vtk data (without extension)"
     stop
  else
     write(*,'("file_base = ", A)') file_base
     write(*,'("filetype  = ", A)') filetype
     write(*,'("tstart    = ", I12)') tstart
     write(*,'("tend      = ", I12)') tend
     write(*,'("file_vtk  = ", A)') file_vtk

     if (trim(filetype) .eq. "primal") then
        n_vertices = 6 ! Hexagonal cells (primal grid)
     elseif (trim(filetype) .eq. "dual") then
        n_vertices = 3 ! Triangular cells (dual grid)
     end if
  end if

  do itime = tstart, tend
     ! Uncompress base file
     write (s_time,'(I3.3)') itime
     command = 'tar xjf ' // trim(file_base) // s_time // '.tbz'
     write(*,*) command
     CALL system(command)

     ! Delete un-needed file
     command = '\rm ' // trim(file_base) // s_time // '00'
     CALL system(command)

     ! Find number of scales
     j = 1 ; jmin = 0
     do
        ! Check if scale j file exists
        write (j_lev,'(I2.2)') j
        filename_in = trim(file_base) // s_time // j_lev
        inquire(FILE=trim(filename_in), EXIST=file_exists)
        
        if (.not.(file_exists)) then
           if (jmin/=0) then
              jmax = j-1
              exit
           end if
        elseif (jmin==0) then
           jmin = j
        end if
        j = j+1
     end do
     write(*,'("Minimum scale = ", i2, " maximum scale = ", i2)') jmin, jmax
     
     ! Loop through all scales
     n_cells = 0
     j = 1
     do j = jmin, jmax
        write (j_lev,'(I2.2)') j
        filename_in = trim(file_base) // s_time // j_lev

        write(*,'("Reading from file ",A)') filename_in
        
        ! First pass to determine number of cells
        open (unit=iunit, file=trim(filename_in), form="formatted")
        
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
           if(allocated(vertices)) deallocate(vertices)
           if(allocated(mask)) deallocate(mask)
           if(allocated(level)) deallocate(level)
           if(allocated(outv)) deallocate(outv)
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

        ! Delete uncompressed file
        command = '\rm ' // trim(filename_in)
        CALL system(command)
     end do
          
     ! Write out data in old vtk format
     write (s_time, '(I3.3)')  itime
     filename_out = trim(file_vtk) // '_' // s_time // '.vtk'

     write (*,'("File ", A," opened to write vtk output field")') trim(filename_out)
     
     open (unit=iunit, file=filename_out, form="unformatted", access='stream',convert='BIG_ENDIAN')
     write(iunit) '# vtk DataFile Version 3.0'//lf
     write(iunit) 'vtk output'//lf                
     write(iunit) 'BINARY'//lf                    
     write(iunit) 'DATASET UNSTRUCTURED_GRID'//lf  

     write(str1(1:12),'(i12)') n_cells*n_vertices
     write(iunit) 'POINTS ' // trim(str1) // ' float'// lf
     
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
           write (iunit) outv(icell,1) * Hdim
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
           write (iunit) outv(icell,3) * Hdim
        end do
        write(iunit) lf

        write(iunit) 'SCALARS KE float'//lf
        write(iunit) 'LOOKUP_TABLE default'//lf
        do icell = 1, n_cells
           write (iunit) outv(icell,4) * Udim**2
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
           write (iunit) outv(icell,1) / Tdim
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
  end do
end program trisk2vtk
 
