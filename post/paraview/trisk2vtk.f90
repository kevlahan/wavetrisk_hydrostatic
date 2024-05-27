program trisk2vtk
  ! Converts trisk data files from wavetrisk code to BINARY .vtk format for paraview
  !
  implicit none
  integer :: fid, i, icell, icoord, ipts, itime, ivert, k, ncells, ncells_adapt, n_vertices, tstart, tend
  integer :: stat
  integer :: j, jmin, jmax, ncells_old, nvar_out, zmax, zmin
  integer, parameter :: iunit=10
  integer, parameter :: vtk_type = 7 ! vtk polygonal data
  integer, dimension(:), allocatable :: tmp, mask, level

  character(2)   :: j_lev
  character(3)   :: zlev
  character(4)   :: s_time
  character(12)  :: str1, str2
  character(6)   :: grid_type
  character(255) :: fmt, nvar, nvert
  character(255) :: arg, bash_cmd, command, filename_in, filename_out, file_base, sim_type
  character(1), parameter :: lf=char(10) ! line feed character

  real(8), dimension (:,:),   allocatable :: tmp_outv, outv
  real(8), dimension (:,:,:), allocatable :: tmp_vertices, vertices

  real(8), parameter :: Hdim = 3.344175893265152e+03, g = 9.80665
  real(8), parameter :: Ldim = 6371e3
  real(8), parameter :: Udim = sqrt(Hdim*g)
  real(8), parameter :: Tdim = Ldim/Udim

  logical :: compressed_file_exists, file_exists

  ! Get input parameters
  call get_command_argument (1, arg); file_base = trim(arg)
  call get_command_argument (2, arg); sim_type  = trim(arg)
  call get_command_argument (3, arg); grid_type = trim(arg)
  call get_command_argument (4, arg); read (arg,'(I12)') tstart
  call get_command_argument (5, arg); read (arg,'(I12)') tend
  call get_command_argument (6, arg); read (arg,'(I12)') jmin
  call get_command_argument (7, arg); read (arg,'(I12)') jmax
  call get_command_argument (8, arg); read (arg,'(I12)') zmin
  call get_command_argument (9, arg); read (arg,'(I12)') zmax
  
  if (grid_type == " ") then
     write (6,'(/,a)') "Converts trisk data files from wavetrisk code to BINARY .vtk format for paraview"
     write (6,'(/,a)')"Usage: trisk2vtk file_base sim_type grid_type tstart tend jmin jmax zmin zmax"
     write (6,'(/,a)')"Example: ./trisk2vtk drake sea hex 0 1 5 7 3 4"
     write (6,'(/,a)')"Reads files drake_003.1.0001.tgz, drake_004.1.0001.tgz, drake_003.1.0002.tgz, drake_004.1.0002.tgz"
     write (6,'(a,/)')"for scale 5 to 7."
     write (6,'(a)') "file_base = name of file"
     write (6,'(a)') "sim_type  = atm (atmosphere, compressible), sea (sea, incompressible), or 2layer (sea, 2 layer)"
     write (6,'(a)') "grid_type = hex (overlapping hexagons on each scale) or tri (non-overlapping triangles of adaptive grid)"
     write (6,'(a)') "tstart    = first file to read"
     write (6,'(a)') "tend      = last file to read"
     write (6,'(a)') "jmin      = minimum scale"
     write (6,'(a)') "jmax      = maximum scale"
     write (6,'(a)') "zmin      = min vertical layer"
     write (6,'(a,/)') "zmax      = max vertical layer"
     stop 
  else
     write (6,'(/,"file_base = ", a)') trim(file_base)
     write (6,'("sim_type  = ", a)') trim(sim_type)
     write (6,'("grid_type = ", a)') trim(grid_type)
     write (6,'("tstart    = ", i4.4)') tstart
     write (6,'("tend      = ", i4.4)') tend
     write (6,'("jmin      = ", i2.2)') jmin
     write (6,'("jmax      = ", i2.2)') jmax
     write (6,'("zmin      = ", i3.3)') zmin
     write (6,'("zmax      = ", i3.3)') zmax

     if (trim(grid_type) == "hex") then     ! hexagonal cells (primal grid)
        n_vertices = 6
     elseif (trim(grid_type) == "tri") then ! triangular cells (dual grid)
        n_vertices = 3 
     end if
     write(nvert,*) 3*n_vertices ! number of vertex coordinates

     if (trim(sim_type) == "2layer" .and. trim(grid_type) == "hex") then
        nvar_out = 11
     else
        nvar_out = 9
     end if
     write (nvar,*) nvar_out ! number of saved variables
  end if

  do itime = tstart, tend
     write (s_time,'(i4.4)') itime

     do k = zmin, zmax
        if (trim(grid_type) == "hex") then
           write (zlev,'(i3.3)') k
           filename_in = trim(file_base)//'_'//zlev//'.1'//s_time
        elseif (trim(grid_type) == "tri") then
           write (zlev,'(i3.3)') k
           filename_in = trim(file_base)//'_'//zlev//'.2'//s_time
        end if

        ! Check if files are compressed
        inquire (FILE=trim(filename_in)//'.tgz', EXIST=compressed_file_exists)
        if (compressed_file_exists) then
           command = 'tar xzf '//trim(filename_in)//'.tgz'
           write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
           call system (trim(bash_cmd))

           ! Delete un-needed file
           if (trim(grid_type) == "hex" .or. trim(grid_type) == "tri" .or. trim(grid_type) == "2layer") then
              write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
              call system (trim(bash_cmd))
           end if
        end if
        write (6,'("Reading from file ",a)') trim(filename_in)
        
        ! Loop through all scales
        ncells = 0
        j = 1
        do j = jmin, jmax
           ! First pass to determine number of cells
           write (j_lev,'(i2.2)') j
           open (iunit, file=trim(filename_in)//j_lev, form="unformatted")

           ncells_old = ncells
           stat = 0
           do while (stat>=0)
              read (iunit, IOSTAT=stat)
              if (stat==0) then
                 ncells = ncells+1
              elseif (stat>0) then
                 write (6,'(a)') "Error reading file"
                 stop
              end if
           end do
           close (iunit)

           write (6,'("  Number of cells at level ", i2.2, " = ", i10)') j, ncells-ncells_old

           if (j > jmin) then ! need to increase size of arrays and keep previous values
              allocate (tmp_vertices(1:ncells,1:n_vertices,1:3))
              tmp_vertices(1:size(vertices,1),:,:) = vertices
              call move_alloc (tmp_vertices, vertices)

              allocate (tmp_outv(1:ncells,1:nvar_out))
              tmp_outv(1:size(outv,1),:) = outv
              call move_alloc (tmp_outv, outv)

              allocate (tmp(1:ncells))
              tmp(1:size(mask)) = mask
              call move_alloc (tmp, mask)

              allocate (tmp(1:ncells))
              tmp(1:size(level)) = level
              call move_alloc (tmp, level)
           else
              if (allocated(vertices)) deallocate (vertices)
              if (allocated(mask))     deallocate (mask)
              if (allocated(level))    deallocate (level)
              if (allocated(level))    deallocate (level)
              if (allocated(outv))     deallocate (outv)
              allocate (vertices(1:ncells,1:n_vertices,1:3), level(1:ncells), mask(1:ncells), outv(1:ncells,1:nvar_out))
           end if

           ! Second pass to read in data
           open (iunit, file=trim(filename_in)//j_lev, form="unformatted")
           if (grid_type == "hex") then
              do icell = ncells_old+1, ncells
                 read (iunit) ((vertices(icell,ivert,icoord),icoord=1,3), ivert = 1, n_vertices), &
                      outv(icell,1:nvar_out), mask(icell), level(icell)
              end do
           elseif (grid_type == "tri") then
              do icell = ncells_old+1, ncells
                 read (iunit) ((vertices(icell,ivert,icoord),icoord=1,3), ivert = 1, n_vertices), &
                      outv(icell,1:nvar_out), mask(icell), level(icell)
              end do
           elseif (grid_type == "2layer") then
              do icell = ncells_old+1, ncells
                 read (iunit) ((vertices(icell,ivert,icoord),icoord=1,3), ivert = 1, n_vertices), &
                      outv(icell,1:nvar_out), mask(icell), level(icell)
              end do
           end if
           close(iunit)
        end do

        ! Delete uncompressed files
        if (compressed_file_exists) then
           command = '\rm ' // trim(filename_in)//"0*"
           write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
           call system (trim(bash_cmd))
        end if

        ! Write out data in old vtk format
        write (6,'("Saving output to file ", a,/)') trim(filename_in)//'.vtk'

        open (unit=iunit, file=trim(filename_in)//'.vtk', form="unformatted", access='stream',convert='BIG_ENDIAN')
        write(iunit) '# vtk DataFile Version 3.0'//lf
        write(iunit) 'vtk output'//lf                
        write(iunit) 'BINARY'//lf                    
        write(iunit) 'DATASET UNSTRUCTURED_GRID'//lf  

        write(str1(1:12),'(i12)') ncells * n_vertices
        write(iunit) 'POINTS ' // trim(str1) // ' double'// lf

        ! Write out vertices
        do icell = 1, ncells
           do ivert = 1, n_vertices
              write(iunit) (vertices(icell,ivert,icoord),icoord=1,3)
           end do
        end do
        write(iunit) lf

        ! Write out cells
        write(str1(1:12),'(i12)') ncells
        write(str2(1:12),'(i12)') ncells*(1+n_vertices)
        write(iunit) 'CELLS '//str1//str2//lf
        do icell = 1, ncells
           write (iunit) n_vertices, ((icell-1)*n_vertices+ivert-1, ivert=1,n_vertices)
        end do
        write(iunit) lf

        ! Write out type of each cell (all polygons!)
        write(iunit) 'CELL_TYPES '//str1//lf
        do icell = 1, ncells
           write (iunit) vtk_type
        end do
        write(iunit) lf

        ! Write out cell data for each cell                        
        write(iunit) 'CELL_DATA '//str1//lf
        
        if (.not. grid_type == "2layer") then ! usual (hex or tri)
           if (sim_type == "atm") then
              write(iunit) 'SCALARS temperature double'//lf
           else
              write(iunit) 'SCALARS density double'//lf
           end if
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,1)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS velocity_zonal double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,2)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS velocity_meridional double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,3)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS OMEGA double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,4)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS vorticity double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,5)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS topography double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,6)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS penalization double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,7)
           end do
           write(iunit) lf

           if (sim_type == "atm") then
              write(iunit) 'SCALARS surf_press double'//lf
           else
              write(iunit) 'SCALARS free_surface_pert double'//lf
           end if
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,8)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS geopot_height double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,9)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS mask int'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) mask(icell)
           end do
           write(iunit) lf
           
           write(iunit) 'SCALARS level int'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) level(icell)
           end do

           close (iunit)
        else ! 2 layer
           write(iunit) 'SCALARS density double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,1)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_velocity_zonal double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,2)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_velocity_meridional double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,3)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_velocity_zonal double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,4)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_velocity_meridional double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,5)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_vort double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,6)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_vort double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,7)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_free_surf double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,8)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_free_surf double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,9)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS topography double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,10)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS land_mass double'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,11)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS mask int'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) mask(icell)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS level int'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) level(icell)
           end do
        end if

        close(iunit)
     end do
  end do
end program

