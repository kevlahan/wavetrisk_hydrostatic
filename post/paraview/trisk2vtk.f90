program trisk2vtk
  ! Converts trisk data files from wavetrisk code to BINARY .vtk format for paraview
  !
  implicit none
  integer, parameter :: VTK_TRIANGLE         = 5
  integer, parameter :: VTK_POLYGON          = 7
  integer, parameter :: VTK_TETRA            = 10
  integer, parameter :: VTK_WEDGE            = 13
  integer, parameter :: VTK_PENTAGONAL_PRISM = 15
  integer, parameter :: VTK_HEXAGONAL_PRISM  = 16

  integer, parameter :: nvar_out = 11 ! number of saved variables
  integer, parameter :: iunit = 10

  integer                                 :: fid, i, icell, icoord, isave, ivert, k, ncells, ncells_adapt, nvertices, tstart, tend
  integer                                 :: stat
  integer                                 :: j, jmin, jmax, ncells_old, zmax, zmin
  integer, dimension(:),      allocatable :: mask, level, vtk_type, tmp

  real(4), dimension (:,:),   allocatable :: tmp_outv, outv
  real(4), dimension (:,:,:), allocatable :: tmp_vertices, vertices

  character(2)    :: j_lev
  character(4)    :: isv
  character(12)   :: str1, str2
  character(6)    :: grid_type
  character(1300) :: arg, bash_cmd, command, file_in, file_out, run_id, sim_type, zlev

  character(1), parameter :: lf=char(10) ! line feed character

  logical :: compressed_file_exists, file_exists

  ! Get input parameters
  call get_command_argument (1, arg); run_id    = trim(arg)
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
     write (6,'(/,a)')"Usage: trisk2vtk run_id sim_type grid_type tstart tend jmin jmax zmin zmax"
     write (6,'(/,a)')"Example: ./trisk2vtk drake sea hex 0 1 5 7 3 4"
     write (6,'(/,a)')"Reads files drake_003.1.0001.tgz, drake_004.1.0001.tgz, drake_003.1.0002.tgz, drake_004.1.0002.tgz"
     write (6,'(a,/)')"for scale 5 to 7."
     write (6,'(a)') "run_id    = base name of files"
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
     write (6,'(/,"run_id = ", a)') trim(run_id)
     write (6,'("sim_type  = ", a)') trim(sim_type)
     write (6,'("grid_type = ", a)') trim(grid_type)
     write (6,'("tstart    = ", i4.4)') tstart
     write (6,'("tend      = ", i4.4)') tend
     write (6,'("jmin      = ", i2.2)') jmin
     write (6,'("jmax      = ", i2.2)') jmax
     write (6,'("zmin      = ", i3.3)') zmin
     write (6,'("zmax      = ", i3.3)') zmax

     if (trim(grid_type) == "hex") then     ! hexagonal cells (primal grid)
        nvertices = 6
     elseif (trim(grid_type) == "tri") then ! triangular cells (dual grid)
        nvertices = 3 
     end if
  end if

  do isave = tstart, tend
     write (isv,'(i4.4)') isave
     file_in = trim(run_id)//"_"//trim(grid_type)//"_"//trim(isv)

     ! Check if files are compressed
     inquire (FILE=trim(file_in)//'.tgz', EXIST=compressed_file_exists)
     if (compressed_file_exists) then
        command = 'tar xf '//trim(file_in)//".tgz"
        write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
        call system (trim(bash_cmd))
     end if
     write (6,'(/,"Reading from archive ",a,/)') trim(file_in)//".tgz"

     do k = zmin, zmax
        write (zlev,'(i3.3)') k
        ! Loop through all scales
        ncells = 0
        do j = jmin, jmax
           ! First pass to determine number of cells
           write (j_lev,'(i2.2)') j

           file_in = trim(run_id)//"_"//trim(grid_type)//"_"//trim(j_lev)//"_"//trim(zlev)//"_"//trim(isv)
           write (6,'("Loading file ",a)') trim (file_in)
           open (iunit, file=trim(file_in), form="unformatted")

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
              allocate (tmp_vertices(1:ncells,1:nvertices,1:3))
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
              if (allocated (vertices)) deallocate (vertices)
              if (allocated (mask))     deallocate (mask)
              if (allocated (level))    deallocate (level)
              if (allocated (level))    deallocate (level)
              if (allocated (outv))     deallocate (outv)
              allocate (vertices(1:ncells,1:nvertices,1:3), level(1:ncells), mask(1:ncells), outv(1:ncells,1:nvar_out))
           end if

           if (allocated(vtk_type)) deallocate (vtk_type)
           allocate (vtk_type(1:ncells)); vtk_type = VTK_POLYGON

           ! Second pass to read in data
           open (iunit, file=trim(file_in), form="unformatted")
           do icell = ncells_old+1, ncells
              read (iunit) ((vertices(icell,ivert,icoord),icoord=1,3), ivert = 1, nvertices), &
                   outv(icell,1:nvar_out), mask(icell), level(icell)
           end do
           close (iunit)
        end do

        ! Write out data in old vtk format
        file_out = trim(run_id)//"_"//trim(grid_type)//"_"//trim(zlev)//"_"//trim(isv)//".vtk"
        write (6,'("Saving output to file ", a,/)') trim (file_out)

        open (unit=iunit, file=trim(file_out), form="unformatted", access='stream',status='replace',convert='BIG_ENDIAN')
        write (iunit) '# vtk DataFile Version 2.0'//lf
        write (iunit) 'vtk output'//lf              
        write (iunit) 'BINARY'//lf                   
        write (iunit) 'DATASET UNSTRUCTURED_GRID'//lf

        ! Write out vertices
        write (str1(1:12),'(i12)') ncells * nvertices
        write (iunit) 'POINTS ' // trim(str1) // ' float'//lf
        do icell = 1, ncells
           do ivert = 1, nvertices
              write (iunit) (vertices(icell,ivert,icoord), icoord=1,3)
           end do
        end do

        ! Write out cells
        write (str1(1:12),'(i12)') ncells
        write (str2(1:12),'(i12)') ncells * (1 + nvertices)
        write (iunit) 'CELLS '//str1//str2//lf
        do icell = 1, ncells
           write (iunit) nvertices, ((icell-1)*nvertices+ivert-1, ivert=1,nvertices)
        end do

        ! Write out type of each cell (all polygons!)
        write (iunit) 'CELL_TYPES '//trim(str1)//lf
        write (iunit) vtk_type

        ! Write out cell data for each cell                        
        write (iunit) 'CELL_DATA '//trim(str1)//lf

        ! Variables
        if (.not. grid_type == "2layer") then ! usual (hex or tri)
           write (iunit) 'SCALARS rho_dz float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,1)

           if (sim_type == "atm") then
              write (iunit) 'SCALARS temperature float'//lf
           else
              write (iunit) 'SCALARS density float'//lf
           end if
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,2)

           write (iunit) 'SCALARS velocity_zonal float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,3)

           write (iunit) 'SCALARS velocity_meridional float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,4)

           write (iunit) 'SCALARS OMEGA float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,5)

           write (iunit) 'SCALARS vorticity float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,6)

           write (iunit) 'SCALARS topography float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,7)

           write (iunit) 'SCALARS penalization float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,8)

           if (sim_type == "atm") then
              write (iunit) 'SCALARS surf_press float'//lf
           else
              write (iunit) 'SCALARS free_surface_pert float'//lf
           end if
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,9)

           write (iunit) 'SCALARS geopot_height float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,10)

           write (iunit) 'SCALARS pressure float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,11)

           write (iunit) 'SCALARS mask int'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) mask

           write (iunit) 'SCALARS level int'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) level
        else ! 2 layer
           write (iunit) 'SCALARS density float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,1)

           write (iunit) 'SCALARS barotropic_velocity_zonal float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,2)

           write (iunit) 'SCALARS barotropic_velocity_meridional float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,3)

           write (iunit) 'SCALARS baroclinic_velocity_zonal float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,4)

           write (iunit) 'SCALARS baroclinic_velocity_meridional float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,5)

           write (iunit) 'SCALARS baroclinic_vort float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,6)

           write (iunit) 'SCALARS barotropic_vort float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,7)

           write (iunit) 'SCALARS baroclinic_free_surf float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,8)

           write (iunit) 'SCALARS barotropic_free_surf float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,9)

           write (iunit) 'SCALARS topography float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,10)

           write (iunit) 'SCALARS land_mass float'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) outv(:,11)

           write (iunit) 'SCALARS mask int'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) mask

           write (iunit) 'SCALARS level int'//lf
           write (iunit) 'LOOKUP_TABLE default'//lf
           write (iunit) level
        end if
        close (iunit)
     end do

     ! Delete uncompressed file for current time
     if (compressed_file_exists) then
        command = "\rm "//trim(run_id)//"_"//trim(grid_type)//"*"//trim(isv)
        write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
        call system (trim(bash_cmd))
     end if
  end do
end program trisk2vtk

