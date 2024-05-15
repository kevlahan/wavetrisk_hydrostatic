program trisk2vtk
  ! Converts standard ASCII trisk data files from wavetrisk code to BINARY .vtk format for paraview
  !
  ! Usage: trisk2vtk file_base file_type tstart tend jmin jmax zmin zmax
  ! Example: trisk2vtk drake primal 0 1 5 7 3 4
  !
  ! Reads files drake_003.1.0001.tgz, drake_004.1.0001.tgz, drake_003.1.0002.tgz, drake_004.1.0002.tgz
  ! for scale 5 to 7.
  !
  ! If zmin < 0 there is no layer label (for backwards compatibility).
  !
  ! file_base = name of file
  ! file_type = primal (hexagons) or multi (triangles)
  ! tstart    = first file to read
  ! tend      = last file to read
  ! jmin      = minimum scale
  ! jmax      = maximum scale
  ! zmin      = min vertical layer
  ! zmax      = max vertical layer

  implicit none
  integer :: fid, i, icell, icoord, ipts, itime, itype, ivert, k, ncells, ncells_adapt, n_vertices, tstart, tend
  integer :: stat
  integer :: j, jmin, jmax, ncells_old, nvar_out, zmax, zmin
  integer, parameter :: iunit=10
  integer, parameter :: vtk_type = 7 ! vtk polygonal data
  integer, dimension(:), allocatable :: tmp, mask, level

  character(2)   :: j_lev
  character(3)   :: zlev
  character(4)   :: s_time
  character(12)  :: str1, str2
  character(6)   :: file_type
  character(255) :: arg, bash_cmd, command, filename_in, filename_out, file_base
  character(1), parameter :: lf=char(10) ! line feed character

  real(4), dimension (:,:),   allocatable :: tmp_outv, outv
  real(4), dimension (:,:,:), allocatable :: tmp_vertices, vertices

  real(4), parameter :: Hdim = 3.344175893265152e+03, g = 9.80665
  real(4), parameter :: Ldim = 6371e3
  real(4), parameter :: Udim = sqrt(Hdim*g)
  real(4), parameter :: Tdim = Ldim/Udim

  logical :: compressed_file_exists, file_exists

  ! Get input parameters
  call get_command_argument(1, arg); file_base = trim(arg)
  call get_command_argument(2, arg); file_type = trim(arg)
  call get_command_argument(3, arg); read (arg,'(I12)') tstart
  call get_command_argument(4, arg); read (arg,'(I12)') tend
  call get_command_argument(5, arg); read (arg,'(I12)') jmin
  call get_command_argument(6, arg); read (arg,'(I12)') jmax
  call get_command_argument(7, arg); read (arg,'(I12)') zmin
  call get_command_argument(8, arg); read (arg,'(I12)') zmax
  
  if (file_type == " ") then
     write (6,'(/,a)') "Converts standard ASCII trisk data files from wavetrisk code to BINARY .vtk format for paraview"
     write (6,'(/,a)')"Usage: trisk2vtk file_base file_type tstart tend jmin jmax zlev"
     write (6,'(/,a)')"Example: ./trisk2vtk drake primal 0 1 5 7 3 4"
     write (6,'(/,a)')"Reads files drake_003.1.0001.tgz, drake_004.1.0001.tgz, drake_003.1.0002.tgz, drake_004.1.0002.tgz"
     write (6,'(a,/)')"for scale 5 to 7."
     write (6,'(a,/)') "If zmin < 0 there is no layer label (for backwards compatibility)."
     write (6,'(a)') "file_base = name of file"
     write (6,'(a)') "file_type = primal (hexagons) or multi (multiscale grid of triangles)"
     write (6,'(a)') "tstart    = first file to read"
     write (6,'(a)') "tend      = last file to read"
     write (6,'(a)') "jmin      = minimum scale"
     write (6,'(a)') "jmax      = maximum scale"
     write (6,'(a)') "zmin      = min vertical layer"
     write (6,'(a,/)') "zmax      = max vertical layer"
     stop
  else
     write (6,'(/,"file_base = ", a)') trim(file_base)
     write (6,'("file_type = ", a)') trim(file_type)
     write (6,'("tstart    = ", i4.4)') tstart
     write (6,'("tend      = ", i4.4)') tend
     write (6,'("jmin      = ", i2.2)') jmin
     write (6,'("jmax      = ", i2.2)') jmax
     write (6,'("zmin      = ", i4.3)') zmin
     if (zmin >= 0) write (6,'("zmax      = ", i3.3)') zmax

     if (trim(file_type) == "primal") then
        itype      = 1 
        n_vertices = 6 ! hexagonal cells (primal grid)
        nvar_out   = 7
     elseif (trim(file_type) == "multi") then
        itype      = 2
        n_vertices = 3 ! triangular cells (multiscale grid)
        nvar_out   = 8
     elseif (trim(file_type) == "2layer") then
        itype      = 1
        zmin       = -1
        n_vertices = 6 ! hexagonal cells (primal grid)
        nvar_out   = 11
     end if

     if (zmin < 0) zmax = zmin 
  end if

  do itime = tstart, tend
     write (s_time,'(I4.4)') itime

     do k = zmin, zmax
        if (itype == 1) then
           if (zmin < 0) then
              filename_in = trim(file_base)//'.1'//s_time
           else
              write (zlev,'(I3.3)') k
              filename_in = trim(file_base)//'_'//zlev//'.1'//s_time
           end if
        elseif (itype == 2) then
           if (zmin < 0) then
              filename_in = trim(file_base)//'.2'//s_time
           else
              write (zlev,'(I3.3)') k
              filename_in = trim(file_base)//'_'//zlev//'.2'//s_time
           end if
        end if

        ! Check if files are compressed
        inquire (FILE=trim(filename_in)//'.tgz', EXIST=compressed_file_exists)
        if (compressed_file_exists) then 
           command = 'tar xzf '//trim(filename_in)//'.tgz'
           write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
           call system (trim(bash_cmd))

           ! Delete un-needed file
           if (trim(file_type) == "primal" .or. trim(file_type) == "2layer") then
              command = '\rm ' // trim(filename_in) // '00'
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
           open (unit=iunit, file=trim(filename_in)//j_lev, form="formatted")

           ncells_old = ncells
           stat = 0
           do while (stat>=0)
              read (iunit, IOSTAT=stat,FMT=*)
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
              if(allocated(vertices)) deallocate(vertices)
              if(allocated(mask)) deallocate(mask)
              if(allocated(level)) deallocate(level)
              if(allocated(level)) deallocate(level)
              if(allocated(outv)) deallocate(outv)
              allocate (vertices(1:ncells,1:n_vertices,1:3), level(1:ncells), mask(1:ncells), outv(1:ncells,1:nvar_out))
           end if

           ! Second pass to read in data
           open (unit=iunit, file=trim(filename_in)//j_lev, form="formatted")
           if (file_type == "primal") then
              do icell = ncells_old+1, ncells
                 read (iunit, fmt='(18(E14.5E2, 1X), 7(E14.5E2, 1X), I3, 1X, I3)') &
                      ((vertices(icell,ivert,icoord),icoord=1,3),ivert=1,n_vertices), &
                      outv(icell,1:nvar_out), mask(icell), level(icell)
              end do
           elseif (file_type == "multi") then
              do icell = ncells_old+1, ncells
                 read (iunit, fmt='(9(e14.5e2, 1x), 8(e14.5e2, 1x), i3)') &
                      ((vertices(icell,ivert,icoord),icoord=1,3),ivert=1,n_vertices), &
                      outv(icell,1:nvar_out), level(icell)
              end do
           elseif (file_type == "2layer") then
              do icell = ncells_old+1, ncells
                 read (iunit, fmt='(18(E14.5E2, 1X), 11(E14.5E2, 1X), I3, 1X, I3)') &
                      ((vertices(icell,ivert,icoord),icoord=1,3),ivert=1,n_vertices), &
                      outv(icell,1:nvar_out), mask(icell), level(icell)
              end do
           end if
           close(iunit)

           ! Delete uncompressed file
           if (compressed_file_exists) then
              command = '\rm ' // trim(filename_in)//j_lev
              write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
              call system (trim(bash_cmd))
           end if
        end do

        ! Write out data in old vtk format
        write (6,'("Saving output to file ", a,/)') trim(filename_in)//'.vtk'

        open (unit=iunit, file=trim(filename_in)//'.vtk', form="unformatted", access='stream',convert='BIG_ENDIAN')
        write(iunit) '# vtk DataFile Version 3.0'//lf
        write(iunit) 'vtk output'//lf                
        write(iunit) 'BINARY'//lf                    
        write(iunit) 'DATASET UNSTRUCTURED_GRID'//lf  

        write(str1(1:12),'(i12)') ncells * n_vertices
        write(iunit) 'POINTS ' // trim(str1) // ' float'// lf

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
        
        if (file_type == "primal") then
           write(iunit) 'SCALARS temperature float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,1)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS velocity_zonal float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,2)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS velocity_meridional float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,3)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS geopot_height float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,4)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS mass float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,5)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS surf_press float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,6)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS vorticity float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,7)
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
        elseif (file_type == "2layer") then
           write(iunit) 'SCALARS density float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,1)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_velocity_zonal float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,2)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_velocity_meridional float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,3)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_velocity_zonal float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,4)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_velocity_meridional float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,5)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_vort float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,6)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_vort float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,7)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS baroclinic_free_surf float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,8)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS barotropic_free_surf float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,9)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS topography float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,10)
           end do
           write(iunit) lf

           write(iunit) 'SCALARS land_mass float'//lf
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
        elseif (file_type == "multi") then
           write(iunit) 'SCALARS temperature float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,1) 
           end do
           write(iunit) lf

           write(iunit) 'SCALARS zonal_velocity float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,2) 
           end do
           write(iunit) lf

           write(iunit) 'SCALARS meridional_velocity float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,3) 
           end do
           write(iunit) lf

           write(iunit) 'SCALARS vorticity float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,4) 
           end do
           write(iunit) lf

           write(iunit) 'SCALARS topography float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,5) 
           end do
           write(iunit) lf

           write(iunit) 'SCALARS OMEGA float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,6) 
           end do
           write(iunit) lf

           write(iunit) 'SCALARS surface_pressure float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,7) 
           end do
           write(iunit) lf
           
           write(iunit) 'SCALARS kinetic_energy float'//lf
           write(iunit) 'LOOKUP_TABLE default'//lf
           do icell = 1, ncells
              write (iunit) outv(icell,8) 
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
end program trisk2vtk

