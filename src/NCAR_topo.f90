module NCAR_topo_mod
  
  ! Routines to assign, load and save NCAR topography data.

  use kind_mod,   only : dp
  use shared_mod, only : N_BDRY, N_GLO_DOMAIN, NONE, max_level, min_level, topo_file, topo_min_level, topo_max_level, z_null

  use arch_mod,       only : barrier, glo_id, n_process, rank
  use comm_mpi_mod,   only : update_bdry
  use domain_mod,     only : Domain, grid, exner_fun, idx, topography, topography_data, mass, scalar, velo1, velo2
  use domain_ops_mod, only : apply_onescale_to_patch
  use geom_mod,       only : dist
  
  implicit none

  private
  public :: assign_NCAR_topo, load_topo, save_topo


  integer, allocatable :: topo_count(:,:)

  
contains

  
  subroutine save_topo
    ! Save topography data
    ! (one file per domain and per level)
    !
    ! !! saves topgraphy data on a non-adaptive grid !!
    !

    implicit none

    integer                   :: d, d_glo, j, l
    character(2)              :: l2
    character(5)              :: d5
    character(9999)           :: filename
    character(:), allocatable :: archive, cmd, files

    call update_bdry (topography, NONE, 966)

    allocate (topo_count(min_level:max_level,1:size(grid))); topo_count = 0

    ! Save a separate file for each domain and each level
    do l = max_level, min_level, -1
       do d = 1, size(grid)
          d_glo = glo_id(rank+1,d) + 1

          write(l2,'(i2.2)') l
          write(d5,'(i5.5)') d_glo

          filename = trim(topo_file) // '.' // l2 // '.' // d5

          open (unit=10, file=filename, form="UNFORMATTED", action='WRITE', status='REPLACE')

          mass   => exner_fun(1)%data(d)%elts
          scalar => topography%data(d)%elts
          velo1  => grid(d)%u_zonal%elts
          velo2  => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (write_topo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (mass, scalar, velo1, velo2)
          close (10)
       end do
    end do

    write (filename, '(a,a)') trim (topo_file), '.count'
    open (unit=10, file=trim (filename), form="UNFORMATTED", action='WRITE', status='REPLACE')
    write (10) topo_count
    close (10)

    ! Compress topography data
    archive = trim (topo_file)//'.tgz'
    write (6, '(a,a)') 'Saving topography file ', archive

    files = trim(topo_file) // '.??.????? ' // trim(filename)

    cmd = 'gtar czf ' // archive // ' ' // files // ' --remove-files'
    call execute_command_line (cmd)
  end subroutine save_topo


  subroutine write_topo (dom, i, j, zlev, offs, dims)
    ! Write out coordinates, topography height, topography gradients and surface pressure
    ! Compute topo_count

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, l

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    l = dom%level%elts(id) ! level

    ! Write out coordinates, topography height and topography gradients
    write (10) dom%node%elts(id), scalar(id), velo1(id), velo2(id), mass(id)

    topo_count(l,d) = topo_count(l,d) + 1
  end subroutine write_topo


  subroutine load_topo
    ! Read topography data from for restart
    ! (one file per domain)
    !
    ! !! assumes topgraphy data was saved on a non-adaptive grid !!
    !

    implicit none

    integer         :: d, d_glo, i, ii, l, r
    character(9999) :: filename
    character(:), allocatable :: archive, cmd, files

    ! Uncompress topography data
    if (rank == 0) then
       archive = trim(topo_file) // '.tgz'
       write (6,'(a,a)') 'Loading topography file ', archive
       cmd = 'gtar xzf ' // archive
       call execute_command_line(cmd)
    end if
    call barrier

    ! Allocate topo_count matrix for all domains on all ranks 
    if (allocated (topo_count)) deallocate (topo_count)
    allocate (topo_count(topo_min_level:topo_max_level,1:N_GLO_DOMAIN))

    do r = 1, n_process

       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call barrier
          cycle 
       end if

       ! Read topography count for all domains
       write (filename, '(a,a)') trim (topo_file), '.count'
       open (unit=10, file=trim (filename), form="UNFORMATTED", action='READ', status='OLD')
       read (10) topo_count
       close (10)

       ! Allocate topgraphy data structures    
       allocate (topography_data(topo_min_level:topo_max_level, 1:size(grid)))
       do l = topo_min_level, topo_max_level
          do d = 1, size(grid)
             d_glo = glo_id(rank+1,d) + 1
             allocate (topography_data(l,d)%node(1:  topo_count(l,d_glo)))
             allocate (topography_data(l,d)%elts(1:4*topo_count(l,d_glo)))
          end do
       end do
    end do

    ! Load topography data
    do l = topo_min_level, topo_max_level
       do d = 1, size(grid)
          d_glo = glo_id(rank+1,d) + 1
          write (filename, '(a,a,i2.2,a,i5.5)') trim (topo_file), '.', l, '.', d_glo
          open (unit=10, file=trim (filename), form="UNFORMATTED", action='READ', status='OLD')
          do ii = 1, topo_count(l,d_glo)
             i = 4*(ii-1) + 1
             read (10) topography_data(l,d)%node(ii), topography_data(l,d)%elts(i:i+3)
          end do
          close (10)
       end do
    end do

    ! Remove temporary files
    call barrier ! do not delete files before everyone has read them
    if (rank == 0) then
       files = trim(topo_file) // '.??.????? ' // trim(topo_file) // '.count'
       cmd   = '\rm -f ' // files
       call execute_command_line (cmd)
    end if
  end subroutine load_topo


  subroutine assign_NCAR_topo (dom, i, j, zlev, offs, dims)
    ! Assign topography data to topography structure for simulation
    ! Sets topography height and surface pressure

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                             :: d, id, ii, jj, l, n_topo
    real(dp), dimension(:), allocatable :: distance

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    l  = dom%level%elts(id)

    n_topo = size (topography_data(l,d)%node); allocate (distance(1:n_topo))
    do ii = 1, n_topo
       distance(ii) = dist (dom%node%elts(id), topography_data(l,d)%node(ii))
    end do
    jj = minloc (distance,1) ; deallocate (distance)

    ii = 4*(jj-1) + 1
    topography%data(d)%elts(id) = topography_data(l,d)%elts(ii)
    dom%surf_press%elts(id)     = topography_data(l,d)%elts(ii+3)
  end subroutine assign_NCAR_topo

  
end module NCAR_topo_mod
