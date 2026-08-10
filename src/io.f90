module io_mod

  use mpi_f08
  use iso_fortran_env, only: int8, int64

  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, Coord, EDGE, N_BDRY,  N_CHDRN, N_VARIABLE, TRIAG, RT, DG, UP, MATH_PI, MULT, &
       cp_idx, istep_cumul, itime, iwrite, &
       level_end, min_level, max_level, n_glo_domain, grav_accel, kappa, p_0, R_d, pressure_save, radius, ref_density, run_id, &
       scalars, threshold, time, topo_file, topo_min_level, topo_max_level, vert_diffuse, zmin, zmax, zlevels, &
       S_MASS, S_TEMP, S_VELO, NONE, LORT, UPLT, z_null

  use arch_mod,         only : barrier, comm, glo_id, n_process, rank, cp_load
  use comm_mod,         only : domain_load
  use comm_mpi_mod,     only : update_bdry
  use domain_mod,       only : init_domain_mod
  use domain_ops_mod,   only : apply_interscale_d, apply_onescale_to_patch, apply_to_pole_d 
  use geom_mod,         only : dist
  use init_mod,         only : dump, load, surf_geopot
  use patch_mod,        only : PATCH_SIZE
  use refine_patch_mod, only : check_child_required, post_refine, refine_patch1, refine_patch2
  use utils_mod,        only : integrate_hex, interp
  use wavelet_mod,      only : Restrict_scalar

  use domain_mod, only : Domain, Float_Field, exner_fun, grid, idx, sol, sol_mean, tke, topography, topography_data, &
       scalar, trend, ke, mass, velo, velo1, velo2, vort, wav_coeff, wav_tke, wc_s

  implicit none 

  private
  public :: dump_adapt_mpi, load_adapt_mpi, read_checkpoint_directory
  public :: assign_NCAR_topo,  kinetic_energy, init_io_mod, load_topo, save_topo, vort_extrema
  public :: topo, pot_energy, total_ke, layer1_ke, layer2_ke, one_layer_ke, barotropic_ke, pot_enstrophy, umag


  ! Magic number identifying a WAVETRISK checkpoint file.
  ! (ASCII encoding of the first eight characters of WAVETRISK)
  integer(int64), parameter :: CP_MAGIC = &
       int(z'5741564554524953', int64)

  ! Version number of the checkpoint file layout.
  integer(int64), parameter :: CP_VERSION = 1_int64

  ! Absolute byte offset of each global domain within the checkpoint file.
  integer(int64), allocatable :: cp_offset(:)

  ! Serialized size (bytes) of each global domain.
  integer(int64), allocatable :: cp_nbytes(:)

  ! Size of the fixed checkpoint header (magic number, version,
  ! number of global domains).
  integer(int64), parameter :: CP_HEADER_BYTES = 24_int64

  ! Byte offset of the domain load array in the checkpoint file.
  integer(int64), parameter :: CP_LOAD_POS = &
       CP_HEADER_BYTES

  ! Byte offset of the domain offset array in the checkpoint file.
  integer(int64), parameter :: CP_OFFSET_POS = &
       CP_LOAD_POS + &
       int(storage_size(0)/8, int64) * int(N_GLO_DOMAIN, int64)

  ! Byte offset of the domain size array in the checkpoint file.
  integer(int64), parameter :: CP_NBYTES_POS = &
       CP_OFFSET_POS + &
       8_int64 * int(N_GLO_DOMAIN, int64)

  ! Byte offset of the serialized domain data.
  integer(int64), parameter :: CP_DATA_POS = &
       CP_NBYTES_POS + &
       8_int64 * int(N_GLO_DOMAIN, int64)


  integer, dimension(:,:), allocatable :: topo_count
  real(dp)                             :: vmin, vmax


contains


  subroutine init_io_mod

    implicit none

    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_domain_mod
    initialized = .true.
  end subroutine init_io_mod


  subroutine vort_extrema (dom, i, j, zlev, offs, dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id, idN, idE
    real(dp) :: vort

    id  = idx(i,   j,   offs, dims)
    idN = idx(i,   j+1, offs, dims)
    idE = idx(i+1, j,   offs, dims)

    if ( dom%mask_e%elts(id*EDGE+DG+1)  >= ADJZONE .or. &
         dom%mask_e%elts(id*EDGE+UP+1)  >= ADJZONE .or. &
         dom%mask_e%elts(idN*EDGE+RT+1) >= ADJZONE) then

       vort = dom%vort%elts(id*TRIAG+UPLT+1)
       vmin = min(vmin, vort)
       vmax = max(vmax, vort)
    end if

    if ( dom%mask_e%elts(id*EDGE+DG+1)  >= ADJZONE .or. &
         dom%mask_e%elts(idE*EDGE+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(id*EDGE+RT+1)  >= ADJZONE) then
       vort = dom%vort%elts(id*TRIAG+LORT+1)
       vmin = min (vmin, vort)
       vmax = max (vmax, vort)
    end if
  end subroutine vort_extrema


  function topo (dom, i, j, zlev, offs, dims) result(val)

    use domain_mod

    implicit none
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    val = topography%data(d)%elts(id)
  end function topo


  function pot_energy (dom, i, j, zlev, offs, dims) result(val)

    use domain_mod

    implicit none
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: id
    real(dp) :: rho_dz

    id = idx (i, j, offs, dims)

    rho_dz = sol_mean(S_MASS,zlev)%data(dom%id+1)%elts(id+1) + sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)
    val = rho_dz**2
  end function pot_energy


  function total_ke (itype) result(val)
    ! Computes total kinetic energy

    implicit none

    character(len=*), intent(in) :: itype
    real(dp)                     :: val

    integer :: k

    val = 0.0_dp
    do k = 1, zlevels
       val = val + integrate_hex (kinetic_energy, k)
    end do
  end function total_ke


  function kinetic_energy (dom, i, j, zlev, offs, dims) result(val)
    ! Kinetic energy u^2/2 at level zlev using approximation to TRiSK formula

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    d  = dom%id + 1
    id = idx (i, j, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    u_prim_RT = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
    u_prim_DG = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
    u_prim_UP = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

    u_dual_RT = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

    u_prim_UP_S  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
    u_prim_DG_SW = sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
    u_prim_RT_W  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

    u_dual_RT_W  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
    u_dual_DG_SW = sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
    u_dual_UP_S  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

    val = (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT + &
         u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
         * dom%areas%elts(id+1)%hex_inv/4.0_dp 
  end function kinetic_energy


  subroutine umag (q)
    ! Evaluate complete velocity trend by adding gradient terms to previously calculated source terms on entire grid
    !
    ! Input:  velocity field q at a single vertical layer
    ! Output: velocity magnitude stored in dom%ke

    implicit none

    type(Float_Field), target, intent(in) :: q

    integer :: k

    integer :: d, p

    do d = 1, size(grid)
       velo => q%data(d)%elts
       ke   => grid(d)%ke%elts
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (cal_umag, grid(d), p-1, k, 0, 1)
       end do
       nullify (ke, velo)
    end do
  end subroutine umag

  
  subroutine cal_umag (dom, i, j, zlev, offs, dims)
    ! Velocity magnitude: sqrt(2*ke) using approximation to TRiSK formula
    ! divide out surface area

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    d  = dom%id + 1
    id = idx (i, j, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    u_prim_RT = velo(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
    u_prim_DG = velo(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
    u_prim_UP = velo(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

    u_dual_RT = velo(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG = velo(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP = velo(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

    u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
    u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
    u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

    u_dual_RT_W  = velo(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
    u_dual_UP_S  = velo(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

    ke(id+1) = sqrt( (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT + &
         u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
         * dom%areas%elts(id+1)%hex_inv/(4.0_dp*MATH_PI*radius**2)) 
  end subroutine cal_umag


  function layer1_ke (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       u_prim_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

       u_prim_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
       u_prim_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

       u_dual_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1))  &
            * (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4.0_dp
    else
       val = 0.0_dp
    end if
  end function layer1_ke


  function layer2_ke (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       u_prim_RT = sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = sol(S_VELO,2)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = sol(S_VELO,2)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = sol(S_VELO,2)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = sol(S_VELO,2)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

       u_prim_UP_S  = sol(S_VELO,2)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
       u_prim_DG_SW = sol(S_VELO,2)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT_W  = sol(S_VELO,2)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

       u_dual_RT_W  = sol(S_VELO,2)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = sol(S_VELO,2)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = sol(S_VELO,2)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = (sol_mean(S_MASS,2)%data(d)%elts(id+1) + sol(S_MASS,2)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       val = 0.0_dp
    end if
  end function layer2_ke


  function one_layer_ke (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       u_prim_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

       u_prim_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
       u_prim_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

       u_dual_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       val = 0.0_dp
    end if
  end function one_layer_ke


  function barotropic_ke (dom, i, j, zlev, offs, dims) result(val)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer                     :: d, id, idE, idNE, idN, idW, idSW, idS
    integer,  dimension(1:EDGE) :: id_edge, id_node
    real(dp)                    :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp)                    :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S
    real(dp), dimension(1:EDGE) :: u

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       id_node = [idE, idNE, idN]
       id_edge = [id,  id,   id ]
       u = barotropic_velo ()
       u_prim_RT = u(1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = u(2) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = u(3) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = u(1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = u(2) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = u(3) * dom%pedlen%elts(EDGE*id+UP+1)

       id_node = [ idW, idSW, idS ]
       id_edge = id_node
       u = barotropic_velo ()
       u_prim_RT_W  = u(1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_DG_SW = u(2) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_S  = u(3) * dom%len%elts(EDGE*idS +UP+1)

       u_dual_RT_W  = u(1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = u(2) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = u(3) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = &
            (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1) + &
            sol_mean(S_MASS,2)%data(d)%elts(id+1) + sol(S_MASS,2)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       val = 0.0_dp
    end if

  contains

    function barotropic_velo ()
      real(dp), dimension(1:EDGE) :: barotropic_velo

      integer                     :: e, id_e, k
      real(dp), dimension(1:EDGE) :: dz

      do e = 1, EDGE
         id_e = EDGE*id_edge(e) + e

         dz = 0.0_dp
         do k = 1, 2
            dz(k) = interp (&
                 sol_mean(S_MASS,k)%data(d)%elts(id+1)         + sol(S_MASS,k)%data(d)%elts(id+1), &
                 sol_mean(S_MASS,k)%data(d)%elts(id_node(e)+1) + sol(S_MASS,k)%data(d)%elts(id_node(e)+1))
         end do

         barotropic_velo(e) = (dz(1)*sol(S_VELO,1)%data(d)%elts(id_e) + dz(2)*sol(S_VELO,2)%data(d)%elts(id_e)) / sum(dz)
      end do

    end function barotropic_velo

  end function barotropic_ke

  
  function pot_enstrophy (dom, i, j, zlev, offs, dims) result(val)
    ! Computes potential enstrophy in two layer mode split case
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, id_i
    real(dp) :: f, h, w

    id = idx (i, j, offs, dims)
    id_i = id + 1

    d = dom%id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       ! Approximate Coriolis term
       f = dom%coriolis%elts(TRIAG*id+LORT+1)/dom%triarea%elts(TRIAG*id+LORT+1)
       ! Total vorticity
       w = dom%ke%elts(id_i)  + f
       ! Height
       if (zlev == 3) then ! barotropic
          h = (sol_mean(S_MASS,1)%data(d)%elts(id_i) + sol(S_MASS,1)%data(d)%elts(id_i) + &
               sol_mean(S_MASS,2)%data(d)%elts(id_i) + sol(S_MASS,2)%data(d)%elts(id_i)) / ref_density
       else ! single layer
          h = (sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)) / ref_density
       end if
       val = 0.5_dp * (w / h)**2
    else
       val = 0.0_dp
    end if
  end function pot_enstrophy


  subroutine dump_adapt_mpi (id)
    ! Save adaptive checkpoint data to one shared MPI-IO file.
    !
    ! Each MPI rank serializes all locally owned domains into one
    ! temporary stream.  The rank-local streams are then written
    ! collectively into one shared checkpoint file.
    !
    ! The checkpoint directory contains, for every global domain:
    !
    !   cp_load   : computational load
    !   cp_offset : absolute byte offset in checkpoint file
    !   cp_nbytes : serialized domain size in bytes
    !
    ! After the shared checkpoint has been closed collectively,
    ! rank zero compresses it using zstd -1 and removes the
    ! uncompressed .bin file.
    !
    ! NOTE: Checkpoint generation modifies the adaptive grid structure
    ! by deleting patches that are not required in the adjacent zone.

    implicit none

    integer, intent(in) :: id

    integer :: d, k, v
    integer :: fid
    integer :: ierr
    integer :: cmdstat, exitstat

    type(MPI_File)   :: fh
    type(MPI_Status) :: status

    integer(int64) :: p0, p1
    integer(int64) :: rank_bytes
    integer(int64) :: rank_disp

    integer,        allocatable :: gid_loc(:)
    integer,        allocatable :: load_loc(:)
    integer(int64), allocatable :: off_loc(:)
    integer(int64), allocatable :: len_loc(:)

    integer(int8), allocatable :: buf(:)

    character(4)              :: cp4
    character(:), allocatable :: filename
    character(:), allocatable :: cmd


    ! ---------------------------------------------------------------
    ! Prepare checkpoint data
    ! ---------------------------------------------------------------

    call update_bdry ( &
         wav_coeff(scalars(1):scalars(2),zmin:zmax), &
         NONE, 964)

    if (vert_diffuse) then
       call update_bdry (wav_tke, NONE, 965)
    end if


    ! Restrict physical solution to coarsest scale.

    do d = 1, size(grid)

       do k = zmin, zmax
          do v = scalars(1), scalars(2)

             scalar => sol(v,k)%data(d)%elts
             wc_s   => wav_coeff(v,k)%data(d)%elts

             call apply_interscale_d ( &
                  Restrict_scalar, &
                  grid(d), &
                  min_level-1, &
                  k, 0, 1)

             nullify (scalar, wc_s)

          end do
       end do

       if (vert_diffuse) then

          do k = 1, zlevels

             scalar => tke(k)%data(d)%elts
             wc_s   => wav_tke(k)%data(d)%elts

             call apply_interscale_d ( &
                  Restrict_scalar, &
                  grid(d), &
                  min_level-1, &
                  k, 0, 1)

             nullify (scalar, wc_s)

          end do

       end if

    end do


    ! ---------------------------------------------------------------
    ! Serialize all locally owned domains into one scratch stream
    ! ---------------------------------------------------------------

    allocate(gid_loc(size(grid)))
    allocate(load_loc(size(grid)))
    allocate(off_loc(size(grid)))
    allocate(len_loc(size(grid)))

    open( &
         newunit=fid, &
         status='scratch', &
         form='unformatted', &
         access='stream', &
         action='readwrite')


    do d = 1, size(grid)

       gid_loc(d)  = glo_id(rank+1,d)
       load_loc(d) = domain_load(grid(d))

       ! Position immediately before this domain.
       inquire(unit=fid, pos=p0)

       ! Complete stage-2 domain serialization.
       call dump_domain(fid, d)

       ! Position immediately after this domain.
       inquire(unit=fid, pos=p1)

       ! Fortran stream positions are one based.
       off_loc(d) = p0 - 1_int64
       len_loc(d) = p1 - p0

    end do


    ! Total rank-local serialized payload size.

    inquire(unit=fid, pos=p1)

    rank_bytes = p1 - 1_int64


    ! ---------------------------------------------------------------
    ! Determine global payload offsets
    ! ---------------------------------------------------------------

    rank_disp = 0_int64

    call MPI_Exscan( &
         rank_bytes, rank_disp, &
         1, MPI_INTEGER8, MPI_SUM, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_Exscan failed"
    end if

    ! MPI_Exscan receive value is undefined on rank zero.
    if (rank == 0) rank_disp = 0_int64


    ! Convert domain offsets from rank-local stream offsets to
    ! absolute byte offsets in the shared checkpoint file.

    off_loc = CP_DATA_POS + rank_disp + off_loc


    ! ---------------------------------------------------------------
    ! Construct global checkpoint directory
    ! ---------------------------------------------------------------

    call build_checkpoint_directory( &
         gid_loc, load_loc, off_loc, len_loc)


    ! build_checkpoint_directory allocates cp_* on every rank but
    ! fills them only on rank zero.  Broadcast the completed
    ! directory so distribute_grid/load_adapt_mpi can use it
    ! immediately without rereading the checkpoint.

    call MPI_Bcast( &
         cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: cp_load broadcast failed"
    end if

    call MPI_Bcast( &
         cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: cp_offset broadcast failed"
    end if

    call MPI_Bcast( &
         cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: cp_nbytes broadcast failed"
    end if


    ! ---------------------------------------------------------------
    ! Copy rank-local scratch stream into memory
    ! ---------------------------------------------------------------

    if (rank_bytes > int(huge(0),int64)) then
       error stop "dump_adapt_mpi: rank checkpoint buffer too large"
    end if

    allocate(buf(int(rank_bytes)))

    if (rank_bytes > 0_int64) then

       flush(fid)

       read(fid, pos=1) buf

    end if

    close(fid)


    ! ---------------------------------------------------------------
    ! Create shared checkpoint file
    ! ---------------------------------------------------------------

    write(cp4,'(i4.4)') id

    filename = &
         trim(run_id) // "_checkpoint_" // cp4 // ".bin"


    call MPI_File_open( &
         comm, &
         filename, &
         MPI_MODE_CREATE + MPI_MODE_WRONLY, &
         MPI_INFO_NULL, &
         fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_open failed"
    end if


    ! Truncate an existing checkpoint with the same name.
    ! MPI_File_set_size is collective.

    call MPI_File_set_size( &
         fh, &
         int(0, MPI_OFFSET_KIND), &
         ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_set_size failed"
    end if


    ! Rank zero writes header and directory.

    if (rank == 0) then
       call write_checkpoint_directory(fh)
    end if


    ! Make sure directory writing has completed before payload writes.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_Barrier failed"
    end if


    ! ---------------------------------------------------------------
    ! Write rank-local payloads collectively
    ! ---------------------------------------------------------------

    call MPI_File_write_at_all( &
         fh, &
         int(CP_DATA_POS + rank_disp, MPI_OFFSET_KIND), &
         buf, &
         int(rank_bytes), &
         MPI_BYTE, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_write_at_all failed"
    end if


    call MPI_File_close(fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_close failed"
    end if


    ! ---------------------------------------------------------------
    ! Compress completed checkpoint
    ! ---------------------------------------------------------------
    !
    ! All MPI ranks must have closed the shared file before rank zero
    ! starts compression.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: pre-zstd barrier failed"
    end if


    if (rank == 0) then

       cmd = 'zstd -q -1 --rm -f "' // trim(filename) // '"'

       call execute_command_line( &
            cmd, &
            exitstat=exitstat, &
            cmdstat=cmdstat)

       if (cmdstat /= 0) then
          error stop "dump_adapt_mpi: failed to execute zstd"
       end if

       if (exitstat /= 0) then
          error stop "dump_adapt_mpi: zstd compression failed"
       end if

    end if


    ! Ensure compression has completed before any rank proceeds.
    ! The checkpoint now exists as filename//".zst".

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: post-zstd barrier failed"
    end if


    ! ---------------------------------------------------------------
    ! Cleanup
    ! ---------------------------------------------------------------

    deallocate(buf)

    deallocate(gid_loc)
    deallocate(load_loc)
    deallocate(off_loc)
    deallocate(len_loc)

  end subroutine dump_adapt_mpi


  subroutine build_checkpoint_directory (gid_loc, load_loc, off_loc, len_loc)

    ! Gather global-domain IDs, loads, offsets and serialized sizes
    ! onto rank zero and construct checkpoint directory arrays in
    ! global-domain order.
    !
    ! cp_load, cp_offset and cp_nbytes are allocated on ALL ranks
    ! here so that dump_adapt_mpi can broadcast them immediately
    ! after this routine returns.

    implicit none

    integer,        intent(in) :: gid_loc(:)
    integer,        intent(in) :: load_loc(:)
    integer(int64), intent(in) :: off_loc(:)
    integer(int64), intent(in) :: len_loc(:)

    integer :: sz
    integer :: n_tot
    integer :: ierr
    integer :: i, r, gid

    integer :: rcounts(n_process)
    integer :: displs(n_process)

    integer, allocatable :: gid_glo(:)
    integer, allocatable :: load_recv(:)

    integer(int64), allocatable :: off_recv(:)
    integer(int64), allocatable :: len_recv(:)

    logical, allocatable :: seen(:)


    ! ---------------------------------------------------------------
    ! Allocate checkpoint directory arrays on every rank
    ! ---------------------------------------------------------------

    if (allocated(cp_load))   deallocate(cp_load)
    if (allocated(cp_offset)) deallocate(cp_offset)
    if (allocated(cp_nbytes)) deallocate(cp_nbytes)

    allocate(cp_load(N_GLO_DOMAIN))
    allocate(cp_offset(N_GLO_DOMAIN))
    allocate(cp_nbytes(N_GLO_DOMAIN))

    ! Initialize for easier debugging.  Only rank zero fills them
    ! before the subsequent broadcasts.

    cp_load   = 0
    cp_offset = 0_int64
    cp_nbytes = 0_int64


    ! ---------------------------------------------------------------
    ! Gather number of domains contributed by each rank
    ! ---------------------------------------------------------------

    sz = size(gid_loc)

    call MPI_Gather( &
         sz, 1, MPI_INTEGER, &
         rcounts, 1, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "build_checkpoint_directory: MPI_Gather failed"
    end if


    ! ---------------------------------------------------------------
    ! Construct Gatherv receive layout
    ! ---------------------------------------------------------------

    if (rank == 0) then

       displs(1) = 0

       do r = 2, n_process
          displs(r) = displs(r-1) + rcounts(r-1)
       end do

       n_tot = sum(rcounts)

       if (n_tot /= N_GLO_DOMAIN) then
          write(*,'(A,I0,A,I0)') &
               "build_checkpoint_directory: expected ", &
               N_GLO_DOMAIN, &
               " domains, received ", &
               n_tot

          error stop
       end if

       allocate(gid_glo(n_tot))
       allocate(load_recv(n_tot))
       allocate(off_recv(n_tot))
       allocate(len_recv(n_tot))

    else

       ! Receive arguments are ignored on non-root ranks.
       rcounts = 0
       displs  = 0

       allocate(gid_glo(0))
       allocate(load_recv(0))
       allocate(off_recv(0))
       allocate(len_recv(0))

    end if


    ! ---------------------------------------------------------------
    ! Gather global-domain IDs
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         gid_loc, sz, MPI_INTEGER, &
         gid_glo, rcounts, displs, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: gid MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Gather domain loads
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         load_loc, sz, MPI_INTEGER, &
         load_recv, rcounts, displs, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: load MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Gather absolute checkpoint offsets
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         off_loc, sz, MPI_INTEGER8, &
         off_recv, rcounts, displs, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: offset MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Gather serialized byte lengths
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         len_loc, sz, MPI_INTEGER8, &
         len_recv, rcounts, displs, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: length MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Construct global-domain ordered directory on rank zero
    ! ---------------------------------------------------------------

    if (rank == 0) then

       allocate(seen(N_GLO_DOMAIN))

       seen = .false.

       do i = 1, n_tot

          gid = gid_glo(i)

          if (gid < 0 .or. gid >= N_GLO_DOMAIN) then

             write(*,'(A,I0)') &
                  "build_checkpoint_directory: invalid gid ", gid

             error stop

          end if


          if (seen(gid+1)) then

             write(*,'(A,I0)') &
                  "build_checkpoint_directory: duplicate gid ", gid

             error stop

          end if


          if (len_recv(i) <= 0_int64) then

             write(*,'(A,I0,A,I0)') &
                  "build_checkpoint_directory: domain ", gid, &
                  " has invalid byte count ", len_recv(i)

             error stop

          end if


          cp_load(gid+1)   = load_recv(i)
          cp_offset(gid+1) = off_recv(i)
          cp_nbytes(gid+1) = len_recv(i)

          seen(gid+1) = .true.

       end do


       if (.not. all(seen)) then
          error stop &
               "build_checkpoint_directory: one or more domains missing"
       end if

       deallocate(seen)

    end if


    ! ---------------------------------------------------------------
    ! Cleanup temporary gather arrays
    ! ---------------------------------------------------------------

    deallocate(gid_glo)
    deallocate(load_recv)
    deallocate(off_recv)
    deallocate(len_recv)

  end subroutine build_checkpoint_directory


  subroutine write_checkpoint_directory (fh)
    implicit none

    type(MPI_File), intent(in) :: fh

    integer :: ierr
    type(MPI_Status) :: status

    integer(int64) :: header(3)

    header(1) = CP_MAGIC
    header(2) = CP_VERSION
    header(3) = int(N_GLO_DOMAIN, int64)

    call MPI_File_write_at( &
         fh, &
         int(0, MPI_OFFSET_KIND), &
         header, 3, MPI_INTEGER8, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: header write failed"
    end if

    call MPI_File_write_at( &
         fh, &
         int(CP_LOAD_POS, MPI_OFFSET_KIND), &
         cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: load write failed"
    end if

    call MPI_File_write_at( &
         fh, &
         int(CP_OFFSET_POS, MPI_OFFSET_KIND), &
         cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: offset write failed"
    end if

    call MPI_File_write_at( &
         fh, &
         int(CP_NBYTES_POS, MPI_OFFSET_KIND), &
         cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: nbytes write failed"
    end if

  end subroutine write_checkpoint_directory


  subroutine read_checkpoint_directory (id)
    implicit none

    integer, intent(in) :: id

    integer :: ierr

    type(MPI_File)   :: fh
    type(MPI_Status) :: status

    integer(int64) :: header(3)

    character(4) :: cp4
    character(:), allocatable :: filename

    write(cp4,'(i4.4)') id

    filename = &
         trim(run_id) // "_checkpoint_" // cp4 // ".bin"

    call MPI_File_open( &
         comm, filename, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "read_checkpoint_directory: MPI_File_open failed"
    end if

    if (allocated(cp_load))   deallocate(cp_load)
    if (allocated(cp_offset)) deallocate(cp_offset)
    if (allocated(cp_nbytes)) deallocate(cp_nbytes)

    allocate(cp_load(N_GLO_DOMAIN))
    allocate(cp_offset(N_GLO_DOMAIN))
    allocate(cp_nbytes(N_GLO_DOMAIN))

    if (rank == 0) then

       call MPI_File_read_at( &
            fh, &
            int(0, MPI_OFFSET_KIND), &
            header, 3, MPI_INTEGER8, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: header read failed"
       end if

       if (header(1) /= CP_MAGIC) then
          error stop "read_checkpoint_directory: bad checkpoint magic"
       end if

       if (header(2) /= CP_VERSION) then
          error stop "read_checkpoint_directory: unsupported version"
       end if

       if (header(3) /= int(N_GLO_DOMAIN, int64)) then
          error stop "read_checkpoint_directory: wrong domain count"
       end if

       call MPI_File_read_at( &
            fh, &
            int(CP_LOAD_POS, MPI_OFFSET_KIND), &
            cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: load read failed"
       end if

       call MPI_File_read_at( &
            fh, &
            int(CP_OFFSET_POS, MPI_OFFSET_KIND), &
            cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: offset read failed"
       end if

       call MPI_File_read_at( &
            fh, &
            int(CP_NBYTES_POS, MPI_OFFSET_KIND), &
            cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: nbytes read failed"
       end if

    end if

    call MPI_File_close(fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "read_checkpoint_directory: MPI_File_close failed"
    end if

    call MPI_Bcast( &
         cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
         0, comm, ierr)

    call MPI_Bcast( &
         cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    call MPI_Bcast( &
         cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

  end subroutine read_checkpoint_directory


  subroutine dump_domain (fid, d)
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: d

    integer :: c, ibeg, iend, j, k
    integer :: l, p_chd, p_lev, p_par, v

    logical, dimension(1:N_CHDRN) :: required

    ! Checkpoint-global/model state.

    write(fid) istep_cumul
    write(fid) time
    write(fid) itime
    write(fid) iwrite
    write(fid) threshold

    call dump(fid)

    ! Coarsest-scale solution.
    call apply_to_pole_d ( &
         write_scalar, grid(d), min_level-1, z_null, fid, .true.)

    p_par = 1

    do k = zmin, zmax
       do v = 1, N_VARIABLE

          ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
          iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1

          write(fid) sol(v,k)%data(d)%elts(ibeg:iend)

       end do
    end do

    if (vert_diffuse) then

       do k = 1, zlevels

          ibeg = grid(d)%patch%elts(p_par+1)%elts_start + 1
          iend = ibeg + PATCH_SIZE**2 - 1

          write(fid) tke(k)%data(d)%elts(ibeg:iend)

       end do

    end if

    ! Finer scales.
    do l = min_level, level_end

       p_lev = 0

       do j = 1, grid(d)%lev(l)%length

          p_par = grid(d)%lev(l)%elts(j)

          ! Do not save children of a deleted parent.
          if (grid(d)%patch%elts(p_par+1)%deleted) then

             do c = 1, N_CHDRN

                p_chd = grid(d)%patch%elts(p_par+1)%children(c)

                if (p_chd > 0) then
                   grid(d)%patch%elts(p_chd+1)%deleted = .true.
                end if

             end do

             cycle

          end if

          ! Wavelet coefficients.
          do k = zmin, zmax
             do v = 1, N_VARIABLE

                ibeg = &
                     MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1

                iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1

                write(fid) &
                     wav_coeff(v,k)%data(d)%elts(ibeg:iend)

             end do
          end do

          if (vert_diffuse) then

             do k = 1, zlevels

                ibeg = &
                     grid(d)%patch%elts(p_par+1)%elts_start + 1

                iend = ibeg + PATCH_SIZE**2 - 1

                write(fid) wav_tke(k)%data(d)%elts(ibeg:iend)

             end do

          end if

          ! Record which children are required.
          do c = 1, N_CHDRN

             p_chd = grid(d)%patch%elts(p_par+1)%children(c)

             if (p_chd > 0) then

                required(c) = &
                     check_child_required(grid(d), p_par, c-1)

                grid(d)%patch%elts(p_chd+1)%deleted = &
                     .not. required(c)

                if (required(c) .and. &
                     p_lev+1 <= size(grid(d)%lev(l+1)%elts)) then

                   p_lev = p_lev + 1
                   grid(d)%lev(l+1)%elts(p_lev) = p_chd

                end if

             else

                required(c) = .false.

             end if

          end do

          ! This used to be written to fid_gr.
          ! It now immediately follows this patch's field data.
          write(fid) required

       end do

       if (l+1 <= max_level) then
          grid(d)%lev(l+1)%length = p_lev
       end if

    end do

  end subroutine dump_domain


  subroutine write_scalar (dom, p, i, j, zlev, offs, dims, fid)
    ! For poles

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: fid, i, j, p, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, k, v

    d = dom%id+1
    id = idx(i, j, offs, dims) + 1

    do k = zmin, zmax
       do v = scalars(1), scalars(2)
          write (fid) sol(v,k)%data(d)%elts(id)
       end do
    end do

    if (vert_diffuse) then
       do k = 1, zlevels
          write (fid) tke(k)%data(d)%elts(id)
       end do
    end if
  end subroutine write_scalar


  subroutine load_adapt_mpi (id)
    ! Read adaptive checkpoint data from one shared MPI-IO file.

    implicit none

    integer, intent(in) :: id

    integer :: d, gid
    integer :: fid
    integer :: ierr

    integer :: cmdstat, exitstat
    integer :: decompress_status
    integer :: delete_unit

    type(MPI_File)   :: fh
    type(MPI_Status) :: status

    integer(int64) :: nbytes
    integer(int64) :: p

    integer(int64), allocatable :: domain_pos(:)
    integer(int8),  allocatable :: buf(:)

    logical :: bin_exists
    logical :: zst_exists
    logical :: temporary_bin

    character(4) :: cp4
    character(:), allocatable :: filename
    character(:), allocatable :: zst_filename
    character(:), allocatable :: cmd


    ! ---------------------------------------------------------------
    ! Validate checkpoint directory
    ! ---------------------------------------------------------------

    if (.not. allocated(cp_offset)) then
       error stop "load_adapt_mpi: cp_offset is not allocated"
    end if

    if (.not. allocated(cp_nbytes)) then
       error stop "load_adapt_mpi: cp_nbytes is not allocated"
    end if

    if (size(cp_offset) /= N_GLO_DOMAIN .or. &
         size(cp_nbytes) /= N_GLO_DOMAIN) then

       error stop "load_adapt_mpi: invalid checkpoint directory"

    end if


    ! ---------------------------------------------------------------
    ! Construct checkpoint filenames
    ! ---------------------------------------------------------------

    write(cp4,'(i4.4)') id

    filename = &
         trim(run_id) // "_checkpoint_" // cp4 // ".bin"

    zst_filename = filename // ".zst"


    ! ---------------------------------------------------------------
    ! Decompress checkpoint if necessary
    ! ---------------------------------------------------------------
    !
    ! If filename already exists, use it directly.
    !
    ! Otherwise rank zero decompresses filename.zst while retaining
    ! the compressed checkpoint.  The resulting .bin is considered
    ! temporary and is deleted after a successful load.

    temporary_bin    = .false.
    decompress_status = 0


    if (rank == 0) then

       inquire(file=filename,     exist=bin_exists)
       inquire(file=zst_filename, exist=zst_exists)


       if (.not. bin_exists) then

          if (.not. zst_exists) then

             decompress_status = 1

          else

             cmd = 'zstd -q -d -k -f "' // trim(zst_filename) // '"'

             call execute_command_line( &
                  cmd, &
                  exitstat=exitstat, &
                  cmdstat=cmdstat)

             if (cmdstat /= 0) then

                decompress_status = 2

             elseif (exitstat /= 0) then

                decompress_status = 3

             else

                ! Verify that decompression actually produced the file.
                inquire(file=filename, exist=bin_exists)

                if (.not. bin_exists) then
                   decompress_status = 4
                else
                   temporary_bin = .true.
                end if

             end if

          end if

       end if

    end if


    ! All ranks need to know whether the decompression succeeded.

    call MPI_Bcast( &
         decompress_status, 1, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: decompression status broadcast failed"
    end if


    if (decompress_status /= 0) then

       select case (decompress_status)

       case (1)
          error stop "load_adapt_mpi: checkpoint .bin and .bin.zst not found"

       case (2)
          error stop "load_adapt_mpi: failed to execute zstd"

       case (3)
          error stop "load_adapt_mpi: zstd decompression failed"

       case (4)
          error stop "load_adapt_mpi: decompressed checkpoint not found"

       case default
          error stop "load_adapt_mpi: checkpoint decompression failed"

       end select

    end if


    ! Broadcast whether the .bin file was created temporarily.
    !
    ! Only rank zero actually needs this information for deletion,
    ! but broadcasting keeps the state consistent on all ranks.

    call MPI_Bcast( &
         temporary_bin, 1, MPI_LOGICAL, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: temporary_bin broadcast failed"
    end if


    ! Ensure that decompression has completed and the file is visible
    ! before any rank attempts MPI_File_open.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: pre-open barrier failed"
    end if


    ! ---------------------------------------------------------------
    ! Open shared checkpoint file
    ! ---------------------------------------------------------------

    call MPI_File_open( &
         comm, filename, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: MPI_File_open failed"
    end if


    ! ---------------------------------------------------------------
    ! Create rank-local scratch stream
    ! ---------------------------------------------------------------

    open( &
         newunit=fid, &
         status="scratch", &
         form="unformatted", &
         access="stream", &
         action="readwrite")


    allocate(domain_pos(size(grid)))


    ! ---------------------------------------------------------------
    ! Read locally owned domains from shared checkpoint
    ! ---------------------------------------------------------------

    do d = 1, size(grid)

       gid = glo_id(rank+1,d)


       if (gid < 0 .or. gid >= N_GLO_DOMAIN) then
          error stop "load_adapt_mpi: invalid global domain ID"
       end if


       nbytes = cp_nbytes(gid+1)


       if (nbytes <= 0_int64) then
          error stop "load_adapt_mpi: invalid domain record size"
       end if


       ! buf() has default-integer bounds, so guard the conversion.
       if (nbytes > int(huge(0),int64)) then
          error stop "load_adapt_mpi: domain record too large"
       end if


       allocate(buf(int(nbytes)))


       call MPI_File_read_at( &
            fh, &
            int(cp_offset(gid+1), MPI_OFFSET_KIND), &
            buf, int(nbytes), MPI_BYTE, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "load_adapt_mpi: domain read failed"
       end if


       inquire(unit=fid, pos=p)

       domain_pos(d) = p


       write(fid) buf


       deallocate(buf)

    end do


    ! ---------------------------------------------------------------
    ! Close shared checkpoint
    ! ---------------------------------------------------------------

    call MPI_File_close(fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: MPI_File_close failed"
    end if


    ! ---------------------------------------------------------------
    ! Reconstruct local domains
    ! ---------------------------------------------------------------

    flush(fid)

    call load_domains_stream(fid, domain_pos)

    close(fid)

    deallocate(domain_pos)


    ! ---------------------------------------------------------------
    ! Remove temporary uncompressed checkpoint
    ! ---------------------------------------------------------------
    !
    ! All ranks must have completely finished reading and reconstructing
    ! their domains before rank zero removes the shared .bin file.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: pre-delete barrier failed"
    end if


    if (rank == 0 .and. temporary_bin) then

       ! Delete using Fortran rather than invoking an external rm command.

       open( &
            newunit=delete_unit, &
            file=filename, &
            status='old', &
            action='read', &
            iostat=ierr)

       if (ierr /= 0) then
          error stop "load_adapt_mpi: unable to open temporary checkpoint for deletion"
       end if

       close(delete_unit, status='delete', iostat=ierr)

       if (ierr /= 0) then
          error stop "load_adapt_mpi: unable to delete temporary checkpoint"
       end if

    end if


    ! Ensure deletion is complete before returning.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: post-delete barrier failed"
    end if


  end subroutine load_adapt_mpi

  subroutine load_domains_stream (fid, domain_pos)
    ! Reconstruct all locally owned domains from a single
    ! rank-local unformatted stream.
    !
    ! domain_pos(d) is the current position in the serialized
    ! record for local domain d.

    use iso_fortran_env, only: int64

    implicit none

    integer,        intent(in)    :: fid
    integer(int64), intent(inout) :: domain_pos(:)

    integer :: c, d
    integer :: ibeg, iend
    integer :: j, k, l
    integer :: old_n_patch
    integer :: p_chd, p_par
    integer :: v

    logical :: first_read

    logical, dimension(N_CHDRN) :: required

    if (size(domain_pos) /= size(grid)) then
       error stop "load_domains_stream: domain_pos has wrong size"
    end if

    !
    ! ------------------------------------------------------------
    ! Coarsest-scale solution
    ! ------------------------------------------------------------
    !
    do d = 1, size(grid)

       !
       ! The first read explicitly seeks to this domain's record.
       ! All subsequent reads for this domain are sequential.
       !
       read(fid, pos=domain_pos(d)) istep_cumul

       read(fid) time
       read(fid) itime
       read(fid) iwrite
       read(fid) threshold

       !
       ! Test-case-specific data. This may be empty.
       !
       call load(fid)

       !
       ! Pole scalar values.
       !
       call apply_to_pole_d ( &
            read_scalar, &
            grid(d), &
            min_level-1, &
            z_null, &
            fid, &
            .true.)

       !
       ! Coarse parent patch.
       !
       p_par = 1

       do k = zmin, zmax
          do v = 1, N_VARIABLE

             ibeg = &
                  MULT(v) * &
                  grid(d)%patch%elts(p_par+1)%elts_start + 1

             iend = &
                  ibeg + MULT(v)*PATCH_SIZE**2 - 1

             read(fid) &
                  sol(v,k)%data(d)%elts(ibeg:iend)

          end do
       end do

       if (vert_diffuse) then

          do k = 1, zlevels

             ibeg = &
                  grid(d)%patch%elts(p_par+1)%elts_start + 1

             iend = ibeg + PATCH_SIZE**2 - 1

             read(fid) &
                  tke(k)%data(d)%elts(ibeg:iend)

          end do

       end if

       !
       ! Save the position of this domain's first finer-scale
       ! wavelet record.
       !
       inquire(unit=fid, pos=domain_pos(d))

    end do

    !
    ! ------------------------------------------------------------
    ! Finer-scale wavelets
    ! ------------------------------------------------------------
    !
    ! Preserve the original level-by-level reconstruction:
    !
    !     process all domains at level l
    !     call post_refine
    !     proceed to newly created level
    !
    l = 1

    do while (level_end > l)

       l = level_end

       if (rank == 0) then
          write(6,'(a,i2)') 'Loading level ', l
       end if

       do d = 1, size(grid)

          old_n_patch = grid(d)%patch%length

          !
          ! The first actual read for this domain at this level
          ! must seek back to domain_pos(d). Thereafter all reads
          ! are sequential.
          !
          first_read = .true.

          do j = 1, grid(d)%lev(l)%length

             p_par = grid(d)%lev(l)%elts(j)

             !
             ! Wavelet coefficients.
             !
             do k = zmin, zmax
                do v = 1, N_VARIABLE

                   ibeg = &
                        MULT(v) * &
                        grid(d)%patch%elts(p_par+1)%elts_start + 1

                   iend = &
                        ibeg + MULT(v)*PATCH_SIZE**2 - 1

                   if (first_read) then

                      read(fid, pos=domain_pos(d)) &
                           wav_coeff(v,k)%data(d)%elts(ibeg:iend)

                      first_read = .false.

                   else

                      read(fid) &
                           wav_coeff(v,k)%data(d)%elts(ibeg:iend)

                   end if

                end do
             end do

             !
             ! Turbulent kinetic energy wavelets.
             !
             if (vert_diffuse) then

                do k = 1, zlevels

                   ibeg = &
                        grid(d)%patch%elts(p_par+1)%elts_start + 1

                   iend = ibeg + PATCH_SIZE**2 - 1

                   !
                   ! Normally first_read is already false because
                   ! N_VARIABLE > 0. This branch keeps the stream
                   ! handling correct even if that ever changes.
                   !
                   if (first_read) then

                      read(fid, pos=domain_pos(d)) &
                           wav_tke(k)%data(d)%elts(ibeg:iend)

                      first_read = .false.

                   else

                      read(fid) &
                           wav_tke(k)%data(d)%elts(ibeg:iend)

                   end if

                end do

             end if

             !
             ! Child-presence information follows this patch's
             ! coefficient data in the stage-2/3 stream.
             !
             if (first_read) then

                read(fid, pos=domain_pos(d)) required
                first_read = .false.

             else

                read(fid) required

             end if

             do c = 1, N_CHDRN

                if (required(c)) then
                   call refine_patch1(grid(d), p_par, c-1)
                end if

             end do

          end do

          !
          ! Only update the domain position if something was
          ! actually read at this level.
          !
          ! This is important for a domain with zero patches at
          ! the current level: its saved stream position must not
          ! accidentally become the position of another domain.
          !
          if (.not. first_read) then
             inquire(unit=fid, pos=domain_pos(d))
          end if

          !
          ! Complete child patch construction exactly as in the
          ! original loader.
          !
          do p_par = 2, old_n_patch

             do c = 1, N_CHDRN

                p_chd = &
                     grid(d)%patch%elts(p_par)%children(c)

                if (p_chd+1 > old_n_patch) then

                   call refine_patch2( &
                        grid(d), &
                        p_par-1, &
                        c-1)

                end if

             end do

          end do

       end do

       !
       ! Must remain outside the domain loop.
       !
       call post_refine

    end do

    !
    ! ------------------------------------------------------------
    ! Boundary state
    ! ------------------------------------------------------------
    !
    sol%bdry_uptodate       = .false.
    wav_coeff%bdry_uptodate = .false.

    if (vert_diffuse) then
       tke%bdry_uptodate     = .false.
       wav_tke%bdry_uptodate = .false.
    end if

  end subroutine load_domains_stream


  subroutine read_scalar (dom, p, i, j, zlev, offs, dims, fid)
    ! For poles

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: fid, i, j, p, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, k, v

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    do k = zmin, zmax
       do v = scalars(1), scalars(2)
          read (fid) sol(v,k)%data(d)%elts(id)
       end do
    end do

    if (vert_diffuse) then
       do k = 1, zlevels
          read (fid) tke(k)%data(d)%elts(id)
       end do
    end if
  end subroutine read_scalar


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
    character(9999    )       :: filename
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

  
end module io_mod
