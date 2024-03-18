module io_mod
  use domain_mod
  use ops_mod
  use smooth_mod
  use adapt_mod
  use utils_mod
  implicit none

  integer, dimension(2,4)                       :: HR_offs
  data                                             HR_offs /0,0, 1,0, 1,1, 0,1/
  
  integer                                       :: ii, N_VAR_OUT
  integer, dimension(:,:), allocatable          :: topo_count
  real(8)                                       :: vmin, vmax
  real(8), dimension(:), allocatable            :: minv, maxv
  integer                                       :: next_fid
  type(Float_Field)                             :: active_level
  type(Topo_Array),  dimension(:,:), allocatable :: topography_data
contains
  subroutine init_io_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_domain_mod
    next_fid = 100

    if (compressible .or. zlevels /= 2) then
       N_VAR_OUT = 7
    else
       N_VAR_OUT = 11
    end if

    allocate (minv(1:N_VAR_OUT), maxv(1:N_VAR_OUT))
    minv = 0d0
    maxv = 0d0
    initialized = .true.
  end subroutine init_io_mod

  integer function get_fid ()
    implicit none

    get_fid  = next_fid
    next_fid = next_fid + 1
  end function get_fid

  subroutine vort_extrema (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idE
    real(8) :: vort

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

  subroutine sum_total_mass (initialgo)
    ! Total mass over all vertical layers
    implicit none
    logical :: initialgo

    integer :: k

    if (initialgo) then
       initotalmass = 0d0
       do k = zmin, zmax
          initotalmass = initotalmass + integrate_hex (full_mass, k, level_start)
       end do
    else
       totalmass = 0d0
       do k = zmin, zmax
          totalmass = totalmass + integrate_hex (full_mass, k, level_start)
       end do
       mass_error = abs (totalmass - initotalmass) / initotalmass
    end if
  end subroutine sum_total_mass

  real(8) function full_mass (dom, i, j, zlev, offs, dims)
    use domain_mod
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    full_mass = sol(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(id+1)
  end function full_mass

  real(8) function topo (dom, i, j, zlev, offs, dims)
    use domain_mod
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    topo = topography%data(d)%elts(id)
  end function topo

  real(8) function pot_energy (dom, i, j, zlev, offs, dims)
    use domain_mod
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: full_mass

    id = idx (i, j, offs, dims)

    full_mass = sol_mean(S_MASS,zlev)%data(dom%id+1)%elts(id+1) + sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)
    pot_energy = full_mass**2
  end function pot_energy

  real(8) function total_ke (itype)
    ! Computes total kinetic energy
    implicit none
    character(*) :: itype

    integer :: k
    
    total_ke = 0d0
    do k = 1, zlevels
       total_ke = total_ke + integrate_hex (kinetic_energy, k)
    end do
  end function total_ke

  subroutine umag (q)
    ! Evaluate complete velocity trend by adding gradient terms to previously calculated source terms on entire grid
    ! Input is velocity field q at a single vertical layer
    implicit none
    type(Float_Field), target :: q
    integer                   :: k

    integer :: d, j, p
    
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
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
   
    integer :: d, id, idW, idSW, idS
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S
    
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
                      * dom%areas%elts(id+1)%hex_inv/(4d0*MATH_PI*radius**2)) 
  end subroutine cal_umag

  real(8) function kinetic_energy (dom, i, j, zlev, offs, dims)
    ! Kinetic energy u^2/2 at level zlev using approximation to TRiSK formula
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
   
    integer :: d, id, idW, idSW, idS
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S
    
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

    kinetic_energy = (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT + &
                      u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
         * dom%areas%elts(id+1)%hex_inv/4d0 
  end function kinetic_energy

  real(8) function layer1_ke (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idW, idSW, idS
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

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

       layer1_ke = (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1))  &
            * (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       layer1_ke = 0d0
    end if
  end function layer1_ke

  real(8) function layer2_ke (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idW, idSW, idS
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

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

       layer2_ke = (sol_mean(S_MASS,2)%data(d)%elts(id+1) + sol(S_MASS,2)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       layer2_ke = 0d0
    end if
  end function layer2_ke

  real(8) function one_layer_ke (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idW, idSW, idS
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

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

       one_layer_ke = (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       one_layer_ke = 0d0
    end if
  end function one_layer_ke

  real(8) function barotropic_ke (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, id, idE, idNE, idN, idW, idSW, idS
    integer, dimension(1:EDGE) :: id_edge, id_node
    real(8)                    :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(8)                    :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S
    real(8), dimension(1:EDGE) :: u

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       id_node = (/ idE, idNE, idN /)
       id_edge = (/ id,  id,   id /)
       u = barotropic_velo ()
       u_prim_RT = u(1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = u(2) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = u(3) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = u(1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = u(2) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = u(3) * dom%pedlen%elts(EDGE*id+UP+1)

       id_node = (/ idW, idSW, idS /)
       id_edge = id_node
       u = barotropic_velo ()
       u_prim_RT_W  = u(1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_DG_SW = u(2) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_S  = u(3) * dom%len%elts(EDGE*idS +UP+1)

       u_dual_RT_W  = u(1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = u(2) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = u(3) * dom%pedlen%elts(EDGE*idS +UP+1)

       barotropic_ke = &
            (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1) + &
            sol_mean(S_MASS,2)%data(d)%elts(id+1) + sol(S_MASS,2)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       barotropic_ke = 0d0
    end if
  contains
    function barotropic_velo ()
      real(8), dimension(1:EDGE) :: barotropic_velo

      integer                    :: e, id_e, k
      real(8), dimension(1:EDGE) :: dz

      do e = 1, EDGE
         id_e = EDGE*id_edge(e) + e

         do k = 1, 2
            dz(k) = interp (sol_mean(S_MASS,k)%data(d)%elts(id+1)         + sol(S_MASS,k)%data(d)%elts(id+1), &
                 sol_mean(S_MASS,k)%data(d)%elts(id_node(e)+1) + sol(S_MASS,k)%data(d)%elts(id_node(e)+1))
         end do

         barotropic_velo(e) = (dz(1)*sol(S_VELO,1)%data(d)%elts(id_e) + dz(2)*sol(S_VELO,2)%data(d)%elts(id_e)) / sum(dz)
      end do
    end function barotropic_velo
  end function barotropic_ke

  real(8) function pot_enstrophy (dom, i, j, zlev, offs, dims)
    ! Computes potential enstrophy in two layer mode split case
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i
    real(8) :: f, h, w

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
       pot_enstrophy = 0.5d0 * (w / h)**2
    else
       pot_enstrophy = 0d0
    end if
  end function pot_enstrophy

   subroutine cal_temp (dom, i, j, zlev, offs, dims)
    ! Compute temperature in compressible case
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: id, d, k
    real(8)                       :: full_mass, full_mass_lower, full_temp
    real(8), dimension(1:zlevels) :: p

    d = dom%id + 1
    id = idx(i, j, offs, dims) + 1

    ! Integrate the pressure upwards
    full_mass = sol_mean(S_MASS,1)%data(d)%elts(id) + sol(S_MASS,1)%data(d)%elts(id)
    p(1) = dom%surf_press%elts(id) - 0.5d0 * grav_accel * full_mass
    do k = 2, zlevels
       full_mass_lower = full_mass
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       p(k) = p(k-1) - grav_accel * interp (full_mass, full_mass_lower)
    end do

    ! Temperature at all vertical levels (saved in exner_fun)
    do k = 1, zlevels
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       full_temp = sol_mean(S_TEMP,k)%data(d)%elts(id) + sol(S_TEMP,k)%data(d)%elts(id)
       
       exner_fun(k)%data(d)%elts(id) = full_temp / full_mass * (p(k)/p_0)**kappa
    end do

    ! Temperature at save levels (saved in trend)
    do k = 1, save_levels
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id) + sol_save(S_MASS,k)%data(d)%elts(id)
       full_temp = sol_mean(S_TEMP,k)%data(d)%elts(id) + sol_save(S_TEMP,k)%data(d)%elts(id)
       
       trend(1,k)%data(d)%elts(id+1) = full_temp / full_mass * (pressure_save(k) / p_0)**kappa
    end do
  end subroutine cal_temp

  subroutine cal_geopot (dom, i, j, zlev, offs, dims)
    ! Compute geopotential in compressible case
    ! Assumes that temperature has already been calculated and stored in exner_fun
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, d, k
    real(8) :: full_mass, pressure_lower, pressure_upper

    d = dom%id + 1
    id = idx(i, j, offs, dims) + 1

    ! Integrate geopotential upwards from surface
    full_mass = sol_mean(S_MASS,1)%data(d)%elts(id) + sol(S_MASS,1)%data(d)%elts(id)

    pressure_lower = dom%surf_press%elts(id)
    pressure_upper = pressure_lower - grav_accel * full_mass
    
    dom%geopot_lower%elts(id) = surf_geopot (d, id) / grav_accel

    k = 1
    do while (pressure_upper > pressure_save(1))
       dom%geopot_lower%elts(id) = dom%geopot_lower%elts(id) + &
            R_d/grav_accel * exner_fun(k)%data(d)%elts(id) * (log(pressure_lower)-log(pressure_upper))

       k = k+1
       full_mass = sol_mean(S_MASS,k+1)%data(d)%elts(id) + sol(S_MASS,k+1)%data(d)%elts(id)

       pressure_lower = pressure_upper
       pressure_upper = pressure_lower - grav_accel * full_mass
    end do

    ! Add additional contribution up to pressure level pressure_save(1)
    dom%geopot_lower%elts(id) = dom%geopot_lower%elts(id) &
         + R_d/grav_accel * exner_fun(k)%data(d)%elts(id) * (log(pressure_lower) - log(pressure_save(1)))
  end subroutine cal_geopot

  subroutine statistics
    ! Calculates zonal statistics
    use domain_mod
    implicit none
    integer :: d, k, p

    if (rank == 0) write (6,'(a)') 'Incrementing zonal averages'

    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))

    do k = 1, zlevels
       do d = 1, size (grid)
          mass   => sol(S_MASS,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          temp   => sol(S_TEMP,k)%data(d)%elts
          velo   => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), p-1, k, 0, 0)
             call apply_onescale_to_patch (cal_pressure,     grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (cal_zonal_avg,    grid(d), p-1, k, 0, 0)
          end do
          nullify (mass, mean_m, temp, velo, velo1, velo2)
       end do
    end do
  end subroutine statistics

  subroutine cal_zonal_avg (dom, i, j, zlev, offs, dims)
    ! Zonal average means and covariances over all checkpoints using stable online algorithm
    ! Uses Welford's stable onlne algorithm
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: bin, d, id_i
    real(8) :: lat, lon, temperature, Tprime, Uprime, Vprime, Tprime_new, Uprime_new, Vprime_new

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    call cart2sph (dom%node%elts(id_i), lon, lat)

    temperature = (temp(id_i)/mass(id_i)) * (dom%press%elts(id_i)/p_0)**kappa

    do bin = 1, nbins
       if (lat*1.8d2/MATH_PI <= bounds(bin)) then
          Nstats(zlev,bin) = Nstats(zlev,bin) + 1

          Tprime = temperature            - zonal_avg(zlev,bin,1)
          Uprime = dom%u_zonal%elts(id_i) - zonal_avg(zlev,bin,3)
          Vprime = dom%v_merid%elts(id_i) - zonal_avg(zlev,bin,4)

          ! Mean values
          zonal_avg(zlev,bin,1) = zonal_avg(zlev,bin,1) + Tprime/Nstats(zlev,bin)
          zonal_avg(zlev,bin,3) = zonal_avg(zlev,bin,3) + Uprime/Nstats(zlev,bin)
          zonal_avg(zlev,bin,4) = zonal_avg(zlev,bin,4) + Vprime/Nstats(zlev,bin)
          zonal_avg(zlev,bin,5) = zonal_avg(zlev,bin,5) +  0.5 * (Uprime**2 + Vprime**2)/Nstats(zlev,bin)

          Tprime_new = temperature            - zonal_avg(zlev,bin,1)
          Uprime_new = dom%u_zonal%elts(id_i) - zonal_avg(zlev,bin,3)
          Vprime_new = dom%v_merid%elts(id_i) - zonal_avg(zlev,bin,4)

          ! Update sums of squares (for variance calculation)

          ! Temperature 
          zonal_avg(zlev,bin,2) = zonal_avg(zlev,bin,2) + Tprime * Tprime_new

          ! Eddy momentum flux (covariance)
          zonal_avg(zlev,bin,6) = zonal_avg(zlev,bin,6) + Uprime * Vprime_new

          ! Zonal velocity variance
          zonal_avg(zlev,bin,7) = zonal_avg(zlev,bin,7) + Uprime * Uprime_new 

          ! Meridional velocity variance
          zonal_avg(zlev,bin,8) = zonal_avg(zlev,bin,8) + Vprime * Vprime_new 

          ! Eddy heat flux (covariance)
          zonal_avg(zlev,bin,9) = zonal_avg(zlev,bin,9) + Vprime * Tprime_new

          exit
       end if
    end do
  end subroutine cal_zonal_avg

  subroutine combine_stats
    ! Uses Chan, Golub and LeVeque (1983) algorithm for partitioned data sets to combine zonal average results from each rank
    !   T.F. Chan, G.H. Golub & R.J. LeVeque (1983):
    !   "Algorithms for computing the sample variance: Analysis and recommendations." The American Statistician 37: 242–247.
#ifdef MPI    
    use mpi
#endif
    implicit none
    integer                                  :: bin, k
    real(8), dimension(nvar_zonal)           :: temp
    integer, dimension(n_process)            :: Nstats_loc
    real(8), dimension(n_process*nvar_zonal) :: zonal_avg_loc

#ifdef MPI
    ! Collect statistics data on rank 0 for combining
    Nstats_glo    = 0d0
    zonal_avg_glo = 0d0
    do k = 1, zlevels
       do bin = 1, nbins
          call MPI_Gather (Nstats(k,bin), 1, MPI_INTEGER, Nstats_loc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

          temp = zonal_avg(k,bin,:)
          call MPI_Gather (temp, nvar_zonal, MPI_DOUBLE_PRECISION, zonal_avg_loc, nvar_zonal, MPI_DOUBLE_PRECISION, &
               0, MPI_COMM_WORLD, ierror)
          
          if (rank == 0) call combine_var
       end do
    end do
#else
    Nstats_glo    = Nstats
    zonal_avg_glo = zonal_avg
#endif
  contains
    subroutine combine_var
      integer :: m, nA, nB, nAB, r
      real(8) :: delta_KE, delta_T, delta_U, delta_V
      real(8) :: mA_T, mB_T, mA_U, mB_U, mA_V, mB_V, mA_zonal, mB_zonal, mA_merid, mB_merid, mA_VT, mB_VT, mA_UV, mB_UV

      do r = 0, n_process-1
         nA = Nstats_glo(k,bin)
         nB = Nstats_loc(r+1)
         if (nB /= 0) then
            m = r*nvar_zonal
            nAB = nA + nB

            delta_T  = zonal_avg_loc(r*nvar_zonal+1) - zonal_avg_glo(k,bin,1)
            delta_U  = zonal_avg_loc(r*nvar_zonal+3) - zonal_avg_glo(k,bin,3)
            delta_V  = zonal_avg_loc(r*nvar_zonal+4) - zonal_avg_glo(k,bin,4)
            delta_KE = zonal_avg_loc(r*nvar_zonal+5) - zonal_avg_glo(k,bin,5)

            mA_T   = zonal_avg_glo(k,bin,2)
            mB_T   = zonal_avg_loc(m+2)

            mA_UV  = zonal_avg_glo(k,bin,6)
            mB_UV  = zonal_avg_loc(m+6)

            mA_zonal = zonal_avg_glo(k,bin,7)
            mB_zonal = zonal_avg_loc(m+7)

            mA_merid = zonal_avg_glo(k,bin,8)
            mB_merid = zonal_avg_loc(m+8)

            mA_VT  = zonal_avg_glo(k,bin,9)
            mB_VT  = zonal_avg_loc(m+9)

            ! Combine means
            zonal_avg_glo(k,bin,1) = zonal_avg_glo(k,bin,1) + delta_T  * nB/nAB
            zonal_avg_glo(k,bin,3) = zonal_avg_glo(k,bin,3) + delta_U  * nB/nAB
            zonal_avg_glo(k,bin,4) = zonal_avg_glo(k,bin,4) + delta_V  * nB/nAB
            zonal_avg_glo(k,bin,5) = zonal_avg_glo(k,bin,5) + delta_KE * nB/nAB

            ! Combine sums of squares (for variances)
            zonal_avg_glo(k,bin,2) = mA_T     + mB_T     + delta_T**2      * nA*nB/nAB ! temperature variance
            zonal_avg_glo(k,bin,6) = mA_UV    + mB_UV    + delta_U*delta_V * nA*nB/nAB ! velocity covariance
            zonal_avg_glo(k,bin,7) = mA_zonal + mB_zonal + delta_U**2      * nA*nB/nAB ! zonal wind variance
            zonal_avg_glo(k,bin,8) = mA_merid + mB_merid + delta_V**2      * nA*nB/nAB ! meridional wind variance
            zonal_avg_glo(k,bin,9) = mA_VT    + mB_VT    + delta_V*delta_T * nA*nB/nAB ! V-T covariance (eddy heat flux)

            ! Update total number of data points
            Nstats_glo(k,bin) = nAB
         end if
      end do
    end subroutine combine_var
  end subroutine combine_stats

  subroutine write_out_stats
    ! Writes out zonal average statistics
    implicit none
    integer            :: ibin, info, k, v
    integer, parameter :: funit = 400
    character(2)       :: var_file
    character(1300)     :: bash_cmd, command

    write (6,'(/,a)') 'Saving statistics'

    ! Find sample covariances from sums of squares
    zonal_avg_glo(:,:,2) = zonal_avg_glo(:,:,2) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,6) = zonal_avg_glo(:,:,6) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,7) = zonal_avg_glo(:,:,7) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,8) = zonal_avg_glo(:,:,8) / (Nstats_glo - 1)
    zonal_avg_glo(:,:,9) = zonal_avg_glo(:,:,9) / (Nstats_glo - 1)

    ! Save number of data points
    write (var_file, '(i2.2)') 00
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="UNFORMATTED", action='WRITE', status='REPLACE')
    write (funit) Nstats_glo
    close (funit)

    ! Save zonal statistics (means and covariances)
    do v = 1, nvar_zonal
       write (var_file, '(i2)') v+10
       open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE', status='REPLACE')
       do k = zlevels,1,-1
          write (funit,'(2047(E15.6, 1X))') zonal_avg_glo(k,:,v)
       end do
       close (funit)
    end do

    ! Longitude values (dummy)
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE', status='REPLACE') 
    write (funit,'(2047(E15.6, 1X))') bounds
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE', status='REPLACE') 
    write (funit,'(2047(E15.6, 1X))') bounds - dbin/2
    close (funit)

    ! Non-dimensional pressure based vertical coordinates p_k/p_s
    Write (var_file, '(i2)') 22
    open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='WRITE', status='REPLACE') 
    write (funit,'(2047(E15.6, 1X))') (0.5*((a_vert(k)+a_vert(k+1))/p_0 + b_vert(k)+b_vert(k+1)), k = zlevels, 1, -1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.3.?? > tmp'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))

    command = 'gtar czf '//trim(run_id)//'.3.tgz -T tmp --remove-files &'
    call system (trim(command), info)
    if (info /= 0) then
       if (rank == 0) write (6,'(a)') "gtar command not present ... aborting"
       call abort
    end if

    command = '\rm -f tmp'
    call system (trim(command))
  end subroutine write_out_stats

  subroutine write_primal (dom, p, i, j, zlev, offs, dims, funit)
    ! Write primal grid for vertical level zlev
    use utils_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, funit
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, idW, idSW, idS, k, outl
    real(8)                       :: full_mass, full_temp, total_depth, total_height
    real(4), dimension(N_VAR_OUT) :: outv

    integer :: idE, idN, idNE
    real(8) :: dz0, dze, r_max
    
    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    
    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       if (allocated(active_level%data)) then ! avoid segfault pre_levelout not used
          outl = nint(active_level%data(d)%elts(id_i))
       else
          outl = 0
       end if

       full_mass = sol(S_MASS,zlev)%data(d)%elts(id_i) + sol_mean(S_MASS,zlev)%data(d)%elts(id_i)
       full_temp = sol(S_TEMP,zlev)%data(d)%elts(id_i) + sol_mean(S_TEMP,zlev)%data(d)%elts(id_i)

       if (compressible .or. zlevels /= 2) then
          if (compressible) then ! temperature in layer zlev
             outv(1) = full_temp / full_mass * (dom%press%elts(id_i)/p_0)**kappa
          else
             outv(1) = ref_density * (1.0_8 - full_temp / full_mass) ! density
          end if

          outv(2) = dom%u_zonal%elts(id_i) * phi_node (d, id_i, zlev)  ! zonal velocity
          outv(3) = dom%v_merid%elts(id_i) * phi_node (d, id_i, zlev)  ! meridional velocity

          if (compressible) then ! topography
             outv(4) = topography%data(d)%elts(id_i)
          else ! vertical velocity
             if (min_level == max_level) then
                outv(4) = trend(S_TEMP,zlev)%data(d)%elts(id_i) 
             else ! average over adjacent hexagons
                outv(4) = ( &
                     trend(S_TEMP,zlev)%data(d)%elts(id_i)   / dom%areas%elts(id_i)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idE+1)  / dom%areas%elts(idE+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idNE+1) / dom%areas%elts(idNE+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idN+1)  / dom%areas%elts(idN+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idW+1)  / dom%areas%elts(idW+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idSW+1) / dom%areas%elts(idSW+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idS+1)  / dom%areas%elts(idS+1)%hex_inv ) / ( &
                     1d0/dom%areas%elts(id_i)%hex_inv + &
                     1d0/dom%areas%elts(idE+1)%hex_inv + 1d0/dom%areas%elts(idNE+1)%hex_inv + 1d0/dom%areas%elts(idN+1)%hex_inv + &
                     1d0/dom%areas%elts(idW+1)%hex_inv + 1d0/dom%areas%elts(idSW+1)%hex_inv + 1d0/dom%areas%elts(idS+1)%hex_inv )
             end if
          end if

          if (compressible .or. .not. penalize) then
              if (min_level == max_level) then
                outv(5) = trend(S_TEMP,zlev)%data(d)%elts(id_i) 
             else ! average over adjacent hexagons
                outv(5) = ( &
                     trend(S_TEMP,zlev)%data(d)%elts(id_i)   / dom%areas%elts(id_i)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idE+1)  / dom%areas%elts(idE+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idNE+1) / dom%areas%elts(idNE+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idN+1)  / dom%areas%elts(idN+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idW+1)  / dom%areas%elts(idW+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idSW+1) / dom%areas%elts(idSW+1)%hex_inv + &
                     trend(S_TEMP,zlev)%data(d)%elts(idS+1)  / dom%areas%elts(idS+1)%hex_inv ) / ( &
                     1d0/dom%areas%elts(id_i)%hex_inv + &
                     1d0/dom%areas%elts(idE+1)%hex_inv + 1d0/dom%areas%elts(idNE+1)%hex_inv + 1d0/dom%areas%elts(idN+1)%hex_inv + &
                     1d0/dom%areas%elts(idW+1)%hex_inv + 1d0/dom%areas%elts(idSW+1)%hex_inv + 1d0/dom%areas%elts(idS+1)%hex_inv )
             end if
          else
             outv(5) = penal_node(zlev)%data(d)%elts(id_i)   ! penalization mask
          end if

          if (compressible) then ! surface pressure
             outv(6) = dom%surf_press%elts(id_i)
          else ! free surface perturbation
             if (mode_split) then 
                outv(6) = sol(S_MASS,zlevels+1)%data(d)%elts(id_i) / phi_node (d, id_i, zlev)
              else
                outv(6) = free_surface (dom, i, j, zlev, offs, dims, sol)
             end if
          end if

          r_max = 0d0
          do k = 1, zlevels
             dz0  = (sol(S_MASS,k)%data(d)%elts(id+1) + sol_mean(S_MASS,k)%data(d)%elts(id+1)) &
                  / porous_density (d, id+1, k)

             dze = (sol(S_MASS,k)%data(d)%elts(idE+1) + sol_mean(S_MASS,k)%data(d)%elts(idE+1)) &
                  / porous_density (d, idE+1, k)
             r_max  = max (r_max, abs (dz0 - dze) / (dz0 + dze))

             dze = (sol(S_MASS,k)%data(d)%elts(idNE+1) + sol_mean(S_MASS,k)%data(d)%elts(idNE+1)) &
                  / porous_density (d, idNE+1, k)
             r_max  = max (r_max, abs (dz0 - dze) / (dz0 + dze))

             dze = (sol(S_MASS,k)%data(d)%elts(idN+1)  + sol_mean(S_MASS,k)%data(d)%elts(idN+1)) &
                  / porous_density(d, idN+1, k)
             r_max  = max (r_max, abs (dz0 - dze) / (dz0 + dze))

             dze = (sol(S_MASS,k)%data(d)%elts(idW+1) + sol_mean(S_MASS,k)%data(d)%elts(idW+1)) &
                  / porous_density (d, idW+1, k)
             r_max  = max (r_max, abs (dz0 - dze) / (dz0 + dze))

             dze = (sol(S_MASS,k)%data(d)%elts(idSW+1) + sol_mean(S_MASS,k)%data(d)%elts(idSW+1)) &
                  / porous_density (d, idSW+1, k)
             r_max  = max (r_max, abs (dz0 - dze) / (dz0 + dze))

             dze = (sol(S_MASS,k)%data(d)%elts(idS+1)  + sol_mean(S_MASS,k)%data(d)%elts(idS+1)) &
                  / porous_density(d, idS+1, k)
             r_max  = max (r_max, abs (dz0 - dze) / (dz0 + dze))
          end do

          outv(7) = dom%ke%elts(id_i)          ! vorticity (at hexagon points)

          write (funit,'(18(E14.5E2, 1X), 7(E14.5E2, 1X), I3, 1X, I3)') &
               dom%ccentre%elts(TRIAG*id  +LORT+1), dom%ccentre%elts(TRIAG*id  +UPLT+1), &
               dom%ccentre%elts(TRIAG*idW +LORT+1), dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
               dom%ccentre%elts(TRIAG*idSW+LORT+1), dom%ccentre%elts(TRIAG*idS +UPLT+1), &
               outv, dom%mask_n%elts(id_i), outl
       else
          ! Density perturbation
          outv(1) = -ref_density * full_temp / full_mass

          ! Barotropic velocity
          outv(2) = dom%u_zonal%elts(id_i)
          outv(3) = dom%v_merid%elts(id_i)

          ! Baroclinic velocity
          outv(4) = trend(S_VELO,1)%data(d)%elts(id_i)
          outv(5) = trend(S_VELO,2)%data(d)%elts(id_i)

          ! Baroclinic vorticity at hexagon points
          outv(6) = dom%geopot_lower%elts(id_i)

          ! Barotropic vorticity at hexagon points
          outv(7) = dom%ke%elts(id_i)

          ! Baroclinic eta
          full_mass = sol_mean(S_MASS,1)%data(d)%elts(id_i) + sol(S_MASS,1)%data(d)%elts(id_i)
          outv(8) = full_mass / (ref_density * phi_node (d, id_i, zlev))

          ! Free surface perturbation (barotropic eta)
          if (mode_split) then 
             outv(9) = sol(S_MASS,zlevels+1)%data(d)%elts(id_i) / phi_node (d, id_i, zlev)
          else
             outv(9) = free_surface (dom, i, j, zlev, offs, dims, sol)
          end if

          ! Topography
          outv(10) = -dom%topo%elts(id_i)

          ! Penalization mask
          outv(11) = penal_node(zlev)%data(d)%elts(id_i)  ! penalization mask

          write (funit,'(18(E14.5E2, 1X), 11(E14.5E2, 1X), I3, 1X, I3)') &
               dom%ccentre%elts(TRIAG*id  +LORT+1), dom%ccentre%elts(TRIAG*id  +UPLT+1), &
               dom%ccentre%elts(TRIAG*idW +LORT+1), dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
               dom%ccentre%elts(TRIAG*idSW+LORT+1), dom%ccentre%elts(TRIAG*idS +UPLT+1), &
               outv, dom%mask_n%elts(id_i), outl
       end if
       where (minv > outv) minv = outv
       where (maxv < outv) maxv = outv
    end if
  end subroutine write_primal

  subroutine write_dual (dom, p, i, j, zlev, offs, dims, funit)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, funit
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                   :: d, id, idE, idN, idNE
    integer, dimension(TRIAG) :: leveldual
    real(8), dimension(TRIAG) :: relvort

    d = dom%id + 1

    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    relvort = get_vort (dom, i, j, offs, dims)

    if (maxval(dom%mask_n%elts((/id, idE, idNE/)+1)) >= ADJZONE) then
       ! avoid segfault if pre_levelout not used
       if (allocated(active_level%data)) leveldual(LORT+1) = maxval(active_level%data(d)%elts((/id, idE, idNE/)+1))

       write (funit,'(9(E14.5E2,1X), E14.5E2, 1X, I3)') dom%node%elts((/id, idE, idNE/)+1), relvort(LORT+1), leveldual(LORT+1)
    end if

    if (maxval(dom%mask_n%elts((/id, idNE, idN/)+1)) >= ADJZONE) then
       ! avoid segfault if pre_levelout not used
       if (allocated(active_level%data)) leveldual(UPLT+1) = maxval(active_level%data(d)%elts((/id, idNE, idN/)+1))

       write (funit,'(9(E14.5E2,1X), E14.5E2, 1X, I3)') dom%node%elts((/id, idNE, idN/)+1), relvort(UPLT+1), leveldual(UPLT+1)
    end if
  end subroutine write_dual

  function get_vort (dom, i, j, offs, dims)
    ! Averages vorticity to get smooth field for visualization
    implicit none
    real(8), dimension(TRIAG)      :: get_vort
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idE, idS, idW

    id  = idx(i,   j,   offs, dims)
    idE = idx(i+1, j,   offs, dims)
    idN = idx(i,   j+1, offs, dims)
    idW = idx(i-1, j,   offs, dims)
    idS = idx(i,   j-1, offs, dims)

    get_vort(UPLT+1) = ( &
         interp(dom%vort%elts(TRIAG*idW+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)+ &
         
         interp(dom%vort%elts(TRIAG*id+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
         
         interp(dom%vort%elts(TRIAG*idN+LORT+1), dom%vort%elts(TRIAG*id+UPLT+1)) &
         *dom%len%elts(EDGE*idN+RT+1)*dom%pedlen%elts(EDGE*idN+RT+1)) &
         
         / (dom%len%elts(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)+ &
         dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)+ &
         dom%len%elts(EDGE*idN+RT+1)*dom%pedlen%elts(EDGE*idN+RT+1))

    get_vort(LORT+1) = ( &
         interp(dom%vort%elts(TRIAG*idS+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)+ &
         
         interp(dom%vort%elts(TRIAG*id+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
         
         interp(dom%vort%elts(TRIAG*idE+UPLT+1), dom%vort%elts(TRIAG*id+LORT+1)) &
         *dom%len%elts(EDGE*idE+UP+1)*dom%pedlen%elts(EDGE*idE+UP+1)) &
         
         / (dom%len%elts(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)+ &
         dom%len%elts(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)+ &
         dom%len%elts(EDGE*idE+UP+1)*dom%pedlen%elts(EDGE*idE+UP+1))
  end function get_vort

  subroutine write_u_wc (dom, p, i, j, offs, dims, fid)
    ! Write wavelet coefficients of velocity
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, info, fid, k

    id = idx(i, j, offs, dims)

    do k = zmin, zmax
       do e = 1, EDGE
          write (fid,*) wav_coeff(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do
  end subroutine write_u_wc

  subroutine dump_adapt_mpi (id, run_id)
    ! Save data in check point files for restart
    ! One file per domain
    ! NOTE: modifies grid structure
    implicit none
    integer      :: id
    character(*) :: run_id

    integer                          :: c, d, ibeg, iend, info, j, k, l, p_chd, p_lev, p_par, r, v
    integer, dimension(1:size(grid)) :: fid_no, fid_gr
    character(9999)                  :: filename_gr, filename_no
    character(9999)                  :: bash_cmd, cmd_archive, cmd_files, command
    logical, dimension(1:N_CHDRN)    :: required
    
    call update_array_bdry (wav_coeff(scalars(1):scalars(2),zmin:zmax), NONE)
    if (vert_diffuse) call update_vector_bdry (wav_tke, NONE, 20)

    do d = 1, size(grid)
       do k = zmin, zmax
          do v = scalars(1), scalars(2)
             scalar => sol(v,k)%data(d)%elts
             wc_s   => wav_coeff(v,k)%data(d)%elts
             call apply_interscale_d (Restrict_scalar, grid(d), min_level-1, k, 0, 1) ! +1 to include poles
             nullify (scalar, wc_s)
          end do
       end do
       if (vert_diffuse) then
          do k = 1, zlevels
             scalar => tke(k)%data(d)%elts
             wc_s   => wav_tke(k)%data(d)%elts
             call apply_interscale_d (Restrict_scalar, grid(d), min_level-1, k, 0, 1) ! +1 to include poles
             nullify (scalar, wc_s)
          end do
       end if
    end do
    
    do r = 1, n_process
#ifdef MPI       
       if (r /= rank+1) then ! write only if our turn, otherwise wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if
#endif       
       do d = 1, size(grid)
          fid_no(d) = id*1000 + 1000000 + glo_id(rank+1,d)
          fid_gr(d) = id*1000 + 3000000 + glo_id(rank+1,d)

          write (filename_no, '(A,A,I4.4,A,I5.5)') trim (run_id), "_coef.", id, "_", glo_id(rank+1,d)
          write (filename_gr, '(A,A,I4.4,A,I5.5)') trim (run_id), "_grid.", id, "_", glo_id(rank+1,d)

          open (unit=fid_no(d), file=trim(filename_no), form="UNFORMATTED", action='WRITE', status='REPLACE')
          open (unit=fid_gr(d), file=trim(filename_gr), form="UNFORMATTED", action='WRITE', status='REPLACE')
       end do
    end do

    do d = 1, size(grid)
       write (fid_no(d)) istep_cumul
       write (fid_no(d)) time
       call dump (fid_no(d))

       ! Write data at coarsest scale (scaling functions)
       call apply_to_pole_d (write_scalar, grid(d), min_level-1, z_null, fid_no(d), .true.)
       p_par = 1
       do k = zmin, zmax
          do v = 1, N_VARIABLE
             ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
             iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
             write (fid_no(d)) sol(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do

       if (vert_diffuse) then
          do k = 1, zlevels
             ibeg = grid(d)%patch%elts(p_par+1)%elts_start + 1
             iend = ibeg + PATCH_SIZE**2 - 1
             write (fid_no(d)) tke(k)%data(d)%elts(ibeg:iend)
          end do
       end if

       ! Write wavelets at finer scales
       do l = min_level, level_end
          p_lev = 0
          do j = 1, grid(d)%lev(l)%length
             p_par = grid(d)%lev(l)%elts(j)
             if (grid(d)%patch%elts(p_par+1)%deleted) then
                do c = 1, N_CHDRN
                   p_chd = grid(d)%patch%elts(p_par+1)%children(c)
                   if (p_chd > 0) grid(d)%patch%elts(p_chd+1)%deleted = .true. 
                end do
                cycle ! No data to write
             end if

             do k = zmin, zmax
                do v = 1, N_VARIABLE
                   ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
                   iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
                   write (fid_no(d)) wav_coeff(v,k)%data(d)%elts(ibeg:iend)
                end do
             end do

             if (vert_diffuse) then
                do k = 1, zlevels
                   ibeg = grid(d)%patch%elts(p_par+1)%elts_start + 1
                   iend = ibeg + PATCH_SIZE**2 - 1
                   write (fid_no(d)) wav_tke(k)%data(d)%elts(ibeg:iend)
                end do
             end if

             ! Record whether patch needs to be refined
             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par+1)%children(c)
                if (p_chd > 0) then
                   required(c) = check_child_required(grid(d), p_par, c-1)
                   grid(d)%patch%elts(p_chd+1)%deleted = .not. required(c) 
                   if (required(c) .and. p_lev+1 <= size(grid(d)%lev(l+1)%elts)) then
                      p_lev = p_lev + 1
                      grid(d)%lev(l+1)%elts(p_lev) = p_chd ! ** grid modified **
                   end if
                else
                   required(c) = .false.
                end if
             end do
             write (fid_gr(d)) required
          end do
          if (l+1 <= max_level) grid(d)%lev(l+1)%length = p_lev ! ** grid modified **
       end do
       close (fid_no(d)); close (fid_gr(d))
    end do

    ! Archive checkpoint (overwriting existing checkpoint if present)
    call barrier ! Make sure all processors have written data
    if (rank == 0) then
       write (cmd_files, '(a,a,i4.4,a,a,a,i4.4)') trim (run_id), '{_grid,_coef}.', cp_idx , '_????? ', &
            trim (run_id), '_conn.', cp_idx
       write (cmd_archive, '(a,i4.4,a)') trim (run_id)//'_checkpoint_' , cp_idx, ".tgz"
       write (command,    '(a,a,a,a,a)') 'gtar cfz ', trim (cmd_archive), ' ', trim (cmd_files), ' --remove-files'
       write (bash_cmd,       '(a,a,a)') 'bash -c "', trim (command), '"'
       call system (trim(bash_cmd), info)
       if (info /= 0) then
          if (rank == 0) write (6,'(a)') "gtar command not present ... aborting"
          call abort
       end if
    end if
    call barrier ! make sure data is archived before restarting
  end subroutine dump_adapt_mpi

  subroutine write_scalar (dom, p, i, j, zlev, offs, dims, fid)
    ! For poles
    implicit none
    type(Domain)                   :: dom
    integer                        :: fid, i, j, p, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

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

  subroutine load_adapt_mpi (id, run_id)
    ! Read data from check point files for restart
    ! One file per domain
    implicit none
    integer      :: id
    character(*) :: run_id
    
    integer                          :: c, d, i, ibeg, iend, j, k, l, old_n_patch, p_chd, p_par, r, v
    integer, dimension(1:size(grid)) :: fid_no, fid_gr
    character(9999)                  :: filename_gr, filename_no
    character(9999)                  :: bash_cmd, cmd_archive, cmd_files, command
    logical, dimension(1:N_CHDRN)    :: required

    do r = 1, n_process
#ifdef MPI
       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if
#endif
       do d = 1, size(grid)
          fid_no(d) = id*1000 + 1000000 + glo_id(rank+1,d)
          fid_gr(d) = id*1000 + 3000000 + glo_id(rank+1,d)

          write (filename_no, '(A,A,I4.4,A,I5.5)') trim (run_id), "_coef.", id, "_", glo_id(rank+1,d)
          write (filename_gr, '(A,A,I4.4,A,I5.5)') trim (run_id), "_grid.", id, "_", glo_id(rank+1,d)

          open (unit=fid_no(d), file=trim(filename_no), form="UNFORMATTED", action='READ', status='OLD')
          open (unit=fid_gr(d), file=trim(filename_gr), form="UNFORMATTED", action='READ', status='OLD')
       end do
    end do
    
    ! Load coarsest scale solution (scaling functions)
    do d = 1, size(grid)
       read (fid_no(d)) istep_cumul
       read (fid_no(d)) time
       call load (fid_no(d))

       call apply_to_pole_d (read_scalar, grid(d), min_level-1, z_null, fid_no(d), .true.)

       p_par = 1
       do k = zmin, zmax
          do v = 1, N_VARIABLE
             ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
             iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
             read (fid_no(d)) sol(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do

       if (vert_diffuse) then
          do k = 1, zlevels
             ibeg = grid(d)%patch%elts(p_par+1)%elts_start + 1
             iend = ibeg + PATCH_SIZE**2 - 1
             read (fid_no(d)) tke(k)%data(d)%elts(ibeg:iend)
          end do
       end if
    end do

    ! Load finer scales (wavelets) if present
    ! (level_end is initially level_start and is incremented by refine_patch1 if children are present)
    l = 1
    do while (level_end > l) ! new level was added -> proceed to it
       l = level_end 
       if (rank == 0) write (6,'(a,i2)') 'Loading level ', l
       do d = 1, size(grid)
          old_n_patch = grid(d)%patch%length
          do j = 1, grid(d)%lev(l)%length
             p_par = grid(d)%lev(l)%elts(j)
             do k = zmin, zmax
                do v = 1, N_VARIABLE
                   ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
                   iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1
                   read (fid_no(d)) wav_coeff(v,k)%data(d)%elts(ibeg:iend)
                end do
             end do

             if (vert_diffuse) then
                do k = 1, zlevels
                   ibeg = grid(d)%patch%elts(p_par+1)%elts_start + 1
                   iend = ibeg + PATCH_SIZE**2 - 1
                   read (fid_no(d)) wav_tke(k)%data(d)%elts(ibeg:iend)
                end do
             end if

             read (fid_gr(d)) required
             do c = 1, N_CHDRN
                if (required(c)) call refine_patch1 (grid(d), p_par, c-1)
             end do
          end do
          do p_par = 2, old_n_patch
             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par)%children(c)
                if (p_chd+1 > old_n_patch) call refine_patch2 (grid(d), p_par-1, c-1)
             end do
          end do
       end do
       call post_refine
    end do
    sol%bdry_uptodate       = .false.
    wav_coeff%bdry_uptodate = .false.
    if (vert_diffuse) then
       tke%bdry_uptodate     = .false.
       wav_tke%bdry_uptodate = .false.
    end if

    do d = 1, size(grid)
       close(fid_no(d)); close(fid_gr(d))
    end do

    ! Delete temporary files
    call barrier ! Do not delete files before everyone has read them
    if (rank == 0) then
       write (cmd_files, '(a,a,i4.4,a,a,a,i4.4)') &
            trim (run_id), '{_grid,_coef}.', cp_idx , '_????? ', trim (run_id), '_conn.', cp_idx
       write (command,    '(a,a)') '\rm ', trim (cmd_files)
       write (bash_cmd, '(a,a,a)') 'bash -c "', trim (command), '"' 
       call system (trim(bash_cmd))
    end if
    call barrier
  end subroutine load_adapt_mpi

  subroutine read_scalar (dom, p, i, j, zlev, offs, dims, fid)
    ! For poles
    implicit none
    type(Domain)                   :: dom
    integer                        :: fid, i, j, p, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

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
    integer         :: d, d_glo, j, l
    character(9999) :: filename
    character(9999) :: bash_cmd, cmd_archive, cmd_files, command

    call update_bdry (topography, NONE)

    allocate (topo_count(min_level:max_level,1:size(grid))); topo_count = 0

    ! Save a separate file for each domain and each level
    do l = max_level, min_level, -1
       do d = 1, size(grid)
          d_glo = glo_id(rank+1,d) + 1
       
          write (filename, '(a,a,i2.2,a,i5.5)') trim (topo_file), '.', l, '.', d_glo           
          open (unit=10, file=trim (filename), form="UNFORMATTED", action='WRITE', status='REPLACE')
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (write_topo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          close (10)
       end do
    end do

    write (filename, '(a,a)') trim (topo_file), '.count'
    open (unit=10, file=trim (filename), form="UNFORMATTED", action='WRITE', status='REPLACE')
    write (10) topo_count
    close (10)
    
    ! Compress topography data
    write (cmd_files,   '(a,a,a,a,a,a,a)') trim (topo_file), '.', '??', '.', '?????', ' ', trim (filename)
    write (cmd_archive,         '(a)') trim (topo_file)//".tgz"
    write (6,             '(/,a,a,/)') 'Saving topography file ', trim (cmd_archive)
    write (command,     '(a,a,a,a,a)') 'gtar czf ', trim (cmd_archive), ' ', trim (cmd_files), ' --remove-files'
    write (bash_cmd,        '(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))
  end subroutine save_topo
  
  subroutine write_topo (dom, i, j, zlev, offs, dims)
    ! Write out level, coordinates and topography height
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, l

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    l = dom%level%elts(id) ! level

    write (10) dom%node%elts(id), topography%data(d)%elts(id)

    topo_count(l,d) = topo_count(l,d) + 1
  end subroutine write_topo

  subroutine load_topo
    ! Read topography data from for restart
    ! (one file per domain)
    !
    ! !! assumes topgraphy data was saved on a non-adaptive grid !!
    !
    implicit none
    integer         :: d, d_glo, l, r
    character(9999) :: filename
    character(9999) :: bash_cmd, cmd_archive, cmd_files, command

    ! Uncompress topography data
    if (rank == 0) then
       write (cmd_archive, '(a)') trim (topo_file)//".tgz"
       write (6,'(/,a,a,/)') 'Loading topography file ', trim (cmd_archive)
       write (command, '(a,a)') 'gtar xzf ', trim (cmd_archive)
       write (bash_cmd,    '(a,a,a)') 'bash -c "', trim (command), '"'
       call system (trim(bash_cmd))
    end if
    call barrier
    
    ! Allocate topo_count matrix for all domains on all ranks 
    if (allocated (topo_count)) deallocate (topo_count)
    allocate (topo_count(min_level:max_level,1:N_GLO_DOMAIN))
    
    do r = 1, n_process
#ifdef MPI
       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if
#endif
       ! Read topography count for all domains
       write (filename, '(a,a)') trim (topo_file), '.count'
       open (unit=10, file=trim (filename), form="UNFORMATTED", action='READ', status='OLD')
       read (10) topo_count
       close (10)

       ! Allocate topgraphy data structures    
       allocate (topography_data(min_level:max_level, 1:size(grid)))
       do l = min_level, max_level
          do d = 1, size(grid)
             d_glo = glo_id(rank+1,d) + 1
             allocate (topography_data(l,d)%node(1:topo_count(l,d_glo)))
             allocate (topography_data(l,d)%elts(1:topo_count(l,d_glo)))
          end do
       end do
    end do

    ! Load topography data
    do l = min_level, max_level
       do d = 1, size(grid)
          d_glo = glo_id(rank+1,d) + 1
          write (filename, '(a,a,i2.2,a,i5.5)') trim (topo_file), '.', l, '.', d_glo
          open (unit=10, file=trim (filename), form="UNFORMATTED", action='READ', status='OLD')

          do ii = 1, topo_count(l,d_glo)
             read (10) topography_data(l,d)%node(ii), topography_data(l,d)%elts(ii)
          end do
          close (10)
       end do
    end do

    ! Remove temporary files
    call barrier ! do not delete files before everyone has read them
    if (rank == 0) then
       write (cmd_files,   '(a,a,a,a,a,a,a,a)') trim (topo_file), '.', '??', '.', '?????', ' ', trim (topo_file)//'.count'
       write (command,     '(a,a)') '\rm ', trim (cmd_files)
       write (bash_cmd,  '(a,a,a)') 'bash -c "', trim (command), '"' 
       call system (trim(bash_cmd))
    end if
  end subroutine load_topo

  subroutine assign_topo (dom, i, j, zlev, offs, dims)
    ! Assign topography data to topography structure for simulation
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                            :: d, id, ii, jj, l, n_topo
    real(8), dimension(:), allocatable :: distance

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1
    l  = dom%level%elts(id)

    n_topo = size (topography_data(l,d)%elts); allocate (distance(1:n_topo))
    do ii = 1, n_topo
       distance(ii) = dist (dom%node%elts(id), topography_data(l,d)%node(ii))
    end do
    jj = minloc (distance,1) ; deallocate (distance)

    topography%data(d)%elts(id) = topography_data(l,d)%elts(jj)
  end subroutine assign_topo

  subroutine proj_xz_plane (cin, cout)
    implicit none
    type(Coord)                        :: cin
    real(8), dimension(2), intent(out) :: cout

    if (cin%y > 0) then
       cout = (/cin%x-radius, cin%z/)
    else
       cout = (/cin%x+radius, cin%z/)
    end if
  end subroutine proj_xz_plane

  subroutine error (msg)
    implicit none
    character(*) :: msg
    write(0,*) "ERROR: ", msg
  end subroutine error

  subroutine read_lonlat_from_binary (arr, n, fid)
    !     Use: real(8) arr(n_lon,n_lat)
    !     call read_lonlat_from_binary(arr(1,1),n_lon*n_lat,fid)
    implicit none
    integer               :: n, fid
    real(8), dimension(n) :: arr

    integer :: i

    read(fid) (arr(i),i=1,n)
  end subroutine read_lonlat_from_binary

  subroutine read_HR_optim_grid
    ! Reads in Heikes & Randall (1995) optimized grid from file in directory grid_HR
    ! Need to provide a symbolic link to grid_HR in working directory
    implicit none
    integer                        :: d_glo, d_HR, d_sub, fid, loz, p, r
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(19+1)                :: filename

    maxerror = 0d0
    l2error = 0d0

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call apply_onescale2 (ccentre, level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -1, 1)
    call apply_onescale (check_d,  level_end-1, z_null,  0, 0)

    l2error = sqrt (sum_real (l2error))
    maxerror = sync_max_real (maxerror)

    if (rank == 0) then
       write (6,'(A)') '-------------------------------------------------------&
            ---------------------------------------------------------------------------'
       write (6,'(A,i2,A,/)') 'Heikes-Randall optimizations of level ', level_start-1, ' grid:'
       write (6,'(A,2(es8.2,A))') 'Grid quality before optimization = ', maxerror, ' m (linf) ', l2error, ' m (l2)'
    end if

    fid = get_fid()
    if (level_start /= level_end) then
       write (0,'(i2,1x,i2)') level_end, level_start
       write (0,'(A)') "Reading HR grid points for level_start not equal to level_end not implemented"
       return
    end if

    do r = 1, n_process
#ifdef MPI
       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if
#endif

       write (filename, '(A,I1)')  "grid_HR/J", level_start-1
       open (unit=fid, file=filename, status='OLD')

       p = 1
       do d_HR = 1, N_ICOSAH_LOZENGE
          loz = dom_id_from_HR_id(d_HR)
          do d_sub = 1, N_SUB_DOM
             d_glo = loz * N_SUB_DOM + sub_dom_id_from_HR_sub_id (d_sub)
             if (owner(d_glo+1) == rank) call get_offs_Domain (grid(loc_id(d_glo+1)+1), p, offs, dims)
             call coord_from_file (d_glo, PATCH_LEVEL, fid, offs, dims, (/0, 0/))
          end do
       end do
       close(fid)
    end do

    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre,    level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,      level_end-1, z_null, -1, 1)
    call apply_onescale2 (check_grid, level_end-1, z_null,  0, 0)

    maxerror = 0d0
    l2error = 0d0

    call comm_nodes3_mpi (get_coord, set_coord, NONE)

    call apply_onescale2 (ccentre, level_end-1, z_null, -2, 1)
    call apply_onescale2 (midpt,   level_end-1, z_null, -1, 1)
    call apply_onescale (check_d,  level_end-1, z_null,  0, 0)

    l2error = sqrt (sum_real(l2error))
    maxerror = sync_max_real (maxerror)
    if (rank == 0) then
       write (6,'(A,2(es8.2,A))') 'Grid quality after optimization  = ', maxerror, ' m (linf) ', l2error, ' m (l2)'
       write (6,'(A)') '(distance between midpoints of primal and dual edges)'
       write (6,'(A)') '-------------------------------------------------------&
     ---------------------------------------------------------------------------'
    end if
  end subroutine read_HR_optim_grid

  integer function dom_id_from_HR_id (d_HR)
    ! d_HR: lozenge id as used by Heikes & Randall (starts from 1)
    ! results: domain id (starts from 0)
    implicit none
    integer :: d_HR

    dom_id_from_HR_id = modulo(d_HR,2)*5 + modulo(d_HR/2-1,5)
  end function dom_id_from_HR_id

  integer function sub_dom_id_from_HR_sub_id (sub_id)
    ! sub_id: lozenge sub id as used by Heikes & Randall (starts from 1)
    ! results: sub domain id (starts from 0)
    implicit none
    integer :: sub_id

    integer :: id, i, j, halv_sub_dom, l, jdiv, idiv

    i = 0
    j = 0
    id = sub_id - 1
    halv_sub_dom = N_SUB_DOM/2
    do l = DOMAIN_LEVEL-1, 0, -1
       jdiv = id/halv_sub_dom
       j = j + jdiv*2**l
       id = modulo (id+4**l,4**(l+1))
       idiv = id/halv_sub_dom
       i = i + idiv*2**l
       halv_sub_dom = halv_sub_dom/4
       id = modulo (id,4**l)
    end do
    sub_dom_id_from_HR_sub_id = j*N_SUB_DOM_PER_DIM + i
  end function sub_dom_id_from_HR_sub_id

  subroutine zrotate (c_in, c_out, angle)
    implicit none
    real(8),      intent(in) :: angle
    type(Coord),  intent(in) :: c_in
    type(Coord), intent(out) :: c_out

    c_out%x =  c_in%x * cos(angle) - c_in%y * sin(angle)
    c_out%y =  c_in%x * sin(angle) + c_in%y * cos(angle)
    c_out%z =  c_in%z
  end subroutine zrotate

  recursive subroutine coord_from_file (d_glo, l, fid, offs, dims, ij0)
    implicit none
    integer,                        intent(in) :: d_glo, l, fid
    integer, dimension(2),          intent(in) :: ij0
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    integer               :: d_loc, k
    integer, dimension(2) :: ij
    type(Coord)           :: node, node_r

    d_loc = loc_id(d_glo+1)
    do k = 1, 4
       ij = ij0 + HR_offs(:,k) * 2**(l-1)
       if (l == 1) then
          if (owner(d_glo+1) == rank) then
             read(fid,*) node
             call zrotate (node, node_r, -0.5d0) ! icosahedron orientation good for tsunami
             grid(d_loc+1)%node%elts(idx(ij(1), ij(2), offs, dims) + 1) = project_on_sphere(node_r)
          else ! if domain is on another process, still read to get to correct position in file
             read(fid,*)
          end if
       else
          call coord_from_file (d_glo, l-1, fid, offs, dims, ij)
       end if
    end do
  end subroutine coord_from_file

  subroutine pre_levelout
    implicit none
    integer :: d, l, max_output_level, num

    ! FIXME cleaner would be to use init_io routine
    call init_Float_Field (active_level, AT_NODE)

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init (active_level%data(d), num)
       active_level%data(d)%elts(1:num) = grid(d)%level%elts(1:num)
    end do

    do l = level_end-1, level_start, -1
       call apply_interscale (restrict_level, l, z_null, 0, 1)
    end do
  end subroutine pre_levelout

  ! now active_level can be used

  subroutine post_levelout
    implicit none
    integer :: d

    do d = 1, size(grid)
       deallocate (active_level%data(d)%elts)
    end do
    deallocate (active_level%data)
  end subroutine post_levelout

  subroutine restrict_level (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: d, id_par, id_chd

    d = dom%id+1

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_chd+1) >= ADJZONE) active_level%data(d)%elts(id_par+1) = active_level%data(d)%elts(id_chd+1)
  end subroutine restrict_level

  subroutine write_and_export (isave)
    use utils_mod
    implicit none
    integer :: isave

    integer        :: d, i, j, k, k1, k2, l, p, u1, u2
    character(3)   :: layer
    character(7)   :: var_file
    character(255) :: prefix

    u1 = 1000000 + 100*isave
    u2 = 2000000 + 100*isave

    if (save_zlev < 0) then ! save all layers
       k1 = 1; k2 = zlevels 
    else
       k1 = save_zlev; k2 = save_zlev
    end if

    if (rank == 0) write(6,'(/,A,i4/)') 'Saving fields ', isave

    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE, 26)

    call pre_levelout

    ! Compute surface pressure
    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))

    do l = level_start, level_end
       minv = 1d63; maxv = -1d63
       if (compressible .or. zlevels /= 2 .or. .not. mode_split) then
          do k = 1, k2
             ! Compute pressure, exner and geopotential
             do d = 1, size(grid)
                mass    => sol(S_MASS,k)%data(d)%elts
                mean_m  => sol_mean(S_MASS,k)%data(d)%elts
                temp    => sol(S_TEMP,k)%data(d)%elts
                mean_t  => sol_mean(S_TEMP,k)%data(d)%elts
                exner   => exner_fun(k)%data(d)%elts
                do j = 1, grid(d)%lev(l)%length
                   call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
                end do
                nullify (mass, mean_m, temp, mean_t, exner)
             end do

             if (save_zlev < 0 .or. k == k2) then ! Compute zonal/meridional velocities and vorticity
                ! Zonal and meridional velocities 
                do d = 1, size(grid)
                   velo  => sol(S_VELO,k)%data(d)%elts
                   velo1 => grid(d)%u_zonal%elts
                   velo2 => grid(d)%v_merid%elts
                   do j = 1, grid(d)%lev(l)%length
                      call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
                   end do
                   nullify (velo, velo1, velo2)
                end do

                ! Vorticity
                do d = 1, size(grid)
                   velo  => sol(S_VELO,k)%data(d)%elts
                   vort  => grid(d)%vort%elts
                   do j = 1, grid(d)%lev(l)%length
                      call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
                   end do
                   call apply_to_penta_d (post_vort, grid(d), l, z_null)
                   nullify (velo, vort)
                end do

                ! Calculate vorticity at hexagon points (stored in ke)
                do d = 1, size(grid)
                   vort => grid(d)%ke%elts
                   do j = 1, grid(d)%lev(l)%length
                      call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
                   end do
                   nullify (vort)
                end do
                
                ! Save fields
                write (layer, '(i3.3)') k
                prefix = trim(run_id)//'_'//layer

                call write_level_mpi (write_primal, u1+l, l, k, .true., prefix)
                do i = 1, N_VAR_OUT
                   minv(i) = -sync_max_real (-minv(i))
                   maxv(i) =  sync_max_real ( maxv(i))
                end do
                if (rank == 0) then
                   write (var_file, '(i7)') u1
                   open (unit=50, file=trim(prefix)//'.'//var_file, status='replace')
                   write (50,'(a, 7(e15.5e2, 1x), i3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
                   write (50,'(a, 7(e15.5e2, 1x), i3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
                   close (50)
                end if
                call write_level_mpi (write_dual, u2+l, l, k, .false., prefix)
             end if
           end do
       else
          ! Barotropic velocity
          do d = 1, size(grid)
             velo => sol(S_VELO,zlevels+1)%data(d)%elts ! barotropic velocity
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (barotropic_velocity, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
             end do
             nullify (velo)
          end do
          do d = 1, size(grid)
             velo  => sol(S_VELO,zlevels+1)%data(d)%elts
             velo1 => grid(d)%u_zonal%elts
             velo2 => grid(d)%v_merid%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (velo, velo1, velo2)
          end do

          ! Barotropic vorticity
          do d = 1, size(grid)
             velo  => sol(S_VELO,zlevels+1)%data(d)%elts
             vort  => grid(d)%vort%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
             end do
             call apply_to_penta_d (post_vort, grid(d), l, z_null)
             nullify (velo, vort)
          end do

          ! Calculate barotropic vorticity at hexagon points 
          do d = 1, size(grid)
             vort => grid(d)%ke%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (vort)
          end do

          ! Baroclinic velocity in top layer
          do d = 1, size(grid)
             velo1  => sol(S_VELO,zlevels+1)%data(d)%elts   ! barotropic velocity
             velo2  => sol(S_VELO,zlevels)%data(d)%elts     ! velocity in top layer
             velo   => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in top layer
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (baroclinic_velocity, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
             end do
             nullify (velo, velo1, velo2)
          end do
          do d = 1, size(grid)
             velo  => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in top layer
             velo1 => trend(S_VELO,1)%data(d)%elts
             velo2 => trend(S_VELO,2)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (velo, velo1, velo2)
          end do

          ! Baroclinic vorticity in top layer
          do d = 1, size(grid)
             velo  => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in top layer
             vort  => grid(d)%vort%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
             end do
             call apply_to_penta_d (post_vort, grid(d), l, z_null)
             nullify (velo, vort)
          end do

          ! Calculate baroclinic vorticity in top layer at hexagon points 
          do d = 1, size(grid)
             vort => grid(d)%geopot_lower%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do
             nullify (vort)
          end do

          ! Save fields
          prefix = trim(run_id)
          call write_level_mpi (write_primal, u1+l, l, z_null, .true., prefix)
          do i = 1, N_VAR_OUT
             minv(i) = -sync_max_real (-minv(i))
             maxv(i) =  sync_max_real ( maxv(i))
          end do
          if (rank == 0) then
             write (var_file, '(i7)') u1
             open (unit=50, file=trim(prefix)//'.'//var_file, status='replace')
             write (50,'(a, 11(e15.5e2, 1x), i3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
             write (50,'(a, 11(e15.5e2, 1x), i3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
             close (50)
          end if
          call write_level_mpi (write_dual, u2+l, l, z_null, .false., prefix)
       end if
    end do
    call post_levelout
    call barrier
    
    if (rank == 0) then
       do k = k1, k2
          write (layer, '(i3.3)') k
          prefix = trim(run_id)//'_'//layer
          call compress_files (isave, prefix)
       end do
    end if
    call barrier
  end subroutine write_and_export

  subroutine barotropic_velocity (dom, i, j, zlev, offs, dims)
    ! Calculate barotropic velocity in two-layer model
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: d, e, id, id_e, id_i, idE, idNE, idN, k
    real(8)                       :: dz0
    real(8), dimension (1:EDGE,2) :: dz

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    do k = 1, 2
       dz0 = sol_mean(S_MASS,k)%data(d)%elts(id+1) + sol(S_MASS,k)%data(d)%elts(id_i)
       dz(RT+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idE+1)  + sol(S_MASS,k)%data(d)%elts(idE+1))
       dz(DG+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idNE+1) + sol(S_MASS,k)%data(d)%elts(idNE+1))
       dz(UP+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idN+1)  + sol(S_MASS,k)%data(d)%elts(idN+1))
    end do

    do e = 1, EDGE
       id_e = EDGE*id + e
       velo(id_e) = (dz(e,1)*sol(S_VELO,1)%data(d)%elts(id_e) + dz(e,2)*sol(S_VELO,2)%data(d)%elts(id_e)) / (dz(e,1) + dz(e,2))
    end do
  end subroutine barotropic_velocity

  subroutine baroclinic_velocity (dom, i, j, zlev, offs, dims)
    ! Calculate baroclinic velocity in top layer
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e

    id = idx (i, j, offs, dims)
    d = dom%id + 1

    do e = 1, EDGE
       id_e = EDGE*id + e
       velo(id_e) = velo2(id_e) - velo1(id_e)
    end do
  end subroutine baroclinic_velocity

  subroutine omega_velocity
    ! Computes vertical velocity in pressure coordinates D_tP = OMEGA [Pa/s]
    ! stored in trend(S_TEMP,1:zlevels)
    implicit none
    integer :: d, j, k, l, p

    call update_array_bdry (sol, NONE)
    
    ! Compute surface pressure
    call cal_surf_press (sol)

    do k = 1, zlevels
       do l = level_end, level_start, -1
          ! Divergence of mass flux at each vertical level
          ! stored in trend(S_MASS,1:zlevels)
          do d = 1, size(grid)
             mass   =>      sol(S_MASS,k)%data(d)%elts
             velo   =>      sol(S_VELO,k)%data(d)%elts
             mean_m => sol_mean(S_MASS,k)%data(d)%elts
             h_flux => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=7)
             end do
             nullify (mass, mean_m, velo)
             if (l < level_end) then
                dscalar => trend(S_MASS,k)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)
                nullify (dscalar)
             end if
             nullify (h_flux)
          end do
          horiz_flux(S_MASS)%bdry_uptodate = .false.
          call update_bdry (horiz_flux(S_MASS), l)
          do d = 1, size(grid)
             dscalar => trend(S_MASS,k)%data(d)%elts
             h_flux  => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do
          
          ! u grad(P_S) at hexagon centres
          ! stored in trend(S_TEMP,1:zlevels)
          do d = 1, size(grid)
             scalar =>         grid(d)%surf_press%elts
             velo   =>      sol(S_VELO,k)%data(d)%elts
             h_flux => horiz_flux(S_TEMP)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=6)
             end do
             if (l < level_end) then
                dscalar => trend(S_TEMP,k)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)
                nullify (dscalar)
             end if
             nullify (h_flux, scalar, velo)
          end do
          horiz_flux(S_TEMP)%bdry_uptodate = .false.
          call update_bdry (horiz_flux(S_TEMP), l)
          do d = 1, size(grid)
             dscalar =>    trend(S_TEMP,k)%data(d)%elts
             h_flux  => horiz_flux(S_TEMP)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do
          
          trend(S_MASS:S_TEMP,k)%bdry_uptodate = .false.
          call update_vector_bdry (trend(S_MASS:S_TEMP,k), l)
       end do
    end do

    ! Compute Omega (velocity flux across interfaces)
    ! stored in trend(S_TEMP,1:zlevels)
    call apply (cal_omega, z_null)

    trend(S_TEMP,1:zlevels)%bdry_uptodate = .false.
    call update_vector_bdry (trend(S_TEMP,:), NONE)
  end subroutine omega_velocity
  
  subroutine cal_omega (dom, i, j, zlev, offs, dims)
    ! Velocity flux across interfaces
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id_i, k
    real(8), dimension(1:zlevels)   :: u_gradP
    real(8), dimension(1:zlevels+1) :: div_mass

    d    = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    ! Vertically integrate div(mass flux) from top to bottom
    ! results at interfaces
    div_mass(zlevels+1) = 0d0 ! zero flux at top boundary
    do k = zlevels, 1, -1
       div_mass(k) = div_mass(k+1) + trend(S_MASS,k)%data(d)%elts(id_i)
    end do

    ! u.gradP at layers
    do k = 1, zlevels
       u_gradP(k) = interp (b_vert(k), b_vert(k+1)) * trend(S_TEMP,k)%data(d)%elts(id_i)
    end do
    
    ! Complete computation of omega
    do k = 1, zlevels
       trend(S_TEMP,k)%data(d)%elts(id_i) = - grav_accel * interp (div_mass(k), div_mass(k+1)) + u_gradP(k) 
    end do
  end subroutine cal_omega

  subroutine vertical_velocity
    ! Computes vertical velocity w [m/s]
    ! stored in trend(S_TEMP,1:zlevels)
    implicit none
    integer :: d, j, k, l, p

    call update_array_bdry (sol, NONE)
    
    ! Compute surface pressure
    call cal_surf_press (sol)

    do k = 1, zlevels
       do l = level_end, level_start, -1
          ! Divergence of mass flux at each vertical level
          ! stored in trend(S_MASS,1:zlevels)
          do d = 1, size(grid)
             mass   =>      sol(S_MASS,k)%data(d)%elts
             velo   =>      sol(S_VELO,k)%data(d)%elts
             mean_m => sol_mean(S_MASS,k)%data(d)%elts
             h_flux => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=7)
             end do
             nullify (mass, mean_m, velo)
             if (l < level_end) then
                dscalar => trend(S_MASS,k)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)
                nullify (dscalar)
             end if
             nullify (h_flux)
          end do
          horiz_flux(S_MASS)%bdry_uptodate = .false.
          call update_bdry (horiz_flux(S_MASS), l)
          do d = 1, size(grid)
             dscalar => trend(S_MASS,k)%data(d)%elts
             h_flux  => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do
       end do
    end do
    trend(S_MASS:S_TEMP,k)%bdry_uptodate = .false.
    call update_vector_bdry (trend(S_MASS:S_TEMP,k), l)


    ! Compute Omega (velocity flux across interfaces)
    ! stored in trend(S_TEMP,1:zlevels)
    call apply (cal_w, z_null)

    trend(S_TEMP,1:zlevels)%bdry_uptodate = .false.
  end subroutine vertical_velocity

  subroutine cal_w (dom, i, j, zlev, offs, dims)
    ! Vertical velocity
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id_i, k
    real(8), dimension(1:zlevels)   :: u_gradP
    real(8), dimension(1:zlevels+1) :: div_mass

    d    = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    ! Vertically integrate div(mass flux) from top to bottom
    !  d_t M = - div_mass(1), where M is the column mass
    ! results at interfaces
    div_mass(zlevels+1) = 0d0 ! zero flux at top boundary
    do k = zlevels, 1, -1
       div_mass(k) = div_mass(k+1) + trend(S_MASS,k)%data(d)%elts(id_i)
    end do
    
    do k = 1, zlevels
       ! Vertical mass flux through layers
       trend(S_TEMP,k)%data(d)%elts(id_i) = - grav_accel * interp (div_mass(k), div_mass(k+1)) &
            - interp (b_vert(k), b_vert(k+1)) * div_mass(1)

       ! Vertical velocity
       trend(S_TEMP,k)%data(d)%elts(id_i) = trend(S_TEMP,k)%data(d)%elts(id_i) / ref_density + proj_vel_vertical ()
    end do
  contains
    real(8) function proj_vel_vertical ()
      ! Computes grad_zonal(z) * u_zonal + grad_merid(z) * u_merid at hexagon centres for vertical velocity computation.
      ! Uses Perot formula as also used for kinetic energy:
      ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_center) ) / cell_area
      implicit none
      integer :: d, id, id_i, idE, idN, idNE, idS, idSW, idW, k

      idE  = idx (i+1, j,   offs, dims)
      idNE = idx (i+1, j+1, offs, dims)
      idN  = idx (i,   j+1, offs, dims)
      idW  = idx (i-1, j,   offs, dims)
      idSW = idx (i-1, j-1, offs, dims)
      idS  = idx (i,   j-1, offs, dims)

      velo => sol(S_VELO,k)%data(d)%elts
    
      proj_vel_vertical =  &
      (vert_vel (i,j,i+1,j,EDGE*id +RT+1) + vert_vel (i+1,j+1,i,j,EDGE*id  +DG+1) + vert_vel (i,j,i,j+1,EDGE*id +UP+1) + &
      (vert_vel (i-1,j,i,j,EDGE*idW+RT+1) + vert_vel (i,j,i-1,j-1,EDGE*idSW+DG+1) + vert_vel (i,j-1,i,j,EDGE*idS+UP+1))) / 6d0
      
      nullify (velo)
    end function proj_vel_vertical

    real(8) function vert_vel (i1, j1, i2, j2, id_e)
      implicit none
      integer :: i1, j1, i2, j2, id_e

      real(8) :: dl, dz

      dz =  z_i (dom, i2, j2, k, offs, dims, sol) - z_i (dom, i1, j1, k, offs, dims, sol)
      dl = dom%len%elts(id_e)

      vert_vel = dz / sqrt (dl**2 + dz**2) * velo(id_e)
    end function vert_vel
  end subroutine cal_w

  subroutine compress_files (isave, run_id)
    implicit none
    integer      :: isave
    character(*) :: run_id

    integer         :: info
    character(4)    :: s_time
    character(1300) :: bash_cmd, command

    write (s_time, '(i4.4)') isave

    command = 'ls -1 '//trim(run_id)//'.1'//s_time//'?? > tmp1'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))

    command = 'gtar cfz '//trim(run_id)//'.1'//s_time//'.tgz -T tmp1 --remove-files'
    call system (trim(command), info)
    if (info /= 0) then
       if (rank == 0) write (6,'(a)') "gtar command not present ... aborting"
       call abort
    end if

    command = 'ls -1 '//trim(run_id)//'.2'//s_time //'?? > tmp2'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))

    command = 'gtar cfz '//trim(run_id)//'.2'//s_time//'.tgz -T tmp2 --remove-files'
    call system (trim(command))

    command = '\rm -f tmp1 tmp2'
    call system (trim(command))
  end subroutine compress_files
end module io_mod
