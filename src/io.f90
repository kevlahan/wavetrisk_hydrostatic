module io_mod
  use kind_mod,   only : dp
  
  use shared_mod, only : ADJZONE, Coord, EDGE, N_BDRY,  N_CHDRN, N_VARIABLE, TRIAG, RT, DG, UP, MATH_PI, MULT, &
       cp_idx, istep_cumul, itime, iwrite, &
       level_end, min_level, max_level, n_glo_domain, grav_accel, kappa, p_0, R_d, pressure_save, radius, ref_density, run_id, &
       scalars, threshold, time, topo_file, topo_min_level, topo_max_level, vert_diffuse, zmin, zmax, zlevels, &
       S_MASS, S_TEMP, S_VELO, NONE, LORT, UPLT, z_null

  use arch_mod,         only : barrier, glo_id, n_process, rank
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
  public :: assign_NCAR_topo, dump_adapt_mpi, kinetic_energy, load_adapt_mpi, init_io_mod, load_topo, save_topo, vort_extrema
  public :: topo, pot_energy, total_ke, layer1_ke, layer2_ke, one_layer_ke, barotropic_ke, pot_enstrophy 
  
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
    ! Save data in check point files for restart
    ! One file per domain
    !
    ! NOTE: modifies adaptive grid structure by deleting any patches that do not contain cells in adjacent zone
    
    implicit none
    
    integer, intent(in) :: id
    
    integer                          :: c, d, ibeg, iend, j, k, l, p_chd, p_lev, p_par, r, v
    integer, dimension(1:size(grid)) :: fid_no, fid_gr
    character(4)                     :: cp4, id4
    character(5)                     :: gid5
    character(:), allocatable        :: filename_no, filename_gr
    character(:), allocatable        :: archive, cmd, files
    logical, dimension(1:N_CHDRN)    :: required

    call update_bdry (wav_coeff(scalars(1):scalars(2),zmin:zmax), NONE, 964)
    if (vert_diffuse) call update_bdry (wav_tke, NONE, 965)

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

    write (id4,'(i4.4)')  id
    
    do r = 1, n_process
#ifdef MPI       
       if (r /= rank+1) then ! write only if our turn, otherwise wait at barrier
          call barrier
          cycle 
       end if
#endif       
       do d = 1, size(grid)
          fid_no(d) = 1000000 + id*1000 + glo_id(rank+1,d)
          fid_gr(d) = 3000000 + id*1000 + glo_id(rank+1,d)

          write (gid5,'(i5.5)') glo_id(rank+1,d)
          
          filename_no = trim(run_id) // '_coef.' // id4 // '_' // gid5
          filename_gr = trim(run_id) // '_grid.' // id4 // '_' // gid5

          open (unit=fid_no(d), file=filename_no, form="UNFORMATTED", action='WRITE', status='REPLACE')
          open (unit=fid_gr(d), file=filename_gr, form="UNFORMATTED", action='WRITE', status='REPLACE')
       end do
    end do

    do d = 1, size(grid)
       write (fid_no(d)) istep_cumul
       write (fid_no(d)) time
       write (fid_no(d)) itime
       write (fid_no(d)) iwrite
       write (fid_no(d)) threshold
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

             ! Do not save any child patches of a deleted parent patch (patches not in adjacent zone)
             if (grid(d)%patch%elts(p_par+1)%deleted) then
                do c = 1, N_CHDRN
                   p_chd = grid(d)%patch%elts(p_par+1)%children(c)
                   if (p_chd > 0) grid(d)%patch%elts(p_chd+1)%deleted = .true. 
                end do
                cycle ! no data to write
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

             ! Record whether patch on finer grid is in adjacent zone (otherwise do not save)
             do c = 1, N_CHDRN
                p_chd = grid(d)%patch%elts(p_par+1)%children(c)
                if (p_chd > 0) then
                   required(c) = check_child_required (grid(d), p_par, c-1)
                   
                   grid(d)%patch%elts(p_chd+1)%deleted = .not. required(c) 
                   
                   if (required(c) .and. p_lev+1 <= size(grid(d)%lev(l+1)%elts)) then
                      p_lev = p_lev + 1
                      grid(d)%lev(l+1)%elts(p_lev) = p_chd 
                   end if
                else
                   required(c) = .false.
                end if
             end do
             write (fid_gr(d)) required ! save whether finer patch is required
          end do
          if (l+1 <= max_level) grid(d)%lev(l+1)%length = p_lev 
       end do
       close (fid_no(d)); close (fid_gr(d))
    end do

    ! Archive checkpoint (overwriting existing checkpoint if present)
    call barrier ! make sure all processors have written data
    if (rank == 0) then
       write (cp4,'(i4.4)') cp_idx
       archive = trim(run_id) // '_checkpoint_' // cp4 // '.tgz'

       files = &
            trim(run_id) // '_grid.' // cp4 // '_????? ' // &
            trim(run_id) // '_coef.' // cp4 // '_????? ' // &
            trim(run_id) // '_conn.' // cp4

       cmd = 'gtar cfz ' // archive // ' ' // files // ' --remove-files'
       call execute_command_line (cmd)
    end if
    call barrier ! make sure data is archived before restarting
  end subroutine dump_adapt_mpi
  

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
    ! Read data from check point files for restart
    ! One file per domain
    
    implicit none
    
    integer, intent(in) :: id
    
    integer                          :: c, d, ibeg, iend, j, k, l, old_n_patch, p_chd, p_par, r, v
    integer, dimension(1:size(grid)) :: fid_no, fid_gr
    Character(4)                     :: cp4, id4
    character(5)                     :: gid5
    character(:), allocatable        :: filename_no, filename_gr
    character(:), allocatable        :: cmd, files
    logical, dimension(1:N_CHDRN)    :: required

    write (id4,'(i4.4)')  id
    
    do r = 1, n_process
#ifdef MPI
       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call barrier
          cycle 
       end if
#endif
       do d = 1, size(grid)
          fid_no(d) = 1000000 + id*1000 + glo_id(rank+1,d)
          fid_gr(d) = 3000000 + id*1000 + glo_id(rank+1,d)

          write (gid5,'(i5.5)') glo_id(rank+1,d)

          filename_no = trim(run_id) // '_coef.' // id4 // '_' // gid5
          filename_gr = trim(run_id) // '_grid.' // id4 // '_' // gid5

          open (unit=fid_no(d), file=filename_no, form="UNFORMATTED", action='READ', status='OLD')
          open (unit=fid_gr(d), file=filename_gr, form="UNFORMATTED", action='READ', status='OLD')
       end do
    end do

    ! Load coarsest scale solution (scaling functions)
    do d = 1, size(grid)
       read (fid_no(d)) istep_cumul
       read (fid_no(d)) time
       read (fid_no(d)) itime
       read (fid_no(d)) iwrite
       read (fid_no(d)) threshold
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
       write (cp4,'(i4.4)') cp_idx

       files = &
            trim(run_id) // '_grid.' // cp4 // '_????? ' // &
            trim(run_id) // '_coef.' // cp4 // '_????? ' // &
            trim(run_id) // '_conn.' // cp4

       cmd = '\rm -f ' // files
       call execute_command_line (cmd)
    end if
  end subroutine load_adapt_mpi
  

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
#ifdef MPI
       if (r /= rank+1) then ! read only if our turn, otherwise wait at barrier
          call barrier
          cycle 
       end if
#endif
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
