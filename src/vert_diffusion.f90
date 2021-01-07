module vert_diffusion_mod
  ! Provides trend routines for vertical diffusion of velocity and buoyancy (temp variable), primarily for ocean models.
  ! (assumes forward Euler time step)
  !
  ! User must supply the following functions in test_case_mod.f90:
  !
  ! viscosity routines
  !        eddy_viscosity, eddy_diffusion
  !
  ! vertical boundary condition routines
  !        top_temp_source, bottom_temp_source
  !        top_velo_source, bottom_velo_source
  use test_case_mod
  implicit none
contains
  subroutine trend_vertical_diffusion (q, dq)
    ! Trend for eddy diffusivity and eddy viscosity 
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq

    integer :: d, k, p

    call update_array_bdry (q, NONE, 27)

    ! Scalars
    do d = 1, size(grid)
       scalar => q(S_MASS,zlevels+1)%data(d)%elts ! free surface
       do k = 1, zlevels
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          mass   =>  q(S_MASS,k)%data(d)%elts
          temp   =>  q(S_TEMP,k)%data(d)%elts
          velo   =>  q(S_VELO,k)%data(d)%elts
          dmass  => dq(S_MASS,k)%data(d)%elts
          dtemp  => dq(S_TEMP,k)%data(d)%elts
          dvelo  => dq(S_VELO,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (trend_scalars, grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_velo,    grid(d), p-1, k, 0, 0)
          end do
          nullify (dmass, dtemp, dvelo, mass, mean_m, mean_t, temp, velo)
       end do
       nullify (scalar)
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_vertical_diffusion

  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Vertical eddy diffusivity of buoyancy
    ! (layer height is not diffused)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: dz_k
    
    id_i = idx (i, j, offs, dims) + 1
      
    dmass(id_i) = 0.0_8

    dz_k = dz_i (dom, i, j, zlev, offs, dims)
    
    if (zlev > 1 .and. zlev < zlevels) then
       dtemp(id_i) = flux_scalar(1) - flux_scalar(-1)
    elseif (zlev == 1) then
       dtemp(id_i) =  flux_scalar( 1) + bottom_temp_source(dom, i, j, zlev, offs, dims)
    elseif (zlev == zlevels) then
       dtemp(id_i) =  top_temp_source(dom, i, j, zlev, offs, dims) - flux_scalar(-1)
    end if
    dtemp(id_i) = porous_density (dom, i, j, zlev, offs, dims) * dtemp(id_i)
  contains
    real(8) function flux_scalar (l)
      ! Computes flux at interface below (l=-1) or above (l=1) vertical level zlev
      implicit none
      integer :: l

      integer :: d
      real(8) :: b_0, b_l, dz_l, eta, mass_0, mass_l, ri, temp_0, temp_l, visc, z

      d = dom%id + 1

      eta = scalar(id_i)
      ri  = richardson (dom, i, j, zlev, offs, dims, l)
      z   = zl_i (dom, i, j, zlev, offs, dims, l)
      visc = eddy_diffusivity (eta, ri, z)

      mass_0 = mean_m(id_i) + mass(id_i)
      temp_0 = mean_t(id_i) + temp(id_i)
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
      temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
      b_l = temp_l / mass_l

      dz_l = 0.5 * (dz_k + dz_i(dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

      flux_scalar = l * visc * (b_l - b_0) / dz_l
    end function flux_scalar
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id
    real(8), dimension(1:EDGE) :: dz_k

    id = idx (i, j, offs, dims)
    
    dz_k = dz_e (dom, i, j, zlev, offs, dims)

    if (zlev > 1 .and. zlev < zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = (flux_velo(1) - flux_velo(-1)) / dz_k
    elseif  (zlev == 1) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = flux_velo(1)/dz_k + bottom_velo_source(dom, i, j, zlev, offs, dims)
    elseif (zlev == zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = top_velo_source(dom, i, j, zlev, offs, dims) - flux_velo(-1)/dz_k 
    end if
  contains
    function flux_velo (l)
      ! Flux at upper interface (l=1) or lower interface (l=-1)
      implicit none
      integer               :: l
      real(8), dimension(3) :: flux_velo

      integer               :: d
      real(8)               :: ri
      real(8), dimension(3) :: dz_l, eta, visc, z
      
      d = dom%id + 1

      eta = eta_e (dom, i, j, zlev, offs, dims)
      ri  = richardson (dom, i, j, zlev, offs, dims, l)
      z   = zl_e  (dom, i, j, zlev, offs, dims, l)
      visc = eddy_viscosity (eta, ri, z)

      dz_l = 0.5 * (dz_k + dz_e (dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

      flux_velo = l * visc * (sol(S_VELO,zlev+l)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - velo(EDGE*id+RT+1:EDGE*id+UP+1)) / dz_l
    end function flux_velo
  end subroutine trend_velo
end module vert_diffusion_mod
