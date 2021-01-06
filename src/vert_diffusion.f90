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

    integer :: d, id, id_i

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1
   
    dmass(id_i) = 0.0_8

    if (zlev > 1 .and. zlev < zlevels) then
       dtemp(id_i) = flux_scalar(1) - flux_scalar(-1)
    elseif (zlev == 1) then
       dtemp(id_i) =  flux_scalar( 1) + bottom_temp_source(dom, i, j, zlev, offs, dims)
    elseif (zlev == zlevels) then
       dtemp(id_i) =  top_temp_source(dom, i, j, zlev, offs, dims) - flux_scalar(-1)
    end if
    dtemp(id_i) = porous_density (d, id, zlev) * dtemp(id_i)
  contains
    real(8) function flux_scalar (itype)
      ! Computes flux at interface below (itype=-1) or above (itype=1) current level
      implicit none
      integer :: itype

      real(8) :: b_0, b_l, dz_l, mass_0, mass_l, temp_0, temp_l

      mass_0 = mean_m(id_i) + mass(id_i)
      temp_0 = mean_t(id_i) + temp(id_i)
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,zlev+itype)%data(d)%elts(id_i) + sol(S_MASS,zlev+itype)%data(d)%elts(id_i)
      temp_l = sol_mean(S_TEMP,zlev+itype)%data(d)%elts(id_i) + sol(S_TEMP,zlev+itype)%data(d)%elts(id_i)
      b_l = temp_l / mass_l

      dz_l = 0.5 * (mass_0 + mass_l) / porous_density (d, id, zlev)

      flux_scalar = itype * eddy_diffusivity (dom, i, j, zlev, offs, dims) * (b_l - b_0) / dz_l
    end function flux_scalar
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id
    real(8), dimension(1:EDGE) :: dz

    id = idx (i, j, offs, dims)
    
    dz = dz_e (dom, i, j, zlev, offs, dims)

    if (zlev > 1 .and. zlev < zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = (flux_velo(1) - flux_velo(-1)) / dz
    elseif  (zlev == 1) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = flux_velo(1)/dz + bottom_velo_source(dom, i, j, zlev, offs, dims)
    elseif (zlev == zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = top_velo_source(dom, i, j, zlev, offs, dims) - flux_velo(-1)/dz 
    end if
  contains
    function flux_velo (itype)
      ! Flux at upper interface (itype=1) or lower interface (itype=-1)
      implicit none
      integer               :: itype
      real(8), dimension(3) :: flux_velo

      integer               :: d
      real(8), dimension(3) ::dz_l
      
      d = dom%id + 1    

      dz_l = 0.5 * (dz + dz_e (dom, i, j, zlev+itype, offs, dims))  ! thickness of layer centred on interface

      flux_velo = itype * eddy_viscosity (dom, i, j, zlev, offs, dims) &
           * (sol(S_VELO,zlev+itype)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - velo(EDGE*id+RT+1:EDGE*id+UP+1)) / dz_l
    end function flux_velo
  end subroutine trend_velo
end module vert_diffusion_mod
