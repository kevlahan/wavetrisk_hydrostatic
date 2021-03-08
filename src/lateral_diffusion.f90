module lateral_diffusion_mod
  use ops_mod
  use lin_solve_mod
  implicit none
contains
  subroutine implicit_lateral_diffusion
    ! Backwards Euler step lateral diffusion of scalars
    use adapt_mod
    implicit none
    integer :: k, v
    
    call update_array_bdry (sol, NONE, 27)

    do k = 1, zlevels
       if (rank==0 .and. log_iter) write (6,'(a,i2,a)') "Level ", k, " lateral diffusion"
       
       if (implicit_diff_sclr) then
          do v = scalars(1), scalars(2)
             call rhs_diffuse (v, k)
             call elliptic_solver (sol(v,k), trend(S_TEMP,1), backwards_euler_diffusion_sclr, "sclr")
          end do
       end if

       if (implicit_diff_divu) then
          call rhs_diffuse (S_VELO, k)
          call elliptic_solver (sol(S_VELO,k), trend(S_VELO,1), backwards_euler_diffusion_divu, "velo")
       end if
    end do
    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE, 27)
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine implicit_lateral_diffusion

  function backwards_euler_diffusion_sclr (q, l)
    ! Calculates linear operator for backwards Euler step for lateral diffusion
    implicit none
    integer                   :: l
    type(Float_Field), target :: backwards_euler_diffusion_sclr, q

    integer :: d, j

    backwards_euler_diffusion_sclr = q

    ! Scalar gradient flux
    do d = 1, size(grid)
       h_flux => horiz_flux(S_MASS)%data(d)%elts
       scalar => q%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
       end do
       nullify (h_flux, scalar)
    end do
    horiz_flux(S_MASS)%bdry_uptodate = .false.
    call update_bdry (horiz_flux(S_MASS), l, 213)

    ! Calculate divergence
    do d = 1, size(grid)
       dscalar => Laplacian_scalar(S_MASS)%data(d)%elts
       h_flux  => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, h_flux)
    end do
    Laplacian_scalar(S_MASS)%bdry_uptodate = .false.
    call update_bdry (Laplacian_scalar(S_MASS), l, 101)

    ! Form complete linear operator 
    do d = 1, size(grid)
       dscalar   => backwards_euler_diffusion_sclr%data(d)%elts
       mass      => q%data(d)%elts
       Laplacian => Laplacian_scalar(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (complete_lo_sclr, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, Laplacian, mass)
    end do
  end function backwards_euler_diffusion_sclr

  subroutine complete_lo_sclr (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id) >= ADJZONE) then
       dscalar(id) = mass(id) + dt*visc_sclr(S_TEMP) * Laplacian(id)
    else
       dscalar(id) = 0.0_8
    end if
  end subroutine complete_lo_sclr

  function backwards_euler_diffusion_divu (q, l)
    ! Calculates linear operator for backwards Euler step for lateral diffusion
    implicit none
    integer                   :: l
    type(Float_Field), target :: backwards_euler_diffusion_divu, q

    integer :: d, j

    backwards_euler_diffusion_divu = q

    ! Divergence of velocity
    do d = 1, size(grid)
       divu => exner_fun(1)%data(d)%elts
       velo => q%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_divu, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (divu, velo)
    end do
    exner_fun(1)%bdry_uptodate = .false.
    call update_bdry (exner_fun(1), l, 101)

    ! Gradient of divergence of velocity
    do d = 1, size(grid)
       divu      => exner_fun(1)%data(d)%elts
       Laplacian => Laplacian_vector(S_ROTU)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (grad_divu, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
       end do
       nullify (divu, Laplacian)
    end do
    Laplacian_vector(S_ROTU)%bdry_uptodate = .false.
    call update_bdry (Laplacian_vector(S_ROTU), l, 101)

    ! Form complete linear operator 
    do d = 1, size(grid)
       dscalar   => backwards_euler_diffusion_divu%data(d)%elts
       velo      => q%data(d)%elts
       Laplacian => Laplacian_vector(S_ROTU)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (complete_lo_divu, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
       end do
       nullify (dscalar, Laplacian, velo)
    end do
  end function backwards_euler_diffusion_divu

   subroutine complete_lo_divu (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e
    
    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) then
          dscalar(id_e) = velo(id_e) - dt*visc_divu * Laplacian(id_e)
       else
          dscalar(id_e) = 0.0_8
       end if
    end do
  end subroutine complete_lo_divu

  subroutine rhs_diffuse (v, k)
    implicit none
    integer :: v, k
    
    integer :: d, j, l
    
    do l = level_end, level_start, -1
       do d = 1, size(grid)
          select case (v)
          case (S_MASS, S_TEMP)
             mass   => sol(v,k)%data(d)%elts
             scalar => trend(S_TEMP,1)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_rhs_diffuse_sclr, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (mass, scalar)
          case (S_VELO)
             velo  => sol(S_VELO,k)%data(d)%elts
             velo1 => trend(S_VELO,1)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_rhs_diffuse_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
             end do
             nullify (velo, velo1)
          end select
       end do
    end do

    select case (v)
    case (S_MASS, S_TEMP)
       trend(S_TEMP,1)%bdry_uptodate = .false.
       call update_bdry (trend(S_TEMP,1), NONE, 600)
    case (S_VELO)
       trend(S_VELO,1)%bdry_uptodate = .false.
       call update_bdry (trend(S_VELO,1), NONE, 600)
    end select
  end subroutine rhs_diffuse

  subroutine cal_rhs_diffuse_sclr (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id) >= ADJZONE) then
       scalar(id) = mass(id)
    else
       scalar(id) = 0.0_8
    end if
  end subroutine cal_rhs_diffuse_sclr
  
  subroutine cal_rhs_diffuse_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e

    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) then
          velo1(id_e) = velo(id_e)
       else
          velo1(id_e) = 0.0_8
       end if
    end do
  end subroutine cal_rhs_diffuse_velo
end module  lateral_diffusion_mod
