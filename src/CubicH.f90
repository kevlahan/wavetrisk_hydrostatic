module CubicH_mod
  ! Computes horizontal pressure gradient using the CubicH algorithm of Shchepetkin and McWilliams (2003)
  ! and adds horizontal gradient of kinetic energy for complete gradient term
  ! (similar to prsgrd32A.F in the ROMS code)
  use utils_mod
  use ops_mod
  implicit none
  real(8), parameter                                     :: eps_small = 1d-33
  real(8), parameter                                     :: OneTwelfth = 1d0/12d0
  type(Float_Field), dimension(:,:), allocatable, target :: q
contains
  subroutine hpg (qq, dq, k)
    ! Compute horizontal pressure gradient at layer k
    ! (need to call pressure_CubicH first to compute pressure over entire grid for all vertical layers)
    implicit none
    integer                                   :: k
    type(Float_Field), dimension(:,:), target :: qq, dq

    integer :: d, l, p

    q = qq ! share variable

    do l = level_end, level_start, -1
       call FC (l, k)
    end do

    ! Compute horizontal pressure gradient at level k over entire grid
    do d = 1, size(grid)
       exner     => exner_fun(k)%data(d)%elts ! pressure
       qe        => grid(d)%qe%elts           ! FC
       ke        => grid(d)%ke%elts           ! kinetic energy
       dvelo     => dq(S_VELO,k)%data(d)%elts
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (du_grad_CubicH, grid(d), p-1, k, 0, 0)
       end do
       nullify (dvelo, exner, ke, qe)
    end do
    dq(S_VELO,k)%bdry_uptodate = .false.
  end subroutine hpg

  subroutine pressure_CubicH (qq)
    ! Compute pressure over entire grid (call before computing horizontal pressure gradient)
    ! pressure is stored in exner_fun(1:zlevels)
    implicit none
    type(Float_Field), dimension(:,:), target :: qq
    
    integer :: d, j, p

    q = qq ! share variable

    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (cal_press_CubicH, grid(d), p-1, z_null, 0, 1)
       end do       
    end do
  end subroutine pressure_CubicH

  subroutine du_grad_CubicH (dom, i, j, zlev, offs, dims)
    ! Computes gradient term using CubicH approximation for horizontal pressure gradient
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, e, id, id_e
    real(8), dimension(1:EDGE) :: gradKE, gradP

    id = idx (i, j, offs, dims)
    d = dom%id + 1

    ! Pressure gradient (not including horizontal gradients of geopotential)
    gradP = gradi_e (exner, dom, i, j, offs, dims)
    
    ! Gradient of kinetic energy
    gradKE = gradi_e (ke, dom, i, j, offs, dims)

    do e = 1, EDGE
       id_e = EDGE*id + e
       dvelo(id_e) = dvelo(id_e)/dom%len%elts(id_e) - gradKE(e) - (gradP(e) + qe(id_e)) / porous_density (d, id+1, zlev)
    end do
  end subroutine du_grad_CubicH

  subroutine cal_press_CubicH (dom, i, j, zlev, offs, dims)
    ! Compute pressure at all vertical layers, stored in exner_fun
    ! (does not require horizontal gradients)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id, k
    real(8)                         :: cff 
    real(8), dimension(1:zlevels)   :: Rk, Zk 
    real(8), dimension(0:zlevels)   :: dRk, dZk

    id = idx (i, j, offs, dims)
    d = dom%id + 1 

    do k = 1, zlevels
       Rk(k)  = density (dom, i, j, k, offs, dims, q)
       Zk(k) =      z_i (dom, i, j, k, offs, dims, q)
    end do

    do k = 1, zlevels-1
       dRk(k) = Rk(k+1) - Rk(k)
       dZk(k) = Zk(k+1) - Zk(k)
    end do

    ! Boundary conditions 
    dRk(0)       = dRk(1)
    dRk(zlevels) = dRk(zlevels-1) 
    dZk(0)       = dZk(1)
    dZk(zlevels) = dZk(zlevels-1)
      
     ! Compute harmonic averages
    do k = zlevels, 1, -1
       dRk(k) = h_av (dRk(k), dRk(k-1))
       dZk(k) = 2d0 * dZk(k) * dZk(k-1) / (dZk(k) + dZk(k-1))
    end do

    ! Compute pressure by integrating down
    exner_fun(zlevels)%data(d)%elts(id+1) = P_top_CubicH ()
    do k = zlevels-1, 1, -1
       exner_fun(k)%data(d)%elts(id+1) = exner_fun(k+1)%data(d)%elts(id+1) &
            + 0.5d0 * grav_accel * ( (Rk(k+1)+Rk(k))  * (Zk(k+1)-Zk(k))  &
            - 0.2d0 * ( (dRk(k+1)-dRk(k)) * ( dZk(k+1)-dZk(k) - OneTwelfth * (dZk(k+1)+dZk(k)) ) &
            - (dZk(k+1)-dZk(k)) * ( dRk(k+1)-dRk(k) - OneTwelfth * (dRk(k+1)+dRk(k)) ) ) )
    end do
  contains
    real(8) function P_top_CubicH ()
      ! Pressure in top layer
      ! (free surface is included in top layer)
      implicit none
      real(8) :: depth, dZk, dRk, Rk_top

      Rk_top = density (dom, i, j, zlevels, offs, dims, q)

      dRk = Rk_top - density (dom, i, j, zlevels-1, offs, dims, q)

      depth = 0.5d0 * dz_i (dom, i, j, zlevels, offs, dims, q)   ! depth of middle of top layer
      dZk   = dz_l (dom, i, j, zlevels-1, offs, dims, q) 

      P_top_CubicH = grav_accel * depth * (Rk_top + 0.5d0 * depth * dRk / dZk)
    end function P_top_CubicH
  end subroutine cal_press_CubicH

  subroutine FC (l, zlev)
    ! Computes variable FC, the geopotential term of the horizontal pressure gradient
    implicit none
    integer :: l, zlev
    integer                   :: d, e, j, p
    ! type(Float_Field), target :: v, D_prime, gradPhi, dZ_prime

    ! ! Initialize local float fields
    ! v        = q(S_VELO,1)
    ! D_prime  = q(S_VELO,1)
    ! dZ_prime = q(S_VELO,1)
    ! gradPhi  = q(S_VELO,1)

    ! do d = 1, size(grid)
    !    ! Compute density differences with results at edges
    !    velo => v%data(d)%elts
    !    do j = 1, grid(d)%lev(l)%length
    !       call apply_onescale_to_patch (cal_diff_density_horiz, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
    !    end do
    !    nullify (velo)
    ! end do
    ! call update_bdry (v, l, 9)

    ! ! Compute horizontal harmonic averages of density differences, D_prime
    ! do d = 1, size(grid)
    !    velo => v%data(d)%elts
    !    dvelo => D_prime%data(d)%elts
    !    do j = 1, grid(d)%lev(l)%length
    !       call apply_onescale_to_patch (harmonic_avgs_horiz, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 1)
    !    end do
    !    nullify (dvelo, velo)
    ! end do
    ! call update_bdry (D_prime, l, 9)

    ! ! Compute z coord differences with results at edges
    !  do d = 1, size(grid)
    !    velo => v%data(d)%elts
    !    do j = 1, grid(d)%lev(l)%length
    !       call apply_onescale_to_patch (cal_dZ_horiz, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
    !    end do
    !    nullify (velo)
    ! end do
    ! call update_bdry (v, l, 9)

    ! ! Compute horizontal harmonic averages of z coord differences, dZ_prime
    ! do d = 1, size(grid)
    !    velo => v%data(d)%elts
    !    dvelo => dZ_prime%data(d)%elts
    !    do j = 1, grid(d)%lev(l)%length
    !       call apply_onescale_to_patch (harmonic_avgs_horiz, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 1)
    !    end do
    !    nullify (dvelo, velo)
    ! end do
    ! call update_bdry (dZ_prime, l, 9)

    ! Compute FC with results stored in qe
    do d = 1, size(grid)
       !velo1 => D_prime%data(d)%elts
       !velo2 => dZ_prime%data(d)%elts
       qe    => grid(d)%qe%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_FC, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
       end do
       !nullify (qe, velo1, velo2)
       nullify (qe)
    end do

    ! Compute FC with results stored in qe
    do d = 1, size(grid)
       qe  => grid(d)%qe%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_FC, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
       end do
       nullify (qe)
    end do
  end subroutine FC

!   subroutine cal_FC (dom, i, j, zlev, offs, dims)
!     ! Computes variable FC (geopotential term of horizontal pressure gradient)
! !!! all horizontal differences should be computed using appropriate restrictions !!!
!     implicit none
!     type(Domain)                   :: dom
!     integer                        :: i, j, zlev
!     integer, dimension(N_BDRY+1)   :: offs
!     integer, dimension(2,N_BDRY+1) :: dims

!     integer                    :: d, e, id, id_i, idE, idNE, idN
!     real(8)                    :: Rho_i, Zi
!     real(8), dimension(1:EDGE) :: dD_prime, dRho, dZ, ddZ_prime, Rho1, Rho_e, sum_D_prime, sum_dZ_prime

!     id = idx (i, j, offs, dims)
!     id_i = id + 1
!     d = dom%id + 1

!     idE  = idx (i+1, j,   offs, dims)
!     idNE = idx (i+1, j+1, offs, dims)
!     idN  = idx (i,   j+1, offs, dims)

!     Rho_i = density (dom, i, j, zlev, offs, dims, q)
!     Zi   =      z_i (dom, i, j, zlev, offs, dims, q)

!     Rho1(RT+1) = density (dom, i+1, j,   zlev, offs, dims, q)
!     Rho1(DG+1) = density (dom, i+1, j+1, zlev, offs, dims, q)
!     Rho1(UP+1) = density (dom, i,   j+1, zlev, offs, dims, q)

!     Rho_e = 0.5d0 * (Rho_i + Rho1)

!     dZ(RT+1) =    z_i (dom, i+1, j,   zlev, offs, dims, q) - Zi
!     dZ(DG+1) = - (z_i (dom, i+1, j+1, zlev, offs, dims, q) - Zi)
!     dZ(UP+1) =    z_i (dom, i,   j+1, zlev, offs, dims, q) - Zi

!     dRho(RT+1) =    Rho1(RT+1) - Rho_i
!     dRho(DG+1) = - (Rho1(DG+1) - Rho_i)
!     dRho(UP+1) =    Rho1(UP+1) - Rho_i

!     sum_D_prime(RT+1) = velo1(EDGE*id+RT+1) + velo1(EDGE*idE+RT+1)
!     sum_D_prime(DG+1) = velo1(EDGE*id+DG+1) + velo1(EDGE*idNE+DG+1)
!     sum_D_prime(UP+1) = velo1(EDGE*id+UP+1) + velo1(EDGE*idN+UP+1)

!     sum_dZ_prime(RT+1) = velo2(EDGE*id+RT+1) + velo2(EDGE*idE+RT+1)
!     sum_dZ_prime(DG+1) = velo2(EDGE*id+DG+1) + velo2(EDGE*idNE+DG+1)
!     sum_dZ_prime(UP+1) = velo2(EDGE*id+UP+1) + velo2(EDGE*idN+UP+1)

!     dD_prime(RT+1) =    velo1(EDGE*idE+RT+1)  - velo1(EDGE*id+RT+1)
!     dD_prime(DG+1) = - (velo1(EDGE*idNE+DG+1) - velo1(EDGE*id+DG+1))
!     dD_prime(UP+1) =    velo1(EDGE*idN+UP+1)  - velo1(EDGE*id+UP+1)

!     ddZ_prime(RT+1) =    velo2(EDGE*idE+RT+1)  - velo2(EDGE*id+RT+1)
!     ddZ_prime(DG+1) = - (velo2(EDGE*idNE+DG+1) - velo2(EDGE*id+DG+1))
!     ddZ_prime(UP+1) =    velo2(EDGE*idN+UP+1)  - velo2(EDGE*id+UP+1)

!     ! qe(EDGE*id+RT+1:EDGE*id+UP+1) = grav_accel * (Rho_e * dZ &
!     !     - 0.1d0 * (dD_prime * (dZ   - OneTwelfth * sum_dZ_prime) &
!     !             - ddZ_prime * (dRho - OneTwelfth * sum_D_prime)) &
!     !             ) / dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1)

!     qe(EDGE*id+RT+1:EDGE*id+UP+1) = grav_accel * Rho_e * dZ / dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1)
!    end subroutine cal_FC

  subroutine cal_FC (dom, i, j, zlev, offs, dims)
    ! Computes variable FC (geopotential term of horizontal pressure gradient)
!!! all horizontal differences should be computed using appropriate restrictions !!!
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                      :: d, e, id, id_e, idE, idNE, idN, idW, idSW, idS, id2E, id2NE, id2N
    real(8), dimension(0:3*EDGE) :: dRx, dZx, Rx, Zx
    real(8), dimension(EDGE,0:1) :: ddRx, ddZx 

    id = idx (i, j, offs, dims)
    d = dom%id + 1
    
    Zx(0)         = z_i (dom, i,   j,   zlev, offs, dims, q)
    Zx(RT+1)      = z_i (dom, i+1, j,   zlev, offs, dims, q)
    Zx(DG+1)      = z_i (dom, i+1, j+1, zlev, offs, dims, q)
    Zx(UP+1)      = z_i (dom, i,   j+1, zlev, offs, dims, q)

    Zx(RT+1+EDGE) = z_i (dom, i-1, j,   zlev, offs, dims, q)
    Zx(DG+1+EDGE) = z_i (dom, i-1, j-1, zlev, offs, dims, q)
    Zx(UP+1+EDGE) = z_i (dom, i,   j-1, zlev, offs, dims, q)

    Zx(RT+1+2*EDGE) = z_i (dom, i+2, j,   zlev, offs, dims, q)
    Zx(DG+1+2*EDGE) = z_i (dom, i-2, j-2, zlev, offs, dims, q)
    Zx(UP+1+2*EDGE) = z_i (dom, i,   j+2, zlev, offs, dims, q)
  
    Rx(0)         = density (dom, i,   j,   zlev, offs, dims, q)
    Rx(RT+1)      = density (dom, i+1, j,   zlev, offs, dims, q)
    Rx(DG+1)      = density (dom, i+1, j+1, zlev, offs, dims, q)
    Rx(UP+1)      = density (dom, i,   j+1, zlev, offs, dims, q)

    Rx(RT+1+EDGE) = density (dom, i-1, j,   zlev, offs, dims, q)
    Rx(DG+1+EDGE) = density (dom, i-1, j-1, zlev, offs, dims, q)
    Rx(UP+1+EDGE) = density (dom, i,   j-1, zlev, offs, dims, q)

    Rx(RT+1+2*EDGE) = density (dom, i+2, j,   zlev, offs, dims, q)
    Rx(DG+1+2*EDGE) = density (dom, i-2, j-2, zlev, offs, dims, q)
    Rx(UP+1+2*EDGE) = density (dom, i,   j+2, zlev, offs, dims, q)

    dRx(RT+1)        =   Rx(RT+1) - Rx(0)
    dRx(DG+1)        = -(Rx(DG+1) - Rx(0))
    dRx(UP+1)        =   Rx(UP+1) - Rx(0)

    dRx(RT+1+EDGE)   = -(Rx(RT+1+EDGE) - Rx(0))
    dRx(DG+1+EDGE)   =   Rx(DG+1+EDGE) - Rx(0)
    dRx(UP+1+EDGE)   = -(Rx(UP+1+EDGE) - Rx(0))

    dRx(RT+1+2*EDGE) = Rx(RT+1+2*EDGE) - Rx(RT+1)
    dRx(DG+1+2*EDGE) = Rx(DG+1+2*EDGE) - Rx(DG+1+EDGE)
    dRx(UP+1+2*EDGE) = Rx(UP+1+2*EDGE) - Rx(UP+1)
    
    dZx(RT+1)        =   Zx(RT+1) - Zx(0)
    dZx(DG+1)        = -(Zx(DG+1) - Zx(0))
    dZx(UP+1)        =   Zx(UP+1) - Zx(0)

    dZx(RT+1+EDGE)   = -(Zx(RT+1+EDGE) - Zx(0))
    dZx(DG+1+EDGE)   =   Zx(DG+1+EDGE) - Zx(0)
    dZx(UP+1+EDGE)   = -(Zx(UP+1+EDGE) - Zx(0))

    dZx(RT+1+2*EDGE) = Zx(RT+1+2*EDGE) - Zx(RT+1)
    dZx(DG+1+2*EDGE) = Zx(DG+1+2*EDGE) - Zx(DG+1+EDGE)
    dZx(UP+1+2*EDGE) = Zx(UP+1+2*EDGE) - Zx(UP+1)

    do e = 1, EDGE
       ddRx(e,0) = h_av (dRx(e), dRx(e+EDGE))
       ddZx(e,0) = h_av (dZx(e), dZx(e+EDGE))
       
       ddRx(e,1) = h_av (dRx(e+EDGE), dRx(e+2*EDGE))
       ddZx(e,1) = h_av (dZx(e+EDGE), dZx(e+2*EDGE))
    end do
    
    do e = 1, EDGE
       id_e = EDGE*id + e
       qe(id_e) = 0.5d0 * grav_accel * ( (Rx(e)+Rx(0)) * dZx(e) &
            - 0.2d0 * ( (ddRx(e,1)-ddRx(e,0)) * ( dZx(e) - OneTwelfth * (ddZx(e,1)+ddZx(e,0)) ) &
                     -  (ddZx(e,1)-ddZx(e,0)) * ( dRx(e) - OneTwelfth * (ddRx(e,1)+ddRx(e,0)) ) ) &
                     ) / dom%len%elts(id_e)
    end do
  end subroutine cal_FC

  Real(8) function h_av (f1, f2)
    ! Harmonic average of f1 and f2
    implicit none
    real(8) :: f1, f2

    real(8) :: cff

    cff = 2d0 * f1 * f2
    if (cff > eps_small) then
       h_av = cff / (f1 + f2)
    else
       h_av = 0d0
    end if
  end function h_av

  subroutine cal_diff_density_horiz (dom, i, j, zlev, offs, dims)
    ! Differences of density, with output at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: Rho_i

    id = idx (i, j, offs, dims)

    Rho_i = density (dom, i, j, zlev, offs, dims, q)

    velo(EDGE*id+RT+1) =    density (dom, i+1, j,   zlev, offs, dims, q) - Rho_i
    velo(EDGE*id+DG+1) = - (density (dom, i+1, j+1, zlev, offs, dims, q) - Rho_i)
    velo(EDGE*id+UP+1) =    density (dom, i,   j+1, zlev, offs, dims, q) - Rho_i
  end subroutine cal_diff_density_horiz
  
  subroutine cal_dZ_horiz (dom, i, j, zlev, offs, dims)
    ! Differences of z coordinate, with output at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: Zi

    id = idx (i, j, offs, dims)

    Zi = z_i (dom, i, j, zlev, offs, dims, q)

    velo(EDGE*id+RT+1) =    z_i (dom, i+1, j,   zlev, offs, dims, q) - Zi
    velo(EDGE*id+DG+1) = - (z_i (dom, i+1, j+1, zlev, offs, dims, q) - Zi)
    velo(EDGE*id+UP+1) =    z_i (dom, i,   j+1, zlev, offs, dims, q) - Zi
  end subroutine cal_dZ_horiz

  subroutine harmonic_avgs_horiz (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id, idS, idSW, idW
    real(8), dimension(1:EDGE) :: cff

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    
    cff(RT+1) = 2d0 * velo(EDGE*id+RT+1) * velo(EDGE*idW+RT+1)
    cff(DG+1) = 2d0 * velo(EDGE*id+DG+1) * velo(EDGE*idSW+DG+1)
    cff(UP+1) = 2d0 * velo(EDGE*id+UP+1) * velo(EDGE*idS+UP+1)

    dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    if (cff(RT+1) > eps_small) dvelo(EDGE*id+RT+1) = cff(RT+1) / (velo(EDGE*id+RT+1) + velo(EDGE*idW+RT+1))
    if (cff(DG+1) > eps_small) dvelo(EDGE*id+DG+1) = cff(DG+1) / (velo(EDGE*id+DG+1) + velo(EDGE*idSW+DG+1))
    if (cff(UP+1) > eps_small) dvelo(EDGE*id+UP+1) = cff(UP+1) / (velo(EDGE*id+UP+1) + velo(EDGE*idS+UP+1))
  end subroutine harmonic_avgs_horiz
end module CubicH_mod
