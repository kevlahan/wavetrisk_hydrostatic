module multi_level_mod
  use ops_mod
  implicit none
contains
  subroutine trend_ml (q, dq)
    ! Compute trends of prognostic variables assuming Lagrangian vertical coordinates
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq

    integer :: k, l, v

    call update_bdry (q, NONE, 967)

    ! Initialize trends
    do k = 1, zlevels
       do v = scalars(1), scalars(2)
          call zero_float_field (dq(v,k), AT_NODE)
       end do
       call zero_float_field (dq(S_VELO,k), AT_EDGE)
    end do

    ! Compute surface pressure on all grids
    call cal_surf_press (q(1:N_VARIABLE,1:zlevels))

    ! Compute each vertical level starting from surface
    do k = 1, zlevels
       if (Laplace_divu /= 0) call cal_divu_ml (q(S_VELO,k))
       if (Laplace_sclr == 2) call cal_Laplacian_scalars (q, k)
       if (Laplace_divu == 2) call cal_Laplacian_divu ! requires divu

       ! Calculate trend on all scales, from fine to coarse
       do l = level_end, level_start, -1
          ! Finish non-blocking communication of dq from level (l+1)
          if (l < level_end) call update_bdry__finish (dq(scalars(1):scalars(2),k),l+1) 

          call basic_operators  (q, dq, k, l)
          call cal_scalar_trend (q, dq, k, l)

          ! Start non-blocking communication of dq for use at next level (l-1)
          if (level_start /= level_end .and. l > level_start) call update_bdry__start (dq(scalars(1):scalars(2),k),l) 

          call velocity_trend_source (q, dq, k, l)
       end do
       call velocity_trend_grad (q, dq, k)
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_ml

  subroutine basic_operators (q, dq, k, l)
    ! Evaluates basic operators on grid level l and computes/restricts Bernoulli, Exner and fluxes
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq
    integer                                                      :: k, l

    integer :: d, j, v

    do d = 1, size(grid)
       mass      => q(S_MASS,k)%data(d)%elts
       temp      => q(S_TEMP,k)%data(d)%elts
       velo      => q(S_VELO,k)%data(d)%elts
       mean_m    => sol_mean(S_MASS,k)%data(d)%elts
       mean_t    => sol_mean(S_TEMP,k)%data(d)%elts
       exner     => exner_fun(k)%data(d)%elts
       bernoulli => grid(d)%bernoulli%elts
       ke        => grid(d)%ke%elts
       vort      => grid(d)%vort%elts
       qe        => grid(d)%qe%elts
       
       ! Compute horizontal fluxes, potential vorticity (qe), Bernoulli, Exner (incompressible case) etc
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          call step1 (dq, q, grid(d), grid(d)%lev(l)%elts(j), k, 0)
       end do
       call apply_to_penta_d (post_step1, grid(d), l, z_null)
       nullify (mass, velo, temp, mean_m, mean_t, ke, qe, vort)

       ! Compute or restrict Bernoulli, Exner and fluxes
       if (l < level_end) then
          scalar => grid(d)%bernoulli%elts
          call cpt_or_restr_scalar (grid(d), l)
          nullify (scalar)

          scalar => exner_fun(k)%data(d)%elts
          call cpt_or_restr_scalar (grid(d), l)
          nullify (scalar)
          
          do v = scalars(1), scalars(2)
             dscalar => dq(v,k)%data(d)%elts
             h_flux  => horiz_flux(v)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) using dscalar (l+1)
             nullify (dscalar, h_flux)
          end do
       end if
       nullify (bernoulli, exner)
    end do
    horiz_flux%bdry_uptodate = .false.
    if (level_start /= level_end) call update_bdry (horiz_flux, l, 968)

    if (Laplace_rotu == 2) call cal_Laplacian_vector_rot (l) ! requires vorticity
  end subroutine basic_operators

  subroutine cal_scalar_trend (q, dq, k, l)
    ! Evaluate scalar trends at level l
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq
    integer                                                      :: k, l

    integer :: d, j, v

    call update_bdry (horiz_flux, l, 969)
    
    do d = 1, size(grid)
       do v = scalars(1), scalars(2)
          dscalar => dq(v,k)%data(d)%elts
          h_flux  => horiz_flux(v)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
    end do
    dq(S_MASS:S_TEMP,k)%bdry_uptodate = .false.
  end subroutine cal_scalar_trend

  subroutine velocity_trend_source (q, dq, k, l)
    ! Evaluate source part of velocity trends at level l
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq
    integer :: k, l

    integer :: d, j

    u_source => du_source
   
    do d = 1, size(grid)
       mass    => q(S_MASS,k)%data(d)%elts
       velo    => q(S_VELO,k)%data(d)%elts
       mean_m  => sol_mean(S_MASS,k)%data(d)%elts
       dvelo   => dq(S_VELO,k)%data(d)%elts
       h_mflux => horiz_flux(S_MASS)%data(d)%elts
       ke      => grid(d)%ke%elts
       qe      => grid(d)%qe%elts
       vort    => grid(d)%vort%elts
       
       if (Laplace_divu == 2) then
          divu => Laplacian_vector(S_DIVU)%data(d)%elts
       else
          divu => grid(d)%divu%elts
       end if

       if (l < level_end) then
          call cpt_or_restr_u_source (grid(d), k, l)
       else
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (u_source, grid(d), grid(d)%lev(level_end)%elts(j), k, 0, 0)
          end do
       end if

       nullify (mass, velo, mean_m, dvelo, h_mflux, divu, ke, qe, vort)
    end do
    dq(S_VELO,k)%bdry_uptodate = .false.

    nullify (u_source)
  end subroutine velocity_trend_source

  subroutine velocity_trend_grad (q, dq, k)
    ! Evaluate complete velocity trend by adding gradient terms to previously calculated source terms on entire grid
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq
    integer                                                      :: k

    integer :: d, p

    do d = 1, size(grid)
       mass      => q(S_MASS,k)%data(d)%elts
       temp      => q(S_TEMP,k)%data(d)%elts
       mean_m    => sol_mean(S_MASS,k)%data(d)%elts
       mean_t    => sol_mean(S_TEMP,k)%data(d)%elts
       dvelo     => dq(S_VELO,k)%data(d)%elts
       exner     => exner_fun(k)%data(d)%elts
       bernoulli => grid(d)%bernoulli%elts
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (du_grad, grid(d), p-1, k, 0, 0)
       end do
       nullify (mass, temp, mean_m, mean_t, dvelo, exner, bernoulli)
    end do
    dq(S_VELO,k)%bdry_uptodate = .false.
  end subroutine velocity_trend_grad
     
  subroutine cal_Laplacian_scalars (q, k)
    ! Computes Laplacian of scalars q, div(grad q)
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q
    integer :: k
    
    integer :: d, j, l, v

    call update_bdry (q(scalars(1):scalars(2),k), NONE, 970)
    
    do l = level_end, level_start, -1
       ! Compute scalar fluxes
       do d = 1, size(grid)
          do v = scalars(1), scalars(2)
             scalar => q(v,k)%data(d)%elts
             h_flux => horiz_flux(v)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
             end do
             nullify (scalar, h_flux)
          end do

          ! Compute or restrict fluxes
          if (l < level_end) then
             do v = scalars(1), scalars(2)
                dscalar => Laplacian_scalar(v)%data(d)%elts
                h_flux  => horiz_flux(v)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) using dscalar (l+1)
                nullify (dscalar, h_flux)
             end do
          end if
       end do
       horiz_flux%bdry_uptodate = .false.
       call update_bdry (horiz_flux, l, 971)

       do d = 1, size(grid)
          do v = scalars(1), scalars(2)
             dscalar => Laplacian_scalar(v)%data(d)%elts
             h_flux  => horiz_flux(v)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do
       end do
       Laplacian_scalar%bdry_uptodate = .false.
       call update_bdry (Laplacian_scalar, l, 972)
    end do
  end subroutine cal_Laplacian_scalars

  subroutine cal_Laplacian_vector_rot (l)
    ! Computes rot(rot(vorticity)) needed for second-order vector Laplacian
    implicit none
    integer :: l
    
    integer :: d, j

    ! Compute rot(vorticity)
    do d = 1, size(grid)
       vort      => grid(d)%vort%elts
       Laplacian => Laplacian_vector(S_ROTU)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_Laplacian_rotu, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
       end do
       nullify (vort, Laplacian)
    end do
    Laplacian_vector(S_ROTU)%bdry_uptodate = .false.
    call update_bdry (Laplacian_vector(S_ROTU), l, 973)
        
    ! Compute rot(rot(vorticity)) using previous result for rot(vorticity)
    !!! grid(d)%vort is now rot(rot(vorticity)), not vorticity !!!
    do d = 1, size(grid)
       velo => Laplacian_vector(S_ROTU)%data(d)%elts
       vort => grid(d)%vort%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=8)
       end do
       call apply_to_penta_d (post_vort, grid(d), l, z_null)
       nullify (velo, vort)
    end do
  end subroutine cal_Laplacian_vector_rot
  
  subroutine cal_Laplacian_divu
    ! Computes Laplacian of divu, div(grad divu)
    implicit none
    integer :: d, j, l, v

    do l = level_end, level_start, -1
       ! Compute scalar fluxes
       do d = 1, size(grid)
          scalar => grid(d)%divu%elts
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
          end do
          nullify (scalar, h_flux)

          ! Compute or restrict fluxes
          if (l < level_end) then
             dscalar => Laplacian_vector(S_DIVU)%data(d)%elts
             h_flux  => horiz_flux(S_MASS)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) using dscalar (l+1)
             nullify (dscalar, h_flux)
          end if
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       call update_bdry (horiz_flux(S_MASS), l, 974)

       do d = 1, size(grid)
          dscalar => Laplacian_vector(S_DIVU)%data(d)%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
       Laplacian_vector(S_DIVU)%bdry_uptodate = .false.
       call update_bdry (Laplacian_vector(S_DIVU), l, 975)
    end do
  end subroutine cal_Laplacian_divu
 
  subroutine cpt_or_restr_scalar (dom, l)
    ! Restrict scalar if possible for grad(scalar) computation
    implicit none
    type(Domain) :: dom
    integer      :: l

    integer                     :: j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .false.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd > 0) restrict(c) = .true.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) then
             call apply_interscale_to_patch3 (scalar_cpt_restr, dom, p_par, c, z_null, 0, 1)
          end if
       end do
    end do
  end subroutine cpt_or_restr_scalar

  subroutine scalar_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) >= RESTRCT) scalar(id_par+1) = scalar(id_chd+1)
  end subroutine scalar_cpt_restr

  subroutine cpt_or_restr_u_source (dom, zlev, l)
    ! Restrict velocity source if possible term u_source(velo)
    ! u_source is a pointer function
    implicit none
    type(Domain) :: dom
    integer      :: zlev, l

    integer :: c, j, p_par, p_chd

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd == 0) call apply_onescale_to_patch (u_source, dom, p_par, zlev, 0, 0)
       end do
       call apply_interscale_to_patch (u_source_cpt_restr, dom, dom%lev(l)%elts(j), zlev, 0, 0)
    end do
  end subroutine cpt_or_restr_u_source

  subroutine u_source_cpt_restr (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                     :: dom
    integer                          :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_par, dims_chd

    integer :: id_par, id_chd, idE_chd, idNE_chd, idN_chd

    id_par = idx (i_par, j_par, offs_par, dims_par)

    id_chd   = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (minval(dom%mask_e%elts(EDGE*id_chd+RT+1:EDGE*id_chd+UP+1)) < ADJZONE) &
         call u_source (dom, i_par, j_par, zlev, offs_par, dims_par)

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) &
         dvelo(EDGE*id_par+RT+1) = dvelo(EDGE*id_chd+RT+1) + dvelo(EDGE*idE_chd+RT+1)

    if (dom%mask_e%elts(EDGE*id_chd+DG+1) >= ADJZONE) &
         dvelo(EDGE*id_par+DG+1) = dvelo(EDGE*id_chd+DG+1) + dvelo(EDGE*idNE_chd+DG+1)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) &
         dvelo(EDGE*id_par+UP+1) = dvelo(EDGE*id_chd+UP+1) + dvelo(EDGE*idN_chd+UP+1)
  end subroutine u_source_cpt_restr

  subroutine cpt_or_restr_flux (dom, l)
    ! Restrict flux if possible for dscalar = div(h_flux) computation
    ! requires dscalar = div(h_flux) in addition to h_flux
    implicit none
    type(Domain) :: dom
    integer      :: l

    integer                     :: j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .false.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd > 0) restrict(c) = .true.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) call apply_interscale_to_patch3 (flux_cpt_restr, dom, p_par, c, z_null, 0, 1)
       end do
    end do
  end subroutine cpt_or_restr_flux

  subroutine flux_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute flux restriction by summing coarse, corrective and small fluxes
    implicit none
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: id_par
    real(8), dimension(4) :: sm_flux

    if (i_chd >= PATCH_SIZE .or. j_chd >= PATCH_SIZE) return

    id_par = idx (i_par, j_par, offs_par, dims_par)
    
    if (maxval(dom%mask_e%elts(EDGE*id_par+RT+1:EDGE*id_par+UP+1)) >= RESTRCT) &
         sm_flux = interp_flux (dom, i_chd, j_chd, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_par+RT+1) >= RESTRCT) h_flux(EDGE*id_par+RT+1) = &
         complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, RT, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_par+DG+1) >= RESTRCT) h_flux(EDGE*id_par+DG+1) = &
         complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, DG, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_par+UP+1) >= RESTRCT) h_flux(EDGE*id_par+UP+1) = &
         complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, UP, offs_chd, dims_chd)
  end subroutine flux_cpt_restr

  function interp_flux (dom, i, j, offs, dims)
    implicit none
    real(8), dimension(4)          :: interp_flux
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer, dimension(20) :: id
    real(8)                :: wgt
    
    call get_indices (dom, i+1, j, RT, offs, dims, id)

    interp_flux(1) = - sum (h_flux(id((/WPM,UZM,VMM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1, j-2, offs, dims)+1)%enc) &
         - sum ((h_flux(id((/VPM,WMMM,UMZ/)+1)+1) - h_flux(id((/UPZ,VPMM,WMM/)+1)+1)) * &
         dom%R_F_wgt%elts(idx(i+1, j-1, offs, dims)+1)%enc) ! UPLT S

    interp_flux(2) = sum (h_flux(id((/WMP,UZP,VPP/)+1)+1)* dom%R_F_wgt%elts(idx(i, j, offs, dims)+1)%enc) &
         + sum ((h_flux(id((/VMP,WPPP,UPZ/)+1)+1) - h_flux(id((/UMZ,VMPP,WPP/)+1)+1))* &
         dom%R_F_wgt%elts(idx(i, j+1, offs, dims)+1)%enc) ! LORT

    call get_indices (dom, i, j+1, UP, offs, dims, id)

    interp_flux(3) = - sum (h_flux(id((/UZM,VMM,WPM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1, j, offs, dims)+1)%enc) &
         - sum((h_flux(id((/WMMM,UMZ,VPM/)+1)+1) - h_flux(id((/VPMM,WMM,UPZ/)+1)+1))* &
         dom%R_F_wgt%elts(idx(i+1, j+1, offs, dims)+1)%enc) ! UPLT

    interp_flux(4) = sum (h_flux(id((/UZP,VPP,WMP/)+1)+1) * dom%R_F_wgt%elts(idx(i-2, j, offs, dims)+1)%enc) &
         + sum ((h_flux(id((/WPPP,UPZ,VMP/)+1)+1) - h_flux(id((/VMPP,WPP,UMZ/)+1)+1))* &
         dom%R_F_wgt%elts(idx(i-2, j+1, offs, dims)+1)%enc) ! LORT W
  end function interp_flux

  real(8) function complete_coarse_flux (sm_flux, dom, i_par, j_par, i_chd, j_chd, e, offs_chd, dims_chd)
    implicit none
    real(8), dimension(4)          :: sm_flux
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    
    integer :: e
    real(8) :: p_flux, c_flux

    if (e == RT) then
       p_flux = part_coarse_flux (dom, i_chd+1, j_chd, RT, offs_chd, dims_chd)
       c_flux = coarse_flux (dom, i_par, j_par, i_chd+1, j_chd, offs_chd, dims_chd, RT)
       complete_coarse_flux = p_flux + c_flux + sm_flux(1) + sm_flux(2)
    elseif (e == DG) then
       p_flux = part_coarse_flux (dom, i_chd+1, j_chd+1, DG, offs_chd, dims_chd)
       c_flux = coarse_flux (dom, i_par, j_par, i_chd+1, j_chd+1, offs_chd, dims_chd, DG)
       complete_coarse_flux = p_flux + c_flux + sm_flux(2) + sm_flux(3)
    elseif (e == UP) then
       p_flux = part_coarse_flux (dom, i_chd, j_chd+1, UP, offs_chd, dims_chd)
       c_flux = coarse_flux (dom, i_par, j_par, i_chd, j_chd+1, offs_chd, dims_chd, UP)
       complete_coarse_flux = p_flux + c_flux + sm_flux(3) + sm_flux(4)
    end if
  end function complete_coarse_flux

  real(8) function part_coarse_flux (dom, i, j, e, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, e
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer, dimension(20) :: id
    real(8), dimension(2)  :: area
    real(8), dimension(4)  :: ol_area

    call get_indices (dom, i, j, e, offs, dims, id)

    area = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%a(1:2)

    ol_area(1:2) = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%split
    ol_area(3:4) = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%a(3:4) - ol_area(1:2)

    area(1) = area(1) + ol_area(1) + ol_area(4)
    area(2) = area(2) + ol_area(2) + ol_area(3)
    area = area / sum(area)

    ol_area(1) = dom%overl_areas%elts(id(PP+1)+1)%split(1)
    ol_area(2) = dom%overl_areas%elts(id(MM+1)+1)%split(2)
    ol_area(3) = dom%overl_areas%elts(id(MP+1)+1)%a(3) - dom%overl_areas%elts(id(MP+1)+1)%split(1)
    ol_area(4) = dom%overl_areas%elts(id(PM+1)+1)%a(4) - dom%overl_areas%elts(id(PM+1)+1)%split(2)

    part_coarse_flux = sum (h_flux(id((/UPZ,UMZ/)+1)+1)*area) - sum (h_flux(id((/VMM,WMP/)+1)+1))*area(2) &
         - sum (h_flux(id((/WPM,VPP/)+1)+1))*area(1) &
         + ol_area(3) * dscalar(id(MP+1)+1) - ol_area(4) * dscalar(id(PM+1)+1)  &
         - ol_area(1) * dscalar(id(PP+1)+1) + ol_area(2) * dscalar(id(MM+1)+1)
  end function part_coarse_flux

  real(8) function coarse_flux (dom, i_par, j_par, i_chd, j_chd, offs_chd, dims_chd, e)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, e
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd

    integer :: id, id_mz, id_pz, id_mp,  id_pp, id_pm, id_mm, id_mm2, id_pm2, id_pp2,  id_mp2

    id_mz = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 1 + 1), offs_chd, dims_chd)
    id_pz = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 4 + 1), offs_chd, dims_chd)
    id_mp = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 0 + 1), offs_chd, dims_chd)
    id_pp = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 5 + 1), offs_chd, dims_chd)
    id_pm = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 3 + 1), offs_chd, dims_chd)
    id_mm = idx2 (i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 2 + 1), offs_chd, dims_chd)

    id_mm2 = idx2 (i_chd, j_chd, bfly_no2(:,1,e+1), offs_chd, dims_chd)
    id_pm2 = idx2 (i_chd, j_chd, bfly_no2(:,2,e+1), offs_chd, dims_chd)
    id_pp2 = idx2 (i_chd, j_chd, bfly_no2(:,3,e+1), offs_chd, dims_chd)
    id_mp2 = idx2 (i_chd, j_chd, bfly_no2(:,4,e+1), offs_chd, dims_chd)

    id = idx (i_chd, j_chd, offs_chd, dims_chd)

    coarse_flux = (dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
         + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
         + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
         + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
         + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
         * (dscalar(id_pz+1) - dscalar(id_mz+1)) + &
         dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
         * 0.5d0 * (dscalar(id_pp2+1) - dscalar(id_mz+1)) + &
         dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
         * 0.5d0 * (dscalar(id_pm2+1) - dscalar(id_mz+1)) + &
         dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
         * 0.5d0 * (dscalar(id_pz+1) - dscalar(id_mp2+1)) + &
         dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
         * 0.5d0 * (dscalar(id_pz+1) - dscalar(id_mm2+1))
  end function coarse_flux

  subroutine cal_divu_ml (q)
    ! Returns flux divergence of velocity in divF using solution q, stored in dom%divu
    implicit none
    type(Float_Field), target :: q
    
    integer :: d, j, l

    call update_bdry (q, NONE, 976)

    do l = level_end, level_start, -1
       ! Calculate velocity flux
       do d = 1, size(grid)
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          velo   => q%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=5)
          end do
          nullify (velo)
          if (l < level_end) then
             dscalar => grid(d)%divu%elts
             call cpt_or_restr_flux (grid(d), l) ! restrict flux if possible
             nullify (dscalar)
          end if
          nullify (h_flux)
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       call update_bdry (horiz_flux(S_MASS), l, 977)

       ! Calculate divergence of vertically integrated velocity flux
       do d = 1, size(grid)
          dscalar => grid(d)%divu%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
    end do
  end subroutine cal_divu_ml

  subroutine get_indices (dom, i, j, e, offs, dims, id)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, e
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer, dimension(20)         :: id
    
    integer, dimension(2)  :: ij_mp, ij_pp, ij_pm, ij_mm

    id(UMZ+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims)
    id(UPZ+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims)
    id(WMP+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims)
    id(VPP+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims)
    id(WPM+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims)
    id(VMM+1) = ed_idx (i, j, hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims)

    ij_mp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 0 + 1)
    id(MP+1) = idx(ij_mp(1), ij_mp(2), offs, dims)

    ij_pp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 5 + 1)
    id(PP+1) = idx(ij_pp(1), ij_pp(2), offs, dims)

    ij_pm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 3 + 1)
    id(PM+1) = idx(ij_pm(1), ij_pm(2), offs, dims)

    ij_mm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 2 + 1)

    id(MM+1)   = idx (ij_mm(1), ij_mm(2), offs, dims)

    id(VMP+1)  = ed_idx (ij_mp(1), ij_mp(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
    id(VMPP+1) = ed_idx (ij_mp(1), ij_mp(2), hex_sides (:, hex_s_offs(e+1) + 1  + 4 + 1), offs, dims)
    id(UZP+1)  = ed_idx (ij_mp(1), ij_mp(2), hex_sides (:, hex_s_offs(e+1) + 0  + 4 + 1), offs, dims)
    id(WPPP+1) = ed_idx (ij_pp(1), ij_pp(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
    id(WPP+1)  = ed_idx (ij_pp(1), ij_pp(2), hex_sides (:, hex_s_offs(e+1) + 1  + 2 + 1), offs, dims)
    id(VPM+1)  = ed_idx (ij_pm(1), ij_pm(2), hex_sides (:, hex_s_offs(e+1) + 1  + 4 + 1), offs, dims)
    id(VPMM+1) = ed_idx (ij_pm(1), ij_pm(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
    id(UZM+1)  = ed_idx (ij_pm(1), ij_pm(2), hex_sides (:,(hex_s_offs(e+1) + 3) - 2 + 1), offs, dims)
    id(WMMM+1) = ed_idx (ij_mm(1), ij_mm(2), hex_sides (:, hex_s_offs(e+1) + 1  + 2 + 1), offs, dims)
    id(WMM+1)  = ed_idx (ij_mm(1), ij_mm(2), hex_sides (:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
  end subroutine get_indices
end module multi_level_mod
