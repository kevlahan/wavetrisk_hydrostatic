module multi_level_mod
  use comm_mod
  use ops_mod
  use wavelet_mod
  use refine_patch_mod
  use comm_mpi_mod
  implicit none
contains
  subroutine trend_ml (q, dq, itype)
    ! Compute trends of prognostic variables assuming Lagrangian vertical coordinates
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q, dq
    integer                                                       :: itype

    integer :: d, j, k, l, p

    call update_array_bdry (q, NONE)

    ! First integrate pressure down across all grid points in order to compute surface pressure
    do k = zlevels, 1, -1
       do d = 1, size(grid)
          mass => q(S_MASS,k)%data(d)%elts
          temp => q(S_TEMP,k)%data(d)%elts

          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (integrate_pressure_down, grid(d), p-1, k, 0, 1)
          end do

          nullify (mass, temp)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate trend on finest scale !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute each vertical level starting from surface
    do k = 1, zlevels
       do d = 1, size(grid)
          mass      => q(S_MASS,k)%data(d)%elts
          temp      => q(S_TEMP,k)%data(d)%elts
          velo      => q(S_VELO,k)%data(d)%elts
          h_mflux   => horiz_flux(S_MASS)%data(d)%elts
          h_tflux   => horiz_flux(S_TEMP)%data(d)%elts
          exner     => exner_fun(k)%data(d)%elts
          bernoulli => grid(d)%bernoulli%elts
          divu      => grid(d)%divu%elts
          vort      => grid(d)%vort%elts
          qe        => grid(d)%qe%elts

          ! Compute pressure, geopotential, Exner (compressible case), specific volume
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(level_end)%elts(j), k, 0, 1)
          end do

          ! Compute horizontal fluxes, potential vorticity (qe), Bernoulli, Exner (incompressible case)
          do j = 1, grid(d)%lev(level_end)%length
             call step1 (grid(d), grid(d)%lev(level_end)%elts(j), k)
          end do
          call apply_to_penta_d (post_step1, grid(d), level_end, k)

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (cal_divu, grid(d), grid(d)%lev(level_end)%elts(j), z_null,  0, 1)
          end do
          
          nullify (mass, velo, temp, h_mflux, h_tflux, bernoulli, divu, exner, vort, qe)
       end do

       if (level_start .lt. level_end) call update_vector_bdry__start (horiz_flux, level_end) ! <= comm flux (Jmax)

       ! Compute scalar trends at finest level
       do d = 1, size(grid)
          dmass   => dq(S_MASS,k)%data(d)%elts
          dtemp   => dq(S_TEMP,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS)%data(d)%elts
          h_tflux => horiz_flux(S_TEMP)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 1)
          end do

          nullify (dmass, dtemp, h_mflux, h_tflux)
       end do

       dq(:,k)%bdry_uptodate    = .False.
       horiz_flux%bdry_uptodate = .False.
       if (level_start .lt. level_end) then
          call update_vector_bdry__finish (horiz_flux, level_end) ! <= finish non-blocking communicate mass flux (Jmax)
          call update_vector_bdry__start (dq(S_MASS:S_TEMP,k), level_end) ! <= start non-blocking communicate dmass (l+1)
       end if

       ! Velocity trend, source part
       do d = 1, size(grid)
          velo    =>  q(S_VELO,k)%data(d)%elts
          dvelo   => dq(S_VELO,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS)%data(d)%elts
          qe      => grid(d)%qe%elts
          divu    => grid(d)%divu%elts
          vort    => grid(d)%vort%elts
          
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (du_source, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 0)
          end do

          nullify (velo, dvelo, h_mflux, divu, qe, vort)
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate trend on coarser scales, from fine to coarse !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do l = level_end-1, level_start, -1
          call update_vector_bdry__finish (dq(S_MASS:S_TEMP,k), l+1) ! <= finish non-blocking communicate dmass (l+1)
          do d = 1, size(grid)
             mass      =>  q(S_MASS,k)%data(d)%elts
             velo      =>  q(S_VELO,k)%data(d)%elts
             temp      =>  q(S_TEMP,k)%data(d)%elts
             dmass     => dq(S_MASS,k)%data(d)%elts
             dtemp     => dq(S_TEMP,k)%data(d)%elts
             exner     => exner_fun(k)%data(d)%elts
             h_mflux   => horiz_flux(S_MASS)%data(d)%elts
             h_tflux   => horiz_flux(S_TEMP)%data(d)%elts
             bernoulli => grid(d)%bernoulli%elts
             divu      => grid(d)%divu%elts
             vort      => grid(d)%vort%elts
             qe        => grid(d)%qe%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do

             do j = 1, grid(d)%lev(l)%length
                call step1 (grid(d), grid(d)%lev(l)%elts(j), k)
             end do
             call apply_to_penta_d (post_step1, grid(d), l, k)

             call cpt_or_restr_Bernoulli_Exner (grid(d), l)

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_divu, grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
             end do
             
             call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) & use dmass (l+1)

             nullify (mass, velo, temp, dmass, dtemp, h_mflux, h_tflux, bernoulli, divu, exner, vort, qe)
          end do

          call update_vector_bdry (horiz_flux, l)

          ! Compute scalar trends at level l
          do d = 1, size(grid)
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS)%data(d)%elts
             h_tflux => horiz_flux(S_TEMP)%data(d)%elts
           
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do

             nullify (dmass, dtemp, h_mflux, h_tflux)
          end do

          dq(S_MASS:S_TEMP,k)%bdry_uptodate = .False.
          if (l .gt. level_start) call update_vector_bdry__start (dq(S_MASS:S_TEMP,k), l)  ! <= start non-blocking communicate dmass (l+1)

          ! Velocity trend, source part
          do d = 1, size(grid)
             velo    => q(S_VELO,k)%data(d)%elts
             dvelo   => dq(S_VELO,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS)%data(d)%elts
             qe      => grid(d)%qe%elts
             divu    => grid(d)%divu%elts
             vort    => grid(d)%vort%elts

             call cpt_or_restr_du_source (grid(d), l)
             nullify (velo, dvelo, h_mflux, divu, qe, vort)
          end do
          dq(S_VELO,k)%bdry_uptodate = .False.
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Evaluate complete velocity trend by adding gradient terms on entire grid at vertical level k !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do d = 1, size(grid)
          dvelo => dq(S_VELO,k)%data(d)%elts
          if (itype.eq.0) then ! Full trend
             mass      =>  q(S_MASS,k)%data(d)%elts
             temp      =>  q(S_TEMP,k)%data(d)%elts
             exner     => exner_fun(k)%data(d)%elts
             bernoulli => grid(d)%bernoulli%elts
          elseif (itype.eq.1) then ! Slow part only
             bernoulli => grid(d)%bern_slow%elts
          end if

          do p = 3, grid(d)%patch%length
             if (itype.eq.0) then ! Full trend
                call apply_onescale_to_patch (du_gradB_gradExn, grid(d), p-1, k, 0, 0)
             elseif (itype.eq.1) then ! Slow part only
                call apply_onescale_to_patch (du_gradB_slow, grid(d),    p-1, k, 0, 0)
             end if
          end do
          nullify (dvelo, bernoulli)
          if (itype.eq.0) nullify (mass, temp, exner)
       end do
    end do
  end subroutine trend_ml

  subroutine trend_fast (q, dq)
    ! Compute fast part of trend of prognostic variables assuming Lagrangian vertical coordinates
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q, dq
    
    integer                                                       :: d, k, p

    call update_array_bdry (q, NONE)

    ! Evaluate complete velocity trend by adding gradient terms on entire grid at vertical level k !
    do k = 1, zlevels
       do d = 1, size(grid)
          mass      =>  q(S_MASS,k)%data(d)%elts
          temp      =>  q(S_TEMP,k)%data(d)%elts
          dvelo     => dq(S_VELO,k)%data(d)%elts
          exner     => exner_fun(k)%data(d)%elts
          bernoulli => bernoulli_fast(k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (du_gradB_gradExn, grid(d), p-1, k, 0, 0)
          end do
          nullify (mass, temp, dvelo, exner, bernoulli)
       end do
    end do
  end subroutine trend_fast
 
  subroutine trend_diffuse (q, dq)
    ! Compute trends of prognostic variables
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q, dq

    integer                                                       :: d, j, k, l, p

    call update_array_bdry (q, NONE)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate trend on finest scale !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute each vertical level starting from surface
    do k = 1, zlevels
       ! Compute scalar fluxes and div, vort
       do d = 1, size(grid)
          mass    => q(S_MASS,k)%data(d)%elts
          temp    => q(S_TEMP,k)%data(d)%elts
          velo    => q(S_VELO,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS)%data(d)%elts
          h_tflux => horiz_flux(S_TEMP)%data(d)%elts
          divu    => grid(d)%divu%elts
          vort    => grid(d)%vort%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (flux_grad_scalar, grid(d), grid(d)%lev(level_end)%elts(j), z_null, -1, 1)
             call apply_onescale_to_patch (cal_divu,         grid(d), grid(d)%lev(level_end)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,         grid(d), grid(d)%lev(level_end)%elts(j), z_null, -1, 0)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_end, k)

          nullify (mass, temp, velo, h_mflux, h_tflux, divu, vort)
       end do

       if (level_start .lt. level_end) call update_vector_bdry__start (horiz_flux, level_end) ! <= comm flux (Jmax)

       ! Compute trends at finest level
       do d = 1, size(grid)
          dmass   => dq(S_MASS,k)%data(d)%elts
          dtemp   => dq(S_TEMP,k)%data(d)%elts
          dvelo   => dq(S_VELO,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS)%data(d)%elts
          h_tflux => horiz_flux(S_TEMP)%data(d)%elts
          divu  => grid(d)%divu%elts
          vort  => grid(d)%vort%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (scalar_trend,      grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 1)
             call apply_onescale_to_patch (du_source_diffuse, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 0)
          end do

          nullify (dmass, dtemp, dvelo, h_mflux, h_tflux, divu, vort)
       end do

       dq(:,k)%bdry_uptodate    = .False.
       horiz_flux%bdry_uptodate = .False.
       if (level_start .lt. level_end) then
          call update_vector_bdry__finish (horiz_flux, level_end) ! <= finish non-blocking communicate mass flux (Jmax)
          call update_vector_bdry__start (dq(S_MASS:S_TEMP,k), level_end) ! <= start non-blocking communicate dmass (l+1)
       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate trend on coarser scales, from fine to coarse !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do l = level_end-1, level_start, -1
          call update_vector_bdry__finish (dq(S_MASS:S_TEMP,k), l+1) ! <= finish non-blocking communicate dmass (l+1)

          ! Compute (or restrict) scalar fluxes and compute div and vort at scale l
          do d = 1, size(grid)
             mass    =>  q(S_MASS,k)%data(d)%elts
             temp    =>  q(S_TEMP,k)%data(d)%elts
             velo    =>  q(S_VELO,k)%data(d)%elts
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS)%data(d)%elts
             h_tflux => horiz_flux(S_TEMP)%data(d)%elts
             divu    => grid(d)%divu%elts
             vort    => grid(d)%vort%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (flux_grad_scalar, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
                call apply_onescale_to_patch (cal_divu,         grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
                call apply_onescale_to_patch (cal_vort,         grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
             end do
             call apply_to_penta_d (post_vort, grid(d), l, k)

             call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) & use dmass (l+1)

             nullify (mass, temp, velo, dmass, dtemp, h_mflux, h_tflux, divu, vort)
          end do

          call update_vector_bdry (horiz_flux, l)

          ! Compute trends at level l
          do d = 1, size(grid)
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             dvelo   => dq(S_VELO,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS)%data(d)%elts
             h_tflux => horiz_flux(S_TEMP)%data(d)%elts
             divu    => grid(d)%divu%elts
             vort    => grid(d)%vort%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             call cpt_or_restr_du_source_diffuse (grid(d), l)

             nullify (dmass, dtemp, dvelo, h_mflux, h_tflux, divu, vort)
          end do

          dq(S_MASS:S_TEMP,k)%bdry_uptodate = .False.
          if (l .gt. level_start) call update_vector_bdry__start (dq(S_MASS:S_TEMP,k), l)  ! <= start non-blocking communicate dmass (l+1)
          dq(S_VELO,k)%bdry_uptodate = .False.
       end do
    end do
  end subroutine trend_diffuse

  subroutine cpt_or_restr_Bernoulli_Exner (dom, l)
    type(Domain) :: dom
    integer      :: l

    integer :: j, p_par, c, p_chd
    logical :: restrict(N_CHDRN)

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .False.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .gt. 0) restrict(c) = .True.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) then
             call apply_interscale_to_patch3 (Bernoulli_Exner_cpt_restr, dom, p_par, c, z_null, 0, 1)
          end if
       end do
    end do
  end subroutine cpt_or_restr_Bernoulli_Exner

  subroutine Bernoulli_Exner_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute or restrict Bernoulli and Exner functions
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) .ge. RESTRCT) then
       bernoulli(id_par+1) = bernoulli(id_chd+1)
       dom%bern_slow%elts(id_par+1) = dom%bern_slow%elts(id_chd+1)
       dom%bern_fast%elts(id_par+1) = dom%bern_fast%elts(id_chd+1)
       exner(id_par+1) = exner(id_chd+1)
    end if
  end subroutine Bernoulli_Exner_cpt_restr

  subroutine cpt_or_restr_du_source (dom, l)
    type(Domain) :: dom
    integer      :: l

    integer :: c, j, p_par, p_chd

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) call apply_onescale_to_patch (du_source, dom, p_par, z_null, 0, 0)
       end do
       call apply_interscale_to_patch (du_source_cpt_restr, dom, dom%lev(l)%elts(j), z_null, 0, 0)
    end do
  end subroutine cpt_or_restr_du_source

  subroutine du_source_cpt_restr (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain)                     :: dom
    integer                          :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_par, dims_chd

    integer :: id_par, id_chd, idE_chd, idNE_chd, idN_chd

    id_par   = idx(i_par,     j_par,     offs_par, dims_par)

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)

    if (minval(dom%mask_e%elts(EDGE*id_chd + RT + 1:EDGE*id_chd + UP + 1)) .lt. ADJZONE) then
       call du_source (dom, i_par, j_par, zlev, offs_par, dims_par)
    end if

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) then
       dvelo(EDGE*id_par+RT+1) = dvelo(EDGE*id_chd+RT+1) + dvelo(EDGE*idE_chd+RT+1)
    end if

    if (dom%mask_e%elts(DG+EDGE*id_chd+1) .ge. ADJZONE) then
       dvelo(EDGE*id_par+DG+1) = dvelo(EDGE*idNE_chd+DG+1) + dvelo(EDGE*id_chd+DG+1)
    end if

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) then
       dvelo(EDGE*id_par+UP+1) = dvelo(EDGE*id_chd+UP+1) + dvelo(EDGE*idN_chd+UP+1)
    end if
  end subroutine du_source_cpt_restr

  subroutine cpt_or_restr_du_source_diffuse (dom, l)
    type(Domain) :: dom
    integer      :: l

    integer :: c, j, p_par, p_chd

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) then
             call apply_onescale_to_patch (du_source_diffuse, dom, p_par, z_null, 0, 0)
          end if
       end do
       call apply_interscale_to_patch (du_source_diffuse_cpt_restr, dom, dom%lev(l)%elts(j), z_null, 0, 0)
    end do
  end subroutine cpt_or_restr_du_source_diffuse

  subroutine du_source_diffuse_cpt_restr (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain)                     :: dom
    integer                          :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_par, dims_chd

    integer :: id_par, id_chd, idE_chd, idNE_chd, idN_chd
    real (8) :: u_prim_RT, u_prim_RT_E, u_prim_DG, u_prim_DG_NE, u_prim_UP, u_prim_UP_N

    id_par   = idx(i_par,     j_par,     offs_par, dims_par)

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)

    if (minval(dom%mask_e%elts(EDGE*id_chd + RT + 1:EDGE*id_chd + UP + 1)) .lt. ADJZONE) then
       call du_source_diffuse (dom, i_par, j_par, zlev, offs_par, dims_par)
    end if

    ! Restriction is on edge integrated velocity
    u_prim_RT    = dvelo(EDGE*id_chd+RT+1)*dom%len%elts(EDGE*id_chd+RT+1)
    u_prim_RT_E  = dvelo(EDGE*idE_chd+RT+1)*dom%len%elts(EDGE*idE_chd+RT+1)
    u_prim_DG    = dvelo(EDGE*id_chd+DG+1)*dom%len%elts(EDGE*id_chd+DG+1)
    u_prim_DG_NE = dvelo(EDGE*idNE_chd+DG+1)*dom%len%elts(EDGE*idNE_chd+DG+1)
    u_prim_UP    = dvelo(EDGE*id_chd+UP+1)*dom%len%elts(EDGE*id_chd+UP+1)
    u_prim_UP_N  = dvelo(EDGE*idN_chd+UP+1)*dom%len%elts(EDGE*idN_chd+UP+1)
    
    if (dom%mask_e%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) &
         dvelo(EDGE*id_par+RT+1) = (u_prim_RT + u_prim_RT_E)/dom%len%elts(EDGE*id_par+RT+1)

    if (dom%mask_e%elts(DG+EDGE*id_chd+1) .ge. ADJZONE) &
         dvelo(EDGE*id_par+DG+1) = (u_prim_DG + u_prim_DG_NE)/dom%len%elts(EDGE*id_par+DG+1)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) &
         dvelo(EDGE*id_par+UP+1) = (u_prim_UP + u_prim_UP_N)/dom%len%elts(EDGE*id_par+UP+1)
  end subroutine du_source_diffuse_cpt_restr

  subroutine cpt_or_restr_flux (dom, l)
    type(Domain) :: dom
    integer      :: l

    integer                     :: j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .False.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .gt. 0) restrict(c) = .True.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) then
             call apply_interscale_to_patch3 (flux_cpt_restr, dom, p_par, c, z_null, 0, 1)
          end if
       end do
    end do
  end subroutine cpt_or_restr_flux

  subroutine flux_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute flux restriction of mass and potential temperature by summing coarse, corrective and small fluxes
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: id_par, id_chd
    real(8), dimension(4) :: sm_flux_m, sm_flux_t

    id_chd   = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par   = idx(i_par, j_par, offs_par, dims_par)

    if (i_chd .ge. PATCH_SIZE .or. j_chd .ge. PATCH_SIZE) return

    call flux_restr (dmass, h_mflux)
    call flux_restr (dtemp, h_tflux)
  contains
    subroutine flux_restr (dscalar, h_flux)
      real(8), dimension(:), pointer :: dscalar, h_flux

      real(8), dimension(4) :: sm_flux

      if (maxval(dom%mask_e%elts(EDGE*id_par+RT+1:EDGE*id_par+UP+1)) .ge. RESTRCT) then
         sm_flux = interp_flux(h_flux, dom, i_chd, j_chd, offs_chd, dims_chd)
      end if

      if (dom%mask_e%elts(EDGE*id_par+RT+1) .ge. RESTRCT) then
         h_flux(EDGE*id_par+RT+1) = &
              complete_coarse_flux (dscalar, h_flux, sm_flux, dom, i_par, j_par, i_chd, j_chd, RT, offs_chd, dims_chd)
      end if

      if (dom%mask_e%elts(EDGE*id_par+DG+1) .ge. RESTRCT) then
         h_flux(EDGE*id_par+DG+1) = &
              complete_coarse_flux (dscalar, h_flux, sm_flux, dom, i_par, j_par, i_chd, j_chd, DG, offs_chd, dims_chd)
      end if

      if (dom%mask_e%elts(EDGE*id_par+UP+1) .ge. RESTRCT) then
         h_flux(EDGE*id_par+UP+1) = &
              complete_coarse_flux (dscalar, h_flux, sm_flux, dom, i_par, j_par, i_chd, j_chd, UP, offs_chd, dims_chd)
      end if
    end subroutine flux_restr

    function complete_coarse_flux (dscalar, flux, sm_flux, dom, i_par, j_par, i_chd, j_chd, e, offs_chd, dims_chd)
      real(8)                        :: complete_coarse_flux
      real(8), dimension(:), pointer :: dscalar, flux
      integer                        :: i_par, j_par, i_chd, j_chd
      integer, dimension(N_BDRY+1)    :: offs_chd
      integer, dimension(2,N_BDRY+1) :: dims_chd

      type(Domain)           :: dom
      integer                :: e
      real(8)                :: p_flux, c_flux
      real(8), dimension(4)  :: sm_flux

      if (e .eq. RT) then
         p_flux = part_coarse_flux (dscalar, flux, dom, i_chd+1, j_chd, RT, offs_chd, dims_chd)
         c_flux = coarse_flux(dscalar, dom, i_par, j_par, i_chd+1, j_chd, RT)
         complete_coarse_flux = p_flux + c_flux + sm_flux(1) + sm_flux(2)
      elseif (e .eq. DG) then
         p_flux = part_coarse_flux (dscalar, flux, dom, i_chd+1, j_chd+1, DG, offs_chd, dims_chd)
         c_flux = coarse_flux(dscalar, dom, i_par, j_par, i_chd+1, j_chd+1, DG)
         complete_coarse_flux = p_flux + c_flux + sm_flux(2) + sm_flux(3)
      elseif (e .eq. UP) then
         p_flux = part_coarse_flux (dscalar, flux, dom, i_chd, j_chd+1, UP, offs_chd, dims_chd)
         c_flux = coarse_flux(dscalar, dom, i_par, j_par, i_chd, j_chd+1, UP)
         complete_coarse_flux = p_flux + c_flux + sm_flux(3) + sm_flux(4)
      end if
    end function complete_coarse_flux

    function coarse_flux (dscalar, dom, i_par, j_par, i_chd, j_chd, e)
      real(8)                        :: coarse_flux
      real(8), dimension(:), pointer :: dscalar
      integer                        :: i_par, j_par, i_chd, j_chd, e

      type(Domain) :: dom
      integer      :: id_mz, id_pz, id_mp,  id_pp, id_pm, id_mm, id_mm2, id_pm2, id_pp2,  id_mp2, id

      id_mz = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 1 + 1), offs_chd, dims_chd)
      id_pz = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 4 + 1), offs_chd, dims_chd)
      id_mp = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 0 + 1), offs_chd, dims_chd)
      id_pp = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 5 + 1), offs_chd, dims_chd)
      id_pm = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 3 + 1), offs_chd, dims_chd)
      id_mm = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 2 + 1), offs_chd, dims_chd)

      id_mm2 = idx2(i_chd, j_chd, bfly_no2(:,1,e+1), offs_chd, dims_chd)
      id_pm2 = idx2(i_chd, j_chd, bfly_no2(:,2,e+1), offs_chd, dims_chd)
      id_pp2 = idx2(i_chd, j_chd, bfly_no2(:,3,e+1), offs_chd, dims_chd)
      id_mp2 = idx2(i_chd, j_chd, bfly_no2(:,4,e+1), offs_chd, dims_chd)

      id = idx(i_chd, j_chd, offs_chd, dims_chd)

      coarse_flux = (dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
           + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
           + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
           + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
           *(dscalar(id_pz+1) - dscalar(id_mz+1)) + &
           dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
           *0.5_8*(dscalar(id_pp2+1) - dscalar(id_mz+1)) + &
           dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           *0.5_8*(dscalar(id_pm2+1) - dscalar(id_mz+1)) + &
           dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
           *0.5_8*(dscalar(id_pz+1) - dscalar(id_mp2+1)) + &
           dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
           *0.5_8*(dscalar(id_pz+1) - dscalar(id_mm2+1))
    end function coarse_flux

    subroutine get_indices (dom, i, j, e, offs, dims, id)
      type(Domain)                   :: dom
      integer                        :: i, j, e
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer, dimension(2)  :: ij_mp, ij_pp, ij_pm, ij_mm
      integer, dimension(20) :: id

      id(UMZ+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims)
      id(UPZ+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims)
      id(WMP+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims)
      id(VPP+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims)
      id(WPM+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims)
      id(VMM+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims)

      ij_mp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 0 + 1)
      id(MP+1) = idx(ij_mp(1), ij_mp(2), offs, dims)
      ij_pp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 5 + 1)
      id(PP+1) = idx(ij_pp(1), ij_pp(2), offs, dims)
      ij_pm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 3 + 1)
      id(PM+1) = idx(ij_pm(1), ij_pm(2), offs, dims)
      ij_mm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 2 + 1)

      id(MM+1)   = idx(ij_mm(1), ij_mm(2), offs, dims)

      id(VMP+1)  = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
      id(VMPP+1) = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1),   offs, dims)
      id(UZP+1)  = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,hex_s_offs(e+1) + 0 + 4 + 1),   offs, dims)
      id(WPPP+1) = ed_idx(ij_pp(1), ij_pp(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
      id(WPP+1)  = ed_idx(ij_pp(1), ij_pp(2), hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1),   offs, dims)
      id(VPM+1)  = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1),   offs, dims)
      id(VPMM+1) = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
      id(UZM+1)  = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,(hex_s_offs(e+1) + 3) - 2 + 1), offs, dims)
      id(WMMM+1) = ed_idx(ij_mm(1), ij_mm(2), hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1),   offs, dims)
      id(WMM+1)  = ed_idx(ij_mm(1), ij_mm(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
    end subroutine get_indices

    function interp_flux (flux, dom, i, j, offs, dims)
      real(8), dimension(4)          :: interp_flux
      real(8), dimension(:), pointer :: flux
      integer                        ::  i, j
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      type(Domain)           :: dom
      real(4)                :: wgt
      integer, dimension(20) :: id

      call get_indices (dom, i+1, j, RT, offs, dims, id)

      interp_flux(1) = - sum(flux(id((/WPM,UZM,VMM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1, j-2, offs, dims)+1)%enc) &
           - sum((flux(id((/VPM,WMMM,UMZ/)+1)+1) -flux(id((/UPZ,VPMM,WMM/)+1)+1)) * &
           dom%R_F_wgt%elts(idx(i+1, j-1, offs, dims)+1)%enc) ! UPLT S

      interp_flux(2) = sum(flux(id((/WMP,UZP,VPP/)+1)+1)* dom%R_F_wgt%elts(idx(i, j, offs, dims)+1)%enc) &
           + sum((flux(id((/VMP,WPPP,UPZ/)+1)+1) - flux(id((/UMZ,VMPP,WPP/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i, j+1, offs, dims)+1)%enc) ! LORT

      call get_indices(dom, i, j+1, UP, offs, dims, id)

      interp_flux(3) = - sum(flux(id((/UZM,VMM,WPM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j, offs, dims)+1)%enc) &
           - sum((flux(id((/WMMM,UMZ,VPM/)+1)+1) - flux(id((/VPMM,WMM,UPZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i+1, j+1, offs, dims)+1)%enc) ! UPLT

      interp_flux(4) = sum(flux(id((/UZP,VPP,WMP/)+1)+1) * dom%R_F_wgt%elts(idx(i-2,j, offs, dims)+1)%enc) &
           + sum((flux(id((/WPPP,UPZ,VMP/)+1)+1) - flux(id((/VMPP,WPP,UMZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i-2, j+1, offs, dims)+1)%enc) ! LORT W
    end function interp_flux

    function part_coarse_flux (dscalar, flux, dom, i, j, e, offs, dims)
      real(8)                        :: part_coarse_flux
      real(8), dimension(:), pointer :: dscalar, flux
      integer                        :: i, j, e
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      type(Domain)           :: dom      
      real(8), dimension(2)  :: area
      real(8), dimension(4)  :: ol_area
      integer, dimension(20) :: id

      call get_indices (dom, i, j, e, offs, dims, id)

      area = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%a(1:2)
      ol_area(1:2) = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%split
      ol_area(3:4) = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%a(3:4) - ol_area(1:2)
      area(1) = area(1) + ol_area(1) + ol_area(4)
      area(2) = area(2) + ol_area(2) + ol_area(3)
      area = area/sum(area)

      ol_area(1) = dom%overl_areas%elts(id(PP+1)+1)%split(1)
      ol_area(2) = dom%overl_areas%elts(id(MM+1)+1)%split(2)
      ol_area(3) = dom%overl_areas%elts(id(MP+1)+1)%a(3) - dom%overl_areas%elts(id(MP+1)+1)%split(1)
      ol_area(4) = dom%overl_areas%elts(id(PM+1)+1)%a(4) - dom%overl_areas%elts(id(PM+1)+1)%split(2)

      part_coarse_flux = &
           sum(flux(id((/UPZ,UMZ/)+1)+1)*area) - sum(flux(id((/VMM,WMP/)+1)+1))*area(2) - sum(flux(id((/WPM,VPP/)+1)+1))*area(1) &
           + ol_area(3)*dscalar(id(MP+1)+1) - ol_area(4)*dscalar(id(PM+1)+1)  &
           - ol_area(1)*dscalar(id(PP+1)+1) + ol_area(2)*dscalar(id(MM+1)+1)
    end function part_coarse_flux
  end subroutine flux_cpt_restr

  subroutine post_refine
    integer :: d, p

    level_end = sync_max(level_end)

    do d = 1, n_domain(rank+1)
       do p = 3, grid(d)%patch%length
          call connect_children(grid(d), p-1)
       end do
    end do

    call comm_patch_conn_mpi

    do d = 1, size(grid)
       call update_comm(grid(d))
    end do

    call comm_communication_mpi
    call comm_nodes9_mpi(get_areas, set_areas, NONE)
    call apply_to_penta(area_post_comm, NONE, z_null)
  end subroutine post_refine

   ! subroutine trend_mass_coord (q, dq)
  !   ! Compute trends of prognostic variables using mass-based coordinate
  !   ! NOT TESTED !!
  !   type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q, dq

  !   type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: vert_flux
  !   integer                                                       :: d, j, k, l, p

  !   call update_array_bdry (q, NONE)

  !   vert_flux = q(S_MASS,:) ! Set correct size of vertical flux (will be over-written)

  !   ! First integrate pressure down across all grid points in order to compute surface pressure
  !   do k = zlevels, 1, -1
  !      do d = 1, size(grid)
  !         mass => q(S_MASS,k)%data(d)%elts
  !         temp => q(S_TEMP,k)%data(d)%elts

  !         do p = 3, grid(d)%patch%length
  !            call apply_onescale_to_patch (integrate_pressure_down, grid(d), p-1, k, 0, 1)
  !         end do

  !         nullify (mass, temp)
  !      end do
  !   end do

  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ! Calculate trend on finest scale !
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !   ! Compute each vertical level starting from surface
  !   do k = 1, zlevels
  !      do d = 1, size(grid)
  !         mass      => q(S_MASS,k)%data(d)%elts
  !         velo      => q(S_VELO,k)%data(d)%elts
  !         temp      => q(S_TEMP,k)%data(d)%elts
  !         h_mflux   => horiz_flux(S_MASS)%data(d)%elts
  !         h_tflux   => horiz_flux(S_TEMP)%data(d)%elts
  !         bernoulli => grid(d)%bernoulli%elts
  !         exner     => grid(d)%exner%elts
  !         vort      => grid(d)%vort%elts
  !         qe        => grid(d)%qe%elts

  !         ! Compute pressure, geopotential, Exner (compressible case), specific volume
  !         do j = 1, grid(d)%lev(level_end)%length
  !            call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(level_end)%elts(j), k, 0, 1)
  !         end do

  !         ! Compute horizontal fluxes, potential vorticity (qe), Bernoulli, Exner (incompressible case)
  !         do j = 1, grid(d)%lev(level_end)%length
  !            call step1 (grid(d), grid(d)%lev(level_end)%elts(j), k)
  !         end do
  !         call apply_to_penta_d (post_step1, grid(d), level_end, k)

  !         nullify (mass, velo, temp, h_mflux, h_tflux, bernoulli, exner, vort, qe)
  !      end do

  !      if (level_start .lt. level_end) call update_vector_bdry__start (horiz_flux, level_end) ! <= comm flux (Jmax)

  !      ! Compute vertically integrated horizontal mass flux (stored in integr_horiz_flux)
  !      if (.not. lagrangian_vertical) then
  !         do d = 1, size(grid)
  !            do j = 1, grid(d)%lev(level_end)%length
  !               call apply_onescale_to_patch (vert_integrated_horiz_flux, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 1)
  !            end do
  !         end do
  !      end if

  !      ! Compute scalar trends at finest level
  !      do d = 1, size(grid)
  !         dmass   => dq(S_MASS,k)%data(d)%elts
  !         dtemp   => dq(S_TEMP,k)%data(d)%elts
  !         mass => q(S_MASS,k)%data(d)%elts
  !         temp => q(S_TEMP,k)%data(d)%elts

  !         if (k<zlevels) then ! Provide temperature and mass at next vertical level
  !            adj_mass_up => q(S_MASS,k+1)%data(d)%elts
  !            adj_temp_up => q(S_TEMP,k+1)%data(d)%elts
  !         end if

  !         ! Use vertically integrated horizontal fluxes
  !         h_mflux => grid(d)%integr_horiz_flux%elts

  !         ! Save vertical mass flux at all vertical levels in vert_flux
  !         v_mflux => vert_flux(k)%data(d)%elts 
  !         h_tflux => horiz_flux(S_TEMP)%data(d)%elts

  !         do j = 1, grid(d)%lev(level_end)%length
  !            call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 1)
  !         end do

  !         nullify (mass, temp, dmass, dtemp, h_mflux, h_tflux, v_mflux, adj_mass_up, adj_temp_up)
  !      end do

  !      dq(:,k)%bdry_uptodate    = .False.
  !      horiz_flux%bdry_uptodate = .False.
  !      if (level_start .lt. level_end) then
  !         call update_vector_bdry__finish (horiz_flux, level_end) ! <= finish non-blocking communicate mass flux (Jmax)
  !         call update_vector_bdry__start (dq(S_MASS:S_TEMP,k), level_end) ! <= start non-blocking communicate dmass (l+1)
  !      end if

  !      ! Velocity trend, source part
  !      do d = 1, size(grid)
  !         velo    =>  q(S_VELO,k)%data(d)%elts
  !         dvelo   => dq(S_VELO,k)%data(d)%elts
  !         h_mflux => horiz_flux(S_MASS)%data(d)%elts
  !         qe      => grid(d)%qe%elts

  !         do j = 1, grid(d)%lev(level_end)%length
  !            call apply_onescale_to_patch (du_source, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 0)
  !         end do

  !         nullify (velo, dvelo, h_mflux, qe)
  !      end do

  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      ! Calculate trend on coarser scales, from fine to coarse !
  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      do l = level_end-1, level_start, -1
  !         call update_vector_bdry__finish (dq(S_MASS:S_TEMP,k), l+1) ! <= finish non-blocking communicate dmass (l+1)
  !         do d = 1, size(grid)
  !            mass      =>  q(S_MASS,k)%data(d)%elts
  !            velo      =>  q(S_VELO,k)%data(d)%elts
  !            temp      =>  q(S_TEMP,k)%data(d)%elts
  !            dmass     =>  dq(S_MASS,k)%data(d)%elts
  !            dtemp     =>  dq(S_TEMP,k)%data(d)%elts
  !            h_mflux   => horiz_flux(S_MASS)%data(d)%elts
  !            h_tflux   => horiz_flux(S_TEMP)%data(d)%elts
  !            bernoulli => grid(d)%bernoulli%elts
  !            exner     => grid(d)%exner%elts
  !            vort      => grid(d)%vort%elts
  !            qe        => grid(d)%qe%elts

  !            do j = 1, grid(d)%lev(l)%length
  !               call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
  !            end do

  !            do j = 1, grid(d)%lev(l)%length
  !               call step1 (grid(d), grid(d)%lev(l)%elts(j), k)
  !            end do
  !            call apply_to_penta_d (post_step1, grid(d), l, k)

  !            call cpt_or_restr_Bernoulli_Exner (grid(d), l)

  !            call cpt_or_restr_flux (grid(d), l)  ! <= compute flux(l) & use dmass (l+1)

  !            nullify (mass, velo, temp, dmass, dtemp, h_mflux, h_tflux, bernoulli, exner, vort, qe)
  !         end do

  !         call update_vector_bdry (horiz_flux, l)

  !         ! Compute vertically integrated horizontal mass flux
  !         do d = 1, size(grid)
  !            do j = 1, grid(d)%lev(l)%length
  !               call apply_onescale_to_patch (vert_integrated_horiz_flux, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
  !            end do
  !         end do

  !         ! Compute scalar trends at level l
  !         do d = 1, size(grid)
  !            dmass   => dq(S_MASS,k)%data(d)%elts
  !            dtemp   => dq(S_TEMP,k)%data(d)%elts
  !            h_mflux => horiz_flux(S_MASS)%data(d)%elts
  !            h_tflux => horiz_flux(S_TEMP)%data(d)%elts
  !            mass => q(S_MASS,k)%data(d)%elts
  !            temp => q(S_TEMP,k)%data(d)%elts
  !            h_mflux => grid(d)%integr_horiz_flux%elts ! Use vertically integrated horizontal fluxes
  !            v_mflux => vert_flux(k)%data(d)%elts ! We save vertical mass flux at all levels at scale l in vert_flux

  !            if (k<zlevels) then ! Provide temperature and mass at next vertical level
  !               adj_mass_up => q(S_MASS,k+1)%data(d)%elts
  !               adj_temp_up => q(S_TEMP,k+1)%data(d)%elts
  !            end if

  !            do j = 1, grid(d)%lev(l)%length
  !               call apply_onescale_to_patch (scalar_trend, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
  !            end do

  !            nullify (dmass, dtemp, h_mflux, h_tflux, mass, temp, v_mflux, adj_mass_up, adj_temp_up)
  !         end do

  !         dq(S_MASS:S_TEMP,k)%bdry_uptodate = .False.
  !         if (l .gt. level_start) call update_vector_bdry__start (dq(S_MASS:S_TEMP,k), l)  ! <= start non-blocking communicate dmass (l+1)

  !         ! Velocity trend, source part
  !         do d = 1, size(grid)
  !            velo    => q(S_VELO,k)%data(d)%elts
  !            dvelo   => dq(S_VELO,k)%data(d)%elts
  !            h_mflux => horiz_flux(S_MASS)%data(d)%elts
  !            qe      => grid(d)%qe%elts
  !            divu    => grid(d)%bernoulli%elts
  !            vort    => grid(d)%vort%elts

  !            call cpt_or_restr_du_source (grid(d), l)
  !            nullify (velo, dvelo, h_mflux, qe, divu, vort)
  !         end do
  !         dq(S_VELO,k)%bdry_uptodate = .False.
  !      end do

  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      ! Evaluate complete velocity trend by adding gradient terms on entire grid at vertical level k !
  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      do d = 1, size(grid)
  !         mass      =>  q(S_MASS,k)%data(d)%elts
  !         temp      =>  q(S_TEMP,k)%data(d)%elts
  !         dvelo     => dq(S_VELO,k)%data(d)%elts
  !         bernoulli => grid(d)%bernoulli%elts
  !         exner     => grid(d)%exner%elts
  !         do p = 3, grid(d)%patch%length
  !            call apply_onescale_to_patch (du_gradB_gradExn, grid(d), p-1, k, 0, 0)
  !         end do
  !         nullify (mass, temp, dvelo, bernoulli, exner)
  !      end do

  !   end do
  ! end subroutine trend_mass_coord
end module multi_level_mod
