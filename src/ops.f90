Module ops_mod
  use domain_mod
  use arch_mod
  use viscous_mod
  implicit none

  real(8) :: totaldmass, totalabsdmass, totaldtemp, totalabsdtemp
  real(8) :: sum_mass, sum_temp
  integer :: tic

contains
  subroutine init_ops_mod()
    logical :: initialized = .False.
    if (initialized) return ! initialize only once
    call init_domain_mod()
    initialized = .True.
  end subroutine init_ops_mod

  subroutine post_step1 (dom, p, c, offs, dims, zlev)
    type(Domain)                 :: dom
    integer                      :: p, c, zlev
    integer, dimension(N_BDRY+1) :: offs
    integer, dimension(2,9)      :: dims

    integer                      :: id, idS, idW, idSW, idN, idE, idNE
    real(8)                      :: pv_SW, pv_W, pv_S, pv_LORT, pv_UPLT, pv_SW_LORT, pv_SW_UPLT, pv
    real(8), dimension(0:N_BDRY) :: full_mass


    if (c .eq. IJMINUS) then
       id   = idx( 0,  0, offs, dims)
       idSW = idx(-1, -1, offs, dims)
       idW  = idx(-1,  0, offs, dims)
       idS  = idx( 0, -1, offs, dims)
       idN  = idx( 0,  1, offs, dims)
       idE  = idx( 1,  0, offs, dims)

       full_mass(0:WEST) = mass((/id,idN,idE,idS,idW/)+1) + mean(S_MASS,zlev)

       dom%vort%elts(TRIAG*idSW+LORT+1) = &
            (velo(EDGE*idW+RT+1)*dom%len%elts(EDGE*idW+RT+1) &
            -velo(EDGE*idSW+1)*dom%len%elts(EDGE*idSW+1) &
            -velo(EDGE*idS+UP+1)*dom%len%elts(EDGE*idS+UP+1))

       pv_SW = (dom%coriolis%elts(TRIAG*idSW+1) + dom%vort%elts(TRIAG*idSW+1))/ &
            (full_mass(WEST)*dom%areas%elts(idW+1)%part(6) &
            + full_mass(0)*sum(dom%areas%elts(id+1)%part(4:5)) &
            + full_mass(SOUTH)*dom%areas%elts(idS+1)%part(3))

       pv_W = (dom%coriolis%elts(TRIAG*idW+LORT+1) + dom%vort%elts(TRIAG*idW+LORT+1)*dom%triarea%elts(TRIAG*idW+LORT+1))/ &
            (full_mass(WEST)*dom%areas%elts(idW+1)%part(1) &
            + full_mass(0)*dom%areas%elts(id+1)%part(3) &
            + full_mass(NORTH)*dom%areas%elts(idN+1)%part(5))

       pv_S = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + dom%vort%elts(TRIAG*idS+UPLT+1)*dom%triarea%elts(TRIAG*idS+UPLT+1))/ &
            (full_mass(SOUTH)*dom%areas%elts(idS+1)%part(2) &
            + full_mass(EAST)*dom%areas%elts(idE+1)%part(4) &
            + full_mass(0)*dom%areas%elts(id+1)%part(6))

       dom%vort%elts(TRIAG*idSW+LORT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)/dom%triarea%elts(LORT+TRIAG*idSW+1)
       dom%vort%elts(TRIAG*idSW+UPLT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)

       qe(EDGE*idW+RT+1) = 0.5_8*(pv_W + pv_SW)
       qe(EDGE*idS+UP+1) = 0.5_8*(pv_S + pv_SW)
    end if

    if (c .eq. IPLUSJMINUS) then
       id   = idx(PATCH_SIZE,    0, offs, dims)
       idSW = idx(PATCH_SIZE-1, -1, offs, dims)
       idS  = idx(PATCH_SIZE,   -1, offs, dims)
       idW  = idx(PATCH_SIZE-1,  0, offs, dims)
       idE  = idx(PATCH_SIZE+1,  0, offs, dims)
       idNE = idx(PATCH_SIZE+1,  1, offs, dims)

       full_mass(0:NORTHEAST) = mass((/id,id,idE,idS,idW,idNE/)+1) + mean(S_MASS,zlev)

       dom%vort%elts(LORT+TRIAG*idSW+1) = - &
            ((velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1) + &
            velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)) - &
            velo(EDGE*id+RT+1)*dom%len%elts(EDGE*id+RT+1))

       pv_SW_LORT = (dom%coriolis%elts(LORT+TRIAG*idSW+1) + &
            dom%vort%elts(LORT+TRIAG*idSW+1))/( &
            full_mass(SOUTHWEST)*dom%areas%elts(idSW+1)%part(1) + &
            full_mass(SOUTH)*dom%areas%elts(idS +1)%part(3) + &
            full_mass(0)*sum(dom%areas%elts(id+1)%part(5:6)))

       pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + dom%vort%elts(TRIAG*id+LORT+1)*dom%triarea%elts(TRIAG*id+LORT+1))/ &
            (full_mass(0)*dom%areas%elts(id  +1)%part(1) &
            + full_mass(EAST)*dom%areas%elts(idE +1)%part(3) &
            + full_mass(NORTHEAST)*dom%areas%elts(idNE+1)%part(5))

       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + &
            dom%vort%elts(TRIAG*idSW+UPLT+1)*dom%triarea%elts(TRIAG*idSW+UPLT+1))/ &
            (full_mass(SOUTHWEST)*dom%areas%elts(idSW+1)%part(2) &
            + full_mass(0)*dom%areas%elts(id  +1)%part(4) &
            + full_mass(WEST)*dom%areas%elts(idW +1)%part(6))

       dom%vort%elts(TRIAG*idSW+LORT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)/dom%triarea%elts(LORT+TRIAG*idSW+1)
       dom%vort%elts(TRIAG*idS +UPLT+1) = dom%vort%elts(LORT+TRIAG*idSW+1)

       pv_S = pv_SW_LORT

      qe(EDGE*id  +RT+1) = 0.5_8*(pv_S + pv_LORT)
      qe(EDGE*idSW+DG+1) = 0.5_8*(pv_SW_LORT + pv_SW_UPLT)
    end if

    if (c .eq. IMINUSJPLUS) then
       id   = idx(0,  PATCH_SIZE,   offs, dims)
       idSW = idx(-1, PATCH_SIZE-1, offs, dims)
       idW  = idx(-1, PATCH_SIZE,   offs, dims)
       idS  = idx(0,  PATCH_SIZE-1, offs, dims)
       idN  = idx(0,  PATCH_SIZE+1, offs, dims)
       idNE = idx(1,  PATCH_SIZE+1, offs, dims)

       full_mass(0:NORTHEAST) = mass((/id,idN,id,idS,idW,idNE/)+1) + mean(S_MASS,zlev)

       dom%vort%elts(TRIAG*idSW+UPLT+1) = &
            - velo(EDGE*id+UP+1)*dom%len%elts(EDGE*id+UP+1) &
            + velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1) &
            + velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)

       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) &
            + dom%vort%elts(TRIAG*idSW+UPLT+1)) &
            /(full_mass(SOUTHWEST)*dom%areas%elts(idSW  +1)%part(2) &
            + full_mass(0)*sum(dom%areas%elts(id+1)%part(3:4)) &
            + full_mass(WEST)*dom%areas%elts(idW+1)%part(6))

       pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + dom%vort%elts(TRIAG*id+UPLT+1)*dom%triarea%elts(TRIAG*id+UPLT+1))/ &
            (full_mass(0)*dom%areas%elts(id+1)%part(2) &
            + full_mass(NORTHEAST)*dom%areas%elts(idNE+1)%part(4) &
            + full_mass(NORTH)*dom%areas%elts(idN+1)%part(6))

       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + &
            dom%vort%elts(TRIAG*idSW+LORT+1)*dom%triarea%elts(TRIAG*idSW+LORT+1))/ &
            (full_mass(SOUTHWEST)*dom%areas%elts(idSW+1)%part(1) &
            + full_mass(SOUTH)*dom%areas%elts(idS+1)%part(3) &
            + full_mass(0)*dom%areas%elts(id+1)%part(5))

       dom%vort%elts(TRIAG*idSW+UPLT+1) = dom%vort%elts(TRIAG*idSW+UPLT+1)/dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       dom%vort%elts(TRIAG*idW +LORT+1) = dom%vort%elts(TRIAG*idSW+UPLT+1)

       pv_W = pv_SW_UPLT

       qe(EDGE*id  +UP+1) = 0.5_8*(pv_W + pv_UPLT)
       qe(EDGE*idSW+DG+1) = 0.5_8*(pv_SW_LORT + pv_SW_UPLT)
    end if

    if (c .eq. IJPLUS) then
       id  = idx(PATCH_SIZE,   PATCH_SIZE,   offs, dims)
       idN = idx(PATCH_SIZE,   PATCH_SIZE+1, offs, dims)
       idE = idx(PATCH_SIZE+1, PATCH_SIZE,   offs, dims)
       idS = idx(PATCH_SIZE,   PATCH_SIZE-1, offs, dims)
       idW = idx(PATCH_SIZE-1, PATCH_SIZE,   offs, dims)

       full_mass(0:WEST) = mass((/id,idN,idE,idS,idW/)+1) + mean(S_MASS,zlev)

       dom%vort%elts(LORT+TRIAG*id+1) = - &
            (velo(EDGE*id +RT+1)*dom%len%elts(EDGE*id+RT+1) - &
            velo(EDGE*idN+RT+1)*dom%len%elts(EDGE*id+DG+1) - &
            velo(EDGE*id +UP+1)*dom%len%elts(EDGE*id+UP+1))

       pv = (dom%coriolis%elts(TRIAG*id+1) + dom%vort%elts(LORT+TRIAG*id+1))/ &          
            (full_mass(EAST)*dom%areas%elts(idE+1)%part(3) + &
            full_mass(0)*sum(dom%areas%elts(id+1)%part(1:2)) + &
            full_mass(NORTH)*dom%areas%elts(idN+1)%part(6))

       pv_W = (dom%coriolis%elts(TRIAG*idW+LORT+1) + dom%vort%elts(TRIAG*idW+LORT+1)*dom%triarea%elts(TRIAG*idW+LORT+1))/ &
            (full_mass(WEST)*dom%areas%elts(idW+1)%part(1) &
            + full_mass(0)*dom%areas%elts(id+1)%part(3) &
            + full_mass(NORTH)*dom%areas%elts(idN+1)%part(5))

       pv_S = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + dom%vort%elts(TRIAG*idS+UPLT+1)*dom%triarea%elts(TRIAG*idS+UPLT+1))/ &
            (full_mass(SOUTH)*dom%areas%elts(idS+1)%part(2) &
            + full_mass(EAST)*dom%areas%elts(idE+1)%part(4) &
            + full_mass(0)*dom%areas%elts(id+1)%part(6))

       dom%vort%elts(LORT+TRIAG*id+1) = dom%vort%elts(LORT+TRIAG*id+1)/dom%triarea%elts(LORT+TRIAG*id+1)
       dom%vort%elts(TRIAG*id+UPLT+1) = dom%vort%elts(LORT+TRIAG*id+1)

       qe(EDGE*id+RT+1) = 0.5_8*(pv + pv_S)
       qe(EDGE*id+UP+1) = 0.5_8*(pv + pv_W)
    end if
  end subroutine post_step1

  subroutine step1 (dom, p, zlev)
    type(Domain) dom
    integer p, zlev
    integer i, j, id, offs(0:N_BDRY), dims(2,N_BDRY)
    integer n, e, s, w, ne, sw
    real(8) u_prim_up, u_dual_up, u_prim_dg, u_dual_dg, u_prim_rt, u_dual_rt
    real(8) u_prim_dn, u_dual_dn, u_prim_sw, u_dual_sw, u_prim_lt, u_dual_lt
    real(8) pv_LORT, pv_UPLT, pv_S, pv_W, vort_W, vort_S, vort_LORT, vort_UPLT
    logical S_bdry, W_bdry
    real(8) full_mass(0:N_BDRY), full_temp(0:N_BDRY)

    call comp_offs3(dom, p, offs, dims)

    S_bdry = (dom%patch%elts(p+1)%neigh(SOUTH) .lt. 0)
    if (S_bdry) S_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(SOUTH)+1)%side .gt. 0)

    W_bdry = (dom%patch%elts(p+1)%neigh(WEST) .lt. 0)
    if (W_bdry) W_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(WEST)+1)%side .gt. 0)

    id = offs(0)

    n = PATCH_SIZE       ; ne = n+1
    w = offs(WEST)       ; e =  +1
    sw = offs(SOUTHWEST) ; s = offs(SOUTH)

    call comput()

    if (W_bdry .or. S_bdry) call comp_ijmin()

    w = -1; sw = s+w
    do id = offs(0)+1, offs(0)+LAST-1
       call comput();
       if (S_bdry) call comp_ijmin()
    end do

    e = offs(EAST); ne = dims(1,EAST) + e
    id = offs(0)+LAST
    call comput()
    if (S_bdry) call comp_ijmin()

    s = -PATCH_SIZE
    do j = 2, PATCH_SIZE-1
       id = offs(0) + PATCH_SIZE*(j-1)
       e = +1; ne = n+e
       w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*(j-1) ! Correct for dimension smaller than patch if boundary
       sw = w-dims(1,WEST)
       call comput()
       if (W_bdry) call comp_ijmin()

       w = -1; sw = s+w; ne = n+e
       do id = offs(0)+PATCH_SIZE*(j-1)+1, offs(0)+PATCH_SIZE*(j-1)+LAST-1
          call comput()
       end do

       id = offs(0)+PATCH_SIZE*j-1
       e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*(j-1); ne = e+dims(1,EAST)
       call comput()
    end do

    id = offs(0)+PATCH_SIZE*LAST
    n = offs(NORTH)
    w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST
    sw = w-dims(1,WEST)
    e = +1; ne = n+e
    call comput()
    if (W_bdry) call comp_ijmin()

    w = -1; sw = s+w
    do id = offs(0)+PATCH_SIZE*LAST+1, offs(0)+PATCH_SIZE*PATCH_SIZE-2
       call comput()
    end do

    id = offs(0)+PATCH_SIZE*PATCH_SIZE-1
    ne = offs(NORTHEAST)
    e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
    call comput()

    if (dom%patch%elts(p+1)%neigh(NORTH) .lt. 0) then ! Neighbour is boundary
       if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(NORTH)+1)%side .gt. 0) then ! Domain boundary
          id = offs(0)+PATCH_SIZE*LAST + offs(NORTH) ! id + n
          s =                                   - offs(NORTH) ! relative to current id
          w = offs(NORTHWEST)                   - offs(NORTH) 
          sw = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST - offs(NORTH)
          e = 1
          n = dims(1,NORTH)
          ne = n+e
          call comput()
          if (W_bdry) call comp_ijmin()

          w = -1
          sw = -offs(NORTH) + w
          do id = offs(0)+PATCH_SIZE*LAST+1+offs(NORTH), offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-2
             call comput()
          end do

          e = offs(NORTHEAST)              - offs(NORTH) 
          ne = e+dims(1,NORTHEAST)
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-1
          call comput()

          ! NORTHEAST point
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTHEAST)-1
          s = -offs(NORTHEAST) + offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
          sw= -offs(NORTHEAST)
          w = -offs(NORTHEAST) + offs(NORTH)
          n = dims(1,NORTHEAST)
          e = 1
          ne = n+e

          ! Avoid doing nodes twice
          if (offs(NORTHEAST).ne.1) call comput()
       end if
    end if

    if (dom%patch%elts(p+1)%neigh(EAST) .lt. 0) then ! neighbour is boundary
       if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(EAST)+1)%side .le. 0) return ! not domain bdry
       id = offs(0) + LAST + offs(EAST) ! id + n
       w  =                - offs(EAST) ! relative to current id
       s  = offs(SOUTHEAST)- offs(EAST) 
       sw = offs(SOUTH)    - offs(EAST)
       n  = dims(1,EAST)
       e  = 1
       ne = n+e
       call comput()
       if (S_bdry) call comp_ijmin()

       s = -dims(1,EAST)
       do j = 1, LAST-1
          id = offs(0) + LAST + offs(EAST) + j*dims(1,EAST)
          w  =                - offs(EAST) + j*(PATCH_SIZE-dims(1,EAST))
          sw = w-PATCH_SIZE
          call comput()
       end do

       id = offs(0) + LAST + offs(EAST) + LAST*dims(1,EAST)
       n  = offs(NORTHEAST) + LAST*(PATCH_SIZE-dims(1,EAST)) - offs(EAST)
       w  =                - offs(EAST) + LAST*(PATCH_SIZE-dims(1,EAST))
       sw = w-PATCH_SIZE
       ne = n+e
       call comput()
    end if
  contains
    subroutine comp_ijmin()
      real(8) vort_SW

      full_mass(SOUTHWEST) = mass(id+sw+1) + mean(S_MASS,zlev)
      full_temp(SOUTHWEST) = temp(id+sw+1) + mean(S_TEMP,zlev)

      vort_SW = - (velo(EDGE*(id+sw)+RT+1)*dom%len%elts(EDGE*(id+sw)+RT+1) + u_prim_sw + u_prim_dn)

      pv_LORT = (dom%coriolis%elts(TRIAG*(id+sw)+LORT+1) + vort_SW)/( &
           full_mass(SOUTHWEST)*dom%areas%elts(id+sw+1)%part(1) + &
           full_mass(SOUTH)*dom%areas%elts(id+s +1)%part(3) + &
           full_mass(0)*dom%areas%elts(id   +1)%part(5))

      vort_SW = u_prim_lt + u_prim_sw + velo(EDGE*(id+sw)+UP+1)*dom%len%elts(EDGE*(id+sw)+UP+1) 

      pv_UPLT = (dom%coriolis%elts(TRIAG*(id+sw)+UPLT+1) + vort_SW)/( &
           full_mass(SOUTHWEST)*dom%areas%elts(id+sw+1)%part(2) + &
           full_mass(0)*dom%areas%elts(id   +1)%part(4) + &
           full_mass(WEST)*dom%areas%elts(id+w +1)%part(6))

      dom%vort%elts(TRIAG*(id+w)+LORT+1) = vort_W/dom%triarea%elts(TRIAG*(id+w)+LORT+1) 
      dom%vort%elts(TRIAG*(id+s)+UPLT+1) = vort_S/dom%triarea%elts(TRIAG*(id+s)+UPLT+1) 
      
      qe(EDGE*(id+ w)+RT+1) = 0.5_8*(pv_W   +pv_UPLT)
      qe(EDGE*(id+sw)+DG+1) = 0.5_8*(pv_LORT+pv_UPLT)
      qe(EDGE*(id+s )+UP+1) = 0.5_8*(pv_LORT+pv_S)
      
      ! Find horizontal fluxes at specific points
      call cal_flux1(h_mflux, full_mass)
      call cal_flux1(h_tflux, full_temp)
    end subroutine comp_ijmin

    subroutine cal_flux1(h_flux, full_scalar)
      real(8), dimension(:), pointer :: h_flux
      real(8), dimension(0:N_BDRY) :: full_scalar

      h_flux(EDGE*(id+s )+UP+1) = 0.5_8*u_dual_dn*(full_scalar(SOUTH) + full_scalar(0))
      h_flux(EDGE*(id+sw)+DG+1) = 0.5_8*u_dual_sw*(full_scalar(0)     + full_scalar(SOUTHWEST))
      h_flux(EDGE*(id+w )+RT+1) = 0.5_8*u_dual_lt*(full_scalar(WEST)  + full_scalar(0))
    end subroutine cal_flux1

    subroutine cal_flux2(h_flux, full_scalar)
      real(8), dimension(:), pointer :: h_flux
      real(8), dimension(0:N_BDRY) :: full_scalar

      ! Find the horizontal mass flux as the velocity multiplied by the mass there
      h_flux(EDGE*id+UP+1) = 0.5_8*u_dual_up*(full_scalar(0)         + full_scalar(NORTH))
      h_flux(EDGE*id+DG+1) = 0.5_8*u_dual_dg*(full_scalar(NORTHEAST) + full_scalar(0))
      h_flux(EDGE*id+RT+1) = 0.5_8*u_dual_rt*(full_scalar(0)         + full_scalar(EAST))
    end subroutine cal_flux2

    subroutine comput()
      ! Computes physical quantities during upward integration
      type (Coord) :: n_e, x_i, vel
      real (8) :: kinetic_energy, Phi_k
     
      ! Find the velocity on primal and dual grid edges, which are equal except for the length of the
      ! side they are on
      u_prim_up = velo(EDGE*id+UP+1)*dom%len%elts(EDGE*id+UP+1)
      u_dual_up = velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
      u_prim_dg = velo(EDGE*id+DG+1)*dom%len%elts(EDGE*id+DG+1)
      u_dual_dg = velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
      u_prim_rt = velo(EDGE*id+RT+1)*dom%len%elts(EDGE*id+RT+1)
      u_dual_rt = velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)

      full_mass(0:NORTHEAST) = mass(id+(/0,n,e,s,w,ne/)+1) + mean(S_MASS,zlev)
      full_temp(0:NORTHEAST) = temp(id+(/0,n,e,s,w,ne/)+1) + mean(S_TEMP,zlev)

      ! Find horizontalfluxes at specific points
      call cal_flux2(h_mflux, full_mass)
      call cal_flux2(h_tflux, full_temp)

      ! Find additional primal and dual velocities: down, southwest (counter-diagonal), left
      u_prim_dn = velo(EDGE*(id+s)+UP+1)*dom%len%elts(EDGE*(id+s)+UP+1)
      u_dual_dn = velo(EDGE*(id+s)+UP+1)*dom%pedlen%elts(EDGE*(id+s)+UP+1)
      u_prim_sw = velo(EDGE*(id+sw)+DG+1)*dom%len%elts(EDGE*(id+sw)+DG+1)
      u_dual_sw = velo(EDGE*(id+sw)+DG+1)*dom%pedlen%elts(EDGE*(id+sw)+DG+1)
      u_prim_lt = velo(EDGE*(id+w)+RT+1)*dom%len%elts(EDGE*(id+w)+RT+1)
      u_dual_lt = velo(EDGE*(id+w)+RT+1)*dom%pedlen%elts(EDGE*(id+w)+RT+1)

      ! Calculate kinetic energy using Perot formula from equation (14) with approximate form (17) in Peixoto (2016)
      ! which gives first order convergence in maximum norm with Heikes-Randall (1995) optimized grids
      
      ! Sum contributions from all six edge velocities reconstructed at hexagon node x_i
      ! u_i = 1/2 A_i sum_e (u_e l_e) d_e n_e, where n_e is the normal vector to the hexagon edge e,
      ! l_e is the length of the hexagon edge (pedlen) and d_e is the length of the triangle edge (len)

      x_i = dom%node%elts(id+1)  ! Coordinate of node
      
      n_e = direction(x_i, dom%node%elts((id+e)+1))
      vel = vec_scale(u_dual_rt * dom%len%elts(EDGE*id+RT+1), n_e)

      n_e = direction(dom%node%elts((id+w)+1), x_i)
      vel = vec_sum(vel, vec_scale(u_dual_lt * dom%len%elts(EDGE*(id+w)+RT+1), n_e))
      
      n_e = direction(dom%node%elts((id+ne)+1), x_i)
      vel = vec_sum(vel, vec_scale(u_dual_dg * dom%len%elts(EDGE*id+DG+1), n_e))

      n_e = direction(x_i, dom%node%elts((id+sw)+1))
      vel = vec_sum(vel, vec_scale(u_dual_sw * dom%len%elts(EDGE*(id+sw)+DG+1), n_e))

      n_e = direction(x_i, dom%node%elts((id+n)+1))
      vel = vec_sum(vel, vec_scale(u_dual_up * dom%len%elts(EDGE*id+UP+1), n_e))

      n_e = direction(dom%node%elts((id+s)+1), x_i)
      vel = vec_sum(vel, vec_scale(u_dual_dn * dom%len%elts(EDGE*(id+s)+UP+1), n_e))

      ! Perot formula (16, 17) of Peixoto (2016)
      kinetic_energy = inner(vel,vel) * dom%areas%elts(id+1)%hex_inv**2 / 8.0_8
      
      ! Formula from TRiSK ... not convergent!
      ! kinetic_energy = &
      !      (u_prim_up*u_dual_up + u_prim_dg*u_dual_dg + u_prim_rt*u_dual_rt + &
      !      u_prim_dn*u_dual_dn + u_prim_sw*u_dual_sw + u_prim_lt*u_dual_lt &
      !      )* (1.0_8/4.0_8)*dom%areas%elts(id+1)%hex_inv

      ! Interpolate geopotential from interfaces to level
      Phi_k =  0.5_8*(dom%geopot%elts(id+1) + dom%geopot%elts(id+1))

      ! Bernoulli function
      if (compressible) then 
         bernoulli(id+1) = kinetic_energy + Phi_k
      else 
         bernoulli(id+1) = kinetic_energy + Phi_k + dom%press%elts(id+1)/ref_density
      end if

      ! Exner function in incompressible case from geopotential
      if (.not. compressible) exner(id+1) = -Phi_k

      if (viscosity .ne. 0.0_8) dom%divu%elts(id+1) = dom%areas%elts(id+1)%hex_inv * &
           (u_dual_up - u_dual_dg + u_dual_rt - u_dual_dn + u_dual_sw - u_dual_lt)

      dom%vort%elts(TRIAG*id+LORT+1) = - (u_prim_rt + u_prim_dg + velo(EDGE*(id+E)+UP+1)*dom%len%elts(EDGE*(id+E)+UP+1))
      dom%vort%elts(TRIAG*id+UPLT+1) =    u_prim_dg + u_prim_up + velo(EDGE*(id+N)+RT+1)*dom%len%elts(EDGE*(id+N)+RT+1) 

      vort_W = - ( &
           u_prim_lt + velo(EDGE*(id+W)+DG+1)*dom%len%elts(EDGE*(id+W)+DG+1) + u_prim_up)
      vort_S = &
           u_prim_rt + velo(EDGE*(id+S)+DG+1)*dom%len%elts(EDGE*(id+S)+DG+1) + u_prim_dn

      pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + dom%vort%elts(TRIAG*id+LORT+1))/( &
           full_mass(0)*dom%areas%elts(id   +1)%part(1) + &
           full_mass(EAST)*dom%areas%elts(id+E +1)%part(3) + &
           full_mass(NORTHEAST)*dom%areas%elts(id+NE+1)%part(5))

      pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + dom%vort%elts(TRIAG*id+UPLT+1))/( &
           full_mass(0)*dom%areas%elts(id   +1)%part(2) + &
           full_mass(NORTHEAST)*dom%areas%elts(id+NE+1)%part(4) + &
           full_mass(NORTH)*dom%areas%elts(id+N +1)%part(6))

      pv_W = (dom%coriolis%elts(TRIAG*(id+W)+LORT+1) + vort_W)/( &
           full_mass(WEST)*dom%areas%elts(id+W+1)%part(1) + &
           full_mass(0)*dom%areas%elts(id  +1)%part(3) + &
           full_mass(NORTH)*dom%areas%elts(id+N+1)%part(5))

      pv_S = (dom%coriolis%elts(TRIAG*(id+S)+UPLT+1) + vort_S)/( &
           full_mass(SOUTH)*dom%areas%elts(id+S+1)%part(2) + &
           full_mass(EAST)*dom%areas%elts(id+E+1)%part(4) + &
           full_mass(0)*dom%areas%elts(id  +1)%part(6))

      dom%vort%elts(TRIAG*id+LORT+1) = dom%vort%elts(TRIAG*id+LORT+1)/dom%triarea%elts(TRIAG*id+LORT+1) 
      dom%vort%elts(TRIAG*id+UPLT+1) = dom%vort%elts(TRIAG*id+UPLT+1)/dom%triarea%elts(TRIAG*id+UPLT+1) 

      qe(EDGE*id+RT+1) = 0.5_8*(pv_S    + pv_LORT)
      qe(EDGE*id+DG+1) = 0.5_8*(pv_UPLT + pv_LORT)
      qe(EDGE*id+UP+1) = 0.5_8*(pv_UPLT + pv_W)
    end subroutine comput
  end subroutine step1

  subroutine interp_vel_hex (dom, i, j, zlev, offs, dims)
    ! Interpolate velocity to hexagon nodes in Cartesian coordinates; uses Perot formula as also used for kinetic energy 
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    type(Coord) :: n_e, vel, x_i
    type(Coord) :: e_zonal, e_merid
    integer     :: id, idN, idE, idNE, idS, idSW, idW
    real(8)     :: u_dual_up, u_dual_dg, u_dual_rt, u_dual_dn, u_dual_sw, u_dual_lt
    real(8)     :: lon, lat

    id   = idx(i,     j,     offs, dims)
    idN  = idx(i, j + 1,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)

    u_dual_up = velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
    u_dual_dg = velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_rt = velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)

    u_dual_dn = velo(EDGE*idS+UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)
    u_dual_sw = velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_lt = velo(EDGE*idW+RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)

    x_i = dom%node%elts(id+1)
    
    n_e = direction(x_i, dom%node%elts(idE+1))
    vel = vec_scale(u_dual_rt * dom%len%elts(EDGE*id+RT+1), n_e)

    n_e = direction(dom%node%elts(idW+1), x_i)
    vel = vec_sum(vel, vec_scale(u_dual_lt * dom%len%elts(EDGE*idW+RT+1), n_e))
    
    n_e = direction(dom%node%elts(idNE+1), x_i)
    vel = vec_sum(vel, vec_scale(u_dual_dg * dom%len%elts(EDGE*id+DG+1), n_e))
    
    n_e = direction(x_i, dom%node%elts(idSW+1))
    vel = vec_sum(vel, vec_scale(u_dual_sw * dom%len%elts(EDGE*idSW+DG+1), n_e))
    
    n_e = direction(x_i, dom%node%elts(idN+1))
    vel = vec_sum(vel, vec_scale(u_dual_up * dom%len%elts(EDGE*id+UP+1), n_e))
    
    n_e = direction(dom%node%elts(idS+1), x_i)
    vel = vec_sum(vel, vec_scale(u_dual_dn * dom%len%elts(EDGE*idS+UP+1), n_e))

    ! Perot formula (16, 17) of Peixoto (2016) gives velocity at hexagonal node in Cartesian coordinates
    vel = vec_scale(0.5_8*dom%areas%elts(id+1)%hex_inv, vel)

    ! Project velocity onto zonal and meridional directions
    call cart2sph (x_i, lon, lat)
    
    e_zonal = Coord (-sin(lon),           cos(lon),           0.0_8)   ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    dom%u_zonal%elts(id+1) = inner(vel, e_zonal)
    dom%v_merid%elts(id+1) = inner(vel, e_merid)
  end subroutine interp_vel_hex
  
  subroutine vert_integrated_horiz_flux (dom, i, j, zlev, offs, dims)
    ! Integrate horizontal fluxes on the three edges vertically 
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id, k

    id   = idx(i, j, offs, dims)

    dom%integr_horiz_flux%elts(EDGE*id+UP+1) = 0.0_8
    dom%integr_horiz_flux%elts(EDGE*id+DG+1) = 0.0_8
    dom%integr_horiz_flux%elts(EDGE*id+RT+1) = 0.0_8
    
    do k = 1, zlevels
       dom%integr_horiz_flux%elts(EDGE*id+UP+1) = dom%integr_horiz_flux%elts(EDGE*id+UP+1) + &
            horiz_flux(S_MASS,k)%data(dom%id+1)%elts(EDGE*id+UP+1)
       
       dom%integr_horiz_flux%elts(EDGE*id+DG+1) = dom%integr_horiz_flux%elts(EDGE*id+DG+1) + &
            horiz_flux(S_MASS,k)%data(dom%id+1)%elts(EDGE*id+DG+1)
       
       dom%integr_horiz_flux%elts(EDGE*id+RT+1) = dom%integr_horiz_flux%elts(EDGE*id+RT+1) + &
            horiz_flux(S_MASS,k)%data(dom%id+1)%elts(EDGE*id+RT+1)
    end do
  end subroutine vert_integrated_horiz_flux

  subroutine compute_vert_flux (dom, i, j, zlev, offs, dims)
    ! Computes vertical mass flux at upper interface of level zlev from mass trend and divergence of horizontal flux
    ! when using mass-based vertical coordinates
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    
    integer :: id, k

    id   = idx(i, j, offs, dims)

    if (zlev.eq.1) then 
       v_mflux(id+1) = 0.0_8                   - dmass(id+1) + horiz_div_flux(h_mflux, dom, i, j, offs, dims, id) ! Flux is zero at surface
    elseif (zlev.eq.zlevels) then ! Flux is zero at upper interface of top level
       v_mflux(id+1) = 0.0_8
    else
       v_mflux(id+1) = dom%adj_mass%elts(id+1) - dmass(id+1) + horiz_div_flux(h_mflux, dom, i, j, offs, dims, id)
    end if

    ! Save current vertical mass flux for lower interface of next vertical level (use adj_mass)
    dom%adj_mass%elts(id+1) = v_mflux(id+1)
  end subroutine compute_vert_flux
  
  subroutine integrate_pressure_up(dom, i, j, zlev, offs, dims)
    ! Integrate pressure (compressible case)/Lagrange multiplier (incompressible case) and geopotential up from surface to top layer
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id
    real(8) :: full_mass, full_temp, full_pot_temp, spec_vol, pert_spec_vol

    id = idx(i, j, offs, dims)

    full_mass     = mass(id+1) + mean(S_MASS,zlev)
    full_temp     = temp(id+1) + mean(S_TEMP,zlev)
    full_pot_temp = full_temp/full_mass

    if (full_mass .lt. 1e-6_8) then
       print *, 'fatal error: a horizontal layer thickness is being squeezed to zero, namely, at zlev=', zlev
       write(6,*) 'mean+perturbation=', full_mass
       stop
    end if

    if (compressible) then ! Compressible case
       if (zlev .eq. 1) then
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*full_mass
       else ! Interpolate mass=rho*dz to lower interface of current level
          dom%press%elts(id+1) = dom%press%elts(id+1) - 0.5_8*grav_accel*(full_mass + dom%adj_mass%elts(id+1))
       end if
       dom%adj_mass%elts(id+1) = full_mass ! Save current mass for pressurce calculation at next vertical level

       if (zlev .eq. zlevels) then !top zlev, purely diagnostic
          if (abs((dom%press%elts(id+1) - 0.5_8*grav_accel*full_mass) - press_infty) .gt. 1e-10_8) then
             write(6,*) 'warning: upward integration of pressure not resulting in zero at top interface'
             write(6,*) 'observed pressure - pressure_infty =', abs((dom%press%elts(id+1) &
                  - 0.5_8*grav_accel*full_mass) - press_infty)
             stop
          end if
       end if

       ! Exner function from pressure
       exner(id+1) = c_p*(dom%press%elts(id+1)/ref_press)**kappa
       
       ! Specific volume alpha = kappa*theta*pi/p
       spec_vol = kappa * full_pot_temp * exner(id+1) / dom%press%elts(id+1)
       pert_spec_vol = spec_vol - mean_spec_vol(zlev)

       ! Find geopotential at upper interface of current level using (18) in DYNAMICO
       if (zlev .eq. 1) then ! Save geopotential at lower interface of level zlev for interpolation in Bernoulli function
          dom%adj_geopot%elts(id+1) = dom%surf_geopot%elts(id+1) 
       else
          dom%adj_geopot%elts(id+1) = dom%geopot%elts(id+1)
       end if
       dom%geopot%elts(id+1) = dom%adj_geopot%elts(id+1) + &
            grav_accel * (full_mass*pert_spec_vol + mass(id+1)*mean_spec_vol(zlev))
       
    else ! Incompressible case
       if (zlev .eq. 1) then 
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*full_temp
       else ! Interpolate to lower interface of current level
          dom%press%elts(id+1) = dom%press%elts(id+1) - 0.5_8*grav_accel*(dom%adj_temp%elts(id+1) + full_temp)
       end if
       dom%adj_mass%elts(id+1) = full_mass

       if (zlev .eq. zlevels) then !top zlev, purely diagnostic
          if (abs(dom%press%elts(id+1)-0.5_8*grav_accel*full_temp - press_infty).gt. 1e-10_8) then
             print *, 'warning: upward integration of Lagrange multiplier not resulting in zero at top interface'
             write(6,'(A,es15.8)') "Pressure at infinity = ", dom%press%elts(id+1)-0.5_8*grav_accel*full_temp
             write(6,'(A,es15.8)') "Press_infty = ", press_infty
              write(6,'(A,es15.8)') "Difference = ", dom%press%elts(id+1)-0.5_8*grav_accel*full_temp - press_infty
             stop
          end if
       end if

       ! Find geopotential at interfaces using (18) in DYNAMICO
       ! Note: since mu is associated with the kinematic mass = inert mass (not the gravitational mass defined by the buyoancy)
       ! we divide by the constant reference density. This is the Boussinesq approximation.
       if (zlev .eq. 1) then ! Save geopotential at lower interface of level zlev for interpolation in Bernoulli function
          dom%adj_geopot%elts(id+1) = dom%surf_geopot%elts(id+1) 
       else
          dom%adj_geopot%elts(id+1) = dom%geopot%elts(id+1)
       end if
       dom%geopot%elts(id+1) = dom%adj_geopot%elts(id+1) + grav_accel* mass(id+1)/ref_density
    end if
  end subroutine integrate_pressure_up

  subroutine integrate_pressure_down(dom, i, j, zlev, offs, dims)
    ! Pressure is computed during downward integration from zlev=zlevels to zlev=1
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id
    real(8) :: full_mass, full_temp

    id   = idx(i, j, offs, dims)

    full_mass = mass(id+1) + mean(S_MASS,zlev)
    full_temp = temp(id+1) + mean(S_TEMP,zlev)

    if (compressible) then !compressible case
       !integrate (or, rather, interpolate) the pressure from top zlev down to bottom zlev; press_infty is user-set
       if (zlev .eq. zlevels) then
          dom%press%elts(id+1) = press_infty + 0.5_8*grav_accel*full_mass
       else ! Interpolate mass to lower interface
          dom%press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*(dom%adj_mass%elts(id+1) + full_mass)
       end if

       !surface pressure is set (even at t=0) from downward numerical integration
       if (zlev .eq. 1) dom%surf_press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*full_mass
    else !incompressible case
       !integrate (or, rather, interpolate) the pressure from top zlev down to bottom zlev; press_infty is user-set
       if (zlev .eq. zlevels) then !top zlev, it is an exception
          dom%press%elts(id+1) = press_infty + 0.5_8*grav_accel*full_temp
       else !other layers equal to half of previous layer and half of current layer
          dom%press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*(dom%adj_temp%elts(id+1) + full_temp)
       end if

       !surface pressure is set (even at t=0) from downward numerical integration
       if (zlev .eq. 1) dom%surf_press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*full_temp
    end if

    !quantities for vertical integration in next zlev
    dom%adj_mass%elts(id+1) = full_mass
    dom%adj_temp%elts(id+1) = full_temp
  end subroutine integrate_pressure_down

  function get_weights(dom, id, offs)
    !find weights for Qperp computation [Aechtner thesis page 44]
    type(Domain) dom
    real(8) get_weights(5)
    real(8) wgt(5)
    integer id, offs, i

    wgt(1) = dom%areas%elts(id+1)%part(1+offs)
    do i = 2, 5
       wgt(i) = wgt(i-1) + dom%areas%elts(id+1)%part(modulo(i+offs-1,6)+1)
    end do
    wgt = 0.5_8 - wgt*dom%areas%elts(id+1)%hex_inv
    get_weights = (/wgt(1), -wgt(2), wgt(3), -wgt(4), wgt(5)/)
  end function get_weights

  subroutine du_Qperp_Enstrophy(dom, i, j, zlev, offs, dims)
    !compute enstrophy-conserving Qperp and store it in dvelo [Aechtner thesis page 44]
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer idNW
    integer idN
    integer idNE
    integer idW
    integer id
    integer idE
    integer idSW
    integer idS
    integer idSE
    real(8) wgt1(5), wgt2(5)

    idNW = idx(i - 1, j + 1, offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idSE = idx(i + 1, j - 1, offs, dims)

    wgt1 = get_weights(dom, id, 0)
    wgt2 = get_weights(dom, idE, 3)

    dvelo(EDGE*id+RT+1) = qe(EDGE*id+RT+1)*( &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)&
         %elts((/ EDGE*id+DG, EDGE*id+UP, EDGE*idW+RT, EDGE*idSW+DG, EDGE*idS+UP /)+1) * wgt1) + &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)&
         %elts((/ EDGE*idS+DG, EDGE*idSE+UP, EDGE*idE+RT, EDGE*idE+DG, EDGE*idE+UP /)+1) * wgt2))

    wgt1 = get_weights(dom, id, 1)
    wgt2 = get_weights(dom, idNE, 4)

    dvelo(EDGE*id+DG+1) = qe(EDGE*id+DG+1)*( &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*id+UP, EDGE*idW+RT, EDGE*idSW+DG, EDGE*idS+UP, EDGE*id+RT/)+1) * wgt1) + &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*idE+UP, EDGE*idNE+RT, DG+EDGE*idNE, EDGE*idNE+UP, EDGE*idN+RT/)+1) * wgt2))

    wgt1 = get_weights(dom, id, 2)
    wgt2 = get_weights(dom, idN, 5)

    dvelo(EDGE*id+UP+1) = qe(EDGE*id+UP+1)*( &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*idW+RT, EDGE*idSW+DG, EDGE*idS+UP, EDGE*id+RT, EDGE*id+DG /)+1) * wgt1) + &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*idN+RT, EDGE*idN+DG, EDGE*idN+UP, EDGE*idNW+RT, EDGE*idW+DG /)+1) * wgt2))
  end subroutine du_Qperp_Enstrophy

  subroutine du_source(dom, i, j, zlev, offs, dims)
    ! Source (non gradient) terms in velocity trend
    ! [Aechtner thesis page 56, Kevlahan, Dubos and Aechtner (2015)]
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id

    id = idx(i, j, offs, dims)

    call du_Qperp(dom, i, j, zlev, offs, dims)

    if (viscosity .ne. 0.0_8) call diff_mom(dom, i, j, zlev, offs, dims)
  end subroutine du_source

  subroutine du_Qperp(dom, i, j, zlev, offs, dims)
    ! Compute energy-conserving Qperp and add it to dvelo [Aechtner thesis page 44]
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id, idNW, idN, idNE, idW, idE, idSW, idS, idSE
    real(8) :: wgt1(5), wgt2(5)

    idNW = idx(i - 1, j + 1, offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idSE = idx(i + 1, j - 1, offs, dims)

    wgt1 = get_weights(dom, id, 0)
    wgt2 = get_weights(dom, idE, 3)

    dvelo(EDGE*id+RT+1) = &
         h_mflux(EDGE*id+DG+1)*0.5_8*(qe(EDGE*id+DG+1) + &
         qe(EDGE*id+RT+1))*wgt1(1) + &
         h_mflux(EDGE*id+UP+1)*0.5_8*(qe(EDGE*id+UP+1) + &
         qe(EDGE*id+RT+1))*wgt1(2) + &
         h_mflux(EDGE*idW+RT+1)*0.5_8*(qe(EDGE*idW+RT+1) + &
         qe(EDGE*id+RT+1))*wgt1(3) + &
         h_mflux(EDGE*idSW+DG+1)*0.5_8*(qe(EDGE*idSW+DG+1) &
         + qe(EDGE*id+RT+1))*wgt1(4) + &
         h_mflux(EDGE*idS+UP+1)*0.5_8*(qe(EDGE*idS+UP+1) + &
         qe(EDGE*id+RT+1))*wgt1(5) + &
         h_mflux(DG+EDGE*idS+1)*0.5_8*(qe(DG+EDGE*idS+1) + &
         qe(EDGE*id+RT+1))*wgt2(1) + &
         h_mflux(EDGE*idSE+UP+1)*0.5_8*(qe(EDGE*idSE+UP+1) &
         + qe(EDGE*id+RT+1))*wgt2(2) + &
         h_mflux(EDGE*idE+RT+1)*0.5_8*(qe(EDGE*idE+RT+1) + &
         qe(EDGE*id+RT+1))*wgt2(3) + &
         h_mflux(DG+EDGE*idE+1)*0.5_8*(qe(DG+EDGE*idE+1) + &
         qe(EDGE*id+RT+1))*wgt2(4) + &
         h_mflux(EDGE*idE+UP+1)*0.5_8*(qe(EDGE*idE+UP+1) + &
         qe(EDGE*id+RT+1))*wgt2(5)

    wgt1 = get_weights(dom, id, 1)
    wgt2 = get_weights(dom, idNE, 4)

    dvelo(EDGE*id+DG+1) = &
         h_mflux(EDGE*id+UP+1)*0.5_8*(qe(EDGE*id+UP+1) + &
         qe(EDGE*id+DG+1))*wgt1(1) + &
         h_mflux(EDGE*idW+RT+1)*0.5_8*(qe(EDGE*idW+RT+1) + &
         qe(EDGE*id+DG+1))*wgt1(2) + &
         h_mflux(EDGE*idSW+DG+1)*0.5_8*(qe(EDGE*idSW+DG+1) &
         + qe(EDGE*id+DG+1))*wgt1(3) + &
         h_mflux(EDGE*idS+UP+1)*0.5_8*(qe(EDGE*idS+UP+1) + &
         qe(EDGE*id+DG+1))*wgt1(4) + &
         h_mflux(EDGE*id+RT+1)*0.5_8*(qe(EDGE*id+RT+1) + &
         qe(EDGE*id+DG+1))*wgt1(5) + &
         h_mflux(EDGE*idE+UP+1)*0.5_8*(qe(EDGE*idE+UP+1) + &
         qe(EDGE*id+DG+1))*wgt2(1) + &
         h_mflux(EDGE*idNE+RT+1)*0.5_8*(qe(EDGE*idNE+RT+1) &
         + qe(EDGE*id+DG+1))*wgt2(2) + &
         h_mflux(EDGE*idNE+DG+1)*0.5_8*(qe(EDGE*idNE+DG+1) &
         + qe(EDGE*id+DG+1))*wgt2(3) + &
         h_mflux(EDGE*idNE+UP+1)*0.5_8*(qe(EDGE*idNE+UP+1) &
         + qe(EDGE*id+DG+1))*wgt2(4) + &
         h_mflux(EDGE*idN+RT+1)*0.5_8*(qe(EDGE*idN+RT+1) + &
         qe(EDGE*id+DG+1))*wgt2(5)

    wgt1 = get_weights(dom, id, 2)
    wgt2 = get_weights(dom, idN, 5)

    dvelo(EDGE*id+UP+1) = &
         h_mflux(EDGE*idW+RT+1)*0.5_8*(qe(EDGE*idW+RT+1) + &
         qe(EDGE*id+UP+1))*wgt1(1) + &
         h_mflux(EDGE*idSW+DG+1)*0.5_8*(qe(EDGE*idSW+DG+1) &
         + qe(EDGE*id+UP+1))*wgt1(2) + &
         h_mflux(EDGE*idS+UP+1)*0.5_8*(qe(EDGE*idS+UP+1) + &
         qe(EDGE*id+UP+1))*wgt1(3) + &
         h_mflux(EDGE*id+RT+1)*0.5_8*(qe(EDGE*id+RT+1) + &
         qe(EDGE*id+UP+1))*wgt1(4) + &
         h_mflux(EDGE*id+DG+1)*0.5_8*(qe(EDGE*id+DG+1) + &
         qe(EDGE*id+UP+1))*wgt1(5) + &
         h_mflux(EDGE*idN+RT+1)*0.5_8*(qe(EDGE*idN+RT+1) + &
         qe(EDGE*id+UP+1))*wgt2(1) + &
         h_mflux(EDGE*idN+DG+1)*0.5_8*(qe(EDGE*idN+DG+1) + &
         qe(EDGE*id+UP+1))*wgt2(2) + &
         h_mflux(EDGE*idN+UP+1)*0.5_8*(qe(EDGE*idN+UP+1) + &
         qe(EDGE*id+UP+1))*wgt2(3) + &
         h_mflux(EDGE*idNW+RT+1)*0.5_8*(qe(EDGE*idNW+RT+1) &
         + qe(EDGE*id+UP+1))*wgt2(4) + &
         h_mflux(EDGE*idW+DG+1)*0.5_8*(qe(EDGE*idW+DG+1) + &
         qe(EDGE*id+UP+1))*wgt2(5)
  end subroutine du_Qperp

  subroutine scalar_trend(dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: v_tflux, full_pot_temp, full_pot_temp_up

    id   = idx(i, j, offs, dims)

    if (lagrangian_vertical) then
       dmass(id+1) = horiz_div_flux(h_mflux, dom, i, j, offs, dims, id)
       dtemp(id+1) = horiz_div_flux(h_tflux, dom, i, j, offs, dims, id) 
    else ! Mass-based vertical coordinates
       
       ! Compute mass trend (mu_t) at level zlev from total mass trend M_t
       ! (our definition of a_vert, b_vert reversed compared with Dynamico paper)
       dmass(id+1) = (b_vert(zlev+1)-b_vert(zlev)) * horiz_div_flux(h_mflux, dom, i, j, offs, dims, id)

       ! Compute horizontal divergence of horizontal temperature flux
       dtemp(id+1) = horiz_div_flux(h_tflux, dom, i, j, offs, dims, id)

       ! Compute vertical mass flux v_mflux at upper interface of level zlev
       call compute_vert_flux (dom, i, j, zlev, offs, dims)

       ! Find non mass-weighted potential temperature at current level
       full_pot_temp = (temp(id+1) + mean(S_TEMP,zlev))/(mass(id+1) + mean(S_MASS,zlev))

       ! Add vertical divergence of vertical potential temperature flux to
       ! trend of mass-weighted potential temperature
       if (zlev.eq.1) then  ! Vertical flux is zero at surface
          ! Potential temperature at next level up
          full_pot_temp_up = (adj_temp_up(id+1) + mean(S_TEMP,zlev))/(adj_mass_up(id+1) + mean(S_MASS,zlev))

          ! Vertical flux of potential temperature at upper interface
          v_tflux = 0.5_8*(full_pot_temp + full_pot_temp_up) * v_mflux(id+1)

          dtemp(id+1) = dtemp(id+1) + v_tflux
          
       elseif (zlev.eq.zlevels) then ! Vertical flux is zero at top interface

          dtemp(id+1) = dtemp(id+1) - dom%adj_vflux%elts(id+1)
          
       else
          ! Potential temperature at next level up
          full_pot_temp_up = (adj_temp_up(id+1) + mean(S_TEMP,zlev))/(adj_mass_up(id+1) + mean(S_MASS,zlev))
          
          ! Vertical flux of potential temperature at upper interface
          v_tflux = 0.5_8*(full_pot_temp + full_pot_temp_up) * v_mflux(id+1)

          dtemp(id+1) = dtemp(id+1) + (v_tflux - dom%adj_vflux%elts(id+1))
          
       end if

       ! Save vertical flux of potential temperature for lower interface of next vertical level
       dom%adj_vflux%elts(id+1) = v_tflux
    end if
  end subroutine scalar_trend

  function horiz_div_flux(h_flux, dom, i, j, offs, dims, id)
    ! Computes negative of divergence of horizontal flux h_flux over hexagons -delta_i(U)
    real(8)                        :: horiz_div_flux
    real(8), dimension(:), pointer :: h_flux
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW

    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)

    horiz_div_flux = -(h_flux(EDGE*id+UP+1)   - h_flux(EDGE*id+DG+1) + &
         h_flux(EDGE*id+RT+1)   - h_flux(EDGE*idS+UP+1) + &
         h_flux(EDGE*idSW+DG+1) - h_flux(EDGE*idW+RT+1)) &
         *dom%areas%elts(id+1)%hex_inv
  end function horiz_div_flux

  subroutine du_gradB_gradExn (dom, i, j, zlev, offs, dims)
    ! Add gradients of Bernoulli and Exner to dvelo [DYNAMICO (23)-(25)]
    ! mass and potential temperature trend is zero
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer               :: e, id, idE, idN, idNE
    real(8)               :: full_pot_temp(0:N_BDRY)
    real(8), dimension(3) :: v_star

    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    ! See DYNAMICO between (23)-(25), geopotential still known from step1_upw
    ! the theta multiplying the Exner gradient is the edge-averaged non-mass-weighted potential temperature

    full_pot_temp(0)         = (temp(id+1)   + mean(S_TEMP,zlev))/(mass(id+1)   + mean(S_MASS,zlev))
    full_pot_temp(NORTH)     = (temp(idN+1)  + mean(S_TEMP,zlev))/(mass(idN+1)  + mean(S_MASS,zlev))
    full_pot_temp(EAST)      = (temp(idE+1)  + mean(S_TEMP,zlev))/(mass(idE+1)  + mean(S_MASS,zlev))
    full_pot_temp(NORTHEAST) = (temp(idNE+1) + mean(S_TEMP,zlev))/(mass(idNE+1) + mean(S_MASS,zlev))

    if (compressible) then
       dvelo(EDGE*id+RT+1) = (dvelo(EDGE*id+RT+1) - (bernoulli(idE+1) - bernoulli(id+1)) &
            - 0.5_8*(full_pot_temp(0)+full_pot_temp(EAST)) * (exner(idE+1) - exner(id+1)))/dom%len%elts(EDGE*id+RT+1)

       dvelo(EDGE*id+DG+1) = (dvelo(EDGE*id+DG+1) - (bernoulli(id+1) - bernoulli(idNE+1)) &
            - 0.5_8*(full_pot_temp(0)+full_pot_temp(NORTHEAST)) * (exner(id+1) - exner(idNE+1)))/dom%len%elts(EDGE*id+DG+1)

       dvelo(EDGE*id+UP+1) = (dvelo(EDGE*id+UP+1) - (bernoulli(idN+1) - bernoulli(id+1)) &
            - 0.5_8*(full_pot_temp(0)+full_pot_temp(NORTH)) * (exner(idN+1) - exner(id+1)))/dom%len%elts(EDGE*id+UP+1)
    else ! Incompressible case
       dvelo(EDGE*id+RT+1) = (dvelo(EDGE*id+RT+1) - (bernoulli(idE+1) - bernoulli(id+1)) &
            - 0.5_8*(2.0_8-full_pot_temp(0)-full_pot_temp(EAST)) * (exner(idE+1) - exner(id+1)))/dom%len%elts(EDGE*id+RT+1)

       dvelo(EDGE*id+DG+1) = (dvelo(EDGE*id+DG+1) - (bernoulli(id+1) - bernoulli(idNE+1)) &
            - 0.5_8*(2.0_8-full_pot_temp(0)-full_pot_temp(NORTHEAST)) * (exner(id+1) - exner(idNE+1)))/dom%len%elts(EDGE*id+DG+1)

       dvelo(EDGE*id+UP+1) = (dvelo(EDGE*id+UP+1) - (bernoulli(idN+1) - bernoulli(id+1)) &
            - 0.5_8*(2.0_8-full_pot_temp(0)-full_pot_temp(NORTH)) * (exner(idN+1) - exner(id+1)))/dom%len%elts(EDGE*id+UP+1)
    end if

    ! Add vertical flux gradient term
    if (.not. lagrangian_vertical) then

       ! Assume free slip boundary conditions at top and bottom of vertical layers (i.e. velocity at top interface and surface equals
       ! velocity at adjacent full level)
       
       if (zlev.eq.1) then ! surface
          ! Horizontal velocities at edges interpolated at upper interface
          do e = 1, 3
             v_star(e) = 0.5_8*(velo(EDGE*id+e) + adj_velo_up(EDGE*id+e))
          end do
          
          dvelo(EDGE*id+RT+1) = dvelo(EDGE*id+RT+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idE+1))  * &
               (v_star(RT+1) - velo(EDGE*id+RT+1))
          
          dvelo(EDGE*id+DG+1) = dvelo(EDGE*id+DG+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idNE+1)) * &
               (v_star(DG+1) - velo(EDGE*id+DG+1))
          
          dvelo(EDGE*id+UP+1) = dvelo(EDGE*id+UP+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idN+1))  * &
               (v_star(UP+1) - velo(EDGE*id+UP+1))

          ! Save horizontal velocities (needed for lower interface at next vertical level)
          do e = 1, 3
             dom%adj_velo%elts(EDGE*id+e) = v_star(e)
          end do
       elseif (zlev.eq.zlevels) then ! top level
          dvelo(EDGE*id+RT+1) = dvelo(EDGE*id+RT+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idE+1))  * &
               (velo(EDGE*id+RT+1) - dom%adj_velo%elts(EDGE*id+RT+1))

          dvelo(EDGE*id+DG+1) = dvelo(EDGE*id+DG+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idNE+1)) * &
               (velo(EDGE*id+DG+1) - dom%adj_velo%elts(EDGE*id+DG+1))
          
          dvelo(EDGE*id+UP+1) = dvelo(EDGE*id+UP+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idN+1))  * &
               (velo(EDGE*id+UP+1) - dom%adj_velo%elts(EDGE*id+UP+1))
       else
          ! Horizontal velocities at edges interpolated at upper interface
          do e = 1, 3
             v_star(e) = 0.5_8*(velo(EDGE*id+e) + adj_velo_up(EDGE*id+e))
          end do
          
          dvelo(EDGE*id+RT+1) = dvelo(EDGE*id+RT+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idE+1)) * &
               (v_star(RT+1) - dom%adj_velo%elts(EDGE*id+RT+1))

          dvelo(EDGE*id+DG+1) = dvelo(EDGE*id+DG+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idNE+1)) * &
               (v_star(DG+1) - dom%adj_velo%elts(EDGE*id+DG+1))
          
          dvelo(EDGE*id+UP+1) = dvelo(EDGE*id+UP+1) + 0.5_8*(dom%vert_velo%elts(id+1) + dom%vert_velo%elts(idN+1)) * &
               (v_star(UP+1) - dom%adj_velo%elts(EDGE*id+UP+1))

          ! Save horizontal velocities interpolated at upper interface (needed for lower interface at next vertical level)
          do e = 1, 3
             dom%adj_velo%elts(EDGE*id+e) = v_star(e)
          end do
       end if
    end if
  end subroutine du_gradB_gradExn

  subroutine interp_vert_velo_at_full_levels (dom, i, j, zlev, offs, dims)
    ! Interpolate non mass-weighted vertical velocity at full levels from vertical velocity at interfaces 
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id
    real (8) :: velo
    
    id   = idx(i, j, offs, dims)

    ! Non-mass weighted vertical velocity at upper interface
    velo = v_mflux(id+1)/(mass(id+1) + mean(S_MASS,zlev))
    
    if (zlev.eq.1) then ! No vertical flux through lower interface of bottom level
       dom%vert_velo%elts(id+1) = 0.5_8 * velo
    else
       dom%vert_velo%elts(id+1) = 0.5_8 * (velo + dom%adj_vflux%elts(id+1))
    end if
    
    ! Save current vertical velocity at upper interface for lower interface of next vertical level (use adj_vflux)
    dom%adj_vflux%elts(id+1) = velo
  end subroutine interp_vert_velo_at_full_levels

  subroutine sum_mass_temp(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id

    id   = idx(i,     j,     offs, dims)

    sum_mass = sum_mass + abs(mass(id+1))
    sum_temp = sum_temp + abs(temp(id+1))
  end subroutine sum_mass_temp

  subroutine comp_offs3(dom, p, offs, dims)
    type(Domain) dom
    integer p 
    integer offs(0:N_BDRY), dims(2,N_BDRY), i, n

    offs(0) = dom%patch%elts(p+1)%elts_start
    do i = 1, N_BDRY
       n = dom%patch%elts(p+1)%neigh(i)
       if (n .gt. 0) then ! regular patch
          offs(i)  = dom%patch%elts(n+1)%elts_start
          dims(:,i) = PATCH_SIZE
       else if (n .lt. 0) then
          offs(i)  = dom%bdry_patch%elts(-n+1)%elts_start
          dims(:,i) = sides_dims(:,abs(dom%bdry_patch%elts(-n+1)%side)+1)
       end if
    end do

    offs(1:N_BDRY)  = offs(1:N_BDRY)  - offs(0) ! make relative
    offs(SOUTH)     = offs(SOUTH)     + dims(1,SOUTH)*(dims(2,SOUTH)-1)
    offs(SOUTHEAST) = offs(SOUTHEAST) + dims(1,SOUTHEAST)*(dims(2,SOUTH)-1)
    offs(WEST)      = offs(WEST)      + (dims(1,WEST)-1)
    offs(NORTHWEST) = offs(NORTHWEST) + (dims(1,NORTHWEST)-1)
    offs(SOUTHWEST) = offs(SOUTHWEST) + dims(1,SOUTHWEST)*dims(2,SOUTHWEST) - 1
    offs(NORTH)     = offs(NORTH)     - PATCH_SIZE*(PATCH_SIZE-1) 
    offs(NORTHWEST) = offs(NORTHWEST) - PATCH_SIZE*(PATCH_SIZE-1) 
    offs(EAST)      = offs(EAST)      - (PATCH_SIZE-1)
    offs(SOUTHEAST) = offs(SOUTHEAST) - (PATCH_SIZE-1)
    offs(NORTHEAST) = offs(NORTHEAST) - (PATCH_SIZE*PATCH_SIZE-1)
  end subroutine comp_offs3

  subroutine sum_dmassdtemp(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idS
    integer idW
    integer idSW

    id   = idx(i,     j,     offs, dims)

    totaldmass = totaldmass + dmass(id+1)/dom%areas%elts(id+1)%hex_inv
    totalabsdmass = totalabsdmass + abs(dmass(id+1)/dom%areas%elts(id+1)%hex_inv)

    totaldtemp = totaldtemp + dtemp(id+1)/dom%areas%elts(id+1)%hex_inv
    totalabsdtemp = totalabsdtemp + abs(dtemp(id+1)/dom%areas%elts(id+1)%hex_inv)
  end subroutine sum_dmassdtemp
end module ops_mod
