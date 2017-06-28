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

  subroutine post_step1(dom, p, c, offs, dims, zlev)
    type(Domain) dom
    integer p
    integer c
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,9) :: dims
    integer id, idS, idW, idSW, idN, idE, idNE
    real(8) pv_SW, pv_W, pv_S, pv_LORT, pv_UPLT, pv_SW_LORT, pv_SW_UPLT, pv
    real(8) phi(0:N_BDRY), full_mass(0:N_BDRY)

    phi(0:NORTHEAST) = 1.0

    if (c .eq. IJMINUS) then
       id   = idx( 0,  0, offs, dims)
       idSW = idx(-1, -1, offs, dims)
       idW  = idx(-1,  0, offs, dims)
       idS  = idx( 0, -1, offs, dims)
       idN  = idx( 0,  1, offs, dims)
       idE  = idx( 1,  0, offs, dims)

       if (penalize) phi(0:WEST) = phi(0:WEST) + alpha_m1*penal%data(dom%id+1)%elts((/id,idN,idE,idS,idW/)+1)

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

       dom%qe%elts(EDGE*idW+RT+1) = 0.5_8*(pv_W + pv_SW)
       dom%qe%elts(EDGE*idS+UP+1) = 0.5_8*(pv_S + pv_SW)
    end if

    if (c .eq. IPLUSJMINUS) then
       id   = idx(PATCH_SIZE,    0, offs, dims)
       idSW = idx(PATCH_SIZE-1, -1, offs, dims)
       idS  = idx(PATCH_SIZE,   -1, offs, dims)
       idW  = idx(PATCH_SIZE-1,  0, offs, dims)
       idE  = idx(PATCH_SIZE+1,  0, offs, dims)
       idNE = idx(PATCH_SIZE+1,  1, offs, dims)

       if (penalize) phi(0:NORTHEAST) = phi(0:NORTHEAST) + alpha_m1*penal%data(dom%id+1)%elts((/id,id,idE,idS,idW,idNE/)+1)

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

       dom%qe%elts(EDGE*id  +RT+1) = 0.5_8*(pv_S + pv_LORT)
       dom%qe%elts(EDGE*idSW+DG+1) = 0.5_8*(pv_SW_LORT + pv_SW_UPLT)
    end if

    if (c .eq. IMINUSJPLUS) then
       id   = idx(0,  PATCH_SIZE,   offs, dims)
       idSW = idx(-1, PATCH_SIZE-1, offs, dims)
       idW  = idx(-1, PATCH_SIZE,   offs, dims)
       idS  = idx(0,  PATCH_SIZE-1, offs, dims)
       idN  = idx(0,  PATCH_SIZE+1, offs, dims)
       idNE = idx(1,  PATCH_SIZE+1, offs, dims)

       if (penalize) phi(0:NORTHEAST) = phi(0:NORTHEAST) + alpha_m1*penal%data(dom%id+1)%elts((/id,idN,id,idS,idW,idNE/)+1)

       full_mass(0:NORTHEAST) = mass((/id,idN,id,idS,idW,idNE/)+1) + mean(S_MASS,zlev)
       
       dom%vort%elts(TRIAG*idSW+UPLT+1) = &
            - velo(EDGE*id+UP+1)*dom%len%elts(EDGE*id+UP+1) &
            + velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1) &
            + velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)

       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) &
            + dom%vort%elts(TRIAG*idSW+UPLT+1)) &
            /(full_mass(SOUTHWEST)*dom%areas%elts(idSW  +1)%part(2) &
            + full_mass(0)*sum(dom%areas%elts(id  +1)%part(3:4)) &
            + full_mass(WEST)*dom%areas%elts(idW +1)%part(6))

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

       dom%qe%elts(EDGE*id  +UP+1) = 0.5_8*(pv_W + pv_UPLT)
       dom%qe%elts(EDGE*idSW+DG+1) = 0.5_8*(pv_SW_LORT + pv_SW_UPLT)
    end if

    if (c .eq. IJPLUS) then
       id  = idx(PATCH_SIZE,   PATCH_SIZE,   offs, dims)
       idN = idx(PATCH_SIZE,   PATCH_SIZE+1, offs, dims)
       idE = idx(PATCH_SIZE+1, PATCH_SIZE,   offs, dims)
       idS = idx(PATCH_SIZE,   PATCH_SIZE-1, offs, dims)
       idW = idx(PATCH_SIZE-1, PATCH_SIZE,   offs, dims)

       if (penalize) phi(0:WEST) = phi(0:WEST) + alpha_m1*penal%data(dom%id+1)%elts((/id,idN,idE,idS,idW/)+1)

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

       dom%qe%elts(EDGE*id+RT+1) = 0.5_8*(pv + pv_S)
       dom%qe%elts(EDGE*id+UP+1) = 0.5_8*(pv + pv_W)
    end if
  end subroutine post_step1

  subroutine step1(dom, p, zlev)
    type(Domain) dom
    integer p, zlev
    integer i, j, id, offs(0:N_BDRY), dims(2,N_BDRY)
    integer n, e, s, w, ne, sw
    real(8) u_prim_up, u_dual_up, u_prim_dg, u_dual_dg, u_prim_rt, u_dual_rt
    real(8) u_prim_dn, u_dual_dn, u_prim_sw, u_dual_sw, u_prim_lt, u_dual_lt
    real(8) pv_LORT, pv_UPLT, pv_S, pv_W, vort_W, vort_S, vort_LORT, vort_UPLT
    logical S_bdry, W_bdry
    real(8) phi(0:N_BDRY), full_mass(0:N_BDRY), full_temp(0:N_BDRY)

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

      phi(SOUTHWEST) = 1.0

      if (penalize) phi(SOUTHWEST) = phi(SOUTHWEST) + alpha_m1*penal%data(dom%id+1)%elts(id+sw+1)

      full_mass(SOUTHWEST) = mass(id+sw+1) + mean(S_MASS,zlev)
      full_temp(SOUTHWEST)  = temp(id+sw+1) + mean(S_TEMP,zlev)
      
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

      dom%qe%elts(EDGE*(id+ w)+RT+1) = 0.5_8*(pv_W   +pv_UPLT)
      dom%qe%elts(EDGE*(id+sw)+DG+1) = 0.5_8*(pv_LORT+pv_UPLT)
      dom%qe%elts(EDGE*(id+s )+UP+1) = 0.5_8*(pv_LORT+pv_S)

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
      h_flux(EDGE*id+UP+1) = u_dual_up*(full_scalar(0)         + full_scalar(NORTH))*0.5_8
      h_flux(EDGE*id+DG+1) = u_dual_dg*(full_scalar(NORTHEAST) + full_scalar(0))*0.5_8
      h_flux(EDGE*id+RT+1) = u_dual_rt*(full_scalar(0)         + full_scalar(EAST))*0.5_8
    end subroutine cal_flux2

    subroutine comput()
      ! Computes physical quantities during upward integration

      ! Find the velocity on primal and dual grid edges, which are equal except for the length of the
      ! side they are on
      u_prim_up = velo(EDGE*id+UP+1)*dom%len%elts(EDGE*id+UP+1)
      u_dual_up = velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
      u_prim_dg = velo(EDGE*id+DG+1)*dom%len%elts(EDGE*id+DG+1)
      u_dual_dg = velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
      u_prim_rt = velo(EDGE*id+RT+1)*dom%len%elts(EDGE*id+RT+1)
      u_dual_rt = velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)

      phi(0:NORTHEAST) = 1.0

      if (penalize) phi(0:NORTHEAST) = phi(0:NORTHEAST) + &
           alpha_m1*penal%data(dom%id+1)%elts(id+(/0,n,e,s,w,ne/)+1)

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
      u_dual_lt = velo(EDGE*(id+W)+RT+1)*dom%pedlen%elts(EDGE*(id+W)+RT+1)

      dom%kin_energy%elts(id+1) = &
           (u_prim_up*u_dual_up + u_prim_dg*u_dual_dg + u_prim_rt*u_dual_rt + &
           u_prim_dn*u_dual_dn + u_prim_sw*u_dual_sw + u_prim_lt*u_dual_lt &
           )* (1.0_8/4.0_8)*dom%areas%elts(id+1)%hex_inv

      if (phi(0) .ne. 0) phi(0) = 1.0_8/phi(0)

      ! Define the Bernoulli function
      if (compressible) then ! Compressible case
         dom%bernoulli%elts(id+1) = dom%kin_energy%elts(id+1) + dom%geopot%elts(id+1)
      else ! Incompressible case
         dom%bernoulli%elts(id+1) = dom%kin_energy%elts(id+1) + dom%geopot%elts(id+1) + dom%press%elts(id+1)
      end if

      if (viscosity .ne. 0) dom%divu%elts(id+1) = dom%areas%elts(id+1)%hex_inv * &
           (u_dual_up - u_dual_dg + u_dual_rt - u_dual_dn + u_dual_sw - u_dual_lt)

      dom%vort%elts(LORT+TRIAG*id+1) = - (u_prim_rt + u_prim_dg + velo(EDGE*(id+E)+UP+1)*dom%len%elts(EDGE*(id+E)+UP+1))
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

      dom%qe%elts(EDGE*id+RT+1) = 0.5_8*(pv_S    + pv_LORT)
      dom%qe%elts(EDGE*id+DG+1) = 0.5_8*(pv_UPLT + pv_LORT)
      dom%qe%elts(EDGE*id+UP+1) = 0.5_8*(pv_UPLT + pv_W)
    end subroutine comput
  end subroutine step1

  subroutine integrate_pressure_up(dom, i, j, zlev, offs, dims)
    !integrate pressure/Lagrange multiplier and pressure quantities upward at all nodes
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id

    id   = idx(i,     j,     offs, dims)

    if (mass(id+1)+mean(S_MASS,zlev) .lt. 1e-4_8) then
       print *, 'fatal error: a horizontal layer thickness is being squeezed to zero, namely, at zlev=', zlev
       write(6,*) mass(id+1), mean(S_MASS,zlev)
       stop
    end if

    if (compressible) then !compressible case
       !integrate the pressure from bottom zlev up to top zlev
       if (zlev .eq. 1) then !bottom zlev, integrate half of a layer further down to the surface
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*mass(id+1)
       else !other layers equal to half of previous layer and half of current layer
          dom%press%elts(id+1) = dom%press%elts(id+1) - 0.5_8*grav_accel*(mass(id+1)+dom%adj_mass%elts(id+1))
       end if

       if (zlev .eq. zlevels) then !top zlev, purely diagnostic
          if (abs((dom%press%elts(id+1) - 0.5_8*grav_accel*mass(id+1)) - press_infty) .gt. 1e-11_8) then
             PRINT *, 'warning: upward integration of pressure not resulting in zero at top interface'
             stop
          end if
       end if

       !compute exner pressure from the pressure
       dom%exner%elts(id+1) = c_p*(dom%press%elts(id+1)/ref_press)**kappa
       !compute the specific volume as kappa*theta*pi/p
       dom%spec_vol%elts(id+1) = kappa*temp(id+1)/mass(id+1)*dom%exner%elts(id+1)/dom%press%elts(id+1)

       !integrate the geopotential; surf_geopot is in shared.f90; (18) and below in DYNAMICO
       if (zlev .eq. 1) then !bottom zlev, integrate half of a layer up from the surface
          dom%geopot%elts(id+1) = dom%surf_geopot%elts(id+1) + 0.5_8*grav_accel*mass(id+1)*dom%spec_vol%elts(id+1)
       else !other layers equal to half of previous layer and half of current layer
          dom%geopot%elts(id+1) = dom%geopot%elts(id+1) + &
               0.5_8*grav_accel*(mass(id+1)*dom%spec_vol%elts(id+1) + dom%adj_mass%elts(id+1)*dom%adj_spec_vol%elts(id+1))
       end if
    else !incompressible case
       !incompressible case: integrate the Lagrange multiplier (saved as pressure) from bottom zlev up to top zlev
       if (zlev .eq. 1) then !bottom zlev, integrate half of a layer further down to the surface
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*temp(id+1)
       else !other layers equal to half of previous layer and half of current layer     
          dom%press%elts(id+1) = dom%press%elts(id+1) - 0.5_8*grav_accel*(dom%adj_temp%elts(id+1) + temp(id+1))
       end if

       if (zlev .eq. zlevels) then !top zlev, purely diagnostic
          if (abs(dom%press%elts(id+1)-0.5_8*grav_accel*temp(id+1) - press_infty).gt. 1e-11_8) then
             print *, 'warning: upward integration of Lagrange multiplier not resulting in zero at top interface'
             stop
          end if
       end if

       !compute the specific volume as 1 divided by the constant density
       dom%spec_vol%elts(id+1) = 1.0_8/cst_density

       !integrate the geopotential; surf_geopot is in shared.f90; (18) and below in DYNAMICO
       if (zlev .eq. 1) then !bottom zlev, integrate half of a layer up from the surface
          dom%geopot%elts(id+1) = dom%surf_geopot%elts(id+1) + 0.5_8*grav_accel*mass(id+1)*dom%spec_vol%elts(id+1)
       else !other layers equal to half of previous layer and half of current layer
          dom%geopot%elts(id+1) = dom%geopot%elts(id+1) + &
               0.5_8*grav_accel*(mass(id+1)*dom%spec_vol%elts(id+1) + dom%adj_mass%elts(id+1)*dom%adj_spec_vol%elts(id+1))
       end if

       !compute exner pressure from the geopotential
       dom%exner%elts(id+1) = -dom%geopot%elts(id+1)
    end if

    !quantities for vertical integration in next zlev
    dom%adj_mass%elts(id+1)     = mass(id+1)
    dom%adj_temp%elts(id+1)     = temp(id+1)
    dom%adj_spec_vol%elts(id+1) = dom%spec_vol%elts(id+1)
  end subroutine integrate_pressure_up

  subroutine integrate_pressure_down(dom, i, j, zlev, offs, dims)
    !pressure is computed here during downward integration from zlev=zlevels to zlev=1
    !INCOMPRESSIBLE CASE TESTED ONLY (JEMF)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id

    id   = idx(i,     j,     offs, dims)

    if (compressible) then !compressible case
       !integrate (or, rather, interpolate) the pressure from top zlev down to bottom zlev; press_infty is user-set
       if (zlev .eq. zlevels) then !top zlev, it is an exception
          dom%press%elts(id+1) = press_infty + 0.5_8*grav_accel*mass(id+1)
       else !other layers equal to half of previous layer and half of current layer
          dom%press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*(mass(id+1) + dom%adj_mass%elts(id+1))
       end if

       !surface pressure is set (even at t=0) from downward numerical integration
       if (zlev .eq. 1) then
          !PRINT *, 'theoretical surface pressure=', dom%surf_press%elts(id+1), &
          !    ', numerical surface pressure=', (dom%press%elts(id+1)+0.5_8*grav_accel*mass(id+1)), &
          !    ', error=', abs(dom%surf_press%elts(id+1)-(dom%pressure%elts(id+1)+0.5_8*grav_accel*mass(id+1)))
          dom%surf_press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*mass(id+1)
       end if
    else !incompressible case
       !integrate (or, rather, interpolate) the pressure from top zlev down to bottom zlev; press_infty is user-set
       if (zlev .eq. zlevels) then !top zlev, it is an exception
          dom%press%elts(id+1) = press_infty + 0.5_8*grav_accel*temp(id+1)
       else !other layers equal to half of previous layer and half of current layer
          dom%press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*(dom%adj_temp%elts(id+1)+ temp(id+1))
       end if

       !surface pressure is set (even at t=0) from downward numerical integration
       if (zlev .eq. 1) then
          !   PRINT *, 'theoretical surface pressure=', dom%surf_press%elts(id+1), &
          !       ', numerical surface pressure=', (dom%press%elts(id+1)+0.5_8*grav_accel*mass(id+1)), &
          !       ', error=', abs(dom%surf_press%elts(id+1)-(dom%press%elts(id+1)+0.5_8*grav_accel*mass(id+1)))
          dom%surf_press%elts(id+1) = dom%press%elts(id+1) + 0.5_8*grav_accel*temp(id+1)
          !PRINT *, 'surf_press', dom%surf_press%elts(id+1)
       end if
    end if

    !quantities for vertical integration in next zlev
    dom%adj_mass%elts(id+1) = mass(id+1)
    dom%adj_temp%elts(id+1) = temp(id+1)
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

    dvelo(EDGE*id+RT+1) = dom%qe%elts(EDGE*id+RT+1)*( &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)&
         %elts((/ EDGE*id+DG, EDGE*id+UP, EDGE*idW+RT, EDGE*idSW+DG, EDGE*idS+UP /)+1) * wgt1) + &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)&
         %elts((/ EDGE*idS+DG, EDGE*idSE+UP, EDGE*idE+RT, EDGE*idE+DG, EDGE*idE+UP /)+1) * wgt2))

    wgt1 = get_weights(dom, id, 1)
    wgt2 = get_weights(dom, idNE, 4)

    dvelo(EDGE*id+DG+1) = dom%qe%elts(EDGE*id+DG+1)*( &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*id+UP, EDGE*idW+RT, EDGE*idSW+DG, EDGE*idS+UP, EDGE*id+RT/)+1) * wgt1) + &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*idE+UP, EDGE*idNE+RT, DG+EDGE*idNE, EDGE*idNE+UP, EDGE*idN+RT/)+1) * wgt2))

    wgt1 = get_weights(dom, id, 2)
    wgt2 = get_weights(dom, idN, 5)

    dvelo(EDGE*id+UP+1) = dom%qe%elts(EDGE*id+UP+1)*( &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*idW+RT, EDGE*idSW+DG, EDGE*idS+UP, EDGE*id+RT, EDGE*id+DG /)+1) * wgt1) + &
         sum(horiz_flux(S_MASS,zlev)%data(dom%id+1)%&
         elts((/ EDGE*idN+RT, EDGE*idN+DG, EDGE*idN+UP, EDGE*idNW+RT, EDGE*idW+DG /)+1) * wgt2))
  end subroutine du_Qperp_Enstrophy

  subroutine du_source(dom, i, j, zlev, offs, dims)
    !add additional effects (viscosity, bottom friction, wind stress, coastal boundaries) to dvelo
    ![Aechtner thesis page 56, Kevlahan, Dubos and Aechtner (2015)]
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, idE, idN, idNE
    real(8), dimension(RT:NODE) :: chi1, phi, hk_phi, kin_energ, u_mag

    id = idx(i, j, offs, dims)

    call du_Qperp(dom, i, j, zlev, offs, dims)

    if (viscosity .ne. 0) call diff_mom(dom, i, j, offs, dims)

    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if (penalize) then
       chi1 = penal%data(dom%id+1)%elts((/idN,idNE,idE,id/)+1)
       phi = 1 + alpha_m1*chi1
       chi1 = -ieta*chi1
       chi1(RT:UP) = 0.5_8*(chi1(NODE)+chi1(RT:UP)) ! interpolate p->u

       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = dvelo(EDGE*id+RT+1:EDGE*id+UP+1) &
            + velo(EDGE*id+RT+1:EDGE*id+UP+1)*dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1)*chi1(RT:UP)
    else
       phi = 1.0_8
    end if

    if (wind_stress .or. bottom_friction) then

       hk_phi = mass((/idN,idNE,idE,id/)+1)

       where (hk_phi .ne. 0) hk_phi = phi/hk_phi

       hk_phi(RT:UP) = 0.5_8*(hk_phi(NODE)+hk_phi(RT:UP)) ! interpolate p->u

       if (wind_stress) then
          dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = dvelo(EDGE*id+RT+1:EDGE*id+UP+1) &
               + dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) &
               *dom%windstress%elts(EDGE*id+RT+1:EDGE*id+UP+1)*hk_phi(RT:UP)
       end if

       if (bottom_friction) then
          kin_energ = dom%kin_energy%elts((/idN,idNE,idE,id/)+1)
          u_mag(RT:UP) = sqrt(2.0_8*0.5_8*(kin_energ(NODE)+kin_energ(RT:UP))) ! interpolate p->u
          dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = dvelo(EDGE*id+RT+1:EDGE*id+UP+1) &
               - dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) &
               *friction_coeff*velo(EDGE*id+RT+1:EDGE*id+UP+1)*u_mag(RT:UP)*hk_phi(RT:UP)
       end if
    end if
  end subroutine du_source

  subroutine du_Qperp(dom, i, j, zlev, offs, dims)
    !compute energy-conserving Qperp and add it to dvelo [Aechtner thesis page 44]
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

    dvelo(EDGE*id+RT+1) = &
         h_mflux(EDGE*id+DG+1)*0.5_8*(dom%qe%elts(EDGE*id+DG+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt1(1) + &
         h_mflux(EDGE*id+UP+1)*0.5_8*(dom%qe%elts(EDGE*id+UP+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt1(2) + &
         h_mflux(EDGE*idW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idW+RT+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt1(3) + &
         h_mflux(EDGE*idSW+DG+1)*0.5_8*(dom%qe%elts(EDGE*idSW+DG+1) &
         + dom%qe%elts(EDGE*id+RT+1))*wgt1(4) + &
         h_mflux(EDGE*idS+UP+1)*0.5_8*(dom%qe%elts(EDGE*idS+UP+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt1(5) + &
         h_mflux(DG+EDGE*idS+1)*0.5_8*(dom%qe%elts(DG+EDGE*idS+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt2(1) + &
         h_mflux(EDGE*idSE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idSE+UP+1) &
         + dom%qe%elts(EDGE*id+RT+1))*wgt2(2) + &
         h_mflux(EDGE*idE+RT+1)*0.5_8*(dom%qe%elts(EDGE*idE+RT+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt2(3) + &
         h_mflux(DG+EDGE*idE+1)*0.5_8*(dom%qe%elts(DG+EDGE*idE+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt2(4) + &
         h_mflux(EDGE*idE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idE+UP+1) + &
         dom%qe%elts(EDGE*id+RT+1))*wgt2(5)

    wgt1 = get_weights(dom, id, 1)
    wgt2 = get_weights(dom, idNE, 4)

    dvelo(EDGE*id+DG+1) = &
         h_mflux(EDGE*id+UP+1)*0.5_8*(dom%qe%elts(EDGE*id+UP+1) + &
         dom%qe%elts(EDGE*id+DG+1))*wgt1(1) + &
         h_mflux(EDGE*idW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idW+RT+1) + &
         dom%qe%elts(EDGE*id+DG+1))*wgt1(2) + &
         h_mflux(EDGE*idSW+DG+1)*0.5_8*(dom%qe%elts(EDGE*idSW+DG+1) &
         + dom%qe%elts(EDGE*id+DG+1))*wgt1(3) + &
         h_mflux(EDGE*idS+UP+1)*0.5_8*(dom%qe%elts(EDGE*idS+UP+1) + &
         dom%qe%elts(EDGE*id+DG+1))*wgt1(4) + &
         h_mflux(EDGE*id+RT+1)*0.5_8*(dom%qe%elts(EDGE*id+RT+1) + &
         dom%qe%elts(EDGE*id+DG+1))*wgt1(5) + &
         h_mflux(EDGE*idE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idE+UP+1) + &
         dom%qe%elts(EDGE*id+DG+1))*wgt2(1) + &
         h_mflux(EDGE*idNE+RT+1)*0.5_8*(dom%qe%elts(EDGE*idNE+RT+1) &
         + dom%qe%elts(EDGE*id+DG+1))*wgt2(2) + &
         h_mflux(EDGE*idNE+DG+1)*0.5_8*(dom%qe%elts(EDGE*idNE+DG+1) &
         + dom%qe%elts(EDGE*id+DG+1))*wgt2(3) + &
         h_mflux(EDGE*idNE+UP+1)*0.5_8*(dom%qe%elts(EDGE*idNE+UP+1) &
         + dom%qe%elts(EDGE*id+DG+1))*wgt2(4) + &
         h_mflux(EDGE*idN+RT+1)*0.5_8*(dom%qe%elts(EDGE*idN+RT+1) + &
         dom%qe%elts(EDGE*id+DG+1))*wgt2(5)

    wgt1 = get_weights(dom, id, 2)
    wgt2 = get_weights(dom, idN, 5)

    dvelo(EDGE*id+UP+1) = &
         h_mflux(EDGE*idW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idW+RT+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt1(1) + &
         h_mflux(EDGE*idSW+DG+1)*0.5_8*(dom%qe%elts(EDGE*idSW+DG+1) &
         + dom%qe%elts(EDGE*id+UP+1))*wgt1(2) + &
         h_mflux(EDGE*idS+UP+1)*0.5_8*(dom%qe%elts(EDGE*idS+UP+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt1(3) + &
         h_mflux(EDGE*id+RT+1)*0.5_8*(dom%qe%elts(EDGE*id+RT+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt1(4) + &
         h_mflux(EDGE*id+DG+1)*0.5_8*(dom%qe%elts(EDGE*id+DG+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt1(5) + &
         h_mflux(EDGE*idN+RT+1)*0.5_8*(dom%qe%elts(EDGE*idN+RT+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt2(1) + &
         h_mflux(EDGE*idN+DG+1)*0.5_8*(dom%qe%elts(EDGE*idN+DG+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt2(2) + &
         h_mflux(EDGE*idN+UP+1)*0.5_8*(dom%qe%elts(EDGE*idN+UP+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt2(3) + &
         h_mflux(EDGE*idNW+RT+1)*0.5_8*(dom%qe%elts(EDGE*idNW+RT+1) &
         + dom%qe%elts(EDGE*id+UP+1))*wgt2(4) + &
         h_mflux(EDGE*idW+DG+1)*0.5_8*(dom%qe%elts(EDGE*idW+DG+1) + &
         dom%qe%elts(EDGE*id+UP+1))*wgt2(5)
  end subroutine du_Qperp

  subroutine masstemp_trend(dom, i, j, zlev, offs, dims)
    type(Domain) :: dom
    integer :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id

    id   = idx(i, j, offs, dims)
    
    dmass(id+1) = scalar_trend(h_mflux, dom, i, j, offs, dims, id)
    dtemp(id+1) = scalar_trend(h_tflux, dom, i, j, offs, dims, id)
  end subroutine masstemp_trend

  function scalar_trend(h_flux, dom, i, j, offs, dims, id)
    real(8) :: scalar_trend
    real(8), dimension(:), pointer :: h_flux
    type(Domain) :: dom
    integer :: i, j, id
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer ::  idS, idW, idSW
   
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    
    scalar_trend = -(h_flux(EDGE*id+UP+1)  - h_flux(EDGE*id+DG+1)   + h_flux(EDGE*id+RT+1) - &
         h_flux(EDGE*idS+UP+1) + h_flux(EDGE*idSW+DG+1) - h_flux(EDGE*idW+RT+1)) &
         *dom%areas%elts(id+1)%hex_inv
  end function scalar_trend

  subroutine du_gradB(dom, i, j, zlev, offs, dims)
    !add gradient of the Bernoulli function to dvelo [Aechtner thesis page 58]
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idE
    integer idN
    integer idNE

    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    dvelo(EDGE*id+RT+1) = (dvelo(EDGE*id+RT+1) - (dom%bernoulli%elts(idE+1) - dom%bernoulli%elts(id+1)))/&
         dom%len%elts(EDGE*id+RT+1)

    dvelo(EDGE*id+DG+1) = (dvelo(EDGE*id+DG+1) - (dom%bernoulli%elts(id+1) - dom%bernoulli%elts(idNE+1)))/&
         dom%len%elts(EDGE*id+DG+1)

    dvelo(EDGE*id+UP+1) = (dvelo(EDGE*id+UP+1) - (dom%bernoulli%elts(idN+1) - dom%bernoulli%elts(id+1)))/&
         dom%len%elts(EDGE*id+UP+1)
  end subroutine du_gradB

  subroutine du_gradB_gradExn(dom, i, j, zlev, offs, dims)
    !add gradient of the Bernoulli and Exner to dvelo [DYNAMICO (23)-(25)]
    !mass and potential temperature trend is zero
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idE
    integer idN
    integer idNE
    real(8) full_pot_temp(0:N_BDRY)

    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    !see DYNAMICO between (23)-(25), geopotential still known from step1_upw
    !the theta multiplying the exner gradient is the edge-averaged non-mass-weighted potential temperature

    full_pot_temp(0)         = (temp(id+1)   + mean(S_TEMP,zlev))/(mass(id+1)   + mean(S_MASS,zlev))
    full_pot_temp(NORTH)     = (temp(idN+1)  + mean(S_TEMP,zlev))/(mass(idN+1)  + mean(S_MASS,zlev))
    full_pot_temp(EAST)      = (temp(idE+1)  + mean(S_TEMP,zlev))/(mass(idE+1)  + mean(S_MASS,zlev))
    full_pot_temp(NORTHEAST) = (temp(idNE+1) + mean(S_TEMP,zlev))/(mass(idNE+1) + mean(S_MASS,zlev))
    
    if (compressible) then
       dvelo(EDGE*id+RT+1) = (dvelo(EDGE*id+RT+1) - &
            (dom%bernoulli%elts(idE+1) - dom%bernoulli%elts(id+1)) &
            - 0.5_8*(full_pot_temp(0)+full_pot_temp(EAST))* &
            (dom%exner%elts(idE+1) - dom%exner%elts(id+1)))/dom%len%elts(EDGE*id+RT+1)

       dvelo(EDGE*id+DG+1) = (dvelo(EDGE*id+DG+1) - &
            (dom%bernoulli%elts(id+1) - dom%bernoulli%elts(idNE+1)) &
            - 0.5_8*(full_pot_temp(0)+full_pot_temp(NORTHEAST))* &
            (dom%exner%elts(id+1) - dom%exner%elts(idNE+1)))/dom%len%elts(EDGE*id+DG+1)

       dvelo(EDGE*id+UP+1) = (dvelo(EDGE*id+UP+1) - &
            (dom%bernoulli%elts(idN+1) - dom%bernoulli%elts(id+1)) &
            - 0.5_8*(full_pot_temp(0)+full_pot_temp(NORTH)) * &
            (dom%exner%elts(idN+1) - dom%exner%elts(id+1)))/dom%len%elts(EDGE*id+UP+1)
    else !incompressible case
       dvelo(EDGE*id+RT+1) = (dvelo(EDGE*id+RT+1) - &
            (dom%bernoulli%elts(idE+1) - dom%bernoulli%elts(id+1)) &
            - 0.5_8*(2.0_8-full_pot_temp(0)-full_pot_temp(EAST))* &
            (dom%exner%elts(idE+1) - dom%exner%elts(id+1)))/dom%len%elts(EDGE*id+RT+1)

       dvelo(EDGE*id+DG+1) = (dvelo(EDGE*id+DG+1) - &
            (dom%bernoulli%elts(id+1) - dom%bernoulli%elts(idNE+1)) &
            - 0.5_8*(2.0_8-full_pot_temp(0)-full_pot_temp(NORTHEAST))* &
            (dom%exner%elts(id+1) - dom%exner%elts(idNE+1)))/dom%len%elts(EDGE*id+DG+1)

       dvelo(EDGE*id+UP+1) = (dvelo(EDGE*id+UP+1) - &
            (dom%bernoulli%elts(idN+1) - dom%bernoulli%elts(id+1)) &
            - 0.5_8*(2.0_8-full_pot_temp(0)-full_pot_temp(NORTH)) * &
            (dom%exner%elts(idN+1) - dom%exner%elts(id+1)))/dom%len%elts(EDGE*id+UP+1)
    end if
  end subroutine du_gradB_gradExn

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
