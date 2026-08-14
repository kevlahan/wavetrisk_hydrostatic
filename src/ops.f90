module ops_mod

  ! Routines used to compute multiscale trend

  use kind_mod,   only : dp
  use shared_mod, only : EDGE, N_BDRY, N_VARIABLE, EAST, NORTH, NORTHEAST, NORTHWEST, SOUTH, SOUTHEAST, SOUTHWEST, WEST, &
       S_MASS, S_VELO, IJMINUS, IJPLUS, IMINUSJPLUS, IPLUSJMINUS, &
       compressible, zlevels, zmax, c_p, eps, grav_accel, &
       kappa, p_0, p_top, ref_density, radius, RT, DG, UP, S_TEMP, TRIAG, LORT, UPLT, TRSK, mode_split, scalars, z_null

  use diagnostics_mod, only : cal_surf_press, cal_vort, div, gradi_e, integrate_pressure_up
  use domain_ops_mod,  only : apply, apply_onescale_to_patch
  use patch_mod,       only : LAST, PATCH_SIZE
  use utils_mod,       only : interp, phi_node, porous_density
  use init_mod,        only : physics_scalar_flux, physics_velo_source, surf_geopot
  
  use domain_mod, only : bernoulli, Domain, Float_Field, grid, sides_dims, exner, exner_fun, h_flux, h_mflux, horiz_flux, ke, &
       mass, mean_m, temp, mean_t, dscalar, dvelo, scalar, sol, sol_mean, topography, Laplacian, qe, velo, vort, idx, id_edge

  implicit none

  private
  public :: post_step1, step1, scalar_trend, du_grad, du_source, Qperp
  public :: comp_offs3

  
contains
  
  
  subroutine step1 (dq, q, dom, p, zlev, itype)
    ! itype = 0 is standard computation of all quantities
    ! itype = 1 computes only scalar flux of pointer scalar
    
    implicit none
    
    integer,           intent(in)              :: itype, p
    integer,           intent(in), optional    :: zlev
    type(Domain),      intent(inout)           :: dom
    type(Float_Field), intent(inout), optional :: dq(1:N_VARIABLE,1:zmax), q(1:N_VARIABLE,1:zmax)

    integer :: j, id,  n, e, s, w, ne, sw, v
    integer :: offs(0:N_BDRY) 
    integer :: dims(2,N_BDRY)

    real(dp) :: u_prim_UP, u_dual_UP, u_prim_DG, u_prim_DG_S, u_prim_DG_W, u_dual_DG, u_prim_RT, u_dual_RT
    real(dp) :: u_prim_UP_S, u_dual_UP_S, u_prim_DG_SW, u_dual_DG_SW, u_prim_RT_W, u_dual_RT_W
    real(dp) :: circ_S_UPLT, circ_W_LORT, pv_LORT, pv_UPLT, pv_S_UPLT, pv_SW_LORT, pv_SW_UPLT, pv_W_LORT

    real(dp), dimension(0:N_BDRY,scalars(1):scalars(2)) :: rho_dz
    real(dp), dimension(1:EDGE)                         :: physics_flux

    logical :: S_bdry, W_bdry

    call comp_offs3 (dom, p, offs, dims)

    S_bdry = (dom%patch%elts(p+1)%neigh(SOUTH) < 0)
    if (S_bdry) S_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(SOUTH)+1)%side > 0)

    W_bdry = (dom%patch%elts(p+1)%neigh(WEST)  < 0)
    if (W_bdry) W_bdry = (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(WEST)+1)%side > 0)

    id = offs(0)

    n = PATCH_SIZE       ; ne = n+1
    w = offs(WEST)       ; e =  +1
    sw = offs(SOUTHWEST) ; s = offs(SOUTH)

    call comput
    if (W_bdry .or. S_bdry) call comp_SW

    w = -1; sw = s+w
    do id = offs(0)+1, offs(0)+LAST-1
       call comput
       if (S_bdry) call comp_SW
    end do

    e = offs(EAST); ne = dims(1,EAST) + e
    id = offs(0)+LAST
    call comput
    if (S_bdry) call comp_SW

    s = -PATCH_SIZE
    do j = 2, PATCH_SIZE-1
       id = offs(0) + PATCH_SIZE*(j-1)
       e = +1; ne = n+e
       w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*(j-1) ! Correct for dimension smaller than patch if boundary
       sw = w-dims(1,WEST)
       call comput
       if (W_bdry) call comp_SW

       w = -1; sw = s+w; ne = n+e
       do id = offs(0)+PATCH_SIZE*(j-1)+1, offs(0)+PATCH_SIZE*(j-1)+LAST-1
          call comput
       end do

       id = offs(0)+PATCH_SIZE*j-1
       e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*(j-1); ne = e+dims(1,EAST)
       call comput
    end do

    id = offs(0)+PATCH_SIZE*LAST
    n = offs(NORTH)
    w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST
    sw = w-dims(1,WEST)
    e = +1; ne = n+e
    call comput
    if (W_bdry) call comp_SW

    w = -1; sw = s+w
    do id = offs(0)+PATCH_SIZE*LAST+1, offs(0)+PATCH_SIZE*PATCH_SIZE-2
       call comput
    end do

    id = offs(0)+PATCH_SIZE*PATCH_SIZE-1
    ne = offs(NORTHEAST)
    e = offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
    call comput

    if (dom%patch%elts(p+1)%neigh(NORTH) < 0) then ! Neighbour is boundary
       if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(NORTH)+1)%side > 0) then ! Domain boundary
          id = offs(0)+PATCH_SIZE*LAST + offs(NORTH) ! id + n
          s =                                   - offs(NORTH) ! relative to current id
          w = offs(NORTHWEST)                   - offs(NORTH) 
          sw = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*LAST - offs(NORTH)
          e = 1
          n = dims(1,NORTH)
          ne = n+e
          call comput
          if (W_bdry) call comp_SW

          w = -1
          sw = -offs(NORTH) + w
          do id = offs(0)+PATCH_SIZE*LAST+1+offs(NORTH), offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-2
             call comput
          end do

          e = offs(NORTHEAST)              - offs(NORTH) 
          ne = e+dims(1,NORTHEAST)
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTH)-1
          call comput

          ! NORTHEAST point
          id = offs(0)+PATCH_SIZE*PATCH_SIZE+offs(NORTHEAST)-1
          s = -offs(NORTHEAST) + offs(EAST) + (dims(1,EAST)-PATCH_SIZE)*LAST
          sw= -offs(NORTHEAST)
          w = -offs(NORTHEAST) + offs(NORTH)
          n = dims(1,NORTHEAST)
          e = 1
          ne = n+e
          call comput
       end if
    end if

    if (dom%patch%elts(p+1)%neigh(EAST) < 0) then ! neighbour is boundary
       if (dom%bdry_patch%elts(-dom%patch%elts(p+1)%neigh(EAST)+1)%side <= 0) return ! not domain bdry
       id = offs(0) + LAST + offs(EAST) ! id + n
       w  =                - offs(EAST) ! relative to current id
       s  = offs(SOUTHEAST)- offs(EAST) 
       sw = offs(SOUTH)    - offs(EAST)
       n  = dims(1,EAST)
       e  = 1
       ne = n+e
       call comput
       if (S_bdry) call comp_SW

       s = -dims(1,EAST)
       do j = 1, LAST-1
          id = offs(0) + LAST + offs(EAST) + j*dims(1,EAST)
          w  =                - offs(EAST) + j*(PATCH_SIZE-dims(1,EAST))
          sw = w-PATCH_SIZE
          call comput
       end do

       id = offs(0) + LAST + offs(EAST) + LAST*dims(1,EAST)
       n  =   offs(NORTHEAST) + LAST*(PATCH_SIZE-dims(1,EAST)) - offs(EAST)
       w  = - offs(EAST)      + LAST*(PATCH_SIZE-dims(1,EAST))
       sw = w-PATCH_SIZE
       ne = n+e
       call comput
    end if
    
  contains
    
    subroutine comput 
      ! Computes physical quantities during upward integration
      
      implicit none
      
      integer  :: idE, idN, idNE, idS, idSW, idW
      integer  :: d, id_i, idE_i, idN_i, idNE_i, idS_i, idW_i, k
      integer  :: id_rhodz(0:NORTHEAST)
      real(dp) :: circ_LORT, circ_UPLT, Phi_k, dz, dz0
      real(dp) :: u_prim_UP_E, u_prim_RT_N
      real(dp) :: csq(0:EDGE) , phi(0:EDGE) 
      real(dp) :: pv_mass(4)

      d = dom%id + 1

      idN  = id+N
      idE  = id+E
      idS  = id+S
      idW  = id+W
      idNE = id+NE
      idSW = id+SW

      id_i   = id+1
      idE_i  = idE+1
      idN_i  = idN+1
      idNE_i = idNE+1
      idS_i  = idS+1
      idW_i  = idW+1

      if (mode_split) then
         phi(0) = phi_node (d, id_i,   zlevels)
         phi(1) = phi_node (d, idE_i,  zlevels)
         phi(2) = phi_node (d, idNE_i, zlevels)
         phi(3) = phi_node (d, idN_i,  zlevels)
      else
         phi = 1.0_dp
      end if

      if (itype == 1) then ! scalar gradient flux         
         h_flux(EDGE*id+RT+1) = (scalar(idE_i) - scalar(id_i))   / dom%len%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
         h_flux(EDGE*id+DG+1) = (scalar(id_i)  - scalar(idNE_i)) / dom%len%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
         h_flux(EDGE*id+UP+1) = (scalar(idN_i) - scalar(id_i))   / dom%len%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)
      elseif (itype == 2) then ! external pressure gradient flux, (H + eta^n) grad(eta^(n+1)) * edge_length, in elliptic operator
         csq(0) = abs (topography%data(d)%elts(id_i))   * phi(0) + mass(id_i)
         csq(1) = abs (topography%data(d)%elts(idE_i))  * phi(1) + mass(idE_i)
         csq(2) = abs (topography%data(d)%elts(idNE_i)) * phi(2) + mass(idNE_i)
         csq(3) = abs (topography%data(d)%elts(idN_i))  * phi(3) + mass(idN_i)
         csq = grav_accel * csq
         
         h_flux(EDGE*id+RT+1) = interp (csq(0), csq(1)) * &
              (scalar(idE_i)/phi(1) - scalar(id_i)/phi(0))   / dom%len%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1) 

         h_flux(EDGE*id+DG+1) = interp (csq(0), csq(2)) * &
              (scalar(id_i)/phi(0)  - scalar(idNE_i)/phi(2)) / dom%len%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1) 

         h_flux(EDGE*id+UP+1) = interp (csq(0), csq(3)) * &
              (scalar(idN_i)/phi(3) - scalar(id_i)/phi(0))   / dom%len%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)
      elseif (itype == 3) then ! external pressure gradient
         h_flux(EDGE*id+RT+1) = grav_accel * (scalar(idE_i)/phi(1) - scalar(id_i)/phi(0))   / dom%len%elts(EDGE*id+RT+1) 
         h_flux(EDGE*id+DG+1) = grav_accel * (scalar(id_i)/phi(0)  - scalar(idNE_i)/phi(2)) / dom%len%elts(EDGE*id+DG+1) 
         h_flux(EDGE*id+UP+1) = grav_accel * (scalar(idN_i)/phi(3) - scalar(id_i)/phi(0))   / dom%len%elts(EDGE*id+UP+1)
      elseif (itype == 4) then ! sum vertical flux
         h_flux(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_dp
         do k = 1, zlevels
            dz0 = q(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)

            dz = interp (dz0, q(S_MASS,k)%data(d)%elts(idE_i) + sol_mean(S_MASS,k)%data(d)%elts(idE_i))
            h_flux(EDGE*id+RT+1) = h_flux(EDGE*id+RT+1) + q(S_VELO,k)%data(d)%elts(EDGE*id+RT+1) * dz

            dz = interp (dz0, q(S_MASS,k)%data(d)%elts(idNE_i) + sol_mean(S_MASS,k)%data(d)%elts(idNE_i))
            h_flux(EDGE*id+DG+1) = h_flux(EDGE*id+DG+1) + q(S_VELO,k)%data(d)%elts(EDGE*id+DG+1) * dz

            dz = interp (dz0, q(S_MASS,k)%data(d)%elts(idN_i) + sol_mean(S_MASS,k)%data(d)%elts(idN_i))
            h_flux(EDGE*id+UP+1) = h_flux(EDGE*id+UP+1) + q(S_VELO,k)%data(d)%elts(EDGE*id+UP+1) * dz
         end do
         h_flux(EDGE*id+RT+1:EDGE*id+UP+1) = h_flux(EDGE*id+RT+1:EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+RT+1:EDGE*id+UP+1) 
      elseif (itype == 5) then ! divu
         h_flux(EDGE*id+RT+1) = velo(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
         h_flux(EDGE*id+DG+1) = velo(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
         h_flux(EDGE*id+UP+1) = velo(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)
      elseif (itype == 6) then ! surface pressure gradient flux (for vertical velocity)
         h_flux(EDGE*id+RT+1) = (scalar(idE_i) - scalar(id_i))   * velo(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1) / 4
         h_flux(EDGE*id+DG+1) = (scalar(id_i)  - scalar(idNE_i)) * velo(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1) / 4
         h_flux(EDGE*id+UP+1) = (scalar(idN_i) - scalar(id_i))   * velo(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1) / 4
      elseif (itype == 7) then ! mass flux (for vertical velocity)
         dz0 = mass(id+1) + mean_m(id+1)

         dz = interp (dz0, mass(idE+1) + mean_m(idE+1))
         h_flux(EDGE*id+RT+1) = velo(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1) * dz

         dz = interp (dz0, mass(idNE+1) + mean_m(idNE+1)) 
         h_flux(EDGE*id+DG+1) = velo(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1) * dz

         dz = interp (dz0, mass(idN+1) + mean_m(idN+1)) 
         h_flux(EDGE*id+UP+1) = velo(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1) * dz

      elseif (itype == 8) then ! vorticity
         u_prim_RT    = velo(EDGE*id  +RT+1) * dom%len%elts(EDGE*id  +RT+1)
         u_prim_RT_N  = velo(EDGE*idN +RT+1) * dom%len%elts(EDGE*idN +RT+1)
         u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
         u_prim_DG    = velo(EDGE*id  +DG+1) * dom%len%elts(EDGE*id  +DG+1)
         u_prim_DG_S  = velo(EDGE*idS +DG+1) * dom%len%elts(EDGE*idS +DG+1)
         u_prim_DG_W  = velo(EDGE*idW +DG+1) * dom%len%elts(EDGE*idW +DG+1)
         u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id  +UP+1)
         u_prim_UP_E  = velo(EDGE*idE +UP+1) * dom%len%elts(EDGE*idE +UP+1)
         u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)

         circ_LORT   =   u_prim_RT   + u_prim_UP_E + u_prim_DG 
         circ_UPLT   = -(u_prim_DG   + u_prim_UP   + u_prim_RT_N)
         circ_W_LORT =   u_prim_RT_W + u_prim_UP   + u_prim_DG_W
         circ_S_UPLT = -(u_prim_RT   + u_prim_DG_S + u_prim_UP_S)

         vort(TRIAG*id +LORT+1) = circ_LORT   / dom%triarea%elts(TRIAG*id +LORT+1) 
         vort(TRIAG*id +UPLT+1) = circ_UPLT   / dom%triarea%elts(TRIAG*id +UPLT+1)
         vort(TRIAG*idW+LORT+1) = circ_W_LORT / dom%triarea%elts(TRIAG*idW+LORT+1)
         vort(TRIAG*idS+UPLT+1) = circ_S_UPLT / dom%triarea%elts(TRIAG*idS+UPLT+1)
      elseif (itype == 0) then ! standard
         id_rhodz = [ id, idN, idE, idS, idW, idNE ]
         do v = scalars(1), scalars(2)
            rho_dz(0:NORTHEAST,v) = q(v,zlev)%data(d)%elts(id_rhodz+1) + sol_mean(v,zlev)%data(d)%elts(id_rhodz+1)
         end do

         u_dual_RT    = velo(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
         u_dual_DG    = velo(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
         u_dual_UP    = velo(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)
         u_dual_RT_W  = velo(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
         u_dual_DG_SW = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
         u_dual_UP_S  = velo(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

         u_prim_RT    = velo(EDGE*id  +RT+1) * dom%len%elts(EDGE*id  +RT+1)
         u_prim_RT_N  = velo(EDGE*idN +RT+1) * dom%len%elts(EDGE*idN +RT+1)
         u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
         u_prim_DG    = velo(EDGE*id  +DG+1) * dom%len%elts(EDGE*id  +DG+1)
         u_prim_DG_S  = velo(EDGE*idS +DG+1) * dom%len%elts(EDGE*idS +DG+1)
         u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
         u_prim_DG_W  = velo(EDGE*idW +DG+1) * dom%len%elts(EDGE*idW +DG+1)
         u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id  +UP+1)
         u_prim_UP_E  = velo(EDGE*idE +UP+1) * dom%len%elts(EDGE*idE +UP+1)
         u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)

         ! Potential vorticity
         circ_LORT   =   u_prim_RT   + u_prim_UP_E + u_prim_DG 
         circ_UPLT   = -(u_prim_DG   + u_prim_UP   + u_prim_RT_N)
         circ_W_LORT =   u_prim_RT_W + u_prim_UP   + u_prim_DG_W
         circ_S_UPLT = -(u_prim_RT   + u_prim_DG_S + u_prim_UP_S)

         pv_mass(1) = &
              rho_dz(0,        S_MASS) * dom%areas%elts(id_i  )%part(1) + &
              rho_dz(EAST,     S_MASS) * dom%areas%elts(idE_i )%part(3) + &
              rho_dz(NORTHEAST,S_MASS) * dom%areas%elts(idNE_i)%part(5)

         pv_mass(2) = &
              rho_dz(0,        S_MASS) * dom%areas%elts(id_i  )%part(2) + &
              rho_dz(NORTHEAST,S_MASS) * dom%areas%elts(idNE_i)%part(4) + &
              rho_dz(NORTH,    S_MASS) * dom%areas%elts(idN_i )%part(6)

         pv_mass(3) = &
              rho_dz(WEST, S_MASS) * dom%areas%elts(idW_i)%part(1) + &
              rho_dz(0,    S_MASS) * dom%areas%elts(id_i )%part(3) + &
              rho_dz(NORTH,S_MASS) * dom%areas%elts(idN_i)%part(5)

         pv_mass(4) = &
              rho_dz(SOUTH,S_MASS)  * dom%areas%elts(idS_i)%part(2) + &
              rho_dz(EAST, S_MASS)  * dom%areas%elts(idE_i)%part(4) + &
              rho_dz(0,    S_MASS)  * dom%areas%elts(id_i )%part(6)

         pv_LORT   = (dom%coriolis%elts(TRIAG*id +LORT+1) + circ_LORT  ) / pv_mass(1)
         pv_UPLT   = (dom%coriolis%elts(TRIAG*id +UPLT+1) + circ_UPLT  ) / pv_mass(2)
         pv_W_LORT = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_W_LORT) / pv_mass(3)
         pv_S_UPLT = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_S_UPLT) / pv_mass(4)

         qe(EDGE*id+RT+1) = interp (pv_S_UPLT, pv_LORT)
         qe(EDGE*id+DG+1) = interp (pv_UPLT,   pv_LORT)
         qe(EDGE*id+UP+1) = interp (pv_UPLT,   pv_W_LORT)

         ! Vorticity (for rotu velocity diffusion)
         vort(TRIAG*id +LORT+1) = circ_LORT   / dom%triarea%elts(TRIAG*id+LORT+1) 
         vort(TRIAG*id +UPLT+1) = circ_UPLT   / dom%triarea%elts(TRIAG*id+UPLT+1)
         vort(TRIAG*idW+LORT+1) = circ_W_LORT / dom%triarea%elts(TRIAG*idW+LORT+1)
         vort(TRIAG*idS+UPLT+1) = circ_S_UPLT / dom%triarea%elts(TRIAG*idS+UPLT+1)

         ! Kinetic energy (TRiSK formula) 
         ke(id_i) = (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT +  &
                     u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
              * dom%areas%elts(id_i)%hex_inv/4

         ! Interpolate geopotential from interfaces to level
         Phi_k = interp (dom%geopot%elts(id_i), dom%geopot_lower%elts(id_i))

         ! Bernoulli function
         if (compressible) then 
            bernoulli(id_i) = ke(id_i) + Phi_k
         else
            bernoulli(id_i) = ke(id_i) + Phi_k + dom%press%elts(id_i) / (ref_density * phi_node (d, id_i, zlev))
         end if

         ! Exner function in incompressible case from geopotential
         if (.not. compressible) exner(id_i) = -Phi_k

         do v = scalars(1), scalars(2)
            physics_flux = physics_scalar_flux (q(:,1:zlevels), dom, id, idE, idNE, idN, v, zlev)

            horiz_flux(v)%data(d)%elts(EDGE*id+RT+1) = u_dual_RT * interp (rho_dz(0,v), rho_dz(EAST,     v)) + physics_flux(RT+1)
            horiz_flux(v)%data(d)%elts(EDGE*id+DG+1) = u_dual_DG * interp (rho_dz(0,v), rho_dz(NORTHEAST,v)) + physics_flux(DG+1)
            horiz_flux(v)%data(d)%elts(EDGE*id+UP+1) = u_dual_UP * interp (rho_dz(0,v), rho_dz(NORTH,    v)) + physics_flux(UP+1)
         end do
      end if
    end subroutine comput
 
    subroutine comp_SW
      
      implicit none
      
      integer  :: d, idS, idSW, idW
      integer  :: id_i, idS_i, idSW_i, idW_i, k
      integer  :: id_rhodz(0:SOUTHWEST)
      real(dp) :: circ_SW_LORT, circ_SW_UPLT, u_prim_RT_SW, u_prim_UP_SW, dz, dz0
      real(dp) :: csq(0:EDGE), phi(0:EDGE)
      real(dp) :: pv_mass(2)

      d = dom%id + 1

      idW  = id+W
      idSW = id+SW
      idS  = id+S

      id_i   = id+1
      idW_i  = idW+1
      idSW_i = idSW+1
      idS_i  = idS+1

      if (mode_split) then
         phi(0) = phi_node (d, id_i,   zlevels)
         phi(1) = phi_node (d, idW_i,  zlevels)
         phi(2) = phi_node (d, idSW_i, zlevels)
         phi(3) = phi_node (d, idS_i,  zlevels)
      else
         phi = 1.0_dp
      end if

      if (itype == 1) then ! scalar gradient flux
         h_flux(EDGE*idW+RT+1)  = -(scalar(idW_i) - scalar(id_i))   / dom%len%elts(EDGE*idW+RT+1)  &
              * dom%pedlen%elts(EDGE*idW+RT+1)
         h_flux(EDGE*idSW+DG+1) = -(scalar(id_i)  - scalar(idSW_i)) / dom%len%elts(EDGE*idSW+DG+1) &
              * dom%pedlen%elts(EDGE*idSW+DG+1)
         h_flux(EDGE*idS+UP+1)  = -(scalar(idS_i) - scalar(id_i))   / dom%len%elts(EDGE*idS+UP+1)  &
              * dom%pedlen%elts(EDGE*idS+UP+1)
      elseif (itype == 2) then ! external pressure gradient flux, (H + eta^n) grad(eta^(n+1)) * edge_length, in elliptic operator
         csq(0) = abs (topography%data(d)%elts(id_i))   * phi(0) + mass(id_i)
         csq(1) = abs (topography%data(d)%elts(idW_i))  * phi(1) + mass(idW_i)
         csq(2) = abs (topography%data(d)%elts(idSW_i)) * phi(2) + mass(idSW_i)
         csq(3) = abs (topography%data(d)%elts(idS_i))  * phi(3) + mass(idS_i)
         csq = grav_accel * csq

         h_flux(EDGE*idW+RT+1)  = - interp (csq(0), csq(1)) * &
              (scalar(idW_i)/phi(1) - scalar(id_i)/phi(0))/dom%len%elts(EDGE*idW+RT+1) * dom%pedlen%elts(EDGE*idW+RT+1)
         
         h_flux(EDGE*idSW+DG+1) = - interp (csq(0), csq(2)) * &
              (scalar(id_i)/phi(0)  - scalar(idSW_i)/phi(2))/dom%len%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1) 

         h_flux(EDGE*idS+UP+1)  = - interp (csq(0), csq(3)) * &
              (scalar(idS_i)/phi(3) - scalar(id_i)/phi(0))/dom%len%elts(EDGE*idS+UP+1) * dom%pedlen%elts(EDGE*idS+UP+1)
      elseif (itype == 3) then ! external pressure gradient
         h_flux(EDGE*idW+RT+1)  = - grav_accel * (scalar(idW_i)/phi(1) - scalar(id_i)/phi(0))   / dom%len%elts(EDGE*idW+RT+1) 
         h_flux(EDGE*idSW+DG+1) = - grav_accel * (scalar(id_i)/phi(0)  - scalar(idSW_i)/phi(2)) / dom%len%elts(EDGE*idSW+DG+1) 
         h_flux(EDGE*idS+UP+1)  = - grav_accel * (scalar(idS_i)/phi(3) - scalar(id_i)/phi(0))   / dom%len%elts(EDGE*idS+UP+1)
      elseif (itype == 4) then ! sum vertical flux
         h_flux(EDGE*idW+RT+1)  = 0.0_dp
         h_flux(EDGE*idSW+DG+1) = 0.0_dp
         h_flux(EDGE*idS+UP+1)  = 0.0_dp
         do k = 1, zlevels
            dz0 = q(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)

            dz = interp (dz0, q(S_MASS,k)%data(d)%elts(idW_i) + sol_mean(S_MASS,k)%data(d)%elts(idW_i))
            h_flux(EDGE*idW+RT+1) = h_flux(EDGE*idW+RT+1) + q(S_VELO,k)%data(d)%elts(EDGE*idW+RT+1) * dz

            dz = interp (dz0, q(S_MASS,k)%data(d)%elts(idSW_i) + sol_mean(S_MASS,k)%data(d)%elts(idSW_i))
            h_flux(EDGE*idSW+DG+1) = h_flux(EDGE*idSW+DG+1) + q(S_VELO,k)%data(d)%elts(EDGE*idSW+DG+1) * dz

            dz = interp (dz0, q(S_MASS,k)%data(d)%elts(idS_i) + sol_mean(S_MASS,k)%data(d)%elts(idS_i))
            h_flux(EDGE*idS+UP+1) = h_flux(EDGE*idS+UP+1) + q(S_VELO,k)%data(d)%elts(EDGE*idS+UP+1) * dz
         end do
         h_flux(EDGE*idW+RT+1)  = h_flux(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1) 
         h_flux(EDGE*idSW+DG+1) = h_flux(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1) 
         h_flux(EDGE*idS+UP+1)  = h_flux(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)  
      elseif (itype == 5) then ! divu
         h_flux(EDGE*idW+RT+1)  = velo(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1) 
         h_flux(EDGE*idSW+DG+1) = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
         h_flux(EDGE*idS+UP+1)  = velo(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)
      elseif (itype == 6) then ! surface pressure gradient flux (for vertical velocity)
         h_flux(EDGE*idW+RT+1)  = - (scalar(idW_i) - scalar(id_i))   * velo(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)  / 4
         h_flux(EDGE*idSW+DG+1) = - (scalar(id_i)  - scalar(idSW_i)) * velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1) / 4
         h_flux(EDGE*idS+UP+1)  = - (scalar(idS_i) - scalar(id_i))   * velo(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)  / 4
      elseif (itype == 7) then ! mass flux for vertical velocity
         dz0 = mass(id+1) + mean_m(id+1)

         dz = interp (dz0, mass(idW+1) + mean_m(idW+1)) 
         h_flux(EDGE*idW+RT+1) = velo(EDGE*idW+RT+1) * dom%pedlen%elts(EDGE*idW+RT+1) * dz

         dz = interp (dz0, mass(idSW+1) + mean_m(idSW+1)) 
         h_flux(EDGE*idSW+DG+1) = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1) * dz

         dz = interp (dz0, mass(idS+1) + mean_m(idS+1))
         h_flux(EDGE*idS+UP+1) = velo(EDGE*idS+UP+1) * dom%pedlen%elts(EDGE*idS+UP+1) * dz
      elseif (itype == 8) then ! vorticity
         u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
         u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1)
         u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
         u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
         u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)

         circ_SW_LORT =   u_prim_RT_SW + u_prim_UP_S  + u_prim_DG_SW
         circ_SW_UPLT = -(u_prim_RT_W  + u_prim_DG_SW + u_prim_UP_SW)

         vort(TRIAG*idSW+LORT+1) = circ_SW_LORT / dom%triarea%elts(TRIAG*idSW+LORT+1)
         vort(TRIAG*idSW+UPLT+1) = circ_SW_UPLT / dom%triarea%elts(TRIAG*idSW+UPLT+1)
      elseif (itype == 0) then ! standard
         id_rhodz = [ id, id, id, idS, idW, id, id, idSW ]
         do v = scalars(1), scalars(2)
            rho_dz(0:SOUTHWEST,v) = q(v,zlev)%data(d)%elts(id_rhodz+1) + sol_mean(v,zlev)%data(d)%elts(id_rhodz+1)
         end do

         u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
         u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
         u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
         u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1)
         u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)

         ! Potential vorticity
         circ_SW_LORT =   u_prim_RT_SW + u_prim_UP_S  + u_prim_DG_SW
         circ_SW_UPLT = -(u_prim_RT_W  + u_prim_DG_SW + u_prim_UP_SW)

         pv_mass(1) = &
              rho_dz(SOUTHWEST,S_MASS) * dom%areas%elts(idSW_i)%part(1) + &
              rho_dz(SOUTH,    S_MASS) * dom%areas%elts(idS_i )%part(3) + &
              rho_dz(0,        S_MASS) * dom%areas%elts(id_i  )%part(5)

         pv_mass(2) = &
              rho_dz(SOUTHWEST,S_MASS) * dom%areas%elts(idSW_i)%part(2) + &
              rho_dz(0,        S_MASS) * dom%areas%elts(id_i  )%part(4) + &
              rho_dz(WEST,     S_MASS) * dom%areas%elts(idW_i )%part(6)

         pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / pv_mass(1)
         pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_SW_UPLT) / pv_mass(2)

         qe(EDGE*idW +RT+1) = interp (pv_W_LORT , pv_SW_UPLT)
         qe(EDGE*idSW+DG+1) = interp (pv_SW_LORT, pv_SW_UPLT)
         qe(EDGE*idS +UP+1) = interp (pv_SW_LORT, pv_S_UPLT)

         ! Vorticity (for rotu velocity diffusion)
         vort(TRIAG*idSW+LORT+1) = circ_SW_LORT / dom%triarea%elts(TRIAG*idSW+LORT+1)
         vort(TRIAG*idSW+UPLT+1) = circ_SW_UPLT / dom%triarea%elts(TRIAG*idSW+UPLT+1)

         ! Scalar fluxes
         do v = scalars(1), scalars(2)
            physics_flux = physics_scalar_flux (q(:,1:zlevels), dom, id, idW, idSW, idS, v, zlev, .true.)

            horiz_flux(v)%data(d)%elts(EDGE*idW+RT +1) = u_dual_RT_W  * interp (rho_dz(0,v), rho_dz(WEST,     v)) &
                 + physics_flux(RT+1)
            
            horiz_flux(v)%data(d)%elts(EDGE*idSW+DG+1) = u_dual_DG_SW * interp (rho_dz(0,v), rho_dz(SOUTHWEST,v)) &
                 + physics_flux(DG+1)
            
            horiz_flux(v)%data(d)%elts(EDGE*idS+UP +1) = u_dual_UP_S  * interp (rho_dz(0,v), rho_dz(SOUTH,    v)) &
                 + physics_flux(UP+1)
         end do
      end if
    end subroutine comp_SW
  end subroutine step1
  

  subroutine post_step1 (dom, p, c, offs, dims, zlev)
    ! Correct values for vorticity and qe at pentagon points
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, c, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id, idS, idW, idSW, idN, idE, idNE
    integer  :: id_rhodz(8)
    real(dp) :: circ_LORT, circ_UPLT, circ_S_UPLT, circ_SW_LORT, circ_SW_UPLT, circ_W_LORT
    real(dp) :: pv_LORT, pv_UPLT, pv_S_UPLT, pv_SW_LORT, pv_SW_UPLT, pv_W_LORT
    real(dp) :: u_prim_RT, u_prim_RT_N, u_prim_RT_SW, u_prim_RT_W
    real(dp) :: u_prim_DG_SW, u_prim_UP, u_prim_UP_S, u_prim_UP_SW
    real(dp) :: pv_mass(3)
    real(dp) :: rho_dz(0:N_BDRY) 

    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case
    if (c == IJMINUS) then
       id   = idx ( 0,  0, offs, dims)
       idN  = idx ( 0,  1, offs, dims)
       idE  = idx ( 1,  0, offs, dims)      
       idS  = idx ( 0, -1, offs, dims)
       idW  = idx (-1,  0, offs, dims)
       idSW = idx (-1, -1, offs, dims)

       id_rhodz = [ id, idN, idE, idS, idW, id, id, idSW ]
       rho_dz(0:SOUTHWEST) = mass(id_rhodz+1) + mean_m(id_rhodz+1)

       u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1) 
       u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)

       circ_SW_LORT = u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW
       circ_W_LORT  = vort(TRIAG*idW+LORT+1) * dom%triarea%elts(TRIAG*idW+LORT+1)
       circ_S_UPLT  = vort(TRIAG*idS+UPLT+1) * dom%triarea%elts(TRIAG*idS+UPLT+1)

       vort(TRIAG*idSW+LORT+1) = circ_SW_LORT / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idSW+UPLT+1) = vort(TRIAG*idSW+LORT+1)

       pv_mass(1) = &
            rho_dz(WEST)  * dom%areas%elts(idW+1)%part(6) + &
            rho_dz(0)     * sum(dom%areas%elts(id+1)%part(4:5)) + &
            rho_dz(SOUTH) * dom%areas%elts(idS+1)%part(3)
      
       pv_mass(2) = &
            rho_dz(WEST)  * dom%areas%elts(idW+1)%part(1) + &
            rho_dz(0)     * dom%areas%elts(id+1)%part(3) + &
            rho_dz(NORTH) * dom%areas%elts(idN+1)%part(5)
       
       pv_mass(3) = &
            rho_dz(SOUTH) * dom%areas%elts(idS+1)%part(2) + &
            rho_dz(EAST)  * dom%areas%elts(idE+1)%part(4) + &
            rho_dz(0)     * dom%areas%elts(id+1)%part(6)

       pv_W_LORT  = (dom%coriolis%elts(TRIAG*idW+LORT +1) + circ_W_LORT ) / pv_mass(1)
       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / pv_mass(2)
       pv_S_UPLT  = (dom%coriolis%elts(TRIAG*idS+UPLT +1) + circ_S_UPLT ) / pv_mass(3)

       pv_SW_UPLT = pv_SW_LORT

       qe(EDGE*idW+RT+1) = interp (pv_W_LORT, pv_SW_UPLT)
       qe(EDGE*idS+UP+1) = interp (pv_S_UPLT, pv_SW_LORT)
    end if

    ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idS+UP+1) = 0 in this case
    if (c == IPLUSJMINUS) then 
       id   = idx (PATCH_SIZE,    0, offs, dims)
       idE  = idx (PATCH_SIZE+1,  0, offs, dims)
       idS  = idx (PATCH_SIZE,   -1, offs, dims)
       idW  = idx (PATCH_SIZE-1,  0, offs, dims)
       idNE = idx (PATCH_SIZE+1,  1, offs, dims)
       idSW = idx (PATCH_SIZE-1, -1, offs, dims)

       id_rhodz = [ id, id, idE, idS, idW, idNE, id, idSW ]
       rho_dz(0:SOUTHWEST) = mass(id_rhodz+1) + mean_m(id_rhodz+1)

       u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT    = velo(EDGE*id  +RT+1) * dom%len%elts(EDGE*id  +RT+1)

       circ_SW_LORT = - u_prim_RT + u_prim_RT_SW + u_prim_DG_SW
       circ_LORT    = vort(TRIAG*id  +LORT+1) * dom%triarea%elts(TRIAG*id  +LORT+1)
       circ_SW_UPLT = vort(TRIAG*idSW+UPLT+1) * dom%triarea%elts(TRIAG*idSW+UPLT+1)
       
       vort(TRIAG*idSW+LORT+1) = circ_SW_LORT / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idS +UPLT+1) = vort(TRIAG*idSW+LORT+1)

       pv_mass(1) = &
            rho_dz(SOUTHWEST) * dom%areas%elts(idSW+1)%part(1) + &
            rho_dz(SOUTH)     * dom%areas%elts(idS +1)%part(3) + &
            rho_dz(0)         * sum(dom%areas%elts(id+1)%part(5:6))

       pv_mass(2) = &
            rho_dz(0)         * dom%areas%elts(id  +1)%part(1) + &
            rho_dz(EAST)      * dom%areas%elts(idE +1)%part(3) + &
            rho_dz(NORTHEAST) * dom%areas%elts(idNE+1)%part(5)

       pv_mass(3) = &
            rho_dz(SOUTHWEST) * dom%areas%elts(idSW+1)%part(2) + &
            rho_dz(0)         * dom%areas%elts(id  +1)%part(4) + &
            rho_dz(WEST)      * dom%areas%elts(idW +1)%part(6)

       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / pv_mass(1)
       pv_LORT    = (dom%coriolis%elts(TRIAG*id+LORT  +1) + circ_LORT   ) / pv_mass(2)
       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_SW_UPLT) / pv_mass(3)

       pv_S_UPLT = pv_SW_LORT

       qe(EDGE*id  +RT+1) = interp (pv_S_UPLT,  pv_LORT)
       qe(EDGE*idSW+DG+1) = interp (pv_SW_LORT, pv_SW_UPLT)
    end if

    ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idW+RT+1) = 0 in this case
    if (c == IMINUSJPLUS) then
       id   = idx ( 0, PATCH_SIZE,   offs, dims)
       idN  = idx ( 0, PATCH_SIZE+1, offs, dims)
       idS  = idx ( 0, PATCH_SIZE-1, offs, dims)
       idW  = idx (-1, PATCH_SIZE,   offs, dims)
       idNE = idx ( 1, PATCH_SIZE+1, offs, dims)
       idSW = idx (-1, PATCH_SIZE-1, offs, dims)

       id_rhodz = [ id, idN, id, idS, idW, idNE, id, idSW ]
       rho_dz(0:SOUTHWEST) = mass(id_rhodz+1) + mean_m(id_rhodz+1)

       u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id  +UP+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)

       circ_SW_UPLT = u_prim_UP - u_prim_DG_SW - u_prim_UP_SW
       circ_UPLT    = vort(TRIAG*id  +UPLT+1) * dom%triarea%elts(TRIAG*id  +UPLT+1)
       circ_SW_LORT = vort(TRIAG*idSW+LORT+1) * dom%triarea%elts(TRIAG*idSW+LORT+1)

       vort(TRIAG*idSW+UPLT+1) = circ_SW_UPLT / dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       vort(TRIAG*idW +LORT+1) = vort(TRIAG*idSW+UPLT+1)

       pv_mass(1) = &
            rho_dz(SOUTHWEST) * dom%areas%elts(idSW+1)%part(2) + &
            rho_dz(0)         * sum(dom%areas%elts(id+1)%part(3:4)) + &
            rho_dz(WEST)      * dom%areas%elts(idW+1)%part(6)

       pv_mass(2) = &
            rho_dz(0)         * dom%areas%elts(id  +1)%part(2) + &
            rho_dz(NORTHEAST) * dom%areas%elts(idNE+1)%part(4) + &
            rho_dz(NORTH)     * dom%areas%elts(idN +1)%part(6)

       pv_mass(3) = &
            rho_dz(SOUTHWEST) * dom%areas%elts(idSW+1)%part(1) + &
            rho_dz(SOUTH)     * dom%areas%elts(idS +1)%part(3) + &
            rho_dz(0)         * dom%areas%elts(id  +1)%part(5)
       
       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_SW_UPLT) / pv_mass(1)
       pv_UPLT    = (dom%coriolis%elts(TRIAG*id+UPLT  +1) + circ_UPLT   ) / pv_mass(2)
       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / pv_mass(3)           

       pv_W_LORT = pv_SW_UPLT

       qe(EDGE*id  +UP+1) = interp (pv_W_LORT,  pv_UPLT)
       qe(EDGE*idSW+DG+1) = interp (pv_SW_LORT, pv_SW_UPLT)
    end if

    ! Parts 1, 2 of hexagon IJPLUS (upper right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*id+DG+1) = 0 in this case
    if (c == IJPLUS) then 
       id  = idx (PATCH_SIZE,   PATCH_SIZE,   offs, dims)
       idN = idx (PATCH_SIZE,   PATCH_SIZE+1, offs, dims)
       idE = idx (PATCH_SIZE+1, PATCH_SIZE,   offs, dims)
       idS = idx (PATCH_SIZE,   PATCH_SIZE-1, offs, dims)
       idW = idx (PATCH_SIZE-1, PATCH_SIZE,   offs, dims)

       id_rhodz(1:5) = [ id, idN, idE, idS, idW ]
       rho_dz(0:WEST) = mass(id_rhodz(1:5)+1) + mean_m(id_rhodz(1:5)+1)

       u_prim_RT   = velo(EDGE*id +RT+1) * dom%len%elts(EDGE*id +RT+1)
       u_prim_RT_N = velo(EDGE*idN+RT+1) * dom%len%elts(EDGE*idN+RT+1)
       u_prim_UP   = velo(EDGE*id +UP+1) * dom%len%elts(EDGE*id +UP+1)

       circ_LORT   = u_prim_RT - u_prim_RT_N - u_prim_UP
       circ_W_LORT = vort(TRIAG*idW+LORT+1) * dom%triarea%elts(TRIAG*idW+LORT+1)
       circ_S_UPLT = vort(TRIAG*idS+UPLT+1) * dom%triarea%elts(TRIAG*idS+UPLT+1)

       vort(TRIAG*id+LORT+1) = circ_LORT / dom%triarea%elts(TRIAG*id+LORT+1)
       vort(TRIAG*id+UPLT+1) = vort(TRIAG*id+LORT+1)

       pv_mass(1) = &
            rho_dz(EAST)  * dom%areas%elts(idE+1)%part(3) + &
            rho_dz(0)     * sum(dom%areas%elts(id+1)%part(1:2)) + &
            rho_dz(NORTH) * dom%areas%elts(idN+1)%part(6)

       pv_mass(2) = &
            rho_dz(WEST)  * dom%areas%elts(idW+1)%part(1) + &
            rho_dz(0)     * dom%areas%elts(id +1)%part(3) + &
            rho_dz(NORTH) * dom%areas%elts(idN+1)%part(5)

       pv_mass(3) = &
            rho_dz(SOUTH) * dom%areas%elts(idS+1)%part(2) + &
            rho_dz(EAST)  * dom%areas%elts(idE+1)%part(4) + &
            rho_dz(0)     * dom%areas%elts(id +1)%part(6)

       pv_LORT   = (dom%coriolis%elts(TRIAG*id+LORT +1) + circ_LORT  ) / pv_mass(1)                   
       pv_W_LORT = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_W_LORT) / pv_mass(2)
       pv_S_UPLT = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_S_UPLT) / pv_mass(3)

       pv_UPLT = pv_LORT

       qe(EDGE*id+RT+1) = interp (pv_LORT, pv_S_UPLT)
       qe(EDGE*id+UP+1) = interp (pv_UPLT, pv_W_LORT)
    end if
  end subroutine post_step1

  
  subroutine scalar_trend (dom, i, j, zlev, offs, dims)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id_i) >= TRSK) then
       dscalar(id_i) = - div (h_flux, dom, i, j, offs, dims)
    else
       dscalar(id_i) = 0.0_dp
    end if
  end subroutine scalar_trend

  
  subroutine du_source (dom, i, j, zlev, offs, dims)
    ! Edge integrated source (non gradient) terms in velocity trend
    ! [Aechtner thesis page 56, Kevlahan, Dubos and Aechtner (2015)]
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id
    integer  :: id_e(EDGE)
    real(dp) :: Qperp_e(EDGE), physics(EDGE)
    
    id   =  idx (i, j, offs, dims)
    id_e = id_edge (id)
    
    if (dom%mask_n%elts(id+1) >= TRSK) then
       ! Calculate Q_perp
       Qperp_e = Qperp (dom, i, j, z_null, offs, dims)

       ! Calculate physics
       physics = physics_velo_source (dom, i, j, zlev, offs, dims)

       ! Trend
       dvelo(id_e) = - Qperp_e + physics * dom%len%elts(id_e)
    else
       dvelo(id_e) = 0.0_dp
    end if
  end subroutine du_source

  
  subroutine du_grad (dom, i, j, zlev, offs, dims)
    ! Add gradients of Bernoulli and Exner to dvelo [DYNAMICO (23)-(25)]
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id, idE, idN, idNE
    integer  :: id_e(EDGE)
    real(dp) :: gradB(EDGE), gradE(EDGE), theta_e(EDGE)
    real(dp) :: rho_dz(0:NORTHEAST), rho_dz_theta(0:NORTHEAST) , theta(0:NORTHEAST) 

    id   = idx (i, j, offs, dims)
    id_e = id_edge (id)
    
    if (dom%mask_n%elts(id+1) >= TRSK) then
       idE  = idx (i+1, j,   offs, dims) 
       idN  = idx (i,   j+1, offs, dims)
       idNE = idx (i+1, j+1, offs, dims)

       rho_dz(0:NORTHEAST)       = mean_m([id,idN,idE,id,id,idNE]+1) + mass([id,idN,idE,id,id,idNE]+1)
       rho_dz_theta(0:NORTHEAST) = mean_t([id,idN,idE,id,id,idNE]+1) + temp([id,idN,idE,id,id,idNE]+1)

       ! See DYNAMICO between (23)-(25), geopotential still known from step1_up
       ! the theta multiplying the Exner gradient is the edge-averaged non-mass-weighted potential temperature
       theta(0)         = rho_dz_theta(0)         / rho_dz(0)
       theta(EAST)      = rho_dz_theta(EAST)      / rho_dz(EAST)
       theta(NORTHEAST) = rho_dz_theta(NORTHEAST) / rho_dz(NORTHEAST)
       theta(NORTH)     = rho_dz_theta(NORTH)     / rho_dz(NORTH)

       ! Interpolate potential temperature to edges
       theta_e(RT+1) = interp (theta(0), theta(EAST))      
       theta_e(DG+1) = interp (theta(0), theta(NORTHEAST)) 
       theta_e(UP+1) = interp (theta(0), theta(NORTH))     

       ! Calculate gradients
       gradB = gradi_e (bernoulli, dom, i, j, offs, dims)
       gradE = gradi_e (exner,     dom, i, j, offs, dims)

       ! Trend
       dvelo(id_e) = dvelo(id_e) / dom%len%elts(id_e) - gradB - theta_e * gradE
    else
       dvelo(id_e) = 0.0_dp
    end if
  end subroutine du_grad

  
  function Qperp (dom, i, j, zlev, offs, dims) result(val)
    ! Compute energy-conserving edge integrated Qperp [Aechtner thesis page 44]
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val(EDGE)
    
    integer  :: id, idNW, idN, idNE, idW, idE, idSW, idS, idSE
    real(dp) :: wgt1(5), wgt2(5)

    id   = idx (i, j, offs, dims)

    idNW = idx (i-1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idSE = idx (i+1, j-1, offs, dims)

    wgt1 = get_weights(dom, id,  0)
    wgt2 = get_weights(dom, idE, 3)

    val(RT+1) = &
         h_mflux(EDGE*id+DG+1)   * interp (qe(EDGE*id+DG+1),   qe(EDGE*id+RT+1))*wgt1(1) + &
         h_mflux(EDGE*id+UP+1)   * interp (qe(EDGE*id+UP+1),   qe(EDGE*id+RT+1))*wgt1(2) + &
         h_mflux(EDGE*idW+RT+1)  * interp (qe(EDGE*idW+RT+1),  qe(EDGE*id+RT+1))*wgt1(3) + &
         h_mflux(EDGE*idSW+DG+1) * interp (qe(EDGE*idSW+DG+1), qe(EDGE*id+RT+1))*wgt1(4) + &
         h_mflux(EDGE*idS+UP+1)  * interp (qe(EDGE*idS+UP+1),  qe(EDGE*id+RT+1))*wgt1(5) + &
         h_mflux(EDGE*idS+DG+1)  * interp (qe(EDGE*idS+DG+1),  qe(EDGE*id+RT+1))*wgt2(1) + &
         h_mflux(EDGE*idSE+UP+1) * interp (qe(EDGE*idSE+UP+1), qe(EDGE*id+RT+1))*wgt2(2) + &
         h_mflux(EDGE*idE+RT+1)  * interp (qe(EDGE*idE+RT+1),  qe(EDGE*id+RT+1))*wgt2(3) + &
         h_mflux(EDGE*idE+DG+1)  * interp (qe(EDGE*idE+DG+1),  qe(EDGE*id+RT+1))*wgt2(4) + &
         h_mflux(EDGE*idE+UP+1)  * interp (qe(EDGE*idE+UP+1),  qe(EDGE*id+RT+1))*wgt2(5)

    wgt1 = get_weights(dom, id,   1)
    wgt2 = get_weights(dom, idNE, 4)

    val(DG+1) = &
         h_mflux(EDGE*id+UP+1)   * interp (qe(EDGE*id+UP+1),   qe(EDGE*id+DG+1))*wgt1(1) + &
         h_mflux(EDGE*idW+RT+1)  * interp (qe(EDGE*idW+RT+1),  qe(EDGE*id+DG+1))*wgt1(2) + &
         h_mflux(EDGE*idSW+DG+1) * interp (qe(EDGE*idSW+DG+1), qe(EDGE*id+DG+1))*wgt1(3) + &
         h_mflux(EDGE*idS+UP+1)  * interp (qe(EDGE*idS+UP+1),  qe(EDGE*id+DG+1))*wgt1(4) + &
         h_mflux(EDGE*id+RT+1)   * interp (qe(EDGE*id+RT+1),   qe(EDGE*id+DG+1))*wgt1(5) + &
         h_mflux(EDGE*idE+UP+1)  * interp (qe(EDGE*idE+UP+1),  qe(EDGE*id+DG+1))*wgt2(1) + &
         h_mflux(EDGE*idNE+RT+1) * interp (qe(EDGE*idNE+RT+1), qe(EDGE*id+DG+1))*wgt2(2) + &
         h_mflux(EDGE*idNE+DG+1) * interp (qe(EDGE*idNE+DG+1), qe(EDGE*id+DG+1))*wgt2(3) + &
         h_mflux(EDGE*idNE+UP+1) * interp (qe(EDGE*idNE+UP+1), qe(EDGE*id+DG+1))*wgt2(4) + &
         h_mflux(EDGE*idN+RT+1)  * interp (qe(EDGE*idN+RT+1),  qe(EDGE*id+DG+1))*wgt2(5)

    wgt1 = get_weights(dom, id,  2)
    wgt2 = get_weights(dom, idN, 5)

    val(UP+1) = &
         h_mflux(EDGE*idW+RT+1)  * interp (qe(EDGE*idW+RT+1),  qe(EDGE*id+UP+1))*wgt1(1) + &
         h_mflux(EDGE*idSW+DG+1) * interp (qe(EDGE*idSW+DG+1), qe(EDGE*id+UP+1))*wgt1(2) + &
         h_mflux(EDGE*idS+UP+1)  * interp (qe(EDGE*idS+UP+1),  qe(EDGE*id+UP+1))*wgt1(3) + &
         h_mflux(EDGE*id+RT+1)   * interp (qe(EDGE*id+RT+1),   qe(EDGE*id+UP+1))*wgt1(4) + &
         h_mflux(EDGE*id+DG+1)   * interp (qe(EDGE*id+DG+1),   qe(EDGE*id+UP+1))*wgt1(5) + &
         h_mflux(EDGE*idN+RT+1)  * interp (qe(EDGE*idN+RT+1),  qe(EDGE*id+UP+1))*wgt2(1) + &
         h_mflux(EDGE*idN+DG+1)  * interp (qe(EDGE*idN+DG+1),  qe(EDGE*id+UP+1))*wgt2(2) + &
         h_mflux(EDGE*idN+UP+1)  * interp (qe(EDGE*idN+UP+1),  qe(EDGE*id+UP+1))*wgt2(3) + &
         h_mflux(EDGE*idNW+RT+1) * interp (qe(EDGE*idNW+RT+1), qe(EDGE*id+UP+1))*wgt2(4) + &
         h_mflux(EDGE*idW+DG+1)  * interp (qe(EDGE*idW+DG+1),  qe(EDGE*id+UP+1))*wgt2(5)
  end function Qperp

  
  function get_weights (dom, id, offs) result(val)
    ! Weights for Qperp computation [Aechtner thesis page 44]
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: id, offs
    real(dp)                    :: val(5)

    integer  :: i
    real(dp) :: wgt(5)

    wgt(1) = dom%areas%elts(id+1)%part(1+offs)
    do i = 2, 5
       wgt(i) = wgt(i-1) + dom%areas%elts(id+1)%part(modulo(i+offs-1,6)+1)
    end do
    wgt = 0.5_dp - wgt*dom%areas%elts(id+1)%hex_inv
    
    val = [ wgt(1), -wgt(2), wgt(3), -wgt(4), wgt(5) ]
  end function get_weights

  
  subroutine comp_offs3 (dom, p, offs, dims)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p
    integer,      intent(out)   :: offs(0:N_BDRY)
    integer,      intent(out)   :: dims(2,N_BDRY)

    integer :: i, n

    offs(0) = dom%patch%elts(p+1)%elts_start
    do i = 1, N_BDRY
       n = dom%patch%elts(p+1)%neigh(i)
       if (n > 0) then ! regular patch
          offs(i)  = dom%patch%elts(n+1)%elts_start
          dims(:,i) = PATCH_SIZE
       elseif (n < 0) then
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
end module ops_mod
