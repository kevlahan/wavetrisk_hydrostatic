module ops_mod
  use test_case_mod
  implicit none
  real(8) :: sum_mass, sum_temp, totaldmass, totalabsdmass, totaldtemp, totalabsdtemp
contains
  subroutine init_ops_mod
    implicit none
    logical :: initialized = .False.
    if (initialized) return ! initialize only once
    call init_domain_mod
    initialized = .true.
  end subroutine init_ops_mod

  subroutine step1 (dom, p, zlev)
    implicit none
    type(Domain) :: dom
    integer      :: p, zlev

    integer :: i, j, id,  n, e, s, w, ne, sw
    integer, dimension(0:N_BDRY) :: offs
    integer, dimension(2,N_BDRY) :: dims

    real(8) :: u_prim_up, u_dual_up, u_prim_dg, u_dual_dg, u_prim_RT, u_dual_rt
    real(8) :: u_prim_UP_S, u_dual_UP_S, u_prim_DG_SW, u_dual_DG_SW, u_prim_RT_W, u_dual_RT_W
    real(8) :: circ_LORT_W, circ_UPLT_S, circ_LORT, circ_UPLT, pv_LORT, pv_UPLT, pv_UPLT_S, pv_LORT_W,  pv_LORT_SW, pv_UPLT_SW

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

    if (W_bdry .or. S_bdry) call comp_ijmin

    w = -1; sw = s+w
    do id = offs(0)+1, offs(0)+LAST-1
       call comput
       if (S_bdry) call comp_ijmin
    end do

    e = offs(EAST); ne = dims(1,EAST) + e
    id = offs(0)+LAST
    call comput
    if (S_bdry) call comp_ijmin

    s = -PATCH_SIZE
    do j = 2, PATCH_SIZE-1
       id = offs(0) + PATCH_SIZE*(j-1)
       e = +1; ne = n+e
       w = offs(WEST) + (dims(1,WEST)-PATCH_SIZE)*(j-1) ! Correct for dimension smaller than patch if boundary
       sw = w-dims(1,WEST)
       call comput
       if (W_bdry) call comp_ijmin

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
    if (W_bdry) call comp_ijmin

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
          if (W_bdry) call comp_ijmin

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
       if (S_bdry) call comp_ijmin

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
    subroutine comp_ijmin
      implicit none
      integer :: idS, idSW, idW
      integer :: id_i, idS_i, idSW_i, idW_i
      real(8) :: circ_LORT_SW, circ_UPLT_SW, u_prim_RT_SW, u_prim_UP_SW
      real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics

      interface
         function physics_scalar_flux (dom, id, idE, idNE, idN, type)
           use domain_mod
           real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
           type(Domain)                             :: dom
           integer                                  :: id, idE, idNE, idN
           logical, optional                        :: type
         end function physics_scalar_flux
      end interface

      idW  = id+W
      idSW = id+SW
      idS  = id+S

      id_i   = id+1
      idW_i  = idW+1
      idSW_i = idSW+1
      idS_i  = idS+1

      u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
      u_prim_UP_SW = velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)

      circ_LORT_SW =   u_prim_RT_SW + u_prim_UP_S  + u_prim_DG_SW
      circ_UPLT_SW = -(u_prim_RT_W  + u_prim_DG_SW + u_prim_UP_SW)

      vort(TRIAG*idW+LORT+1) = circ_LORT_W/dom%triarea%elts(TRIAG*idW+LORT+1) 
      vort(TRIAG*idS+UPLT+1) = circ_UPLT_S/dom%triarea%elts(TRIAG*idS+UPLT+1) 

      pv_LORT_SW = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_LORT_SW)/( &
           mass(idSW_i)*dom%areas%elts(idSW_i)%part(1) + &
           mass(idS_i)*dom%areas%elts(idS_i)%part(3) + &
           mass(id_i)*dom%areas%elts(id_i)%part(5))

      pv_UPLT_SW = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_UPLT_SW)/( &
           mass(idSW_i)*dom%areas%elts(idSW_i)%part(2) + &
           mass(id_i)*dom%areas%elts(id_i)%part(4) + &
           mass(idW_i)*dom%areas%elts(idW_i)%part(6))

      qe(EDGE*idW+RT+1)  = interp (pv_LORT_W , pv_UPLT_SW)
      qe(EDGE*idSW+DG+1) = interp (pv_LORT_SW, pv_UPLT_SW)
      qe(EDGE*idS+UP+1)  = interp (pv_LORT_SW, pv_UPLT_S)

      ! Mass and temperature fluxes
      physics = physics_scalar_flux (dom, id, idW, idSW, idS, .true.)

      h_mflux(EDGE*idW+RT+1)  = u_dual_RT_W  * interp (mass(id_i), mass(idW_i))  + physics(S_MASS,RT+1)
      h_mflux(EDGE*idSW+DG+1) = u_dual_DG_SW * interp (mass(id_i), mass(idSW_i)) + physics(S_MASS,DG+1)
      h_mflux(EDGE*idS+UP+1)  = u_dual_UP_S  * interp (mass(id_i), mass(idS_i))  + physics(S_MASS,UP+1)

      h_tflux(EDGE*idW+RT+1)  = u_dual_RT_W  * interp (temp(id_i), temp(idW_i))  + physics(S_TEMP,RT+1)
      h_tflux(EDGE*idSW+DG+1) = u_dual_DG_SW * interp (temp(id_i), temp(idSW_i)) + physics(S_TEMP,DG+1)
      h_tflux(EDGE*idS+UP+1)  = u_dual_UP_S  * interp (temp(id_i), temp(idS_i))  + physics(S_TEMP,UP+1)
    end subroutine comp_ijmin

    subroutine comput
      ! Computes physical quantities during upward integration
      implicit none
      type (Coord)                             :: x_e, x_i, vel
      integer                                  :: idE, idN, idNE, idS, idSW, idW
      integer                                  :: id_i, idE_i, idN_i, idNE_i, idS_i, idSW_i, idW_i
      real(8)                                  :: kinetic_energy, Phi_k, circ_LORT, circ_UPLT
      real(8)                                  :: u_prim_UP_E, u_prim_RT_N, u_prim_DG_W, u_prim_DG_S
      real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics
      type (Coord), dimension(6)               :: hex_nodes

      interface
         function physics_scalar_flux (dom, id, idE, idNE, idN, type)
           use domain_mod
           real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
           type(Domain)                             :: dom
           integer                                  :: id, idE, idNE, idN
           logical, optional                        :: type
         end function physics_scalar_flux
      end interface
      
      idE  = id+E
      idN  = id+N
      idNE = id+NE
      idS  = id+S
      idSW = id+SW
      idW  = id+W

      id_i   = id+1
      idE_i  = idE+1
      idN_i  = idN+1
      idNE_i = idNE+1
      idS_i  = idS+1
      idSW_i = idSW+1
      idW_i  = idW+1

      ! Find the velocity on primal and dual grids
      u_prim_RT    = velo(EDGE*id  +RT+1)*dom%len%elts(EDGE*id+RT+1)
      u_prim_UP    = velo(EDGE*id  +UP+1)*dom%len%elts(EDGE*id+UP+1)
      u_prim_DG    = velo(EDGE*id  +DG+1)*dom%len%elts(EDGE*id+DG+1)
      u_prim_RT_W  = velo(EDGE*idW +RT+1)*dom%len%elts(EDGE*idW+RT+1)
      u_prim_UP_S  = velo(EDGE*idS +UP+1)*dom%len%elts(EDGE*idS+UP+1)
      u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
      u_prim_UP_E  = velo(EDGE*idE +UP+1)*dom%len%elts(EDGE*idE+UP+1)
      u_prim_RT_N  = velo(EDGE*idN +RT+1)*dom%len%elts(EDGE*idN+RT+1)
      u_prim_DG_W  = velo(EDGE*idW +DG+1)*dom%len%elts(EDGE*idW+DG+1)
      u_prim_DG_S  = velo(EDGE*idS +DG+1)*dom%len%elts(EDGE*idS+DG+1)

      u_dual_RT    = velo(EDGE*id  +RT+1)*dom%pedlen%elts(EDGE*id+RT+1)
      u_dual_UP    = velo(EDGE*id  +UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
      u_dual_DG    = velo(EDGE*id  +DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
      u_dual_RT_W  = velo(EDGE*idW +RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)
      u_dual_UP_S  = velo(EDGE*idS +UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)
      u_dual_DG_SW = velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1)

      ! Formula from TRiSK 
      kinetic_energy = (u_prim_UP*u_dual_UP + u_prim_DG*u_dual_DG + u_prim_RT*u_dual_RT + &
           u_prim_UP_S*u_dual_UP_S + u_prim_DG_SW*u_dual_DG_SW + u_prim_RT_W*u_dual_RT_W) * dom%areas%elts(id_i)%hex_inv/4

      ! Interpolate geopotential from interfaces to level
      Phi_k = interp (dom%geopot%elts(id_i), dom%geopot_lower%elts(id_i))

      ! Bernoulli function
      if (compressible) then 
         bernoulli(id_i) = kinetic_energy + Phi_k
      else 
         bernoulli(id_i) = kinetic_energy + Phi_k + dom%press%elts(id_i)/ref_density
      end if

      ! Exner function in incompressible case from geopotential
      if (.not. compressible) exner(id_i) = -Phi_k

      ! Calculate div(u) for velocity diffusion
      if (Laplace_order /= 0) &
           divu(id_i) = (u_dual_RT-u_dual_RT_W + u_dual_DG_SW-u_dual_DG + u_dual_UP-u_dual_UP_S) * dom%areas%elts(id_i)%hex_inv 

      circ_LORT   =   u_prim_RT    + u_prim_UP_E + u_prim_DG 
      circ_UPLT   = -(u_prim_DG    + u_prim_UP   + u_prim_RT_N)
      circ_LORT_W =   u_prim_RT_W  + u_prim_UP   + u_prim_DG_W
      circ_UPLT_S = -(u_prim_RT    + u_prim_DG_S + u_prim_UP_S)

      vort(TRIAG*id+LORT+1) = circ_LORT/dom%triarea%elts(TRIAG*id+LORT+1) 
      vort(TRIAG*id+UPLT+1) = circ_UPLT/dom%triarea%elts(TRIAG*id+UPLT+1)

      pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT)/( &
           mass(id_i)*dom%areas%elts(id_i)%part(1) + &
           mass(idE_i)*dom%areas%elts(idE_i)%part(3) + &
           mass(idNE_i)*dom%areas%elts(idNE_i)%part(5))

      pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + circ_UPLT)/( &
           mass(id_i)*dom%areas%elts(id_i)%part(2) + &
           mass(idNE_i)*dom%areas%elts(idNE_i)%part(4) + &
           mass(idN_i)*dom%areas%elts(idN_i)%part(6))

      pv_LORT_W = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_LORT_W)/( &
           mass(idW_i)*dom%areas%elts(idW_i)%part(1) + &
           mass(id_i)*dom%areas%elts(id_i)%part(3) + &
           mass(idN_i)*dom%areas%elts(idN_i)%part(5))

      pv_UPLT_S = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_UPLT_S)/( &
           mass(idS_i)*dom%areas%elts(idS_i)%part(2) + &
           mass(idE_i)*dom%areas%elts(idE_i)%part(4) + &
           mass(id_i)*dom%areas%elts(id_i)%part(6))

      qe(EDGE*id+RT+1) = interp (pv_UPLT_S, pv_LORT)
      qe(EDGE*id+DG+1) = interp (pv_UPLT,   pv_LORT)
      qe(EDGE*id+UP+1) = interp (pv_UPLT,   pv_LORT_W)

      ! Mass and temperature fluxes
      physics = physics_scalar_flux (dom, id, idE, idNE, idN)

      h_mflux(EDGE*id+RT+1) = u_dual_RT * interp (mass(id_i), mass(idE_i))  + physics(S_MASS,RT+1)
      h_mflux(EDGE*id+DG+1) = u_dual_DG * interp (mass(id_i), mass(idNE_i)) + physics(S_MASS,DG+1)
      h_mflux(EDGE*id+UP+1) = u_dual_UP * interp (mass(id_i), mass(idN_i))  + physics(S_MASS,UP+1)

      h_tflux(EDGE*id+RT+1) = u_dual_RT * interp (temp(id_i), temp(idE_i))  + physics(S_TEMP,RT+1)
      h_tflux(EDGE*id+DG+1) = u_dual_DG * interp (temp(id_i), temp(idNE_i)) + physics(S_TEMP,DG+1)
      h_tflux(EDGE*id+UP+1) = u_dual_UP * interp (temp(id_i), temp(idN_i))  + physics(S_TEMP,UP+1)
    end subroutine comput
  end subroutine step1

  subroutine post_step1 (dom, p, c, offs, dims, zlev)
    ! Correct values for vorticity and qe at pentagon points
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer  :: id, idS, idW, idSW, idN, idE, idNE
    real(8)  :: pv_LORT_W, pv_UPLT_S, pv_LORT, pv_UPLT, pv_LORT_SW, pv_UPLT_SW, pv 
    real(8)  :: circ_LORT, circ_LORT_SW, circ_UPLT_SW, circ_LORT_W, circ_UPLT, circ_UPLT_S
    real(8)  :: u_prim_RT, u_prim_RT_N, u_prim_RT_SW, u_prim_RT_W, u_prim_DG_SW, u_prim_UP, u_prim_UP_S, u_prim_UP_SW

    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case
    if (c == IJMINUS) then
       id   = idx( 0,  0, offs, dims)
       idSW = idx(-1, -1, offs, dims)
       idW  = idx(-1,  0, offs, dims)
       idS  = idx( 0, -1, offs, dims)
       idN  = idx( 0,  1, offs, dims)
       idE  = idx( 1,  0, offs, dims)

       u_prim_RT_W  = velo(EDGE*idW +RT+1)*dom%len%elts(EDGE*idW +RT+1)
       u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1) 
       u_prim_UP_S  = velo(EDGE*idS +UP+1)*dom%len%elts(EDGE*idS +UP+1)

       circ_LORT_SW = u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW
       circ_LORT_W  = vort(TRIAG*idW+LORT+1)*dom%triarea%elts(TRIAG*idW+LORT+1)
       circ_UPLT_S  = vort(TRIAG*idS+UPLT+1)*dom%triarea%elts(TRIAG*idS+UPLT+1)

       vort(TRIAG*idSW+LORT+1) = circ_LORT_SW/dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idSW+UPLT+1) = vort(TRIAG*idSW+LORT+1)

       pv_LORT_SW = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_LORT_SW)/ &
            (mass(idW+1)*dom%areas%elts(idW+1)%part(6) &
            + mass(id+1)*sum(dom%areas%elts(id+1)%part(4:5)) &
            + mass(idS+1)*dom%areas%elts(idS+1)%part(3))

       pv_LORT_W = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_LORT_W)/ &
            (mass(idW+1)*dom%areas%elts(idW+1)%part(1) &
            + mass(id+1)*dom%areas%elts(id+1)%part(3) &
            + mass(idN+1)*dom%areas%elts(idN+1)%part(5))

       pv_UPLT_S = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_UPLT_S)/ &
            (mass(idS+1)*dom%areas%elts(idS+1)%part(2) &
            + mass(idE+1)*dom%areas%elts(idE+1)%part(4) &
            + mass(id+1)*dom%areas%elts(id+1)%part(6))

       pv_UPLT_SW = pv_LORT_SW

       qe(EDGE*idW+RT+1) = interp (pv_LORT_W, pv_UPLT_SW)
       qe(EDGE*idS+UP+1) = interp (pv_UPLT_S, pv_LORT_SW)
    end if

    ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idS+UP+1) = 0 in this case
    if (c == IPLUSJMINUS) then 
       id   = idx(PATCH_SIZE,    0, offs, dims)
       idSW = idx(PATCH_SIZE-1, -1, offs, dims)
       idS  = idx(PATCH_SIZE,   -1, offs, dims)
       idW  = idx(PATCH_SIZE-1,  0, offs, dims)
       idE  = idx(PATCH_SIZE+1,  0, offs, dims)
       idNE = idx(PATCH_SIZE+1,  1, offs, dims)

       u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT    = velo(EDGE*id  +RT+1)*dom%len%elts(EDGE*id  +RT+1)

       circ_LORT_SW = - u_prim_RT + u_prim_RT_SW + u_prim_DG_SW 
       circ_LORT    = vort(TRIAG*id+LORT+1)*dom%triarea%elts(TRIAG*id+LORT+1)
       circ_UPLT_SW = vort(TRIAG*idSW+UPLT+1)*dom%triarea%elts(TRIAG*idSW+UPLT+1)

       vort(TRIAG*idSW+LORT+1) = circ_LORT_SW/dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idS +UPLT+1) = vort(TRIAG*idSW+LORT+1)

       pv_LORT_SW = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_LORT_SW)/ &
            (mass(idSW+1)*dom%areas%elts(idSW+1)%part(1) + &
            mass(idS+1)*dom%areas%elts(idS+1)%part(3) + &
            mass(id+1)*sum(dom%areas%elts(id+1)%part(5:6)))

       pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT)/ &
            (mass(id+1)*dom%areas%elts(id+1)%part(1) &
            + mass(idE+1)*dom%areas%elts(idE+1)%part(3) &
            + mass(idNE+1)*dom%areas%elts(idNE+1)%part(5))

       pv_UPLT_SW = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_UPLT_SW)/ &
            (mass(idSW+1)*dom%areas%elts(idSW+1)%part(2) &
            + mass(id+1)*dom%areas%elts(id+1)%part(4) &
            + mass(idW+1)*dom%areas%elts(idW+1)%part(6))

       pv_UPLT_S = pv_LORT_SW

       qe(EDGE*id  +RT+1) = interp (pv_UPLT_S,  pv_LORT)
       qe(EDGE*idSW+DG+1) = interp (pv_LORT_SW, pv_UPLT_SW)
    end if

    ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idW+RT+1) = 0 in this case
    if (c == IMINUSJPLUS) then
       id   = idx( 0, PATCH_SIZE,   offs, dims)
       idSW = idx(-1, PATCH_SIZE-1, offs, dims)
       idW  = idx(-1, PATCH_SIZE,   offs, dims)
       idS  = idx( 0, PATCH_SIZE-1, offs, dims)
       idN  = idx( 0, PATCH_SIZE+1, offs, dims)
       idNE = idx( 1, PATCH_SIZE+1, offs, dims)

       u_prim_UP    = velo(EDGE*id  +UP+1)*dom%len%elts(EDGE*id  +UP+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_SW = velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)

       circ_UPLT_SW = u_prim_UP - u_prim_DG_SW - u_prim_UP_SW
       circ_UPLT    = vort(TRIAG*id+UPLT+1)*dom%triarea%elts(TRIAG*id+UPLT+1)
       circ_LORT_SW = vort(TRIAG*idSW+LORT+1)*dom%triarea%elts(TRIAG*idSW+LORT+1)

       vort(TRIAG*idSW+UPLT+1) = circ_UPLT_SW/dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       vort(TRIAG*idW +LORT+1) = vort(TRIAG*idSW+UPLT+1)

       pv_UPLT_SW = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_UPLT_SW)/ &
            (mass(idSW+1)*dom%areas%elts(idSW+1)%part(2) &
            + mass(id+1)*sum(dom%areas%elts(id+1)%part(3:4)) &
            + mass(idW+1)*dom%areas%elts(idW+1)%part(6))

       pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + circ_UPLT)/ &
            (mass(id+1)*dom%areas%elts(id+1)%part(2) &
            + mass(idNE+1)*dom%areas%elts(idNE+1)%part(4) &
            + mass(idN+1)*dom%areas%elts(idN+1)%part(6))

       pv_LORT_SW = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_LORT_SW)/ &
            (mass(idSW+1)*dom%areas%elts(idSW+1)%part(1) &
            + mass(idS+1)*dom%areas%elts(idS+1)%part(3) &
            + mass(id+1)*dom%areas%elts(id+1)%part(5))

       pv_LORT_W = pv_UPLT_SW

       qe(EDGE*id  +UP+1) = interp (pv_LORT_W,  pv_UPLT)
       qe(EDGE*idSW+DG+1) = interp (pv_LORT_SW, pv_UPLT_SW)
    end if

    ! Parts 1, 2 of hexagon IJPLUS (upper right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*id+DG+1) = 0 in this case
    if (c == IJPLUS) then 
       id  = idx(PATCH_SIZE,   PATCH_SIZE,   offs, dims)
       idN = idx(PATCH_SIZE,   PATCH_SIZE+1, offs, dims)
       idE = idx(PATCH_SIZE+1, PATCH_SIZE,   offs, dims)
       idS = idx(PATCH_SIZE,   PATCH_SIZE-1, offs, dims)
       idW = idx(PATCH_SIZE-1, PATCH_SIZE,   offs, dims)

       u_prim_RT   = velo(EDGE*id +RT+1)*dom%len%elts(EDGE*id +RT+1)
       u_prim_RT_N = velo(EDGE*idN+RT+1)*dom%len%elts(EDGE*idN+DG+1)
       u_prim_UP   = velo(EDGE*id +UP+1)*dom%len%elts(EDGE*id +UP+1)

       circ_LORT   = u_prim_RT - u_prim_RT_N - u_prim_UP
       circ_LORT_W = vort(TRIAG*idW+LORT+1)*dom%triarea%elts(TRIAG*idW+LORT+1)
       circ_UPLT_S = vort(TRIAG*idS+UPLT+1)*dom%triarea%elts(TRIAG*idS+UPLT+1)

       vort(TRIAG*id+LORT+1) = circ_LORT/dom%triarea%elts(TRIAG*id+LORT+1)
       vort(TRIAG*id+UPLT+1) = vort(TRIAG*id+LORT+1)

       pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT)/ &          
            (mass(idE+1)*dom%areas%elts(idE+1)%part(3) + &
            mass(id+1)*sum(dom%areas%elts(id+1)%part(1:2)) + &
            mass(idN+1)*dom%areas%elts(idN+1)%part(6))

       pv_LORT_W = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_LORT_W)/ &
            (mass(idW+1)*dom%areas%elts(idW+1)%part(1) &
            + mass(id+1)*dom%areas%elts(id+1)%part(3) &
            + mass(idN+1)*dom%areas%elts(idN+1)%part(5))

       pv_UPLT_S = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_UPLT_S)/ &
            (mass(idS+1)*dom%areas%elts(idS+1)%part(2) &
            + mass(idE+1)*dom%areas%elts(idE+1)%part(4) &
            + mass(id+1)*dom%areas%elts(id+1)%part(6))

       pv_UPLT = pv_LORT

       qe(EDGE*id+RT+1) = interp (pv_LORT, pv_UPLT_S)
       qe(EDGE*id+UP+1) = interp (pv_UPLT, pv_LORT_W)
    end if
  end subroutine post_step1

  subroutine post_vort (dom, p, c, offs, dims, zlev)
    ! Correct values for vorticity at pentagon points (used in diffusion time step)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idN, idS, idSW, idW
    real(8) :: circ_LORT, circ_LORT_SW, circ_UPLT_SW
    real(8) :: u_prim_DG_SW, u_prim_RT, u_prim_RT_N, u_prim_RT_SW,  u_prim_RT_W, u_prim_UP, u_prim_UP_S, u_prim_UP_SW

    if (c == IJMINUS) then ! Parts 4, 5 of hexagon IJMINUS (SW corner of lozenge) combined to form pentagon
       id   = idx( 0,  0, offs, dims)
       idS  = idx( 0, -1, offs, dims)
       idSW = idx(-1, -1, offs, dims)
       idW  = idx(-1,  0, offs, dims)

       u_prim_RT_W  = velo(EDGE*idW +RT+1)*dom%len%elts(EDGE*idW+RT+1)
       u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
       u_prim_UP_S  = velo(EDGE*idS +UP+1)*dom%len%elts(EDGE*idS+UP+1)

       circ_LORT_SW = u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW

       vort(TRIAG*idSW+LORT+1) = circ_LORT_SW/dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idSW+UPLT+1) = vort(TRIAG*idSW+LORT+1)
    end if

    if (c == IPLUSJMINUS) then ! Parts 5, 6 of hexagon IPLUSJMINUS (SE corner of lozenge) combined to form pentagon
       id   = idx(PATCH_SIZE,    0, offs, dims)
       idS  = idx(PATCH_SIZE,   -1, offs, dims)
       idSW = idx(PATCH_SIZE-1, -1, offs, dims)

       u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT    = velo(EDGE*id+RT+1)  *dom%len%elts(EDGE*id+RT+1)

       circ_LORT_SW = - u_prim_RT + u_prim_RT_SW + u_prim_DG_SW 

       vort(TRIAG*idSW+LORT+1) = circ_LORT_SW/dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idS +UPLT+1) = vort(TRIAG*idSW+LORT+1)
    end if

    if (c == IMINUSJPLUS) then ! Parts 3, 4 of hexagon IMINUSJPLUS (NW corner of lozenge) combined to form pentagon
       id   = idx( 0, PATCH_SIZE,   offs, dims)
       idSW = idx(-1, PATCH_SIZE-1, offs, dims)
       idW  = idx(-1, PATCH_SIZE,   offs, dims)

       u_prim_UP    = velo(EDGE*id+UP+1)  *dom%len%elts(EDGE*id+UP+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_SW = velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)

       circ_UPLT_SW = u_prim_UP - u_prim_DG_SW - u_prim_UP_SW

       vort(TRIAG*idSW+UPLT+1) = circ_UPLT_SW/dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       vort(TRIAG*idW +LORT+1) = vort(TRIAG*idSW+UPLT+1)
    end if

    if (c == IJPLUS) then ! Parts 1, 2 of hexagon IJPLUS (NE corner of lozenge) combined to form pentagon
       id  = idx(PATCH_SIZE,   PATCH_SIZE,   offs, dims)
       idN = idx(PATCH_SIZE,   PATCH_SIZE+1, offs, dims)

       u_prim_RT   = velo(EDGE*id +RT+1)*dom%len%elts(EDGE*id+RT+1)
       u_prim_RT_N = velo(EDGE*idN+RT+1)*dom%len%elts(EDGE*idN+RT+1)
       u_prim_UP   = velo(EDGE*id +UP+1)*dom%len%elts(EDGE*id+UP+1)

       circ_LORT   = u_prim_RT - u_prim_RT_N - u_prim_UP

       vort(TRIAG*id+LORT+1) = circ_LORT/dom%triarea%elts(TRIAG*id+LORT+1)
       vort(TRIAG*id+UPLT+1) = vort(TRIAG*id+LORT+1)
    end if
  end subroutine post_vort

  subroutine interp_vel_hex (dom, i, j, zlev, offs, dims)
    ! Interpolate velocity to hexagon nodes in Cartesian coordinates; uses Perot formula as also used for kinetic energy
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    type(Coord) :: vel, x_e, x_i
    type(Coord) :: e_zonal, e_merid
    integer     :: id, idN, idE, idNE, idS, idSW, idW
    real(8)     :: lon, lat, u_dual_RT, u_dual_UP, u_dual_DG, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)

    ! Find the velocity on primal and dual grid edges, which are equal except for the length of the
    ! side they are on
    u_dual_RT    = velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_UP    = velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
    u_dual_DG    = velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_RT_W  = velo(EDGE*idW+RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_UP_S  = velo(EDGE*idS+UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1)

    x_i = dom%node%elts(id+1)

    ! Perot formula (15, 16) of Peixoto (2016) for velocity at hexagonal node from velocities at six adjacent edges
    x_e = dom%midpt%elts(EDGE*id+RT+1)
    vel = vec_scale(u_dual_RT, vec_minus(x_e,x_i))

    x_e = dom%midpt%elts(EDGE*idW+RT+1)
    vel = vec_plus(vel, vec_scale(u_dual_RT_W,  vec_minus(x_i,x_e)))

    x_e = dom%midpt%elts(EDGE*id+DG+1)
    vel = vec_plus(vel, vec_scale(u_dual_DG,    vec_minus(x_i,x_e)))

    x_e = dom%midpt%elts(EDGE*idSW+DG+1)
    vel = vec_plus(vel, vec_scale(u_dual_DG_SW, vec_minus(x_e,x_i)))

    x_e = dom%midpt%elts(EDGE*id+UP+1)
    vel = vec_plus(vel, vec_scale(u_dual_UP,    vec_minus(x_e,x_i)))

    x_e = dom%midpt%elts(EDGE*idS+UP+1)
    vel = vec_plus(vel, vec_scale(u_dual_UP_S,  vec_minus(x_i,x_e)))

    vel = vec_scale(dom%areas%elts(id+1)%hex_inv, vel) ! construct velocity at hexagonal node

    ! Project velocity onto zonal and meridional directions
    call cart2sph (x_i, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),             0.0_8) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    dom%u_zonal%elts(id+1) = inner(vel, e_zonal)
    dom%v_merid%elts(id+1) = inner(vel, e_merid)
  end subroutine interp_vel_hex

  subroutine cal_surf_press (q)
    implicit none
    ! Compute surface pressure
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q

    integer :: d, k, mass_type, p

    if (compressible) then
       mass_type = S_MASS
    else
       mass_type = S_TEMP
    end if

    do d = 1, size(grid)
       grid(d)%surf_press%elts = 0.0_8
       do k = 1, zlevels
          mass => q(mass_type,k)%data(d)%elts
          do p = 2, grid(d)%patch%length
             call apply_onescale_to_patch (column_mass, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass)
       end do
       grid(d)%surf_press%elts = grav_accel*grid(d)%surf_press%elts + p_top
    end do
  contains
    subroutine column_mass (dom, i, j, zlev, offs, dims)
      ! Sum up total mass over column id
      implicit none
      type (Domain)                  :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i

      id_i = idx(i, j, offs, dims) + 1

      dom%surf_press%elts(id_i) = dom%surf_press%elts(id_i) + mass(id_i)
    end subroutine column_mass
  end subroutine cal_surf_press

  subroutine integrate_pressure_up (dom, i, j, zlev, offs, dims)
    ! Integrate pressure (compressible case)/Lagrange multiplier (incompressible case) and geopotential up from surface to top layer
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: p_lower, p_upper

    id_i = idx (i, j, offs, dims) + 1

    if (compressible) then ! Compressible case
       if (zlev /= 1) then
          p_lower = dom%press_lower%elts(id_i)
       else 
          p_lower = dom%surf_press%elts(id_i)
       end if
       p_upper = p_lower - grav_accel*mass(id_i)
       dom%press%elts(id_i) = interp (p_lower, p_upper)
       dom%press_lower%elts(id_i) = p_upper

       ! Find geopotential at upper interface of current level using (18) in DYNAMICO
       if (zlev /= 1) then 
          dom%geopot_lower%elts(id_i) = dom%geopot%elts(id_i) 
       else
          dom%geopot_lower%elts(id_i) = surf_geopot (dom%node%elts(id_i))
       end if
       exner(id_i) = c_p * (dom%press%elts(id_i)/p_0)**kappa
       dom%geopot%elts(id_i) = dom%geopot_lower%elts(id_i) + grav_accel*kappa*temp(id_i)*exner(id_i)/dom%press%elts(id_i)
    else ! Incompressible case
       if (zlev == 1) then 
          dom%press%elts(id_i) = dom%surf_press%elts(id_i) - 0.5*grav_accel*temp(id_i)
       else ! Interpolate to lower interface of current level
          dom%press%elts(id_i) = dom%press%elts(id_i) - grav_accel*interp (dom%adj_temp%elts(id_i), temp(id_i))
       end if
       dom%adj_temp%elts(id_i) = temp(id_i)

       if (zlev == zlevels) then !top zlev, purely diagnostic
          if (abs(dom%press%elts(id_i)-0.5*grav_accel*temp(id_i) - p_top)> 1d-10) then
             print *, 'warning: upward integration of Lagrange multiplier not resulting in zero at top interface'
             write(6,'(A,es15.8)') "Pressure at infinity = ", dom%press%elts(id_i)-0.5*grav_accel*temp(id_i)
             write(6,'(A,es15.8)') "p_top = ", p_top
             write(6,'(A,es15.8)') "Difference = ", dom%press%elts(id_i)-0.5*grav_accel*temp(id_i) - p_top
             stop
          end if
       end if

       ! Find geopotential at interfaces using (18) in DYNAMICO
       ! Note: since mu is associated with the kinematic mass = inert mass (not the gravitational mass defined by the buyoancy)
       ! we divide by the constant reference density. This is the Boussinesq approximation.
       if (zlev == 1) then ! Save geopotential at lower interface of level zlev for interpolation in Bernoulli function
          dom%geopot_lower%elts(id_i) = surf_geopot (dom%node%elts(id_i))
       else
          dom%geopot_lower%elts(id_i) = dom%geopot%elts(id_i)
       end if
       dom%geopot%elts(id_i) = dom%geopot_lower%elts(id_i) + grav_accel*mass(id_i)/ref_density
    end if
  end subroutine integrate_pressure_up

  subroutine cal_pressure (dom, i, j, zlev, offs, dims)
    ! Integrate pressure up from surface to top layer
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: p_lower, p_upper

    id_i = idx(i, j, offs, dims)+1

    if (compressible) then ! Compressible case
        if (zlev /= 1) then
          p_lower = dom%press_lower%elts(id_i)
       else 
          p_lower = dom%surf_press%elts(id_i)
       end if
       p_upper = p_lower - grav_accel*mass(id_i)
       dom%press%elts(id_i) = interp (p_lower, p_upper)
       dom%press_lower%elts(id_i) = p_upper
    else ! Incompressible case
       if (zlev == 1) then 
          dom%press%elts(id_i) = dom%surf_press%elts(id_i) - 0.5*grav_accel*temp(id_i)
       else ! Interpolate to lower interface of current level
          dom%press%elts(id_i) = dom%press%elts(id_i) - grav_accel*interp (dom%adj_temp%elts(id_i), temp(id_i))
       end if
       dom%adj_temp%elts(id_i) = temp(id_i)
    end if
  end subroutine cal_pressure

  subroutine du_source (dom, i, j, zlev, offs, dims)
    ! Edge integrated source (non gradient) terms in velocity trend
    ! [Aechtner thesis page 56, Kevlahan, Dubos and Aechtner (2015)]
    implicit none
    type(Domain),                   intent(in) :: dom
    integer,                        intent(in) :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    interface
       function physics_velo_source (dom, i, j, zlev, offs, dims)
         use domain_mod
         real(8), dimension(1:EDGE)     :: physics_velo_source
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function physics_velo_source
    end interface

    integer                :: id, id_i
    real(8), dimension (3) :: Qperp_e, physics

    id = idx(i, j, offs, dims)
    id_i = id+1

    ! Calculate Q_perp
    Qperp_e = Qperp (dom, i, j, z_null, offs, dims)

    ! Calculate physics
    physics = physics_velo_source (dom, i, j, zlev, offs, dims)

    dvelo(EDGE*id+1:EDGE*id_i) = - Qperp_e + physics*dom%len%elts(EDGE*id+1:EDGE*id_i)
  end subroutine du_source

  function Qperp (dom, i, j, zlev, offs, dims)
    ! Compute energy-conserving edge integrated Qperp [Aechtner thesis page 44]
    implicit none
    real(8), dimension(3)          :: Qperp
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer               :: id, idNW, idN, idNE, idW, idE, idSW, idS, idSE
    real(8), dimension(5) :: wgt1, wgt2

    id   = idx(i, j, offs, dims)

    idNW = idx(i-1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idSE = idx(i+1, j-1, offs, dims)

    wgt1 = get_weights(dom, id,  0)
    wgt2 = get_weights(dom, idE, 3)

    Qperp(RT+1) = &
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

    Qperp(DG+1) = &
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

    Qperp(UP+1) = &
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

  function Qperp_Gassmann (dom, i, j, zlev, offs, dims)
    ! Compute energy-conserving edge integrated Qperp using Gassmann (2018) formula
    implicit none
    real(8), dimension(3)          :: Qperp_Gassmann
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer               :: id, idNW, idN, idNE, idW, idE, idSW, idS, idSE
    real(8), dimension(5) :: wgt1, wgt2

    id   = idx(i, j, offs, dims)

    idNW = idx(i-1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idSE = idx(i+1, j-1, offs, dims)

    ! RT edge
    wgt1 = get_weights(dom, id,  0)
    wgt2 = get_weights(dom, idE, 3)

    Qperp_Gassmann(RT+1) = &
         ! Adjacent neighbour edges (Gassmann rule 1)
         h_mflux(EDGE*id +DG+1) * qe(EDGE*idE+UP+1)*wgt1(1) + &
         h_mflux(EDGE*idS+UP+1) * qe(EDGE*idS+DG+1)*wgt1(5) + &
         h_mflux(EDGE*idS+DG+1) * qe(EDGE*idS+UP+1)*wgt2(1) + &
         h_mflux(EDGE*idE+UP+1) * qe(EDGE*id +DG+1)*wgt2(5) + &
         
         ! Second neighbour edges (Gassmann rule 2)
         h_mflux(EDGE*id  +UP+1) * qe(EDGE*id +DG+1)*wgt1(2) + &
         h_mflux(EDGE*idSW+DG+1) * qe(EDGE*idS+UP+1)*wgt1(4) + &
         h_mflux(EDGE*idSE+UP+1) * qe(EDGE*idS+DG+1)*wgt2(2) + &
         h_mflux(EDGE*idE +DG+1) * qe(EDGE*idE+UP+1)*wgt2(4)

    ! ! Third neighbour edges (Gassmann rule 3 = TRSK)
    ! h_mflux(EDGE*idW+RT+1)  * interp (qe(EDGE*idW+RT+1), qe(EDGE*id+RT+1))*wgt1(3) + &
    ! h_mflux(EDGE*idE+RT+1)  * interp (qe(EDGE*idE+RT+1), qe(EDGE*id+RT+1))*wgt2(3)

    if (dom%pedlen%elts(EDGE*idSW+DG+1)/=0.0_8) then ! Hexagon, third neighbour edge (Gassmann rule 3)
       Qperp_Gassmann(RT+1) = Qperp_Gassmann(RT+1) + h_mflux(EDGE*idW+RT+1)*interp (qe(EDGE*idW+RT+1),qe(EDGE*id+RT+1))*wgt1(3)
    else ! Pentagon, second neighbour edge (Gassmann rule 2)
       Qperp_Gassmann(RT+1) = Qperp_Gassmann(RT+1) + h_mflux(EDGE*idW+RT+1)*qe(EDGE*idS+UP+1)*wgt1(3)
    end if

    if (dom%pedlen%elts(EDGE*idSE+UP+1)/=0.0_8) then ! Hexagon, third neighbour edge (Gassmann rule 3)
       Qperp_Gassmann(RT+1) = Qperp_Gassmann(RT+1) + h_mflux(EDGE*idE+RT+1)*interp (qe(EDGE*idE+RT+1),qe(EDGE*id+RT+1))*wgt2(3)
    else ! Pentagon, second neighbour edge (Gassmann rule 2)
       Qperp_Gassmann(RT+1) = Qperp_Gassmann(RT+1) + h_mflux(EDGE*idE+RT+1)*qe(EDGE*idS+DG+1)*wgt2(3)
    end if

    ! DG edge - no modification required for pentagons since TRSK rule edges have zero flux in pentagon case
    wgt1 = get_weights(dom, id,   1)
    wgt2 = get_weights(dom, idNE, 4)

    Qperp_Gassmann(DG+1) = &
                                ! Adjacent neighbour edges (Gassmann rule 1)
         h_mflux(EDGE*id +UP+1) * qe(EDGE*idN+RT+1)*wgt1(1) + &
         h_mflux(EDGE*id +RT+1) * qe(EDGE*idE+UP+1)*wgt1(5) + &
         h_mflux(EDGE*idE+UP+1) * qe(EDGE*id +RT+1)*wgt2(1) + &
         h_mflux(EDGE*idN+RT+1) * qe(EDGE*id +UP+1)*wgt2(5) + &
         
                                ! Second neighbour edges (Gassmann rule 2)
         h_mflux(EDGE*idW +RT+1) * qe(EDGE*id +UP+1)*wgt1(2) + &
         h_mflux(EDGE*idS +UP+1) * qe(EDGE*id +RT+1)*wgt1(4) + &
         h_mflux(EDGE*idNE+RT+1) * qe(EDGE*idE+UP+1)*wgt2(2) + &
         h_mflux(EDGE*idNE+UP+1) * qe(EDGE*idN+RT+1)*wgt2(4) + &
         
                                ! Third neighbour edges (Gassmann rule 3 = TRSK)
         h_mflux(EDGE*idSW+DG+1) * interp (qe(EDGE*idSW+DG+1), qe(EDGE*id+DG+1))*wgt1(3) + &
         h_mflux(EDGE*idNE+DG+1) * interp (qe(EDGE*idNE+DG+1), qe(EDGE*id+DG+1))*wgt2(3)

    ! UP edge
    wgt1 = get_weights(dom, id,  2)
    wgt2 = get_weights(dom, idN, 5)

    Qperp_Gassmann(UP+1) = &
         ! Adjacent neighbour edges (Gassmann rule 1)
         h_mflux(EDGE*idW+RT+1)  * qe(EDGE*idW+DG+1)*wgt1(1) + &
         h_mflux(EDGE*id +DG+1)  * qe(EDGE*idN+RT+1)*wgt1(5) + &
         h_mflux(EDGE*idN+RT+1)  * qe(EDGE*id +DG+1)*wgt2(1) + &
         h_mflux(EDGE*idW+DG+1)  * qe(EDGE*idW+RT+1)*wgt2(5) + &
         
         ! Second neighbour edges (Gassmann rule 2)
         h_mflux(EDGE*idSW+DG+1) * qe(EDGE*idW+RT+1)*wgt1(2) + &
         h_mflux(EDGE*id+RT+1)   * qe(EDGE*id +DG+1)*wgt1(4) + &
         h_mflux(EDGE*idN+DG+1)  * qe(EDGE*idN+RT+1)*wgt2(2) + &         
         h_mflux(EDGE*idNW+RT+1) * qe(EDGE*idW+DG+1)*wgt2(4)

    ! ! Third neighbour edges (Gassmann rule 3 = TRSK)
    ! h_mflux(EDGE*idS+UP+1)  * interp (qe(EDGE*idS+UP+1),  qe(EDGE*id+UP+1))*wgt1(3) + &
    ! h_mflux(EDGE*idN+UP+1)  * interp (qe(EDGE*idN+UP+1),  qe(EDGE*id+UP+1))*wgt2(3)

    if (dom%pedlen%elts(EDGE*idSW+DG+1)/=0.0_8) then ! Hexagon, third neighbour edge (Gassmann rule 3 = TRSK)
       Qperp_Gassmann(UP+1) = Qperp_Gassmann(UP+1) + h_mflux(EDGE*idS+UP+1)*interp (qe(EDGE*idS+UP+1),qe(EDGE*id+UP+1))*wgt1(3)
    else ! Pentagon, second neighbour edge (Gassmann rule 2)
       Qperp_Gassmann(UP+1) = Qperp_Gassmann(UP+1) + h_mflux(EDGE*idS+UP+1)*qe(EDGE*idW+RT+1)*wgt1(3)
    end if

    if (dom%pedlen%elts(EDGE*idNW+RT+1)/=0.0_8) then ! Hexagon, third neighbour edge (Gassmann rule 3)
       Qperp_Gassmann(UP+1) = Qperp_Gassmann(UP+1) + h_mflux(EDGE*idN+UP+1)*interp (qe(EDGE*idN+UP+1),qe(EDGE*id+UP+1))*wgt2(3)
    else ! Pentagon, second neighbour edge (Gassmann rule 2)
       Qperp_Gassmann(UP+1) = Qperp_Gassmann(UP+1) + h_mflux(EDGE*idN+UP+1)*qe(EDGE*idW+DG+1)*wgt2(3)
    end if
  end function Qperp_Gassmann

  function get_weights(dom, id, offs)
    ! Weights for Qperp computation [Aechtner thesis page 44]
    implicit none
    type(Domain) :: dom
    integer      :: id, offs

    integer               :: i
    real(8), dimension(5) :: get_weights, wgt

    wgt(1) = dom%areas%elts(id+1)%part(1+offs)
    do i = 2, 5
       wgt(i) = wgt(i-1) + dom%areas%elts(id+1)%part(modulo(i+offs-1,6)+1)
    end do
    wgt = 0.5_8 - wgt*dom%areas%elts(id+1)%hex_inv
    get_weights = (/wgt(1), -wgt(2), wgt(3), -wgt(4), wgt(5)/)
  end function get_weights

  subroutine scalar_trend (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                           :: id_i
    real(8), dimension(S_MASS:S_TEMP) :: physics

    interface
       function physics_scalar_source (dom, i, j, zlev, offs, dims)
         use domain_mod
         real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
         type(Domain)                      :: dom
         integer                           :: i, j, zlev
         integer, dimension(N_BDRY+1)      :: offs
         integer, dimension(2,N_BDRY+1)    :: dims
       end function physics_scalar_source
    end interface

    id_i = idx(i, j, offs, dims) + 1

    physics = physics_scalar_source (dom, i, j, zlev, offs, dims)

    dmass(id_i) = - div (h_mflux, dom, i, j, offs, dims) + physics(S_MASS)
    dtemp(id_i) = - div (h_tflux, dom, i, j, offs, dims) + physics(S_TEMP)
  end subroutine scalar_trend

  subroutine cal_Laplacian_scalar (dom, i, j, zlev, offs, dims)
    ! Calculate divergence of gradient of scalar
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idW, idSW, idS
    integer :: id_i, idE_i, idNE_i, idN_i, idW_i, idSW_i, idS_i
    
    id   = idx(i, j, offs, dims)
    id_i = id+1
    
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)

    idE_i  = idE+1
    idNE_i = idNE+1
    idN_i  = idN+1
    idW_i  = idW+1
    idSW_i = idSW+1
    idS_i  = idS+1

    Laplacian(id_i) = div_grad (grad_flux())
  contains
    function grad_flux()
      ! Calculates gradient flux
      implicit none
      real(8), dimension(6) :: grad_flux

      grad_flux(1) =  (sclr(idE_i) - sclr(id_i))  /dom%len%elts(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
      grad_flux(2) =  (sclr(id_i)  - sclr(idNE_i))/dom%len%elts(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
      grad_flux(3) =  (sclr(idN_i) - sclr(id_i))  /dom%len%elts(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)
      grad_flux(4) = -(sclr(idW_i) - sclr(id_i))  /dom%len%elts(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
      grad_flux(5) = -(sclr(id_i)  - sclr(idSW_i))/dom%len%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
      grad_flux(6) = -(sclr(idS_i) - sclr(id_i))  /dom%len%elts(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)
    end function grad_flux

    real(8) function div_grad (grad)
      ! Calculates divergence of gradient flux
      implicit none
      real(8), dimension(6) :: grad

      div_grad = (grad(1)-grad(4) + grad(5)-grad(2) + grad(3)-grad(6)) * dom%areas%elts(id_i)%hex_inv
    end function div_grad
  end subroutine cal_Laplacian_scalar

  subroutine cal_Laplacian_rotu (dom, i, j, zlev, offs, dims)
    ! Curl of vorticity given at triangle circumcentres x_v, i.e. rotational part of vector Laplacian
    ! output is at edges x_e
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW

    id   = idx(i,   j,   offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)

    Laplacian(EDGE*id+RT+1) = -(vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1))/dom%pedlen%elts(EDGE*id+RT+1)
    Laplacian(EDGE*id+DG+1) = -(vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+DG+1)
    Laplacian(EDGE*id+UP+1) = -(vort(TRIAG*idW+LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+UP+1)
  end subroutine cal_Laplacian_rotu

  function gradi_e (scalar, dom, i, j, offs, dims)
    ! Gradient of a scalar at nodes x_i
    ! output is at edges
    ! If type = .true. then compute the gradient at the southwest edges of the hexagon
    implicit none
    real(8), dimension(3)          :: gradi_e
    real(8), dimension(:), pointer :: scalar
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idN, idNE

    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    gradi_e(RT+1) = (scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*id+RT+1)
    gradi_e(DG+1) = (scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*id+DG+1)
    gradi_e(UP+1) = (scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*id+UP+1)
  end function gradi_e

  function curlv_e (curl, dom, i, j, offs, dims)
    ! Curl of vorticity given at triangle circumcentres x_v, rot(rot(u))
    ! output is at edges x_e
    implicit none
    real(8), dimension(:), pointer :: curl
    real(8), dimension(3)          :: curlv_e
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW

    id   = idx(i,   j,   offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)

    curlv_e(RT+1) = (curl(TRIAG*id +LORT+1) - curl(TRIAG*idS+UPLT+1))/dom%pedlen%elts(EDGE*id+RT+1)
    curlv_e(DG+1) = (curl(TRIAG*id +LORT+1) - curl(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+DG+1)
    curlv_e(UP+1) = (curl(TRIAG*idW+LORT+1) - curl(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+UP+1)
  end function curlv_e

  function div (hflux, dom, i, j, offs, dims)
    ! Divergence at nodes x_i given horizontal fluxes at edges x_e
    implicit none
    real(8)                         :: div
    real(8), dimension(:), pointer  :: hflux
    type(Domain)                    :: dom
    integer                         :: i, j
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: id, idW, idS, idSW

    id   = idx(i,   j,   offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)

    div = (hflux(EDGE*id+RT+1)-hflux(EDGE*idW+RT+1) + hflux(EDGE*idSW+DG+1)-hflux(EDGE*id+DG+1) &
         + hflux(EDGE*id+UP+1)-hflux(EDGE*idS+UP+1)) * dom%areas%elts(id+1)%hex_inv
  end function div

  real(8) function interp (e1, e2)
    ! Centred average interpolation of quantities e1 and e2
    implicit none
    real(8) :: e1, e2

    interp = 0.5 * (e1 + e2)
  end function interp

  subroutine du_grad (dom, i, j, zlev, offs, dims)
    ! Add gradients of Bernoulli and Exner to dvelo [DYNAMICO (23)-(25)]
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer                      :: id, id_i, idE_i, idN_i, idNE_i
    real(8), dimension(3)        :: gradB, gradE, theta_e
    real(8), dimension(0:N_BDRY) :: theta

    id = idx(i, j, offs, dims)
    id_i = id+1

    idE_i  = idx(i+1, j,   offs, dims)+1
    idN_i  = idx(i,   j+1, offs, dims)+1
    idNE_i = idx(i+1, j+1, offs, dims)+1

    ! See DYNAMICO between (23)-(25), geopotential still known from step1_upw
    ! the theta multiplying the Exner gradient is the edge-averaged non-mass-weighted potential temperature
    theta(0)         = temp(id_i)/mass(id_i)
    theta(NORTH)     = temp(idN_i)/mass(idN_i)
    theta(EAST)      = temp(idE_i)/mass(idE_i)
    theta(NORTHEAST) = temp(idNE_i)/mass(idNE_i)

    ! Interpolate potential temperature to edges
    if (compressible) then
       theta_e(1) = interp (theta(0), theta(EAST))
       theta_e(2) = interp (theta(0), theta(NORTHEAST))
       theta_e(3) = interp (theta(0), theta(NORTH))
    else
       theta_e(1) = interp (1.0_8-theta(0), 1.0_8-theta(EAST))
       theta_e(2) = interp (1.0_8-theta(0), 1.0_8-theta(NORTHEAST)) 
       theta_e(3) = interp (1.0_8-theta(0), 1.0_8-theta(NORTH)) 
    end if

    ! Calculate gradients
    gradB = gradi_e (bernoulli, dom, i, j, offs, dims)
    gradE = gradi_e (exner,     dom, i, j, offs, dims)

    ! Update velocity trend (source dvelo calculated was edge integrated)
    dvelo(EDGE*id+1:EDGE*id_i) = dvelo(EDGE*id+1:EDGE*id_i)/dom%len%elts(EDGE*id+1:EDGE*id_i) - gradB - theta_e*gradE
  end subroutine du_grad

  subroutine cal_divu (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i, idS, idW, idSW
    real(8) :: u_dual_RT, u_dual_RT_W, u_dual_DG_SW, u_dual_DG, u_dual_UP, u_dual_UP_S

    id   = idx(i,   j,   offs, dims)
    id_i = id+1
    
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)

    u_dual_RT    = velo(EDGE*id  +RT+1)*dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_RT_W  = velo(EDGE*idW +RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_DG    = velo(EDGE*id  +DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    = velo(EDGE*id  +UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
    u_dual_UP_S  = velo(EDGE*idS +UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)

    divu(id_i) =  (u_dual_RT-u_dual_RT_W + u_dual_DG_SW-u_dual_DG + u_dual_UP-u_dual_UP_S) * dom%areas%elts(id_i)%hex_inv
  end subroutine cal_divu

  subroutine cal_vort (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idN
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_UP_E, u_prim_RT_N

    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)

    u_prim_RT   = velo(EDGE*id+RT +1)*dom%len%elts(EDGE*id+RT+1)
    u_prim_DG   = velo(EDGE*id+DG +1)*dom%len%elts(EDGE*id+DG+1)
    u_prim_UP   = velo(EDGE*id+UP +1)*dom%len%elts(EDGE*id+UP+1)
    u_prim_UP_E = velo(EDGE*idE+UP+1)*dom%len%elts(EDGE*idE+UP+1)
    u_prim_RT_N = velo(EDGE*idN+RT+1)*dom%len%elts(EDGE*idN+RT+1)

    vort(TRIAG*id+LORT+1) =   (u_prim_RT + u_prim_UP_E + u_prim_DG)  /dom%triarea%elts(TRIAG*id+LORT+1) 
    vort(TRIAG*id+UPLT+1) = - (u_prim_DG + u_prim_UP   + u_prim_RT_N)/dom%triarea%elts(TRIAG*id+UPLT+1)
  end subroutine cal_vort

  subroutine sum_mass_temp (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id   = idx(i, j, offs, dims)

    sum_mass = sum_mass + abs(mass(id+1))
    sum_temp = sum_temp + abs(temp(id+1))
  end subroutine sum_mass_temp

  subroutine sum_dmassdtemp (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer :: id_i

    id_i = idx(i, j, offs, dims)+1

    totaldmass = totaldmass + dmass(id_i)/dom%areas%elts(id_i)%hex_inv
    totalabsdmass = totalabsdmass + abs (dmass(id_i)/dom%areas%elts(id_i)%hex_inv)

    totaldtemp = totaldtemp + dtemp(id_i)/dom%areas%elts(id_i)%hex_inv
    totalabsdtemp = totalabsdtemp + abs (dtemp(id_i)/dom%areas%elts(id_i)%hex_inv)
  end subroutine sum_dmassdtemp

  subroutine comp_offs3 (dom, p, offs, dims)
    implicit none
    type(Domain)                 :: dom
    integer                      :: p 
    integer, dimension(0:N_BDRY) :: offs
    integer, dimension(2,N_BDRY) :: dims

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

  subroutine evals_diffusion
    ! Estimates largest eigenvalues of Laplacian diffusion operators by power iteration
    use wavelet_mod
    implicit none

    integer                                                       :: d, iter, j, k
    integer, parameter                                            :: iter_max = 100, seed = 86456
    real(8), dimension(3)                                         :: err, eval, eval_old, inner_prod
    real(8), dimension(S_MASS:S_VELO,1:zlevels)                   :: lnorm
    real(8), parameter                                            :: err_max = 1d-4

    character(3), parameter :: order = "2"

    if (rank == 0) write (6,'(/,A)') "Finding diffusion lengthscales by power iteration:"
    
    ! Initialize variables to random values and normalize
    call init_rand

    iter = 1; err = 1d16; eval = 0.0_8
    do while (maxval (err) > err_max .and. iter <= iter_max)
       ! Apply Laplacian operators
       call Ax

       ! Normalize new eigenvectors
       call apply_normalize

       ! Find Rayleigh quotient approximation to largest eigenvalues
       eval_old = eval
       call Ray_quotient

       ! Find relative change in eigenvalues
       err = abs (eval - eval_old) / abs (eval)
       
       iter = iter + 1
    end do

    if (rank == 0) then
       if (iter > iter_max) then
          write (6,'(2(A,es8.2),A,i3,A)') "Warning: eigenvalue error ",  &
               maxval (err), " not converged to specified error ", err_max, " after ", iter_max, " iterations"
       else
          write (6,'(A,es8.2,A,i3,A)') "Eigenvalues converged to relative error ", maxval (err), " after ", iter-1, " iterations"
       end if
    end if

    ! Find diffusion length scales
    L_diffusion(1) = 1/sqrt(-eval(1))
    L_diffusion(2:3) = 1/sqrt(-eval(2:3))
    if (rank == 0) write (6,'(3(A,es8.2,1x),/)') &
         "dx_scalar = ", MATH_PI*L_diffusion(1), "dx_divu = ",MATH_PI*L_diffusion(2),"dx_rotu = ", MATH_PI*L_diffusion(3)
  contains
    subroutine init_rand
      ! Applies random initial conditions
      implicit none
      integer :: i, m
      
      call random_seed (size=m)
      call random_seed (put=(/(i,i=1,m)/))

      call apply_onescale (init_rand_scalar, level_start, z_null, 0, 1)
      call apply_onescale (init_rand_velo,   level_start, z_null, 0, 0)

      call update_evec
      call apply_normalize
    end subroutine init_rand

    subroutine init_rand_scalar (dom, i, j, zlev, offs, dims)
      ! Initializes mass to a random number on [-1,1) 
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      integer :: d, id_i
      real(8) :: harvest

      d  = dom%id+1
      id_i = idx(i, j, offs, dims)+1

      call random_number (harvest)
      trend(S_MASS,1)%data(d)%elts(id_i) = harvest - 0.5_8
    end subroutine init_rand_scalar

    subroutine init_rand_velo (dom, i, j, zlev, offs, dims)
      ! Initializes velocities to a uniform random number on [-1, 1)
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      integer :: d, e, id, id_e
      real(8) :: harvest

      d  = dom%id+1
      id = idx(i, j, offs, dims)

      do e = 1, EDGE
         id_e = EDGE*id+e
         call random_number (harvest)
         trend(S_VELO,1)%data(d)%elts(id_e) = harvest - 0.5_8
         
         call random_number (harvest)
         trend(S_VELO,2)%data(d)%elts(id_e) = harvest - 0.5_8
      end do
    end subroutine init_rand_velo

    subroutine apply_normalize
      ! L2 normalization of eigenvectors
      implicit none

      ! Find norms
      call cal_lnorm (trend, order, lnorm)

      call apply_onescale (normalize_scalar, level_start, z_null, 0, 1)
      call apply_onescale (normalize_velo,   level_start, z_null, 0, 0)

      call update_evec
    end subroutine apply_normalize

    subroutine normalize_scalar (dom, i, j, zlev, offs, dims)
      ! Normalizes the mass eigenvector
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      integer :: d, id_i, v

      d  = dom%id+1
      id_i = idx(i, j, offs, dims)+1

      trend(S_MASS,1)%data(d)%elts(id_i) = trend(S_MASS,1)%data(d)%elts(id_i)/lnorm(S_MASS,1)
    end subroutine normalize_scalar

    subroutine normalize_velo (dom, i, j, zlev, offs, dims)
      ! Normalizes the velocity eigenvectors
      implicit none
      type (Domain)                   :: dom
      integer                         :: i, j, zlev
      integer, dimension (N_BDRY+1)   :: offs
      integer, dimension (2,N_BDRY+1) :: dims

      integer :: d, e, id, id_e

      d  = dom%id+1
      id = idx(i, j, offs, dims)
      do e = 1, EDGE
         id_e = EDGE*id+e
         trend(S_VELO,1)%data(d)%elts(id_e) = trend(S_VELO,1)%data(d)%elts(id_e)/lnorm(S_VELO,1)
         trend(S_VELO,2)%data(d)%elts(id_e) = trend(S_VELO,2)%data(d)%elts(id_e)/lnorm(S_VELO,2)
      end do
    end subroutine normalize_velo

    subroutine Ax
      ! Applies Laplacian operators to previous eigenvectors to find new eigenvectors
      implicit none
      integer :: d, j

      ! Find div(u), rot(u)
      do d = 1, size(grid) 
         velo => trend(S_VELO,1)%data(d)%elts
         divu => grid(d)%divu%elts
         do j = 1, grid(d)%lev(level_start)%length
            call apply_onescale_to_patch (cal_divu, grid(d), grid(d)%lev(level_start)%elts(j), z_null, 0, 1)
         end do
         nullify (velo, divu)

         velo => trend(S_VELO,2)%data(d)%elts
         vort => grid(d)%vort%elts
         do j = 1, grid(d)%lev(level_start)%length
            call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(level_start)%elts(j), z_null, -1, 0)
         end do
         call apply_to_penta_d (post_vort, grid(d), level_start, z_null)
         nullify (velo, vort)
      end do

      ! Find Laplacians
      do d = 1, size(grid)
         sclr => trend(S_MASS,1)%data(d)%elts
         velo => trend(S_VELO,1)%data(d)%elts
         divu => grid(d)%divu%elts
         Laplacian => Laplacian_scalar(S_MASS)%data(d)%elts
         do j = 1, grid(d)%lev(level_start)%length
            call apply_onescale_to_patch (cal_Laplacian_scalar, grid(d), grid(d)%lev(level_start)%elts(j), z_null, 0, 1)
            call apply_onescale_to_patch (cal_Laplacian_divu,   grid(d), grid(d)%lev(level_start)%elts(j), z_null, 0, 0)
         end do
         nullify (sclr, velo, divu, Laplacian)
         
         velo => trend(S_VELO,2)%data(d)%elts
         vort => grid(d)%vort%elts
         do j = 1, grid(d)%lev(level_start)%length
            call apply_onescale_to_patch (cal_Laplacian_rotu, grid(d), grid(d)%lev(level_start)%elts(j), z_null, 0, 0)
         end do
         nullify (velo, vort)
      end do

      ! Update scalar eigenvector
      trend(S_MASS,1) = Laplacian_scalar(S_MASS)
      
      ! Update boundaries
      call update_evec
    end subroutine Ax

    subroutine cal_Laplacian_divu (dom, i, j, zlev, offs, dims)
      ! Calculate divergence part of Laplacian
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer               :: id, id_i
      real(8), dimension(3) :: grad_divu

      id = idx(i, j, offs, dims)
      id_i = id+1
      
      grad_divu = gradi_e (divu, dom, i, j, offs, dims) 

      velo(EDGE*id+1:EDGE*id_i) = grad_divu
    end subroutine cal_Laplacian_divu

    subroutine cal_Laplacian_rotu (dom, i, j, zlev, offs, dims)
      ! Calculate rotational part of Laplacian
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer               :: id, id_i
      real(8), dimension(3) :: curl_rotu

      curl_rotu = curlv_e (vort, dom, i, j, offs, dims)

      id = idx(i, j, offs, dims)
      id_i = id+1

      velo(EDGE*id+1:EDGE*id_i) = - curl_rotu
    end subroutine cal_Laplacian_rotu

    subroutine Ray_quotient
      implicit none
      integer :: ii

      ! Save current eigenvectors, x (already normalized)
      trend(S_MASS,2) = trend(S_MASS,1)
      trend(S_VELO,3) = trend(S_VELO,1)
      trend(S_VELO,4) = trend(S_VELO,2)

      ! Apply Laplacian operators to find Ax
      call Ax

      ! Rayleigh quotient approximation to eigenvalue <Ax,x>/<x,x> (x is already normalized)
      inner_prod = 0.0_8
      call apply_onescale (cal_inner_prod_scalar, level_start, z_null, 0, 1)
      call apply_onescale (cal_inner_prod_velo,   level_start, z_null, 0, 0)
      do ii = 1, 3
         eval(ii) = sum_real (inner_prod(ii))
      end do
    end subroutine Ray_quotient

    subroutine cal_inner_prod_scalar (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id_i, v

      d = dom%id+1
      id_i = idx(i, j, offs, dims)+1

      inner_prod = inner_prod + trend(S_MASS,2)%data(d)%elts(id_i)*trend(S_MASS,1)%data(d)%elts(id_i)
    end subroutine cal_inner_prod_scalar

    subroutine cal_inner_prod_velo (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, e, id, id_e

      d = dom%id+1
      id = idx(i, j, offs, dims)

      do e = 1, EDGE
         id_e  = EDGE*id+e
         inner_prod(2) = inner_prod(2) + trend(S_VELO,3)%data(d)%elts(id_e)*trend(S_VELO,1)%data(d)%elts(id_e)
         inner_prod(3) = inner_prod(3) + trend(S_VELO,4)%data(d)%elts(id_e)*trend(S_VELO,2)%data(d)%elts(id_e)
      end do
    end subroutine cal_inner_prod_velo

    subroutine update_evec
      implicit none

      trend(S_MASS,1)%bdry_uptodate = .false.
      trend(S_VELO,1:2)%bdry_uptodate = .false.
      call update_bdry (trend(S_MASS,1), level_start)
      call update_vector_bdry (trend(S_VELO,1:2), level_start)
    end subroutine update_evec
  end subroutine evals_diffusion
end module ops_mod
