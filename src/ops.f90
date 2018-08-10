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
      
      u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
      u_prim_UP_SW = velo(EDGE*idSW+UP+1)*dom%len%elts(EDGE*idSW+UP+1)

      circ_LORT_SW =   u_prim_RT_SW + u_prim_UP_S  + u_prim_DG_SW
      circ_UPLT_SW = -(u_prim_RT_W  + u_prim_DG_SW + u_prim_UP_SW)

      vort(TRIAG*idW+LORT+1) = circ_LORT_W/dom%triarea%elts(TRIAG*idW+LORT+1) 
      vort(TRIAG*idS+UPLT+1) = circ_UPLT_S/dom%triarea%elts(TRIAG*idS+UPLT+1) 

      pv_LORT_SW = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_LORT_SW)/( &
           mass(idSW+1)*dom%areas%elts(idSW+1)%part(1) + &
           mass(idS+1)*dom%areas%elts(idS+1)%part(3) + &
           mass(id+1)*dom%areas%elts(id+1)%part(5))

      pv_UPLT_SW = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_UPLT_SW)/( &
           mass(idSW+1)*dom%areas%elts(idSW+1)%part(2) + &
           mass(id+1)*dom%areas%elts(id+1)%part(4) + &
           mass(idW+1)*dom%areas%elts(idW+1)%part(6))

      qe(EDGE*idW+RT+1)  = interp (pv_LORT_W , pv_UPLT_SW)
      qe(EDGE*idSW+DG+1) = interp (pv_LORT_SW, pv_UPLT_SW)
      qe(EDGE*idS+UP+1)  = interp (pv_LORT_SW, pv_UPLT_S)

      ! Mass and temperature fluxes
      physics = physics_scalar_flux (dom, id, idW, idSW, idS, .true.)

      h_mflux(EDGE*idW+RT+1)  = u_dual_RT_W  * interp (mass(id+1), mass(idW+1))  + physics(S_MASS,RT+1)
      h_mflux(EDGE*idSW+DG+1) = u_dual_DG_SW * interp (mass(id+1), mass(idSW+1)) + physics(S_MASS,DG+1)
      h_mflux(EDGE*idS+UP+1)  = u_dual_UP_S  * interp (mass(id+1), mass(idS+1))  + physics(S_MASS,UP+1)

      h_tflux(EDGE*idW+RT+1)  = u_dual_RT_W  * interp (temp(id+1), temp(idW+1))  + physics(S_TEMP,RT+1)
      h_tflux(EDGE*idSW+DG+1) = u_dual_DG_SW * interp (temp(id+1), temp(idSW+1)) + physics(S_TEMP,DG+1)
      h_tflux(EDGE*idS+UP+1)  = u_dual_UP_S  * interp (temp(id+1), temp(idS+1))  + physics(S_TEMP,UP+1)
    end subroutine comp_ijmin

    subroutine comput
      ! Computes physical quantities during upward integration
      implicit none
      type (Coord)                             :: x_e, x_i, vel
      integer                                  :: idE, idN, idNE, idS, idSW, idW
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
           u_prim_UP_S*u_dual_UP_S + u_prim_DG_SW*u_dual_DG_SW + u_prim_RT_W*u_dual_RT_W) * dom%areas%elts(id+1)%hex_inv/4.0_8

      ! Interpolate geopotential from interfaces to level
      Phi_k = interp (dom%geopot%elts(id+1), dom%adj_geopot%elts(id+1))

      ! Bernoulli function
      if (compressible) then 
         bernoulli(id+1) = kinetic_energy + Phi_k
      else 
         bernoulli(id+1) = kinetic_energy + Phi_k + dom%press%elts(id+1)/ref_density
      end if

      ! Exner function in incompressible case from geopotential
      if (.not. compressible) exner(id+1) = -Phi_k

      ! Calculate div(u) for velocity diffusion
      divu(id+1) = (u_dual_RT-u_dual_RT_W + u_dual_DG_SW-u_dual_DG + u_dual_UP-u_dual_UP_S) * dom%areas%elts(id+1)%hex_inv 

      circ_LORT   =   u_prim_RT    + u_prim_UP_E + u_prim_DG 
      circ_UPLT   = -(u_prim_DG    + u_prim_UP   + u_prim_RT_N)
      circ_LORT_W =   u_prim_RT_W  + u_prim_UP   + u_prim_DG_W
      circ_UPLT_S = -(u_prim_RT    + u_prim_DG_S + u_prim_UP_S)

      vort(TRIAG*id+LORT+1) = circ_LORT/dom%triarea%elts(TRIAG*id+LORT+1) 
      vort(TRIAG*id+UPLT+1) = circ_UPLT/dom%triarea%elts(TRIAG*id+UPLT+1)

      pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT)/( &
           mass(id+1)*dom%areas%elts(id+1)%part(1) + &
           mass(idE+1)*dom%areas%elts(idE+1)%part(3) + &
           mass(idNE+1)*dom%areas%elts(idNE+1)%part(5))

      pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + circ_UPLT)/( &
           mass(id+1)*dom%areas%elts(id+1)%part(2) + &
           mass(idNE+1)*dom%areas%elts(idNE+1)%part(4) + &
           mass(idN+1)*dom%areas%elts(idN+1)%part(6))

      pv_LORT_W = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_LORT_W)/( &
           mass(idW+1)*dom%areas%elts(idW+1)%part(1) + &
           mass(id+1)*dom%areas%elts(id+1)%part(3) + &
           mass(idN+1)*dom%areas%elts(idN+1)%part(5))

      pv_UPLT_S = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_UPLT_S)/( &
           mass(idS+1)*dom%areas%elts(idS+1)%part(2) + &
           mass(idE+1)*dom%areas%elts(idE+1)%part(4) + &
           mass(id+1)*dom%areas%elts(id+1)%part(6))

      qe(EDGE*id+RT+1) = interp (pv_UPLT_S, pv_LORT)
      qe(EDGE*id+DG+1) = interp (pv_UPLT,   pv_LORT)
      qe(EDGE*id+UP+1) = interp (pv_UPLT,   pv_LORT_W)

      ! Mass and temperature fluxes
      physics = physics_scalar_flux (dom, id, idE, idNE, idN)

      h_mflux(EDGE*id+RT+1) = u_dual_RT * interp (mass(id+1), mass(idE+1))  + physics(S_MASS,RT+1)
      h_mflux(EDGE*id+DG+1) = u_dual_DG * interp (mass(id+1), mass(idNE+1)) + physics(S_MASS,DG+1)
      h_mflux(EDGE*id+UP+1) = u_dual_UP * interp (mass(id+1), mass(idN+1))  + physics(S_MASS,UP+1)

      h_tflux(EDGE*id+RT+1) = u_dual_RT * interp (temp(id+1), temp(idE+1))  + physics(S_TEMP,RT+1)
      h_tflux(EDGE*id+DG+1) = u_dual_DG * interp (temp(id+1), temp(idNE+1)) + physics(S_TEMP,DG+1)
      h_tflux(EDGE*id+UP+1) = u_dual_UP * interp (temp(id+1), temp(idN+1))  + physics(S_TEMP,UP+1)
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

       pv_LORT_SW = (dom%coriolis%elts(TRIAG*idSW+1) + circ_LORT_SW)/ &
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

       pv_LORT = (dom%coriolis%elts(TRIAG*id+1) + circ_LORT)/ &          
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
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (column_mass, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass)
       end do
       grid(d)%surf_press%elts = grav_accel*grid(d)%surf_press%elts + press_infty
    end do
  end subroutine cal_surf_press

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

  subroutine integrate_pressure_up (dom, i, j, zlev, offs, dims)
    ! Integrate pressure (compressible case)/Lagrange multiplier (incompressible case) and geopotential up from surface to top layer
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: err, spec_vol

    id = idx(i, j, offs, dims)

    if (compressible) then ! Compressible case
       if (zlev == 1) then
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*mass(id+1)
       else ! Interpolate mass=rho*dz to lower interface of current level
          dom%press%elts(id+1) = dom%press%elts(id+1) - grav_accel*interp (mass(id+1), dom%adj_mass%elts(id+1))
       end if
       dom%adj_mass%elts(id+1) = mass(id+1) ! Save current mass for pressure calculation at next vertical level

       if (zlev == zlevels) then !top zlev, purely diagnostic
          err = abs((dom%press%elts(id+1) - 0.5_8*grav_accel*mass(id+1)) - press_infty)/dom%surf_press%elts(id+1)
          if (err > 1d-10) then
             write(6,*) 'Warning: upward integration of pressure not resulting in zero at top interface'
             write(6,*) '(observed pressure - pressure_infty)/P_Surface =', err
             stop
          end if
       end if

       ! Exner function from pressure
       exner(id+1) = c_p*(dom%press%elts(id+1)/ref_press)**kappa

       ! Specific volume alpha = kappa*theta*pi/p
       spec_vol = kappa * temp(id+1)/mass(id+1) * exner(id+1) / dom%press%elts(id+1)

       ! Find geopotential at upper interface of current level using (18) in DYNAMICO
       if (zlev == 1) then ! Save geopotential at lower interface of level zlev for interpolation in Bernoulli function
          dom%adj_geopot%elts(id+1) = surf_geopot (dom%node%elts(id+1))
       else
          dom%adj_geopot%elts(id+1) = dom%geopot%elts(id+1)
       end if
       dom%geopot%elts(id+1) = dom%adj_geopot%elts(id+1) + grav_accel*mass(id+1)*spec_vol
    else ! Incompressible case
       if (zlev == 1) then 
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*temp(id+1)
       else ! Interpolate to lower interface of current level
          dom%press%elts(id+1) = dom%press%elts(id+1) - grav_accel*interp (dom%adj_temp%elts(id+1), temp(id+1))
       end if
       dom%adj_temp%elts(id+1) = temp(id+1)

       if (zlev == zlevels) then !top zlev, purely diagnostic
          if (abs(dom%press%elts(id+1)-0.5_8*grav_accel*temp(id+1) - press_infty)> 1e-10_8) then
             print *, 'warning: upward integration of Lagrange multiplier not resulting in zero at top interface'
             write(6,'(A,es15.8)') "Pressure at infinity = ", dom%press%elts(id+1)-0.5_8*grav_accel*temp(id+1)
             write(6,'(A,es15.8)') "Press_infty = ", press_infty
             write(6,'(A,es15.8)') "Difference = ", dom%press%elts(id+1)-0.5_8*grav_accel*temp(id+1) - press_infty
             stop
          end if
       end if

       ! Find geopotential at interfaces using (18) in DYNAMICO
       ! Note: since mu is associated with the kinematic mass = inert mass (not the gravitational mass defined by the buyoancy)
       ! we divide by the constant reference density. This is the Boussinesq approximation.
       if (zlev == 1) then ! Save geopotential at lower interface of level zlev for interpolation in Bernoulli function
          dom%adj_geopot%elts(id+1) = surf_geopot (dom%node%elts(id+1))
       else
          dom%adj_geopot%elts(id+1) = dom%geopot%elts(id+1)
       end if
       dom%geopot%elts(id+1) = dom%adj_geopot%elts(id+1) + grav_accel*mass(id+1)/ref_density
    end if
  end subroutine integrate_pressure_up

  subroutine cal_pressure (dom, i, j, zlev, offs, dims)
    ! Integrate pressure up from surface to top layer
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)

    if (compressible) then ! Compressible case
       if (zlev == 1) then
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*mass(id+1)
       else ! Interpolate mass=rho*dz to lower interface of current level
          dom%press%elts(id+1) = dom%press%elts(id+1) - grav_accel*interp (mass(id+1), dom%adj_mass%elts(id+1))
       end if
       dom%adj_mass%elts(id+1) = mass(id+1) ! Save current mass for pressure calculation at next vertical level
    else ! Incompressible case
       if (zlev == 1) then 
          dom%press%elts(id+1) = dom%surf_press%elts(id+1) - 0.5_8*grav_accel*temp(id+1)
       else ! Interpolate to lower interface of current level
          dom%press%elts(id+1) = dom%press%elts(id+1) - grav_accel*interp (dom%adj_temp%elts(id+1), temp(id+1))
       end if
       dom%adj_temp%elts(id+1) = temp(id+1)
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

    integer                :: e, id
    real(8), dimension (3) :: Qperp_e, physics

    id = idx(i, j, offs, dims)

    ! Calculate Q_perp
    Qperp_e = Qperp (dom, i, j, z_null, offs, dims)

    ! Calculate physics
    physics = physics_velo_source (dom, i, j, zlev, offs, dims)

    do e = 1, EDGE 
       dvelo(EDGE*id+e) = - Qperp_e(e) + physics(e)*dom%len%elts(EDGE*id+e)
    end do
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

    integer                           :: id
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
   
    physics = physics_scalar_source (dom, i, j, zlev, offs, dims)

    id = idx(i, j, offs, dims)

    dmass(id+1) = - div (h_mflux, dom, i, j, offs, dims) + physics(S_MASS)
    dtemp(id+1) = - div (h_tflux, dom, i, j, offs, dims) + physics(S_TEMP)
  end subroutine scalar_trend

  subroutine cal_Laplacian_scalar (dom, i, j, zlev, offs, dims)
    ! Calculate divergence of gradient of scalar
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idW, idSW, idS

    id = idx(i, j, offs, dims)

    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idS  = idx(i,   j-1, offs, dims)

    Laplacian_scalar(S_MASS)%data(dom%id+1)%elts(id+1) = div_grad (grad_flux (mass))
    Laplacian_scalar(S_TEMP)%data(dom%id+1)%elts(id+1) = div_grad (grad_flux (temp))
  contains
    function grad_flux (scalar)
      ! Calculates gradient flux
      implicit none
      real(8), dimension(6) :: grad_flux
      real(8), dimension(:) :: scalar

      grad_flux(1) =  (scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
      grad_flux(2) =  (scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
      grad_flux(3) =  (scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)
      grad_flux(4) = -(scalar(idW+1) - scalar(id+1))  /dom%len%elts(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
      grad_flux(5) = -(scalar(id+1)  - scalar(idSW+1))/dom%len%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
      grad_flux(6) = -(scalar(idS+1) - scalar(id+1))  /dom%len%elts(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)
    end function grad_flux

    real(8) function div_grad (grad)
      ! Calculates divergence of gradient flux
      implicit none
      real(8), dimension(6) :: grad

      div_grad = (grad(1)-grad(4) + grad(5)-grad(2) + grad(3)-grad(6)) * dom%areas%elts(id+1)%hex_inv
    end function div_grad
  end subroutine cal_Laplacian_scalar

  subroutine cal_Laplacian_u (dom, i, j, zlev, offs, dims)
    ! Calculate Laplacian(u) = grad(divu) - curl(vort)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer               :: d, e, id
    real(8), dimension(3) :: grad_divu, curl_rotu

    grad_divu = gradi_e (divu, dom, i, j, offs, dims)
    curl_rotu = curlv_e (vort, dom, i, j, offs, dims)

    id = idx(i, j, offs, dims)
    do e = 1, EDGE
       Laplacian_u%data(dom%id+1)%elts(EDGE*id+e) = grad_divu(e) - curl_rotu(e)
    end do
  end subroutine cal_Laplacian_u

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

  function curlv_e (curlu, dom, i, j, offs, dims)
    ! Curl of vorticity given at triangle circumcentres x_v used in calculating curl(curl(u))
    ! output is at edges x_e
    implicit none
    real(8), dimension(3)          :: curlv_e
    real(8), dimension(:), pointer :: curlu
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW

    id   = idx(i,   j,   offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)

    curlv_e(RT+1) = (dom%vort%elts(TRIAG*id +LORT+1) - dom%vort%elts(TRIAG*idS+UPLT+1))/dom%pedlen%elts(EDGE*id+RT+1)
    curlv_e(DG+1) = (dom%vort%elts(TRIAG*id +LORT+1) - dom%vort%elts(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+DG+1)
    curlv_e(UP+1) = (dom%vort%elts(TRIAG*idW+LORT+1) - dom%vort%elts(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+UP+1)
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

    interp = (e1 + e2)/2.0_8
  end function interp

  subroutine du_grad (dom, i, j, zlev, offs, dims)
    ! Add gradients of Bernoulli and Exner to dvelo [DYNAMICO (23)-(25)]
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer                      :: e, id, idE, idN, idNE
    real(8), dimension(3)        :: gradB, gradE, theta_e
    real(8), dimension(0:N_BDRY) :: theta

    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    ! See DYNAMICO between (23)-(25), geopotential still known from step1_upw
    ! the theta multiplying the Exner gradient is the edge-averaged non-mass-weighted potential temperature
    theta(0)         = temp(id+1)/mass(id+1)
    theta(NORTH)     = temp(idN+1)/mass(idN+1)
    theta(EAST)      = temp(idE+1)/mass(idE+1)
    theta(NORTHEAST) = temp(idNE+1)/mass(idNE+1)

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
    do e = 1, EDGE
       dvelo(EDGE*id+e) = dvelo(EDGE*id+e)/dom%len%elts(EDGE*id+e) - gradB(e) - theta_e(e)*gradE(e)
    end do
  end subroutine du_grad

  subroutine cal_divu (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW
    real(8) :: u_dual_RT, u_dual_RT_W, u_dual_DG_SW, u_dual_DG, u_dual_UP, u_dual_UP_S

    id   = idx(i,   j,   offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)

    u_dual_RT    = velo(EDGE*id  +RT+1)*dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_RT_W  = velo(EDGE*idW +RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_DG    = velo(EDGE*id  +DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    = velo(EDGE*id  +UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
    u_dual_UP_S  = velo(EDGE*idS +UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)

    divu(id+1) =  (u_dual_RT-u_dual_RT_W + u_dual_DG_SW-u_dual_DG + u_dual_UP-u_dual_UP_S) * dom%areas%elts(id+1)%hex_inv
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

    id = idx(i, j, offs, dims)

    totaldmass = totaldmass + dmass(id+1)/dom%areas%elts(id+1)%hex_inv
    totalabsdmass = totalabsdmass + abs(dmass(id+1)/dom%areas%elts(id+1)%hex_inv)

    totaldtemp = totaldtemp + dtemp(id+1)/dom%areas%elts(id+1)%hex_inv
    totalabsdtemp = totalabsdtemp + abs(dtemp(id+1)/dom%areas%elts(id+1)%hex_inv)
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
end module ops_mod
