module ops_mod
  use test_case_mod
  implicit none
contains
  subroutine init_ops_mod
    implicit none
    logical :: initialized = .False.
    if (initialized) return ! initialize only once
    call init_domain_mod
    initialized = .true.
  end subroutine init_ops_mod

  subroutine step1 (dq, q, dom, p, zlev)
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: dq, q
    type(Domain) :: dom
    integer      :: p, zlev

    integer                      :: i, j, id,  n, e, s, w, ne, sw, v
    integer, dimension(0:N_BDRY) :: offs
    integer, dimension(2,N_BDRY) :: dims

    real(8) :: u_prim_UP, u_dual_UP, u_prim_DG, u_prim_DG_S, u_prim_DG_W, u_dual_DG, u_prim_RT, u_dual_RT
    real(8) :: u_prim_UP_S, u_dual_UP_S, u_prim_DG_SW, u_dual_DG_SW, u_prim_RT_W, u_dual_RT_W
    real(8) :: circ_LORT, circ_UPLT, circ_S_UPLT, circ_W_LORT, pv_LORT, pv_UPLT, pv_S_UPLT, pv_SW_LORT, pv_SW_UPLT, pv_W_LORT
    
    real(8), dimension(0:N_BDRY,scalars(1):scalars(2)) :: full
    real(8), dimension(1:EDGE)                         :: physics_flux

    logical :: S_bdry, W_bdry
    
    interface
       function physics_scalar_flux (q, dom, id, idE, idNE, idN, v, zlev, type)
         import
         implicit none
         real(8), dimension(1:EDGE)                           :: physics_scalar_flux
         type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
         type(domain)                                         :: dom
         integer                                              :: d, id, idE, idNE, idN, v, zlev
         logical, optional                                    :: type
       end function physics_scalar_flux
    end interface

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
      integer                                          :: idE, idN, idNE, idS, idSW, idW
      integer                                          :: d, id_i, idE_i, idN_i, idNE_i, idS_i, idW_i
      real(8)                                          :: circ_LORT, circ_UPLT, Phi_k
      real(8)                                          :: u_prim_UP_E, u_prim_RT_N
      real(8), dimension(scalars(1):scalars(2))        :: physics_source
      type (Coord), dimension(6)                       :: hex_nodes
      
      interface
         function physics_scalar_source (q, id, zlev)
           import
           implicit none
           real(8), dimension(scalars(1):scalars(2))            :: physics_scalar_source
           integer                                              :: id, zlev
           type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
         end function physics_scalar_source
      end interface

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

      do v = scalars(1), scalars(2)
         full(0:NORTHEAST,v) = q(v,zlev)%data(d)%elts((/id,idN,idE,idS,idW,idNE/)+1) &
                      + sol_mean(v,zlev)%data(d)%elts((/id,idN,idE,idS,idW,idNE/)+1)
      end do

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

      u_dual_RT    = velo(EDGE*id  +RT+1) * dom%pedlen%elts(EDGE*id  +RT+1)
      u_dual_RT_W  = velo(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
      u_dual_DG_SW = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
      u_dual_DG    = velo(EDGE*id  +DG+1) * dom%pedlen%elts(EDGE*id  +DG+1)
      u_dual_UP    = velo(EDGE*id  +UP+1) * dom%pedlen%elts(EDGE*id  +UP+1)
      u_dual_UP_S  = velo(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

      ! Potential vorticity
      circ_LORT   =   u_prim_RT   + u_prim_UP_E + u_prim_DG 
      circ_UPLT   = -(u_prim_DG   + u_prim_UP   + u_prim_RT_N)
      circ_W_LORT =   u_prim_RT_W + u_prim_UP   + u_prim_DG_W
      circ_S_UPLT = -(u_prim_RT   + u_prim_DG_S + u_prim_UP_S)

      pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT) / &
           (full(0,S_MASS)         * dom%areas%elts(id_i)%part(1) + &
            full(EAST,S_MASS)      * dom%areas%elts(idE_i)%part(3) + &
            full(NORTHEAST,S_MASS) * dom%areas%elts(idNE_i)%part(5))

      pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + circ_UPLT) / &
           (full(0,S_MASS)         * dom%areas%elts(id_i)%part(2) + &
            full(NORTHEAST,S_MASS) * dom%areas%elts(idNE_i)%part(4) + &
            full(NORTH,S_MASS)     * dom%areas%elts(idN_i)%part(6))

      pv_W_LORT = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_W_LORT) /  &
           (full(WEST,S_MASS)  * dom%areas%elts(idW_i)%part(1) + &
            full(0,S_MASS)     * dom%areas%elts(id_i)%part(3) + &
            full(NORTH,S_MASS) * dom%areas%elts(idN_i)%part(5))

      pv_S_UPLT = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_S_UPLT) / &
           (full(SOUTH,S_MASS)  * dom%areas%elts(idS_i)%part(2) + &
            full(EAST,S_MASS)   * dom%areas%elts(idE_i)%part(4) + &
            full(0,S_MASS)      * dom%areas%elts(id_i)%part(6))

      qe(EDGE*id+RT+1) = interp (pv_S_UPLT, pv_LORT)
      qe(EDGE*id+DG+1) = interp (pv_UPLT,   pv_LORT)
      qe(EDGE*id+UP+1) = interp (pv_UPLT,   pv_W_LORT)

      ! Vorticity (for velocity diffusion)
      vort(TRIAG*id+LORT+1) = circ_LORT / dom%triarea%elts(TRIAG*id+LORT+1) 
      vort(TRIAG*id+UPLT+1) = circ_UPLT / dom%triarea%elts(TRIAG*id+UPLT+1)

      ! Velocity divergence (for velocity diffusion)
      if (Laplace_order /= 0) &
           divu(id_i) = (u_dual_RT-u_dual_RT_W + u_dual_DG_SW-u_dual_DG + u_dual_UP-u_dual_UP_S) * dom%areas%elts(id_i)%hex_inv

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
         bernoulli(id_i) = ke(id_i) + Phi_k + dom%press%elts(id_i) / (ref_density * porosity (d, id_i, zlev))
      end if

      ! Exner function in incompressible case from geopotential
      if (.not. compressible) exner(id_i) = -Phi_k

      ! Mass and temperature fluxes
      physics_source = physics_scalar_source (q, id, zlev)

      do v = scalars(1), scalars(2)
         physics_flux = physics_scalar_flux (q, dom, id, idE, idNE, idN, v, zlev)

         horiz_flux(v)%data(d)%elts(EDGE*id+RT+1) = u_dual_RT * interp (full(0,v), full(EAST,v))      + physics_flux(RT+1)
         horiz_flux(v)%data(d)%elts(EDGE*id+DG+1) = u_dual_DG * interp (full(0,v), full(NORTHEAST,v)) + physics_flux(DG+1)
         horiz_flux(v)%data(d)%elts(EDGE*id+UP+1) = u_dual_UP * interp (full(0,v), full(NORTH,v))     + physics_flux(UP+1)

         dq(v,zlev)%data(d)%elts(id+1) = physics_source(v)
      end do
    end subroutine comput

    subroutine comp_SW
      implicit none
      integer :: d, idS, idSW, idW
      integer :: id_i, idS_i, idSW_i, idW_i
      real(8) :: circ_SW_LORT, circ_SW_UPLT, u_prim_RT_SW, u_prim_UP_SW

      d = dom%id + 1

      idW  = id+W
      idSW = id+SW
      idS  = id+S

      id_i   = id+1
      idW_i  = idW+1
      idSW_i = idSW+1
      idS_i  = idS+1

      do v = scalars(1), scalars(2)
         full(0:SOUTHWEST,v) = q(v,zlev)%data(d)%elts((/id,id,id,idS,idW,id,id,idSW/)+1) &
                      + sol_mean(v,zlev)%data(d)%elts((/id,id,id,idS,idW,id,id,idSW/)+1)
      end do

      u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1)
      u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
      u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
      u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
      u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)   

      ! Potential vorticity
      circ_SW_LORT =   u_prim_RT_SW + u_prim_UP_S  + u_prim_DG_SW
      circ_SW_UPLT = -(u_prim_RT_W  + u_prim_DG_SW + u_prim_UP_SW)

      pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / &
           (full(SOUTHWEST,S_MASS) * dom%areas%elts(idSW_i)%part(1) + &
            full(SOUTH,S_MASS)     * dom%areas%elts(idS_i)%part(3) + &
            full(0,S_MASS)         * dom%areas%elts(id_i)%part(5))

      pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_SW_UPLT) / &
           (full(SOUTHWEST,S_MASS) * dom%areas%elts(idSW_i)%part(2) + &
            full(0,S_MASS)         * dom%areas%elts(id_i)%part(4) + &
            full(WEST,S_MASS)      * dom%areas%elts(idW_i)%part(6))

      qe(EDGE*idW +RT+1) = interp (pv_W_LORT , pv_SW_UPLT)
      qe(EDGE*idSW+DG+1) = interp (pv_SW_LORT, pv_SW_UPLT)
      qe(EDGE*idS +UP+1) = interp (pv_SW_LORT, pv_S_UPLT)

      ! Vorticity (for velocity diffusion)
      vort(TRIAG*idW+LORT+1) = circ_W_LORT / dom%triarea%elts(TRIAG*idW+LORT+1) 
      vort(TRIAG*idS+UPLT+1) = circ_S_UPLT / dom%triarea%elts(TRIAG*idS+UPLT+1)
      
      ! Scalar fluxes
       do v = scalars(1), scalars(2)
         physics_flux = physics_scalar_flux (q, dom, id, idW, idSW, idS, v, zlev, .true.)
          
         horiz_flux(v)%data(d)%elts(EDGE*idW+RT+1)  = u_dual_RT_W  * interp (full(0,v), full(WEST,v))      + physics_flux(RT+1)
         horiz_flux(v)%data(d)%elts(EDGE*idSW+DG+1) = u_dual_DG_SW * interp (full(0,v), full(SOUTHWEST,v)) + physics_flux(DG+1)
         horiz_flux(v)%data(d)%elts(EDGE*idS+UP+1)  = u_dual_UP_S  * interp (full(0,v), full(SOUTH,v))     + physics_flux(UP+1)
      end do
    end subroutine comp_SW
  end subroutine step1

  subroutine post_step1 (dom, p, c, offs, dims, zlev)
    ! Correct values for vorticity and qe at pentagon points
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW, idN, idE, idNE
    real(8) :: circ_LORT, circ_UPLT, circ_S_UPLT, circ_SW_LORT, circ_SW_UPLT, circ_W_LORT
    real(8) :: pv, pv_LORT, pv_UPLT, pv_S_UPLT, pv_SW_LORT, pv_SW_UPLT, pv_W_LORT
    real(8) :: u_prim_RT, u_prim_RT_N, u_prim_RT_SW, u_prim_RT_W, u_prim_DG_SW, u_prim_UP, u_prim_UP_S, u_prim_UP_SW
    real(8), dimension(0:N_BDRY) :: full_mass

    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case
    if (c == IJMINUS) then
       id   = idx ( 0,  0, offs, dims)
       idN  = idx ( 0,  1, offs, dims)
       idE  = idx ( 1,  0, offs, dims)      
       idS  = idx ( 0, -1, offs, dims)
       idW  = idx (-1,  0, offs, dims)
       idSW = idx (-1, -1, offs, dims)

       full_mass(0:SOUTHWEST) = mass((/id,idN,idE,idS,idW,id,id,idSW/)+1) + mean_m((/id,idN,idE,idS,idW,id,id,idSW/)+1)

       u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1) 
       u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)

       circ_SW_LORT = u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW
       circ_W_LORT  = vort(TRIAG*idW+LORT+1) * dom%triarea%elts(TRIAG*idW+LORT+1)
       circ_S_UPLT  = vort(TRIAG*idS+UPLT+1) * dom%triarea%elts(TRIAG*idS+UPLT+1)

       vort(TRIAG*idSW+LORT+1) = circ_SW_LORT / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idSW+UPLT+1) = vort(TRIAG*idSW+LORT+1)

       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / &
            (full_mass(WEST)  * dom%areas%elts(idW+1)%part(6) + &
            full_mass(0)     * sum(dom%areas%elts(id+1)%part(4:5)) + &
            full_mass(SOUTH) * dom%areas%elts(idS+1)%part(3))

       pv_W_LORT = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_W_LORT) / &
            (full_mass(WEST)  * dom%areas%elts(idW+1)%part(1) + &
            full_mass(0)     * dom%areas%elts(id+1)%part(3) + &
            full_mass(NORTH) * dom%areas%elts(idN+1)%part(5))

       pv_S_UPLT = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_S_UPLT) / &
            (full_mass(SOUTH) * dom%areas%elts(idS+1)%part(2) + &
            full_mass(EAST)  * dom%areas%elts(idE+1)%part(4) + &
            full_mass(0)     * dom%areas%elts(id+1)%part(6))

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

       full_mass(0:SOUTHWEST) = mass((/id,id,idE,idS,idW,idNE,id,idSW/)+1) + mean_m((/id,id,idE,idS,idW,idNE,id,idSW/)+1)

       u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT    = velo(EDGE*id  +RT+1) * dom%len%elts(EDGE*id  +RT+1)

       circ_SW_LORT = - u_prim_RT + u_prim_RT_SW + u_prim_DG_SW
       circ_LORT    = vort(TRIAG*id  +LORT+1) * dom%triarea%elts(TRIAG*id  +LORT+1)
       circ_SW_UPLT = vort(TRIAG*idSW+UPLT+1) * dom%triarea%elts(TRIAG*idSW+UPLT+1)

       vort(TRIAG*idSW+LORT+1) = circ_SW_LORT / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idS +UPLT+1) = vort(TRIAG*idSW+LORT+1)

       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / &
            (full_mass(SOUTHWEST) * dom%areas%elts(idSW+1)%part(1) + &
            full_mass(SOUTH)     * dom%areas%elts(idS +1)%part(3) + &
            full_mass(0)         * sum(dom%areas%elts(id+1)%part(5:6)))

       pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT) / &
            (full_mass(0)         * dom%areas%elts(id  +1)%part(1) + &
            full_mass(EAST)      * dom%areas%elts(idE +1)%part(3) + &
            full_mass(NORTHEAST) * dom%areas%elts(idNE+1)%part(5))

       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_SW_UPLT) / &
            (full_mass(SOUTHWEST) * dom%areas%elts(idSW+1)%part(2) + &
            full_mass(0)         * dom%areas%elts(id  +1)%part(4) + &
            full_mass(WEST)      * dom%areas%elts(idW +1)%part(6))

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

       full_mass(0:SOUTHWEST) = mass((/id,idN,id,idS,idW,idNE,id,idSW/)+1) + mean_m((/id,idN,id,idS,idW,idNE,id,idSW/)+1)

       u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id  +UP+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)

       circ_SW_UPLT = u_prim_UP - u_prim_DG_SW - u_prim_UP_SW
       circ_UPLT    = vort(TRIAG*id  +UPLT+1) * dom%triarea%elts(TRIAG*id  +UPLT+1)
       circ_SW_LORT = vort(TRIAG*idSW+LORT+1) * dom%triarea%elts(TRIAG*idSW+LORT+1)

       vort(TRIAG*idSW+UPLT+1) = circ_SW_UPLT / dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       vort(TRIAG*idW +LORT+1) = vort(TRIAG*idSW+UPLT+1)

       pv_SW_UPLT = (dom%coriolis%elts(TRIAG*idSW+UPLT+1) + circ_SW_UPLT) / &
            (full_mass(SOUTHWEST) * dom%areas%elts(idSW+1)%part(2) + &
            full_mass(0)         * sum(dom%areas%elts(id+1)%part(3:4)) + &
            full_mass(WEST)      * dom%areas%elts(idW+1)%part(6))

       pv_UPLT = (dom%coriolis%elts(TRIAG*id+UPLT+1) + circ_UPLT) / &
            (full_mass(0)         * dom%areas%elts(id  +1)%part(2) + &
            full_mass(NORTHEAST) * dom%areas%elts(idNE+1)%part(4) + &
            full_mass(NORTH)     * dom%areas%elts(idN +1)%part(6))

       pv_SW_LORT = (dom%coriolis%elts(TRIAG*idSW+LORT+1) + circ_SW_LORT) / &
            (full_mass(SOUTHWEST) * dom%areas%elts(idSW+1)%part(1) + &
            full_mass(SOUTH)     * dom%areas%elts(idS +1)%part(3) + &
            full_mass(0)         * dom%areas%elts(id  +1)%part(5))

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

       full_mass(0:WEST) = mass((/id,idN,idE,idS,idW/)+1) + mean_m((/id,idN,idE,idS,idW/)+1)

       u_prim_RT   = velo(EDGE*id +RT+1) * dom%len%elts(EDGE*id +RT+1)
       u_prim_RT_N = velo(EDGE*idN+RT+1) * dom%len%elts(EDGE*idN+RT+1)
       u_prim_UP   = velo(EDGE*id +UP+1) * dom%len%elts(EDGE*id +UP+1)

       circ_LORT   = u_prim_RT - u_prim_RT_N - u_prim_UP
       circ_W_LORT = vort(TRIAG*idW+LORT+1) * dom%triarea%elts(TRIAG*idW+LORT+1)
       circ_S_UPLT = vort(TRIAG*idS+UPLT+1) * dom%triarea%elts(TRIAG*idS+UPLT+1)

       vort(TRIAG*id+LORT+1) = circ_LORT / dom%triarea%elts(TRIAG*id+LORT+1)
       vort(TRIAG*id+UPLT+1) = vort(TRIAG*id+LORT+1)

       pv_LORT = (dom%coriolis%elts(TRIAG*id+LORT+1) + circ_LORT) / &          
            (full_mass(EAST)  * dom%areas%elts(idE+1)%part(3) + &
            full_mass(0)     * sum(dom%areas%elts(id+1)%part(1:2)) + &
            full_mass(NORTH) * dom%areas%elts(idN+1)%part(6))

       pv_W_LORT = (dom%coriolis%elts(TRIAG*idW+LORT+1) + circ_W_LORT) / &
            (full_mass(WEST)  * dom%areas%elts(idW+1)%part(1) + &
            full_mass(0)     * dom%areas%elts(id +1)%part(3) + &
            full_mass(NORTH) * dom%areas%elts(idN+1)%part(5))

       pv_S_UPLT = (dom%coriolis%elts(TRIAG*idS+UPLT+1) + circ_S_UPLT) / &
            (full_mass(SOUTH) * dom%areas%elts(idS+1)%part(2) + &
            full_mass(EAST)  * dom%areas%elts(idE+1)%part(4) + &
            full_mass(0)     * dom%areas%elts(id +1)%part(6))

       pv_UPLT = pv_LORT

       qe(EDGE*id+RT+1) = interp (pv_LORT, pv_S_UPLT)
       qe(EDGE*id+UP+1) = interp (pv_UPLT, pv_W_LORT)
    end if
  end subroutine post_step1

  subroutine scalar_trend (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                           :: id_i

    id_i = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       dscalar(id_i) = - div (h_flux, dom, i, j, offs, dims) + dscalar(id_i) ! dscalar already contains any source terms
    else
       dscalar(id_i) = 0.0_8
    end if
  end subroutine scalar_trend

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
         import
         real(8), dimension(1:EDGE)     :: physics_velo_source
         type(domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function physics_velo_source
    end interface

    integer                :: e, id
    real(8), dimension (3) :: Qperp_e, physics

    id = idx (i, j, offs, dims)

    if (maxval (dom%mask_e%elts(EDGE*id+RT+1:EDGE*id+UP+1)) >= ADJZONE) then
       ! Calculate Q_perp
       Qperp_e = Qperp (dom, i, j, z_null, offs, dims)
       
       ! Calculate physics
       physics = physics_velo_source (dom, i, j, zlev, offs, dims)

       ! Trend
       do e = 1, EDGE
          if (dom%mask_e%elts(EDGE*id+e) >= ADJZONE) then
             dvelo(EDGE*id+e) = - Qperp_e(e) + physics(e) * dom%len%elts(EDGE*id+e)
          else
             dvelo(EDGE*id+e) = 0.0_8
          end if
       end do
    else
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
    end if
  end subroutine du_source

  subroutine du_grad (dom, i, j, zlev, offs, dims)
    ! Add gradients of Bernoulli and Exner to dvelo [DYNAMICO (23)-(25)]
    implicit none
    type(Domain)                     :: dom
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer                         :: e, id, idE, idN, idNE
    real(8), dimension(1:EDGE)      :: gradB, gradE, theta_e
    real(8), dimension(0:NORTHEAST) :: full_mass, full_temp, theta

    id = idx (i, j, offs, dims)
    
    if (maxval (dom%mask_e%elts(EDGE*id+RT+1:EDGE*id+UP+1)) >= ADJZONE) then
       idE  = idx (i+1, j,   offs, dims) 
       idN  = idx (i,   j+1, offs, dims)
       idNE = idx (i+1, j+1, offs, dims)

       full_mass(0:NORTHEAST) = mass((/id,idN,idE,id,id,idNE/)+1) + mean_m((/id,idN,idE,id,id,idNE/)+1)
       full_temp(0:NORTHEAST) = temp((/id,idN,idE,id,id,idNE/)+1) + mean_t((/id,idN,idE,id,id,idNE/)+1)

       ! See DYNAMICO between (23)-(25), geopotential still known from step1_up
       ! the theta multiplying the Exner gradient is the edge-averaged non-mass-weighted potential temperature
       theta(0)         = full_temp(0)         / full_mass(0)
       theta(EAST)      = full_temp(EAST)      / full_mass(EAST)
       theta(NORTHEAST) = full_temp(NORTHEAST) / full_mass(NORTHEAST)
       theta(NORTH)     = full_temp(NORTH)     / full_mass(NORTH)

       ! Interpolate potential temperature to edges
       theta_e(RT+1) = interp (theta(0), theta(EAST))
       theta_e(DG+1) = interp (theta(0), theta(NORTHEAST))
       theta_e(UP+1) = interp (theta(0), theta(NORTH))

       ! Incompressible: theta is normalized density perturbation, want theta_e = (rho0-rho)/rho0
       if (.not. compressible) theta_e = 1.0_8 - theta_e

       ! Calculate gradients
       gradB = gradi_e (bernoulli, dom, i, j, offs, dims)
       gradE = gradi_e (exner,     dom, i, j, offs, dims)

       ! Trend
       do e = 1, EDGE
          if (dom%mask_e%elts(EDGE*id+e) >= ADJZONE) then
             dvelo(EDGE*id+e) = dvelo(EDGE*id+e)/dom%len%elts(EDGE*id+e) - gradB(e) - theta_e(e)*gradE(e)
          else
             dvelo(EDGE*id+e) = 0.0_8
          end if
       end do
    else
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
    end if
  end subroutine du_grad

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

    id   = idx (i, j, offs, dims)

    idNW = idx (i-1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idSE = idx (i+1, j-1, offs, dims)

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

  function get_weights (dom, id, offs)
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

  subroutine cal_surf_press (q)
    implicit none
    ! Compute surface pressure and save in press_lower for upward integration
    ! Set geopotential to surface geopotential for upward integration
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q

    integer :: d, k, mass_type, p

    if (compressible) then
       mass_type = S_MASS
    else
       mass_type = S_TEMP
    end if

    call apply (set_surf_geopot, z_null)

    do d = 1, size(grid)
       grid(d)%surf_press%elts = 0.0_8
       do k = 1, zlevels
          mass   => q(mass_type,k)%data(d)%elts
          mean_m => sol_mean(mass_type,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (column_mass, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass, mean_m)
       end do
       grid(d)%surf_press%elts = grav_accel*grid(d)%surf_press%elts + p_top
       grid(d)%press_lower%elts = grid(d)%surf_press%elts
    end do
  end subroutine cal_surf_press

  subroutine column_mass (dom, i, j, zlev, offs, dims)
    ! Sum up mass
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: full_mass

    id_i = idx (i, j, offs, dims) + 1

    full_mass = mass(id_i) + mean_m(id_i)
    dom%surf_press%elts(id_i) = dom%surf_press%elts(id_i) + full_mass
  end subroutine column_mass

  subroutine set_surf_geopot (dom, i, j, zlev, offs, dims)
    ! Set initial geopotential to surface geopotential (negative for incompressible ocean flows)
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1

    if (compressible) then
       dom%geopot%elts(id_i) = surf_geopot (dom%node%elts(id_i))
    else
       dom%geopot%elts(id_i) = grav_accel * dom%topo%elts(id_i)
    end if
  end subroutine set_surf_geopot

  subroutine integrate_pressure_up (dom, i, j, zlev, offs, dims)
    ! Integrate pressure (compressible case)/Lagrange multiplier (incompressible case) and geopotential up from surface to top layer
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i
    real(8) :: full_mass, full_temp, p_upper

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    dom%geopot_lower%elts(id_i) = dom%geopot%elts(id_i)
    if (compressible) then ! compressible case
       p_upper = dom%press_lower%elts(id_i) - grav_accel*mass(id_i)
       dom%press%elts(id_i) = interp (dom%press_lower%elts(id_i), p_upper)

       exner(id_i) = c_p * (dom%press%elts(id_i)/p_0)**kappa

       dom%geopot%elts(id_i) = dom%geopot_lower%elts(id_i) + grav_accel*kappa*temp(id_i)*exner(id_i)/dom%press%elts(id_i)
    else ! incompressible case
       full_mass = mass(id_i) + mean_m(id_i)
       full_temp = temp(id_i) + mean_t(id_i)

       p_upper = dom%press_lower%elts(id_i) - grav_accel*full_temp
       dom%press%elts(id_i) = interp (dom%press_lower%elts(id_i), p_upper)

       dom%geopot%elts(id_i) = dom%geopot_lower%elts(id_i) + grav_accel*full_mass / (ref_density * porosity (d, id_i, zlev))
    end if
    dom%press_lower%elts(id_i) = p_upper
  end subroutine integrate_pressure_up

  subroutine cal_pressure (dom, i, j, zlev, offs, dims)
    ! Integrate pressure up from surface to top layer
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: full_mass, full_temp, p_upper

    id_i = idx (i, j, offs, dims) + 1

    if (compressible) then ! Compressible case
       full_mass = mass(id_i) + mean_m(id_i)
       p_upper = dom%press_lower%elts(id_i) - grav_accel * full_mass
    else ! Incompressible case
       full_temp = temp(id_i) + mean_t(id_i)
       p_upper = dom%press_lower%elts(id_i) - grav_accel * full_temp
    end if
    dom%press%elts(id_i) = interp (dom%press_lower%elts(id_i), p_upper)
    dom%press_lower%elts(id_i) = p_upper
  end subroutine cal_pressure

  subroutine post_vort (dom, p, c, offs, dims, zlev)
    ! Correct values for vorticity and qe at pentagon points
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW, idN
    real(8) :: u_prim_RT, u_prim_RT_N, u_prim_RT_SW, u_prim_RT_W, u_prim_DG_SW, u_prim_UP, u_prim_UP_S, u_prim_UP_SW

    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case
    if (c == IJMINUS) then
       idS  = idx ( 0, -1, offs, dims)
       idSW = idx (-1, -1, offs, dims)
       idW  = idx (-1,  0, offs, dims)

       u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1) 
       u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)

       vort(TRIAG*idSW+LORT+1) = (u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW) / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idSW+UPLT+1) = vort(TRIAG*idSW+LORT+1)
    end if

    ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idS+UP+1) = 0 in this case
    if (c == IPLUSJMINUS) then 
       id   = idx (PATCH_SIZE,    0, offs, dims)
       idS  = idx (PATCH_SIZE,   -1, offs, dims)
       idSW = idx (PATCH_SIZE-1, -1, offs, dims)

       u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT    = velo(EDGE*id  +RT+1)*dom%len%elts(EDGE*id  +RT+1)

       vort(TRIAG*idSW+LORT+1) = (- u_prim_RT + u_prim_RT_SW + u_prim_DG_SW) / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idS +UPLT+1) = vort(TRIAG*idSW+LORT+1)
    end if

    ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idW+RT+1) = 0 in this case
    if (c == IMINUSJPLUS) then
       id   = idx ( 0, PATCH_SIZE,   offs, dims)
       idSW = idx (-1, PATCH_SIZE-1, offs, dims)
       idW  = idx (-1, PATCH_SIZE,   offs, dims)

       u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id  +UP+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)

       vort(TRIAG*idSW+UPLT+1) = (u_prim_UP - u_prim_DG_SW - u_prim_UP_SW) / dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       vort(TRIAG*idW +LORT+1) = vort(TRIAG*idSW+UPLT+1)
    end if

    ! Parts 1, 2 of hexagon IJPLUS (upper right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*id+DG+1) = 0 in this case
    if (c == IJPLUS) then 
       id  = idx (PATCH_SIZE, PATCH_SIZE,   offs, dims)
       idN = idx (PATCH_SIZE, PATCH_SIZE+1, offs, dims)

       u_prim_RT   = velo(EDGE*id +RT+1) * dom%len%elts(EDGE*id +RT+1)
       u_prim_RT_N = velo(EDGE*idN+RT+1) * dom%len%elts(EDGE*idN+RT+1)
       u_prim_UP   = velo(EDGE*id +UP+1) * dom%len%elts(EDGE*id +UP+1)

       vort(TRIAG*id+LORT+1) = (u_prim_RT - u_prim_RT_N - u_prim_UP) / dom%triarea%elts(TRIAG*id+LORT+1)
       vort(TRIAG*id+UPLT+1) = vort(TRIAG*id+LORT+1)
    end if
  end subroutine post_vort

  subroutine cal_vort (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idN
    real(8) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_UP_E, u_prim_RT_N

    id  = idx (i,   j,   offs, dims)
    idE = idx (i+1, j,   offs, dims)
    idN = idx (i,   j+1, offs, dims)

    u_prim_RT   = velo(EDGE*id+RT +1) * dom%len%elts(EDGE*id +RT+1)
    u_prim_DG   = velo(EDGE*id+DG +1) * dom%len%elts(EDGE*id +DG+1)
    u_prim_UP   = velo(EDGE*id+UP +1) * dom%len%elts(EDGE*id +UP+1)
    u_prim_UP_E = velo(EDGE*idE+UP+1) * dom%len%elts(EDGE*idE+UP+1)
    u_prim_RT_N = velo(EDGE*idN+RT+1) * dom%len%elts(EDGE*idN+RT+1)

    vort(TRIAG*id+LORT+1) =   (u_prim_RT + u_prim_UP_E + u_prim_DG  ) / dom%triarea%elts(TRIAG*id+LORT+1)
    vort(TRIAG*id+UPLT+1) = - (u_prim_DG + u_prim_UP   + u_prim_RT_N) / dom%triarea%elts(TRIAG*id+UPLT+1)
  end subroutine cal_vort

  subroutine interp_vel_hex (dom, i, j, zlev, offs, dims)
    ! Interpolate velocity to hexagon nodes in Cartesian coordinates; uses Perot formula as also used for kinetic energy
    ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_center) ) / cell_area
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    type(Coord) :: vel, x_e, x_i
    type(Coord) :: e_zonal, e_merid
    integer     :: id, idN, idE, idNE, idS, idSW, idW
    real(8)     :: lon, lat, u_dual_RT, u_dual_UP, u_dual_DG, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    ! Fluxes normal to hexagon edges
    u_dual_RT    =  velo(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG    = -velo(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    =  velo(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)

    u_dual_RT_W  = -velo(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_DG_SW =  velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_UP_S  = -velo(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)

    ! Coordinate of hexagon centre (circumcentre)
    x_i = dom%node%elts(id+1)

    ! Sum over 6 hexagon edges
    vel = Coord (0.0_8, 0.0_8, 0.0_8)
    
    x_e = dom%midpt%elts(EDGE*id+RT+1)
    vel = vec_plus (vel, vec_scale (u_dual_RT,    vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*idW+RT+1)
    vel = vec_plus (vel, vec_scale (u_dual_RT_W,  vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*id+DG+1)
    vel = vec_plus (vel, vec_scale (u_dual_DG,    vec_minus(x_e, x_i)))
 
    x_e = dom%midpt%elts(EDGE*idSW+DG+1)
    vel = vec_plus (vel, vec_scale (u_dual_DG_SW, vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*id+UP+1)
    vel = vec_plus (vel, vec_scale (u_dual_UP,    vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*idS+UP+1)
    vel = vec_plus (vel, vec_scale (u_dual_UP_S,  vec_minus(x_e, x_i)))

    vel = vec_scale (dom%areas%elts(id+1)%hex_inv, vel) ! construct velocity at hexagonal node

    ! Project velocity onto zonal and meridional directions
    call cart2sph (x_i, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),             0.0_8) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    dom%u_zonal%elts(id+1) = inner (vel, e_zonal)
    dom%v_merid%elts(id+1) = inner (vel, e_merid)
  end subroutine interp_vel_hex

  subroutine cal_Laplacian_scalar (dom, i, j, zlev, offs, dims)
    ! Calculate divergence of gradient of scalar
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idW, idSW, idS
    integer :: id_i, idE_i, idNE_i, idN_i, idW_i, idSW_i, idS_i

    id   = idx (i, j, offs, dims)
    id_i = id+1

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

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

      grad_flux(1) =  (scalar(idE_i) - scalar(id_i))  /dom%len%elts(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
      grad_flux(2) =  (scalar(id_i)  - scalar(idNE_i))/dom%len%elts(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
      grad_flux(3) =  (scalar(idN_i) - scalar(id_i))  /dom%len%elts(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)
      grad_flux(4) = -(scalar(idW_i) - scalar(id_i))  /dom%len%elts(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
      grad_flux(5) = -(scalar(id_i)  - scalar(idSW_i))/dom%len%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
      grad_flux(6) = -(scalar(idS_i) - scalar(id_i))  /dom%len%elts(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)
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

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)

    Laplacian(EDGE*id+RT+1) = -(vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1))/dom%pedlen%elts(EDGE*id+RT+1)

    if (dom%pedlen%elts(EDGE*id+DG+1) /= 0.0_8) then
       Laplacian(EDGE*id+DG+1) = -(vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+DG+1)
    else
       Laplacian(EDGE*id+DG+1) = 0.0_8
    end if

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

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

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

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)

    curlv_e(RT+1) = (curl(TRIAG*id +LORT+1) - curl(TRIAG*idS+UPLT+1)) / dom%pedlen%elts(EDGE*id+RT+1)
    curlv_e(DG+1) = (curl(TRIAG*id +LORT+1) - curl(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+DG+1)
    curlv_e(UP+1) = (curl(TRIAG*idW+LORT+1) - curl(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+UP+1)
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

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)

    div = (hflux(EDGE*id+RT+1)-hflux(EDGE*idW+RT+1) + hflux(EDGE*idSW+DG+1)-hflux(EDGE*id+DG+1) &
         + hflux(EDGE*id+UP+1)-hflux(EDGE*idS+UP+1)) * dom%areas%elts(id+1)%hex_inv
  end function div

  real(8) function interp (e1, e2)
    ! Centred average interpolation of quantities e1 and e2
    implicit none
    real(8) :: e1, e2

    interp = 0.5 * (e1 + e2)
  end function interp

  subroutine cal_divu (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i, idS, idW, idSW
    real(8) :: u_dual_RT, u_dual_RT_W, u_dual_DG_SW, u_dual_DG, u_dual_UP, u_dual_UP_S

    id   = idx (i,   j,   offs, dims)
    id_i = id+1

    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)

    u_dual_RT    = velo(EDGE*id  +RT+1)*dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_RT_W  = velo(EDGE*idW +RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_DG    = velo(EDGE*id  +DG+1)*dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    = velo(EDGE*id  +UP+1)*dom%pedlen%elts(EDGE*id+UP+1)
    u_dual_UP_S  = velo(EDGE*idS +UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)

    divu(id_i) =  (u_dual_RT-u_dual_RT_W + u_dual_DG_SW-u_dual_DG + u_dual_UP-u_dual_UP_S) * dom%areas%elts(id_i)%hex_inv
  end subroutine cal_divu

  real(8) function porosity (d, id_i, zlev)
    ! Returns porosity at position given by (d, id_i, zlev)
    implicit none
    integer :: d, id_i, zlev

    porosity = 1.0_8 + (alpha - 1.0_8) * penal_node(zlev)%data(d)%elts(id_i)
  end function porosity

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
