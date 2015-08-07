module init_mod
  use geom_mod
  use domain_mod
  use arch_mod
  implicit none
  real(8), parameter :: YANGLE = 0.0_8

contains
  subroutine init_init_mod()
    logical :: initialized = .False.
    if (initialized) return ! initialize only once
    call init_sphere_mod()
    call init_domain_mod()
    call init_arch_mod()
    initialized = .True.
  end subroutine init_init_mod

  subroutine init_grid()
    integer k, v
    integer d
    integer s
    integer b
    integer p
    integer i_
    integer loz

    allocate (grid(n_domain(rank+1)))
    allocate (sol(S_MASS:S_VELO,1:zlevels), trend(S_MASS:S_VELO,1:zlevels), wav_coeff(S_MASS:S_VELO, 1:zlevels))
    allocate (horiz_massflux(1:zlevels))

    do k = 1, zlevels
       do v = S_MASS, S_VELO
          call init_Float_Field(sol(v,k), v)
          call init_Float_Field(trend(v,k), v)
       end do
       call init_Float_Field(horiz_massflux(k), S_VELO)
    end do

    if (penalize) call init_Float_Field(penal, S_MASS)

    do d = 1, n_domain(rank+1)
       call init_Domain(grid(d))
       do k = 1, zlevels
          call init(sol(S_MASS,k)%data(d), 1)
          call init(sol(S_VELO,k)%data(d), EDGE)
       end do
    end do
    
    !  initializes grid for icosahedron
    do d = 1, n_domain(rank+1)

       grid(d)%id = d-1
       
       do s = 1, N_BDRY
          b = add_bdry_patch_Domain(grid(d), s)
       end do
       
       grid(d)%penta = .False.

       p = add_patch_Domain(grid(d), min_level-1)

       grid(d)%patch%elts(p+1)%children = 0
       grid(d)%patch%elts(p+1)%neigh    = (/ ( i_ , i_ = -1, -N_BDRY, -1 ) /)
    end do

    do loz = 1, N_ICOSAH_LOZANGE
       call set_penta(N_SUB_DOM*(loz - 1),                                 SOUTHWEST)
       call set_penta(N_SUB_DOM*(loz - 1) + N_SUB_DOM_PER_DIM - 1,         SOUTHEAST)
       call set_penta(N_SUB_DOM*(loz - 1) + N_SUB_DOM - N_SUB_DOM_PER_DIM, NORTHWEST)
       call set_penta(N_SUB_DOM*(loz - 1) + N_SUB_DOM - 1,                 NORTHEAST)
    end do

    call init_connections()
    call init_coordinates()

    do d = 1, n_domain(rank+1)
       if (penalize) call init(penal%data(d), grid(d)%node%length)
    end do
  end subroutine init_grid

  subroutine init_coordinates()
    real(8), dimension(4) :: lat
    real(8), dimension(10) :: lon
    integer ii, jj
    integer loz
    type(Coord) ne, se, sw, nw
    type(Coord), dimension(2,2) :: cnr
    integer j, i
    integer d_glo
    integer d
    integer k

    lat = (/-MATH_PI/2.0_8, -atan(1.0_8/2.0_8), atan(1.0_8/2.0_8), &
         MATH_PI/2.0_8/)
    lon = (/ ((MATH_PI*k)/5.0_8, k = 0, 10-1) /)

    do ii = 1, 2
       do jj = 1, 5
          loz = 5*ii - 5 + (jj - 1)

          ne = sph2cart(lon(modulo(ii + 2*jj - 2, 10) + 1), lat(ii+1))
          se = sph2cart(lon(ii+2*jj-2),                     lat(ii)  )
          sw = sph2cart(lon(modulo(ii + 2*jj - 4, 10) + 1), lat(ii+1))
          nw = sph2cart(lon(ii+2*jj-2),                     lat(ii+2))

          call yrotate(nw, cnr(1,2), YANGLE); call yrotate(ne, cnr(2,2), YANGLE)
          call yrotate(sw, cnr(1,1), YANGLE); call yrotate(se, cnr(2,1), YANGLE)

          do j = 1, N_SUB_DOM_PER_DIM
             do i = 1, N_SUB_DOM_PER_DIM
                d_glo = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 1) + (i - 1)

                if (.not. owner(d_glo+1) .eq. rank) cycle

                d = loc_id(d_glo+1)
                call assign_coord(grid(d+1), 1, get_J0_coord(i, j, &
                     DOMAIN_LEVEL), get_J0_coord(i, j - 1, &
                     DOMAIN_LEVEL), get_J0_coord(i - 1, j - 1, &
                     DOMAIN_LEVEL), get_J0_coord(i - 1, j, &
                     DOMAIN_LEVEL))
             end do
          end do
       end do
    end do

  contains

    subroutine  yrotate(c_in, c_out, angle)
      real(8), intent(in) :: angle
      type(Coord), intent(in) :: c_in
      type(Coord), intent(out) :: c_out

      c_out%x =  c_in%x*cos(angle) - c_in%z*sin(angle)
      c_out%y =  c_in%y
      c_out%z =  c_in%x*sin(angle) + c_in%z*cos(angle)
    end subroutine yrotate

    type(Coord) recursive function get_J0_coord(i, j, l) result(c)
      integer i
      integer j
      integer l

      if (l .gt. 0) then
         c = mid_pt( get_J0_coord(i/2,       j/2,       l-1), &
                     get_J0_coord((i + 1)/2, (j + 1)/2, l-1) )
         return
      else
         c = cnr(i+1,j+1)
         return
      end if
    end function get_J0_coord

  end subroutine init_coordinates

  subroutine ccentre_penta(dom, p)
    type(Domain) dom
    integer p
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,9) :: dims
    type(Coord) pt3
    integer id_mm
    integer id_0m1
    integer id_m10
    integer id_0m2
    integer id_m20
    integer id_m12
    integer id_m21
    integer id, idN, idNE, idW, idE, idS
    integer i, j

    call get_offs_Domain(dom, p, offs, dims)

    if (is_penta(dom, p, IJMINUS - 1)) then
       id_mm  = idx(-1, -1, offs, dims)
       id_0m1 = idx( 0, -1, offs, dims)
       id_m10 = idx(-1,  0, offs, dims)

       dom%ccentre%elts(TRIAG*id_mm+LORT+1) = &
            circumcentre(dom%node%elts(idx(0,0,offs,dims)+1), dom%node%elts(id_0m1+1), dom%node%elts(id_m10+1))
       
       dom%ccentre%elts(TRIAG*id_mm+UPLT+1) = dom%ccentre%elts(TRIAG*id_mm+LORT+1)
       
       id_0m2 = idx( 0, -2, offs, dims)
       id_m20 = idx(-2,  0, offs, dims)
       
       pt3 = mid_pt(dom%node%elts(id_0m2+1), dom%node%elts(id_m20+1))
       
       id_m12 = idx(-1, -2, offs, dims)
       id_m21 = idx(-2, -1, offs, dims)
       
       dom%ccentre%elts(TRIAG*id_m12+LORT+1) = &
            circumcentre(dom%node%elts(id_0m1+1), dom%node%elts(id_0m2+1), pt3)

       dom%ccentre%elts(TRIAG*id_m21+UPLT+1) = &
            circumcentre(pt3, dom%node%elts(id_m20+1), dom%node%elts(id_m10+1))
    end if
    
    if (is_penta(dom, p, IPLUSJMINUS - 1)) then
       i = PATCH_SIZE
       j = -1

       id   = idx(i,     j,     offs, dims)
       idN  = idx(i,     j + 1, offs, dims)
       idNE = idx(i + 1, j + 1, offs, dims)
       idW  = idx(i - 1, j,     offs, dims)

       dom%ccentre%elts(TRIAG*idW+LORT+1) = &
            circumcentre(dom%node%elts(idW+1), dom%node%elts(idN+1), dom%node%elts(idNE+1))
       
       dom%ccentre%elts(TRIAG*id+UPLT+1) = dom%ccentre%elts(TRIAG*idW+LORT+1)
    end if
    
    if (is_penta(dom, p, IMINUSJPLUS - 1)) then
       i = -1
       j = PATCH_SIZE

       id   = idx(i,     j,     offs, dims)
       idE  = idx(i + 1, j,     offs, dims)
       idNE = idx(i + 1, j + 1, offs, dims)
       idS  = idx(i,     j - 1, offs, dims)
       
       dom%ccentre%elts(TRIAG*idS+UPLT+1) = &
            circumcentre(dom%node%elts(idS+1), dom%node%elts(idNE+1), dom%node%elts(idE+1))
       
       dom%ccentre%elts(TRIAG*id+LORT+1) = dom%ccentre%elts(TRIAG*idS+UPLT+1)
    end if
    
    if (is_penta(dom, p, IJPLUS - 1)) then
       i = PATCH_SIZE
       j = PATCH_SIZE

       id  = idx(i,     j,     offs, dims)
       idN = idx(i,     j + 1, offs, dims)
       idE = idx(i + 1, j,     offs, dims)

       dom%ccentre%elts(TRIAG*id+LORT+1) = &
            circumcentre(dom%node%elts(id+1), dom%node%elts(idN+1), dom%node%elts(idE+1))
       
       dom%ccentre%elts(TRIAG*id+UPLT+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
    end if
  end subroutine ccentre_penta

  subroutine lengths(dom, p, i, j, zlev, offs, dims)
    type(Domain) dom
    integer p, i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, idS, idW, idN, idE, idNE

    id   = idx(i,     j,     offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if (j .ge. PATCH_SIZE + 1) then
       dom%len%elts(EDGE*id+RT+1) = dist(dom%node%elts(id+1), dom%node%elts(idE+1))
       return
    end if
    
    if (i .ge. PATCH_SIZE + 1) then
       dom%len%elts(EDGE*id+UP+1) = dist(dom%node%elts(id+1), dom%node%elts(idN+1))
       return
    end if
    
    dom%len%elts(EDGE*id+RT+1)    = dist(dom%node%elts(id+1), dom%node%elts(idE+1))
    dom%len%elts(EDGE*id+DG+1)    = dist(dom%node%elts(idNE+1), dom%node%elts(id+1))
    dom%len%elts(EDGE*id+UP+1)    = dist(dom%node%elts(id+1), dom%node%elts(idN+1))
    dom%pedlen%elts(EDGE*id+RT+1) = dist(dom%ccentre%elts(TRIAG*idS+UPLT+1), dom%ccentre%elts(LORT+TRIAG*id+1))
    dom%pedlen%elts(EDGE*id+DG+1) = dist(dom%ccentre%elts(TRIAG*id+UPLT+1), dom%ccentre%elts(LORT+TRIAG*id+1))
    dom%pedlen%elts(EDGE*id+UP+1) = dist(dom%ccentre%elts(TRIAG*id+UPLT+1), dom%ccentre%elts(LORT+TRIAG*idW+1))
    
    if (i .eq. -1 .and. j .eq. -1 .and. is_penta(dom, p, IJMINUS - 1)) then
       dom%pedlen%elts(EDGE*id+DG+1) = 0
       dom%len%elts(EDGE*id+1) = dist(dom%node%elts(idE+1), dom%node%elts(idN+1))
    end if
    
    if (i .eq. PATCH_SIZE .and. j .eq. PATCH_SIZE .and. is_penta(dom, p, IJPLUS - 1)) then
       dom%len%elts(EDGE*id+DG+1) = dist(dom%node%elts(idE+1), dom%node%elts(idN+1))
    end if
  end subroutine lengths

  subroutine init_geometry()
    integer d
    integer i
    
    do d = 1, size(grid)
       call init(grid(d)%ccentre, grid(d)%node%length*TRIAG)

       do i = 1, grid(d)%node%length*TRIAG
          call init_Coord(grid(d)%ccentre%elts(i), 0.0_8, 0.0_8, 0.0_8)
       end do

       call init(grid(d)%midpt, grid(d)%node%length*EDGE)

       do i = 1, grid(d)%node%length*EDGE
          call init_Coord(grid(d)%midpt%elts(i), 0.0_8, 0.0_8, 0.0_8)
       end do

       call init(grid(d)%areas,      grid(d)%node%length)
       call init(grid(d)%pedlen,     grid(d)%node%length*EDGE)
       call init(grid(d)%len,        grid(d)%node%length*EDGE)
       call init(grid(d)%triarea,    grid(d)%node%length*TRIAG)
       call init(grid(d)%coriolis,    grid(d)%node%length*TRIAG)
       call init(grid(d)%windstress, grid(d)%node%length*EDGE)
       call init(grid(d)%topo,       grid(d)%node%length)
    end do
  end subroutine init_geometry

  subroutine precompute_geometry()
    integer d, k
    
    call apply_onescale2(ccentre, min_level-1, z_null, -2, 1)

    do d = 1, n_domain(rank+1)
       call ccentre_penta(grid(d), 1)
    end do

    call apply_onescale2(midpt,      min_level-1, z_null, -1, 2)
    call apply_onescale2(cpt_areas,  min_level-1, z_null, -1, 2)
    call apply_onescale2(lengths,    min_level-1, z_null, -1, 2)
    
    call apply_onescale(cpt_triarea, min_level-1, z_null, -1, 1)
    call apply_onescale(coriolis,    min_level-1, z_null, -1, 1)
    
    do k = 1, zlevels
       do d = 1, size(grid)
          call init(trend(S_MASS,k)%data(d),   grid(d)%node%length)
          call init(trend(S_VELO,k)%data(d),   grid(d)%node%length*EDGE); trend(S_VELO,k)%data(d)%elts = 0
          call init(horiz_massflux(k)%data(d), grid(d)%node%length*EDGE)
       end do
    end do

    do d = 1, size(grid)
       call init(grid(d)%bernoulli, grid(d)%node%length)
       call init(grid(d)%kin_energy, grid(d)%node%length)
       call init(grid(d)%qe,   grid(d)%node%length*EDGE)
       call init(grid(d)%divu, grid(d)%node%length)
       call init(grid(d)%vort, grid(d)%node%length*TRIAG)
    end do

  end subroutine precompute_geometry

  subroutine init_connections()
    logical pole_assigned(2)
    integer, dimension(2) :: neigh_over_pole
    integer ii, jj
    integer loz
    integer i, j
    integer split
    integer d_glo
    integer d
    integer k
    integer rot
    integer ngb_loz
    integer s
    
    pole_assigned = .False.

    do ii = 1, 2
       do jj = 1, 5
          loz = 5*ii - 5 + (jj - 1)
          do i = 1, N_SUB_DOM_PER_DIM
             do j = 1, N_SUB_DOM_PER_DIM
                split = N_SUB_DOM_PER_DIM*(j - 1) + (i - 1)
                d_glo = N_SUB_DOM*loz + split
                
                if (.not. owner(d_glo+1) .eq. rank) then
                   ! check if pole master on other rank
                   if (ii - 1 .eq. 1 .and. i - 1 .eq. 0 .and. j - 1 .eq. N_SUB_DOM_PER_DIM - 1) pole_assigned(1) = .True.
                   if (ii - 1 .eq. 0 .and. j - 1 .eq. 0 .and. i - 1 .eq. N_SUB_DOM_PER_DIM - 1) pole_assigned(2) = .True.
                   cycle
                end if
                
                d = loc_id(d_glo+1)
                grid(d+1)%neigh = NONE

                neigh_over_pole = (/5*ii - 5 + modulo(jj + 1, 5), 5*ii - 5 + modulo(jj - 3, 5)/)
                
                if (ii - 1 .eq. 1 .and. i - 1 .eq. 0 .and. j - 1 .eq. N_SUB_DOM_PER_DIM - 1) then
                   grid(d+1)%neigh(NORTHWEST) = POLE
                   grid(d+1)%neigh_over_pole = N_SUB_DOM*neigh_over_pole + split
                   
                   if (.not. pole_assigned(1)) then
                      grid(d+1)%pole_master(2) = .True.
                      pole_assigned(1) = .True.
                   end if
                else
                   if (ii - 1 .eq. 0 .and. j - 1 .eq. 0 .and. i - 1 .eq. N_SUB_DOM_PER_DIM - 1) then
                      grid(d+1)%neigh(SOUTHEAST) = POLE
                      grid(d+1)%neigh_over_pole = N_SUB_DOM*neigh_over_pole + split
                      if (.not. pole_assigned(2)) then
                         grid(d+1)%pole_master(1) = .True.
                         pole_assigned(2) = .True.
                      end if
                   end if
                end if
                
                call init(grid(d+1)%neigh_pa_over_pole, min_level*2)
                
                grid(d+1)%neigh_pa_over_pole%elts(min_level*2-1:min_level*2) = 1
                grid(d+1)%neigh_rot = 0
                
                grid(d+1)%neigh(NORTH) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*j + (i - 1)
                grid(d+1)%neigh(EAST)  = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 1) + i - 1 + 1
                grid(d+1)%neigh(SOUTH) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 2) + (i - 1)
                grid(d+1)%neigh(WEST) = (N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 1) + i - 1) - 1
                
                if (i .lt. N_SUB_DOM_PER_DIM) then
                   if (j .lt. N_SUB_DOM_PER_DIM) then
                      grid(d+1)%neigh(IJPLUS) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*j + i - 1 + 1
                   end if
                   
                   if (j - 1 .gt. 0) then
                      grid(d+1)%neigh(IPLUSJMINUS) = N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 2) + i - 1 + 1
                   end if
                end if
                
                if (i - 1 .gt. 0) then
                   if (j .lt. N_SUB_DOM_PER_DIM) then
                      grid(d+1)%neigh(IMINUSJPLUS) = (N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*j + i - 1) - 1
                   end if
                   
                   if (j - 1 .gt. 0) then
                      grid(d+1)%neigh(IJMINUS) = (N_SUB_DOM*loz + N_SUB_DOM_PER_DIM*(j - 2) + i - 1) - 1
                   end if
                end if
                
                if (j .eq. N_SUB_DOM_PER_DIM) then
                   rot = ii - 1
                   ngb_loz = (5 + modulo(ii + jj - 2, 5))*N_SUB_DOM
                   call set_dom_neigh(d, NORTH, ngb_loz, i - 1, 0, NORTH, rot)
                   
                   if (i .lt. N_SUB_DOM_PER_DIM) then
                      call set_dom_neigh(d, IJPLUS, ngb_loz, i, 0, NORTH, rot)
                   end if
                   
                   if (i - 1 .gt. 0) then
                      call set_dom_neigh(d, IMINUSJPLUS, ngb_loz, i - 2, 0, NORTH, rot)
                   end if
                end if
                
                if (i .eq. N_SUB_DOM_PER_DIM) then
                   rot = 1 - (ii - 1)
                   ngb_loz = (0 + modulo(jj, 5))*N_SUB_DOM
                   call set_dom_neigh(d, EAST, ngb_loz, 0, j - 1, EAST, rot)
                   
                   if (j .lt. N_SUB_DOM_PER_DIM) then
                      call set_dom_neigh(d, IJPLUS, ngb_loz, 0, j, EAST, rot)
                   end if
                   
                   if (j - 1 .gt. 0) then
                      call set_dom_neigh(d, IPLUSJMINUS, ngb_loz, 0, j - 2, EAST, rot)
                   end if
                end if
                
                if (j - 1 .eq. 0) then
                   rot = 1 - (ii - 1)
                   ngb_loz = (0 + modulo(ii + jj - 3, 5))*N_SUB_DOM
                   
                   call set_dom_neigh(d, SOUTH, ngb_loz, i - 1, N_SUB_DOM_PER_DIM - 1, SOUTH, rot)
                   
                   if (i .lt. N_SUB_DOM_PER_DIM) then
                      call set_dom_neigh(d, IPLUSJMINUS, ngb_loz, i, N_SUB_DOM_PER_DIM - 1, SOUTH, rot)
                   end if
                   
                   if (i - 1 .gt. 0) then
                      call set_dom_neigh(d, IJMINUS, ngb_loz, i - 2, N_SUB_DOM_PER_DIM - 1, SOUTH, rot)
                   end if
                end if
                
                if (i - 1 .eq. 0) then
                   rot = ii - 1
                   ngb_loz = (5 + modulo(jj - 2, 5))*N_SUB_DOM
                   
                   call set_dom_neigh(d, WEST, ngb_loz, N_SUB_DOM_PER_DIM - 1, j - 1, WEST, rot)
                   
                   if (j .lt. N_SUB_DOM_PER_DIM) then
                      call set_dom_neigh(d, IMINUSJPLUS, ngb_loz, N_SUB_DOM_PER_DIM - 1, j, WEST, rot)
                   end if
                   
                   if (j - 1 .gt. 0) then
                      call set_dom_neigh(d, IJMINUS, ngb_loz, N_SUB_DOM_PER_DIM - 1, j - 2, WEST, rot)
                   end if
                end if
                do s = 1, N_BDRY
                   grid(d+1)%bdry_patch%elts(s+1)%neigh = 1
                end do
             end do
          end do
       end do
    end do
  end subroutine init_connections

  subroutine set_penta(d_glo, side)
    integer d_glo
    integer side
    
    if (owner(d_glo+1) .eq. rank) then
       grid(loc_id(d_glo+1) + 1)%penta(side) = .True.
    end if
  end subroutine set_penta

  subroutine assign_coord(dom, p, ne, se, sw, nw)
    type(Domain) dom
    integer p
    type(Coord) ne
    type(Coord) se
    type(Coord) sw
    type(Coord) nw
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer l
    integer pc_incr
    integer pf_offs
    integer j
    integer i

    call get_offs_Domain(dom, p, offs, dims)

    dom%node%elts(offs(NORTHEAST+1) + 1) = project_on_sphere(ne)
    dom%node%elts(offs(EAST+1) + 1)      = project_on_sphere(se)
    dom%node%elts(offs(1) + 1)           = project_on_sphere(sw)
    dom%node%elts(offs(NORTH+1) + 1)     = project_on_sphere(nw)

    do l = 1, PATCH_LEVEL
       pc_incr = 2**(PATCH_LEVEL - l + 1)
       pf_offs = pc_incr/2
       
       do j = 1, PATCH_SIZE, pc_incr
          do i = 1, PATCH_SIZE, pc_incr
             dom%node%elts(idx(i + pf_offs - 1, j - 1, offs, dims) + 1) = &
                  mid_pt(dom%node%elts(idx(i - 1,           j - 1, offs, dims) + 1), &
                         dom%node%elts(idx(i + pc_incr - 1, j - 1, offs, dims) + 1))
             
             dom%node%elts(idx(i - 1, j + pf_offs - 1, offs, dims) + 1) = &
                  mid_pt(dom%node%elts(idx(i - 1, j - 1,           offs, dims) + 1), &
                         dom%node%elts(idx(i - 1, j + pc_incr - 1, offs, dims) + 1))
             
             dom%node%elts(idx(i + pf_offs - 1, j + pf_offs - 1, offs, dims) + 1) = &
                  mid_pt(dom%node%elts(idx(i + pc_incr - 1, j + pc_incr - 1, offs, dims) + 1), &
                         dom%node%elts(idx(i - 1,           j - 1,           offs, dims) + 1))
          end do
          
          i = PATCH_SIZE
          dom%node%elts(idx(i, j + pf_offs - 1, offs, dims) + 1) = &
               mid_pt(dom%node%elts(idx(i, j - 1,           offs, dims) + 1), &
                      dom%node%elts(idx(i, j + pc_incr - 1, offs, dims) + 1))
       end do
       
       j = PATCH_SIZE
       do i = 1, PATCH_SIZE, pc_incr
          dom%node%elts(idx(i + pf_offs - 1, j, offs, dims) + 1) = &
               mid_pt(dom%node%elts(idx(i - 1,           j, offs, dims) + 1), &
                      dom%node%elts(idx(i + pc_incr - 1, j, offs, dims) + 1))
       end do
    end do
  end subroutine assign_coord

  subroutine cpt_areas(dom, p, i, j, zlev, offs, dims)
    type(Domain) dom
    integer p
    integer i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer idNW, idN, idNE, idW, id, idE, idSW, idS, idSE
    type(Areas) area
    
    idNW = idx(i - 1, j + 1, offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idSE = idx(i + 1, j - 1, offs, dims)
    
    call init_Areas(area, dom%node%elts(id+1), &
         
         (/ dom%ccentre%elts(TRIAG*id  +LORT+1), &
            dom%ccentre%elts(TRIAG*id  +UPLT+1), &
            dom%ccentre%elts(TRIAG*idW +LORT+1), &
            dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
            dom%ccentre%elts(TRIAG*idSW+LORT+1), &
            dom%ccentre%elts(TRIAG*idS +UPLT+1) /), &
            
            (/ dom%midpt%elts(EDGE*id  +RT+1), &
               dom%midpt%elts(EDGE*id  +DG+1), &
               dom%midpt%elts(EDGE*id  +UP+1), &
               dom%midpt%elts(EDGE*idW +RT+1), &
               dom%midpt%elts(EDGE*idSW+DG+1), &
               dom%midpt%elts(EDGE*idS +UP+1)/) )
    
    if (j .ge. PATCH_SIZE + 1) then
       dom%areas%elts(id+1)%part(4:6) = area%part(4:6)
       return
    end if

    if (i .ge. PATCH_SIZE + 1) then
       dom%areas%elts(id+1)%part(3:5) = area%part(3:5)
       return
    end if

    dom%areas%elts(id+1) = area
  end subroutine cpt_areas

  integer function sub_dom_id(i, j, s, rot)
    integer i, j
    integer s
    integer rot
    integer, dimension(2) :: ij

    ij = (/i, j/)
    if (rot .eq. 1) then
       ij(modulo(s, 2) + 1) = N_SUB_DOM_PER_DIM - 1 - ij(modulo(s, 2) + 1)
       ij = (/ij(2), ij(1)/)
    end if

    sub_dom_id = ij(2)*N_SUB_DOM_PER_DIM + ij(1)
  end function sub_dom_id

  subroutine ccentre(dom, p, i, j, zlev, offs, dims)
    type(Domain) dom
    integer p
    integer i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, idN, idE, idNE, idS, idW

    id   = idx(i, j, offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)

    dom%ccentre%elts(TRIAG*id+LORT+1) = circumcentre(dom%node%elts(id+1), dom%node%elts(idNE+1), dom%node%elts(idE+1))
    dom%ccentre%elts(TRIAG*id+UPLT+1) = circumcentre(dom%node%elts(id+1), dom%node%elts(idN+1),  dom%node%elts(idNE+1))
  end subroutine ccentre

  subroutine cpt_triarea(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, idN, idE, idNE
    
    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    
    dom%triarea%elts(TRIAG*id+LORT+1) = dom%areas%elts(id+1)%part(1) + dom%areas%elts(idE+1)%part(3) &
         + dom%areas%elts(idNE+1)%part(5)
    dom%triarea%elts(TRIAG*id+UPLT+1) = dom%areas%elts(id+1)%part(2) + dom%areas%elts(idNE+1)%part(4) &
         + dom%areas%elts(idN+1)%part(6)
  end subroutine cpt_triarea

  subroutine coriolis(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, idN, idE, idNE

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    
    dom%coriolis%elts(TRIAG*id+LORT+1) = - &
         (2.0_8*dom%ccentre%elts(LORT+TRIAG*id+1)%z)/radius*omega* &
         (dom%areas%elts(id+1)%part(1) + dom%areas%elts(idE+1)%part(3) + dom%areas%elts(idNE+1)%part(5))
    
    dom%coriolis%elts(TRIAG*id+UPLT+1) = - &
         (2.0_8*dom%ccentre%elts(TRIAG*id+UPLT+1)%z)/radius*omega* &
         (dom%areas%elts(id+1)%part(2) + dom%areas%elts(idNE+1)%part(4) + dom%areas%elts(idN+1)%part(6))
  end subroutine coriolis

  subroutine set_dom_neigh(d, s, ngb_loz, i, j, s1, rot)
    integer d
    integer s
    integer ngb_loz
    integer i
    integer j
    integer s1
    integer rot

    call set_neigh_Domain(grid(d+1), s, ngb_loz + sub_dom_id(i, j, s1 - 1, rot), rot)
  end subroutine set_dom_neigh

  subroutine set_level(dom, p, i, j, zlev, offs, dims)
    type(Domain) dom
    integer p
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    
    dom%level%elts(idx(i,j,offs,dims)+1) = dom%patch%elts(p+1)%level
  end subroutine set_level

  subroutine midpt(dom, p, i, j, zlev, offs, dims)
    type(Domain) dom
    integer p
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idN
    integer idE
    integer idNE

    id   = idx(i,     j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if (j .ge. PATCH_SIZE + 1) then
       dom%midpt%elts(EDGE*id+RT+1) = mid_pt(dom%node%elts(id+1), dom%node%elts(idE+1))
       return
    end if
    
    if (i .ge. PATCH_SIZE + 1) then
       dom%midpt%elts(EDGE*id+UP+1) = mid_pt(dom%node%elts(id+1), dom%node%elts(idN+1))
       return
    end if

    if (j .eq. -1) then
       if (i .eq. -1 .and. is_penta(dom, p, IJMINUS - 1)) then
          dom%midpt%elts(EDGE*id+DG+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
          return
       else
          if (i .eq. PATCH_SIZE .and. is_penta(dom, p, IPLUSJMINUS - 1)) then
             dom%midpt%elts(EDGE*id+UP+1) = dom%ccentre%elts(TRIAG*id+UPLT+1)
             return
          end if
       end if
    else
       if (j .eq. PATCH_SIZE .and. is_penta(dom, p, IMINUSJPLUS - 1)) then
          if (i .eq. -1) then
             dom%midpt%elts(EDGE*id+RT+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
             return
          end if
       end if
    end if
    
    dom%midpt%elts(EDGE*id+RT+1) = mid_pt(dom%node%elts(id+1),   dom%node%elts(idE+1))
    dom%midpt%elts(EDGE*id+DG+1) = mid_pt(dom%node%elts(idNE+1), dom%node%elts(id+1))
    dom%midpt%elts(EDGE*id+UP+1) = mid_pt(dom%node%elts(id+1),   dom%node%elts(idN+1))
    
    if (j .eq. PATCH_SIZE) then
       if (i .eq. PATCH_SIZE .and. is_penta(dom, p, IJPLUS - 1)) then
          dom%midpt%elts(EDGE*id+DG+1) = dom%ccentre%elts(TRIAG*id+LORT+1)
       end if
    end if
  end subroutine midpt
  
end module init_mod
