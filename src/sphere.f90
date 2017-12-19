module geom_mod
  use param_mod
  use shared_mod
  implicit none

  type Coord
     real(8) x
     real(8) y
     real(8) z
  end type Coord

  type Areas
     real(8), dimension(6) :: part
     real(8) :: hex_inv
  end type Areas

  type(Coord), parameter :: ORIGIN = Coord(0.0_8, 0.0_8, 0.0_8)

  integer, parameter :: N_GLO_DOMAIN = N_ICOSAH_LOZANGE*N_SUB_DOM

contains
  subroutine init_sphere_mod
    ! if needed in future
  end subroutine init_sphere_mod

  type(Coord) function direction (init, term)
    type(Coord) init
    type(Coord) term
    type(Coord) v

    v = vector(init, term)
    direction = normalize_Coord(v)
  end function direction

  function vec_plus (v1, v2)
    type(Coord) :: vec_plus
    type(Coord) :: v1, v2

    vec_plus = Coord (v1%x+v2%x, v1%y+v2%y, v1%z+v2%z)
  end function vec_plus

  function vec_minus (v1, v2)
    type(Coord) :: vec_minus
    type(Coord) :: v1, v2

    vec_minus = Coord (v1%x-v2%x, v1%y-v2%y, v1%z-v2%z)
  end function vec_minus

  type(Coord) function vec_scale(alpha, v)
    real(8) :: alpha
    type(Coord) :: v

    vec_scale = Coord (alpha*v%x, alpha*v%y, alpha*v%z)
  end function vec_scale

  real(8) function dist(p, q)
    type(Coord) p
    type(Coord) q

    dist = asin(sqrt((p%y*q%z - p%z*q%y)**2 + (p%z*q%x - p%x*q%z)**2 + (p%x*q%y - p%y*q%x)**2)/radius**2)*radius
  end function dist

  type(Coord) function sph2cart(lon, lat)
    real(8) lon
    real(8) lat
    sph2cart = Coord(cos(lon)*cos(lat), sin(lon)*cos(lat), sin(lat))
  end function sph2cart

  type(Coord) function cross(u, v)
    type(Coord) u
    type(Coord) v
    cross = Coord(u%y*v%z - u%z*v%y, u%z*v%x - u%x*v%z, u%x*v%y - u%y*v%x)
  end function cross

  type(Coord) function project_on_sphere(p)
    type(Coord) p
    real(8) nrm
    nrm = sqrt(p%x**2 + p%y**2 + p%z**2)
    p%x = p%x*(radius/nrm)
    p%y = p%y*(radius/nrm)
    p%z = p%z*(radius/nrm)
    project_on_sphere = p
  end function project_on_sphere

  subroutine arc_inters (arc1_no1, arc1_no2, arc2_no1, arc2_no2, inters_pt, does_inters, troubles)
    type(Coord) :: arc1_no1, arc1_no2, arc2_no1, arc2_no2, inters_pt, neg_int_pt, normal1, normal2
    
    real(8) :: inpr
    logical :: does_inters, troubles
    
    inters_pt = arc2_no2
    does_inters = .True.
    troubles = .False.

    if (norm(vector(arc1_no2, arc2_no2)) .lt. eps()) return

    normal1 = cross(arc1_no1, arc1_no2)
    inpr = inner(normal1, arc2_no1)*inner(normal1, arc2_no2)
    if (inpr .gt. 0.0) then
       if (inpr .lt. (eps()*radius**2)**2) troubles = .True.
       does_inters = .False.
       return
    end if

    normal2 = cross(arc2_no1, arc2_no2)
    inpr = inner(normal2, arc1_no1)*inner(normal2, arc1_no2)

    if (inpr .gt. 0.0) then
       if (inpr .lt. (eps()*radius**2)**2) troubles = .True.
       does_inters = .False.
       return
    end if

    inters_pt = project_on_sphere(cross(normal1, normal2))
    call init_Coord(neg_int_pt, -inters_pt%x, -inters_pt%y, -inters_pt%z)

    if (norm(vector(neg_int_pt, arc1_no1)) .lt. norm(vector(inters_pt, &
         arc1_no1))) then
       inters_pt = neg_int_pt
    end if

    does_inters = .True.
  end subroutine arc_inters

  function vector (init, term)
    type(Coord) :: vector
    type(Coord) :: init, term

    vector = Coord(term%x - init%x, term%y - init%y, term%z - init%z)
  end function vector

  function inner (u, v)
    real(8)     :: inner
    type(Coord) :: u, v

    inner = u%x*v%x + u%y*v%y + u%z*v%z
  end function inner

  function triarea (A, B, C)
    real(8) :: triarea
    type(Coord) :: A, B, C

    real(8) :: ab, ac, bc, s, t

    ab = distn(A, B)
    ac = distn(A, C)
    bc = distn(B, C)

    s = (ab + ac + bc)*0.5_8

    t = tan(0.5_8*s)*tan(-0.5_8*ab + 0.5_8*s)*tan(-0.5_8*ac + &
         0.5_8*s)*tan(-0.5_8*bc + 0.5_8*s)

    if (t .lt. 1.0d-64) then
       triarea = 0.0
       return
    end if

    triarea = 4*radius**2*atan(sqrt(t))
  end function triarea

  function distn (p, q)
    real(8)     :: distn
    type(Coord) :: p, q

    real(8) :: sindist

    sindist = (1.0_8/radius)**2*sqrt((p%y*q%z - p%z*q%y)**2 + (p%z*q%x - &
         p%x*q%z)**2 + (p%x*q%y - p%y*q%x)**2)

    if (sindist .gt. 1) then
       distn = asin(1.0_8)
       return
    end if

    distn = asin(sindist)
  end function distn

  subroutine cart2sph(c, lon, lat)
    type(Coord) :: c
    real(8)     :: lat, lon

    lat = asin (c%z/radius)
    lon = atan2 (c%y, c%x)
  end subroutine cart2sph

  function circumcentre(A, B, C)
    type(Coord) :: circumcentre
    type(Coord) :: A, B, C

    type(Coord) :: centre

    centre = cross(Coord(A%x - B%x, A%y - B%y, A%z - B%z), Coord(C%x - B%x, &
         C%y - B%y, C%z - B%z))

    if (norm(centre) .lt. eps()) then
       circumcentre = centre
       return
    end if

    circumcentre = project_on_sphere(centre)
  end function circumcentre

  function norm (c)
    real(8)     :: norm
    type(Coord) :: c

    norm = sqrt(c%x**2 + c%y**2 + c%z**2)
  end function norm

  function mid_pt(p, q)
    type(Coord) :: mid_pt
    type(Coord) :: p, q

    mid_pt = project_on_sphere(Coord(p%x + q%x, p%y + q%y, p%z + q%z))
  end function mid_pt

  function normalize_Coord (self)
    type(Coord) :: normalize_Coord
    type(Coord) :: self
    
    real(8) :: nrm

    nrm = sqrt(self%x**2 + self%y**2 + self%z**2)
    normalize_Coord = Coord(self%x/nrm, self%y/nrm, self%z/nrm)
  end function normalize_Coord

  subroutine init_Coord (self, x, y, z)
    type(Coord) :: self
    real(8)     :: x, y, z

    self%x = x
    self%y = y
    self%z = z
  end subroutine init_Coord

  subroutine init_Areas (self, centre, corners, midpts)
    type(Areas)               :: self
    type(Coord)               :: centre
    type(Coord), dimension(6) :: corners, midpts
    
    integer :: i
    real(8) :: area

    do i = 1, 6
       self%part(i) = triarea(centre, corners(i), midpts(i)) + triarea(centre, corners(i), midpts(modulo(i,6)+1))
    end do
    area = sum(self%part)

    if (area .eq. 0.0) then ! Avoid overflow for unused zero area hexagons (points at origin and 10 lozenge vertices)
       self%hex_inv = 1.0_8
    else
       self%hex_inv = 1.0_8/area
    end if
  end subroutine init_Areas

  function proj_vel (vel_fun, ep1, ep2)
    ! Finds velocity in direction from points ep1 to ep2 at mid-point of this vector
    ! given a function for zonal u and meridional v velocities as a function of longitude and latitude
    real(8) :: proj_vel
    external    :: vel_fun
    type(Coord) :: ep1, ep2
    
    type(Coord)           :: co, e_zonal, e_merid, vel
    real(8)               :: lon, lat, u_zonal, v_merid

    co = mid_pt(ep1, ep2)

    ! Find longitude and latitude coordinates of point co
    call cart2sph(co, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),             0.0_8) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    ! Function returning zonal and meridional velocities given longitude and latitude
    call vel_fun(lon, lat, u_zonal, v_merid)

    ! Velocity vector in Cartesian coordinates
    vel = vec_plus(vec_scale(u_zonal,e_zonal), vec_scale(v_merid,e_merid))

    ! Project velocity vector on direction given by points ep1, ep2
    proj_vel = inner(direction(ep1, ep2), vel)
  end function proj_vel

  function proj_vel_eta (vel_fun, ep1, ep2, eta_z)
    !extention of proj_vel that allows for another parameter eta_z to be passed in vel_fun
    real(8)     :: proj_vel_eta
    external    :: vel_fun
    type(Coord) :: co, ep1, ep2
    real(8)     :: eta_z

    real(8)               :: lon, lat, u, v
    real(8), dimension(3) :: e_lat, e_lon, vel

    co = mid_pt(ep1, ep2)

    call cart2sph(co, lon, lat)

    e_lat = (/-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)/)
    e_lon = (/-sin(lon), cos(lon), 0.0_8/)

    call vel_fun(lon, lat, u, v, eta_z)

    vel = e_lat*v + e_lon*u
    proj_vel_eta = inner(direction(ep1, ep2), Coord(vel(1), vel(2), vel(3)))
  end function proj_vel_eta

end module geom_mod
