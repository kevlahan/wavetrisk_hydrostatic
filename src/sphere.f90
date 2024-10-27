module geom_mod
  use param_mod
  use shared_mod
  use coord_arithmetic_mod
  implicit none
contains
  type(Coord) function direction (init, term)
    implicit none
    type(Coord) :: init, term, v

    v = vector (init, term)
    direction = normalize_Coord (v)
  end function direction

  real(8) function dist (p, q)
    ! Geodesic distance between points on the sphere with coordinates p and q
    implicit none
    type(Coord) :: p, q

    dist = radius * asin (sqrt ((p%y * q%z - p%z * q%y)**2 + (p%z * q%x - p%x * q%z)**2 + (p%x * q%y - p%y * q%x)**2)/radius**2)
  end function dist

  real(8) function dist_sph (lon1, lat1, lon2, lat2)
    ! Distance between points on the sphere angular coordinates (lat1, lon1) and (lat2, lon2)
    implicit none

    real(8) :: lat1, lat2, lon1, lon2

    type(Coord) :: x1, x2
    
    x1 = sph2cart (lon1, lat1)
    x2 = sph2cart (lon2, lat2)

    dist_sph = dist (x1, x2)
  end function dist_sph

  real(8) function geodesic (p, q)
    ! Great circle (minimum) distance between points with coordinates p and q
    implicit none
    type (Coord) :: p, q

    real(8) :: lat1, lat2, lon1, lon2

    call cart2sph (p, lon1, lat1)
    call cart2sph (q, lon2, lat2)
    geodesic = radius * acos (sin (lat1) * sin (lat2) + cos (lat1) * cos (lat2) * cos (lon2-lon1))
  end function geodesic

  subroutine cart2sph (c, lon, lat)
    ! Angular coordinates (in radians) of a point with coordinates c on the sphere
    implicit none
    type(Coord) :: c
    real(8)     :: lat, lon

    lat = asin (c%z/radius)
    lon = atan2 (c%y, c%x)
  end subroutine cart2sph

  type(Coord) function sph2cart (lon, lat)
    ! Cartesian coordinates of point with longitude lon and latitude lat on the sphere
    implicit none
    real(8) :: lon, lat
    
    sph2cart = radius * Coord (cos(lon)*cos(lat), sin(lon)*cos(lat), sin(lat))
  end function sph2cart

  type(Coord) function project_on_sphere (p)
    implicit none
    type(Coord) :: p
    real(8)     :: nrm
    
    nrm = sqrt (p%x**2 + p%y**2 + p%z**2)
    p%x = p%x*(radius/nrm)
    p%y = p%y*(radius/nrm)
    p%z = p%z*(radius/nrm)
    project_on_sphere = p
  end function project_on_sphere

  subroutine arc_inters (arc1_no1, arc1_no2, arc2_no1, arc2_no2, inters_pt, does_inters, troubles)
    implicit none
    type(Coord) :: arc1_no1, arc1_no2, arc2_no1, arc2_no2, inters_pt, neg_int_pt, normal1, normal2
    
    real(8) :: inpr
    logical :: does_inters, troubles
    
    inters_pt = arc2_no2
    does_inters = .true.
    troubles = .false.

    if (norm(vector(arc1_no2, arc2_no2)) < eps()) return

    normal1 = cross (arc1_no1, arc1_no2)
    inpr = inner (normal1, arc2_no1) * inner(normal1, arc2_no2)
    if (inpr > 0d0) then
       if (inpr < (eps()*radius**2)**2) troubles = .true.
       does_inters = .false.
       return
    end if

    normal2 = cross (arc2_no1, arc2_no2)
    inpr = inner (normal2, arc1_no1) * inner (normal2, arc1_no2)

    if (inpr > 0d0) then
       if (inpr < (eps()*radius**2)**2) troubles = .true.
       does_inters = .false.
       return
    end if

    inters_pt = project_on_sphere (cross(normal1, normal2))
    call init_Coord (neg_int_pt, -inters_pt%x, -inters_pt%y, -inters_pt%z)

    if (norm(vector(neg_int_pt, arc1_no1)) < norm (vector (inters_pt, arc1_no1))) then
       inters_pt = neg_int_pt
    end if

    does_inters = .true.
  end subroutine arc_inters

  type(Coord) function vector (init, term)
    implicit none
    type(Coord) :: init, term

    vector = Coord (term%x - init%x, term%y - init%y, term%z - init%z)
  end function vector
  
  real(8) function inner (u, v)
    implicit none
    type(Coord) :: u, v

    inner = u%x*v%x + u%y*v%y + u%z*v%z
  end function inner

  type(Coord) function cross(u, v)
    implicit none
    type(Coord) :: u, v

    cross = Coord (u%y*v%z - u%z*v%y, u%z*v%x - u%x*v%z, u%x*v%y - u%y*v%x)
  end function cross

  real(8) function triarea (A, B, C)
    implicit none
    type(Coord) :: A, B, C

    real(8) :: ab, ac, bc, s, t

    ab = distn (A, B)
    ac = distn (A, C)
    bc = distn (B, C)

    s = (ab + ac + bc)/2

    t = tan(0.5d0*s) * tan ((s-ab)/2d0) * tan ((s-ac)/2d0) * tan ((s-bc)/2d0) 

    if (t < 1d-64) then
       triarea = 0d0
       return
    end if

    triarea = 4d0*radius**2 * atan (sqrt (t))
  end function triarea

  real(8) function distn (p, q)
    implicit none
    type(Coord) :: p, q

    real(8) :: sindist

    sindist = (1d0/radius)**2 * sqrt ((p%y*q%z - p%z*q%y)**2 + (p%z*q%x - p%x*q%z)**2 + (p%x*q%y - p%y*q%x)**2)

    if (sindist > 1d0) then
       distn = asin (1d0)
       return
    end if
    distn = asin (sindist)
  end function distn

  type(Coord) function circumcentre (A, B, C)
    implicit none
    type(Coord) :: A, B, C

    type(Coord) :: centre

    centre = cross (Coord(A%x - B%x, A%y - B%y, A%z - B%z), Coord(C%x - B%x, C%y - B%y, C%z - B%z))

    if (norm(centre) < eps()) then
       circumcentre = centre
       return
    end if
    circumcentre = project_on_sphere (centre)
  end function circumcentre
  
  type(Coord) function centroid (points, n)
    ! Computes centroid of polygon given coordinates for its n nodes
    ! Simple area-weighted average (second-order accurate, stable)
    implicit none
    integer                   :: n
    type(Coord), dimension(n) :: points

    integer     :: i, j
    type(Coord) :: cc
    real(8)     :: area

    ! Arithmetic mean used as center
    cc = points(1)
    do i = 2, n
       cc = cc + points(i)
    end do
    cc = cc / 6d0
    
    centroid = ORIGIN
    do i = 1, n
       j = mod(i,n)+1
       area = triarea (cc, points(i), points(j))
       centroid = centroid + area * (cc + points(i) + points(j))
    end do
    centroid = project_on_sphere (centroid/6d0)
  end function centroid

  real(8) function norm (c)
    implicit none
    type(Coord) :: c

    norm = sqrt (c%x**2 + c%y**2 + c%z**2)
  end function norm

  type(Coord) function mid_pt(p, q)
    implicit none
    type(Coord) :: p, q

    mid_pt = project_on_sphere (Coord (p%x + q%x, p%y + q%y, p%z + q%z))
  end function mid_pt

  type(Coord) function normalize_Coord (self)
    implicit none
    type(Coord) :: self
    
    real(8) :: nrm

    nrm = sqrt (self%x**2 + self%y**2 + self%z**2)
    if(nrm >= eps()) then
       normalize_Coord = Coord (self%x/nrm, self%y/nrm, self%z/nrm)
    else
       normalize_Coord = ORIGIN
    end if
  end function normalize_Coord

  subroutine init_Coord (self, x, y, z)
    implicit none
    type(Coord) :: self
    real(8)     :: x, y, z

    self%x = x
    self%y = y
    self%z = z
  end subroutine init_Coord

  subroutine init_Areas (self, centre, corners, midpts)
    implicit none
    type(Areas)               :: self
    type(Coord)               :: centre
    type(Coord), dimension(6) :: corners, midpts
    
    integer :: i

    do i = 1, 6
       self%part(i) = triarea (centre, corners(i), midpts(i)) + triarea (centre, corners(i), midpts(modulo(i,6)+1))
    end do
    self%hex_inv = 1d0
    if (sum(self%part) /= 0d0) self%hex_inv = 1d0 / sum (self%part)
  end subroutine init_Areas

  subroutine wrap_lonlat (lat, lon)
    ! Wraps longitude and latitude onto [-pi,pi] and [-pi/2,pi/2]
    implicit none
    real(8) :: lat, lon

    if (lat > MATH_PI/2d0) then
       lat =  MATH_PI/2d0 - mod (lat, MATH_PI/2d0)
       lon = lon + MATH_PI
    elseif (lat < -MATH_PI/2d0) then 
       lat = -MATH_PI/2d0 - mod (lat, MATH_PI/2d0)
       lon = lon + MATH_PI
    end if

    if (lon > MATH_PI) then
       lon = -MATH_PI + mod (lon, MATH_PI)
    elseif (lon < - MATH_PI) then
       lon =  MATH_PI + mod (lon, MATH_PI)
    end if
  end subroutine wrap_lonlat

  real(8) function proj_vel (vel_fun, ep1, ep2)
    ! Finds velocity in direction from points ep1 to ep2 at mid-point of this vector
    ! given a function for zonal u and meridional v velocities as a function of longitude and latitude
    implicit none
    external    :: vel_fun
    type(Coord) :: ep1, ep2
    
    type(Coord) :: co, e_zonal, e_merid, vel
    real(8)     :: lon, lat, u_zonal, v_merid

    co = mid_pt (ep1, ep2)

    ! Find longitude and latitude coordinates of point co
    call cart2sph(co, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),               0d0) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    ! Function returning zonal and meridional velocities given longitude and latitude
    call vel_fun (lon, lat, u_zonal, v_merid)

    ! Velocity vector in Cartesian coordinates
    vel = u_zonal * e_zonal + v_merid * e_merid

    ! Project velocity vector on direction given by points ep1, ep2
    proj_vel = inner (direction (ep1, ep2), vel)
  end function proj_vel

  real(8) function hex_area_avg (l)
    ! Average are of hexagonal cells at level l.
    ! (area of sphere divided by number of hexagonal/penatagonal cells)
    implicit none
    integer, intent(in) :: l

    hex_area_avg = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**dble(l) + 2d0)
  end function hex_area_avg
end module geom_mod



