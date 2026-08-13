module geom_mod
  use kind_mod,   only : dp
  use shared_mod, only : Areas, Coord, eps, MATH_PI, ORIGIN, radius
  
  use coord_arithmetic_mod
  
  implicit none

  private
  public :: arc_intersect_test, cart2sph, centroid, circumcentre, cross, direction, dist, dist_sph
  public :: geodesic, triarea, init_Areas, init_Coord, inner, min_dist
  public :: norm, normalize_coord, mid_pt, number_hex, proj_vel, proj_xz_plane, project_on_sphere, sph2cart, vector, wrap_lonlat

  
  abstract interface
     
     subroutine vel_lonlat (lon, lat, u, v)
       use kind_mod, only : dp
       implicit none
       real(dp), intent(in)  :: lon, lat
       real(dp), intent(out) :: u, v
     end subroutine vel_lonlat
     
  end interface

  
contains

  
  function direction (init, term) result(val)
    implicit none
    type(Coord), intent(in) :: init, term
    type(Coord)             :: val
    
    type(Coord) :: v

    v = vector (init, term)
    val = normalize_Coord (v)
  end function direction


  function dist (p, q) result(val)
    ! Geodesic distance between points on the sphere with coordinates p and q
    implicit none
    type(Coord), intent(in) :: p, q
    real(dp)                :: val

    val = radius * asin (sqrt ((p%y * q%z - p%z * q%y)**2 + (p%z * q%x - p%x * q%z)**2 + (p%x * q%y - p%y * q%x)**2)/radius**2)
  end function dist

  
  subroutine min_dist (p, q, dmin, imin)
    ! Minimum distance between a point p and a vector of points q in R3
    implicit none
    type(Coord), intent(in)  :: p
    type(Coord), intent(in)  :: q(:)
    integer,     intent(out) :: imin
    real(dp),    intent(out) :: dmin

    integer               :: i, n
    real(dp), allocatable :: diff_pq(:)

    n = size (q)
    allocate (diff_pq(n))

    diff_pq = 0.0_dp
    do i = 1, n
       diff_pq(i) = sqrt ((p%x - q(i)%x)**2 + (p%y - q(i)%y)**2 + (p%z - q(i)%z)**2)
    end do

    dmin = minval (diff_pq)
    imin = minloc (diff_pq, 1)
  end subroutine min_dist

  
  function dist_sph (lon1, lat1, lon2, lat2) result(val)
    ! Distance between points on the sphere angular coordinates (lat1, lon1) and (lat2, lon2)
    implicit none
    real(dp), intent(in) :: lat1, lat2, lon1, lon2
    real(dp)             :: val

    type(Coord) :: x1, x2
    
    x1 = sph2cart (lon1, lat1)
    x2 = sph2cart (lon2, lat2)

    val = dist (x1, x2)
  end function dist_sph


  function geodesic (p, q) result(val)
    ! Great circle (minimum) distance between points with coordinates p and q
    implicit none
    type (Coord), intent(in) :: p, q
    real(dp)                 :: val

    real(dp) :: lat1, lat2, lon1, lon2

    call cart2sph (p, lon1, lat1)
    call cart2sph (q, lon2, lat2)
    
    val = radius * acos (sin (lat1) * sin (lat2) + cos (lat1) * cos (lat2) * cos (lon2-lon1))
  end function geodesic

  
  subroutine proj_xz_plane (cin, cout)
    implicit none
    type(Coord), intent(in)  :: cin
    real(dp),    intent(out) :: cout(2)

    if (cin%y > 0) then
       cout = [ cin%x-radius, cin%z ]
    else
       cout = [ cin%x+radius, cin%z ]
    end if
  end subroutine proj_xz_plane

  
  subroutine cart2sph (c, lon, lat)
    ! Angular coordinates (in radians) of a point with coordinates c on the sphere
    implicit none
    type(Coord), intent(in)  :: c
    real(dp),    intent(out) :: lat, lon

    lat = asin (c%z/radius)
    lon = atan2 (c%y, c%x)
  end subroutine cart2sph
  
  
  function sph2cart (lon, lat) result(val)
    ! Cartesian coordinates of point with longitude lon and latitude lat on the sphere
    implicit none
    real(dp), intent(in) :: lon, lat
    type(Coord)          :: val

    val = radius * Coord (cos(lon)*cos(lat), sin(lon)*cos(lat), sin(lat))
  end function sph2cart

  
  function project_on_sphere (p) result(val)
    implicit none
    type(Coord), intent(in) :: p
    type(Coord)             :: val

    real(dp) :: r
    
    r = norm (p)

    if (r < eps(radius)) then
       val = ORIGIN
    else
       val = radius * p / r
    end if
  end function project_on_sphere

  
  subroutine arc_intersect_test (arc1_node1, arc1_node2, arc2_node1, arc2_node2, &
       intersection_pt, intersects, degenerate)
    ! Determines whether two great-circle arcs on sphere intersect.
    !
    ! Uses signed orientation (same-side) tests with respect to planes of two great circles.
    !
    ! If arcs intersect, corresponding spherical intersection point is computed from cross product
    ! of two great-circle normals and projected onto sphere.
    !
    ! Near-degenerate/tangent configurations are flagged as `degenerate` to be corrected.
    implicit none

    type(Coord), intent(in)  :: arc1_node1, arc1_node2
    type(Coord), intent(in)  :: arc2_node1, arc2_node2

    type(Coord), intent(out) :: intersection_pt ! intersection point 
    logical,     intent(out) :: intersects      ! flags intersections
    logical,     intent(out) :: degenerate      ! flags degenerate and near-zero cases
    
    real(dp)    :: side_product ! Product of signed distances/orientations relative to great-circle plane with normal given by
                                ! cross product of arc end points.
                                ! Measures whether two endpoints lie on same side or opposite sides of great-circle plane:
                                ! > 0 same side
                                ! < 0 opposite side
                                ! ~ 0 one point nearly on plane
    
    real(dp)    :: dneg, dpos, tol, nx
    
    type(Coord) :: neg_intersection_pt, normal1, normal2, x
    
    intersects = .true.
    degenerate = .false.

    ! Coincident arc endpoints case
    tol = 100.0_dp * eps(radius)
    if (norm(vector(arc1_node2, arc2_node2)) < tol) then
       intersection_pt = arc2_node2
       return
    end if

    ! Empirical scale-aware tolerance for side_product (smaller than worst case radius**6)
    tol = 100.0_dp * eps(radius**4)

    ! Test end points of arc2 relative to great-circle plane of arc1
    normal1 = cross(arc1_node1, arc1_node2)
    side_product = inner(normal1, arc2_node1) &
                 * inner(normal1, arc2_node2)
    if (side_product > tol) then
       intersects = .false.
       return
    elseif (abs(side_product) <= tol) then
       degenerate = .true.
    end if

    ! Test end points of arc1 relative to great-circle plane of arc2
    normal2 = cross(arc2_node1, arc2_node2)
    side_product = inner(normal2, arc1_node1) &
                 * inner(normal2, arc1_node2)
    if (side_product > tol) then
       intersects = .false.
       return
    elseif (abs(side_product) <= tol) then
       degenerate = .true.
    end if

    ! Intersection of the two great circles
    x  = cross(normal1, normal2)
    nx = norm(x)
    tol = 100.0_dp * eps(radius**4)
    if (nx <= tol) then ! nearly parallel great circles: intersection direction is ill-conditioned
       intersection_pt = arc2_node2
       intersects = .false.
       degenerate = .true.
       return
    end if

    ! Choose antipodal intersection point closer to arc1_node1
    intersection_pt = project_on_sphere (x)
    neg_intersection_pt = (-1.0_dp) * intersection_pt

    ! Distances from two antipodal intersection candidates to reference arc endpoint.
    dneg = norm(vector(neg_intersection_pt, arc1_node1))
    dpos = norm(vector(intersection_pt,     arc1_node1))

    tol = 100.0_dp * eps(radius)
    if (dneg < dpos - tol) then
       intersection_pt = neg_intersection_pt
    elseif (abs(dneg - dpos) <= tol) then
       degenerate = .true.
    end if

    intersects = .true.

  end subroutine arc_intersect_test

  
  function vector (init, term) result(val)
    implicit none
    type(Coord), intent(in) :: init, term
    type(Coord)             :: val

    val = Coord (term%x - init%x, term%y - init%y, term%z - init%z)
  end function vector

  
  function inner (u, v) result(val)
    implicit none
    type(Coord), intent(in) :: u, v
    real(dp)                :: val

    val = u%x*v%x + u%y*v%y + u%z*v%z
  end function inner

  
  function cross (u, v) result(val)
    implicit none
    type(Coord), intent(in) :: u, v
    type(Coord)             :: val

    val = Coord (u%y*v%z - u%z*v%y, u%z*v%x - u%x*v%z, u%x*v%y - u%y*v%x)
  end function cross

  
  function triarea (A, B, C) result(val)
    implicit none
    type(Coord), intent(in) :: A, B, C
    real(dp)                :: val

    real(dp) :: ab, ac, bc, s, t

    ab = distn (A, B)
    ac = distn (A, C)
    bc = distn (B, C)

    s = (ab + ac + bc)/2

    t = tan(0.5*s) * tan ((s-ab)/2) * tan ((s-ac)/2) * tan ((s-bc)/2) 

    if (t < 1e-64_dp) then
       val = 0.0_dp
       return
    end if

    val = 4*radius**2 * atan (sqrt (t))
  end function triarea

  
  function distn (p, q) result(val)
    implicit none
    type(Coord), intent(in) :: p, q
    real(dp)                :: val

    real(dp) :: sindist

    sindist = (1/radius)**2 * sqrt ((p%y*q%z - p%z*q%y)**2 + (p%z*q%x - p%x*q%z)**2 + (p%x*q%y - p%y*q%x)**2)

    if (sindist > 1.0_dp) then
       val = asin (1.0_dp)
       return
    end if
    val = asin (sindist)
  end function distn


  function circumcentre (A, B, C) result(val)
    implicit none
    type(Coord), intent(in) :: A, B, C
    type(Coord)             :: val

    type(Coord) :: centre

    centre = cross (Coord(A%x - B%x, A%y - B%y, A%z - B%z), Coord(C%x - B%x, C%y - B%y, C%z - B%z))

    if (norm(centre) < eps(radius)) then
       val = centre
       return
    end if
    val = project_on_sphere (centre)
  end function circumcentre

  
 function centroid (points, n) result(val)
    ! Computes centroid of polygon given coordinates for its n nodes
    ! Simple area-weighted average (second-order accurate, stable)
    implicit none
    integer,                   intent(in) :: n
    type(Coord), dimension(n), intent(in) :: points
    type(Coord)                           :: val

    integer     :: i, j
    type(Coord) :: cc
    real(dp)    :: area

    ! Arithmetic mean used as center
    cc = points(1)
    do i = 2, n
       cc = cc + points(i)
    end do
    cc = cc / 6.0_dp
    
    val = ORIGIN
    do i = 1, n
       j = mod(i,n)+1
       area = triarea (cc, points(i), points(j))
       val = val + area * (cc + points(i) + points(j))
    end do
    val = project_on_sphere (val/6.0_dp)
  end function centroid


  function norm (c) result(val)
    implicit none
    type(Coord), intent(in) :: c
    real(dp)                :: val

    val = sqrt (c%x**2 + c%y**2 + c%z**2)
  end function norm

  
  function mid_pt (p, q) result(val)
    implicit none
    type(Coord), intent(in) :: p, q
    type(Coord)             :: val

    val = project_on_sphere (Coord (p%x + q%x, p%y + q%y, p%z + q%z))
  end function mid_pt

  
  function normalize_Coord (self) result(val)
    implicit none
    type(Coord), intent(in) :: self
    type(Coord)             :: val
    
    real(dp) :: nrm

    nrm = sqrt (self%x**2 + self%y**2 + self%z**2) + eps (radius)
    if(nrm >= eps(radius)) then
       val = Coord (self%x/nrm, self%y/nrm, self%z/nrm)
    else
       val = ORIGIN
    end if
  end function normalize_Coord

  
  subroutine init_Coord (self, x, y, z)
    implicit none
    real(dp),    intent(in)  :: x, y, z
    type(Coord), intent(out) :: self
    
    self%x = x
    self%y = y
    self%z = z
  end subroutine init_Coord

  
  subroutine init_Areas (self, centre, corners, midpts)
    implicit none
    type(Coord), intent(in)  :: centre
    type(Coord), intent(in)  :: corners(6), midpts(6)
    type(Areas), intent(out) :: self
    
    integer :: i

    do i = 1, 6
       self%part(i) = triarea (centre, corners(i), midpts(i)) + triarea (centre, corners(i), midpts(modulo(i,6)+1))
    end do
    self%hex_inv = 1.0_dp
    if (sum(self%part) > eps(radius)**2) self%hex_inv = 1.0_dp / sum (self%part)
  end subroutine init_Areas

  
  subroutine wrap_lonlat (lat, lon)
    ! Wraps longitude and latitude onto [-pi,pi] and [-pi/2,pi/2]
    implicit none
    real(dp), intent(inout) :: lat, lon

    if (lat > MATH_PI/2) then
       lat =  MATH_PI/2 - mod (lat, MATH_PI/2)
       lon = lon + MATH_PI
    elseif (lat < -MATH_PI/2) then 
       lat = -MATH_PI/2 - mod (lat, MATH_PI/2)
       lon = lon + MATH_PI
    end if

    if (lon > MATH_PI) then
       lon = -MATH_PI + mod (lon, MATH_PI)
    elseif (lon < - MATH_PI) then
       lon =  MATH_PI + mod (lon, MATH_PI)
    end if
  end subroutine wrap_lonlat

  
  function proj_vel (velocity, ep1, ep2) result(val)
    ! Finds velocity in direction from points ep1 to ep2 at mid-point of this vector
    ! given a function for zonal u and meridional v velocities as a function of longitude and latitude
    implicit none
    type(Coord), intent(in) :: ep1, ep2
    real(dp)                :: val
    
    type(Coord)           :: co, e_zonal, e_merid, vel
    real(dp)              :: lon, lat, u_zonal, v_merid
    procedure(vel_lonlat) :: velocity

    co = mid_pt (ep1, ep2)

    ! Find longitude and latitude coordinates of point co
    call cart2sph(co, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),            0.0_dp) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    ! Function returning zonal and meridional velocities given longitude and latitude
    call velocity (lon, lat, u_zonal, v_merid)
    
    ! Velocity vector in Cartesian coordinates
    vel = u_zonal * e_zonal + v_merid * e_merid

    ! Project velocity vector on direction given by points ep1, ep2
    val = inner (direction (ep1, ep2), vel)
  end function proj_vel

  
  function number_hex (l) result(val)
    ! Number of hexagonal/pentagonal cells for level l
    integer, intent(in) :: l
    integer             :: val

    val = 10 * 4**l + 2
  end function number_hex

  
end module geom_mod



