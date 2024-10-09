module comgeomfi
  implicit none
  save
  integer                         :: ngridmax, nlayermx, nsoilmx
  real, dimension(:), allocatable :: long, lati, sinlon, coslon, sinlat, coslat
contains
  subroutine init_comgeomfi (klon, klev, longitude, latitude)
    integer,               intent(in) :: klon, klev
    real, dimension(klon), intent(in) :: longitude, latitude ! in radians

    ngridmax = klon
    nlayermx = klev

    allocate (long(klon), lati(klon))
    allocate (sinlon(klon), coslon(klon), sinlat(klon), coslat(klon))

    long  =  longitude
    lati  =  latitude

    sinlat = sin (lati)
    coslat = cos (lati)
    sinlon = sin (long)
    coslon = cos (long)
  end subroutine init_comgeomfi
end module comgeomfi
