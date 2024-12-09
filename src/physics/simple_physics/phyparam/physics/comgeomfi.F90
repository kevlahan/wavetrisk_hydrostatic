module comgeomfi
  implicit none
  save
  integer                         :: ngridmax, nlayermx, nsoilmx
  real, dimension(:), allocatable :: long, lati, SinLon, CosLon, SinLat, CosLat
contains
  subroutine init_comgeomfi (klon, klev, longitude, latitude)
    integer,               intent(in) :: klon, klev
    real, dimension(klon), intent(in) :: longitude, latitude ! in radians

    ngridmax = klon
    nlayermx = klev

    allocate (long(klon), lati(klon))
    allocate (SinLon(klon), CosLon(klon), SinLat(klon), CosLat(klon))

    long = longitude
    lati = latitude

    SinLat = sin (lati)
    CosLat = cos (lati)
    SinLon = sin (long)
    CosLon = cos (long)
  end subroutine init_comgeomfi
end module comgeomfi
