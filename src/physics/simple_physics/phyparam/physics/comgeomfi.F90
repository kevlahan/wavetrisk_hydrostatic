MODULE comgeomfi
  IMPLICIT NONE
  SAVE

  REAL, ALLOCATABLE :: long(:), lati(:), sinlon(:), coslon(:), sinlat(:), coslat(:)
  INTEGER :: ngridmax, nlayermx, nsoilmx
  !$OMP THREADPRIVATE(long,lati,sinlon,coslon,sinlat,coslat,totarea)
  !$OMP THREADPRIVATE(ngridmax,nlayermx,nsoilmx)
  !$acc declare create(sinlon, coslon, sinlat, coslat)

CONTAINS

  SUBROUTINE init_comgeomfi(klon, klev, longitude, latitude, soil_levs)
    INTEGER, INTENT(IN) :: klon, klev
    REAL, INTENT(IN) :: longitude(klon), latitude(klon) ! in radians
    INTEGER, INTENT(IN), OPTIONAL :: soil_levs ! number of soil levels
    
    ! Check if soil_levs is present and set number soil layers nsoilmx (if not set default)
    if (present(soil_levs)) then
       if (soil_levs == 0 .or. soil_levs == 1) then
          ! Make sure soil_levs is not 0 or 1 (soil model not set up for this)
          PRINT*, 'STOP!! the soil model can not be set to have ', soil_levs, ' layer(s)!'
          STOP
       end if 
       nsoilmx = soil_levs
    else
        nsoilmx = 10
    end if
    
    ngridmax=klon
    nlayermx=klev
    allocate(long(klon))
    allocate(lati(klon))
    allocate(sinlon(klon))
    allocate(coslon(klon))
    allocate(sinlat(klon))
    allocate(coslat(klon))
    long(:) = longitude(:)
    lati(:) = latitude(:)
    sinlat(:)=sin(lati(:))
    coslat(:)=cos(lati(:))
    sinlon(:)=sin(long(:))
    coslon(:)=cos(long(:))
    !$acc update device(sinlon, coslon, sinlat, coslat)
  END SUBROUTINE init_comgeomfi

END MODULE comgeomfi
