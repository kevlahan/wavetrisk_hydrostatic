! $Id: physiq.F 1565 2011-08-31 12:53:29Z jghattas $
MODULE physiq_mod
  USE dimphy, only : klon,klev

  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL, PARAMETER :: oneday = 86400. ! hard-coded 

  INTEGER :: itau=0 ! counter to count number of calls to physics
!$OMP THREADPRIVATE(itau)
  INTEGER :: nid_hist ! NetCDF ID of output file
!$OMP THREADPRIVATE(nid_hist)

  INTEGER :: iwrite_phys ! output every iwrite_phys physics step

  PUBLIC :: physiq

CONTAINS

  SUBROUTINE physiq (nlon,nlev, &
       &            debut,lafin,pdtphys, &
       &            paprs,pplay,pphi,pphis,presnivs, &
       &            u,v,t,qx, &
       &            flxmass_w, &
       &            d_u, d_v, d_t, d_qx, d_ps)
    
    USE infotrac_phy, only : nqtot
    USE iophys,       ONLY : iophys_ini
    USE phyparam_mod, ONLY : phyparam
    !
    ! Routine argument:
    !
    integer,intent(in) :: nlon ! number of atmospheric colums
    integer,intent(in) :: nlev ! number of vertical levels (should be =klev)
    logical,intent(in) :: debut ! signals first call to physics
    logical,intent(in) :: lafin ! signals last call to physics
    real,intent(in) :: pdtphys ! physics time step (s)
    real,intent(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
    real,intent(in) :: pplay(klon,klev) ! mid-layer pressure (Pa)
    real,intent(in) :: pphi(klon,klev) ! geopotential at mid-layer
    real,intent(in) :: pphis(klon) ! surface geopotential
    real,intent(in) :: presnivs(klev) ! pseudo-pressure (Pa) of mid-layers
    real,intent(in) :: u(klon,klev) ! eastward zonal wind (m/s)
    real,intent(in) :: v(klon,klev) ! northward meridional wind (m/s)
    real,intent(in) :: t(klon,klev) ! temperature (K)
    real,intent(in) :: qx(klon,klev,nqtot) ! tracers (.../kg_air)
    real,intent(in) :: flxmass_w(klon,klev) ! vertical mass flux
    real,intent(out) :: d_u(klon,klev) ! physics tendency on u (m/s/s)
    real,intent(out) :: d_v(klon,klev) ! physics tendency on v (m/s/s)
    real,intent(out) :: d_t(klon,klev) ! physics tendency on t (K/s)
    real,intent(out) :: d_qx(klon,klev,nqtot) ! physics tendency on tracers
    real,intent(out) :: d_ps(klon) ! physics tendency on surface pressure
    
    logical, save :: first=.true.
    !$OMP THREADPRIVATE(first)
    
    ! For I/Os
    REAL,SAVE :: rjourvrai=0.
    REAL,SAVE :: gmtime=0.
    !$OMP THREADPRIVATE(rjourvrai,gmtime)
        
    ! initializations
    if (debut) CALL init_physiq_mod(pdtphys, presnivs) ! Things to do only for the first call to physics 
    
    ! increment local time counter itau
    itau=itau+1
    gmtime=gmtime+pdtphys/oneday
    IF (gmtime>=1.) THEN
       gmtime=gmtime-1.
       rjourvrai=rjourvrai+1.
    ENDIF
    
    IF(debut) CALL iophys_ini('phys.nc    ',presnivs) ! calls iotd_ini

    CALL phyparam(klon,klev,       &
         debut,lafin,              &
         rjourvrai,gmtime,pdtphys, &
         paprs,pplay,pphi,         &
         u,v,t,                    &
         d_u,d_v,d_t,d_ps) 

    d_qx(:,:,:)=0.

    IF(lafin) THEN
       call iotd_fin
       PRINT*,'Ecriture du fichier de reinitialiastion de la physique'
!       write(75)  tsurf,tsoil FIXME
    ENDIF

    
    print*,'PHYDEV: itau=',itau

    CALL write_hist(paprs, t, u, v)      ! write some outputs

    ! if lastcall, then it is time to write "restartphy.nc" file
    IF (lafin) CALL phyredem("restartphy.nc")
    
  END SUBROUTINE physiq

  SUBROUTINE init_physiq_mod(dtime, presnivs)
    USE iophy, only : histbeg_phy
    USE phys_state_var_mod, only : phys_state_var_init
    USE ioipsl, only : getin,histvert,histdef,histend,ymds2ju
    USE mod_phys_lmdz_para, only : jj_nb
    USE mod_grid_phy_lmdz, ONLY: nbp_lon,nbp_lat
    REAL, INTENT(IN) :: dtime
    REAL, INTENT(IN) :: presnivs(klev) ! pseudo-pressure (Pa) of mid-layers

    ! For I/Os
    INTEGER, PARAMETER :: itau0=0, annee0=1979, month=1, dayref=1
    REAL, PARAMETER :: hour=0.0
    INTEGER :: nid_hori    ! NetCDF ID of horizontal coordinate
    INTEGER :: nid_vert ! NetCDF ID of vertical coordinate
    REAL :: t_ops ! frequency of the IOIPSL operations (eg average over...)
    REAL :: t_wrt ! frequency of the IOIPSL outputs
    REAL :: zjulian

    ! load initial conditions for physics (including the grid)
    call phys_state_var_init() ! some initializations, required before calling phyetat0
    call phyetat0("startphy.nc")
    
    ! Initialize outputs:

    ! NB: getin() is not threadsafe; only one thread should call it.
    !$OMP MASTER
    iwrite_phys=1 !default: output every physics timestep
    call getin("iwrite_phys",iwrite_phys)
    !$OMP END MASTER
    !$OMP BARRIER

    t_ops=dtime*iwrite_phys ! frequency of the IOIPSL operation
    t_wrt=dtime*iwrite_phys ! frequency of the outputs in the file

    ! compute zjulian for annee0=1979 and month=1 dayref=1 and hour=0.0
    CALL ymds2ju(annee0, month, dayref, hour, zjulian)
#ifndef CPP_IOIPSL_NO_OUTPUT
    ! Initialize IOIPSL output file
    call histbeg_phy("histins.nc",itau0,zjulian,dtime,nid_hori,nid_hist)
#endif
    
    !$OMP MASTER
    
#ifndef CPP_IOIPSL_NO_OUTPUT 
    ! IOIPSL
    ! define vertical coordinate
    call histvert(nid_hist,"presnivs","Vertical levels","Pa",klev, &
         presnivs,nid_vert,'down')
    ! define variables which will be written in "histins.nc" file
    call histdef(nid_hist,'temperature','Atmospheric temperature','K', &
         nbp_lon,jj_nb,nid_hori,klev,1,klev,nid_vert,32, &
         'inst(X)',t_ops,t_wrt)
    call histdef(nid_hist,'u','Eastward Zonal Wind','m/s', &
         nbp_lon,jj_nb,nid_hori,klev,1,klev,nid_vert,32, &
         'inst(X)',t_ops,t_wrt)
    call histdef(nid_hist,'v','Northward Meridional Wind','m/s', &
         nbp_lon,jj_nb,nid_hori,klev,1,klev,nid_vert,32, &
         'inst(X)',t_ops,t_wrt)
    call histdef(nid_hist,'ps','Surface Pressure','Pa', &
         nbp_lon,jj_nb,nid_hori,1,1,1,nid_vert,32, &
         'inst(X)',t_ops,t_wrt)
    ! end definition sequence
    call histend(nid_hist)
#endif
    
#ifdef CPP_XIOS
    !XIOS
    ! Declare available vertical axes to be used in output files:    
    CALL wxios_add_vaxis("presnivs", klev, presnivs)
    
    ! Declare calendar and time step
    CALL wxios_set_cal(dtime,"earth_360d",1,1,1,0.0,1,1,1,0.0)
    
    !Finalize the context:
    CALL wxios_closedef()
#endif
    !$OMP END MASTER
    !$OMP BARRIER
  END SUBROUTINE init_physiq_mod
  
  SUBROUTINE write_hist(paprs, t, u, v)
    USE iophy, only : histwrite_phy
#ifdef CPP_XIOS
    USE xios, ONLY: xios_update_calendar
    USE wxios, only: wxios_add_vaxis, wxios_set_cal, wxios_closedef
    USE iophy, ONLY: histwrite_phy
#endif
    REAL, INTENT(IN) :: paprs(klon,klev+1), & ! interlayer pressure (Pa)
         &              t(klon,klev),       & ! temperature (K)
         &              u(klon,klev),       & ! eastward zonal wind (m/s)
         &              v(klon,klev)          ! northward meridional wind (m/s)
    ! output using IOIPSL
#ifndef CPP_IOIPSL_NO_OUTPUT 
    if (modulo(itau,iwrite_phys)==0) then
       call histwrite_phy(nid_hist,.false.,"temperature",itau,t)
       call histwrite_phy(nid_hist,.false.,"u",itau,u)
       call histwrite_phy(nid_hist,.false.,"v",itau,v)
       call histwrite_phy(nid_hist,.false.,"ps",itau,paprs(:,1))
    endif
#endif
    
    ! output using XIOS
#ifdef CPP_XIOS
    !$OMP MASTER
    !Increment XIOS time
    CALL xios_update_calendar(itau)
    !$OMP END MASTER
    !$OMP BARRIER
    
    !Send fields to XIOS: (NB these fields must also be defined as
    ! <field id="..." /> in iodef.xml to be correctly used
    CALL histwrite_phy("temperature",t)
    CALL histwrite_phy("temp_newton",temp_newton)
    CALL histwrite_phy("u",u)
    CALL histwrite_phy("v",v)
    CALL histwrite_phy("ps",paprs(:,1))
#endif
    
  END SUBROUTINE write_hist
  
END MODULE physiq_mod
