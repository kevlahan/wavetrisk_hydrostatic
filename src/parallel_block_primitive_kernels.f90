module parallel_block_primitive_kernels_mod
  ! Explicit compact operand kernels. No Domain storage, field pointers, MPI,
  ! callbacks or module numerical state. Scheduling/ownership is separate.
  use kind_mod, only: dp
  implicit none
  private
  integer, parameter :: EDGE=3, TRIAG=2, RT=0, DG=1, UP=2, LORT=0, UPLT=1
  integer, parameter :: NORTH=1,EAST=2,SOUTH=3,WEST=4,NORTHEAST=5,SOUTHWEST=7
  integer, parameter :: IJMINUS=7,IPLUSJMINUS=6,IMINUSJPLUS=8,IJPLUS=5
  public :: primitive_regular, primitive_southwest, primitive_pentagon
  public :: pressure_layer, surface_pressure, scalar_diffusion_flux, velocity_diffusion_source
contains
  pure real(dp) function interp(a,b) result(v)
    real(dp), intent(in) :: a,b
    v=0.5*(a+b)
  end function

  pure real(dp) function surface_pressure(mass,mean_mass,gravity,top_pressure) result(pressure)
    real(dp), intent(in) :: mass(:),mean_mass(:),gravity,top_pressure
    integer :: k
    pressure=0.0_dp
    do k=1,size(mass)
       pressure=pressure+(mean_mass(k)+mass(k))
    end do
    pressure=gravity*pressure+top_pressure
  end function

  pure subroutine pressure_layer(mass,mean_mass,temp,mean_temp,gravity,cp,p0,kappa, &
       pressure_lower,geopot,pressure,exner,geopot_lower)
    real(dp), intent(in) :: mass,mean_mass,temp,mean_temp,gravity,cp,p0,kappa
    real(dp), intent(inout) :: pressure_lower,geopot
    real(dp), intent(out) :: pressure,exner,geopot_lower
    real(dp) :: rho_dz,rho_dz_theta,p_upper
    rho_dz=mean_mass+mass
    rho_dz_theta=mean_temp+temp
    geopot_lower=geopot
    p_upper=pressure_lower-gravity*rho_dz
    pressure=interp(pressure_lower,p_upper)
    exner=cp*(pressure/p0)**kappa
    geopot=geopot_lower+gravity*kappa*rho_dz_theta*exner/pressure
    pressure_lower=p_upper
  end subroutine

  pure function scalar_diffusion_flux(values,length,pedlen,order,nu) result(flux)
    real(dp), intent(in) :: values(4),length(3),pedlen(3),nu
    integer, intent(in) :: order
    real(dp) :: flux(3),grad(3)
    flux=0.0_dp
    if(order==0) return
    grad(1)=(values(2)-values(1))/length(1)
    grad(2)=(values(1)-values(3))/length(2)
    grad(3)=(values(4)-values(1))/length(3)
    flux=(-1)**order*nu*grad*pedlen
  end function

  pure function velocity_diffusion_source(divergence,vorticity,length,pedlen,div_order,rot_order, &
       div_nu,rot_nu) result(source)
    real(dp), intent(in) :: divergence(4),vorticity(4),length(3),pedlen(3),div_nu,rot_nu
    integer, intent(in) :: div_order,rot_order
    real(dp) :: source(3),grad(3),curl(3)
    source=0.0_dp
    if(div_order/=0) then
       grad(1)=(divergence(2)-divergence(1))/length(1)
       grad(2)=(divergence(1)-divergence(3))/length(2)
       grad(3)=(divergence(4)-divergence(1))/length(3)
       source=(-1)**(div_order-1)*div_nu*grad
    end if
    if(rot_order/=0) then
       ! C lower, C upper, S upper, W lower. Preserve callback expression order.
       curl(1)=(vorticity(1)-vorticity(3))/pedlen(1)
       curl(2)=(vorticity(1)-vorticity(2))/pedlen(2)
       curl(3)=(vorticity(4)-vorticity(2))/pedlen(3)
       source=source-(-1)**(rot_order-1)*rot_nu*curl
    end if
  end function

  pure subroutine primitive_regular(rho,velo,length,pedlen,area,hex_inv,triarea,coriolis, &
       geopot,geopot_lower,qe,vort,kinetic,bernoulli,pv_context,dual)
    real(dp), intent(in) :: rho(0:7),velo(3,0:7),length(3,0:7),pedlen(3,0:7),area(6,0:7)
    real(dp), intent(in) :: hex_inv,triarea(2,0:7),coriolis(2,0:7),geopot,geopot_lower
    real(dp), intent(inout) :: qe(3,0:7),vort(2,0:7)
    real(dp), intent(out) :: kinetic,bernoulli,pv_context(2),dual(3)
    real(dp) :: u_prim_UP,u_dual_UP,u_prim_DG,u_prim_DG_S,u_prim_DG_W,u_dual_DG,u_prim_RT,u_dual_RT
    real(dp) :: u_prim_UP_S,u_dual_UP_S,u_prim_DG_SW,u_dual_DG_SW,u_prim_RT_W,u_dual_RT_W
    real(dp) :: u_prim_UP_E,u_prim_RT_N
    real(dp) :: circ_LORT,circ_UPLT,circ_S_UPLT,circ_W_LORT
    real(dp) :: pv_LORT,pv_UPLT,pv_S_UPLT,pv_W_LORT,pv_mass(4),Phi_k
         u_dual_RT    = velo(RT+1,0)   * pedlen(RT+1,0)
         u_dual_DG    = velo(DG+1,0)   * pedlen(DG+1,0)
         u_dual_UP    = velo(UP+1,0)   * pedlen(UP+1,0)
         u_dual_RT_W  = velo(RT+1,WEST) * pedlen(RT+1,WEST)
         u_dual_DG_SW = velo(DG+1,SOUTHWEST) * pedlen(DG+1,SOUTHWEST)
         u_dual_UP_S  = velo(UP+1,SOUTH) * pedlen(UP+1,SOUTH)

         u_prim_RT    = velo(RT+1,0) * length(RT+1,0)
         u_prim_RT_N  = velo(RT+1,NORTH) * length(RT+1,NORTH)
         u_prim_RT_W  = velo(RT+1,WEST) * length(RT+1,WEST)
         u_prim_DG    = velo(DG+1,0) * length(DG+1,0)
         u_prim_DG_S  = velo(DG+1,SOUTH) * length(DG+1,SOUTH)
         u_prim_DG_SW = velo(DG+1,SOUTHWEST) * length(DG+1,SOUTHWEST)
         u_prim_DG_W  = velo(DG+1,WEST) * length(DG+1,WEST)
         u_prim_UP    = velo(UP+1,0) * length(UP+1,0)
         u_prim_UP_E  = velo(UP+1,EAST) * length(UP+1,EAST)
         u_prim_UP_S  = velo(UP+1,SOUTH) * length(UP+1,SOUTH)

         ! Potential vorticity
         circ_LORT   =   u_prim_RT   + u_prim_UP_E + u_prim_DG
         circ_UPLT   = -(u_prim_DG   + u_prim_UP   + u_prim_RT_N)
         circ_W_LORT =   u_prim_RT_W + u_prim_UP   + u_prim_DG_W
         circ_S_UPLT = -(u_prim_RT   + u_prim_DG_S + u_prim_UP_S)

         pv_mass(1) = &
              rho(0) * area(1,0) + &
              rho(EAST) * area(3,EAST) + &
              rho(NORTHEAST) * area(5,NORTHEAST)

         pv_mass(2) = &
              rho(0) * area(2,0) + &
              rho(NORTHEAST) * area(4,NORTHEAST) + &
              rho(NORTH) * area(6,NORTH)

         pv_mass(3) = &
              rho(WEST) * area(1,WEST) + &
              rho(0) * area(3,0) + &
              rho(NORTH) * area(5,NORTH)

         pv_mass(4) = &
              rho(SOUTH)  * area(2,SOUTH) + &
              rho(EAST)  * area(4,EAST) + &
              rho(0)  * area(6,0)

         pv_LORT   = (coriolis(LORT+1,0) + circ_LORT  ) / pv_mass(1)
         pv_UPLT   = (coriolis(UPLT+1,0) + circ_UPLT  ) / pv_mass(2)
         pv_W_LORT = (coriolis(LORT+1,WEST) + circ_W_LORT) / pv_mass(3)
         pv_S_UPLT = (coriolis(UPLT+1,SOUTH) + circ_S_UPLT) / pv_mass(4)

         qe(RT+1,0) = interp (pv_S_UPLT, pv_LORT)
         qe(DG+1,0) = interp (pv_UPLT,   pv_LORT)
         qe(UP+1,0) = interp (pv_UPLT,   pv_W_LORT)

         ! Vorticity (for rotu velocity diffusion)
         vort(LORT+1,0) = circ_LORT   / triarea(LORT+1,0)
         vort(UPLT+1,0) = circ_UPLT   / triarea(UPLT+1,0)
         vort(LORT+1,WEST) = circ_W_LORT / triarea(LORT+1,WEST)
         vort(UPLT+1,SOUTH) = circ_S_UPLT / triarea(UPLT+1,SOUTH)

         ! Kinetic energy (TRiSK formula)
         kinetic = (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT +  &
                     u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
              * hex_inv/4

         ! Interpolate geopotential from interfaces to level
         Phi_k = interp (geopot, geopot_lower)

         ! Bernoulli function
         bernoulli=kinetic+Phi_k
         pv_context=[pv_W_LORT,pv_S_UPLT]
         dual=[u_dual_RT,u_dual_DG,u_dual_UP]
  end subroutine

  pure subroutine primitive_southwest(rho,velo,length,area,triarea,coriolis,pv_context,qe,vort)
    real(dp), intent(in) :: rho(0:7),velo(3,0:7),length(3,0:7),area(6,0:7)
    real(dp), intent(in) :: triarea(2,0:7),coriolis(2,0:7),pv_context(2)
    real(dp), intent(inout) :: qe(3,0:7),vort(2,0:7)

    real(dp) :: u_prim_UP_S,u_prim_DG_SW,u_prim_RT_W
    real(dp) :: u_prim_RT_SW,u_prim_UP_SW
    real(dp) :: circ_SW_LORT,circ_SW_UPLT
    real(dp) :: pv_S_UPLT,pv_SW_LORT,pv_SW_UPLT,pv_W_LORT,pv_mass(4)
    pv_W_LORT=pv_context(1)
    pv_S_UPLT=pv_context(2)
         u_prim_RT_W  = velo(RT+1,WEST) * length(RT+1,WEST)
         u_prim_DG_SW = velo(DG+1,SOUTHWEST) * length(DG+1,SOUTHWEST)
         u_prim_UP_S  = velo(UP+1,SOUTH) * length(UP+1,SOUTH)
         u_prim_RT_SW = velo(RT+1,SOUTHWEST) * length(RT+1,SOUTHWEST)
         u_prim_UP_SW = velo(UP+1,SOUTHWEST) * length(UP+1,SOUTHWEST)

         ! Potential vorticity
         circ_SW_LORT =   u_prim_RT_SW + u_prim_UP_S  + u_prim_DG_SW
         circ_SW_UPLT = -(u_prim_RT_W  + u_prim_DG_SW + u_prim_UP_SW)

         pv_mass(1) = &
              rho(SOUTHWEST) * area(1,SOUTHWEST) + &
              rho(SOUTH) * area(3,SOUTH) + &
              rho(0) * area(5,0)

         pv_mass(2) = &
              rho(SOUTHWEST) * area(2,SOUTHWEST) + &
              rho(0) * area(4,0) + &
              rho(WEST) * area(6,WEST)

         pv_SW_LORT = (coriolis(LORT+1,SOUTHWEST) + circ_SW_LORT) / pv_mass(1)
         pv_SW_UPLT = (coriolis(UPLT+1,SOUTHWEST) + circ_SW_UPLT) / pv_mass(2)

         qe(RT+1,WEST) = interp (pv_W_LORT , pv_SW_UPLT)
         qe(DG+1,SOUTHWEST) = interp (pv_SW_LORT, pv_SW_UPLT)
         qe(UP+1,SOUTH) = interp (pv_SW_LORT, pv_S_UPLT)

         ! Vorticity (for rotu velocity diffusion)
         vort(LORT+1,SOUTHWEST) = circ_SW_LORT / triarea(LORT+1,SOUTHWEST)
         vort(UPLT+1,SOUTHWEST) = circ_SW_UPLT / triarea(UPLT+1,SOUTHWEST)

  end subroutine

  pure subroutine primitive_pentagon(c,rho,velo,length,area,triarea,coriolis,qe,vort)
    integer, intent(in) :: c
    real(dp), intent(in) :: rho(0:7),velo(3,0:7),length(3,0:7),area(6,0:7)
    real(dp), intent(in) :: triarea(2,0:7),coriolis(2,0:7)
    real(dp), intent(inout) :: qe(3,0:7),vort(2,0:7)
    real(dp) :: u_prim_UP,u_prim_RT
    real(dp) :: u_prim_UP_S,u_prim_DG_SW,u_prim_RT_W
    real(dp) :: u_prim_RT_N,u_prim_RT_SW,u_prim_UP_SW
    real(dp) :: circ_LORT,circ_UPLT,circ_S_UPLT,circ_SW_LORT,circ_SW_UPLT,circ_W_LORT
    real(dp) :: pv_LORT,pv_UPLT,pv_S_UPLT,pv_SW_LORT,pv_SW_UPLT,pv_W_LORT,pv_mass(4)
    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*SOUTHWEST+DG+1) = 0 in this case
    if (c == IJMINUS) then

       u_prim_RT_W  = velo(RT+1,WEST) * length(RT+1,WEST)
       u_prim_RT_SW = velo(RT+1,SOUTHWEST) * length(RT+1,SOUTHWEST)
       u_prim_UP_S  = velo(UP+1,SOUTH) * length(UP+1,SOUTH)

       circ_SW_LORT = u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW
       circ_W_LORT  = vort(LORT+1,WEST) * triarea(LORT+1,WEST)
       circ_S_UPLT  = vort(UPLT+1,SOUTH) * triarea(UPLT+1,SOUTH)

       vort(LORT+1,SOUTHWEST) = circ_SW_LORT / triarea(LORT+1,SOUTHWEST)
       vort(UPLT+1,SOUTHWEST) = vort(LORT+1,SOUTHWEST)

       pv_mass(1) = &
            rho(WEST)  * area(6,WEST) + &
            rho(0)     * sum(area(4:5,0)) + &
            rho(SOUTH) * area(3,SOUTH)

       pv_mass(2) = &
            rho(WEST)  * area(1,WEST) + &
            rho(0)     * area(3,0) + &
            rho(NORTH) * area(5,NORTH)

       pv_mass(3) = &
            rho(SOUTH) * area(2,SOUTH) + &
            rho(EAST)  * area(4,EAST) + &
            rho(0)     * area(6,0)

       pv_W_LORT  = (coriolis(LORT+1,WEST) + circ_W_LORT ) / pv_mass(1)
       pv_SW_LORT = (coriolis(LORT+1,SOUTHWEST) + circ_SW_LORT) / pv_mass(2)
       pv_S_UPLT  = (coriolis(UPLT+1,SOUTH) + circ_S_UPLT ) / pv_mass(3)

       pv_SW_UPLT = pv_SW_LORT

       qe(RT+1,WEST) = interp (pv_W_LORT, pv_SW_UPLT)
       qe(UP+1,SOUTH) = interp (pv_S_UPLT, pv_SW_LORT)
    end if

    ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*SOUTH+UP+1) = 0 in this case
    if (c == IPLUSJMINUS) then

       u_prim_RT_SW = velo(RT+1,SOUTHWEST) * length(RT+1,SOUTHWEST)
       u_prim_DG_SW = velo(DG+1,SOUTHWEST) * length(DG+1,SOUTHWEST)
       u_prim_RT    = velo(RT+1,0) * length(RT+1,0)

       circ_SW_LORT = - u_prim_RT + u_prim_RT_SW + u_prim_DG_SW
       circ_LORT    = vort(LORT+1,0) * triarea(LORT+1,0)
       circ_SW_UPLT = vort(UPLT+1,SOUTHWEST) * triarea(UPLT+1,SOUTHWEST)

       vort(LORT+1,SOUTHWEST) = circ_SW_LORT / triarea(LORT+1,SOUTHWEST)
       vort(UPLT+1,SOUTH) = vort(LORT+1,SOUTHWEST)

       pv_mass(1) = &
            rho(SOUTHWEST) * area(1,SOUTHWEST) + &
            rho(SOUTH)     * area(3,SOUTH) + &
            rho(0)         * sum(area(5:6,0))

       pv_mass(2) = &
            rho(0)         * area(1,0) + &
            rho(EAST)      * area(3,EAST) + &
            rho(NORTHEAST) * area(5,NORTHEAST)

       pv_mass(3) = &
            rho(SOUTHWEST) * area(2,SOUTHWEST) + &
            rho(0)         * area(4,0) + &
            rho(WEST)      * area(6,WEST)

       pv_SW_LORT = (coriolis(LORT+1,SOUTHWEST) + circ_SW_LORT) / pv_mass(1)
       pv_LORT    = (coriolis(LORT+1,0) + circ_LORT   ) / pv_mass(2)
       pv_SW_UPLT = (coriolis(UPLT+1,SOUTHWEST) + circ_SW_UPLT) / pv_mass(3)

       pv_S_UPLT = pv_SW_LORT

       qe(RT+1,0) = interp (pv_S_UPLT,  pv_LORT)
       qe(DG+1,SOUTHWEST) = interp (pv_SW_LORT, pv_SW_UPLT)
    end if

    ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*WEST+RT+1) = 0 in this case
    if (c == IMINUSJPLUS) then

       u_prim_UP    = velo(UP+1,0) * length(UP+1,0)
       u_prim_DG_SW = velo(DG+1,SOUTHWEST) * length(DG+1,SOUTHWEST)
       u_prim_UP_SW = velo(UP+1,SOUTHWEST) * length(UP+1,SOUTHWEST)

       circ_SW_UPLT = u_prim_UP - u_prim_DG_SW - u_prim_UP_SW
       circ_UPLT    = vort(UPLT+1,0) * triarea(UPLT+1,0)
       circ_SW_LORT = vort(LORT+1,SOUTHWEST) * triarea(LORT+1,SOUTHWEST)

       vort(UPLT+1,SOUTHWEST) = circ_SW_UPLT / triarea(UPLT+1,SOUTHWEST)
       vort(LORT+1,WEST) = vort(UPLT+1,SOUTHWEST)

       pv_mass(1) = &
            rho(SOUTHWEST) * area(2,SOUTHWEST) + &
            rho(0)         * sum(area(3:4,0)) + &
            rho(WEST)      * area(6,WEST)

       pv_mass(2) = &
            rho(0)         * area(2,0) + &
            rho(NORTHEAST) * area(4,NORTHEAST) + &
            rho(NORTH)     * area(6,NORTH)

       pv_mass(3) = &
            rho(SOUTHWEST) * area(1,SOUTHWEST) + &
            rho(SOUTH)     * area(3,SOUTH) + &
            rho(0)         * area(5,0)

       pv_SW_UPLT = (coriolis(UPLT+1,SOUTHWEST) + circ_SW_UPLT) / pv_mass(1)
       pv_UPLT    = (coriolis(UPLT+1,0) + circ_UPLT   ) / pv_mass(2)
       pv_SW_LORT = (coriolis(LORT+1,SOUTHWEST) + circ_SW_LORT) / pv_mass(3)

       pv_W_LORT = pv_SW_UPLT

       qe(UP+1,0) = interp (pv_W_LORT,  pv_UPLT)
       qe(DG+1,SOUTHWEST) = interp (pv_SW_LORT, pv_SW_UPLT)
    end if

    ! Parts 1, 2 of hexagon IJPLUS (upper right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*0+DG+1) = 0 in this case
    if (c == IJPLUS) then

       u_prim_RT   = velo(RT+1,0) * length(RT+1,0)
       u_prim_RT_N = velo(RT+1,NORTH) * length(RT+1,NORTH)
       u_prim_UP   = velo(UP+1,0) * length(UP+1,0)

       circ_LORT   = u_prim_RT - u_prim_RT_N - u_prim_UP
       circ_W_LORT = vort(LORT+1,WEST) * triarea(LORT+1,WEST)
       circ_S_UPLT = vort(UPLT+1,SOUTH) * triarea(UPLT+1,SOUTH)

       vort(LORT+1,0) = circ_LORT / triarea(LORT+1,0)
       vort(UPLT+1,0) = vort(LORT+1,0)

       pv_mass(1) = &
            rho(EAST)  * area(3,EAST) + &
            rho(0)     * sum(area(1:2,0)) + &
            rho(NORTH) * area(6,NORTH)

       pv_mass(2) = &
            rho(WEST)  * area(1,WEST) + &
            rho(0)     * area(3,0) + &
            rho(NORTH) * area(5,NORTH)

       pv_mass(3) = &
            rho(SOUTH) * area(2,SOUTH) + &
            rho(EAST)  * area(4,EAST) + &
            rho(0)     * area(6,0)

       pv_LORT   = (coriolis(LORT+1,0) + circ_LORT  ) / pv_mass(1)
       pv_W_LORT = (coriolis(LORT+1,WEST) + circ_W_LORT) / pv_mass(2)
       pv_S_UPLT = (coriolis(UPLT+1,SOUTH) + circ_S_UPLT) / pv_mass(3)

       pv_UPLT = pv_LORT

       qe(RT+1,0) = interp (pv_LORT, pv_S_UPLT)
       qe(UP+1,0) = interp (pv_UPLT, pv_W_LORT)
    end if
  end subroutine
end module
