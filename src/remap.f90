module remap_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use comm_mpi_mod
  implicit none

contains

  subroutine remap_column(dom, p, i, j, k, offs, dims, fid)
    !remap the Lagrangian coordinate to prevent layer squeezing
    !inputs k and fid are unused
    !tested only for zlevels>7
    type(Domain) dom
    integer p
    integer i, j, k, m, kb
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,9) :: dims
    integer fid
    integer id
    integer e
    integer d
    real(8) integrated_mass(zlevels), integrated_temp(zlevels), pressure(zlevels)
    real(8) integrated_velo(zlevels,EDGE)
    real(8) new_mass(zlevels), new_temp(zlevels), new_velo(zlevels,EDGE) !all are integrated quantities
    real(8) top_press, bottom_press, centre_press
    integer stencil(7)

    d = dom%id + 1
    id   = idx(i,     j,     offs, dims)

    !integrate all quantities vertically downward
    integrated_mass(1)=sol(S_MASS,zlevels)%data(dom%id+1)%elts(id+1)
    integrated_temp(1)=sol(S_TEMP,zlevels)%data(dom%id+1)%elts(id+1)
    do e = 1, EDGE
       integrated_velo(1,e)=sol(S_VELO,zlevels)%data(dom%id+1)%elts(EDGE*id+e)
    end do
    do kb = 2 , zlevels
       integrated_mass(kb)=integrated_mass(kb-1)+sol(S_MASS,zlevels-kb+1)%data(dom%id+1)%elts(id+1)
       integrated_temp(kb)=integrated_temp(kb-1)+sol(S_TEMP,zlevels-kb+1)%data(dom%id+1)%elts(id+1)
       do e = 1, EDGE
          integrated_velo(kb,e)=integrated_velo(kb-1,e)+sol(S_VELO,zlevels-kb+1)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do

    !calculate current pressure distribution as it will be the independent coordinate
    pressure(1)=0.5_8*grav_accel*sol(S_MASS,zlevels)%data(dom%id+1)%elts(id+1)
    do kb = 2 , zlevels
       pressure(kb)=pressure(kb-1)+ 0.5_8*grav_accel* &
            (sol(S_MASS,zlevels-kb+1)%data(dom%id+1)%elts(id+1)+ &
            sol(S_MASS,zlevels-kb+2)%data(dom%id+1)%elts(id+1))
    end do
    PRINT *, 'pressure is', pressure(1:zlevels)

    !interpolate using the moving stencil (note that in case of extreme shifting of the layers, we may extrapolate)
    do kb = 1, zlevels
       !top_press=a_vert(kb+1)*ref_press+b_vert(kb+1)*dom%surf_press%elts(id+1) !should be ref_press_t JEMF
       !bottom_press=a_vert(kb)*ref_press+b_vert(kb)*dom%surf_press%elts(id+1)
       centre_press=0.5_8 +1.0_8*(kb-1) !(top_press+bottom_press) !JEMF temporarily disabled for tenlayergauss
       if (kb .le. 4) then
          stencil=(/ (m, m = 1,7) /)
       else if (kb .ge. (zlevels-3)) then
          stencil=(/ (m, m = zlevels-6,zlevels) /)
       else
          stencil=(/ (m, m = kb-3, kb+3) /)
       end if
       new_mass(kb)=seven_point_interp(pressure(stencil), integrated_mass(stencil), centre_press)
       new_temp(kb)=seven_point_interp(pressure(stencil), integrated_temp(stencil), centre_press)
       do e = 1, EDGE
          new_velo(kb,e)=seven_point_interp(pressure(stencil), integrated_velo(stencil,e), centre_press)
       end do
    end do

    !assign values to mass, temp and velo field
    sol(S_MASS,zlevels)%data(dom%id+1)%elts(id+1)=new_mass(1)
    sol(S_TEMP,zlevels)%data(dom%id+1)%elts(id+1)=new_temp(1)
    do e = 1, EDGE
       sol(S_VELO,zlevels)%data(dom%id+1)%elts(EDGE*id+e)=new_velo(1,e)
    end do
    do kb = 2 , zlevels
       sol(S_MASS,zlevels-kb+1)%data(dom%id+1)%elts(id+1)=new_mass(kb)-new_mass(kb-1)
       sol(S_TEMP,zlevels-kb+1)%data(dom%id+1)%elts(id+1)=new_temp(kb)-new_temp(kb-1)
       do e = 1, EDGE
          sol(S_VELO,zlevels-kb+1)%data(dom%id+1)%elts(EDGE*id+e)=new_velo(kb,e)-new_velo(kb-1,e)
       end do
    end do
  end subroutine remap_column

  function seven_point_interp(xv, yv, xd)
    !seven point interpolation as in Yang (2001)
    !uses the seven values in xv and yv to calculate the value yd at xd
    real(8) xv(0:6), yv(0:6), xd, seven_point_interp
    real(8) interp_diff(0:6,0:6) !second index describes order of finite difference
    integer mi, ni

    !construct interpolating polynomial by calculating Newton differences
    interp_diff(0,0:6)=yv(0:6) !zeroth order finite differences for x_0 thru x_6

    do mi = 1, 6
       do ni = 0, 6-mi
          interp_diff(mi,ni)=(interp_diff(mi-1,ni+1)-interp_diff(mi-1,ni))/(xv(ni+mi)-xv(ni))
       end do
    end do

    !evaluate the polynomial using Horner's algorithm
    seven_point_interp=interp_diff(6,0)

    do mi = 5, 0, -1
       seven_point_interp=interp_diff(mi,0)+(xd-xv(mi))*seven_point_interp
    end do
  end function seven_point_interp

  subroutine remap_coordinates
    integer l

    PRINT *, 'we are remapping the coordinates'

    do l = level_start, level_end
       call write_level_mpi(remap_column, 0, l, 0, .True.) !we are not really writing; both 0 spots are not used
    end do

    call barrier
  end subroutine remap_coordinates
end module remap_mod
