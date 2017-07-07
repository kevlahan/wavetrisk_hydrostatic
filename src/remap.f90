module remap_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use comm_mpi_mod
  implicit none

contains

  subroutine remap_column(dom, p, i, j, k, offs, dims, fid)
    !remap the Lagrangian coordinate with the goal of preventing layer squeezing
    !the inputs k and fid are unused
    !tested only for zlevels>7
    !we interpolate on the full (not the perturbation) quantities
    type(Domain) dom
    integer p
    integer i, j, k, m, kb
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,9) :: dims
    integer fid
    integer id
    integer e
    integer d
    real(8) pressure(zlevels+1), layer_pressure
    real(8) integrated_velo(zlevels+1,EDGE), integrated_mass(zlevels+1), integrated_temp(zlevels+1)
    real(8) new_mass(zlevels+1), new_temp(zlevels+1), new_velo(zlevels+1,EDGE) !all are integrated quantities
    integer stencil(7)

    d = dom%id + 1
    id   = idx(i,     j,     offs, dims)

    !integrate all quantities vertically downward from the top, all quantities located at interfaces
    integrated_mass(1) = 0.0_8
    integrated_temp(1) = 0.0_8
    do e = 1, EDGE
       integrated_velo(1,e) = 0.0_8
    end do
    do kb = 2, zlevels + 1
       integrated_mass(kb) = integrated_mass(kb-1) + &
            sol(S_MASS,zlevels-kb+2)%data(dom%id+1)%elts(id+1) + mean(S_MASS,zlevels-kb+2)
       integrated_temp(kb) = integrated_temp(kb-1) + &
            sol(S_TEMP,zlevels-kb+2)%data(dom%id+1)%elts(id+1) + mean(S_TEMP,zlevels-kb+2)
       do e = 1, EDGE
          integrated_velo(kb,e) = integrated_velo(kb-1,e) + &
               sol(S_VELO,zlevels-kb+2)%data(dom%id+1)%elts(EDGE*id+e) + mean(S_VELO,zlevels-kb+2)

          !if the mean velocity is non-zero, this routine will not work properly
          if (abs(mean(S_VELO,zlevels-kb+2)).gt.1.0e-14) then
            write(6,*) 'fatal error: mean velocity is not zero'
            stop
          end if
       end do
    end do

    !calculate current pressure distribution (at the interfaces) as it will be the independent coordinate
    pressure(1) = press_infty
    do kb = 2, zlevels + 1
       pressure(kb) = pressure(kb-1) + grav_accel * &
            (sol(S_MASS,zlevels-kb+2)%data(dom%id+1)%elts(id+1)+mean(S_MASS,zlevels-kb+2))
    end do

    !interpolate using the moving stencil (note that in case of extreme shifting of the layers, we may extrapolate)
    !again new quantities are computed top-down
    do kb = 1, zlevels + 1
       if (allocated(a_vert).and.allocated(b_vert)) then !a_vert and b_vert are allocated so they will be used
          layer_pressure=a_vert(kb)*ref_press+b_vert(kb)*pressure(zlevels+1) !should be ref_press_t JEMF
       else  !layers will be equi-distributed
          layer_pressure=press_infty+(kb-1.0_8)*(pressure(zlevels+1)-press_infty)/zlevels !JEMF check ref_press or infty
       end if

       if (kb .le. 4) then
          stencil = (/ (m, m = 1, 7) /)
       else if (kb .ge. (zlevels-3)) then
          stencil = (/ (m, m = zlevels-5, zlevels+1) /)
       else
          stencil = (/ (m, m = kb-3, kb+3) /)
       end if
       new_mass(kb) = seven_point_interp(pressure(stencil), integrated_mass(stencil), layer_pressure)
       new_temp(kb) = seven_point_interp(pressure(stencil), integrated_temp(stencil), layer_pressure)
       do e = 1, EDGE
          new_velo(kb,e) = seven_point_interp(pressure(stencil), integrated_velo(stencil,e), layer_pressure)
       end do
    end do

    !temporarily assign full quantities to mass, temp and velo field (instead of perturbation quantities)
    do kb = 1 , zlevels
       sol(S_MASS,kb)%data(dom%id+1)%elts(id+1) = new_mass(zlevels-kb+2) - new_mass(zlevels-kb+1)
       sol(S_TEMP,kb)%data(dom%id+1)%elts(id+1) = new_temp(zlevels-kb+2) - new_temp(zlevels-kb+1)
       do e = 1, EDGE
          sol(S_VELO,kb)%data(dom%id+1)%elts(EDGE*id+e) = new_velo(zlevels-kb+2,e) - new_velo(zlevels-kb+1,e)
       end do
    end do

    !now assign perturbation quantities to mass, temp and velo field
    do kb = 1 , zlevels
       sol(S_MASS,kb)%data(dom%id+1)%elts(id+1) = sol(S_MASS,kb)%data(dom%id+1)%elts(id+1) - mean(S_MASS,kb)
       sol(S_TEMP,kb)%data(dom%id+1)%elts(id+1) = sol(S_TEMP,kb)%data(dom%id+1)%elts(id+1) - mean(S_TEMP,kb)
       do e = 1, EDGE
          sol(S_VELO,kb)%data(dom%id+1)%elts(EDGE*id+e) = sol(S_VELO,kb)%data(dom%id+1)%elts(EDGE*id+e) - &
               mean(S_VELO,kb)
       end do
    end do
  end subroutine remap_column

  function seven_point_interp(xv, yv, xd)
    !seven point interpolation as in Yang (2001)
    !uses the seven values in xv and yv to calculate the value yd at xd
    !this routine has been verified to work properly
    real(8) xv(0:6), yv(0:6), xd, seven_point_interp
    real(8) interp_diff(0:6,0:6) !second index describes order of finite difference
    integer mi, ni

    !construct interpolating polynomial by calculating Newton differences
    interp_diff(0,0:6) = yv(0:6) !zeroth order finite differences for x_0 thru x_6

    do mi = 1, 6
       do ni = 0, 6-mi
          interp_diff(mi,ni) = (interp_diff(mi-1,ni+1)-interp_diff(mi-1,ni))/(xv(ni+mi)-xv(ni))
       end do
    end do

    !evaluate the polynomial using Horner's algorithm
    seven_point_interp = interp_diff(6,0)

    do mi = 5, 0, -1
       seven_point_interp = interp_diff(mi,0)+(xd-xv(mi))*seven_point_interp
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
