module remap_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use adapt_mod
  use smooth_mod
  use comm_mpi_mod
  implicit none

contains

  subroutine remap_lagrangian(dom, p, i, j, k, offs, dims, fid)
    !remap the Lagrangian coordinate to prevent layer squeezing
    !input k is unused, as is input fid
    !tested only for zlevels>7
    type(Domain) dom
    integer p
    integer i, j, k, m
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,9) :: dims
    integer fid
    integer id
    integer e
    integer d
    real(8) integrated_mass(zlevels), integrated_temp(zlevels), pressure(zlevels)
    real(8) integrated_velo(zlevels,EDGE)
    real(8) new_mass(zlevels), new_temp(zlevels), new_velo(zlevels,EDGE) !all are integrated quantities
    real(8) dep_stencil(7), indep_stencil(7), top_press, bottom_press, centre_press
    integer stencil(7)

    d = dom%id + 1

    id   = idx(i,     j,     offs, dims)

    !integrate all quantities vertically downward
    integrated_mass(1)=sol(S_MASS,zlevels)%data(dom%id+1)%elts(id+1)
    integrated_temp(1)=sol(S_TEMP,zlevels)%data(dom%id+1)%elts(id+1)
    do e = 1, EDGE
       integrated_velo(1,e)=sol(S_VELO,zlevels)%data(dom%id+1)%elts(EDGE*id+e)
    end do
    do k = 2 , zlevels
       integrated_mass(k)=integrated_mass(k-1)+sol(S_MASS,zlevels-k+1)%data(dom%id+1)%elts(id+1)
       integrated_temp(k)=integrated_temp(k-1)+sol(S_TEMP,zlevels-k+1)%data(dom%id+1)%elts(id+1)
       do e = 1, EDGE
          integrated_velo(k,e)=integrated_velo(k-1,e)+sol(S_VELO,zlevels-k+1)%data(dom%id+1)%elts(EDGE*id+e)
       end do
    end do

    !calculate pressure as it will be the independent coordinate
    pressure(1)=0.5_8*grav_accel*sol(S_MASS,zlevels)%data(dom%id+1)%elts(id+1)
    do k = 2 , zlevels
       pressure(k)=integrated_mass(k-1)+ 0.5_8*grav_accel* &
            (sol(S_MASS,zlevels-k+1)%data(dom%id+1)%elts(id+1)+ &
            sol(S_MASS,zlevels-k+2)%data(dom%id+1)%elts(id+1))
    end do

    !interpolate using the moving stencil
    do k = 1, zlevels
       top_press=a_vert(k+1)*ref_press+b_vert(k+1)*dom%surf_press%elts(id+1) !should be ref_press_t JEMF
       bottom_press=a_vert(k)*ref_press+b_vert(k)*dom%surf_press%elts(id+1)
       centre_press=0.5_8*(top_press+bottom_press)
       if (k .le. 4) then
          stencil=(/ (m, m = 1,7) /)
       else if (k .ge. (zlevels-3)) then
          stencil=(/ (m, m = zlevels-6,zlevels) /)
       else
          stencil=(/ (m, m = k-3, k+3) /)
       end if
       new_mass(k)=seven_point_interp(pressure(stencil), integrated_mass(stencil), centre_press)
       new_temp(k)=seven_point_interp(pressure(stencil), integrated_temp(stencil), centre_press)
       do e = 1, EDGE
          new_velo(k,e)=seven_point_interp(pressure(stencil), integrated_velo(stencil,e), centre_press)
       end do
    end do
  end subroutine remap_lagrangian

  function seven_point_interp(xv, yv, xd)
    !seven point interpolation as in Yang (2001)
    !uses the seven values in xv and yv to calculate the value yd at xd
    real(8) xv(7), yv(7), xd, seven_point_interp, interp_diff(7,7)
    integer mi, ni

    !construct interpolating polynomial by calculating Newton differences
    interp_diff(1,:)=yv

    do mi = 2, 7
       do ni = 1, mi
          interp_diff(mi,ni)=(interp_diff(mi-1,ni+1)-interp_diff(mi-1,ni))/(xv(ni+1)-xv(ni))
       end do
    end do

    !evaluate the polynomial using Horner's algorithm
    seven_point_interp=interp_diff(7,1)

    do mi = 6, 1, -1
       seven_point_interp=interp_diff(mi,1)+(xd-xv(mi))*seven_point_interp
    end do
  end function seven_point_interp
end module remap_mod
