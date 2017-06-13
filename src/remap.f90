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
    type(Domain) dom
    integer p
    integer i, j, k, m
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,9) :: dims
    integer fid
    integer id, idW, idSW, idS
    integer d, outl
    real(8) integrated_mass(zlevels), pressure(zlevels)

    d = dom%id + 1

    id   = idx(i,     j,     offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)

    !integrate all quantities vertically downward
    integrated_mass(1)=0.5_8*(sol(S_MASS,zlevels)%data(dom%id+1)%elts(id+1)+mean(S_MASS,zlevels))
    pressure(1)=grav_accel*integrated_mass(1)
    do k = 2 , zlevels
       integrated_mass(k)=integrated_mass(k-1)+ 0.5_8* &
            (sol(S_MASS,zlevels-k+1)%data(dom%id+1)%elts(id+1)+mean(S_MASS,zlevels-k+1)+ &
            sol(S_MASS,zlevels-k+2)%data(dom%id+1)%elts(id+1)+mean(S_MASS,zlevels-k+2))
       pressure(k)=grav_accel*integrated_mass(k)
    end do
  end subroutine remap_lagrangian

  subroutine seven_point_interpolation(xv, yv, xd, yd)
    !seven point interpolation as in Yang (2001)
    !uses the seven values in xv and yv to calculate the value yd at xd
    real(8) xv(7), yv(7), xd, yd, interp_diff(7,7)
    integer mi, ni

    !construct interpolating polynomial by calculating Newton differences
    interp_diff(1,:)=yv

    do mi = 2, 7
       do ni = 1, mi
          interp_diff(mi,ni)=(interp_diff(mi-1,ni+1)-interp_diff(mi-1,ni))/(xv(ni+1)-xv(ni))
       end do
    end do

    !evaluate the polynomial using Horner's algorithm
    yd=interp_diff(7,1)

    do mi = 6, 1, -1
       yd=interp_diff(mi,1)+(xd-xv(mi))*yd
    end do
  end subroutine seven_point_interpolation
end module remap_mod
