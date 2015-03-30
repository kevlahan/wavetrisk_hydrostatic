module time_integr_mod
  use dyn_arrays
  use comm_mod
  use ops_mod
  use multi_level_mod
  use io_mod
  use wavelet_mod
  use adapt_mod
  use arch_mod
  implicit none
  type(Float_Field), dimension(S_HEIGHT:S_VELO) :: q1, q2, q3, q4, dq1
contains
  subroutine init_time_integr_mod()
      logical :: initialized = .False.
      integer v, d
      if (initialized) return ! initialize only once
      call init_comm_mod()
      call init_ops_mod()
      call init_multi_level_mod()
      dt = 30.
      initialized = .True.
  end subroutine init_time_integr_mod
  subroutine RK_sub_step4(sol1, sol2, sol3, sol4, trend1, trend2, alpha, dt, dest)
      real(8) :: alpha(4), dt(2)
      type(Float_Field), dimension(S_HEIGHT:S_VELO) :: sol1, sol2, sol3, sol4
      type(Float_Field), dimension(S_HEIGHT:S_VELO)  :: trend1, trend2
      type(Float_Field), dimension(S_HEIGHT:S_VELO), intent(inout) :: dest
      integer v, s, t, d, start
      
      do d = 1, n_domain(rank+1)
          do v = S_HEIGHT, S_VELO
              start = (1+2*(v-1))*grid(d)%patch%elts(2+1)%elts_start ! start of second level
              dest(v)%data(d)%elts(start+1:dest(v)%data(d)%length) = &
                   alpha(1)*sol1(v)%data(d)%elts(start+1:sol1(v)%data(d)%length) &
                 + alpha(2)*sol2(v)%data(d)%elts(start+1:sol2(v)%data(d)%length) &
                 + alpha(3)*sol3(v)%data(d)%elts(start+1:sol3(v)%data(d)%length) &
                 + alpha(4)*sol4(v)%data(d)%elts(start+1:sol4(v)%data(d)%length) &
                 + dt(1)*trend1(v)%data(d)%elts(start+1:trend1(v)%data(d)%length) &
                 + dt(2)*trend2(v)%data(d)%elts(start+1:trend2(v)%data(d)%length)
          end do
      end do
      dest(:)%bdry_uptodate = .False.
  end subroutine


  subroutine RK_sub_step1(sols, trends, alpha, dt, dest)
      real(8) :: alpha, dt
      type(Float_Field), dimension(S_HEIGHT:S_VELO) :: sols
      type(Float_Field), dimension(S_HEIGHT:S_VELO)  :: trends
      type(Float_Field), dimension(S_HEIGHT:S_VELO), intent(inout) :: dest
      integer v, s, t, d, start
      
     do d = 1, n_domain(rank+1)
      do v = S_HEIGHT, S_VELO
          start = (1+2*(v-1))*grid(d)%patch%elts(2+1)%elts_start ! start of second level
          dest(v)%data(d)%elts(start+1:dest(v)%data(d)%length) = &
                   alpha*sols(v)%data(d)%elts(start+1:sols(v)%data(d)%length) &
                 + dt*trends(v)%data(d)%elts(start+1:trends(v)%data(d)%length)
      end do
    end do
    dest(:)%bdry_uptodate = .False.
  end subroutine
  subroutine RK_sub_step2(sol1, sol2, trends, alpha, dt, dest)
      real(8) :: alpha(2), dt
      type(Float_Field), dimension(S_HEIGHT:S_VELO) :: sol1, sol2
      type(Float_Field), dimension(S_HEIGHT:S_VELO)  :: trends
      type(Float_Field), dimension(S_HEIGHT:S_VELO), intent(inout) :: dest
      integer v, s, t, d, start
      
      do d = 1, n_domain(rank+1)
        do v = S_HEIGHT, S_VELO
          start = (1+2*(v-1))*grid(d)%patch%elts(2+1)%elts_start ! start of second level
          dest(v)%data(d)%elts(start+1:dest(v)%data(d)%length) = &
                   alpha(1)*sol1(v)%data(d)%elts(start+1:sol1(v)%data(d)%length) &
                 + alpha(2)*sol2(v)%data(d)%elts(start+1:sol2(v)%data(d)%length) &
                 + dt*trends(v)%data(d)%elts(start+1:trends(v)%data(d)%length)
        end do
      end do
      dest(:)%bdry_uptodate = .False.
  end subroutine

  subroutine init_RK_mem()
      integer d, v
      do v = S_HEIGHT, S_VELO
          call init_Float_Field(q1(v), v)
          call init_Float_Field(q2(v), v)
          call init_Float_Field(q3(v), v)
          call init_Float_Field(q4(v), v)
          call init_Float_Field(dq1(v), v)
      end do
      do d = 1, n_domain(rank+1)
          do v = S_HEIGHT, S_VELO
              call init(q1(v)%data(d), sol(v)%data(d)%length); q1(v)%data(d)%elts = dble(2-v)
              call init(q2(v)%data(d), sol(v)%data(d)%length); q2(v)%data(d)%elts = dble(2-v)
              call init(q3(v)%data(d), sol(v)%data(d)%length); q3(v)%data(d)%elts = dble(2-v)
              call init(q4(v)%data(d), sol(v)%data(d)%length); q4(v)%data(d)%elts = dble(2-v)
              call init(dq1(v)%data(d), sol(v)%data(d)%length); dq1(v)%data(d)%elts = dble(2-v)
          end do
      end do
  end subroutine
  subroutine manage_RK_mem()
      integer d, v, n_new
      do d = 1, n_domain(rank+1)
          do v = S_HEIGHT, S_VELO
              n_new = sol(v)%data(d)%length - q1(v)%data(d)%length
              if (n_new .gt. 0) then
                  call extend(q1(v)%data(d), n_new, dble(2-v))
                  call extend(q2(v)%data(d), n_new, dble(2-v))
                  call extend(q3(v)%data(d), n_new, dble(2-v))
                  call extend(q4(v)%data(d), n_new, dble(2-v))
                  call extend(dq1(v)%data(d), n_new, dble(2-v))
              end if
          end do
      end do
  end subroutine

  subroutine RK45_opt()
      real(8), dimension(5,5) :: alpha
      real(8), dimension(5,5) :: beta
      integer d, v
      alpha = reshape((/1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.444370494067_8, 0.555629505933_8, 0.0_8, 0.0_8, 0.0_8,  &
!     alpha = reshape((/1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
!                       0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8,  &
          0.620101851385_8, 0.0_8, 0.379898148615_8, 0.0_8, 0.0_8, 0.178079954108_8, 0.0_8, 0.0_8, 0.821920045892_8, 0.0_8,  &
          0.00683325884039_8, 0.0_8, 0.51723167209_8, 0.127598311333_8, 0.348336757737_8/), (/5, 5/))
      beta = reshape((/0.391752227004_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.36841059263_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
!     beta = reshape((/1.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
!                      0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8,  &
                       0.251891774247_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
                       0.544974750212_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0846041633821_8,  &
          0.226007483194_8/), (/5, 5/))

!     do d = 1, size(grid)
!         call associate_Float_Array(q0(HEIGHT,d), grid(d)%height)
!         call associate_Float_Array(q0(VELO,d), grid(d)%velo)
!         q0(HEIGHT:VELO,d) = (/grid(d)%height, grid(d)%velo/)
!     end do
      call manage_RK_mem()

      call trend_ml(sol, trend)
!     call RK_sub_step(sol, trend, (/alpha(1,1)/), dt*(/beta(1,1)/), q1)
!Euler
!     call RK_sub_step1(sol, trend, alpha(1,1), dt, sol)
!     call WT_after_step(sol, level_start-1)
!     return
!end Euler
!Heun
!     call RK_sub_step1(sol, trend, alpha(1,1), dt, q1)
!     call WT_after_step(q1)
!     call trend_ml(q1, trend)
!     call RK_sub_step2(sol, q1, trend, (/0.5_8, 0.5_8/), dt*0.5_8, sol)
!     call WT_after_step(sol, level_start-1)
!     return
!end Heun
      call RK_sub_step1(sol, trend, alpha(1,1), dt*beta(1,1), q1)
      call WT_after_step(q1)

      call trend_ml(q1, trend)
      call RK_sub_step2(sol, q1, trend, alpha(1:2,2), dt*beta(2,2), q2)
      call WT_after_step(q2)

      call trend_ml(q2, trend)
      call RK_sub_step2(sol, q2, trend, (/alpha(1,3), alpha(3,3)/), dt*beta(3,3), q3)
      call WT_after_step(q3)

      call trend_ml(q3, trend)
      call RK_sub_step2(sol, q3, trend, (/alpha(1,4), alpha(4,4)/), dt*beta(4,4), q4)
      call WT_after_step(q4)

      call trend_ml(q4, dq1)
      call RK_sub_step4(sol, q2, q3, q4, trend, dq1, (/alpha(1,5), alpha(3:5,5)/), &
                           dt*beta(4:5,5), sol)
      call WT_after_step(sol, level_start-1)
  end subroutine
  subroutine WT_after_step(q, l_start0)
      type(Float_Field), target :: q(S_HEIGHT:S_VELO)
      integer, optional :: l_start0
      integer l, d, k, l_start
      !  everything needed in terms of forward and backward wavelet transform
      !         after one time step (e.g. RK sub-step)
      !         A) compute height wavelets and perform backwards transform
      !            to conserve mass
      !         B) interpolate values for next step

      if (present(l_start0)) then
          l_start = l_start0
          if (max_level .gt. min_level) call apply_interscale(restrict_u, level_start-1, 0, 0)
      else
          l_start = level_start
      end if
      call update_bdry(q(S_HEIGHT), NONE)
      call update_bdry(q(S_VELO), NONE)
      do l = l_start, level_end-1
          do d = 1, n_domain(rank+1)
              wc_u => wav_coeff(S_VELO)%data(d)%elts
              velo => q(S_VELO)%data(d)%elts
              wc_h => wav_coeff(S_HEIGHT)%data(d)%elts
              height => q(S_HEIGHT)%data(d)%elts
              call apply_interscale_d(cpt_velo_wc, grid(d), l, 0, 0)
              call apply_interscale_d(cpt_height_wc, grid(d), l, 0, 0)
              call apply_to_penta_d(cpt_vel_wc_penta, grid(d), l)
          end do
          wav_coeff(:)%bdry_uptodate = .False.
      end do

      do l = level_start+1, level_end
          do d = 1, n_domain(rank+1)
              do k = 1, grid(d)%lev(l)%length
                  call apply_onescale_to_patch(compress, grid(d), grid(d)%lev(l)%elts(k), 0, 1)
              end do
          end do
          wav_coeff(:)%bdry_uptodate = .False.
      end do

      call invers_wavelet_transform(q)
  end subroutine
end module time_integr_mod
