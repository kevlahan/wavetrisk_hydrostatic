module template_mod
  use main_mod

  implicit none

contains
  subroutine apply_initial_conditions()
  end subroutine

  subroutine set_thresholds()
      toll_height = 
      toll_velo   = 
  end subroutine
end module

program template
  use main_mod
  use template_mod
  implicit none

  time_end = 

  call init_main_mod()

  call initialize(apply_initial_conditions, 1, set_thresholds, default_dump, default_load)
  do while (time .lt. time_end)
      call time_step(time_end, aligned)
  end do
  call finalize()
end program
