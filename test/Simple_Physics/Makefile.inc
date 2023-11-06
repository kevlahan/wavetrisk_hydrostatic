SRC = $(PARAM).f90 shared.f90 coord_arithmetic.f90 sphere.f90  patch.f90 dyn_array.f90 \
base_$(PROC).f90 dgesv.f90 dgtsv.f90 spline.f90 domain.f90 init.f90 comm.f90 comm_$(PROC).f90 utils.f90 projection.f90 equation_of_state.f90 \
wavelet.f90 lnorms.f90 mask.f90 refine_patch.f90 smooth.f90 ops.f90 multi_level.f90 adapt.f90 lin_solve.f90  \
barotropic_2d.f90 time_integr.f90 vert_diffusion.f90 io.f90 remap.f90 main.f90 test_case_module.f90 \
init_physics_module.f90 physics_call_module.f90 phys_processing.f90 test.f90
