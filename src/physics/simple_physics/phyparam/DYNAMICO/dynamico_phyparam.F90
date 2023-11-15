PROGRAM icosa_gcm
  USE icosa_init_mod, ONLY : icosa_init
  USE icosa_phyparam_mod, ONLY : init_physics, physics
  USE physics_mod, ONLY : init_physics_plugin, physics_plugin

  init_physics_plugin => init_physics
  physics_plugin => physics

  CALL icosa_init

END PROGRAM icosa_gcm

SUBROUTINE initialize_external_physics
END SUBROUTINE initialize_external_physics

SUBROUTINE external_physics
END SUBROUTINE external_physics
