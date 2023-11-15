MODULE profiling

  IMPLICIT NONE
  SAVE
  PRIVATE

  INTEGER :: id_phyparam

  INTERFACE  ! Explicit interfaces for plugins

     ! Plugin that registers a string identifier and returns an integer id to be used in enter/exit
     SUBROUTINE plugin_profile_register(thename, id)
       CHARACTER(*), INTENT(IN) :: thename
       INTEGER, INTENT(OUT) :: id
     END SUBROUTINE plugin_profile_register

     ! Plugin to be called when entering/exiting a region identified by its id
     SUBROUTINE plugin_profile_enter_exit(id)
       INTEGER, INTENT(IN) :: id
     END SUBROUTINE plugin_profile_enter_exit

  END INTERFACE

  ! This module provides default implementations of plugins that do nothing
  ! but the top-level driver is welcome to override them.

  ! Note F2003/F2008: pgfortran (F2003) accepts to initialize pointers only to NULL()
  ! => plugins are initialzed to NULL() and set to default values in flush_log and log_gridpoint

  PUBLIC :: profile_register, profile_enter, profile_exit, &
       id_phyparam

  PROCEDURE(plugin_profile_register),   POINTER, PUBLIC :: profile_register_plugin => NULL()
  PROCEDURE(plugin_profile_enter_exit), POINTER, PUBLIC :: profile_enter_plugin    => NULL(), &
       &                                                   profile_exit_plugin     => NULL()

CONTAINS

  SUBROUTINE profile_register(thename, id)
    CHARACTER(*), INTENT(IN) :: thename
    INTEGER, INTENT(OUT) :: id
    IF(ASSOCIATED(profile_register_plugin)) CALL profile_register_plugin(thename, id)
    !    PRINT *, 'profile_register'
  END SUBROUTINE profile_register

  SUBROUTINE profile_enter(id)
    INTEGER, INTENT(IN) :: id
    IF(ASSOCIATED(profile_enter_plugin)) CALL profile_enter_plugin(id)
    !    PRINT *, 'profile_enter'
  END SUBROUTINE profile_enter

  SUBROUTINE profile_exit(id)
    INTEGER, INTENT(IN) :: id
    IF(ASSOCIATED(profile_exit_plugin)) CALL profile_exit_plugin(id)
    !    PRINT *, 'profile_exit'
  END SUBROUTINE profile_exit

END MODULE profiling
