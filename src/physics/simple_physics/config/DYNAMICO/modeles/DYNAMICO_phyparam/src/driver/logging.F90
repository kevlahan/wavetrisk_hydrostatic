MODULE logging

  ! see also use_logging.h
  ! macro LOGBUF accumulates log output into logging_buffer
  IMPLICIT NONE
  SAVE
  PRIVATE

  INTERFACE  ! Explicit interfaces for plugins
     ! Plugin that typically prints all lines in the loggin buffer 'buf' and prepends tags (log level, timestamp, ...)
     SUBROUTINE plugin(lev, taglen, tag, buflen, bufsize, buf) BIND(C)
       USE, INTRINSIC :: iso_c_binding, ONLY : c_char, c_null_char, c_int
       INTEGER(c_int), INTENT(IN), VALUE :: lev, taglen, buflen, bufsize
       CHARACTER(KIND=c_char), INTENT(IN) :: tag(taglen), buf(buflen, bufsize)
     END SUBROUTINE plugin

     ! Plugin that writes into string 'line' information about the gridpoint of index 'index'
     SUBROUTINE plugin_log_gridpoint(index, line_len, line) BIND(C)
       USE, INTRINSIC :: iso_c_binding, ONLY : c_char, c_null_char, c_int
       INTEGER(c_int), INTENT(IN), VALUE :: index, line_len ! index of gridpoint, LEN(line)
       CHARACTER(KIND=c_char), INTENT(OUT) :: line(line_len)
     END SUBROUTINE plugin_log_gridpoint

  END INTERFACE


  ! This module provides a default implementations of plugins but the top-level driver is welcome to override them.
  ! Note F2003/F2008: pgfortran (F2003) accepts to initialize pointers only to NULL()
  ! => plugins are initialzed to NULL() and set to default values in flush_log and log_gridpoint
  PUBLIC :: flush_plugin, log_gridpoint_plugin
  PROCEDURE(plugin), POINTER :: flush_plugin => NULL()
  PROCEDURE(plugin_log_gridpoint), POINTER :: log_gridpoint_plugin => NULL()

  INTEGER, PARAMETER :: linesize=10000, logging_bufsize=100
  CHARACTER(linesize) :: logging_buf(logging_bufsize)

  INTEGER :: logging_lineno=0

  ! messages with a log level > max_log_level are not printed
  INTEGER, PARAMETER, PUBLIC :: log_level_fatal=1, log_level_error=2, log_level_warn=3, log_level_info=4, log_level_dbg=5
  INTEGER, PUBLIC :: max_log_level = log_level_info
  CHARACTER(3), DIMENSION(log_level_dbg), PUBLIC :: dbtag = (/ 'FAT', 'ERR', 'WRN', 'INF', 'DBG' /)

  PUBLIC :: logging_buf, logging_bufsize, logging_lineno, flush_log, log_gridpoint, &
       missing_plugin

CONTAINS

  SUBROUTINE set_plugins(flush_plugin_c) BIND(C, name='phyparam_set_plugins_logging')
    !$cython header void phyparam_set_plugins_logging(void *);
    USE, INTRINSIC :: ISO_C_BINDING
    TYPE(C_FUNPTR), INTENT(IN), VALUE :: flush_plugin_c
    CALL C_F_PROCPOINTER(flush_plugin_c, flush_plugin)
  END SUBROUTINE set_plugins

  SUBROUTINE missing_plugin(name, mod)
    CHARACTER(*), INTENT(IN) :: name, mod
    WRITE(logging_buf(1),*) 'WARNING : plugin ', name, ' not provided by the driver program'
    WRITE(logging_buf(2),*) '        see ', mod
    logging_lineno=2
    CALL flush_log(log_level_warn, 'missing_plugin')
  END SUBROUTINE missing_plugin

  !-------------------------------------------------------------------------------------------------

  SUBROUTINE flush_log(lev,tag)
    INTEGER, INTENT(IN) :: lev
    CHARACTER(*), INTENT(IN) :: tag
    IF(.NOT.ASSOCIATED(flush_plugin)) THEN
       flush_plugin => default_flush_plugin
       WRITE(*,*) 'WARNING : plugin flush_plugin not provided by the driver program'
       WRITE(*,*) '        see logging.F90'
    END IF
    IF(logging_lineno>0 .AND. lev<=max_log_level) &
         CALL flush_plugin(lev, LEN(tag), TRIM(tag), linesize, logging_lineno, logging_buf)
    logging_lineno=0
  END SUBROUTINE flush_log

  SUBROUTINE default_flush_plugin(lev, taglen, tag, buflen, bufsize, buf) BIND(C)
    USE, INTRINSIC :: iso_c_binding, ONLY : c_char, c_null_char, c_int
    INTEGER(c_int), INTENT(IN), VALUE :: lev, taglen, buflen, bufsize
    CHARACTER(KIND=c_char), INTENT(IN) :: tag(taglen), buf(buflen, bufsize)
    CHARACTER(100) :: prefix
    CHARACTER(buflen+1) :: line
    INTEGER :: i
    WRITE(prefix,*) '[', dbtag(lev), ' ', tag, ']'
    DO i=1, bufsize
       WRITE(line,*) buf(:,i)
       WRITE(*,*) TRIM(prefix), TRIM(line)
    END DO
  END SUBROUTINE default_flush_plugin

  !-------------------------------------------------------------------------------------------------

  SUBROUTINE log_gridpoint(index)
    INTEGER, INTENT(IN) :: index
    logging_lineno = logging_lineno+1
    IF(.NOT.ASSOCIATED(log_gridpoint_plugin)) THEN
       log_gridpoint_plugin => default_log_gridpoint
       WRITE(*,*) 'WARNING : plugin log_gridpoint_plugin not provided by the driver program'
       WRITE(*,*) '        see logging.F90'
    END IF
    CALL log_gridpoint_plugin(index, linesize, logging_buf(logging_lineno))
  END SUBROUTINE log_gridpoint

  SUBROUTINE default_log_gridpoint(index, line_len, line) BIND(C)
    USE, INTRINSIC :: iso_c_binding, ONLY : c_char, c_null_char, c_int
    INTEGER(c_int), INTENT(IN), VALUE :: index, line_len ! index of gridpoint, LEN(line)
    CHARACTER(KIND=c_char), INTENT(OUT) :: line(line_len)
    line=' '
  END SUBROUTINE default_log_gridpoint

END MODULE logging
