  USE logging
#define WRITELOG(junk, fmt) logging_lineno = MIN(logging_bufsize, logging_lineno+1) ; WRITE(logging_buf(logging_lineno), fmt)
#define LOG_WARN(tag) CALL flush_log(log_level_warn, tag)
#define LOG_INFO(tag) CALL flush_log(log_level_info, tag)
#define LOG_DBG(tag) CALL flush_log(log_level_dbg, tag)
