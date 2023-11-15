#---------------------------------- logging plugin -------------------------------#

cdef to_str(ptrlen, char[] ptr, start=0):
     return str(bytes([ ptr[i+start*ptrlen] for  i in range(ptrlen)]), encoding='utf8').rstrip()

"""
  SUBROUTINE default_flush_plugin(lev, taglen, tag, buflen, bufsize, buf)
    INTEGER(c_int), INTENT(IN), VALUE :: lev, taglen, buflen, bufsize
    CHARACTER(KIND=c_char), INTENT(IN) :: tag(taglen), buf(buflen, bufsize)
"""

cdef flush_plugin_wrap(int lev, int taglen, char[] tag, int buflen, int bufsize, char[] buf):
    flush_plugin(to_str(taglen, tag), [to_str(buflen, buf, i) for i in range(bufsize) ] )

def default_flush_plugin(tag, buf):
    if tag != 'missing_plugin' :
        n = len(buf)
        for i,line in enumerate(buf):
            print( '[%s %2.2d/%2.2d] %s'%(tag, i+1, n, line) )

def set_flush(plugin):
    flush_plugin = plugin

#----------------------------- plug plugins in ----------------------------------#

flush_plugin = default_flush_plugin

phy.phyparam_set_plugins_logging(&flush_plugin_wrap)
