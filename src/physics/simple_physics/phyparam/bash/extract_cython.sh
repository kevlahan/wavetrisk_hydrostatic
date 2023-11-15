#!/bin/bash

function extract_header()
{
    grep header tmp/cython | sed -E 's/.*cython.*header//g' | awk '{print "  " $0}'
}

function extract_wrapper()
{
    grep wrapper tmp/cython | sed -n -e 's/^.* wrapper //p'
}


#	grep -v extern python/phyparam.hpxd > python/phyparam.h
#	cat python/phyparam.hpxd | tr ';' ' ' > python/phyparam.pxd
#	cd python ; rm -rf build *.c ; python setup.py build_ext --inplace

#cat tmp/all_code | extract_header

grep '!$cython' tmp/all_code > tmp/cython

extract_header > python/phyparam.h

cat <<EOF > python/phyparam.pxd
cdef extern from "phyparam.h" :
$(cat python/phyparam.h | tr -d ';')
EOF

cat <<EOF > python/phyparam.pyx
# cython: language_level=3

cimport phyparam as phy

#------------------------ wrappers  extracted from Fortran --------------------------#

$(extract_wrapper)

#-------------------------- plugins from python/plugins.pyx -------------------------#

$(cat python/plugins.pyx)
EOF
