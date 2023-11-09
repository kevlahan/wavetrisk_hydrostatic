#!/bin/bash

# NOTE : sed -i seems to work differently on GNU and OSX systems
# a temporary file is used instead

which emacs > /dev/null 2>&1 || exit 0

TMP=$(mktemp)

for x in $* ; do
    # indent source file using emacs in batch mode
    # following standard settings of emacs f90-mode
    emacs -batch $x --eval '(f90-mode)' --eval '(indent-region (point-min) (point-max) nil)' --eval '(delete-trailing-whitespace)' -f save-buffer 2>/dev/null
    # remove blank lines at end of file
    sed  -e :a -e '/^\n*$/{$d;N;};/\n$/ba' $x  > $TMP
    diff -q $TMP $x || cp -f $TMP $x
done

rm -f $TMP
