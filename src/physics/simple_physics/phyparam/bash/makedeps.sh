#!/bin/bash

# input file $1 is output of bash/concatenate_all_code.sh, see phyparam/Makefile

function print_deps()
{
    tr ',' ' ' | awk '{print "obj/" $1 ".o : " "obj/" $4 ".o" }'
}

grep '\bUSE\b' $1 | grep -i -v 'INTRINSIC' | grep -i -v '\bEND\b' | print_deps
