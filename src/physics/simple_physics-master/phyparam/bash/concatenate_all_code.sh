#!/bin/bash

# Usage : $0 name1 name2 ...

for src in $* ; do
    awk "{ print \"$src \" NR \" \" \$0; }" physics/${src}.F90
done
