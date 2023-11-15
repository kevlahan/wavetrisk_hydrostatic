#!/bin/bash

# some quotes in French-language comments break strict conformance to F2003
# this script removes them but could affect real strings ending with a single l or d after a space
# => check with svn diff after applying

for x in $* ; do 
    sed -i -e "s/ l'/ l /g" $x
    sed -i -e "s/ d'/ d /g" $x
done
