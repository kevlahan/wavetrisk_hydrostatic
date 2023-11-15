#!/bin/bash

echo <<EOF
usage : replace_print file1 file2 ...
This script replaces occurences of PRINT and WRITE with the macro WRITELOG
User is asked for confirmation before overwriting original file.
This script is basic and you should check results as corner cases are not handled.
EOF

function replace()
{
    FILE=$1.new
    cp -f $1 $FILE
# WRITE => WRITELOG
    sed -i -e "s/\bwrite(\b/WRITELOG(/g" $FILE
    sed -i -e "s/\bWRITE(\b/WRITELOG(/g" $FILE
# print * => WRITELOG(*,*)
    sed -i -e "s/\bprint\*\,/WRITELOG(*,*) /g" $FILE
    sed -i -e "s/\bPRINT\*\,/WRITELOG(*,*) /g" $FILE
    diff $1 $FILE
    cp -i -u $FILE $1
    rm -f $FILE
}

for FILE in $* ; do
    replace $FILE
done
