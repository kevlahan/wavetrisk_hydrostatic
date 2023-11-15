#!/bin/bash

echo <<EOF
usage : replace_keyword file1 file2 ...
This script replaces lowercase Fortran keywords by uppercase
User is asked for confirmation before overwriting original file.
This script is basic and you should check results as corner cases are not handled.
EOF

function replace_key()
{
    LOW="$1"
    UP="$2"
    sed -i -e "s/\b${LOW}\b/${UP}/g" $3
}

function replace()
{
    ORIG=$1
    NEW=$1.new
    cp -f $ORIG $NEW
# enddo => ENDDO
    replace_key subroutine SUBROUTINE $NEW
    replace_key 'do i=' 'DO i=' $NEW
    replace_key 'do k=' 'DO k=' $NEW
    replace_key 'end do' 'END DO' $NEW
    replace_key 'parameter' 'PARAMETER' $NEW
    replace_key 'integer' 'INTEGER' $NEW
    replace_key 'real' 'REAL' $NEW
    replace_key 'then' 'THEN' $NEW
    replace_key 'else' 'ELSE' $NEW
    replace_key 'endif' 'END IF' $NEW
    replace_key 'end if' 'END IF' $NEW
    diff $NEW $ORIG
    cp -i -u $NEW $ORIG
    rm -f $NEW
}

for FILE in $* ; do
    replace $FILE
done
