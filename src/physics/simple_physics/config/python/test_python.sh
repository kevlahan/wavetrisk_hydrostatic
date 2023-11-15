#!/bin/env bash

# modules for camelot :
#     module load python/3.6-anaconda50gnu/7.2.0

cd ../../phyparam
make py
cd -

rm -f lib phyparam*.so
ln -s ../../phyparam/python/phyparam*.so .
ln -s ../../phyparam/lib .

LD_LIBRARY_PATH=lib python test_python.py
