#!/bin/bash


# Compiles the omni-compiler suite, especially the Fortran front-end F_Front 
# which is used to parse the physics code and (later) check that coding standards are followed.

# F_Front will be installed in the directory from which this script has been invoked. 
# Before executing this script, cd to the directory where you want to install.
# e.g. > cd $HOME/local . Then F_Front will be in $HOME/local/bin . 

# There are chances this script fails on OSX if the standard gcc is used. You can try with a more recent gcc, which is easy to 
# install through Homebrew. The XML library libxml2 will be required too :
#   * install gcc and libxml2 : > brew install gcc libxml2
#   * check : > which gfortran
#   * cd to simple_physics/bash/omnicompiler-1.3.2
#   * run the configure script with :
#          > CC=/usr/local/bin/gcc-9 ./configure --prefix=$HOME/local --with-libxml2=/usr/local/Cellar/libxml2/2.9.10/
#    (your version numbers may be different
#   * make all install
#   * even if the full build fails, you may find a working F_Front in xcodeml-tools/F-FrontEnd/src :
#          > xcodeml-tools/F-FrontEnd/src/F_Front --help
#   * copy F_Front where desired
 
PREFIX=$(pwd)
ROOT=$(dirname $0)
ROOT=$(cd -P $ROOT ; pwd)

OMNI=omnicompiler-1.3.2

#echo "prefix : $PREFIX"
#echo "ROOT   : $ROOT"
cd $ROOT
rm -rf ${OMNI}
tar xjf ${OMNI}.tar.bz2
cd ${OMNI}
./configure --prefix=$PREFIX
make -j 8
make install
ls -l $PREFIX/bin
