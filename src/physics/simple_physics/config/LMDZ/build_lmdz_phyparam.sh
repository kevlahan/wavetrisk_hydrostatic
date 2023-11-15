#!/bin/bash

function cmd_install_lmdz()
{
# Installation du modele en mode sequentiel
    echo "cmd_install_lmdz"
    wget -N http://www.lmd.jussieu.fr/~lmdz/pub/install_lmdz.sh
    export LANG=C # fixes issue with sed on MaxOSX
    sed -e 's/veget=1/veget=0/g' install_lmdz.sh > install_lmdz_patched.sh
    sed -i .bak -e 's/makelmdz_fcm/echo makelmdz_fcm/g'  install_lmdz_patched.sh
    rm -f *.bak
    chmod +x install_lmdz_patched.sh
    echo "Watch $PWD/install_lmdz_patched.log"
    ./install_lmdz_patched.sh -parallel none -v $version >install_lmdz_patched.log 2>&1
}

function cmd_patch_lmdz()
{
# Modification du code source pour prendre en compte la physique
# a 20 parametres
    echo "cmd_patch_lmdz"
    cd $LMDZ/libf
    rm -rf phyparam dynphy_lonlat/phyparam
    mkdir phyparam dynphy_lonlat/phyparam
    cd phyparam
    ln -s ../phydev/* .
    rm -f physiq_mod.F90 # we have our own one in dynphy_lonlat
    ln -sf $ROOT/phyparam/param/* .
    ln -sf $ROOT/phyparam/physics/* .
    cd ../dynphy_lonlat/phyparam
    ln -s ../phydev/* .
    ln -sf $ROOT/phyparam/dynphy_lonlat/* .
    cd $LMDZ
    echo "./makelmdz_fcm $* -rrtm false  -v false -arch local -j 8 -p param -d 32x32x39 gcm" > compile.sh
    chmod +x compile.sh
    ./compile.sh
}

function cmd_full()
{
    cd $LMDZ
    echo "./makelmdz_fcm -rrtm false  -v false -arch local -j 8 -p param -d 32x32x39 -full gcm" > compile.sh
    chmod +x compile.sh
    ./compile.sh
}

function cmd_()
{
    rm -rf LMDZ$version install_lmdz.*
    cmd_install_lmdz
    cmd_patch_lmdz
    echo "Now cd TEST_PARAM and execute ./gcm.e"
}

# On peut choisir la version de LMDZ a insitaller
version=20191106.trunk
LMDZ=$PWD/LMDZ$version/modipsl/modeles/LMDZ
ROOT=$(cd -P ../.. ; pwd)

CMD=$1
shift
cmd_$CMD $*

