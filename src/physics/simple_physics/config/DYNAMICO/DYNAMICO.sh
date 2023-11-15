#!/bin/bash

cat <<EOF 
Usage : 
$0 install [ --branch devel2master ]
   download XIOS and DYNAMICO ; optionnally select specific git branch for DYNAMICO (defaults to trunk2master)
$0 avail
   find out possible values for ARCH parameter
$0 build_nophys ARCH [options]
   build XIOS and DYNAMICO with --arch=ARCH ; [options] (e.g. -debug) are passed to make_icosa.
$0 build ARCH [options]
   build XIOS, DYNAMICO and simple physics with --arch=ARCH
EOF

function cmd_install()
{
    echo "Downloading XIOS and DYNAMICO from forge.ipsl.fr via http"
    rm -rf modeles/XIOS modeles/DYNAMICO
    svn co "http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk" -r 2109 modeles/XIOS
    cp -pr patch/XIOS/arch/* modeles/XIOS/arch/
    git clone https://gitlab.in2p3.fr/ipsl/projets/dynamico/dynamico.git $* modeles/DYNAMICO
}

function cmd_avail()
{
    cp -pru patch/XIOS/arch/* modeles/XIOS/arch/
    cp -pru patch/DYNAMICO/arch/* modeles/DYNAMICO/arch/    

    echo "Available architectures values : "
    cd modeles/DYNAMICO/arch
    ls *.fcm
    echo "Each file arch-ARCH.fcm corresponds to a value of argument ARCH."
}

function cmd_build_nophys()
{
    cmd_avail
    
    ARCH=$1 ; shift
    echo "$0 build_nophys $ARCH $*" > rebuild_DYNAMICO.sh
    chmod a+x rebuild_DYNAMICO.sh

    cd -P $ROOT/modeles/XIOS
    echo "./make_xios --arch $ARCH --job 16" > rebuild
    chmod a+x rebuild
    ./rebuild | tee build.log || exit -1

    cd -P $ROOT/modeles/DYNAMICO
    git pull --ff-only
    ./make_icosa -arch $ARCH -parallel mpi -job 8 $* || exit -1
    cd $ROOT
}

function cmd_build_all()
{
    cmd_build_nophys $* || exit -1
    echo "$0 build_all $*" > rebuild_DYNAMICO.sh
    chmod a+x rebuild_DYNAMICO.sh

    cd -P $ROOT/modeles/DYNAMICO_phyparam
    ./make_dynamico_phyparam

    cd $ROOT/TEST_devel
    cp $ROOT/modeles/DYNAMICO/xml/* .
    cp xml/* .
}

function cmd_build()
{
    cmd_build_all $* -with_xios 
}

function cmd_srun()
{
    ACCOUNT=$1 ; shift
    set +x
    srun --partition=compil --ntasks=1 --cpus-per-task=16 --hint=nomultithread  --account=$ACCOUNT ./DYNAMICO.sh $*
}

function cmd_()
{
    echo
}

which svn 2>/dev/null || echo "We need svn to check out XIOS and DYNAMICO."
which svn 2>/dev/null || exit

cd -P $(dirname $0)
ROOT=$PWD

CMD=$1 ; shift
cmd_$CMD $*
