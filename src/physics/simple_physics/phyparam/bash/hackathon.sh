#!/bin/bash

function main()
{
    rm -rf driver/driver_$1.exe outputs_$1/
    make -j driver/driver.exe || exit -1
    driver/driver.exe > outputs/log || exit -1
    mv outputs outputs_$1
    mv driver/driver.exe driver/driver_$1.exe    
}

function cmd_cpu()
{
    make clean
    F90FLAGS="-i4 -r8 -gopt -traceback -fast" main cpu || exit -1
    make clean
}

function main_gpu()
{
    export LDFLAGS="-L $CUDA_HOME/lib64/ -lnvToolsExt"
    export F90FLAGS="-i4 -r8 -gopt -traceback -fast -acc -Minfo=accel -Minline $*"
    main gpu || exit -1
    mkdir -p outputs
    for x in pt pu pv; do
	file="output_$x"
	echo "Verifying outputs/$file"
	diff outputs_cpu/$file outputs_gpu/$file
    done
}

function cmd_verif()
{
    # see https://developer.nvidia.com/blog/detecting-divergence-using-pcast-to-compare-gpu-to-cpu-results/
    export PCAST_COMPARE=rel=7
    main_gpu -ta=tesla:cc70,autocompare || exit -1
}

function cmd_gpu()
{
    main_gpu -ta=tesla:cc70 || exit -1
}

function cmd_()
{
    true
}

module purge
module load nvidia-compilers/21.3 cuda
export F90=nvfortran

cat <<EOF
Usage :
* FIRST compile and run without GPU, once :
  $0 cpu
* compile and run with GPU+autocompare, as many times as needed :
  $0 verif
* compile and run with GPU, as many times as needed :
  $0 gpu
EOF

cd -P $(dirname $0)/..

cmd_$1
