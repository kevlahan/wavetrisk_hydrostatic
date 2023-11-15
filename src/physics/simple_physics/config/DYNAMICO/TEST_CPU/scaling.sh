#!/bin/bash

# these commands differ from their equivalent in TEST_GPU
function cmd_build()
{
    cd ..
    ./DYNAMICO.sh build_all X64_JEANZAY -no_io
}

function cmd_submit()
{
    for RUNDEF_NBP in $* ; do
	export RUNDEF_NBP
	sbatch job_X64_JEANZAY_scaling.sh
    done
}

# to avoid duplication the commands below are redirected to the same script in TEST_GPU

function redir_gpu()
{
    ../TEST_GPU/scaling.sh $*
}

function cmd_watch()
{
    redir_gpu watch $*
}

function cmd_grep()
{
    redir_gpu grep $*
}

function cmd_()
{
    redir_gpu
}

CMD=$1 ; shift
cmd_$CMD $*
