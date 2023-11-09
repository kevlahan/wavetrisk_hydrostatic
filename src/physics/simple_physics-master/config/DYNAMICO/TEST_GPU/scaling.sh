#!/bin/bash

function cmd_build()
{
    cd ..
    ./DYNAMICO.sh build_all JEANZAY_NVIDIA_ACC -no_io
}

function cmd_submit()
{
    for RUNDEF_NBP in $* ; do
	export RUNDEF_NBP
	sbatch job_JEANZAY_ACC_scaling.sh
    done
}

function cmd_watch()
{
    tail -f rundir.$1/gcm.log
}

function cmd_grep()
{
    for RUNDEF_NBP in $* ; do
	LOG=rundir.$RUNDEF_NBP/gcm.log
	for str in 'Whole ' 'Throughput' 'It No' 'MPI_waitall' 'dyn ' 'phys ' 'phyparam '  ; do 
	    grep "$str" $LOG | tail -n 1
	done
    done
}

function cmd_()
{
    cat <<EOF
Usage : 
      $0 build
      $0 submit NBP1 NBP2 NBP3 ...
      $0 watch NBP
      $0 grep NBP1 NBP2 NBP3 ...
EOF
}

CMD=$1 ; shift
cmd_$CMD $*
