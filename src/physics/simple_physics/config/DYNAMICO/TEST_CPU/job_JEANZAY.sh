#!/bin/bash
## Request name
#SBATCH --job-name=SimPhy
## Number of tasks (=MPI processes) to use
#SBATCH --ntasks=40
## Number of OpenMP threads
#SBATCH --cpus-per-task=1

## computing project
#SBATCH -A wuu@cpu
## QoS for short jobs
#SBATCH --qos=qos_cpu-dev
## Elapsed time limit HH:MM:SS
#SBATCH --time=00:15:00

# do not use hyperthreading
#SBATCH --hint=nomultithread
# standard outputs
#SBATCH --output=DYNAMICO%j.out
#SBATCH --error=DYNAMICO%j.out

export OMP_NUM_THREADS=1
# OpenMP binding
export OMP_PLACES=cores

# stack
export OMP_STACKSIZE=128M
ulimit -s unlimited

# move to submission directory
cd ${SLURM_SUBMIT_DIR}

# load the same modules as during compilation
source ../modeles/DYNAMICO/arch.env
module list

# cleanup execution directory and run
rm -rf rundir
mkdir rundir
cp *.def *.xml rundir/
cd rundir
date > gcm.log
ulimit -s unlimited
# srun ../../modeles/DYNAMICO/bin/icosa_gcm.exe >> gcm.log
srun ../../modeles/DYNAMICO_phyparam/bin/dynamico_phyparam.exe >> gcm.log
date >> gcm.log
