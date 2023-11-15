#!/bin/bash
## Request name
#SBATCH --job-name=DCMIP41_mpi

# perime #SBATCH --partition=gpu_p1          # partition GPU choisie
# cf http://www.idris.fr/eng/jean-zay/gpu/jean-zay-gpu-exec_partition_slurm-eng.html

#SBATCH --nodes=1                  # nombre de noeud
#SBATCH --ntasks=4                 # nombre de tache MPI (= nombre de GPU ici)
#SBATCH --ntasks-per-node=4        # nombre de tache MPI par noeud (= nombre de GPU ici)
#SBATCH --gres=gpu:4               # nombre de GPU par noeud
#SBATCH --cpus-per-task=10         # nombre de coeurs CPU par tache

## computing project
#SBATCH -A stu@gpu
## Elapsed time limit HH:MM:SS
#SBATCH --time=00:30:00

# do not use hyperthreading
#SBATCH --hint=nomultithread
# standard outputs
#SBATCH --output=DYNAMICO%j.out
#SBATCH --error=DYNAMICO%j.out

set +x

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

# set up execution directory
RUNDIR=rundir.$RUNDEF_NBP
rm -rf $RUNDIR
mkdir $RUNDIR
cp *.xml $RUNDIR
cp run.$RUNDEF_NBP.def $RUNDIR/run.def 
# and run
cd $RUNDIR

pwd
ls-lrth
grep -A 10 Resolution gcm.log

echo "Run started : $(date)" | tee gcm.log
ulimit -s unlimited
srun ../../modeles/DYNAMICO_phyparam/bin/dynamico_phyparam.exe >> gcm.log
echo "Finished : $(date)" >> gcm.log
