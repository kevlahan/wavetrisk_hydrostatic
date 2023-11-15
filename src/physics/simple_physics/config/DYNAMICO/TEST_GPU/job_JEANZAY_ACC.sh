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
#SBATCH -A wuu@gpu
## Elapsed time limit HH:MM:SS
#SBATCH --time=00:10:00

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
module load nvidia-nsight-systems/2021.1.1
module list

# set up execution directory
rm -rf rundir
mkdir rundir
cp *.def *.xml rundir/

# and run
cd rundir
export TMPDIR=$JOBSCRATCH
ln -s $JOBSCRATCH /tmp/nvidia

echo "Run started : $(date)" > gcm.log
ulimit -s unlimited
srun --unbuffered nsys profile -t nvtx,openacc -b dwarf -o profile_${CI_COMMIT_SHORT_SHA}_%q{SLURM_PROCID} ../../modeles/DYNAMICO_phyparam/bin/dynamico_phyparam.exe >> gcm.log

echo "Run finished at $(date), now collecting profiling data (takes a couple of minutes)" >> gcm.log
nsys stats profile_${CI_COMMIT_SHORT_SHA}_0.qdrep >> gcm.log
echo "Finished : $(date)" >> gcm.log
