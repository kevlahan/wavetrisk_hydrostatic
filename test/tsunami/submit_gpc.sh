#!/bin/bash
#PBS -l nodes=8:ppn=8:m32g,walltime=24:00:00
#PBS -N tsunami
cd $PBS_O_WORKDIR
module load openmpi
mpirun -np 64 ./tsunami
