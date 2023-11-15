#!/bin/bash
## Request name
#MSUB -r DYNAMICO_phyparam_mpi
## Number of tasks (=MPI processes) to use
#MSUB -n 40
## Number of OpenMP threads
#MSUB -c 1
## Elapsed time limit in seconds
#MSUB -T 300
# account, partition (Xeon/KNL)
#MSUB -A gen0239
#MSUB -q skylake
# Mount  on compute node
#MSUB -m work
## Quality of Service required (long [3 days], normal [1 day], test [30 min])
#MSUB -Q test

export OMP_NUM_THREADS=1
export OMP_STACKSIZE=128M

cd ${BRIDGE_MSUB_PWD} 

# module c++/gnu is needed to link with XIOS
# module feature/openmpi/mpi_compiler/intel is made necessary by c++/gnu

module purge
module load feature/openmpi/mpi_compiler/intel
module load c++/gnu/7.3.0
module load intel/17.0.6.256
module load mpi/openmpi/2.0.4
module load flavor/hdf5/parallel
module load netcdf-fortran/4.4.4
module load mkl/17.0.6.256
module load hdf5/1.8.20
module load parmetis/4.0.3

rm -rf gcm.log logs *.nc netcdf

date > gcm.log
ulimit -s unlimited
ccc_mprun ./gcm.exe >> gcm.log
date >> gcm.log

mkdir -p netcdf
cp gcm.log *.def netcdf
mv *.nc netcdf

mkdir -p logs
cp *.xml logs
mv xios_client_*.err xios_client_*.out gcm.log logs
