#!/bin/bash
#PBS -N DYNAMICO_phyparam_mpi
#PBS -q std
#PBS -n
#PBS -l nodes=1:ppn=60
#PBS -l walltime=04:00:00
#PBS -l mem=120gb
#PBS -l vmem=120gb

#####comment:PBS -l mem=31922327552,ncpus=60
#Jobs start in the HOME directory, cd to submitted directory
cd "$PBS_O_WORKDIR"

############################################################

echo ------------------JOB SUMMARY--------------------------
echo 'PBS OUT: Job allocated on: '${NCPU}' cpu(s)'
echo 'PBS OUT: Job running on the following node(s): '
cat nodelist.txt
echo -------------------------------------------------------
echo PBS OUT: Job running on: $PBS_O_HOST
echo PBS OUT: Originating queue: $PBS_O_QUEUE
echo PBS OUT: Executing queue: $PBS_QUEUE
echo PBS OUT: Work directory: $PBS_O_WORKDIR
echo PBS OUT: Execution mode: $PBS_ENVIRONMENT
echo PBS OUT: Job identifier: $PBS_JOBID
echo PBS OUT: Job name: $PBS_JOBNAME
echo PBS OUT: Node file: $PBS_NODEFILE
echo PBS OUT: Current home directory: $PBS_O_HOME
echo PBS OUT: PATH= $PBS_O_PATH
echo -------------------------------------------------------


export OMP_NUM_THREADS=1
export OMP_STACKSIZE=128M
ulimit -s unlimited

module purge
module load gnu/4.9.3 
module load intel/15.0.6.233
module load openmpi/1.6.5-ifort
module load hdf5/1.8.18-parallel-ifort
module load netcdf4/4.4.1.1-parallel-ifort
module list

rm -rf gcm.log logs *.nc netcdf
ls -lrth
which mpirun
mpirun -np 40 -cpus-per-proc $OMP_NUM_THREADS -bycore -bind-to-core -report-bindings ./gcm.exe > gcm.log 2>&1
date >> gcm.log

mkdir -p netcdf
cp gcm.log *.def netcdf
mv *.nc netcdf

mkdir -p logs
cp *.xml logs
mv xios_client_*.err xios_client_*.out gcm.log logs
