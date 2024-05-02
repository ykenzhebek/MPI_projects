#!/bin/bash

#PBS -q S3792025
#PBS -l walltime=00:09:00
#PBS -l select=2:ncpus=8:mpiprocs=2:mem=8gb
#PBS -m n

cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')

echo "Number of MPI process: $MPI_NP"
echo 'File $PBS_NODEFILE:'
cat $PBS_NODEFILE
echo

mpirun -machinefile $PBS_NODEFILE -np $MPI_NP ./hybrid.o 4