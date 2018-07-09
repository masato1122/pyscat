#!/bin/sh
#QSUB -queue F4cpu
#QSUB -node 1
#QSUB -mpi 2
#QSUB -omp 12
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=05:00:00
#PBS -N imp007

cd $PBS_O_WORKDIR
nprocs=2
mpirun -np ${nprocs} vasp > log.txt

