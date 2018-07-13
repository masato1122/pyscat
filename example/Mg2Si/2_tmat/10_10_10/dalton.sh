#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=1:cnode012
#PBS -j oe
#PBS -N 10

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
nprocs=12
MPIRUN=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpirun

python cal_rscat.py

