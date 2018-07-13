#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -N Mg2Si_grn

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
nprocs=12
MPIRUN=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpirun

python test_green.py

