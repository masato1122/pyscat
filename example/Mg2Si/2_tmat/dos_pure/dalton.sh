#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=1:cnode012
#PBS -j oe
#PBS -N dos_pure

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
MPIRUN=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpirun

python cal_rscat.py

