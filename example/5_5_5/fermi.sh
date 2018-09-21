#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=24
#PBS -j oe
#PBS -N test

export LANG=C
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
nprocs=24

cd $PBS_O_WORKDIR
mpirun -np ${nprocs} vasp > log.txt

