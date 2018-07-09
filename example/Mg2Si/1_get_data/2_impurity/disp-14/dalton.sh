#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=16
#PBS -j oe
#PBS -N imp014

export OMP_NUM_THREADS=16
cd $PBS_O_WORKDIR
nprocs=1
MPIRUN=/opt/intel/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpirun


#----- RUN the JOB -----#
$MPIRUN -np ${nprocs} vasp > log.txt

