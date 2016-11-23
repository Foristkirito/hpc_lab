#!/bin/bash
#SBATCH -J mpi_test
#SBATCH -n 3
#SBATCH -N 3
#SBATCH -c 4
#SBATCH -w cn02,cn03,cn05
#SBATCH -t 00:30:00
#SBATCH -o out
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./stencil 100 100 100 10 
