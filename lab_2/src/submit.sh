#!/bin/bash
#SBATCH -J foristkirito
#SBATCH -n 8
#SBATCH -N 8
#SBATCH -c 10
#SBATCH -w cn03,cn06,cn07,cn08,cn09,cn10,cn11,cn12
#SBATCH -t 00:1:00
#SBATCH -o out_9
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./stencil 1200 1200 1200 100
