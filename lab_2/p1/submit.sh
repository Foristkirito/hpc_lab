#!/bin/bash
#SBATCH -J foristkirito
#SBATCH -n 9
#SBATCH -N 9
#SBATCH -c 4
#SBATCH -w cn04,cn05,cn06,cn07,cn08,cn09,cn10,cn11,cn12
#SBATCH -t 00:1:00
#SBATCH -o out_9
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./stencil 1000 1000 1000 10
