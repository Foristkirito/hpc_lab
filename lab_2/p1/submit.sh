#!/bin/bash
#SBATCH -J foristkirito
#SBATCH -n 7
#SBATCH -N 7
#SBATCH -c 24
#SBATCH -w cn02,cn03,cn04,cn05,cn06,cn07,cn08
#SBATCH -t 00:1:00
#SBATCH -o out_9
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./stencil 50 50 50 100
