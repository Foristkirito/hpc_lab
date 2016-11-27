#!/bin/bash
#SBATCH -J foristkirito
#SBATCH -n 2
#SBATCH -N 2
#SBATCH -c 4
#SBATCH -w cn08,cn09
#SBATCH -t 00:1:00
#SBATCH -o out_9
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./stencil 200 100 100 1
