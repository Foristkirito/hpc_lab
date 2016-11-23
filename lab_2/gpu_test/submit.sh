#!/bin/bash
#SBATCH -J gpu_test
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -w cn15
#SBATCH -t 00:30:00
#SBATCH -o out
#SBATCH -p gpu
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact ./a.out
