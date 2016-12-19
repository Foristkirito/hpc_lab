#!/bin/bash
#SBATCH -J foristkirito
#SBATCH -n 10
#SBATCH -N 10
#SBATCH -c 2
#SBATCH -w cn03,cn04,cn05,cn06,cn07,cn08,cn09,cn10,cn11,cn12
#SBATCH -t 00:1:00
#SBATCH -o out_9
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./main 2 5 1 360 180 38 ./data/case_1bin/data_A_v1.bin ./data/case_1bin/data_x0_v1.bin ./data/case_1bin/data_b_v1.bin