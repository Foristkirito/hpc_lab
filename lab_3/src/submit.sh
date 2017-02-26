#!/bin/bash
#SBATCH -J foristkirito
#SBATCH -n 8
#SBATCH -N 8
#SBATCH -c 24
#SBATCH -w cn02,cn03,cn07,cn08,cn09,cn10,cn11,cn12
#SBATCH -t 00:1:00
#SBATCH -o out_9_case2
#SBATCH -p batch
unset I_MPI_PMI_LIBRARY
mpiexec.hydra -bootstrap slurm -l \
  -genv KMP_AFFINITY compact \
 ./main 4 2 1 720 360 38 ./data/case_2bin/data_A_v1.bin ./data/case_2bin/data_x0_v1.bin ./data/case_2bin/data_b_v1.bin