#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv); //Init MPI
  int myrank, nprocs, nthreads, i, j;
  int cpus[100];
  char hostname[100];
  gethostname(hostname, 100); //Get my host name
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //Get Process ID
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Get # of Processes
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
    cpus[omp_get_thread_num()] = sched_getcpu();
  }
  for(i = 0; i < nprocs; i ++) {
    if(i == myrank) {
      printf("%s has %d threads\n", hostname, nthreads);
      for(j = 0; j < nthreads; j ++)
        printf("Need core %d\n", cpus[j]);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return 0;
}
