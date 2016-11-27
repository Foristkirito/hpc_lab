#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

__global__ void test(float *A)
{
  int id = threadIdx.x + blockIdx.x*blockDim.x;
  if(id == 0)
    A[0] = 233333;
  return;
}
int main(int argc, char **argv)
{
  char hostname[100];
  gethostname(hostname, 100); //Get my host name
  cout << hostname << " Hello World" << endl;
  float *A = NULL;
  float *B = NULL;
  A = (float*)malloc(sizeof(float));
  cudaMalloc(&B, sizeof(float));
  test<<<1,1>>>(B);
  cudaMemcpy(A, B, sizeof(float), cudaMemcpyDeviceToHost);
  cout << A[0] << endl;
  free(A);
  cudaFree(B);
  return 0;
}
