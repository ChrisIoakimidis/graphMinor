/*
This main function is used to test the Graph Minor function that computes the
minor of a Graph G given its cluster vector 'c' and its adjacency matrix A in
sparse form, specifically in CSR, since conversion between COO and CSR is
relatively fast.
*/
#include "graphMinor.h"
#include "myMatrix.h"
#include "myStructs.h"
#include <omp.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#define L_ID 1000

double get_time(struct timeval a, struct timeval b);

int main(int argc, char *argv[]) {
  // Tests
  COO A_coo;
  CSR A_csr;
  struct timeval a, b;
  struct timezone tz;
  A_coo = mm2coo(argv[1], 0, 0);
  A_csr = coo2csr(A_coo);
  CSR M;
  int N = A_csr.N;
  int nz = A_csr.nz;
  // Random id vector
  int *c = parse_cluster(argv[1], A_csr.N);
  printf("Data:\nN: %d\nnz: %d\n", N, nz);
  gettimeofday(&a, &tz);
  M = graphMinor(A_csr, c);
  gettimeofday(&b, &tz);

  double t = get_time(a, b);
  printf("Time: %f\n", t);
  printf("NZ: %d\n", M.nz);

  freeCOO(A_coo);
  freeCSR(A_csr);
  free(c);
  return 0;
}

double get_time(struct timeval a, struct timeval b) {
  return (double)((b.tv_sec - a.tv_sec) * 1000000 + b.tv_usec - a.tv_usec) /
         1000000;
}
