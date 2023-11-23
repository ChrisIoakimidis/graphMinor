// This script is used to calculate the time needed for the graph Minor
// algorithm
// Argv[1] contains the matrix name (no suffix)
// Output P A.N A.nz M.N M.nz t_mix t_max t_median
#include "graphMinor.h"
#include "myMatrix.h"
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#define TEST_NUM 7

double get_time(struct timeval a, struct timeval b);

void insert_sort(double *a, int N);

int main(int argc, char *argv[]) {
  struct timeval a, b;
  struct timezone tz;
  // Get problem data
  COO A_coo = mm2coo(argv[1], 0, 0);
  CSR A_csr = coo2csr(A_coo);
  int *c = parse_cluster(argv[1], A_csr.N);
  CSR M;
  // Get times
  double t[TEST_NUM];
  for (int k = 0; k < TEST_NUM; k++) {
    gettimeofday(&a, &tz);
    M = graphMinor(A_csr, c);
    gettimeofday(&b, &tz);
    t[k] = get_time(a, b);
  }
  // Sort times
  insert_sort(t, TEST_NUM);
  // Prints the line that is pipelined to a txt file
  printf("%d %d %d %d %d %f %f %f\n", P, A_csr.N, A_csr.nz, M.N, M.nz, t[0],
         t[TEST_NUM - 1], t[TEST_NUM / 2]);
}

double get_time(struct timeval a, struct timeval b) {
  return (double)((b.tv_sec - a.tv_sec) * 1000000 + b.tv_usec - a.tv_usec) /
         1000000;
}

void insert_sort(double *a, int N) {
  int k, l;
  double temp;
  for (k = 1; k < N; k++) {
    temp = a[k];
    l = k - 1;
    while (l >= 0 && a[l] > temp) {
      a[l + 1] = a[l];
      l--;
    }
    a[l + 1] = temp;
  }
}
