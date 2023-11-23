#include "graphMinor.h"
#include "myMatrix.h"
#include "myStructs.h"
#include <pthread.h>
#include <stdio.h>

CSR graphMinor(CSR A, int *c) {
  CSR M;
  P_ARGS args[P];
  pthread_t threads[P];
  for (int k = 0; k < P; k++) {
    args[k].id = k;
    args[k].cpu_count = P;
    args[k].c = c;
    args[k].A = &A;
    args[k].mode = 1;
    pthread_create(&threads[k], NULL, parallel_mode, &args[k]);
  }
  for (int k = 0; k < P; k++) {
    pthread_join(threads[k], NULL);
  }
  // Find L
  M.N = 0;
  for (int k = 0; k < P; k++) {
    if (args[k].local_max > M.N)
      M.N = args[k].local_max;
  }
  M.N++;
  // Allocate count array and Row pointers
  Row_pointer *r_pointers = (Row_pointer *)malloc(M.N * sizeof(Row_pointer));
  for (int k = 0; k < M.N; k++) {
    r_pointers[k] = r_pointer_init(k, P);
  }
  // Count count
  for (int k = 0; k < P; k++) {
    args[k].r_pointers = r_pointers;
    args[k].mode = 2;
    pthread_create(&threads[k], NULL, parallel_mode, &args[k]);
  }
  for (int k = 0; k < P; k++) {
    pthread_join(threads[k], NULL);
  }
  int *minor_count_array = (int *)calloc(M.N, sizeof(int));
  for (int k = 0; k < P; k++) {
    args[k].M = &M;
    args[k].minor_count_array = minor_count_array;
    args[k].mode = 3;
    pthread_create(&threads[k], NULL, parallel_mode, &args[k]);
  }
  for (int k = 0; k < P; k++) {
    pthread_join(threads[k], NULL);
  }
  // Create M.ii
  M.ii = (int *)malloc((M.N + 1) * sizeof(int));
  for (int k = 0; k < M.N; k++) {
    M.ii[k + 1] = M.ii[k] + minor_count_array[k];
  }
  M.nz = M.ii[M.N];
  M.jv = (JV *)malloc(M.nz * sizeof(JV));

  for (int k = 0; k < P; k++) {
    args[k].mode = 4;
    pthread_create(&threads[k], NULL, parallel_mode, &args[k]);
  }

  for (int k = 0; k < M.N; k++)
    r_pointer_free(&r_pointers[k]);
  free(r_pointers);
  return M;
}

void *parallel_mode(void *args) {
  P_ARGS *a = (P_ARGS *)args;
  int id = a->id;
  int cpu_count = a->cpu_count;
  switch (a->mode) {
  case 1:
    a->local_max = find_c_max(a->c, a->A->N, id, cpu_count);
    break;
  case 2:
    count_rows(a->A, a->c, a->r_pointers, a->id, a->cpu_count);
    break;
  case 3:
    compress_all_rows(a->A, a->c, a->r_pointers, a->M->N, a->minor_count_array,
                      a->id, a->cpu_count);
    break;
  case 4:
    for (int k = id; k < a->M->N; k += a->cpu_count) {
      cRow_transfer(&a->r_pointers[k].compressed, a->M);
    }
    break;
  }

  return NULL;
}

int find_c_max(int *c, int N, int id, int cpu_count) {
  int max = 0;
  for (int k = id; k < N; k += cpu_count) {
    if (max < c[k])
      max = c[k];
  }
  return max;
}

void count_rows(CSR *A, int *c, Row_pointer *r_pointers, int id,
                int cpu_count) {
  // Assumes I'm thread id and there are cpu_count threads in total
  for (int k = id; k < A->N; k += cpu_count) {
    // If row not empty
    if (A->ii[k] < A->ii[k + 1]) {
      list_insert(&r_pointers[c[k]].indices[id], k);
    }
  }
}

void compress_all_rows(CSR *A, int *c, Row_pointer *r, int L,
                       int *minor_count_array, int id, int cpu_count) {
  FRow fRow = fRow_init(L, -1, c);
  // Iterating through whole row pointer array
  for (int k = id; k < L; k += cpu_count) {
    r_pointer_to_compress(A, c, &r[k], L, minor_count_array, &fRow);
  }
  fRow_free(&fRow);
}

// From A.jv to this.compressed
void r_pointer_to_compress(CSR *A, int *c, Row_pointer *r, int L,
                           int *minor_count_array, FRow *fRow) {
  fRow->indices.cnt = 0;
  fRow->minor_i = r->minor_i;
  // Iterate through all counts of all cpus
  for (int k = 0; k < r->core_num; k++) {
    for (int l = 0; l < r->indices[k].cnt; l++) {
      for (int m = A->ii[r->indices[k].data[l]];
           m < A->ii[r->indices[k].data[l] + 1]; m++) {
        fRow_insert(fRow, A->jv[m]);
      }
    }
  }
  // Compress fRow to the cRow
  r->compressed = fRow_compress(fRow);
  // Set count of minor count array
  minor_count_array[r->minor_i] = r->compressed.size;
}

int *parse_cluster(char *filename, int N) {
  // Init output
  int *c = (int *)malloc(N * sizeof(int));
  char *target = (char *)malloc(150 * sizeof(char));
  sprintf(target, "%s%s%s%s", MM_PATH, CLUSTER_PATH, filename, ".txt");
  FILE *file = fopen(target, "r");
  if (file == NULL) {
    printf("Error in reading cluster\n");
    exit(1);
  }
  int temp;
  int cnt = 0;
  while (!feof(file)) {
    fscanf(file, "%d\n", &temp);
    c[cnt++] = temp - 1;
  }
  return c;
}
