#include "graphMinor.h"
#include "myMatrix.h"
#include "myStructs.h"
#include <omp.h>

CSR graphMinor(CSR A, int *c) {
  CSR M;
  // Get size of minor
  M.N = find_c_max(c, A.N);
  // Allocate r_pointers
  Row_pointer *r_pointers = (Row_pointer *)malloc(M.N * sizeof(Row_pointer));
#pragma omp parallel for
  for (int k = 0; k < M.N; k++) {
    r_pointers[k] = r_pointer_init(k, P, A.N);
  }
  // Count rows
  count_rows(A, c, r_pointers, P);
  // Init count array
  int *minor_count_array = (int *)calloc(M.N, sizeof(int));
  // Compress all rows
  compress_all(A, c, r_pointers, M.N, minor_count_array);
  // Create M.ii
  M.ii = (int *)malloc((M.N + 1) * sizeof(int));
  for (int k = 0; k < M.N; k++) {
    M.ii[k + 1] = M.ii[k] + minor_count_array[k];
  }
  M.nz = M.ii[M.N];
  M.jv = (JV *)malloc(M.nz * sizeof(JV));
// Transfer all compressed rows to M
#pragma omp parallel for
  for (int k = 0; k < M.N; k++) {
    cRow_transfer(&r_pointers[k].compressed, &M);
  }
  // Free useless things
#pragma omp parallel for
  for (int k = 0; k < M.N; k++) {
    r_pointer_free(&r_pointers[k]);
  }
  free(r_pointers);
  return M;
}

int find_c_max(int *c, int N) {
  int max = 0;
  int local_max;
#pragma omp parallel num_threads(P) private(local_max)
  {
    local_max = 0;
#pragma omp for
    for (int k = 0; k < N; k++) {
      if (local_max < c[k])
        local_max = c[k];
    }
#pragma omp critical
    {
      if (max < local_max) {
        max = local_max;
      }
    }
  }
  return max + 1;
}

void count_rows(CSR A, int *c, Row_pointer *r_pointers, int cpu_count) {
  // Assumes I'm thread is and there are cpu_count threads in total
  int id;
#pragma omp parallel num_threads(P) private(id)
  {
    id = omp_get_thread_num();
#pragma omp for
    for (int k = 0; k < A.N; k++) {
      // If row not empty
      if (A.ii[k] < A.ii[k + 1]) {
        list_insert(&r_pointers[c[k]].indices[id], k);
      }
    }
  }
}
// Parallel compressing
void compress_all(CSR A, int *c, Row_pointer *r, int L,
                  int *minor_count_array) {
  FRow fRow;
#pragma omp parallel num_threads(P) private(fRow)
  {
    fRow = fRow_init(L, -1, c);
#pragma omp for
    for (int k = 0; k < L; k++) {
      r_pointer_to_compress(A, c, &r[k], L, minor_count_array, &fRow);
    }
  }
}

// From A.jv to this.compressed
void r_pointer_to_compress(CSR A, int *c, Row_pointer *r, int L,
                           int *minor_count_array, FRow *fRow) {
  fRow->minor_i = r->minor_i;
  fRow->indices.cnt = 0;
  // Iterate through all counts of all cpus
  for (int k = 0; k < r->core_num; k++) {
    for (int l = 0; l < r->indices[k].cnt; l++) {
      for (int m = A.ii[r->indices[k].data[l]];
           m < A.ii[r->indices[k].data[l] + 1]; m++) {
        fRow_insert(fRow, A.jv[m]);
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
