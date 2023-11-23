#include "graphMinor.h"
#include "myMatrix.h"
#include "myStructs.h"
#include <cilk/cilk.h>

CSR graphMinor(CSR A, int *c) {
  CSR M;
  // Get size of minor
  M.N = find_c_max(c, A.N);
  // Allocate r_pointers
  Row_pointer *r_pointers = (Row_pointer *)malloc(M.N * sizeof(Row_pointer));
  cilk_scope {
    cilk_for(int k = 0; k < M.N; k++) {
      r_pointers[k] = r_pointer_init(k, P, A.N);
    }
  }
  cilk_sync;
  // Count rows
  count_rows(A, c, r_pointers, P);
  // Init count array
  int *minor_count_array = (int *)calloc(M.N, sizeof(int));
  // Compress all rows
  compress_all(A, c, r_pointers, M.N, minor_count_array);
  cilk_sync;
  // Create M.ii
  M.ii = (int *)malloc((M.N + 1) * sizeof(int));
  M.ii[0] = 0;
  for (int k = 0; k < M.N; k++) {
    M.ii[k + 1] = M.ii[k] + minor_count_array[k];
  }
  cilk_sync;
  M.nz = M.ii[M.N];
  M.jv = (JV *)malloc(M.nz * sizeof(JV));
  // Transfer all compressed rows to M
  cilk_scope {
    cilk_for(int k = 0; k < M.N; k++) {
      cRow_transfer(&r_pointers[k].compressed, &M);
    }
  }
  // Free useless things
  free(r_pointers);
  return M;
}

int find_c_max(int *c, int N) {
  int max = 0;
  int local_max[P];
  cilk_scope {
    cilk_for(int id = 0; id < P; id++) {
      local_max[id] = 0;
      for (int k = id; k < N; k += P)
        if (local_max[id] < c[k])
          local_max[id] = c[k];
    }
  }
  for (int id = 0; id < P; id++)
    if (local_max[id] > max)
      max = local_max[id];
  return max + 1;
}

void count_rows(CSR A, int *c, Row_pointer *r_pointers, int cpu_count) {
  // Assumes I'm thread id and there are cpu_count threads in total
  cilk_scope {
    cilk_for(int id = 0; id < P; id++) {
      for (int k = id; k < A.N; k += P) {
        // If row not empty
        if (A.ii[k] < A.ii[k + 1]) {
          list_insert(&r_pointers[c[k]].indices[id], k);
        }
      }
    }
  }
}
// Parallel compressing
void compress_all(CSR A, int *c, Row_pointer *r, int L,
                  int *minor_count_array) {
  FRow fRow[P];
  cilk_scope {
    cilk_for(int id = 0; id < P; id++) {
      fRow[id] = fRow_init(L, -1, c);
      for (int k = id; k < L; k += P) {
        r_pointer_to_compress(A, c, &r[k], L, minor_count_array, &fRow[id]);
      }
      fRow_free(&fRow[id]);
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
