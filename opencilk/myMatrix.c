// This fils implements the functions for matrix manipulation, conversion
// I/O to files etc
#include "myMatrix.h"
#include "mmio.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>

/* Parses the matrix (in mm format) and shuffles it (and saves it if prompted)
 * Arguments:
 * source: Name of matrix (without extensions)
 * has_perm: If true, permutate the matrix according to its permutation file
 * save: If true, save the matrix in mm format in the UNSORTED section
 */
COO mm2coo(char *source, bool has_perm, bool save) {
  // Allocate size and return pointer
  Edge *A_el;
  int M, N, nz;
  int l = strlen(source);
  MM_typecode matcode;
  int *perm_indices;

  // Define paths to permutation and mtx
  char *mat = (char *)malloc(200 * sizeof(char));
  char *perm = (char *)malloc(200 * sizeof(char));
  // If permutation is required, access file from sorted list
  if (has_perm) {
    sprintf(mat, "%s%s%s%s", MM_PATH, READ_PATH_S, source, ".mtx");
    sprintf(perm, "%s%s%s%s", MM_PATH, READ_PATH_S, source, ".permutation.txt");
  } else {
    sprintf(mat, "%s%s%s%s", MM_PATH, READ_PATH_U, source, ".mtx");
    sprintf(perm, "%s%s%s%s", MM_PATH, READ_PATH_U, source, ".permutation.txt");
  }
  // Read permutation
  errno = 0;
  FILE *matFile = fopen(mat, "r");
  if (matFile == NULL) {
    printf("ERROR IN READING\n");
    printf("Error: %d\n", errno);
    exit(errno);
  }
  // Check if matrix is readable
  if (mm_read_banner(matFile, &matcode)) {
    printf("Could not process Martix banner.\n");
    exit(1);
  }
  if (!(mm_is_matrix(matcode) && mm_is_coordinate(matcode))) {
    printf("Matrix is not in the desired form of MCIG.\n");
    exit(1);
  }
  // Get sizes
  if (mm_read_mtx_crd_size(matFile, &M, &N, &nz))
    exit(1);
  if (mm_is_symmetric(matcode))
    nz *= 2;
  if (M != N) {
    printf("Matrix is not square, extending smaller dimension to make it "
           "square.\n");
    if (M > N) {
      N = M;
    } else {
      M = N;
    }
  }
  // If a permutation has been set do it
  if (has_perm) {
    // Read file
    FILE *permFile = fopen(perm, "r");
    perm_indices = (int *)malloc(nz * sizeof(int));
    int buffer;
    int k = 0;
    // Parse list into array
    while (!feof(permFile)) {
      fscanf(permFile, " %d ", &buffer);
      perm_indices[k++] = buffer;
    }
    fclose(permFile);

    // Allocate Matrix
    A_el = (Edge *)malloc(nz * sizeof(Edge));
    // Parse matrix to COO pointer
    int i, j, v;
    k = 0;
    // Store in the random indices determined by permutation.
    while (!feof(matFile)) {
      fscanf(matFile, "%d %d %d", &i, &j, &v);
      A_el[perm_indices[k]].i = i - 1;
      A_el[perm_indices[k]].j = j - 1;
      A_el[perm_indices[k++]].v = v;
    }
    free(perm_indices);
  } else {
    A_el = (Edge *)malloc(nz * sizeof(Edge));
    // Parse matrix to COO pointer
    int i, j, v;
    int k = 0;
    // Store in the random indices determined by permutation.
    while (!feof(matFile)) {
      fscanf(matFile, "%d %d %d", &i, &j, &v);
      if (mm_is_pattern(matcode))
        v = 1;
      A_el[k].i = i - 1;
      A_el[k].j = j - 1;
      A_el[k++].v = v;
      if (mm_is_symmetric(matcode)) {
        A_el[k].j = i - 1;
        A_el[k].i = j - 1;
        A_el[k++].v = v;
      }
    }
  }
  // Save if specified (only makes sense for having a perm so save to unsorted)
  if (save) {
    char *dest = (char *)malloc(200 * sizeof(char));
    sprintf(dest, "%s%s%s%s", MM_PATH, READ_PATH_U, source, ".mtx");
    FILE *destFile = fopen(dest, "w");
    mm_write_banner(destFile, matcode);
    fprintf(destFile, "%%\n");
    mm_write_mtx_crd_size(destFile, M, N, nz);
    for (int k = 0; k < nz; k++) {
      fprintf(destFile, "%d %d %d\n", A_el[k].i, A_el[k].j, A_el[k].v);
    }
  }
  fclose(matFile);
  free(perm);
  free(mat);
  COO A;
  A.N = N;
  A.nz = nz;
  A.el = A_el;
  return A;
}

void freeCOO(COO A) {
  // Firstly free edge matrix then itself
  free(A.el);
}

void freeCSR(CSR A) {
  free(A.ii);
  free(A.jv);
}

/* Converts a matrix given in COO format to CSR format (sequencially)
 *  Arguments:
 *  Edge: COO Matrix
 *  Returns:
 *  CSR* pointer to CSR form (MUST BE FREED BY CALLER)
 */
CSR coo2csr(COO A_coo) {
  CSR A;
  // Pass size and nz
  int N = A.N = A_coo.N;
  int nz = A.nz = A_coo.nz;
  A.ii = malloc((N + 1) * sizeof(int));
  A.jv = malloc(nz * sizeof(JV));
  Edge *el = A_coo.el; // Allocate counter for index array
  int *i_cnt = (int *)calloc(N, sizeof(int));
  // Count all rows
  for (int k = 0; k < nz; k++) {
    i_cnt[el[k].i]++;
  }
  // Prefix sum copied to ii data
  A.ii[0] = 0;
  for (int k = 0; k < N; k++) {
    A.ii[k + 1] = A.ii[k] + i_cnt[k];
  }
  // Dump data from COO to CSR
  int i, cnt;
  for (int k = 0; k < nz; k++) {
    i = el[k].i;
    cnt = i_cnt[i];
    if (cnt > 0) {
      A.jv[A.ii[i] + cnt - 1].j = el[k].j;
      A.jv[A.ii[i] + cnt - 1].v = el[k].v;
      i_cnt[i]--;
    }
  }
  free(i_cnt);
  return A;
}
