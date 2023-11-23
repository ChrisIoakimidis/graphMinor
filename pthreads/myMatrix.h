#ifndef MYMATRIX_H
#define MYMATRIX_H

// Define location of matrix data
#define MM_PATH "/home/chris/Documents/pads/mm/"
#define READ_PATH_S "read/sorted/"
#define CLUSTER_PATH "data/clusters/"
#define READ_PATH_U "data/datasets/"

#include "mmio.h"
#include <stdbool.h>
#include <string.h>

// Define matrix and point structures
// For COO: Edge [(i,j,v) pair]
typedef struct {
  int i;
  int j;
  int v;
} Edge;

// For CSR: jv
typedef struct {
  int j;
  int v;
} JV;

// Row of Matrix
typedef struct {
  int i;
  JV *jvs;
  int cnt; // ammount of jvs
} ROW;

// Actual CSR matrix (index array and Col containing j's and v's)
typedef struct {
  int N;
  int nz;
  int *ii;
  JV *jv;
} CSR;

// COO Matrix
typedef struct {
  int N;
  int nz;
  Edge *el;
} COO;

// Free functions for dynamically allocated matrices
void freeCOO(COO A);
void freeCSR(CSR A);

// This function parses the .mtx format using mmio (and a possible permutation
// matrix) and returns the COO form of the matrix
COO mm2coo(char *source, bool has_perm, bool save);

// Convert a matrix from COO to CSR (parallelization is assumed after
// this has been done)
CSR coo2csr(COO A_coo);

#endif // !MYMATRIX_H
