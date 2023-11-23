#ifndef GRAPH_MINOR_H
#define GRAPH_MINOR_H

#include "myMatrix.h"
#include "myStructs.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define P 128

/* Main function that does 3 things:
 * 1. Divides A.ii to P groups and for each group counts number of elements that
 * map to each row and stores it in its specific row_pointer
 *   2. Divide the row_pointers (there are L of them) and for each portion
 *      compress the elements
 *   3. Create M.ii's count array using the row_pointers compressed vectors
 *   4. Transfer all the vectors to the M.jv
 */
CSR graphMinor(CSR A, int *c);

int find_c_max(int *c, int N);

// The ii array is scanned and all positions are mapped to each row pointer
void count_rows(CSR A, int *c, Row_pointer *r_pointers, int cpu_count);

// Compress the elements of the pointers given to the structures' cRow
// L is the dimension of the minor array
// Also sets the count of the row for the count array of M
void r_pointer_to_compress(CSR A, int *c, Row_pointer *r, int L,
                           int *minor_count_array, FRow *fRow);

void compress_all(CSR A, int *c, Row_pointer *r, int L, int *minor_count_array);

// Just a function to parse the cluster vectors
// Returns pointer to cluster vector
// FILENAME WITHOUT EXTENSION
int *parse_cluster(char *filename, int N);

#endif
