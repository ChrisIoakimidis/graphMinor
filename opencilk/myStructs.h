#ifndef MY_STRUCTS_H
#define MY_STRUCTS_H

#define LIST_INIT_SIZE 2
#include "myMatrix.h"

// Dynamic list
typedef struct {
  int size;
  int *data;
  int cnt;
} List;

// Compressed row (includes only non zero)
typedef struct {
  int minor_i; // row id;
  int size;    // Size of vector
  JV *data;    // data array
} CRow;

// Full row vector for each row to be compressed
typedef struct {
  int minor_i;  // row id;
  int size;     // vector size
  JV *data;     // data array
  List indices; // array of indices of non zero elements
  int *c;       // cluster vector for abstraction reasons
} FRow;

// For each minor row, you point to all the positions
// in the ii vector (so rows, that compress into the target row)
typedef struct {
  int minor_i;     // id of target row
  int core_num;    // # of threads for this process
  List *indices;   // a list of indices (1 for each core);
  CRow compressed; // compressed row to put all the data in once done
} Row_pointer;

// List methods

// Initialize the list to size 16 with no elements
List list_init();

// Frees the list
void list_free(List *list);

// Inserts element into list and increases size if needed
void list_insert(List *list, int element);

// cRow methods

// Initializes an empty row
CRow cRow_init(int size, int minor_i);

// Frees cRow
void cRow_free(CRow *cRow);

// Copy data to minor array, free cRow
void cRow_transfer(CRow *cRow, CSR *A);

// FRow methods

// Initializes fRow (all mapping to c is done inside the structure)
FRow fRow_init(int size, int minor_i, int *c);

// Free fRow
void fRow_free(FRow *fRow);

// Insert element to fRow (takes as input the original array's jv)
void fRow_insert(FRow *fRow, JV element);

// Compress to CRow and return the pointer to it. Free fRow
CRow fRow_compress(FRow *fRow);

// Row pointer methods

// Initializes the struct
Row_pointer r_pointer_init(int minor_i, int core_num, int N);

// Frees the pointer
void r_pointer_free(Row_pointer *r_pointer);

#endif
