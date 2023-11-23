#include "myStructs.h"

// List methods
List list_init() {
  List out;
  out.cnt = 0;
  out.size = LIST_INIT_SIZE;
  out.data = (int *)malloc(out.size * sizeof(int));
  return out;
}

void list_free(List *list) {
  free(list->data);
  list->cnt = -1;  // mark it as dead
  list->size = -1; // mark it as dead
}

void list_insert(List *list, int element) {
  if (list == NULL) {
    printf("ERROR\n");
    getchar();
  }
  // If at size limit, double
  if (list->cnt + 1 == list->size) {
    int *temp = (int *)malloc(2 * list->size * sizeof(int));
    for (int k = 0; k < list->cnt; k++) {
      temp[k] = list->data[k];
    }
    free(list->data);
    list->data = temp;
    list->size = 2 * list->size;
  }
  // Add element to list
  list->data[list->cnt] = element;
  list->cnt = list->cnt + 1;
}

// CRow methods

CRow cRow_init(int size, int minor_i) {
  CRow out;
  out.minor_i = minor_i;
  out.data = (JV *)malloc(size * sizeof(JV));
  out.size = size;
  return out;
}

void cRow_free(CRow *cRow) {
  free(cRow->data);
  cRow->size = -1;
}

void cRow_transfer(CRow *cRow, CSR *A) {
  // Transfer
  int m = 0;
  for (int k = A->ii[cRow->minor_i]; k < A->ii[cRow->minor_i + 1]; k++) {
    A->jv[k] = cRow->data[m++];
  }
  // Free cRow
  cRow_free(cRow);
}

// FRow methods

FRow fRow_init(int size, int minor_i, int *c) {
  FRow out;
  out.size = size;
  out.indices = list_init();
  out.minor_i = minor_i;
  out.c = c;
  out.data = (JV *)malloc(size * sizeof(JV));
  for (int k = 0; k < size; k++)
    out.data[k].j = -1;
  return out;
}

void fRow_free(FRow *fRow) {
  fRow->size = -1;
  free(fRow->data);
  list_free(&fRow->indices);
}

void fRow_insert(FRow *fRow, JV element) {
  // Keep the val change the j
  JV temp = element;
  temp.j = fRow->c[element.j];
  // If it doesn't exist enter it
  if (fRow->data[temp.j].j == -1) {
    list_insert(&fRow->indices, temp.j);
    fRow->data[temp.j] = temp;
  } else {
    fRow->data[temp.j].v += temp.v;
  }
}

CRow fRow_compress(FRow *fRow) {
  CRow out = cRow_init(fRow->indices.cnt, fRow->minor_i);
  for (int k = 0; k < fRow->indices.cnt; k++) {
    out.data[k] = fRow->data[fRow->indices.data[k]];
    fRow->data[fRow->indices.data[k]].j = -1;
  }
  return out;
}

// R pointer methods

Row_pointer r_pointer_init(int minor_i, int core_num) {
  Row_pointer out;
  out.minor_i = minor_i;
  out.core_num = core_num;
  out.indices = (List *)malloc(out.core_num * sizeof(List));
  for (int k = 0; k < core_num; k++) {
    out.indices[k] = list_init();
  }
  return out;
}

// Assume compressed row has been transferd so it's already free
void r_pointer_free(Row_pointer *r_pointer) {
  for (int k = 0; k < r_pointer->core_num; k++) {
    list_free(&r_pointer->indices[k]);
  }
  free(r_pointer->indices);
}
