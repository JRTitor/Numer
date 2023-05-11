#include "s21_matrix.h"
// variant 23
int main() {
  double eps = 0.01;
  matrix_t A, eigenvalues;
  int iter = 0;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1.;
  A.matrix[0][1] = 5.;
  A.matrix[0][2] = -6.;

  A.matrix[1][0] = 9.;
  A.matrix[1][1] = -7.;
  A.matrix[1][2] = -9.;

  A.matrix[2][0] = 6.;
  A.matrix[2][1] = -1.;
  A.matrix[2][2] = -9.;

  printf("matrix A:\n");
  Print_mat(&A);
  printf("eps = %g\n", eps);
  eigenvalue(&A, &eigenvalues, eps, &iter);
  printf("ITERATIONS: %d\n", iter);
  printf("eigenvalus\n");
  Print_mat(&eigenvalues);

  s21_remove_matrix(&A);
  s21_remove_matrix(&eigenvalues);
}

