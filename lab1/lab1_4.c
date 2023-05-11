#include "s21_matrix.h"
// variant 23
int main() {
  double eps = 0.1;
  int iter = 0;
  matrix_t A, euginvalues, euginvalue_vec;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 9.;
  A.matrix[0][1] = -5.;
  A.matrix[0][2] = -6.;

  A.matrix[1][0] = -5.;
  A.matrix[1][1] = 1.;
  A.matrix[1][2] = -8.;

  A.matrix[2][0] = -6.;
  A.matrix[2][1] = -8.;
  A.matrix[2][2] = -3.;
  printf("matrix A:\n");
  Print_mat(&A);
  printf("eps = %g\n", eps);
  rotate(&A, &euginvalue_vec, &euginvalues ,eps, &iter);
  printf("iterations: %d\n", iter);
  printf("euginvalues\n");
  Print_mat(&euginvalues);
  printf("euginvalue_vectors\n");
  Print_mat(&euginvalue_vec);
  s21_remove_matrix(&A);
  s21_remove_matrix(&euginvalue_vec);
  s21_remove_matrix(&euginvalues);
}

