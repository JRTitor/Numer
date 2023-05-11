#include "s21_matrix.h"
// variant 23
int main () {
  matrix_t A, B, X;
  double eps = 0.00001;
  int iter = 0;

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = -24.;
  A.matrix[0][1] = -6.;
  A.matrix[0][2] = 4.;
  A.matrix[0][3] = 7.;

  A.matrix[1][0] = -8.;
  A.matrix[1][1] = 21.;
  A.matrix[1][2] = 4.;
  A.matrix[1][3] = -2.;

  A.matrix[2][0] = 6.;
  A.matrix[2][1] = 6.;
  A.matrix[2][2] = 16.;
  A.matrix[2][3] = 0.;

  A.matrix[3][0] = -7.;
  A.matrix[3][1] = -7.;
  A.matrix[3][2] = 5.;
  A.matrix[3][3] = 24.;

  s21_create_matrix(4, 1, &B);
  B.matrix[0][0] = 130.;
  B.matrix[1][0] = 139.;
  B.matrix[2][0] = -84.;
  B.matrix[3][0] = -165.;

  printf("A:\n");
  Print_mat(&A);
  printf("B:\n");
  Print_mat(&B);
  printf("eps = %g\n", eps);

  printf("iteration:\n");
  iteration(&A, &B, eps, &X, &iter);
  printf("ITERATIONS: %d\n", iter);
  printf("X = \n");
  Print_mat(&X);
  s21_remove_matrix(&X);

  printf("seidel:\n");
  seidel(&A, &B, eps, &X, &iter);
  printf("ITERATIONS: %d\n", iter);
  printf("X = \n");
  Print_mat(&X);
  s21_remove_matrix(&X);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

