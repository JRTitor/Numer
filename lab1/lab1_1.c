#include "s21_matrix.h"

// variant 23

int main() {
  matrix_t A, B, X;

  s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = -7;
  A.matrix[0][2] = 8;
  A.matrix[0][3] = -4;

  A.matrix[1][0] = 0;
  A.matrix[1][1] = -1;
  A.matrix[1][2] = 4;
  A.matrix[1][3] = -1;

  A.matrix[2][0] = 3;
  A.matrix[2][1] = -4;
  A.matrix[2][2] = 2;
  A.matrix[2][3] = -1;

  A.matrix[3][0] = -9;
  A.matrix[3][1] = 1;
  A.matrix[3][2] = -4;
  A.matrix[3][3] = -6;

  s21_create_matrix(1, 4, &B);
  B.matrix[0][0] = 57;
  B.matrix[0][1] = 24;
  B.matrix[0][2] = 28;
  B.matrix[0][3] = 12;
  printf("A:\n");
  Print_mat(&A);
  printf("B^T:\n");
  Print_mat(&B);

  printf("det(A) = %g\n",det_(&A));
  
  printf("X^T = \n");
  Solve_compact(&A, &B, &X);
  Print_mat(&X);
  s21_remove_matrix(&X);

  printf("A^(-1):\n");
  Inverse_c(&A, &X);
  Print_mat(&X);
  s21_remove_matrix(&X);

  s21_remove_matrix(&B);
  s21_remove_matrix(&A);
  
}

