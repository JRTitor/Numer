#include "s21_matrix.h"
// variant 23
int main () {
    matrix_t A, B, X;
    s21_create_matrix(5, 3, &A);
    s21_create_matrix(1, 5, &B);
    A.matrix[0][0] = 0;
    A.matrix[0][1] = 7;
    A.matrix[0][2] = -5;

    A.matrix[1][0] = -6;
    A.matrix[1][1] = 19;
    A.matrix[1][2] = -9;

    A.matrix[2][0] = 6;
    A.matrix[2][1] = -18;
    A.matrix[2][2] = 7;


    A.matrix[3][0] = -7;
    A.matrix[3][1] = -11;
    A.matrix[3][2] = -2;

    A.matrix[4][0] = 5;
    A.matrix[4][1] = -7;
    A.matrix[4][2] = 0;

    B.matrix[0][0] = 38;
    B.matrix[0][1] = 14;
    B.matrix[0][2] = -45;
    B.matrix[0][3] = 30;
    B.matrix[0][4] = 48;
    printf("\nA with 3 diagonals = \n");
    Print_mat(&A);
    printf("\nB^T = \n");
    Print_mat(&B);
    
    QP(&A, &B, &X);
    printf("\nX^T = \n");
    Print_mat(&X);
    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
    s21_remove_matrix(&X);
}

