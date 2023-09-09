#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>


#define SUCCESS 1
#define FAILURE 0
#define OK 0
#define INCORRECT 1
#define UNABLE 2
#define s21_NAN 0.0 / 0.0
#define s21_INF 1.0 / 0.0
#define M_PI 3.14159265358979323846

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
void s21_copy_without_lines(int a, int b, matrix_t *A, matrix_t *tmp);
void LU(matrix_t* A, matrix_t* LU, int show_lu);
void Print_mat(matrix_t* A);
void Solve(matrix_t *A, matrix_t *B, matrix_t *X, int show_lu);
double det(matrix_t *A);
void Inverse(matrix_t *A, matrix_t *res);
void QP(matrix_t *A, matrix_t *B, matrix_t *X);
void noralize(matrix_t *A, matrix_t *B, matrix_t *a, matrix_t *b);
void iteration(matrix_t *A, matrix_t *B, double eps, matrix_t *res, int *iter);
void seidel(matrix_t *A, matrix_t *B, double eps, matrix_t *res, int *iter);
void rotate(matrix_t *A, matrix_t *eigenvalue_vec, matrix_t *eigenvalue, double eps, int *iter);
void household(matrix_t *A, matrix_t *Q, matrix_t *R);
void eigenvalue(matrix_t *A, matrix_t *eigenvalue,double eps, int *iter);
void Solve_compact(matrix_t *A, matrix_t *B, matrix_t *X);
double det_(matrix_t *A);
void Inverse_c(matrix_t *A, matrix_t *res);
