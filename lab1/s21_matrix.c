#include "s21_matrix.h"

int incorrect(matrix_t *A) {
  int res = OK;
  if (!A || A->rows <= 0 || A->columns <= 0 || abs(A->rows) == s21_INF ||
      abs(A->columns) == s21_INF || A->rows == s21_NAN || A->columns == s21_NAN)
    res = INCORRECT;
  return res;
}
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = OK;
  if (!result || rows < 0 || columns < 0) {
    res = INCORRECT;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
    }
  }
  return res;
}
void s21_remove_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    free(A->matrix[i]);
  }
  free(A->matrix);
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (incorrect(A) || incorrect(B)) {
    res = INCORRECT;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    res = UNABLE;
  } else {
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-6) {
          res = FAILURE;
          break;
        }
      }
      if (res == FAILURE) break;
    }
  }
  return res;
}
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (incorrect(A) || incorrect(B)) {
    res = INCORRECT;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    res = UNABLE;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return res;
}
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (incorrect(A) || incorrect(B)) {
    res = INCORRECT;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    res = UNABLE;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return res;
}
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = OK;
  if (incorrect(A)) {
    res = INCORRECT;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return res;
}
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (incorrect(A)) {
    res = INCORRECT;
  } else if (A->columns != B->rows) {
    res = UNABLE;
  } else {
    s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < B->columns; j++) {
        double sum = 0;
        for (int k = 0; k < A->columns; k++) {
          sum += (A->matrix[i][k]) * (B->matrix[k][j]);
        }
        result->matrix[i][j] = sum;
      }
    }
  }
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (incorrect(A)) {
    res = INCORRECT;
  } else {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return res;
}

void s21_copy_without_lines(int a, int b, matrix_t *A, matrix_t *tmp) {
  for (int i = 0, k = 0; i < A->rows; i++) {
    int l = 0;
    for (int j = 0; j < A->columns; j++) {
      if (i == a || j == b) continue;
      tmp->matrix[k][l] = A->matrix[i][j];
      l++;
    }
    if (l != 0) k++;
  }
}
int s21_determinant(matrix_t *A, double *result) {
  int res = OK;
  if (incorrect(A)) {
    res = INCORRECT;
  } else if (A->rows != A->columns) {
    res = UNABLE;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];

    } else {
      double det = 0.;
      double one = 1.;
      for (int i = 0; i < A->columns; i++) {
        double tmp_det = 0.;
        if (A->rows > 1.) {
          matrix_t tmp;
          s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
          s21_copy_without_lines(0, i, A, &tmp);
          s21_determinant(&tmp, &tmp_det);
          det += A->matrix[0][i] * one * tmp_det;
          s21_remove_matrix(&tmp);
        } else {
          det += A->matrix[0][i] * one;
        }
        one = -one;
      }
      *result = det;
    }
  }
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (incorrect(A)) {
    res = INCORRECT;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    matrix_t tmp;
    double det;

    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
        s21_copy_without_lines(i, j, A, &tmp);
        s21_determinant(&tmp, &det);
        result->matrix[i][j] = det * pow(-1, i + j);
        s21_remove_matrix(&tmp);
      }
    }
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = OK;
  double det = 0;
  s21_determinant(A, &det);
  if (incorrect(A)) {
    res = INCORRECT;
  } else if (det == 0) {
    res = UNABLE;
  } else {
    s21_calc_complements(A, result);
    matrix_t tmp;
    s21_transpose(result, &tmp);
    s21_remove_matrix(result);
    s21_mult_number(&tmp, 1 / det, result);
    s21_remove_matrix(&tmp);
  }
  return res;
}

void Identity(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    for(int j = 0; j < A->columns; j++) {
      A->matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
}

void LU (matrix_t* A, matrix_t* L, matrix_t* U) {
  if (A->rows != A->columns) {
    printf("\nMATRIX IS NOT SQUARE\n");
    return;
  }
  s21_create_matrix(A->rows, A->columns, L);
  s21_create_matrix(A->rows, A->columns, U);
  s21_copy_without_lines(-1, -1, A, U);
  Identity(L);

  int n = A->rows;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      L->matrix[j][i] = U->matrix[j][i] / U->matrix[i][i];
    }
  }

  for (int k = 1; k < n; k++) {
    for (int i = k - 1; i < n; i++) {
      for (int j = i; j < n; j++) {
        L->matrix[j][i] = U->matrix[j][i] / U->matrix[i][i];
      }
    }
    for (int i = k; i < n; i++) {
      for (int j = k - 1; j < n; j++) {
        U->matrix[i][j] = U->matrix[i][j] - L->matrix[i][k-1] * U->matrix[k - 1][j];
      }
    }
  }
}

void Print_mat(matrix_t* A) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (j != 0) printf(" ");
      printf("%lf", A->matrix[i][j]);
    }
    printf("\n");
  }
}

void Solve(matrix_t *A, matrix_t *B, matrix_t *X) {
  matrix_t tmp;
  s21_create_matrix(B->rows, B->columns, &tmp);
  s21_create_matrix(B->rows, B->columns, X);
  int n = A->columns;

  for (int i = 0; i < n; i++) {
    tmp.matrix[0][i] = 0;
    X->matrix[0][i] = 0;
  }
  matrix_t L, U;
  LU(A, &L, &U);

  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < i; ++j) {
      if (i == j) continue;
      s += tmp.matrix[0][j] * L.matrix[i][j];
    }
    tmp.matrix[0][i] = B->matrix[0][i] - s;
  }
  
  for (int i = n - 1; i >= 0; --i) {
    double s = tmp.matrix[0][i];
    for (int j = n - 1; j >= i; --j) {
      if (i == j) continue;
      s -=   U.matrix[i][j] * X->matrix[0][j];
    }
    X->matrix[0][i] = s / U.matrix[i][i];
  }

  s21_remove_matrix(&tmp);
  s21_remove_matrix(&L);
  s21_remove_matrix(&U);
}

double det(matrix_t *A) {
  matrix_t L, U;
  double det = 1;
  LU(A, &L, &U);
  for (int i = 0; i < A->rows; ++i) {
    det *= U.matrix[i][i];
    det *= L.matrix[i][i];
  }
  s21_remove_matrix(&L);
  s21_remove_matrix(&U);
  return det;
}


void Inverse(matrix_t *A, matrix_t *res) {
  matrix_t inv, I, tmp, L, U;
  int n = A->rows;
  s21_create_matrix(A->rows, A->columns, &tmp);
  s21_create_matrix(A->rows, A->columns, &inv);
  s21_create_matrix(A->rows, A->columns, &I);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      I.matrix[i][j] = 0;
      tmp.matrix[i][j] = 0;
      inv.matrix[i][j] = 0;
      if (i == j) I.matrix[i][j] = 1;
    }
  }
  LU(A, &L, &U);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double s = 0;
      for (int k = 0; k < j; ++k) {
        s += tmp.matrix[k][i] * L.matrix[j][k];
      }
      tmp.matrix[j][i] = (I.matrix[j][i] - s) / L.matrix[j][j];
    }
  } 
  for (int i = 0; i < n; ++i) {
    for (int j = n - 1; j >= 0; --j) {
      double s = 0;
      for (int k = n - 1; k >= j; --k) {
        s += inv.matrix[k][i] * U.matrix[j][k];
      }
      inv.matrix[j][i] = (tmp.matrix[j][i] - s) / U.matrix[j][j];
    }
  }
  s21_create_matrix(inv.rows, inv.columns, res);
  s21_copy_without_lines(-1, -1, &inv, res);

  s21_remove_matrix(&tmp);
  s21_remove_matrix(&I);
  s21_remove_matrix(&inv);
  s21_remove_matrix(&L);
  s21_remove_matrix(&U);
}

void QP(matrix_t *A, matrix_t *B, matrix_t *X) {
  matrix_t Q, P;
  int n = A->rows;
  s21_create_matrix(1, n, &Q);
  s21_create_matrix(1, n, &P);
  s21_create_matrix(1, n, X);
  for (int i = 0; i < n; ++i) {
    Q.matrix[0][i] = 0;
    P.matrix[0][i] = 0;
    X->matrix[0][i] = 0;
  }
  for (int i = 0; i < n; ++i) {
    double tmp_q = 0, tmp_p = 0;
    if (i) {
      tmp_q = Q.matrix[0][i - 1];
      tmp_p = P.matrix[0][i - 1];
    }
    P.matrix[0][i] = -1. * A->matrix[i][2] / (A->matrix[i][0] * tmp_p + A->matrix[i][1]);
    Q.matrix[0][i] = (B->matrix[0][i] - A->matrix[i][0] * tmp_q) / (A->matrix[i][0] * tmp_p + A->matrix[i][1]);
  }

  X->matrix[0][n - 1] = Q.matrix[0][n - 1];

  for (int i = n - 2; i >= 0; --i) {
    X->matrix[0][i] = P.matrix[0][i] * X->matrix[0][i + 1] + Q.matrix[0][i];
  }

  s21_remove_matrix(&Q);
  s21_remove_matrix(&P);
}

void noralize(matrix_t *A, matrix_t *B, matrix_t *a, matrix_t *b) {
  for(int i = 0; i < A->rows; i++) {
    b->matrix[i][0] = B->matrix[i][0] / A->matrix[i][i];
    for (int j = 0; j < A->columns; j++) {
      a->matrix[i][j] = (i == j) ? 0 : - A->matrix[i][j] / A->matrix[i][i];
    }
  }
}

double norm(matrix_t *A) {
  double s = 0;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      s = pow(A->matrix[i][j], 2);
    }
  }
  return sqrt(s);
}

void iteration(matrix_t *A, matrix_t *B, double eps, matrix_t *res, int *iter) {
  if (det(A) == 0) {
    printf("DETERMINANT A == 0\n");
    return;
  }
  matrix_t X, a, b, tmp;

  s21_create_matrix(A->rows, A->columns, &a);
  s21_create_matrix(A->rows, 1, &b);
  s21_create_matrix(A->rows, 1, &X);
  s21_create_matrix(A->rows, 1, res);

  noralize(A, B, &a, &b);
  s21_copy_without_lines(-1, -1, &b, &X);
  *iter = 0;
  for (double error = 10. * eps; error > eps; (*iter)++) {
    s21_remove_matrix(res);
    s21_mult_matrix(&a, &X, &tmp);
    s21_sum_matrix(&tmp, &b, res);
    s21_remove_matrix(&tmp);

    error = 0;
    s21_sub_matrix(res, &X, &tmp);
    error = norm(&tmp);
    s21_remove_matrix(&tmp);

    s21_copy_without_lines(-1, -1, res, &X);
  }

  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&X);
}

void noralize_seidel(matrix_t *a, matrix_t *BB, matrix_t *CC) {
  for (int i = 0; i < a->rows; ++i) {
    for (int j = 0; j < a->columns; ++j) {
      BB->matrix[i][j] = (i > j) ? a->matrix[i][j] : 0. ;
      CC->matrix[i][j] = (j >= i) ? a->matrix[i][j] : 0. ;
    }
  }
}


void seidel(matrix_t *A, matrix_t *B, double eps, matrix_t *res, int *iter) {
  if (det(A) == 0) {
    printf("DETERMINANT A == 0\n");
    return;
  }
  matrix_t X, a, b, BB, CC, I, tmp, E_B_inv, tmp2, tmp3;

  s21_create_matrix(A->rows, A->columns, &a);
  s21_create_matrix(A->rows, A->columns, &BB);
  s21_create_matrix(A->rows, A->columns, &CC);
  s21_create_matrix(A->rows, A->columns, &I);
  s21_create_matrix(A->rows, 1, &b);
  s21_create_matrix(A->rows, 1, &X);
  s21_create_matrix(A->rows, 1, res);

  noralize(A, B, &a, &b);
  Identity(&I);
  noralize_seidel(&a, &BB, &CC);
  s21_copy_without_lines(-1, -1, &b, &X);
  *iter = 0;

  s21_sub_matrix(&I, &BB, &tmp);
  Inverse(&tmp, &E_B_inv);
  s21_remove_matrix(&tmp);
  s21_remove_matrix(&I);
  s21_mult_matrix(&E_B_inv, &CC, &tmp);
  s21_mult_matrix(&E_B_inv, &b, &tmp2);

  s21_remove_matrix(&CC);
  s21_remove_matrix(&BB);
  
  for (double error = 10. * eps; error > eps; (*iter)++) {
    s21_remove_matrix(res);
    s21_remove_matrix(&E_B_inv);

    s21_mult_matrix(&tmp, &X, &E_B_inv);
    s21_sum_matrix(&E_B_inv, &tmp2, res);


    error = 0;
    s21_sub_matrix(res, &X, &tmp3);
    error =  norm(&tmp3);
    s21_remove_matrix(&tmp3);
    s21_copy_without_lines(-1, -1, res, &X);
  }
  s21_remove_matrix(&E_B_inv);
  s21_remove_matrix(&tmp);
  s21_remove_matrix(&tmp2);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&X);
}

int symmetry_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; ++i) {
    for (int j = 0; j < A->columns; j++) {
      if (A->matrix[i][j] != A->matrix[j][i]) return FAILURE;
    }
  }
  return SUCCESS;
}
double sum_end(matrix_t *A) {
  double s = 0.;
  for (int i = 0; i < A->rows; i++) {
    for (int j = i+1; j < A->columns; j++) {
      s += pow(A->matrix[i][j], 2.);
    }
  }
  return sqrt(s);
}
void get_index_max_off_diagonal(matrix_t *A, int *i, int *j) {
  for (int k = 0; k < A->rows; ++k) {
    for (int m = k+1; m < A->columns; ++m) {
      if (k != m && fabs(A->matrix[k][m]) > fabs(A->matrix[*i][*j])) {
        *i = k;
        *j = m;
      }
    }
  }
}




void rotate(matrix_t *A, matrix_t *eigenvalue_vec, matrix_t *eigenvalue, double eps, int *iter) {
  if (A->rows != A->columns) {
    printf("matrix is not square\n");
    return;
  }
  
  if (!symmetry_matrix(A)) {
    printf("matrix is not symmetric\n");
    return;
  }
  
  matrix_t U, U_T, A_tmp, tmp;
  double phi = 0;
  
  
  s21_create_matrix(A->rows, A->columns, &A_tmp);
  s21_create_matrix(A->rows, A->columns, eigenvalue_vec);
  Identity(eigenvalue_vec);
  
  s21_copy_without_lines(-1, -1, A, &A_tmp);
  s21_create_matrix(A->rows, A->columns, &U);
  for (;sum_end(&A_tmp) > eps;  (*iter)++) {
    
    int i = 0, j = 1;
    get_index_max_off_diagonal(&A_tmp, &i, &j);
    
    if (fabs(A_tmp.matrix[i][i] - A_tmp.matrix[j][j]) < 1e-8) {
      phi = M_PI/4;
    } else {
      phi = 0.5 * atan(2.0 * A_tmp.matrix[i][j] / (A_tmp.matrix[j][j] - A_tmp.matrix[i][i]));
    }
    Identity(&U);

    U.matrix[i][i] = cos(phi);
    U.matrix[i][j] = -sin(phi);
    U.matrix[j][i] = sin(phi);
    U.matrix[j][j] = cos(phi);
    
    s21_transpose(&U, &U_T);
    s21_mult_matrix(&U_T, &A_tmp, &tmp);
    s21_remove_matrix(&A_tmp);
    s21_mult_matrix(&tmp, &U, &A_tmp);
    s21_remove_matrix(&tmp);
    s21_remove_matrix(&U_T);

    s21_mult_matrix(eigenvalue_vec, &U, &tmp);
    s21_copy_without_lines(-1, -1, &tmp, eigenvalue_vec);
    s21_remove_matrix(&tmp);
  
  }
  s21_remove_matrix(&U);
  s21_create_matrix(1, A->rows, eigenvalue);
  for (int i = 0; i < A->rows; i++) {
    eigenvalue->matrix[0][i] =  A_tmp.matrix[i][i];
  }

  s21_transpose(eigenvalue_vec, &tmp);
  s21_copy_without_lines(-1, -1, &tmp, eigenvalue_vec);
  s21_remove_matrix(&tmp);
  s21_remove_matrix(&A_tmp);
}


double sign(double t) {
  return (double)(t > 0);
}
double norm_col(matrix_t *A, int j) {
  double s = 0;
  for (int i = j; i < A->rows; i++) {
    s += pow(A->matrix[i][j], 2);
  }
  return sqrt(s);
}
void vector_fill_zeros_till_index(matrix_t *V, int i) {
  for (int j = 0; j < i; ++j) V->matrix[j][0] = 0.;
}
void householder(matrix_t* A, matrix_t* Q, matrix_t* R) {
  matrix_t A_tmp, V, V_T, V_res, ONE, H, I;
  if (A->rows != A->columns) {
    printf("matrix is not square\n");
    return;
  }

  s21_create_matrix(A->rows, A->columns, &A_tmp);
  s21_create_matrix(A->rows, A->columns, Q);
  s21_create_matrix(A->rows, A->columns, R);
  Identity(Q);
  s21_create_matrix(A->rows, A->columns, &I);
  Identity(&I);
  s21_copy_without_lines(-1, -1, A, &A_tmp);
  for (int i = 0; i < A->rows; i++) {
    s21_create_matrix(A->rows, 1, &V);
    vector_fill_zeros_till_index(&V, i);
    V.matrix[i][0] = A_tmp.matrix[i][i] + sign(A_tmp.matrix[i][i]) * norm_col(&A_tmp, i);
    for (int j = i + 1; j < A->rows; j++) {
      V.matrix[j][0] = A_tmp.matrix[j][i];
    }


    s21_transpose(&V, &V_T);
    s21_mult_matrix(&V, &V_T, &V_res);
    s21_mult_matrix(&V_T, &V, &ONE);
    s21_remove_matrix(&V_T);
    s21_remove_matrix(&V);

    double num = 2. / ONE.matrix[0][0];
    s21_remove_matrix(&ONE);


    s21_mult_number(&V_res, num, &V_T);
    s21_remove_matrix(&V_res);
    s21_sub_matrix(&I, &V_T, &H);
    s21_remove_matrix(&V_T);
    s21_mult_matrix(&H, &A_tmp, &V);
    s21_copy_without_lines(-1, -1, &V, &A_tmp);
    s21_remove_matrix(&V);


    s21_mult_matrix(Q, &H, &ONE);
    s21_copy_without_lines(-1, -1, &ONE, Q);
    s21_copy_without_lines(-1, -1, &A_tmp, R);
    s21_remove_matrix(&H);
    s21_remove_matrix(&ONE);
  }
  s21_remove_matrix(&I);
  s21_remove_matrix(&A_tmp);
}

double off_diagonal_norm(matrix_t *A) {
  double s = 0;
  for (int i = 1; i < A->rows; i++) {
    s+= pow(A->matrix[i][0], 2);
  }
  return sqrt(s);
}
void eigenvalue(matrix_t *A, matrix_t *eigenvalue, double eps, int *iter) {
  if (A->rows != A->columns) {
    printf("matrix is not square\n");
    return;
  }
  matrix_t Q, R, A_tmp;
  s21_create_matrix(A->rows, A->columns, &A_tmp);
  s21_copy_without_lines(-1, -1, A, &A_tmp);
  double error = eps * 10.;
  *iter = 0;
  for(; error > eps; (*iter)++) {
    householder(&A_tmp, &Q, &R);
    s21_remove_matrix(&A_tmp);
    s21_mult_matrix(&R, &Q, &A_tmp);
    error = off_diagonal_norm(&A_tmp);
    s21_remove_matrix(&Q);
    s21_remove_matrix(&R);
  }
  s21_create_matrix(A->rows, 1, eigenvalue);
  for (int i = 0; i < A->rows; i++) {
    eigenvalue->matrix[i][0] = A_tmp.matrix[i][i];
  }
  s21_remove_matrix(&A_tmp);
}

