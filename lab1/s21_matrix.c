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



void Identity(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    for(int j = 0; j < A->columns; j++) {
      A->matrix[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
}

void switch_lines(matrix_t *LU, int i, int s) {
  double tmp = 0;
  for (int j = 0; j < LU->columns; j++) {
    tmp = LU->matrix[i][j];
    LU->matrix[i][j] = LU->matrix[s][j];
    LU->matrix[s][j] = tmp;
  }
}
int check_LU(matrix_t* A) {
  int n = A->rows;
  for (int i = 0; i < n - 1; i++) {
    int tmp = i + 1;
    while (A->matrix[i][i] == 0) {
      if (tmp == n) {
        printf("LU is not supported for this matrix");
        return 0;
      } 
      switch_lines(A, i, tmp);
      tmp++;
    }
  }
  return 1;
}


void LU(matrix_t *A, matrix_t *Lower, int show_lu) {
  if (A->rows != A->columns) {
    printf("\nMATRIX IS NOT SQUARE\n");
    return;
  }
  if (!check_LU(A)) {
    _Exit(0);
  }

  int n = A->columns;
  double mu = 0;
  matrix_t M, M2, Upper, tmp_A2, tmp_LU;
  s21_create_matrix(n, n, Lower);
  s21_create_matrix(n, n, &Upper);

  s21_copy_without_lines(-1, -1, A, &Upper);
  Identity(Lower);
  s21_create_matrix(n, n, &M);
  s21_create_matrix(n, n, &M2);
  Identity(&M);
  Identity(&M2);

  for (int j = 0; j < n - 1; j++) {
    Identity(&M);
    Identity(&M2);
    for(int i = j + 1; i < n; i++) {
      mu = Upper.matrix[i][j] / Upper.matrix[j][j];
      M.matrix[i][j] = -mu;
      M2.matrix[i][j] = mu;
    }
    
    s21_mult_matrix(&M, &Upper, &tmp_A2);
    s21_copy_without_lines(-1, -1, &tmp_A2, &Upper);
    s21_remove_matrix(&tmp_A2);

    s21_mult_matrix(Lower, &M2, &tmp_LU);
    s21_copy_without_lines(-1, -1, &tmp_LU, Lower);
    s21_remove_matrix(&tmp_LU);
  }
  // Putting Lower and Upper together
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j ++) {
      Upper.matrix[i][j] = Lower->matrix[i][j];
    }
  }
  s21_copy_without_lines(-1, -1, &Upper, Lower);
  s21_remove_matrix(&Upper);
  if (show_lu == 1) {
    printf("LU разложение матрицы:\n");
    Print_mat(Lower);
  }
}

void Print_mat(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (j != 0) printf(" ");
      printf("%lf", A->matrix[i][j]);
    }
    printf("\n");
  }
}

void forward_sub(matrix_t *L, matrix_t *B, matrix_t *X) {
  int n = L->columns;
  double tmp = 0;
  s21_create_matrix(B->rows, B->columns, X);
  fill_zero(X);
  X->matrix[0][0] = B->matrix[0][0];
  for (int i = 1; i < n; i++) {
    tmp = B->matrix[0][i];
    for (int j = 0; j < i; j++) {
      tmp -= L->matrix[i][j] * X->matrix[0][j];
    }
    X->matrix[0][i] = tmp;
  }
}

void backward_sub(matrix_t *U, matrix_t *B, matrix_t *X) {
  int n = U->columns;
  double tmp = 0;
  s21_create_matrix(B->rows, B->columns, X);
  fill_zero(X);
  for (int i = n - 1; i > -1; i--) {
    tmp = B->matrix[0][i];
    for (int j = i + 1; j < n; j++) {
      tmp -= U->matrix[i][j] * X->matrix[0][j];
    }
    X->matrix[0][i] = tmp / U->matrix[i][i];
  }
}

void Solve(matrix_t *A, matrix_t *B, matrix_t *X, int show_lu) {
  matrix_t tmp, LU_m;
  LU(A, &LU_m, show_lu);
  forward_sub(&LU_m, B, &tmp);
  backward_sub(&LU_m, &tmp, X);
  s21_remove_matrix(&tmp);
  s21_remove_matrix(&LU_m);
}



double det(matrix_t *A, int show_lu) {
  matrix_t LU_m;
  double det = 1.;
  LU(A, &LU_m, show_lu);
  for (int i = 0; i < A->rows; i++) {
    det *= LU_m.matrix[i][i];
  }
  s21_remove_matrix(&LU_m);
  return det;
}


void Inverse(matrix_t *A, matrix_t *res, int show_lu) {
  if (A->rows != A->columns) {
    printf("\nMATRIX IS NOT SQUARE\n");
    return;
  }
  matrix_t inv, I, tmp, LU_m;
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
  LU(A, &LU_m, show_lu);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double s = 0;
      for (int k = 0; k < j; ++k) {
        s += tmp.matrix[k][i] * LU_m.matrix[j][k];
      }
      tmp.matrix[j][i] = (I.matrix[j][i] - s);
    }
  } 
  for (int i = 0; i < n; ++i) {
    for (int j = n - 1; j >= 0; --j) {
      double s = 0;
      for (int k = n - 1; k >= j; --k) {
        s += inv.matrix[k][i] * LU_m.matrix[j][k];
      }
      inv.matrix[j][i] = (tmp.matrix[j][i] - s) / LU_m.matrix[j][j];
    }
  }
  s21_create_matrix(inv.rows, inv.columns, res);
  s21_copy_without_lines(-1, -1, &inv, res);

  s21_remove_matrix(&tmp);
  s21_remove_matrix(&I);
  s21_remove_matrix(&inv);
  s21_remove_matrix(&LU_m);
}

void fill_zero(matrix_t *A) {
  for (int i = 0;  i < A->rows; i++) {
    for (int j = 0;  j < A->rows; j++) {
      A->matrix[i][j] = 0.;
    }
  }
}

// void Inverse(matrix_t *A, matrix_t *res, int show_lu) {
//   if (A->rows != A->columns) {
//     printf("\nMATRIX IS NOT SQUARE\n");
//     return;
//   }

//   matrix_t x1, res1;
//   s21_create_matrix(A->rows, A->columns, res);
//   s21_create_matrix(1, A->columns, &x1);
//   fill_zero(&x1);
//   for (int i = 0;  i < A->rows; i++) {
//     if (i > 0)
//       x1.matrix[0][i - 1] = 0.;
//     x1.matrix[0][i] = 1.;
//     Solve(A, &x1, &res1, show_lu);
//     for (int j = 0;  j < A->rows; j++) {
//       res->matrix[j][i] = res1.matrix[0][j];
//     }
//     s21_remove_matrix(&res1);
//   }
//   s21_remove_matrix(&x1);
// }

void is_tridiagonal(matrix_t *A) {

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (i == j || i == j - 1 || i == j + 1) {
        if (A->matrix[i][j] == 0) {
          printf("\nMatrix is not tridiagonal\n");
          return;
        }
      } else if (A->matrix[i][j] != 0) {
        printf("\nMatrix is not tridiagonal\n");
        return;
      }
    }
  }
}

void stable_qr(matrix_t *A) {
  int n = A->rows;
  int flag = 0;
  for (int i = 0; i < n - 1; i++) {
    if ((i > 0) && (i < (n - 1)) && (A->matrix[i][0] == 0 || A->matrix[i][1] == 0 || A->matrix[i][2] == 0)) {
      printf("\nUNSTABLE\n");
      return;
    }
    if (fabs(A->matrix[i][1]) < fabs(A->matrix[i][0]) + fabs(A->matrix[i][2])) {
      printf("\nUNSTABLE\n");
      return;
    } if (fabs(A->matrix[i][1]) > fabs(A->matrix[i][0]) + fabs(A->matrix[i][2])) { flag++;}
  }
  if (!flag) {
    printf("\nUNSTABLE\n");
  }
}

void QP(matrix_t *A, matrix_t *B, matrix_t *X) {
  stable_qr(A);
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

double norm(matrix_t *A ) {
  double s = 0;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      s = pow(A->matrix[i][j], 2);
    }
  }
  return sqrt(s);
}

void norm_m(matrix_t *A) {
  double m = -1.;
  for (int i = 0; i < A->rows; i++) {
    double tmp = 0;
    for (int j = 0; j < A->columns; j++) {
      tmp += fabs(A->matrix[i][j]);
      m = (m < tmp) ? tmp : m; 
    }
  }
  if (m < 1)  printf("\nmethod converges\n");
  else  printf("\nmethod dont converges\n");
}

void iteration(matrix_t *A, matrix_t *B, double eps, matrix_t *res, int *iter) {
  if (det(A, 0) == 0) {
    printf("DETERMINANT A == 0\n");
    return;
  }
  matrix_t X, a, b, tmp;

  s21_create_matrix(A->rows, A->columns, &a);
  s21_create_matrix(A->rows, 1, &b);
  s21_create_matrix(A->rows, 1, &X);
  s21_create_matrix(A->rows, 1, res);

  noralize(A, B, &a, &b);
  norm_m(&a);
  s21_copy_without_lines(-1, -1, &b, &X);
  *iter = 0;
  for (double error = 10. * eps; error > eps; (*iter)++) {
    s21_remove_matrix(res);
    s21_mult_matrix(&a, &X, &tmp);
    s21_sum_matrix(&tmp, &b, res);
    s21_remove_matrix(&tmp);


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
  if (det(A, 0) == 0) {
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
  norm_m(&a);
  Identity(&I);
  noralize_seidel(&a, &BB, &CC);
  s21_copy_without_lines(-1, -1, &b, &X);
  *iter = 0;

  s21_sub_matrix(&I, &BB, &tmp);
  Inverse(&tmp, &E_B_inv, 0);
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


