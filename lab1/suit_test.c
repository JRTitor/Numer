#include <check.h>

#include "s21_matrix.h"

START_TEST(add_test) {
  matrix_t a, b, c, d;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &b);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  b.matrix[0][0] = 3;
  b.matrix[0][1] = 2;
  b.matrix[0][2] = 1;
  b.matrix[1][0] = 3;
  b.matrix[1][1] = 2;
  b.matrix[1][2] = 1;
  b.matrix[2][0] = 3;
  b.matrix[2][1] = 2;
  b.matrix[2][2] = 1;
  d.matrix[0][0] = 4;
  d.matrix[0][1] = 4;
  d.matrix[0][2] = 4;
  d.matrix[1][0] = 4;
  d.matrix[1][1] = 4;
  d.matrix[1][2] = 4;
  d.matrix[2][0] = 4;
  d.matrix[2][1] = 4;
  d.matrix[2][2] = 4;
  int res = s21_sum_matrix(&a, &b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), SUCCESS);
  else
    ck_assert_int_eq(res, OK);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(sub_test) {
  matrix_t a, b, c, d;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &b);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  b.matrix[0][0] = 3;
  b.matrix[0][1] = 2;
  b.matrix[0][2] = 1;
  b.matrix[1][0] = 3;
  b.matrix[1][1] = 2;
  b.matrix[1][2] = 1;
  b.matrix[2][0] = 3;
  b.matrix[2][1] = 2;
  b.matrix[2][2] = 1;
  d.matrix[0][0] = -2;
  d.matrix[0][1] = 0;
  d.matrix[0][2] = 2;
  d.matrix[1][0] = -2;
  d.matrix[1][1] = 0;
  d.matrix[1][2] = 2;
  d.matrix[2][0] = -2;
  d.matrix[2][1] = 0;
  d.matrix[2][2] = 2;
  int res = s21_sub_matrix(&a, &b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), SUCCESS);
  else
    ck_assert_int_eq(res, OK);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(mul_num) {
  matrix_t a, c, d;
  double b = 2;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  d.matrix[0][0] = 2;
  d.matrix[0][1] = 4;
  d.matrix[0][2] = 6;
  d.matrix[1][0] = 2;
  d.matrix[1][1] = 4;
  d.matrix[1][2] = 6;
  d.matrix[2][0] = 2;
  d.matrix[2][1] = 4;
  d.matrix[2][2] = 6;
  int res = s21_mult_number(&a, b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), SUCCESS);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(mul_matrix) {
  matrix_t a, b, c, d;
  s21_create_matrix(2, 3, &a);
  s21_create_matrix(3, 4, &b);
  s21_create_matrix(2, 4, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  b.matrix[0][0] = 3;
  b.matrix[0][1] = 2;
  b.matrix[0][2] = 1;
  b.matrix[0][3] = 1;
  b.matrix[1][0] = 3;
  b.matrix[1][1] = 2;
  b.matrix[1][2] = 1;
  b.matrix[1][3] = 1;
  b.matrix[2][0] = 3;
  b.matrix[2][1] = 2;
  b.matrix[2][2] = 1;
  b.matrix[2][3] = 1;
  d.matrix[0][0] = 18;
  d.matrix[0][1] = 12;
  d.matrix[0][2] = 6;
  d.matrix[0][3] = 6;
  d.matrix[1][0] = 18;
  d.matrix[1][1] = 12;
  d.matrix[1][2] = 6;
  d.matrix[1][3] = 6;
  int res = s21_mult_matrix(&a, &b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), SUCCESS);
  else
    ck_assert_int_eq(res, OK);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(trans_test) {
  matrix_t a, c, d;
  s21_create_matrix(2, 3, &a);
  s21_create_matrix(3, 2, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  d.matrix[0][0] = 1;
  d.matrix[0][1] = 1;
  d.matrix[1][0] = 2;
  d.matrix[1][1] = 2;
  d.matrix[2][0] = 3;
  d.matrix[2][1] = 3;
  int res = s21_transpose(&a, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), SUCCESS);
  else
    ck_assert_int_eq(res, OK);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(inv_test) {
  matrix_t a, c, d;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 5;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 7;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  d.matrix[0][0] = -1.0 / 6;
  d.matrix[0][1] = 0;
  d.matrix[0][2] = 1.0 / 6;
  d.matrix[1][0] = -1.0 / 3;
  d.matrix[1][1] = 1.0 / 3;
  d.matrix[1][2] = 0;
  d.matrix[2][0] = 11.0 / 18;
  d.matrix[2][1] = -2.0 / 9;
  d.matrix[2][2] = -1.0 / 18;
  int res = s21_inverse_matrix(&a, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), SUCCESS);
  else
    ck_assert_int_eq(res, OK);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(eq_test) {
  matrix_t a, c;
  s21_create_matrix(2, 2, &a);
  s21_create_matrix(2, 2, &c);
  a.matrix[0][0] = NAN;
  a.matrix[0][1] = 0;
  a.matrix[1][0] = 0;
  a.matrix[1][1] = 0;
  c.matrix[0][0] = NAN;
  c.matrix[0][1] = 1;
  c.matrix[1][0] = 0;
  c.matrix[1][1] = 0;
  int res = s21_eq_matrix(&a, &c);
  ck_assert_int_eq(res, FAILURE);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(det_test_0) {
  matrix_t a;
  double b;
  s21_create_matrix(5, 5, &a);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[0][3] = 4;
  a.matrix[0][4] = 5;
  a.matrix[1][0] = 6;
  a.matrix[1][1] = 7;
  a.matrix[1][2] = 8;
  a.matrix[1][3] = 9;
  a.matrix[1][4] = -10;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  a.matrix[2][3] = 4;
  a.matrix[2][4] = 5;
  a.matrix[3][0] = 6;
  a.matrix[3][1] = 7;
  a.matrix[3][2] = 8;
  a.matrix[3][3] = 9;
  a.matrix[3][4] = -10;
  a.matrix[4][0] = 1;
  a.matrix[4][1] = 2;
  a.matrix[4][2] = 3;
  a.matrix[4][3] = 4;
  a.matrix[4][4] = 5;
  s21_determinant(&a, &b);
  ck_assert_double_eq_tol(b, 0, 1e-6);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(det_test_1) {
  matrix_t a;
  double b;
  s21_create_matrix(5, 5, &a);
  a.matrix[0][0] = 55;
  a.matrix[0][1] = 54;
  a.matrix[0][2] = 53;
  a.matrix[0][3] = -52;
  a.matrix[0][4] = 51;
  a.matrix[1][0] = 50;
  a.matrix[1][1] = 49;
  a.matrix[1][2] = 48;
  a.matrix[1][3] = 47;
  a.matrix[1][4] = -10;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  a.matrix[2][3] = 4;
  a.matrix[2][4] = -7;
  a.matrix[3][0] = 6;
  a.matrix[3][1] = 7;
  a.matrix[3][2] = 8;
  a.matrix[3][3] = 9;
  a.matrix[3][4] = 11;
  a.matrix[4][0] = 13;
  a.matrix[4][1] = 17;
  a.matrix[4][2] = 19;
  a.matrix[4][3] = 11;
  a.matrix[4][4] = 7;
  s21_determinant(&a, &b);
  ck_assert_double_eq_tol(b, 208624, 1e-6);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(inverse_matrix_1) {
  matrix_t A, B, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  s21_inverse_matrix(&A, &B);
  C.matrix[0][0] = 1;
  C.matrix[0][1] = -1;
  C.matrix[0][2] = 1;
  C.matrix[1][0] = -38;
  C.matrix[1][1] = 41;
  C.matrix[1][2] = -34;
  C.matrix[2][0] = 27;
  C.matrix[2][1] = -29;
  C.matrix[2][2] = 24;
  ck_assert_int_eq(s21_eq_matrix(&B, &C), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(inverse_matrix_2) {
  matrix_t A, B, C;
  s21_create_matrix(4, 4, &A);
  s21_create_matrix(4, 4, &C);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = 4;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[1][3] = 1;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  A.matrix[2][3] = 3;
  A.matrix[3][0] = 3;
  A.matrix[3][1] = 7;
  A.matrix[3][2] = 2;
  A.matrix[3][3] = 4;
  s21_inverse_matrix(&A, &B);
  C.matrix[0][0] = -0.125;
  C.matrix[0][1] = 0.19166666;
  C.matrix[0][2] = 0.325;
  C.matrix[0][3] = -0.16666666;
  C.matrix[1][0] = -0.125;
  C.matrix[1][1] = 0.05833333;
  C.matrix[1][2] = -0.075;
  C.matrix[1][3] = 0.16666666;
  C.matrix[2][0] = 0.125;
  C.matrix[2][1] = 0.675;
  C.matrix[2][2] = 0.275;
  C.matrix[2][3] = -0.5;
  C.matrix[3][0] = 0.25;
  C.matrix[3][1] = -0.58333333;
  C.matrix[3][2] = -0.25;
  C.matrix[3][3] = 0.33333333;
  ck_assert_int_eq(s21_eq_matrix(&B, &C), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(inverse_matrix_3) {
  matrix_t A, B;
  s21_create_matrix(3, 4, &A);
  ck_assert_int_eq(s21_inverse_matrix(&A, &B), UNABLE);
  s21_remove_matrix(&A);
}
END_TEST
START_TEST(incorrect_add) {
  matrix_t A, B, C;
  s21_create_matrix(3, 4, &A);
  s21_create_matrix(4, 3, &B);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &C), UNABLE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST
START_TEST(incorrect_sub) {
  matrix_t A, B, C;
  s21_create_matrix(7, 4, &A);
  s21_create_matrix(4, 9, &B);
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &C), UNABLE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST
START_TEST(correct_mul) {
  matrix_t A, B, C;
  s21_create_matrix(7, 4, &A);
  s21_create_matrix(4, 9, &B);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &C), OK);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST
START_TEST(incorrect_mul) {
  matrix_t A, B, C;
  s21_create_matrix(7, 4, &A);
  s21_create_matrix(11, 9, &B);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &C), UNABLE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST
START_TEST(incorrect_det) {
  matrix_t A;
  double d;
  s21_create_matrix(7, 4, &A);
  ck_assert_int_eq(s21_determinant(&A, &d), UNABLE);
  s21_remove_matrix(&A);
}
END_TEST
START_TEST(incorrect_inverse) {
  matrix_t A, B;
  s21_create_matrix(7, 4, &A);
  ck_assert_int_eq(s21_inverse_matrix(&A, &B), UNABLE);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(zero_det) {
  matrix_t A;
  double d = 0;
  s21_create_matrix(3, 3, &A);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A.matrix[i][j] = 0.;
    }
  }
  ck_assert_int_eq(s21_determinant(&A, &d), OK);
  ck_assert_int_eq(d, 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(random_inverse) {
  matrix_t A, B, C, D;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &D);

  A.matrix[0][0] = 21.0;
  A.matrix[0][1] = 99.3;
  A.matrix[0][2] = 11.78;
  A.matrix[1][0] = -122.1;
  A.matrix[1][1] = 0.93;
  A.matrix[1][2] = 37.1;
  A.matrix[2][0] = -17.9;
  A.matrix[2][1] = 21.28;
  A.matrix[2][2] = -67.34;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j)
        D.matrix[i][j] = 1;
      else
        D.matrix[i][j] = 0;
    }
  }
  ck_assert_int_eq(s21_inverse_matrix(&A, &B), OK);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &C), OK);
  ck_assert_int_eq(s21_eq_matrix(&D, &C), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
  s21_remove_matrix(&D);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("s21_matrix");
  TCase *tc1_1 = tcase_create("all");
  SRunner *sr = srunner_create(s1);
  int nf;

  suite_add_tcase(s1, tc1_1);
  tcase_add_test(tc1_1, add_test);
  tcase_add_test(tc1_1, sub_test);
  tcase_add_test(tc1_1, mul_num);
  tcase_add_test(tc1_1, mul_matrix);
  tcase_add_test(tc1_1, trans_test);
  tcase_add_test(tc1_1, inv_test);
  tcase_add_test(tc1_1, eq_test);
  tcase_add_test(tc1_1, det_test_0);
  tcase_add_test(tc1_1, det_test_1);
  tcase_add_test(tc1_1, inverse_matrix_1);
  tcase_add_test(tc1_1, inverse_matrix_2);
  tcase_add_test(tc1_1, inverse_matrix_3);
  tcase_add_test(tc1_1, incorrect_add);
  tcase_add_test(tc1_1, incorrect_sub);
  tcase_add_test(tc1_1, incorrect_mul);
  tcase_add_test(tc1_1, correct_mul);
  tcase_add_test(tc1_1, incorrect_det);
  tcase_add_test(tc1_1, incorrect_inverse);
  tcase_add_test(tc1_1, zero_det);
  tcase_add_test(tc1_1, random_inverse);

  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}
