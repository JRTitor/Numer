#include"labs.h"

void load_lab(const int l) {
    if (l == 1) lab1();
    if (l == 2) lab2();
    if (l == 3) lab3();
    if (l == 4) lab4();
    if (l == 5) lab5();
}

int defaul_or_input() {
    int stop = 1;
    int choice = 0;
    while (stop) {
        printf("Введите 1 для демонстации лаботаторной работы с дефолтными параметрами \nВведите 2 для ввода данных\n");
        if (scanf("%d", &choice) == 1 && (choice == 1 || choice == 2)) {
            break;
        } else {
            printf("Вам требуется ввести 1 или 2 \n");
        }
        
    }
    return choice;
}

void lab1() {
    int choice = 0;
    choice = defaul_or_input();

    if (choice == 1) {
        lab1_default();
    }
    if (choice == 2) {
        lab1_input();
    }

}



void get_row_col(matrix_t *A) {
    int row = 0;
    int col = 0;
    printf("Введите размерность киличество строк и столбцов через пробел  (минимальная размерность: 1  1)\n");
    while (scanf("%d %d", &row, &col) != 2) {
        if (row  > 0 && col > 0) {
            break;
        }
        printf("Некоррктная размерность\n");
    }
    s21_create_matrix(row, col, A);
}

void fill_values(matrix_t *A) {
    double tmp = 0;
    fill_zero(A);
    printf("Введите матрицу\n");
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
            Print_mat(A);
            printf("Введите элемент %d %d", i+1, j+1);
            while (scanf("%lf", &tmp) != 1) {
                printf("Введите корректно число  \n ");
            }
            A->matrix[i][j] = tmp;

        }
    }
}

void Constractor(matrix_t *A) {
    get_row_col(A);
    fill_values(A);
}

void lab1_input() {
    matrix_t A, B, X;
    int stop = 1;

    while (stop) {
        printf("Введите матрицу A\n");
        Constractor(&A);
        printf("Введите транспонированную B (строка вместо столбца)\n");
        Constractor(&B);

        if (A.rows == B.columns) {
            break;
        }
        printf("Ввод некорректен\n");
        s21_remove_matrix(&A);
        s21_remove_matrix(&B);
    }
    
    printf("A:\n");
    Print_mat(&A);
    printf("B^T:\n");
    Print_mat(&B);

    printf("det(A) = %lf\n",det(&A, 0));
    
    printf("X^T = \n");
    Solve(&A, &B, &X, 0);
    Print_mat(&X);
    s21_remove_matrix(&X);

    printf("A^(-1):\n");
    Inverse(&A, &X, 0);
    Print_mat(&X);
    s21_remove_matrix(&X);

    s21_remove_matrix(&B);
    s21_remove_matrix(&A);
}

void lab1_default() {
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

    printf("det(A) = %lf\n",det(&A, 0));
    
    printf("X^T = \n");
    Solve(&A, &B, &X, 0);
    Print_mat(&X);
    s21_remove_matrix(&X);

    printf("A^(-1):\n");
    Inverse(&A, &X, 0);
    Print_mat(&X);
    s21_remove_matrix(&X);

    s21_remove_matrix(&B);
    s21_remove_matrix(&A);
}

void lab2() {
    int choice = 0;
    choice = defaul_or_input();

    if (choice == 1) {
        lab2_default();
    }
    if (choice == 2) {
        lab2_input();
    }
}

void lab2_default() {
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

void lab2_input() {
    matrix_t A, B, X;

    int stop = 1;

    while (stop) {
        printf("Введите матрицу A в трёхдиагональном виде (её размерность может быть, например 5 3)\n");
        Constractor(&A);
        printf("Введите транспонированную B (строка вместо столбца)\n");
        Constractor(&B);

        if (A.rows == B.columns) {
            break;
        }
        printf("Ввод некорректен\n");
        s21_remove_matrix(&A);
        s21_remove_matrix(&B);
    }


    printf("\nA в рёхдмагональном виде = \n");
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

void lab3() {
    int choice = 0;
    choice = defaul_or_input();

    if (choice == 1) {
        lab3_default();
    }
    if (choice == 2) {
        lab3_input();
    }
}

void lab3_default() {
    matrix_t A, B, X;
    double eps = 0.01;
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
    printf("eps = %lf\n", eps);

    printf("Метод инетраций:\n");
    iteration(&A, &B, eps, &X, &iter);
    printf("количество инетраций: %d\n", iter);
    printf("X = \n");
    Print_mat(&X);
    s21_remove_matrix(&X);

    printf("Метод зейделя:\n");
    seidel(&A, &B, eps, &X, &iter);
    printf("количество инетраций: %d\n", iter);
    printf("X = \n");
    Print_mat(&X);
    s21_remove_matrix(&X);

    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
}

void lab3_input() {
    matrix_t A, B, X;
    double eps = 0.01;
    int iter = 0;

    while (1) {
        printf("Введите eps");
        while (scanf("%lf", &eps) != 1 || eps < 0) {
            printf("eps не некорректен, возможно отрицателен попробуйте снова, желательно сделать его меньше 1");
        }
        
        printf("Введите матрицу A\n");
        Constractor(&A);
        printf("Введите B (столбец)\n");
        Constractor(&B);

        if (A.rows == B.columns) {
            break;
        }
        printf("Ввод некорректен\n");
        s21_remove_matrix(&A);
        s21_remove_matrix(&B);
    }

    printf("A:\n");
    Print_mat(&A);
    printf("B:\n");
    Print_mat(&B);
    printf("eps = %lf\n", eps);

    printf("Метод инетраций:\n");
    iteration(&A, &B, eps, &X, &iter);
    printf("количество инетраций: %d\n", iter);
    printf("X = \n");
    Print_mat(&X);
    s21_remove_matrix(&X);

    printf("Метод зейделя:\n");
    seidel(&A, &B, eps, &X, &iter);
    printf("количество инетраций: %d\n", iter);
    printf("X = \n");
    Print_mat(&X);
    s21_remove_matrix(&X);

    s21_remove_matrix(&A);
    s21_remove_matrix(&B);
}

void lab4() {
    int choice = 0;
    choice = defaul_or_input();

    if (choice == 1) {
        lab4_default();
    }
    if (choice == 2) {
        lab4_input();
    }

}

void lab4_default() {
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
    printf("eps = %lf\n", eps);
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

void lab4_input() {
    matrix_t A, euginvalues, euginvalue_vec;
    double eps = 0.1;
    int iter = 0;

    printf("Введите eps");
    while (scanf("%lf", &eps) != 1 || eps < 0) {
        printf("eps не некорректен, возможно отрицателен попробуйте снова, желательно сделать его меньше 1");
    }
        
    printf("Введите матрицу A\n");
    Constractor(&A);

    

    printf("Матрица A:\n");
    Print_mat(&A);
    printf("eps = %g\n", eps);
    rotate(&A, &euginvalue_vec, &euginvalues ,eps, &iter);
    printf("инерации: %d\n", iter);
    printf("собственные значения\n");
    Print_mat(&euginvalues);
    printf("векторы собственных значений\n");
    Print_mat(&euginvalue_vec);
    s21_remove_matrix(&A);
    s21_remove_matrix(&euginvalue_vec);
    s21_remove_matrix(&euginvalues);
}

void lab5() {
    int choice = 0;
    choice = defaul_or_input();

    if (choice == 1) {
        lab5_default();
    }
    if (choice == 2) {
        lab5_input();
    }
}

void lab5_input() {
    double eps = 0.01;
    matrix_t A, eigenvalues;
    int iter = 0;


    printf("Введите eps");
    while (scanf("%lf",&eps) != 1 || eps < 0) {
        printf("eps не некорректен, возможно отрицателен попробуйте снова, желательно сделать его меньше 1");
    }
        
    printf("Введите матрицу A\n");
    Constractor(&A);


    printf("матрица A:\n");
    Print_mat(&A);
    printf("eps = %g\n", eps);
    eigenvalue(&A, &eigenvalues, eps, &iter);
    printf("Итерации: %d\n", iter);
    printf("Собственные значения\n");
    Print_mat(&eigenvalues);

    s21_remove_matrix(&A);
    s21_remove_matrix(&eigenvalues);
}

void lab5_default() {
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

    printf("матрица A:\n");
    Print_mat(&A);
    printf("eps = %g\n", eps);
    eigenvalue(&A, &eigenvalues, eps, &iter);
    printf("Итерации: %d\n", iter);
    printf("Собственные значения\n");
    Print_mat(&eigenvalues);

    s21_remove_matrix(&A);
    s21_remove_matrix(&eigenvalues);
}


