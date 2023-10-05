CC = gcc 
FLAGS = -Werror -Wextra -Wall -std=c11 -pedantic -g
GCOVFLAGS =  --coverage
SUITE_CASES_C = test/suite_*.c
S21_MAT_C = s21_*.c
S21_MAT_O = s21_*.o
OS = $(shell uname -s)
LINUX = -lcheck -lm -lrt -lpthread -lsubunit

s21_matrix.o:  s21_matrix.h
		$(CC) $(FLAGS) -c $(S21_MAT_C)

s21_matrix.cov:  s21_matrix.h
		$(CC) $(FLAGS) -c $(S21_MAT_C) $(GCOVFLAGS)

s21_matrix.a: s21_matrix.o
		ar rcs s21_matrix.a $(S21_MAT_O)

s21_matrix.a_cov: s21_matrix.cov
		ar rcs s21_matrix.a $(S21_MAT_O)

menu: s21_matrix.a
	$(CC) -Werror -Wextra -Wall -std=c11 labs.c menu.c -lm -L. -l:s21_matrix.a -o interface
	./interface


test: s21_matrix.a_cov suit_test.c
ifeq ($(OS), Darwin)
		$(CC) -Werror -Wextra -Wall -std=c11    suit_test.c  s21_matrix.a -lcheck -o s21_test $(GCOVFLAGS)
else
		$(CC)   -Werror -Wextra -Wall -std=c11  suit_test.c  -L. -l:s21_matrix.a -o s21_test  $(GCOVFLAGS) $(LINUX)
endif
		./s21_test

gcov_report: test
		lcov -t "gcov_report" -o s21_gcov_report.info -c -d .
		genhtml -o s21_report s21_gcov_report.info
		open ./s21_report/index.html
lab1_1: clean s21_matrix.a
	@echo
	@echo
	@echo "lab1_1"
	gcc -Werror -Wextra -Wall -std=c11 lab1_1.c -lm -L. -l:s21_matrix.a -o 1.out
	./1.out

lab1_2: clean s21_matrix.a
	@echo
	@echo
	@echo "lab1_2"
	gcc -Werror -Wextra -Wall -std=c11 lab1_2.c -lm -L. -l:s21_matrix.a -o 2.out
	./2.out

lab1_3: clean s21_matrix.a
	@echo
	@echo
	@echo "lab1_3"
	gcc -Werror -Wextra -Wall -std=c11 lab1_3.c -lm -L. -l:s21_matrix.a -o 3.out
	./3.out

lab1_4: clean s21_matrix.a
	@echo
	@echo
	@echo "lab1_4"
	gcc -Werror -Wextra -Wall -std=c11 lab1_4.c -lm -L. -l:s21_matrix.a -o 4.out
	./4.out

lab1_5: clean s21_matrix.a
	@echo
	@echo
	@echo "lab1_5"
	gcc -Werror -Wextra -Wall -std=c11 lab1_5.c -lm -L. -l:s21_matrix.a -o 5.out
	./5.out


clean : 
		rm -rf s21_matrix.a *.o *.a *.gcno  *.gcda *.info report.out *.gcov  s21_test s21_report *.out RES* 

