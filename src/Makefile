CC = gcc -lstdc++ -std=c++17 -g -Wall -Werror -Wextra
GCOVFLAGS = -fprofile-arcs -ftest-coverage
SOURSES = *.cc
OBJECT = $(patsubst %*.cc, %*.o,  ${SOURSES})
TEST_FLAGS =-lgtest -lgmock -pthread
USER=develop
COMMIT=Dev

all: clean s21_matrix_oop.a

s21_matrix_oop.a:$(SOURSES)
	$(CC) -c $(SOURSES)
	@ar rcs s21_matrix_oop.a $(OBJECT)

test: clean s21_matrix_oop.a
	@$(CC) -lgtest s21_matrix_oop.a tests.cc -o unit_test
	@./unit_test

gcov_report: clean
	$(CC) --coverage *.cc -o unit_test -lgtest 
	chmod +x *
	./unit_test
	lcov -t "test" -o unit_test.info --no-external -c -d .
	genhtml -o report/ unit_test.info
	open ./report/index.html

leak:
	@$(CC) tests.cc s21_matrix_oop.cc -lgtest -o leaks -g
	CK_FORK=no leaks --atExit -- ./leaks
	@rm ./leaks

style_check:
	clang-format -style=Google -n *.cc *.h

style_correct:
	clang-format -style=Google -i *cc *.h

git: git_add git_commit git_push

git_add:
	git add .

git_commit:
	git commit -m "$(COMMIT)"

git_push:
	git push origin $(USER)

clean:
	@rm -rf ./*.o ./*.a ./a.out ./GcovReport ./*.gcno ./*.gcda ./report ./*.info ./decimal_o ./*.dSYM ./unit_test ./CPPLINT*


