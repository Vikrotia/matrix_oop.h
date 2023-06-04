#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

class S21MatrixTest : public testing::Test {
 protected:
  S21Matrix *matrix_1x1;
  S21Matrix *matrix_3x3;
  S21Matrix *matrix_2x3;

  void SetUp() {
    matrix_1x1 = new S21Matrix(1, 1);
    matrix_1x1->operator()(0, 0) = 10;

    matrix_3x3 = new S21Matrix(3, 3);
    matrix_3x3->operator()(0, 0) = 1;
    matrix_3x3->operator()(0, 1) = 0;
    matrix_3x3->operator()(0, 2) = -2;
    matrix_3x3->operator()(1, 0) = 3;
    matrix_3x3->operator()(1, 1) = 2;
    matrix_3x3->operator()(1, 2) = 1;
    matrix_3x3->operator()(2, 0) = 1;
    matrix_3x3->operator()(2, 1) = 2;
    matrix_3x3->operator()(2, 2) = -2;

    matrix_2x3 = new S21Matrix(2, 3);
    matrix_2x3->operator()(0, 0) = 0.0;
    matrix_2x3->operator()(0, 1) = 1.0;
    matrix_2x3->operator()(0, 2) = 2.0;
    matrix_2x3->operator()(1, 0) = 3.0;
    matrix_2x3->operator()(1, 1) = 4.0;
    matrix_2x3->operator()(1, 2) = 5.0;
  }

  void TearDown() {
    delete matrix_1x1;
    delete matrix_3x3;
    delete matrix_2x3;
  }
};

TEST_F(S21MatrixTest, ExceptionConstructor1) {
  EXPECT_THROW(S21Matrix matrix(-2, 8), std::out_of_range);
}
TEST_F(S21MatrixTest, ExceptionConstructor2) {
  EXPECT_THROW(S21Matrix matrix(0, 0), std::out_of_range);
}

TEST_F(S21MatrixTest, Constructor) {
  S21Matrix matrix;
  EXPECT_EQ(0, matrix.get_rows());
  EXPECT_EQ(0, matrix.get_cols());
}

TEST_F(S21MatrixTest, ConstructorWithParameters) {
  S21Matrix matrix(10, 3);
  EXPECT_EQ(10, matrix.get_rows());
  EXPECT_EQ(3, matrix.get_cols());
}

TEST_F(S21MatrixTest, ExceptionSetCols) {
  S21Matrix matrix(10, 3);
  EXPECT_THROW(matrix.set_cols(-10);, std::out_of_range);
}

TEST_F(S21MatrixTest, ExceptionSetRows) {
  S21Matrix matrix(10, 3);
  EXPECT_THROW(matrix.set_rows(-3);, std::out_of_range);
}

TEST_F(S21MatrixTest, SetGet1) {
  S21Matrix matrix(10, 10);
  S21Matrix result(1, 1);
  matrix.set_rows(1);
  matrix.set_cols(1);
  EXPECT_EQ(1, result.get_rows());
  EXPECT_EQ(1, result.get_cols());
}

TEST_F(S21MatrixTest, ConstructorCopy1) {
  S21Matrix matrix(*matrix_2x3);
  EXPECT_EQ(2, matrix.get_rows());
  EXPECT_EQ(3, matrix.get_cols());
  EXPECT_EQ(2, matrix_2x3->get_rows());
  EXPECT_EQ(3, matrix_2x3->get_cols());
  EXPECT_EQ(1, *matrix_2x3 == matrix);
}

TEST_F(S21MatrixTest, ConstructorCopy2) {
  S21Matrix matrix1;
  S21Matrix matrix(matrix1);
  EXPECT_EQ(0, matrix.get_rows());
  EXPECT_EQ(0, matrix.get_cols());
  EXPECT_EQ(0, matrix1.get_rows());
  EXPECT_EQ(0, matrix1.get_cols());
  EXPECT_EQ(1, matrix1 == matrix);
}

TEST_F(S21MatrixTest, ConstructorMove) {
  S21Matrix matrix1(*matrix_2x3);
  S21Matrix matrix(std::move(*matrix_2x3));
  EXPECT_EQ(2, matrix.get_rows());
  EXPECT_EQ(3, matrix.get_cols());
  EXPECT_EQ(1, matrix1 == matrix);
  EXPECT_EQ(0, matrix_2x3->get_rows());
  EXPECT_EQ(0, matrix_2x3->get_cols());
}

TEST_F(S21MatrixTest, EqualsOperator1) {
  S21Matrix matrix1(4, 3), matrix2(4, 3);

  EXPECT_EQ(0, *matrix_2x3 == *matrix_3x3);
  EXPECT_EQ(0, *matrix_2x3 == matrix2);
  EXPECT_EQ(1, matrix1 == matrix2);
}

TEST_F(S21MatrixTest, EqualsOperator2) {
  S21Matrix matrix1(4, 3);
  matrix1 = matrix1;
  S21Matrix result(4, 3);
  result = matrix1;
  EXPECT_EQ(1, result == matrix1);
}

TEST_F(S21MatrixTest, ExceptionSumTest) {
  S21Matrix test_result1(3, 3);
  S21Matrix test_result2(2, 3);
  EXPECT_THROW(test_result1.SumMatrix(test_result2), std::invalid_argument);
}

TEST_F(S21MatrixTest, SumOperator) {
  S21Matrix result(2, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 2.0;
  result(0, 2) = 4.0;
  result(1, 0) = 6.0;
  result(1, 1) = 8.0;
  result(1, 2) = 10.0;

  S21Matrix test_result(3, 3);
  test_result = *matrix_2x3 + *matrix_2x3;
  EXPECT_EQ(1, result == test_result);
}

TEST_F(S21MatrixTest, SumOverloading) {
  S21Matrix result(2, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 2.0;
  result(0, 2) = 4.0;
  result(1, 0) = 6.0;
  result(1, 1) = 8.0;
  result(1, 2) = 10.0;

  *matrix_2x3 += *matrix_2x3;
  EXPECT_EQ(1, result == *matrix_2x3);
}

TEST_F(S21MatrixTest, ExceptionSubTest) {
  S21Matrix test_result1(3, 3);
  S21Matrix test_result2(2, 3);
  EXPECT_THROW(test_result1.SubMatrix(test_result2), std::invalid_argument);
}

TEST_F(S21MatrixTest, SubOperator) {
  S21Matrix result(3, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 0.0;
  result(0, 2) = 0.0;
  result(1, 0) = 0.0;
  result(1, 1) = 0.0;
  result(1, 2) = 0.0;

  S21Matrix result_test(3, 3);
  result_test = *matrix_3x3 - *matrix_3x3;
  EXPECT_EQ(1, result == result_test);
}

TEST_F(S21MatrixTest, SubOverloading) {
  matrix_2x3->operator()(0, 0) = 1;
  S21Matrix result(2, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 0.0;
  result(0, 2) = 0.0;
  result(1, 0) = 0.0;
  result(1, 1) = 0.0;
  result(1, 2) = 0.0;

  *matrix_2x3 -= *matrix_2x3;
  EXPECT_EQ(1, result == *matrix_2x3);
}

TEST_F(S21MatrixTest, ExceptionMulMatrixTest) {
  S21Matrix test_result1(3, 3);
  S21Matrix test_result2(2, 3);
  EXPECT_THROW(test_result1.MulMatrix(test_result2), std::invalid_argument);
}

TEST_F(S21MatrixTest, MulMatrixOperator) {
  S21Matrix result(2, 3);
  result(0, 0) = 5.0;
  result(0, 1) = 6.0;
  result(0, 2) = -3.0;
  result(1, 0) = 20.0;
  result(1, 1) = 18.0;
  result(1, 2) = -12.0;

  S21Matrix result_test(2, 3);
  result_test = *matrix_2x3 * (*matrix_3x3);
  EXPECT_EQ(1, result == result_test);
}

TEST_F(S21MatrixTest, MulMatrixOverloading) {
  S21Matrix result(2, 3);
  result(0, 0) = 5.0;
  result(0, 1) = 6.0;
  result(0, 2) = -3.0;
  result(1, 0) = 20.0;
  result(1, 1) = 18.0;
  result(1, 2) = -12.0;

  *matrix_2x3 *= *matrix_3x3;
  EXPECT_EQ(1, result == *matrix_2x3);
}

TEST_F(S21MatrixTest, MulNumOperator) {
  S21Matrix result(2, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 5.6;
  result(0, 2) = 11.2;
  result(1, 0) = 16.8;
  result(1, 1) = 22.4;
  result(1, 2) = 28.0;

  S21Matrix result_test(2, 3);
  result_test = *matrix_2x3 * 5.6;
  EXPECT_EQ(1, result == result_test);
}

TEST_F(S21MatrixTest, MulNumOverloading) {
  S21Matrix result(2, 3);
  result(0, 0) = 0.0;
  result(0, 1) = 0.5;
  result(0, 2) = 1.0;
  result(1, 0) = 1.5;
  result(1, 1) = 2.0;
  result(1, 2) = 2.5;

  *matrix_2x3 *= 0.5;
  EXPECT_EQ(1, result == *matrix_2x3);
}

TEST_F(S21MatrixTest, ExceptionOperatorBrascet1) {
  S21Matrix result(2, 3);
  EXPECT_THROW(result(-4, 2) = 1, std::out_of_range);
}
TEST_F(S21MatrixTest, ExceptionOperatorBrascet2) {
  S21Matrix result(3, 3);
  EXPECT_THROW(result(3, 3) = 1, std::out_of_range);
}

TEST_F(S21MatrixTest, ExceptionOperatorBrascetConst) {
  const S21Matrix result(2, 3);
  double elem = 3.0;
  EXPECT_THROW(elem = result(-4, 2), std::out_of_range);
  EXPECT_EQ(1, 3.0 == elem);
}

TEST_F(S21MatrixTest, OperatorBrascetConst) {
  double elem = 3.0;
  const S21Matrix matrix_const(*matrix_2x3);
  elem = matrix_const(0, 0);
  EXPECT_EQ(1, 0.0 == elem);
}

TEST_F(S21MatrixTest, Transpose) {
  S21Matrix result(3, 2);
  result(0, 0) = 0.0;
  result(0, 1) = 3.0;
  result(1, 0) = 1.0;
  result(1, 1) = 4.0;
  result(2, 0) = 2.0;
  result(2, 1) = 5.0;

  S21Matrix test_result(3, 2);
  test_result = matrix_2x3->Transpose();
  EXPECT_EQ(1, result == test_result);
}

TEST_F(S21MatrixTest, ExceptionCalcCompliments) {
  S21Matrix result(2, 3);
  EXPECT_THROW(result.CalcComplements(), std::out_of_range);
}

TEST_F(S21MatrixTest, CalcComplements1) {
  S21Matrix result(1, 1);
  result(0, 0) = 10;
  S21Matrix test_result(1, 1);
  test_result = matrix_1x1->CalcComplements();
  EXPECT_EQ(1, result == test_result);
}

TEST_F(S21MatrixTest, CalcComplements2) {
  S21Matrix result(3, 3);
  result(0, 0) = -6.0;
  result(0, 1) = 7.0;
  result(0, 2) = 4.0;
  result(1, 0) = -4.0;
  result(1, 1) = 0.0;
  result(1, 2) = -2.0;
  result(2, 0) = 4.0;
  result(2, 1) = -7.0;
  result(2, 2) = 2.0;
  S21Matrix test_result(3, 3);
  test_result = matrix_3x3->CalcComplements();
  EXPECT_EQ(1, result == test_result);
}

TEST_F(S21MatrixTest, ExceptionDeterminant) {
  S21Matrix result(2, 3);
  EXPECT_THROW(result.Determinant(), std::out_of_range);
}

TEST_F(S21MatrixTest, Determinant1) {
  EXPECT_EQ(1, -14.0 == matrix_3x3->Determinant());
}

TEST_F(S21MatrixTest, Determinant2) {
  EXPECT_EQ(1, 10.0 == matrix_1x1->Determinant());
}

TEST_F(S21MatrixTest, ExceptionInverseMatrix) {
  S21Matrix result(3, 3);
  result(0, 0) = 1.0;
  result(0, 1) = 2.0;
  result(0, 2) = 3.0;
  result(1, 0) = 4.0;
  result(1, 1) = 5.0;
  result(1, 2) = 6.0;
  result(2, 0) = 7.0;
  result(2, 1) = 8.0;
  result(2, 2) = 9.0;
  EXPECT_THROW(result.InverseMatrix(), std::out_of_range);
}

TEST_F(S21MatrixTest, InverseMatrix1) {
  S21Matrix result(1, 1);
  result(0, 0) = 0.1;
  EXPECT_EQ(1, result == matrix_1x1->InverseMatrix());
}

TEST_F(S21MatrixTest, InverseMatrix2) {
  S21Matrix matrix_2x2(2, 2);
  matrix_2x2(0, 0) = 3.0;
  matrix_2x2(0, 1) = 4.0;
  matrix_2x2(1, 0) = 5.0;
  matrix_2x2(1, 1) = 7.0;

  S21Matrix result(2, 2);
  result(0, 0) = 7.0;
  result(0, 1) = -4.0;
  result(1, 0) = -5.0;
  result(1, 1) = 3.0;

  EXPECT_EQ(1, result == matrix_2x2.InverseMatrix());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
