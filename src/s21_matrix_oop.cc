#include "s21_matrix_oop.h"

// CONSTRUCTORS & DESTRUCTORS
S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::out_of_range(
        "Matrix size cannot contain negative values and zero.");
  }
  MemoryAllocation();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  MemoryAllocation();
  for (int row = 0; row < rows_; ++row) {
    std::memcpy(matrix_[row], other.matrix_[row], other.cols_ * sizeof(double));
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    for (int row = 0; row < rows_; ++row) {
      delete[] matrix_[row];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

void S21Matrix::MemoryAllocation() {
  matrix_ = new double*[rows_]();
  for (int row = 0; row < rows_; row++) {
    matrix_[row] = new double[cols_]();
  }
}

// ACCESSORS & MUTATORS
int S21Matrix::get_rows() const { return rows_; }

int S21Matrix::get_cols() const { return cols_; }

void S21Matrix::set_cols(int cols) {
  if (cols <= 0) {
    throw std::out_of_range(
        "Matrix size cannot contain negative values and zero.");
  }
  S21Matrix res_matrix(rows_, cols);
  for (int row = 0; row < rows_; ++row) {
    std::memcpy(res_matrix.matrix_[row], matrix_[row],
                res_matrix.cols_ * sizeof(double));
  }
  *this = res_matrix;
}

void S21Matrix::set_rows(int rows) {
  if (rows <= 0) {
    throw std::out_of_range(
        "Matrix size cannot contain negative values and zero.");
  }
  S21Matrix res_matrix(rows, cols_);
  for (int row = 0; row < res_matrix.rows_ && row < rows_; ++row) {
    std::memcpy(res_matrix.matrix_[row], matrix_[row], cols_ * sizeof(double));
  }
  *this = res_matrix;
}

// SIMPLE OPERATIONS
bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool result = true;
  if (matrix_ != other.matrix_) {
    if (rows_ == other.rows_ && cols_ == other.cols_) {
      for (int i = 0; i < rows_ && result; ++i) {
        for (int j = 0; j < cols_ && result; ++j) {
          if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) result = false;
        }
      }
    } else {
      result = false;
    }
  }
  return result;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("Incompatible matrix sizes for addition.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("Incompatible matrix sizes for subtraction.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument(
        "Matrix dimensions are not compatible for multiplication.");
  }
  S21Matrix res_matrix(rows_, other.cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      double sum_elem = 0;
      for (int k = 0; k < other.rows_; ++k) {
        sum_elem += matrix_[i][k] * other.matrix_[k][j];
      }
      res_matrix.matrix_[i][j] = sum_elem;
    }
  }
  *this = res_matrix;
}

// OVERLOADING OPERATORS
S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  res_matrix.SumMatrix(other);
  return res_matrix;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  res_matrix.SubMatrix(other);
  return res_matrix;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  res_matrix.MulMatrix(other);
  return res_matrix;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res_matrix(*this);
  res_matrix.MulNumber(num);
  return res_matrix;
}

bool S21Matrix::operator==(const S21Matrix& other) {
  S21Matrix res_matrix(*this);
  return res_matrix.EqMatrix(other);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;
  }
  this->~S21Matrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  MemoryAllocation();
  CopyElem(other);
  return *this;
}

void S21Matrix::CopyElem(const S21Matrix& other) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    this->~S21Matrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
    throw std::out_of_range("Incorrect index values.");
  }
  return matrix_[i][j];
}

double S21Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
    throw std::out_of_range("Incorrect index values.");
  }
  return matrix_[i][j];
}

// OPERATIONS
S21Matrix S21Matrix::Transpose() {
  S21Matrix res_matrix(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      res_matrix.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res_matrix;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::out_of_range(
        "Incompatible matrix sizes to search for a determinant.");
  }
  S21Matrix res_matrix(rows_, cols_);
  if (rows_ == 1) {
    res_matrix.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int row = 0; row < rows_; ++row) {
      for (int col = 0; col < cols_; ++col) {
        res_matrix.matrix_[row][col] =
            pow(-1, row + col) * Minor(row, col).Determinant();
      }
    }
  }
  return res_matrix;
}

S21Matrix S21Matrix::Minor(int deleted_row, int deleted_col) {
  S21Matrix res_matrix(rows_ - 1, cols_ - 1);
  int new_row = 0;
  for (int row = 0; row < rows_; ++row) {
    if (row != deleted_row) {
      int new_col = 0;
      for (int col = 0; col < cols_; ++col) {
        if (col != deleted_col) {
          res_matrix.matrix_[new_row][new_col] = matrix_[row][col];
          ++new_col;
        }
      }
      ++new_row;
    }
  }
  return res_matrix;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::out_of_range(
        "Incompatible matrix sizes to search for a determinant.");
  }
  double determinant = 0;
  if (rows_ == 1) {
    determinant = matrix_[0][0];
  } else if (rows_ == 2) {
    determinant = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int col = 0; col < cols_; ++col) {
      determinant +=
          pow(-1, col + 2) * matrix_[0][col] * Minor(0, col).Determinant();
    }
  }
  return determinant;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determinant = Determinant();
  if (fabs(determinant) < EPS || rows_ != cols_) {
    throw std::out_of_range(
        "Incompatible matrix sizes to search inverse matrix.");
  }
  S21Matrix res_matrix(rows_, cols_);
  if (rows_ == 1) {
    res_matrix.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    res_matrix = CalcComplements().Transpose();
    res_matrix.MulNumber(1 / determinant);
  }
  return res_matrix;
}
