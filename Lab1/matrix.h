#pragma once

#include <random>
#include<cstdlib>
#include <ctime>
#include <stdexcept>

using namespace std;

template<typename T>
class Matrix {
private:
	size_t _rows, _cols;
	T** _matrix;

public:
	Matrix(int rows, int cols, T val = 0) {
		if (rows <= 0 || cols <= 0) throw invalid_argument("Size of matrix must be greater than zero.");

		_rows = rows;
		_cols = cols;

		_matrix = new T*[_rows];
		for (int i = 0; i < _rows; ++i) {
			_matrix[i] = new T[_cols];
			for (int j = 0; j < _cols; ++j) {
				_matrix[i][j] = val;
			}
		}
	}

	Matrix(size_t rows, size_t cols, T max_val, T min_val){
		if (rows <= 0 || cols <= 0) throw invalid_argument("Size of matrix must be greater than zero.");

		_rows = rows;
		_cols = cols;

		srand(time(0));
		_matrix = new T * [_rows];
		for (int i = 0; i < rows; ++i) {
			_matrix[i] = new T[_cols];
			for (int j = 0; j < cols; ++j) {
				_matrix[i][j] = ((T)rand() / RAND_MAX) * (max_val - min_val) + min_val;
			}
		}
	}

	Matrix(const Matrix& other) {
		_rows = other._rows;
		_cols = other._cols;

		_matrix = new T * [_rows];
		for (int i = 0; i < _rows; ++i) {
			_matrix[i] = new T[_cols];
			for (int j = 0; j < _cols; ++j) {
				_matrix[i][j] = other._matrix[i][j];
			}
		}
	}

	~Matrix(){
		for (int i = 0; i < _rows; ++i) {
			delete[] _matrix[i];
		}
		delete[] _matrix;
	}

	size_t rows() const {
		return _rows;
	}

	size_t cols() const {
		return _cols;
	}

	T& operator()(size_t row, size_t col) {
		if (row < 0 || row >= _rows || col < 0 || col >= _cols) throw out_of_range("Invalid index!");
		return _matrix[row][col];
	}

	T operator()(size_t row, size_t col) const {
		if (row < 0 || row >= _rows || col < 0 || col >= _cols) throw out_of_range("Invalid index!");
		return _matrix[row][col];
	}

	Matrix operator+(const Matrix& other) {
		if (_rows != other._rows || _cols != other._cols) throw invalid_argument("Sizes of matrixes must match.");
		Matrix res_m(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res_m(i,j) = _matrix[i][j] + other._matrix[i][j];
			}
		}
		return res_m;
	}

	Matrix operator-(const Matrix& other) {
		if (_rows != other._rows || _cols != other._cols) throw invalid_argument("Sizes of matrixes must match.");
		Matrix res_m(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res_m(i, j) = _matrix[i][j] - other._matrix[i][j];
			}
		}
		return res_m;
	}

	Matrix operator*(const Matrix& other) {
		if (_cols != other._rows) throw std::invalid_argument("The number of columns of the first matrix "
																"must match the number of rows of the second matrix");
		Matrix res_m(_rows, other._cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < other._cols; ++j) {
				for (int k = 0; k < _cols; ++k) {
					res_m(i, j) += _matrix[i][k] * other._matrix[k][j];
				}	
			}
		}
		return res_m;

	}

	Matrix operator*(const T val) {
		Matrix res_m(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res_m(i, j) = _matrix[i][j] * val;
			}
		}
		return res_m;
	}

};


template <typename T>
ostream& operator<<(ostream& os, const Matrix<T>& m) {
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.cols(); ++j) {
			os << m(i, j) << '\t';
		}
		os << '\n';
	}
	return os;
}

template <typename T>
Matrix<T> operator*(const T val, Matrix <T>& m) {
	return m * val;
}