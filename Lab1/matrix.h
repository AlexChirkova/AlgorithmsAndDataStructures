#pragma once

#include <random>
#include<cstdlib>
#include <ctime>
#include <stdexcept>
#include <iomanip>
#include <complex>

using namespace std;

static constexpr double EPSILON = 1e-5;

template<typename T>
class Matrix {
private:
	size_t _rows, _cols;
	T** _matrix;

public:
	Matrix(size_t rows, size_t cols, T val = 0) {
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

	Matrix(size_t rows, size_t cols, T min_val, T max_val){
		if (rows <= 0 || cols <= 0) throw invalid_argument("Size of matrix must be greater than zero.");

		if ((min_val) > (max_val)) swap(min_val, max_val);

		_rows = rows;
		_cols = cols;

		srand(static_cast<unsigned int>(time(0)));
		random_device rd;
		mt19937 generator(rd());
		uniform_real_distribution<T> distribution(min_val, max_val);
		
		_matrix = new T * [_rows];
		for (int i = 0; i < rows; ++i) {
			_matrix[i] = new T[_cols];
			for (int j = 0; j < cols; ++j) {
				_matrix[i][j] = static_cast<T>(distribution(generator));
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

	Matrix operator/(const T val) {
		Matrix res_m(_rows, _cols);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				res_m(i, j) = _matrix[i][j] / val;
			}
		}
		return res_m;
	}

	bool operator==(const Matrix& other) const{
		if (_rows != other._rows || _cols != other._cols) return false;

		size_t size = _rows * _cols;
		for (size_t i = 0; i < size; ++i){
			if (std::abs(_matrix[i] - other._matrix[i]) > EPSILON) return false;
		}
		return true;
	}

	bool operator!=(const Matrix& other) const
	{
		return !(*this == other);
	}

	T trace() {
		if (_rows != _cols) throw std::invalid_argument("The trace can only be calculated for a square matrix");
		
		T trace = 0;
		for (int i = 0; i < _rows; ++i) {
			trace += _matrix[i][i];
		}
		return trace;
	}

	Matrix& downTr() {
		/*if (_rows != _cols) throw std::invalid_argument("Matrix mast be square to downTr it.");

		for (int k = int(_rows)-1; k >=0; --k) {
			T point = _matrix[k][k];
			if (point == 0) throw runtime_error("The zero element on the main diagonal.");

			for (int i = k - 1; i >= 0; --i) {
				T ratio = _matrix[i][k] / _matrix[k][k];
				for (int j = k; j >= 0; --j) {
					_matrix[i][j] -= ratio * _matrix[k][j];
				}
			}
		}
		return *this;*/
		
		if (_rows != _cols) throw std::invalid_argument("Only square matrices are reduced to a triangular form!");

		size_t n = _rows;
		for (size_t i = n; i > 0; --i)
		{
			size_t row = i - 1;

			if (std::abs(_matrix[row][row]) < EPSILON)
			{
				bool swapped = false;
				for (size_t j = row - 1; j < n; --j)
				{
					if (std::abs(_matrix[j][row]) > EPSILON)
					{
						for (size_t k = 0; k < n; ++k)
						{
							std::swap(_matrix[row][k], _matrix[j][k]);
						}
						swapped = true;
						break;
					}
				}
				if (!swapped)
				{
					throw std::runtime_error("Matrix is singular and cannot be reduced!");
				}
			}

			for (size_t j = 0; j < row; ++j)
			{
				T multiplier = _matrix[j][row] / _matrix[row][row];
				for (size_t k = 0; k < n; ++k)
				{
					_matrix[j][k] -= multiplier * _matrix[row][k];
				}

			}
			
		}
		return *this;
	}

};


template <typename T>
ostream& operator<<(ostream& os, const Matrix<T>& m) {
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.cols(); ++j) {
			os << m(i, j) << setw(10);
		}
		os << '\n';
	}
	return os;
}

template <typename T>
Matrix<T> operator*(const T val, Matrix <T>& m) {
	return m * val;
}