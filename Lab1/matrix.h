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
	Matrix(size_t rows, size_t cols, T val = 0) {
		if (r <= 0 || c <= 0) throw invalid_argument("Size of matrix must be greater than zero.");

		_rows = rows;
		_cols = cols;

		_matrix = new T[_rows];
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
		_matrix = new T[_rows];
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

	T& operator()(size_t r, size_t c) {
		if (r < 0 || r >= _rows || c < 0 || c >= _cols) throw out_of_range("Invalid index!");
		return &_matrix[r][c];
	}

	Matrix operator+(const Matrix& other) {
		if (_rows != other._rows || _cols != other._cols) throw invalid_argument("Sizes of matrixes must match.");
		Matrix sum_m(_rows, _cols, 0);
		for (int i = 0; i < _rows; ++i) {
			for (int j = 0; j < _cols; ++j) {
				sum_m[i][j] = _matrix[i][j] + other._matrix[i][j];
			}
		}
		return sum_m;
	}




};