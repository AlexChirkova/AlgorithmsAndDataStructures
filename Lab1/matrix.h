#include <random>
#include<cstdlib>
#include <ctime>


template<typename T>
class Matrix {
private:
	size_t _rows, _cols;
	T _matrix[_rows][_cols];

public:
	Matrix(size_t r, size_t c, T val): _rows(r), _cols(c) {
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				_matrix[i][j] = val;
			}
		}
	}

	Matrix(size_t r, size_t c, T max_val, T min_val): _rows(r), _cols(c) {
		srand(time(0));
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				_matrix[i][j] = ((T)rand() / RAND_MAX) * (max_val - min_val) + min_val;
			}
		}
	}

};