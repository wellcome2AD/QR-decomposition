#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

#include "TVector.h"

size_t thread_num = 8;

template <typename T> class TMatrix {
public:
	TMatrix<T>() : _matrix() {}

	TMatrix<T>(ptrdiff_t n, ptrdiff_t m, T def_value = T()) : _matrix(n, TVector<T>(m, def_value)) {}

	TMatrix<T>(std::initializer_list<TVector<T>> l) : _matrix(l) {}

	TMatrix<T>(const TMatrix<T>& other) : _matrix(other._matrix) {}

	TMatrix<T> operator-(const TMatrix<T>& m) const {
		TMatrix<T> res(*this);
		auto N = _matrix.Size();
		for (auto i = 0; i < N; ++i) {
			for (auto j = 0; j < N; ++j) {
				res._matrix[i][j] -= m._matrix[i][j];
			}
		}
		return res;
	}

	TMatrix<T> operator*(const TMatrix<T>& m) const {
		assert(Size() == m.Size());
		auto N = Size();
		TMatrix<T> res(N, N, 0);
		int s = 50;
#pragma omp parallel for num_threads(thread_num)
		for (int i = 0; i < N; i++) {
			for (int jj = 0; jj < N; jj += s) {
				for (int kk = 0; kk < N; kk += s) {
					for (int j = jj; j < ((jj + s) > N ? N : (jj + s)); j++) {
						T temp = 0;
						for (int k = kk; k < ((kk + s) > N ? N : (kk + s)); k++) {
							temp += _matrix[i][k] * m[k][j];
						}
						res[i][j] += temp;
					}
				}
			}
		}
		return res;
	}

	TMatrix<T>& operator=(const TMatrix<T>& other) {
		_matrix = other._matrix;
		return *this;
	}

	bool operator==(const TMatrix<T>& m) const {
		auto N = _matrix.Size();
		auto M = m._matrix[0].Size();
		const float eps = 0.000001;
		for (auto i = 0; i < N; ++i) {
			for (auto j = 0; j < M; ++j) {
				if (abs(_matrix[i][j] - m._matrix[i][j]) > eps) {
					return false;
				}
			}
		}
		return true;
	}

	bool operator!=(const TMatrix<T>& v) const {
		return !(*this == v);
	}


	auto Size() const {
		return _matrix.Size();
	}

	TVector<T>& operator[](ptrdiff_t index) {
		return _matrix[index];
	}

	const TVector<T>& operator[](ptrdiff_t index) const {
		return _matrix[index];
	}

	friend std::ostream& operator<<(std::ostream& os, const TMatrix<T>& m) {
		os << "[";
		auto N = m._matrix.Size();
		auto M = m._matrix[0].Size();
		for (auto i = 0; i < N; ++i) {
			std::cout << std::left << std::setw(5);
			for (auto j = 0; j < M - 1; ++j) {
				os << m._matrix[i][j] << ", ";
			}
			os << m._matrix[i][M - 1];
			if (i == N - 1) {
				os << "]";
			}
			os << std::endl;
		}
		return os;
	}


private:
	TVector<TVector<T>> _matrix;
};

template <typename T>
TMatrix<T> generate_matrix(int N, int M) {
	TMatrix<T> res(N, M);
	std::srand(std::time(0));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < M; ++j) {
			res[i][j] = std::rand() % 20 - 10;
		}
	}
	return res;
}
