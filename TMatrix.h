#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

#include "TVector.h"

template <typename T> class TMatrix {
public:
	TMatrix<T>() : _matrix() {}

	TMatrix<T>(ptrdiff_t n, ptrdiff_t m, T def_value = T()) : _matrix(n, TVector<T>(m, def_value)) {}

	TMatrix<T>(std::initializer_list<TVector<T>> l) : _matrix(l) {}

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

	TMatrix<T> operator*(T a) const {
		auto N = _matrix.Size();
		auto M = _matrix[0].Size();
		TMatrix<T> res(N, M, 0);
		for (auto i = 0; i < N; ++i) {
			for (auto j = 0; j < M; ++j) {
				res[i][j] = _matrix[i][j] * a;
			}
		}
		return res;
	}

	TVector<T> operator*(const TVector<T>& v) const {
		auto N = _matrix.Size();
		TVector<T> res(N, 0);
		for (auto i = 0; i < N; ++i) {
			for (auto j = 0; j < N; ++j) {
				res[i] += _matrix[i][j] * v[j];
			}
		}
		return res;
	}

	TMatrix<T> operator*(const TMatrix<T>& m) const {
		auto N = _matrix.Size();
		auto M = m._matrix[0].Size();
		TMatrix<T> res(N, M, 0);
		for (auto i = 0; i < N; ++i) {
			for (auto j = 0; j < M; ++j) {
				for (auto k = 0; k < N; ++k) {
					res[i][j] += _matrix[i][k] * m._matrix[k][j];
				}
			}
		}
		return res;
	}

	void Transpone() {
		const auto N = _matrix.Size();
		const auto M = _matrix[0].Size();
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < M; j++) {
				auto temp = _matrix[i][j];
				_matrix[i][j] = _matrix[j][i];				
				_matrix[j][i] = temp;
			}
		}
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
TMatrix<T> operator*(TVector<T> v1, TVector<T> v2) {
	auto N = v1.Size();
	auto M = v2.Size();
	TMatrix<T> res(N, N, 0);
	for (auto i = 0; i < N; ++i) {
		for (auto j = 0; j < M; ++j) {
			res[i][j] = v1[i] * v2[j];
		}
	}
	return res;
}
