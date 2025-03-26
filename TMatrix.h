#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

#include "TVector.h"

template <typename T> class TMatrix {
public:
    TMatrix<T>(ptrdiff_t n, ptrdiff_t m, T def_value = T()) : _matrix(n, TVector<T>(m def_value)) {}
	
	TMatrix<T>(std::initializer_list<TVector<T>> l) : _matrix(l) {}

    TMatrix<T> operator-(TMatrix<T>& m) {
        TMatrix<T> res(*this);
        auto N = _matrix.Size();
        for (auto i = 0; i < N; ++i) {
            for (auto j = 0; j < N; ++j) {
                res._matrix[i][j] -= m._matrix[i][j];
            }
        }
        return res;
    }

    TVector<T> operator*(TVector<T> v) {
        auto N = _matrix.size();
        TVector<T> res(N, 0);
        for (auto i = 0; i < N; ++i) {
            for (auto j = 0; j < N; ++j) {
                res[i] += _matrix[i][j] * v[j];
            }
        }
        return res;
    }

    TMatrix<T> operator*(TMatrix<T> m) {
        auto N = _matrix.size();
        auto M = m._matrix[0].size();
        TMatrix<T> res(N, TVector<T>(M, 0));
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
        auto rows = _matrix.Size();
        auto cols = _matrix[0].Size();
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                _matrix[i][j] = mas[j][i];
            }
        }
    }
    
	auto Size() const {
		return _matrix.size();
	}

	T& operator[](ptrdiff_t index) {
		return _matrix[index];
	}

	const T& operator[](ptrdiff_t index) const {
		return _matrix[index];
	}

	template <typename T>
	friend std::ostream& operator<<(std::ostream& os, const TMatrix<T>& m) {
		os << "[";
        auto N = m._matrix.size();
        auto M = m._matrix[0].size();
		for (auto i = 0; i < N; ++i) {
            for (auto j = 0; j < M - 1; ++j) {
			    os << m._matrix[i][j] << ", ";
            }
            if (i == N - 1) {
                os <<  m._matrix[i] << std::endl;
            } else {
		        os << vector._vector.size() - 1 << "]" << std::endl;
            }
		}
		return os;
	}


private:
    TVector<TVector<T>> _matrix;
}
