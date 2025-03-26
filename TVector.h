#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

#include "TMatrix.h"

template <typename T> class TVector {
public:
	TVector<T>(ptrdiff_t size, T def_value = T()) : _vector(size, def_value) {}
	
	TVector<T>(std::initializer_list<T> l) : _vector(l) {}

	TVector<T>& operator*(T a) {
		for (auto i = 0; i < _vector.size(); ++i) {
			_vector[i] = _vector[i] * a;
		}
		return *this;
	}

	TMatrix<T> operator*(TVector<T> v) {
		auto N = _vector.size();
		TMatrix<T> res(N, N, 0);
		for (auto i = 0; i < N; ++i) {
			for (auto j = 0; j < N; ++j) {
				res[i][j] = _vector[i] * v._vector[j];
			}
		}
		return res;
	}

	auto Size() const {
		return _vector.size();
	}

	T& operator[](ptrdiff_t index) {
		return _vector[index];
	}

	const T& operator[](ptrdiff_t index) const {
		return _vector[index];
	}

	template <typename T>
	friend std::ostream& operator<<(std::ostream& os, const TVector<T>& vector) {
		os << "[";
		auto N = vector._vector.size();
		for (auto i = 0; i < N - 1; ++i) {
			os << vector._vector[i];
			os << ", ";
		}
		os << vector._vector[N-1] << "]" << std::endl;
		return os;
	}

private:
	std::vector<T> _vector;
};
