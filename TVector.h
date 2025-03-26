#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

typedef std::vector<float> fvector;
typedef std::vector<std::vector<float>> fmatrix;

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
		for (auto i = 0; i < vector._vector.size() - 1; ++i) {
			os << vector._vector[i];
			os << ", ";
		}
		os << vector._vector.size() - 1 << "]" << std::endl;
		return os;
	}

private:
	std::vector<T> _vector;
};

template <typename T>
TVector<T> vec_mult(TVector<T> v, float a) {
	TVector<T> res(v.Size());
	for (auto i = 0; i < v.Size(); ++i) {
		res[i] = v[i] * a;
	}
	return res;
}

template <typename T>
fmatrix vec_mult(TVector<T> v1, TVector<T> v2) {
	auto N = v1.Size();
	fmatrix res(N, std::vector<float>(N));
	for (auto i = 0; i < N; ++i) {
		for (auto j = 0; j < N; ++j) {
			res[i][j] = v1[i] * v2[j];
		}
	}
	return res;
}
