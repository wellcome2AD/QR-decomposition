#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

template <typename T> class TVector {
public:
	TVector<T>() : _vector() {}

	TVector<T>(ptrdiff_t size, T def_value = T()) : _vector(size, def_value) {}

	TVector<T>(std::initializer_list<T> l) : _vector(l) {}

	TVector<T>& operator*(T a) {
#pragma omp parallel for num_treads(thread_num)
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

	bool operator==(const TVector<T>& v) const {
		const float eps = 0.000001;
		for (auto i = 0; i < _vector.size(); ++i) {
			if (abs(_vector[i] - v._vector[i]) > eps) {
				return false;
			}
		}
		return true;
	}

	bool operator!=(const TVector<T>& v) const {
		return !(*this == v);
	}

	friend std::ostream& operator<<(std::ostream& os, const TVector<T>& vector) {
		os << "[";
		auto N = vector._vector.size();
		for (auto i = 0; i < N - 1; ++i) {
			os << vector._vector[i];
			os << ", ";
		}
		os << vector._vector[N - 1] << "]";
		return os;
	}

private:
	std::vector<T> _vector;
};
