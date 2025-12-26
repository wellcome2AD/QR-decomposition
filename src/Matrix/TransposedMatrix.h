#pragma once

#include <initializer_list>
#include <ostream>
#include <vector>

template <typename T> class TransposedMatrix {
public:
	TransposedMatrix() : _matrix(0) {}

	TransposedMatrix(ptrdiff_t n, ptrdiff_t m, T def_value = T()) : _n(n), _m(m), _matrix(n* m, def_value) {}

	TransposedMatrix(std::initializer_list<std::vector<std::vector<T>>> l)
	{
		const auto& rows = *l.begin();
		_n = rows.size();
		_m = rows[0].size();
		_matrix = std::vector<T>(_n * _m);
#pragma omp parallel for num_threads(thread_num) if (_n >= 1000)
		for (int j = 0; j < _m; ++j)
		{
			for (size_t i = 0; i < _n; ++i)
			{
				_matrix[j * _m + i] = rows[i][j];
			}
		}
	}

	TransposedMatrix(const std::vector<std::vector<T>>& other)
	{
		_n = other.size();
		_m = other[0].size();
		_matrix = std::vector<T>(_n * _m);
#pragma omp parallel for num_threads(thread_num) if (_n >= 1000)
		for (int j = 0; j < _m; ++j)
		{
			for (size_t i = 0; i < _n; ++i)
			{
				_matrix[j * _m + i] = other[i][j];
			}
		}
	}

	TransposedMatrix(const TransposedMatrix & other) : _matrix(other._matrix) {}

	~TransposedMatrix() = default;

	size_t GetN() const
	{
		return _n;
	}

	size_t GetM() const
	{
		return _m;
	}

	T& At(size_t row, size_t col)
	{
		return _matrix[col * _m + row];
	}

	const T& At(size_t row, size_t col) const
	{
		return _matrix[col * _m + row];
	}

	std::vector<std::vector<T>> Transpose() const
	{
		auto res = std::vector<std::vector<T>>(_n, std::vector<T>(_m, 0));
		for (size_t j = 0; j < _m; ++j)
		{
			for (size_t i = 0; i < _n; ++i)
			{
				res[i][j] = _matrix[j * _m + i];
			}
		}
		return res;
	}

	friend std::ostream& operator<<(std::ostream& os, const TransposedMatrix<T>& m) {
		os << "[\n  ";
		auto N = m._matrix.GetN();
		auto M = m._matrix.GetM();
		for (auto i = 0; i < N; ++i) {
			std::cout << std::left << std::setw(7);
			for (auto j = 0; j < M - 1; ++j) {
				os << m._matrix.At(i, j) << ", ";
			}
			os << m._matrix.At(i, M - 1);
			if (i == N - 1) {
				os << "\n]";
			}
			os << "\n  ";
		}
		return os;
	}

private:
	std::vector<T> _matrix;
	size_t _n, _m;
};
