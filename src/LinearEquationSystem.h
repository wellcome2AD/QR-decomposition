#pragma once

#include <cstdlib>
#include <ctime>

#include "TMatrix.h"

template <typename T>
void generate_linear_equation_system(TMatrix<T>& A, TVector<T>& b, TVector<T>& x) {
	const size_t N = A.Size();
	std::srand(std::time(0));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			A[i][j] = std::rand() % 20 - 10;
			x[i] = std::rand() % 20 - 10;
		}
	}

	b = substitution(A, x);
}

template <typename T>
TVector<T> substitution(const TMatrix<T>& A, const TVector<T>& x) {
	TVector<T> res = A * x;
	return res;
}


template <typename T>
TVector<T> back_substitution(const TMatrix<T>& Q, const TMatrix<T>& R, const TVector<T>& b)
{
	const size_t N = R.Size();
	TVector<T> res(N);

	auto Q_copy = Q;
	Q_copy.Transpone();

	auto Q_invert_b = Q_copy * b;

	res[N - 1] = Q_invert_b[N - 1] / R[N - 1][N - 1];
	for (ptrdiff_t i = N - 2; i >= 0; --i)
	{
		float sum = 0;
		for (size_t j = i + 1; j < N; ++j)
		{
			sum += R[i][j] * res[j];
		}
		res[i] = (Q_invert_b[i] - sum) / R[i][i];
	}

	return res;
}
