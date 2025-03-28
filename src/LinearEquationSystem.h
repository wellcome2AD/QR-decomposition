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
