#pragma once

#include "TMatrix.h"

float sign(float x) {
	if (x > 0.0) return 1.0;
	if (x < 0.0) return -1.0;
	return 1.0;
}

template <typename T>
TMatrix<T> count_H(float beta, float mu, TVector<T> w) {
	auto N = w.Size();
	TMatrix<T> E(N, N);
	for (auto i = 0; i < N; ++i) {
		E[i][i] = 1.0;
	}
	auto wt = w;
	return E - w * wt * 2.0;
}

template <typename T>
TVector<T> back_substitution(TMatrix<T> Q, TMatrix<T> A, TVector<T> b) {
	const auto N = A.size();
	TVector<T> res(N);

	Q.Transpone();
	std::cout << "Q_invert\n" << Q << std::endl;

	Q_invert  = Q_invert * b;
	printf("Q_invert_b\n");
	std::cout << Q_invert_b << std::endl;
	std::cout << std::endl;

	res[N - 1] = Q_invert_b[N - 1] / A[N - 1][N - 1];
	for (auto i = N - 2; i >= 0; --i) {
		float sum = 0;
		for (auto j = i + 1; j < N; ++j) {
			sum += A[i][j] * res[j];
		}
		res[i] = (Q_invert_b[i] - sum) / A[i][i];
	}

	return res;
}

template <typename T>
TVector<T> QR_decomposition(const TMatrix<T>& A, const TVector<T>& b) {
	const auto N = A.Size();
	const auto A_copy = A;
	TMatrix<float> Q;
	bool isFirst = true;
	for (auto k = 0; k < N - 1; ++k) { // k-ый шаг алгоритма
		TVector<float> w(N);
		float sum_by_k_col = 0.0;
		for (auto str = k; str < N; ++str) { // нужен не весь столбец, а только начиная с k-ой строки
			w[str] = A[str][k];
			sum_by_k_col += A[str][k] * A[str][k];
		}

		auto beta = sign(-A[k][k]) * sqrt(sum_by_k_col);
		printf("beta%ld = %e\n", k, beta);

		auto mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * A[k][k]);
		printf("mu%ld = %e\n", k, mu);

		w[k] -= beta;
		w = w * mu;
		for (auto index = 0; index < k; ++index) { // первые k позиций - нулевые
			w[index] = 0;
		}

		auto H = count_H(beta, mu, w);
		printf("H%td\n", k + 1);
		std::cout << H << std::endl;

		if (isFirst) {
			isFirst = false;
			Q = H;
		}
		else {
			Q = Q * H;
		}

		// умножить Hk на A, получим новую Ak
		A = H * A;
		printf("A%td\n", k + 1);
		print_matr(A);
		std::cout << std::endl;
	}

	std::cout << "Q\n" << Q << std::endl;

	return back_substitution<float>(Q, A, b);
}
