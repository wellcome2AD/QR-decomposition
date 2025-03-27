#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

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
	TVector<T> wt = w;
	return E - (w * wt * 2.0);
}

template <typename T>
TVector<T> back_substitution(TMatrix<T> Q, TMatrix<T> A, TVector<T> b) {
	const int N = A.Size();
	TVector<T> res(N);

	Q.Transpone();
	std::cout << "Q_invert\n" << Q << std::endl;

	auto Q_invert_b = Q * b;
	printf("Q_invert_b\n");
	std::cout << Q_invert_b << std::endl;

	res[N - 1] = Q_invert_b[N - 1] / A[N - 1][N - 1];
	for (ptrdiff_t i = N - 2; i >= 0; --i) {
		float sum = 0;
		for (auto j = i + 1; j < N; ++j) {
			sum += A[i][j] * res[j];
		}
		res[i] = (Q_invert_b[i] - sum) / A[i][i];
	}

	return res;
}

int main() {
	const int N = 3;
	TMatrix<float> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} };
	const TMatrix<float> A_copy = A;
	const TVector<float> b = { 1, 8, 5 };
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
		std::cout << A << std::endl;
	}

	std::cout << "Q\n" << Q << std::endl;

	TVector<float> x = back_substitution<float>(Q, A, b);

	printf("x\n");
	std::cout << x;

	// проверка
	const float eps = 0.000001;
	TVector<float> A_x = A_copy * x;
	for (auto i = 0; i < A_x.Size(); ++i) {
		float bi = b[i];
		if (abs(bi - A_x[i]) >= eps) {
			printf("error in %I32d: %f != %f\n", i, A_x[i], b[i]);
		}
	}

	return 0;
}
