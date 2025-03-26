#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "TVector.h"

float sign(float x) {
	if (x > 0.0) return 1.0;
	if (x < 0.0) return -1.0;
	return 1.0;
}

fmatrix matr_diff(fmatrix m1, fmatrix m2) {
	auto N = m1.size();
	for (auto i = 0; i < N; ++i) {
		for (auto j = 0; j < N; ++j) {
			m1[i][j] -= m2[i][j];
		}
	}
	return m1;
}

template <typename T>
TVector<T> matr_mult(fmatrix m, TVector<T> v) {
	auto N = m.size();
	TVector<T> res(N, 0);
	for (auto i = 0; i < N; ++i) {
		for (auto j = 0; j < N; ++j) {
			res[i] += m[i][j] * v[j];
		}
	}
	return res;
}

fmatrix matr_mult(fmatrix m1, fmatrix m2) {
	auto N = m1.size();
	auto M = m2[0].size();
	fmatrix res(N, std::vector<float>(M, 0));
	for (auto i = 0; i < N; ++i) {
		for (auto j = 0; j < M; ++j) {
			for (auto k = 0; k < N; ++k) {
				res[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return res;
}

void print_matr(fmatrix m) {
	auto N = m.size();
	auto M = m[0].size();
	for (auto i = 0; i < N; ++i) {
		std::cout << std::left << std::setw(15);
		for (auto j = 0; j < M; ++j) {
			std::cout << m[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

fmatrix transpone(fmatrix mas) {
	auto rows = mas.size();
	auto cols = mas[0].size();
	fmatrix res(rows, fvector(cols, 0));
	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++) {
			res[i][j] = mas[j][i];
		}
	}
	return res;
}

template <typename T>
fmatrix count_H(float beta, float mu, TVector<T> w) {
	auto N = w.Size();
	fmatrix E(N, std::vector<float>(N, 0));
	for (auto i = 0; i < N; ++i) {
		E[i][i] = 1;
	}
	return matr_diff(E, vec_mult<T>(vec_mult<T>(w, 2.0), w));
}

template <typename T>
TVector<T> back_substitution(fmatrix Q, fmatrix A, TVector<T> b) {
	const auto N = A.size();
	TVector<T> res(N);

	fmatrix Q_invert = transpone(Q);
	printf("Q_invert\n");
	print_matr(Q_invert);
	std::cout << std::endl;

	TVector<T> Q_invert_b = matr_mult<T>(Q_invert, b);
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

int main() {
	const int N = 3;
	fmatrix A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} };
	const fmatrix A_copy = A;
	const TVector<float> b = { 1, 8, 5 };
	fmatrix Q;
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

		fmatrix H = count_H(beta, mu, w);
		printf("H%td\n", k + 1);
		print_matr(H);
		std::cout << std::endl;

		if (isFirst) {
			isFirst = false;
			Q = H;
		}
		else {
			Q = matr_mult(Q, H);
		}

		// умножить Hk на A, получим новую Ak
		A = matr_mult(H, A);
		printf("A%td\n", k + 1);
		print_matr(A);
		std::cout << std::endl;
	}

	printf("Q\n");
	print_matr(Q);
	std::cout << std::endl;

	TVector<float> x = back_substitution<float>(Q, A, b);

	printf("x\n");
	std::cout << x;

	// проверка
	const float eps = 0.000001;
	TVector<float> A_x = matr_mult<float>(A_copy, x);
	for (auto i = 0; i < A_x.Size(); ++i) {
		float bi = b[i];
		if (abs(bi - A_x[i]) >= eps) {
			printf("error in %td: %f != %f\n", i, A_x[i], b[i]);
		}
	}

	return 0;
}
