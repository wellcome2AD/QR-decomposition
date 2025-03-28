#pragma once

#include "TMatrix.h"

float sign(float x)
{
	if (x > 0.0)
		return 1.0;
	if (x < 0.0)
		return -1.0;
	return 1.0;
}

template <typename T>
float count_beta(const TMatrix &A, int k_step)
{
	const auto N = A.size();
	float sum_by_k_col = 0.0;
	for (auto str = k_step; str < N; ++str)
	{ // нужен не весь столбец, а только начиная с k-ой строки
		sum_by_k_col += A[str][k_step] * A[str][k_step];
	}
	auto beta = sign(-A[k_step][k_step]) * sqrt(sum_by_k_col);
	printf("beta%ld = %e\n", k_step, beta);
	return beta;
}

float count_mu(float beta, int k_step)
{
	auto mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * A[k_step][k_step]);
	printf("mu%ld = %e\n", k_step, mu);
	return mu;
}

TVector<T> count_w(const TMatrix &A, int k_step, float beta, float mu)
{
	for (auto str = k_step; str < N; ++str)
	{ // нужен не весь столбец, а только начиная с k-ой строки
		w[str] = A[str][k_step];
	}

	w[k_step] -= beta;
	w = w * mu;
	for (auto index = 0; index < k_step; ++index)
	{ // первые k позиций - нулевые
		w[index] = 0;
	}
	
	return w;
}

template <typename T>
TMatrix<T> count_H(float beta, float mu, TVector<T> w)
{
	auto N = w.Size();
	TMatrix<T> E(N, N);
	for (auto i = 0; i < N; ++i)
	{
		E[i][i] = 1.0;
	}
	auto wt = w;
	auto H = E - w * wt * 2.0;
	return H;
}

TMatrix<T>& count_Q(const TMatrix<T>& H, const TMatrix<T>& Q) {
	if (isFirst)
		{
			isFirst = false;
			Q = H;
		}
		else
		{
			Q = Q * H;
		}

	return Q;
}

template <typename T>
TVector<T> back_substitution(const TMatrix<T>& Q, const TMatrix<T>& A, const TVector<T>& b)
{
	const auto N = A.size();
	TVector<T> res(N);

	Q.Transpone();
	std::cout << "Q_invert\n"
			  << Q << std::endl;

	Q_invert = Q_invert * b;
	printf("Q_invert_b\n");
	std::cout << Q_invert_b << std::endl;
	std::cout << std::endl;

	res[N - 1] = Q_invert_b[N - 1] / A[N - 1][N - 1];
	for (auto i = N - 2; i >= 0; --i)
	{
		float sum = 0;
		for (auto j = i + 1; j < N; ++j)
		{
			sum += A[i][j] * res[j];
		}
		res[i] = (Q_invert_b[i] - sum) / A[i][i];
	}

	return res;
}

template <typename T>
TVector<T> QR_decomposition(const TMatrix<T> &A, const TVector<T> &b)
{
	const auto N = A.Size();
	const auto A_copy = A;
	TMatrix<float> Q;
	bool isFirst = true;
	for (auto k = 0; k < N - 1; ++k)
	{ // k-ый шаг алгоритма
		auto beta = count_bet();
		auto mu = count_mu();
		auto w = count_w();
		auto H = count_H(beta, mu, w);
		printf("H%td\n", k_step + 1);
		std::cout << H << std::endl;
		Q = count_Q();

		// умножить Hk на A, получим новую Ak
		A = H * A;
		printf("A%td\n", k + 1);
		print_matr(A);
		std::cout << std::endl;
	}

	std::cout << "Q\n" << Q << std::endl;

	return back_substitution<float>(Q, A, b);
}
