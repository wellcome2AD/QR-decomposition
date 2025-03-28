#pragma once

#include "TMatrix.h"

inline static float sign(float x)
{
	if (x > 0.0)
		return 1.0;
	if (x < 0.0)
		return -1.0;
	return 1.0;
}

template <typename T>
static float count_beta(const TMatrix<T>& Ak, int k_step)
{
	const size_t N = Ak.Size();
	float sum_by_k_col = 0.0;
	for (size_t str = k_step; str < N; ++str)
	{ // нужен не весь столбец, а только начиная с k-ой строки
		sum_by_k_col += Ak[str][k_step] * Ak[str][k_step];
	}
	auto beta = sign(-Ak[k_step][k_step]) * sqrt(sum_by_k_col);
	// printf("beta%ld = %e\n", k_step, beta);
	return beta;
}

template <typename T>
static float count_mu(float beta, int k_step, const TMatrix<T>& Ak)
{
	auto mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * Ak[k_step][k_step]);
	// printf("mu%ld = %e\n", k_step, mu);
	return mu;
}

template <typename T>
static TVector<T> count_w(const TMatrix<T>& Ak, int k_step, float beta, float mu)
{
	const size_t N = Ak.Size();
	TVector<T> w(N, 0);
	for (size_t str = k_step; str < N; ++str)
	{ // нужен не весь столбец, а только начиная с k-ой строки
		w[str] = Ak[str][k_step];
	}

	w[k_step] -= beta;
	w = w * mu;

	return w;
}

template <typename T>
static TMatrix<T> count_H(float beta, float mu, TVector<T> w)
{
	size_t N = w.Size();
	TMatrix<T> E(N, N);
	for (size_t i = 0; i < N; ++i)
	{
		E[i][i] = 1.0;
	}
	auto H = E - w * w * 2.0;
	return H;
}

template <typename T>
static TMatrix<T>& count_Q(const TMatrix<T>& H, TMatrix<T>& Q) {
	bool isFirst = Q.IsZero();
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
static TVector<T> back_substitution(const TMatrix<T>& Q, const TMatrix<T>& AH, const TVector<T>& b)
{
	const size_t N = AH.Size();
	TVector<T> res(N);

	auto Q_copy = Q;
	Q_copy.Transpone();
	/*std::cout << "Q_invert\n"
		<< Q_copy << std::endl;*/

	auto Q_invert_b = Q_copy * b;
	/*printf("Q_invert_b\n");
	std::cout << Q_invert_b << std::endl;*/

	res[N - 1] = Q_invert_b[N - 1] / AH[N - 1][N - 1];
	for (ptrdiff_t i = N - 2; i >= 0; --i)
	{
		float sum = 0;
		for (size_t j = i + 1; j < N; ++j)
		{
			sum += AH[i][j] * res[j];
		}
		res[i] = (Q_invert_b[i] - sum) / AH[i][i];
	}

	return res;
}

template <typename T>
TVector<T> QR_decomposition(const TMatrix<T>& A, const TVector<T>& b)
{
	const size_t N = A.Size();
	auto A_copy = A;
	TMatrix<float> Q;
	bool isFirst = true;
	for (size_t k = 0; k < N - 1; ++k)
	{ // k-ый шаг алгоритма
		auto beta = count_beta(A_copy, k);
		auto mu = count_mu(beta, k, A_copy);
		auto w = count_w(A_copy, k, beta, mu);

		auto H = count_H(beta, mu, w);
		/*printf("H%td\n", k + 1);
		std::cout << H << std::endl;*/

		Q = count_Q(H, Q);

		// умножить Hk на A, получим новую Ak
		A_copy = H * A_copy;
		/*printf("A%td\n", k + 1);
		std::cout << A_copy << std::endl;*/
	}

	// std::cout << "Q\n" << Q << std::endl;

	return back_substitution<float>(Q, A_copy, b);
}
