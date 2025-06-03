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
	return beta;
}

template <typename T>
static float count_mu(float beta, int k_step, const TMatrix<T>& Ak)
{
	auto mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * Ak[k_step][k_step]);
	return mu;
}

template <typename T>
static TVector<T> count_w(const TMatrix<T>& Ak, int k_step, float beta, float mu)
{
	const size_t N = Ak.Size();
	TVector<T> w(N);
	for (size_t str = 0; str < N; ++str)
	{ // нужен не весь столбец, а только начиная с k-ой строки
		w[str] = str < k_step ? 0 : Ak[str][k_step] * mu;
	}

	w[k_step] -= beta * mu;
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
static void count_Q(const TMatrix<T>& H, TMatrix<T>& Q) {
	bool isFirst = Q.IsZero();
	if (isFirst)
	{
		Q = H;
	}
	else
	{
		Q = Q * H;
	}

	return;
}

template <typename T>
void QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R)
{
	const size_t N = A.Size();
	R = A;
	bool isFirst = true;
	for (size_t k = 0; k < N - 1; ++k)
	{
		auto beta = count_beta(R, k);
		auto mu = count_mu(beta, k, R);
		auto w = count_w(R, k, beta, mu);

		auto H = count_H(beta, mu, w);
		count_Q(H, Q);

		R = H * R;
	}

	return;
}
