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
static double count_beta(const TMatrix<T>& Ak, int k_step)
{
	const ptrdiff_t N = Ak.Size();
	double sum_by_k_col = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum_by_k_col)
	for (int str = k_step; str < N; ++str) {
		sum_by_k_col += Ak[str][k_step] * Ak[str][k_step];
	}
	auto beta = sign(-Ak[k_step][k_step]) * sqrt(sum_by_k_col);
	return beta;
}

template <typename T>
static double count_mu(float beta, int k_step, const TMatrix<T>& Ak)
{
	double mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * Ak[k_step][k_step]);
	return mu;
}

template <typename T>
static TVector<T> count_w(const TMatrix<T>& Ak, int k_step, float beta, float mu)
{
	const ptrdiff_t N = Ak.Size();
	TVector<T> w(N);
#pragma omp parallel for num_threads(thread_num)
	for (int str = 0; str < N; ++str) {
		w[str] = str < k_step ? 0 : Ak[str][k_step] * mu;
	}
	w[k_step] -= beta * mu;
	return w;
}

template <typename T>
static TMatrix<T> count_H(float beta, float mu, TVector<T> w)
{
	ptrdiff_t N = w.Size();
	TMatrix<T> H(N, N, 0);
#pragma omp parallel for num_threads(thread_num)
	for (int i = 0; i < N; ++i) {
		H[i][i] = 1.0;
		for (int j = 0; j < N; ++j) {
			H[i][j] -= w[i] * w[j] * 2.0;
		}
	}
	return H;
}

template <typename T>
void Household_QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R)
{
	const ptrdiff_t N = A.Size();
	R = A;
	bool isFirst = true;
	for (ptrdiff_t k = 0; k < N - 1; ++k) {
		auto beta = count_beta(R, k);
		auto mu = count_mu(beta, k, R);
		auto w = count_w(R, k, beta, mu);

		auto H = count_H(beta, mu, w);

		if (k == 0) {
			Q = H;
		}
		else {
			Q = Q * H;
		}

		R = H * R;
	}

	return;
}

template <typename T>
void Household_QR_decomposition_experimental(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R)
{
	const ptrdiff_t N = A.Size();
	R = A;
	TMatrix<T> all_w(N, N);
	// вычисление R
	for (ptrdiff_t k = 0; k < N - 1; ++k) {
		auto beta = count_beta(R, k);
		auto mu = count_mu(beta, k, R);
		auto w = count_w(R, k, beta, mu);

		// H * R = (I - 2 * w * wᵀ) * R = R - 2 * w * (wT * R)
		// формируем (wT * R)
		TMatrix<T> y(1, N, 0);
		for (ptrdiff_t i = 0; i < N; ++i)
		{
			for (ptrdiff_t j = 0; j < N; ++j)
			{
				y[0][i] += w[j] * R[j][i];
			}
		}

		// y = 2 * y
		for (ptrdiff_t i = 0; i < N; ++i)
		{
			y[0][i] *= 2.0;
		}

		// w * y = матрица размером n×n, где (i,j)-й элемент = w[i] * y[j]
		// w это вектор nx1, y это вектор 1xn
		TMatrix<T> wy(N, N);
		for (ptrdiff_t i = 0; i < N; ++i)
		{
			for (ptrdiff_t j = 0; j < N; ++j)
			{
				wy[i][j] = w[i] * y[0][j];
			}
		}

		// R = R - w * y   (w * y — внешнее произведение, outer product)
		R = R - wy;

		// сохраним в столбце k вектор wk с элемента k+1
		for (auto i = 0; i < N; ++i)
		{
			all_w[i][k] = w[i];
		}
	}

	Q = TMatrix<T>(N, N, 0);
	for (auto i = 0; i < N; i++)
	{
		Q[i][i] = 1;
	}
	// теперь для вычисления Q достанем вектора wk
	for (ptrdiff_t k = N - 2; k >= 0; --k) {
		double sum = 0.0;
		for (ptrdiff_t i = 0; i < N; ++i) {
			sum += all_w[i][k] * all_w[i][k];
		}

		// Q = H_{N-2} * Q = Q - 2 * w * (wT * Q)
		TMatrix<T> wTQ(1, N, 0);
		for (ptrdiff_t i = k; i < N; ++i)
		{
			for (ptrdiff_t j = k; j < N; ++j) // todo: нужно оптимизировать, умножать не полностью; видимо, с k-го элемента
			{
				wTQ[0][i] += all_w[j][k] * Q[j][i];
			}
		}

		// w = 2 * w
		for (ptrdiff_t i = k; i < N; ++i)
		{
			all_w[i][k] *= 2;
		}

		// w * wTQ = матрица размером n×n, где (i,j)-й элемент = w[i] * y[j]
		// w это вектор nx1, wTQ это вектор 1xn
		TMatrix<T> wwTQ(N, N);
		for (ptrdiff_t i = 0; i < N; ++i)
		{
			for (ptrdiff_t j = 0; j < N; ++j)
			{
				wwTQ[i][j] = all_w[i][k] * wTQ[0][j];
			}
		}

		// теперь вычитаем из Q полученную матрицу wwTQ
		for (ptrdiff_t i = 0; i < N; ++i)
		{
			for (ptrdiff_t j = 0; j < N; ++j)
			{
				Q[i][j] -= wwTQ[i][j];
			}
		}
	}
}
