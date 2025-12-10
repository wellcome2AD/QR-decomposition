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
	TMatrix<T> all_w(N, N, 0);
	// вычисление R
	for (ptrdiff_t k = 0; k < N - 1; ++k) {
		double tau = 0.0, sigma = 0.0;
		TVector<T> w(N, 0);
		// вычисление w, tau, sigma
		// w = x - sigma*(1,0...0)
		// tau = 2 / (wTw)
		// sigma = -sign(x0) * ||x||
		[&R, k, N, &w, &tau, &sigma]() {
			double sum_sq = 0.0;
			for (ptrdiff_t str = k; str < N; ++str) {
				sum_sq += R[str][k] * R[str][k];
			}

			// sigma = -sign(x0) * ||x||
			double norm_x = sqrt(sum_sq);
			sigma = -sign(R[k][k]) * norm_x;

			// w = x - sigma*(1,0...0)
			w[k] = R[k][k] - sigma;
			for (ptrdiff_t str = k + 1; str < N; ++str) {
				w[str] = R[str][k];
			}

			// tau = 2 / (wTw)
			double denominator = 0;
			for (ptrdiff_t i = k; i < N; ++i) {
				denominator += w[i] * w[i];
			}

			if (denominator != 0) {
				tau = 2.0 / denominator;
			}
		}();

		R[k][k] = sigma;
		for (ptrdiff_t col = k + 1; col < N; ++col) {
			double dot = 0;
			for (ptrdiff_t row = k; row < N; ++row) {
				dot += w[row] * R[row][col];
			}

			double scale = tau * dot;
			for (ptrdiff_t row = k; row < N; ++row) {
				R[row][col] -= w[row] * scale;
			}
		}

		// зануляем столбец под диагональю
		for (ptrdiff_t i = k + 1; i < N; ++i) {
			R[i][k] = 0;
		}

		// сохраняем w для Q (масштабированную)
		if (w[k] != 0) {
			for (ptrdiff_t i = k; i < N; ++i) {
				all_w[i][k] = w[i] / w[k];
			}
		}
	}

	// вычисление Q
	Q = TMatrix<T>(N, N, 0);
	for (auto i = 0; i < N; i++) {
		Q[i][i] = 1;
	}
	// применяем отражения в обратном порядке
	for (ptrdiff_t k = N - 2; k >= 0; --k) {
		// берём вектор w из all_w
		int m = N - k;  // размер w
		TVector<T> w(m);
		w[0] = 1.0;
		for (ptrdiff_t i = 1; i < m; ++i) {
			w[i] = all_w[k + i][k];
		}

		// tau = 2 / (wTw)
		double w_norm_sq = 0;
		for (ptrdiff_t i = 0; i < m; ++i) {
			w_norm_sq += w[i] * w[i];
		}
		double tau = 2.0 / w_norm_sq;

		// H = I - tau * w * wT к Q[k:N, :]
		for (ptrdiff_t col = 0; col < N; ++col) {
			// dot = wT * Q[k:N, col]
			double dot = 0;
			for (ptrdiff_t i = 0; i < m; ++i) {
				dot += w[i] * Q[k + i][col];
			}
			dot *= tau;

			// Q[k:N, col] -= dot * w
			for (ptrdiff_t i = 0; i < m; ++i) {
				Q[k + i][col] -= dot * w[i];
			}
		}
	}
}
