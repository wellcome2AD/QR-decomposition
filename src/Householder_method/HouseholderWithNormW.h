#pragma once

#include "../TMatrix.h"
#include "../math.h"

template <typename T>
class HouseholderMethodWithNormW : public IHouseholderMethodSolver<T> {
	// третья версия, на основе второй. w вычисляется проще для уменьшения роста ошибки, нормированный w сохраняется в R
	virtual void QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R) override
	{
		const ptrdiff_t N = A.Size();
		R = A;
		// вычисление R
		for (ptrdiff_t k = 0; k < N - 1; ++k) {
			double tau = 0.0, sigma = 0.0;
			TVector<T> w(N, 0);
			// вычисление w, tau, sigma
			// w = x - sigma*(1,0...0)
			// tau = 2 / (wTw)
			// sigma = -sign(x0) * ||x||
			double sum_sq = 0.0;
#pragma omp parallel for reduction(+:sum_sq)
			for (ptrdiff_t str = k; str < N; ++str) {
				sum_sq += R[str][k] * R[str][k];
			}

			double norm_x = sqrt(sum_sq);
			sigma = -sign(R[k][k]) * norm_x;

			w[k] = R[k][k] - sigma;
#pragma omp parallel for // if(N-k-1 > 1000)
			for (ptrdiff_t str = k + 1; str < N; ++str) {
				w[str] = R[str][k];
			}

			double denominator = 0;
#pragma omp parallel for reduction(+:denominator)
			for (ptrdiff_t i = k; i < N; ++i) {
				denominator += w[i] * w[i];
			}

			if (denominator != 0) {
				tau = 2.0 / denominator;
			}

			// сохраняем w для Q (нормированный, чтобы w[0]=1), теперь можно хранить в R под диагональю
			if (w[k] != 0) {
				double scale = 1.0 / w[k];
#pragma omp parallel for
				for (ptrdiff_t i = k + 1; i < N; ++i) {
					R[i][k] = w[i] * scale;
				}
			}

			R[k][k] = sigma;

			// применяем отражения к правой части
#pragma omp parallel for
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
		}

		// вычисление Q
		Q = TMatrix<T>(N, N, 0);
#pragma omp parallel for
		for (ptrdiff_t i = 0; i < N; ++i) {
			Q[i][i] = 1;
		}

		// применяем отражения в обратном порядке
		for (ptrdiff_t k = N - 2; k >= 0; --k) {
			// берём вектор w из R
			int m = N - k;  // размер w
			TVector<T> w(m);
			w[0] = 1.0;
			for (ptrdiff_t i = 1; i < m; ++i) {
				w[i] = R[k + i][k];
				R[k + i][k] = 0;
			}

			// tau = 2 / (wTw)
			double w_norm_sq = 1.0; // w[0] = 1
			for (ptrdiff_t i = 1; i < m; ++i) {
				w_norm_sq += w[i] * w[i];
			}
			double tau = 2.0 / w_norm_sq;

			// H = I - tau * w * wT к Q[k:N, :]
#pragma omp parallel for
			for (ptrdiff_t col = 0; col < N; ++col) {
				// dot = wT * Q[k:N, col]
				double dot = Q[k][col]; // w[0] = 1
				for (ptrdiff_t i = 1; i < m; ++i) {
					dot += w[i] * Q[k + i][col];
				}
				dot *= tau;

				// Q[k:N, col] -= dot * w
				Q[k][col] -= dot;  // w[0] = 1
				for (ptrdiff_t i = 1; i < m; ++i) {
					Q[k + i][col] -= dot * w[i];
				}
			}
		}
	}
};
