#pragma once

#include <vector>

#include "../IQRSolver.h"
#include "../Matrix/TransposedMatrix.h"

template <typename T>
class HouseholderWithNormWInplace : public IQRSolver<T> {
public:
	// третья версия, на основе второй. w вычисляется проще для уменьшения роста ошибки, нормированный w сохраняется в R
	virtual void QR_decomposition(
		const std::vector<std::vector<T>>& A,
		std::vector<std::vector<T>>& Q,
		std::vector<std::vector<T>>& R) override
	{
		const size_t N = A.size();
		TransposedMatrix<T> A_T(A), R_T(A);
		TransposedMatrix<T> Q_T(N, N, 0.0);

		// вычисление R
		for (ptrdiff_t k = 0; k < N - 1; ++k) {
			double tau = 0.0, sigma = 0.0;

			// вычисление w, tau, sigma
			// sigma = -sign(Akk) * ||Ak||
			// w = Ak - sigma*(1,0...0)
			// tau = 2 / (wTw)

			double sum_sq = 0.0;
			for (ptrdiff_t row = k; row < N; ++row) {
				T val = R_T.At(row, k);
				sum_sq += static_cast<double>(val) * static_cast<double>(val);
			}

			double norm = sqrt(sum_sq);
			sigma = -sign(R_T.At(k, k)) * norm;

			// вычисляем w, храним без k нулей в начале
			std::vector<T> w(N - k);
			for (ptrdiff_t i = 1; i < N - k; ++i) {
				ptrdiff_t row = k + i;
				w[i] = R_T.At(row, k);
			}
			w[0] = R_T.At(k, k) - sigma;

			double denominator = w[0] * w[0];
			for (ptrdiff_t i = 1; i < N - k; ++i) {
				denominator += w[i] * w[i];
			}

			if (denominator != 0) {
				tau = 2.0 / denominator;
			}

			// сохраняем нормированный w для Q (w[0]=1) в R под диагональю
			if (w[0] != 0) {
				double scale = 1.0 / w[0];
				for (ptrdiff_t i = 1; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					R_T.At(row, k) = w[i] * scale;
				}
			}
			R_T.At(k, k) = sigma;

			// применяем отражения к правой части
#pragma omp parallel for if(N >= 1000)
			for (ptrdiff_t col = k + 1; col < N; ++col) {
				double wTR = 0;
				// Вычисляем wTR = wT * столбец col
				for (ptrdiff_t i = 0; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					wTR += w[i] * R_T.At(row, col);
				}

				double scale = tau * wTR;

				// обновляем столбец col: R[:,col] -= w * scale
				for (ptrdiff_t i = 0; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					R_T.At(row, col) -= w[i] * scale;
				}
			}
		}

		// вычисление Q
		for (ptrdiff_t i = 0; i < N; ++i) {
			Q_T.At(i, i) = 1.0;
		}

		// применяем отражения в обратном порядке
		for (ptrdiff_t k = N - 2; k >= 0; --k) {
			// берём вектор w из R
			ptrdiff_t m = N - k;  // размер w
			std::vector<T> w(m);
			w[0] = 1.0;

			for (ptrdiff_t i = 1; i < m; ++i) {
				ptrdiff_t row = k + i;
				w[i] = R_T.At(row, k);
				R_T.At(row, k) = 0;  // очищаем
			}

			// tau = 2 / (wTw)
			double w_norm_sq = 1.0; // w[0] = 1
			for (ptrdiff_t i = 1; i < m; ++i) {
				w_norm_sq += w[i] * w[i];
			}
			double tau = 2.0 / w_norm_sq;

			// H = I - tau * w * wT к Q[k:N, :]
#pragma omp parallel for if(N >= 1000)
			for (ptrdiff_t col = 0; col < N; ++col) {
				// wTQ = wT * Q[k:N, col]
				double wTQ = Q_T.At(k, col); // w[0] = 1

				for (ptrdiff_t i = 1; i < m; ++i) {
					ptrdiff_t row = k + i;
					wTQ += w[i] * Q_T.At(row, col);
				}
				wTQ *= tau;

				// Q[k:N, col] -= wTQ * w
				Q_T.At(k, col) -= wTQ;  // w[0] = 1

				for (ptrdiff_t i = 1; i < m; ++i) {
					ptrdiff_t row = k + i;
					Q_T.At(row, col) -= wTQ * w[i];
				}
			}
		}

		// перевести в обычную матрицу
		Q = Q_T.Transpose();
		R = R_T.Transpose();
	}
private:
	static double sign(double x)
	{
		if (x > 0.0)
			return 1.0;
		if (x < 0.0)
			return -1.0;
		return 1.0;
	}
};
