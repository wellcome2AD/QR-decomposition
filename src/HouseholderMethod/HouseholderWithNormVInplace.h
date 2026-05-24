#pragma once

#include <vector>

#include "../IQRSolver.h"
#include "../Matrix/TransposedMatrix.h"

template <typename T>
class HouseholderWithNormVInplace : public IQRSolver<T> {
public:
	// третья версия, на основе второй. v вычисляется проще для уменьшения роста ошибки, нормированный v сохраняется в R
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

			// вычисление v, tau, sigma
			// sigma = -sign(Akk) * ||Ak||
			// v = Ak - sigma*(1,0...0)
			// tau = 2 / (vTv)

			double sum = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum)
			for (ptrdiff_t i = k; i < N; ++i) {
				sum += R_T.At(i, k) * R_T.At(i, k);
			}
			double norm = sqrt(sum);
			sigma = -sign(R_T.At(k, k)) * norm;

			// вычисляем v, храним без k нулей в начале
			std::vector<T> v(N - k);
			for (ptrdiff_t i = 1; i < N - k; ++i) {
				ptrdiff_t row = k + i;
				v[i] = R_T.At(row, k);
			}
			v[0] = R_T.At(k, k) - sigma;

			double denominator = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:denominator)
			for (ptrdiff_t i = 0; i < N - k; ++i) {
				denominator += v[i] * v[i];
			}

			if (denominator != 0) {
				tau = 2.0 / denominator;
			}

			// сохраняем нормированный v для Q (v[0]=1) в R под диагональю
			if (v[0] != 0) {
				double scale = 1.0 / v[0];
#pragma omp parallel for num_threads(thread_num)
				for (ptrdiff_t i = 1; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					R_T.At(row, k) = v[i] * scale;
				}
			}
			R_T.At(k, k) = sigma;

			// применяем отражения к правой части
#pragma omp parallel for num_threads(thread_num)
			for (ptrdiff_t col = k + 1; col < N; ++col) {
				double vTR = 0;
				// Вычисляем vTR = vT * столбец col
				for (ptrdiff_t i = 0; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					vTR += v[i] * R_T.At(row, col);
				}

				double scale = tau * vTR;

				// обновляем столбец col: R[:,col] -= v * scale
				for (ptrdiff_t i = 0; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					R_T.At(row, col) -= v[i] * scale;
				}
			}
		}

		// вычисление Q
		for (ptrdiff_t i = 0; i < N; ++i) {
			Q_T.At(i, i) = 1.0;
		}

		// применяем отражения в обратном порядке
		for (ptrdiff_t k = N - 2; k >= 0; --k) {
			// берём вектор v из R
			ptrdiff_t m = N - k;  // размер v
			std::vector<T> v(m);
			v[0] = 1.0;

			for (ptrdiff_t i = 1; i < m; ++i) {
				ptrdiff_t row = k + i;
				v[i] = R_T.At(row, k);
				R_T.At(row, k) = 0;  // очищаем
			}

			// tau = 2 / (vTv)
			double sum = 1.0; // v[0] = 1
#pragma omp parallel for num_threads(thread_num) reduction(+:sum)
			for (ptrdiff_t i = 1; i < m; ++i) {
				sum += v[i] * v[i];
			}
			double tau = 2.0 / sum;

			// H = I - tau * v * vT к Q[k:N, :]
#pragma omp parallel for num_threads(thread_num)
			for (ptrdiff_t j = 0; j < N; ++j) {
				// vTQ = vT * Q[k:N, col]
				double vTQ = Q_T.At(k, j); // v[0] = 1

				for (ptrdiff_t i = 1; i < m; ++i) {
					ptrdiff_t row = k + i;
					vTQ += v[i] * Q_T.At(row, j);
				}
				vTQ *= tau;

				// Q[k][j] -= vTQ * v
				Q_T.At(k, j) -= vTQ;  // v[0] = 1

				for (ptrdiff_t i = 1; i < m; ++i) {
					ptrdiff_t row = k + i;
					Q_T.At(row, j) -= vTQ * v[i];
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
