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
		std::vector<T> R_T(N * N, 0.0);
		std::vector<T> Q_T(N * N, 0.0);

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for(int i = 0; i < N; i++)
		{
			Q_T[i * N + i] = 1.0;
			for (int j = 0; j < N; j++)
			{
				R_T[j * N + i] = A[i][j];
			}
		}

		std::vector<T> v(N);
		// вычисление R
		for (ptrdiff_t k = 0; k < N - 1; ++k) {
			double tau = 0.0, sigma = 0.0;

			// вычисление v, tau, sigma
			// sigma = -sign(Akk) * ||Ak||
			// v = Ak - sigma*(1,0...0)
			// tau = 2 / (vTv)

			double sum = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum) if (N >= 1000)
			for (ptrdiff_t i = k; i < N; ++i) {
				sum += R_T[k * N + i] * R_T[k * N + i];
			}
			double norm = sqrt(sum);
			sigma = -sign(R_T[k * N + k]) * norm;

			// вычисляем v, храним без k нулей в начале
			for (ptrdiff_t i = 1; i < N - k; ++i) {
				ptrdiff_t row = k + i;
				v[i] = R_T[k * N + row];
			}
			v[0] = R_T[k * N + k] - sigma;

			double denominator = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:denominator) if (N >= 1000)
			for (ptrdiff_t i = 0; i < N - k; ++i) {
				denominator += v[i] * v[i];
			}

			if (denominator != 0) {
				tau = 2.0 / denominator;
			}

			// сохраняем нормированный v для Q (v[0]=1) в R под диагональю
			if (v[0] != 0) {
				double scale = 1.0 / v[0];
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (ptrdiff_t i = 1; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					R_T[k * N + row] = v[i] * scale;
				}
			}
			R_T[k * N + k] = sigma;

			// применяем отражения к правой части
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
			for (ptrdiff_t col = k + 1; col < N; ++col) {
				double vTR = 0;
				// Вычисляем vTR = vT * столбец col
				for (ptrdiff_t i = 0; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					vTR += v[i] * R_T[col * N + row];
				}

				double scale = tau * vTR;

				// обновляем столбец col: R[:,col] -= v * scale
				for (ptrdiff_t i = 0; i < N - k; ++i) {
					ptrdiff_t row = k + i;
					R_T[col * N + row] -= v[i] * scale;
				}
			}
		}

		// вычисление Q
		// применяем отражения в обратном порядке
		for (ptrdiff_t k = N - 2; k >= 0; --k) {
			// берём вектор v из R
			ptrdiff_t m = N - k;  // размер v
			std::vector<T> v(m);
			v[0] = 1.0;

			for (ptrdiff_t i = 1; i < m; ++i) {
				ptrdiff_t row = k + i;
				v[i] = R_T[k * N + row];
				       R_T[k * N + row] = 0;  // очищаем
			}

			// tau = 2 / (vTv)
			double sum = 1.0; // v[0] = 1
#pragma omp parallel for num_threads(thread_num) reduction(+:sum) if (N >= 1000)
			for (ptrdiff_t i = 1; i < m; ++i) {
				sum += v[i] * v[i];
			}
			double tau = 2.0 / sum;

			// H = I - tau * v * vT к Q[k:N, :]
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
			for (ptrdiff_t j = 0; j < N; ++j) {
				// vTQ = vT * Q[k:N, col]
				double vTQ = Q_T[j * N + k]; // v[0] = 1

				for (ptrdiff_t i = 1; i < m; ++i) {
					ptrdiff_t row = k + i;
					vTQ += v[i] * Q_T[j * N + row];
				}
				vTQ *= tau;

				// Q[k][j] -= vTQ * v
				Q_T[j * N + k] -= vTQ;  // v[0] = 1

				for (ptrdiff_t i = 1; i < m; ++i) {
					ptrdiff_t row = k + i;
					Q_T[j * N + row] -= vTQ * v[i];
				}
			}
		}

		// перевести в обычную матрицу
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				R[i][j] = R_T[j * N + i];
				Q[i][j] = Q_T[j * N + i];
			}
		}
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
