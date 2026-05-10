#pragma once

#include <vector>
#include <cmath>

#include "../IQRSolver.h"
#include "math.h"

template <typename T>
class GivensMethodBasic : public IQRSolver<T> {
public:
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A,
		std::vector<std::vector<T>>& Q,
		std::vector<std::vector<T>>& R) override
	{
		auto N = A.size();
		// матрица Q будет хранится транспонированной
		Q = std::vector<std::vector<T>>(N, std::vector<T>(N, 0));
		R = A;

		for (int i = 0; i < N; ++i)
			Q[i][i] = T(1);

		for (int j = 0; j < N - 1; ++j) {
			for (int i = j + 1; i < N; ++i) {
				// приведение к double для точных проверок и вычислений
				double Rjj = static_cast<double>(R[j][j]);
				double Rij = static_cast<double>(R[i][j]);

				// проверки малости (пороги для double)
				if (std::abs(Rij) < 1e-11)
					continue;

				double sqrt_val = std::sqrt(Rjj * Rjj + Rij * Rij);
				if (std::abs(sqrt_val) < 1e-11)
					continue;

				double c = Rjj / sqrt_val;
				double s = -Rij / sqrt_val;

				// вычисление матрицы R путём вращения строк
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (int k = j; k < N; ++k) {
					T temp = static_cast<T>(R[j][k] * c - R[i][k] * s);
					R[i][k] = static_cast<T>(R[j][k] * s + R[i][k] * c);
					R[j][k] = temp;
				}

				//std::cout << "expected R, rotation (" << i << "," << j << "):";
				//printMatrix(R);
				//std::cout << std::endl;

				// сохранение коэффициентов вращения
				double tau = s / (1.0 + c);
				R[i][j] = static_cast<T>(tau);
			}
		}

		for (int j = 0; j < N - 1; ++j)
		{
			for (int i = j + 1; i < N; ++i)
			{
				if (std::abs(R[i][j]) < 1e-11) continue;
				// восстановление коэффициентов вращения
				double tau = static_cast<double>(R[i][j]);
				double tau2 = tau * tau;
				double denom = 1.0 + tau2;
				double c = (1.0 - tau2) / denom;
				double s = 2.0 * tau / denom;

				// вычисление матрицы Q путём вращения столбцов
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (int k = 0; k < N; ++k)
				{// обращение к элементам транспонированной матрицы Q
					T temp = static_cast<T>(c * Q[j][k] - s * Q[i][k]);
					Q[i][k] = static_cast<T>(s * Q[j][k] + c * Q[i][k]);
					Q[j][k] = temp;
				}
				R[i][j] = T(0);
			}
		}
		// транспонирование
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < i; ++j)
				std::swap(Q[i][j], Q[j][i]);
		}
	}
};
