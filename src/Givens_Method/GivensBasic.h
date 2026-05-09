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

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (int k = j; k < N; ++k) {
					T temp = static_cast<T>(R[j][k] * c - R[i][k] * s);
					R[i][k] = static_cast<T>(R[j][k] * s + R[i][k] * c);
					R[j][k] = temp;
				}

				double tau = s / (1.0 + c);
				R[i][j] = static_cast<T>(tau);
			}
		}

		for (int j = 0; j < N - 1; ++j)
		{
			for (int i = j + 1; i < N; ++i)
			{
				if (std::abs(R[i][j]) < 1e-11) continue;

				double tau = static_cast<double>(R[i][j]);
				double tau2 = tau * tau;
				double denom = 1.0 + tau2;
				double c = (1.0 - tau2) / denom;
				double s = 2.0 * tau / denom;

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (int k = 0; k < N; ++k)
				{
					T temp = static_cast<T>(c * Q[k][j] - s * Q[k][i]);
					Q[k][i] = static_cast<T>(s * Q[k][j] + c * Q[k][i]);
					Q[k][j] = temp;
				}
				R[i][j] = T(0);
			}
		}
	}
};
