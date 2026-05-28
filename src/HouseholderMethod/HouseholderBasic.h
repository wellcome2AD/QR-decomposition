#pragma once

#include <vector>

#include "../IQRSolver.h"
#include "Matrix/MatrixOperations.h"

template <typename T>
class HouseholderMethodBasic : public IQRSolver<T> {
	// первая, самая неоптимальная версия. на каждом шаге алгоритма выполняется два матричных умножения,
	// используется много лишней памяти
public:
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		const ptrdiff_t N = A.size();
		R = A;
		Q = std::vector<std::vector<T>>(N, std::vector<T>(N, 0));
		for (int i = 0; i < N; ++i) {
			Q[i][i] = 1.0;
		}

		for (ptrdiff_t k = 0; k < N - 1; ++k) {
			// count beta
			double sum = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum)
			for (ptrdiff_t i = k; i < N; ++i) {
				sum += R[i][k] * R[i][k];
			}
			double beta = sign(-R[k][k]) * sqrt(sum);

			// count mu
			double mu = 1.0 / sqrt(2.0 * beta * (beta - R[k][k]));

			// count v
			std::vector<double> v(N, 0);
#pragma omp parallel for num_threads(thread_num)
			for (int i = k; i < N; ++i) {
				v[i] = R[i][k] * mu;
			}
			v[k] -= beta * mu;

			// count H
			std::vector<std::vector<T>> H(N, std::vector<T>(N, 0));
#pragma omp parallel for num_threads(thread_num)
			for (int i = 0; i < N; ++i) {
				H[i][i] = 1.0;
				for (int j = 0; j < N; ++j) {
					H[i][j] -= v[i] * v[j] * 2.0;
				}
			}

			R = multiplyMatrix(H, R);
			Q = multiplyMatrix(Q, H);
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
