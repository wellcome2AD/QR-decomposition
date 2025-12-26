#pragma once

#include <vector>

#include "IHouseholderMethodSolver.h"
#include "math.h"
#include "matrix.h"

template <typename T>
class HouseholderMethodBasic : public IHouseholderMethodSolver<T> {
	// первая, самая неоптимальная версия. на каждом шаге алгоритма выполняется два матричных умножения,
	// используется много лишней памяти
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		const ptrdiff_t N = A.size();
		R = A;
		bool isFirst = true;
		for (ptrdiff_t k = 0; k < N - 1; ++k) {
			double beta = count_beta(R, k);
			double mu = count_mu(beta, k, R);
			auto w = count_w(R, k, beta, mu);

			auto H = count_H(beta, mu, w);

			if (k == 0) {
				Q = H;
			}
			else {
				Q = multiplyMatrix(Q, H);
			}

			R = multiplyMatrix(H, R);
		}
	}

private:
	double count_beta(const std::vector<std::vector<T>>& Ak, int k_step)
	{
		const ptrdiff_t N = Ak.size();
		double sum_by_k_col = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum_by_k_col)
		for (int str = k_step; str < N; ++str) {
			sum_by_k_col += Ak[str][k_step] * Ak[str][k_step];
		}
		double beta = sign(-Ak[k_step][k_step]) * sqrt(sum_by_k_col);
		return beta;
	}

	double count_mu(double beta, int k_step, const std::vector<std::vector<T>>& Ak)
	{
		double mu = 1.0 / sqrt(2.0 * beta * beta - 2 * beta * double(Ak[k_step][k_step]));
		return mu;
	}

	std::vector<T> count_w(const std::vector<std::vector<T>>& Ak, int k_step, double beta, double mu)
	{
		const ptrdiff_t N = Ak.size();
		std::vector<T> w(N);
#pragma omp parallel for num_threads(thread_num)
		for (int str = 0; str < N; ++str) {
			w[str] = str < k_step ? 0 : Ak[str][k_step] * mu;
		}
		w[k_step] -= beta * mu;
		return w;
	}

	 std::vector<std::vector<T>> count_H(double beta, double mu, std::vector<T> w)
	{
		ptrdiff_t N = w.size();
		 std::vector<std::vector<T>> H(N, std::vector<T>(N, 0));
#pragma omp parallel for num_threads(thread_num)
		for (int i = 0; i < N; ++i) {
			H[i][i] = 1.0;
			for (int j = 0; j < N; ++j) {
				H[i][j] -= w[i] * w[j] * 2.0;
			}
		}
		return H;
	}
};
