#pragma once

#include "..\TMatrix.h"

template <typename T>
class HouseholderMethodWithoutMatrixMults : public IHouseholderMethodSolver<T> {
	// вторая версия, избавлена от недостатка первой -- двух матричных умножений на каждой итерации
	virtual void QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R) override
	{
		const ptrdiff_t N = A.Size();
		R = A;
		TMatrix<T> all_w(N, N);
		// вычисление R
		for (ptrdiff_t k = 0; k < N - 1; ++k) {
			double beta = count_beta(R, k);
			double mu = count_mu(beta, k, R);
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
			for (ptrdiff_t i = 0; i < N; ++i)
			{
				all_w[i][k] = w[i];
			}
		}

		Q = TMatrix<T>(N, N, 0);
		for (ptrdiff_t i = 0; i < N; i++)
		{
			Q[i][i] = 1;
		}
		// теперь для вычисления Q достанем вектора wk
		for (ptrdiff_t k = N - 2; k >= 0; --k) {
			double sum = 0.0;
			for (ptrdiff_t i = 0; i < N; ++i) {
				sum += all_w[i][k] * all_w[i][k];
			}

			// Q = H * Q = Q - 2 * w * (wT * Q)
			TMatrix<T> wTQ(1, N, 0);
			for (ptrdiff_t i = k; i < N; ++i)
			{
				for (ptrdiff_t j = k; j < N; ++j)
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

private:
	double count_beta(const TMatrix<T>& Ak, int k_step)
	{
		const ptrdiff_t N = Ak.Size();
		double sum_by_k_col = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum_by_k_col)
		for (int str = k_step; str < N; ++str) {
			sum_by_k_col += Ak[str][k_step] * Ak[str][k_step];
		}
		double beta = sign(-Ak[k_step][k_step]) * sqrt(sum_by_k_col);
		return beta;
	}

	double count_mu(double beta, int k_step, const TMatrix<T>& Ak)
	{
		double mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * Ak[k_step][k_step]);
		return mu;
	}

	TVector<T> count_w(const TMatrix<T>& Ak, int k_step, double beta, double mu)
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
};
