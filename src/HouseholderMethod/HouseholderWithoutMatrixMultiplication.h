#pragma once

#include <vector>

#include "../IQRSolver.h"
#include "../Matrix/MatrixOperations.h"

template <typename T>
class HouseholderMethodWithoutMatrixMults : public IQRSolver<T> {
public:
	// вторая версия, избавлена от недостатка первой -- двух матричных умножений на каждой итерации
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		const size_t N = A.size();
		// Конвертируем входные данные в плоские массивы
		std::vector<T> A_flat(N * N), Q_flat;
		for (ptrdiff_t j = 0; j < N; ++j) {
			for (ptrdiff_t i = 0; i < N; ++i) {
				A_flat[i * N + j] = A[i][j];
			}
		}

		auto R_flat = A_flat;
		std::vector<T> all_w(N * (N - 1), 0.0); // все вектора Хаусхолдера будут хранится в одной матрице подряд, без незначащих нулей

		// вычисление R
		// H * R = (I - 2 * w * wT) * R = R - 2 * w * (wT * R)
		for (ptrdiff_t k = 0; k < N - 1; ++k)
		{
			T* w_k = &all_w[k * N]; // часть матрицы с вектором wk
			ptrdiff_t w_length = N - k; // длина вектора wk
			ptrdiff_t kN = k * N; // предвычисленное смещение для k строки

			T* R_col_k = &R_flat[kN + k]; // 0 элемент k столбца матрицы R, чтобы получить следующий -- прибавить строку
			double beta = 0.0;
			// count beta
			{
				double sum_by_k_col = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum_by_k_col) if(N >= 1000)
				for (int i = 0; i < w_length; ++i) {
					T val = R_col_k[i * N];
					sum_by_k_col += static_cast<double>(val) * static_cast<double>(val);
				}
				if (fabs(sum_by_k_col) < std::numeric_limits<T>::epsilon()) // eсли знаменатель близок к 0, вектор отражения не нужен
				{
					continue;
				}
				double norm_x = sqrt(sum_by_k_col);
				beta = (R_flat[kN + k] >= 0.0) ? -norm_x : norm_x;
			}

			double mu = 0.0;
			// count mu
			{
				double denominator = 2.0 * beta * (beta - R_flat[kN + k]);
				if (fabs(denominator) < std::numeric_limits<T>::epsilon()) // eсли знаменатель близок к 0, вектор отражения не нужен
				{
					continue;
				}
				mu = 1.0 / sqrt(denominator);
			}

			// count w
#pragma omp parallel for num_threads(thread_num) if(N >= 1000)
			for (int i = 0; i < w_length; ++i) {
				w_k[i] = R_col_k[i * N] * mu;
			}
			w_k[0] -= beta * mu;

			// отражение
			// R - 2 * w * (wT * R)
			std::vector<double> _2wTR(N - k); // предварительно вычисляем wTR для каждого столбца
#pragma omp parallel for num_threads(thread_num) if(N >= 1000)
			for (ptrdiff_t j = k; j < N; ++j) {
				T* R_col_j = &R_flat[kN + j]; // начальный элемент для умножения из j-ого столбца R, прибавить строку
				double wTR = 0.0;
				for (int i = 0; i < w_length; ++i) {
					wTR += w_k[i] * R_col_j[i * N];
				}
				_2wTR[j - k] = 2.0 * wTR;
			}

			// обновление R
			// R - w * 2wTR
#pragma omp parallel for num_threads(thread_num) if(N >= 1000)
			for (ptrdiff_t j = k; j < N; ++j) {
				for (int i = 0; i < w_length; ++i) {
					T* R_col_j = &R_flat[kN + j];
					R_col_j[i * N] -= w_k[i] * _2wTR[j - k];
				}
			}
		}

		// вычисление Q
		// начать с Q_{N-1} = H_{N-2} * I = I - 2 * w_{N-2} * (w_{N-2}T * I) и домножать слева на матрицы отражения
		// Q = H_{1} * ... * (H_{k} * Q_{k-1} ) = H_{1} * ... * (Q_{k-1} - 2 * w_{k} * (w_{k}T * Q_{k-1})
		Q_flat.assign(N * N, 0.0); // создание единичной матрицы
#pragma omp parallel for num_threads(thread_num) if(N >= 1000)
		for (ptrdiff_t i = 0; i < N; ++i) {
			Q_flat[i * N + i] = 1.0;
		}

		std::vector<T> _2wTQ(N, 0.0);
		for (ptrdiff_t k = N - 2; k >= 0; --k)
		{
			T* w_k = &all_w[k * N]; // взятие вектора wk
			ptrdiff_t w_length = N - k; // длина вектора wk

			// 2wTQ параллельно по столбцам
#pragma omp parallel for num_threads(thread_num) if(N >= 1000)
			for (ptrdiff_t j = 0; j < N; ++j) {
				double sum = 0.0;
				for (int i = 0; i < w_length; ++i) {
					ptrdiff_t row = k + i;
					sum += w_k[i] * Q_flat[row * N + j];
				}
				_2wTQ[j] = 2.0 * sum;
			}

			// Q - w * 2wTQ параллельно по строкам
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
			for (int i = 0; i < w_length; ++i) {
				ptrdiff_t row = k + i;
				double wi = w_k[i];
				T* Q_row = &Q_flat[row * N];

				for (ptrdiff_t j = 0; j < N; ++j) {
					Q_row[j] -= wi * _2wTQ[j];
				}
			}
		}

#pragma omp parallel for num_threads(thread_num) collapse(2)
		for (ptrdiff_t i = 0; i < N; ++i) {
			for (ptrdiff_t j = 0; j < N; ++j) {
				R[i][j] = R_flat[i * N + j];
				Q[i][j] = Q_flat[i * N + j];
			}
		}
	}
};
