#pragma once

#include <vector>

#include "../IQRSolver.h"

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
		std::vector<T> all_v(N * (N - 1), 0.0); // все вектора Хаусхолдера будут хранится в одной матрице подряд, без незначащих нулей

		// вычисление R
		// H * R = (I - 2 * v * v^T) * R = R - 2 * v * (v^T * R)
		for (ptrdiff_t k = 0; k < N - 1; ++k)
		{
			T* v_k = &all_v[k * N]; // часть матрицы с вектором vk
			ptrdiff_t v_length = N - k; // длина вектора vk
			ptrdiff_t kN = k * N; // предвычисленное смещение для k строки

			T* R_col_k = &R_flat[kN + k]; // 0 элемент k столбца матрицы R, чтобы получить следующий -- прибавить строку

			double beta = 0.0;
			// count beta
			double sum = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum)
			for (int i = 0; i < v_length; ++i) {
				sum += R_col_k[i * N] * R_col_k[i * N];
			}
			if (fabs(sum) < std::numeric_limits<T>::epsilon()) // eсли знаменатель близок к 0, вектор отражения не нужен
			{
				continue;
			}
			double norm_x = sqrt(sum);
			beta = (R_flat[kN + k] >= 0.0) ? -norm_x : norm_x;

			double mu = 0.0;
			// count mu
			double denominator = 2.0 * beta * (beta - R_flat[kN + k]);
			if (fabs(denominator) < std::numeric_limits<T>::epsilon()) // eсли знаменатель близок к 0, вектор отражения не нужен
			{
				continue;
			}
			mu = 1.0 / sqrt(denominator);

			// count v
#pragma omp parallel for num_threads(thread_num)
			for (int i = 0; i < v_length; ++i) {
				v_k[i] = R_col_k[i * N] * mu;
			}
			v_k[0] -= beta * mu;

			// отражение
			// R - 2 * v * (v^T * R)
			std::vector<double> _2vTR(N - k); // предварительно вычисляем v^TR для каждого столбца
#pragma omp parallel for num_threads(thread_num)
			for (ptrdiff_t j = k; j < N; ++j) {
				T* R_col_j = &R_flat[kN + j]; // начальный элемент для умножения из j-ого столбца R, прибавить строку
				double vTR = 0.0;
				for (int i = 0; i < v_length; ++i) {
					vTR += v_k[i] * R_col_j[i * N];
				}
				_2vTR[j - k] = 2.0 * vTR;
			}

			// обновление R
			// R - v * 2v^TR
#pragma omp parallel for num_threads(thread_num)
			for (ptrdiff_t j = k; j < N; ++j) {
				for (int i = 0; i < v_length; ++i) {
					T* R_col_j = &R_flat[kN + j];
					R_col_j[i * N] -= v_k[i] * _2vTR[j - k];
				}
			}
		}

		// вычисление Q
		// начать с Q_{N-1} = H_{N-2} * I = I - 2 * v_{N-2} * (v_{N-2}^T * I) и домножать слева на матрицы отражения
		// Q = H_{1} * ... * (H_{k} * Q_{k-1} ) = H_{1} * ... * (Q_{k-1} - 2 * v_{k} * (v_{k}^T * Q_{k-1})
		Q_flat.assign(N * N, 0.0); // создание единичной матрицы
#pragma omp parallel for num_threads(thread_num)
		for (ptrdiff_t i = 0; i < N; ++i) {
			Q_flat[i * N + i] = 1.0;
		}

		std::vector<T> _2vTQ(N, 0.0);
		for (ptrdiff_t k = N - 2; k >= 0; --k)
		{
			T* v_k = &all_v[k * N]; // взятие вектора vk
			ptrdiff_t v_length = N - k; // длина вектора vk

			// 2v^TQ параллельно по столбцам
#pragma omp parallel for num_threads(thread_num)
			for (ptrdiff_t j = 0; j < N; ++j) {
				double sum = 0.0;
				for (int i = 0; i < v_length; ++i) {
					ptrdiff_t row = k + i;
					sum += v_k[i] * Q_flat[row * N + j];
				}
				_2vTQ[j] = 2.0 * sum;
			}

			// Q - v * 2v^TQ параллельно по строкам
#pragma omp parallel for num_threads(thread_num)
			for (int i = 0; i < v_length; ++i) {
				ptrdiff_t row = k + i;
				double vi = v_k[i];
				T* Q_row = &Q_flat[row * N];

				for (ptrdiff_t j = 0; j < N; ++j) {
					Q_row[j] -= vi * _2vTQ[j];
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
