#pragma once

#include <vector>
#include <iostream>
#include <cmath>

#include "IGivensMethodSolver.h"
#include "math.h"

template <typename T>
class GivensMethodBasic : public IGivensMethodSolver<T> {
public:
	// первая, неоптимальная версия
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		auto N = A.size();
		// Конвертируем входные данные в плоские массивы
		std::vector<T> A_flat(N * N);
//#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				A_flat[i * N + j] = A[i][j];
			}
		}
		auto R_flat = A_flat;

		TransposedMatrix<T> Q_T_flat(N, N, 0.0);
//#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (auto i = 0; i < N; i++)
		{
			Q_T_flat.At(i, i) = 1.0;
		}

		for (auto j = 0; j < N - 1; j++)
		{ // зануляем весь столбец j
			for (auto i = j + 1; i < N; i++)
			{ // c j+1 строк
				double Rjj = R_flat[j * N + j];
				double Rij = R_flat[i * N + j];

				if (std::abs(Rij) < 1e-11 * std::max(std::abs(Rjj), 1.0)) {
					continue; // относительные чисал
				}

				if (std::abs(Rij) < 1e-11 * std::abs(Rjj)) { // R[i][j] пренебрежимо мал по сравнению с R[j][j]					
					continue;
				}

				if (std::abs(Rij) < 1e-11) { // уже 0
					continue;
				}

				auto sqrt = std::sqrt(Rjj * Rjj + Rij * Rij);
				if (std::abs(sqrt) < 1e-11) { // для корректировки ошибки игнорируется малый знаменатель
					continue;
				}
				double c = Rjj / sqrt, s = -Rij / sqrt;

//#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (int k = j; k < N; k++)
				{ // вращаем только две строки -- i и j
					auto temp = R_flat[j * N + k] * c - R_flat[i * N + k] * s;
					R_flat[i * N + k] = R_flat[j * N + k] * s + R_flat[i * N + k] * c;
					R_flat[j * N + k] = temp;
				}

//#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (auto k = 0; k < N; k++)
				{ // меняются только два столбца -- i и j
					auto temp = c * Q_T_flat.At(k, j) - s * Q_T_flat.At(k, i);
					Q_T_flat.At(k, i) = s * Q_T_flat.At(k, j) + c * Q_T_flat.At(k, i);
					Q_T_flat.At(k, j) = temp;
				}
			}
		}

		R = Q = std::vector<std::vector<T>>(N, std::vector<T>(N));
#pragma omp parallel for num_threads(thread_num) collapse(2)
		for (int j = 0; j < N; ++j)
		{
			for (int i = 0; i < N; ++i)
			{
				R[i][j] = R_flat[i * N + j];
				Q[i][j] = Q_T_flat.At(i, j);
			}
		}
	}
};
