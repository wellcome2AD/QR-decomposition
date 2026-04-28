#pragma once

#include <vector>
#include <iostream>
#include <cmath>

#include "../IQRSolver.h"
#include "math.h"

template <typename T>
class GivensMethodBasic : public IQRSolver<T> {
public:
	// первая, неоптимальная версия
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		auto N = A.size();
		Q = std::vector<std::vector<T>>(N, std::vector<T>(N, 0));
		R = A;
		for (int i = 0; i < N; ++i)
		{
			Q[i][i] = 1.0;
		}

		for (auto j = 0; j < N - 1; j++)
		{ // зануляем весь столбец j
			for (auto i = j + 1; i < N; i++)
			{ // под главной диагональю
				double Rjj = R[j][j];
				double Rij = R[i][j];

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

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
				for (auto k = 0; k < N; k++)
				{
					if (k >= j)
					{// меняются только две строки -- i и j
						auto temp = R[j][k] * c - R[i][k] * s;
						R[i][k] = R[j][k] * s + R[i][k] * c;
						R[j][k] = temp;
					}

					// меняются только два столбца -- i и j
					auto temp = c * Q[k][j] - s * Q[k][i];
					Q[k][i] = s * Q[k][j] + c * Q[k][i];
					Q[k][j] = temp;
				}
			}
		}
	}
};
