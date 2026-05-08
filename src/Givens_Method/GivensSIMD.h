#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <memory>
#include <new>

#include "../IQRSolver.h"
#include "math.h"

static constexpr auto c_align = 64;

template <typename T>
class GivensMethodBasic : public IQRSolver<T> {
public:
	// первая, неоптимальная версия
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		auto N = A.size();
		auto R_aligned = static_cast<T *>(::operator new[](N * N * sizeof(T), std::align_val_t(c_align)));
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				R_aligned[i * N + j] = A[i][j];
			}
		}

		auto Q_aligned = static_cast<T*>(::operator new[](N* N * sizeof(T), std::align_val_t(c_align)));
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i)
		{
			Q_aligned[i * N + i] = 1.0;
		}

		auto&& Q_T = Q_aligned;

		for (auto j = 0; j < N - 1; j++)
		{
			for (auto i = j + 1; i < N; i++)
			{
				double Rjj = R[j][j];
				double Rij = R[i][j];

				double denom = std::max(std::abs(Rjj), 1.0);
				if (std::abs(Rij) < 1e-11 * denom || std::abs(Rij) < 1e-11 * std::abs(Rjj))
					continue;

				double sqrt_val = std::sqrt(Rjj * Rjj + Rij * Rij);
				if (std::abs(sqrt_val) < 1e-11)
					continue;

				double c = Rjj / sqrt_val, s = -Rij / sqrt_val;

#pragma omp simd aligned(R_aligned, Q_T: 64)
				for (auto k = 0; k < N; k++) {
					T Rj =  R_aligned[j * N + k];
					T Ri = R_aligned[i * N + k];
					T Qj = Q_T[j * N + k];
					T Qi = Q_T[i * N + k];

					R_aligned[j * N + k] = Rj * c - Ri * s;
					R_aligned[i * N + k] = Rj * s + Ri * c;
					Q_T[j * N + k] = Qj * c - Qi * s;
					Q_T[i * N + k] = Qj * s + Qi * c;
				}
			}
		}

		R = Q = std::vector<std::vector<T>>(N, std::vector<T>(N));
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				R[i][j] = R_aligned[i * N + j];
				Q[i][j] = Q_T[j * N + i];
			}
		}

		::operator delete[](R_aligned, std::align_val_t(c_align));
		::operator delete[](Q_aligned, std::align_val_t(c_align));
	}
};
