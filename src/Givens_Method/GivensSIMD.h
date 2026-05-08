#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <memory>
#include <new>
#include <algorithm>

#include "../IQRSolver.h"
#include "math.h"

static constexpr std::align_val_t c_align = std::align_val_t(64);

template <typename T>
class GivensMethodSIMD : public IQRSolver<T> {
public:
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		const auto N = A.size();

		// выделение памяти с выравниванием строк
		constexpr size_t elem_size = sizeof(T);
		constexpr size_t align_mod = 64 / elem_size;

		size_t stride = N; // каждая строка будет начинаться с выровненного адреса
		size_t rem = N % align_mod;
		if (rem != 0) stride += align_mod - rem; // stride кратен 8

		// обёртки для автоматического освобождения памяти
		auto deleter = [](T* p) { ::operator delete[](p, c_align); };
		std::unique_ptr<T[], decltype(deleter)> R_data{
			static_cast<T*>(::operator new[](N* stride* elem_size, c_align)), deleter };
		std::unique_ptr<T[], decltype(deleter)> Q_data{
			static_cast<T*>(::operator new[](N* stride* elem_size, c_align)), deleter };

#define R(i, j) R_data[(i)*stride + (j)] // функция доступа: R_data[i * stride + j]
#define Q_T(i, j) Q_data[(i)*stride + (j)] // Q в транспонированном виде

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i) {
			for (size_t j = 0; j < N; ++j) {
				R(i, j) = A[i][j];
			}
			Q_T(i, i) = 1.0;
		}

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i) {
			for (size_t j = 0; j < N; ++j) {
				if (i != j) Q_T(i, j) = 0.0;
			}
		}

		for (int j = 0; j < N - 1; ++j) {
			for (int i = j + 1; i < N; ++i) {
				const double Rjj = R(j, j);
				const double Rij = R(i, j);

				const double denom = std::max(std::abs(Rjj), 1.0);
				if (std::abs(Rij) < 1e-11 * denom ||
					std::abs(Rij) < 1e-11 * std::abs(Rjj))
					continue;

				const double sqrt_val = std::sqrt(Rjj * Rjj + Rij * Rij);
				if (std::abs(sqrt_val) < 1e-11)
					continue;

				const double c = Rjj / sqrt_val;
				const double s = -Rij / sqrt_val;

				// локальные указатели на строки с выравниванием
				T* __restrict Rj_ptr = &R_data[j * stride];
				T* __restrict Ri_ptr = &R_data[i * stride];
				T* __restrict Qj_ptr = &Q_data[j * stride];
				T* __restrict Qi_ptr = &Q_data[i * stride];

#pragma omp simd aligned(Rj_ptr, Ri_ptr, Qj_ptr, Qi_ptr : 64)
				for (size_t k = 0; k < stride; ++k) {
					const T Rj = Rj_ptr[k];
					const T Ri = Ri_ptr[k];
					const T Qj = Qj_ptr[k];
					const T Qi = Qi_ptr[k];

					Rj_ptr[k] = Rj * c - Ri * s;
					Ri_ptr[k] = Rj * s + Ri * c;
					Qj_ptr[k] = Qj * c - Qi * s;
					Qi_ptr[k] = Qj * s + Qi * c;
				}
			}
		}

		R = Q = std::vector<std::vector<T>>(N, std::vector<T>(N));
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i) {
			for (size_t j = 0; j < N; ++j) {
				R[i][j] = R(i, j);
				Q[i][j] = Q_T(j, i);
			}
		}

#undef R
#undef Q_T
	}
};
