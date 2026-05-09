#pragma once

#include <immintrin.h>   // AVX2
#include <vector>
#include <cmath>
#include <algorithm>

#include "../IQRSolver.h"
#include "math.h"

template <typename T>
class GivensMethodSIMD : public IQRSolver<T> {
public:
    virtual void QR_decomposition(const std::vector<std::vector<T>>& A,
        std::vector<std::vector<T>>& Q,
        std::vector<std::vector<T>>& R) override
    {
        const auto N = A.size();
        if (N == 0) return;

        // Плоские массивы для R и Q_T (транспонированная Q)
        size_t total = static_cast<size_t>(N) * N;
        auto R_data = std::make_unique<T[]>(total);
        auto Q_data = std::make_unique<T[]>(total);   // Q_T

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                R_data[i * N + j] = A[i][j];
                Q_data[i * N + j] = (i == j) ? T(1) : T(0);
            }
        }

        for (int j = 0; j < N - 1; ++j) {
            for (int i = j + 1; i < N; ++i) {
                const double Rjj = R_data[j * N + j];
                const double Rij = R_data[i * N + j];

 
                if (std::abs(Rij) < 1e-11)
                    continue;

                const double sqrt_val = std::sqrt(Rjj * Rjj + Rij * Rij);
                if (std::abs(sqrt_val) < 1e-11)
                    continue;

                const double c = Rjj / sqrt_val;
                const double s = -Rij / sqrt_val;

                // Указатели на начала строк
                T* Rj_ptr = &R_data[j * N];
                T* Ri_ptr = &R_data[i * N];
                T* Qj_ptr = &Q_data[j * N];
                T* Qi_ptr = &Q_data[i * N];

                // Векторный цикл (по 4 double, AVX2)
                // Если N не кратен 4, оставшийся хвост обработаем после цикла
                size_t k = 0;
                for (; k + 3 < N; k += 4) {
                    __m256d Rj = _mm256_loadu_pd(&Rj_ptr[k]);
                    __m256d Ri = _mm256_loadu_pd(&Ri_ptr[k]);
                    __m256d Qj = _mm256_loadu_pd(&Qj_ptr[k]);
                    __m256d Qi = _mm256_loadu_pd(&Qi_ptr[k]);

                    __m256d c_vec = _mm256_set1_pd(c);
                    __m256d s_vec = _mm256_set1_pd(s);

                    // Обновление R и Q (FMA)
                    __m256d new_Rj = _mm256_fmsub_pd(Rj, c_vec, _mm256_mul_pd(Ri, s_vec));
                    __m256d new_Ri = _mm256_fmadd_pd(Rj, s_vec, _mm256_mul_pd(Ri, c_vec));
                    __m256d new_Qj = _mm256_fmsub_pd(Qj, c_vec, _mm256_mul_pd(Qi, s_vec));
                    __m256d new_Qi = _mm256_fmadd_pd(Qj, s_vec, _mm256_mul_pd(Qi, c_vec));

                    _mm256_storeu_pd(&Rj_ptr[k], new_Rj);
                    _mm256_storeu_pd(&Ri_ptr[k], new_Ri);
                    _mm256_storeu_pd(&Qj_ptr[k], new_Qj);
                    _mm256_storeu_pd(&Qi_ptr[k], new_Qi);
                }

                // Обработка остатка (если N не кратно 4)
                for (; k < N; ++k) {
                    double Rj_val = Rj_ptr[k];
                    double Ri_val = Ri_ptr[k];
                    double Qj_val = Qj_ptr[k];
                    double Qi_val = Qi_ptr[k];

                    Rj_ptr[k] = Rj_val * c - Ri_val * s;
                    Ri_ptr[k] = Rj_val * s + Ri_val * c;
                    Qj_ptr[k] = Qj_val * c - Qi_val * s;
                    Qi_ptr[k] = Qj_val * s + Qi_val * c;
                }
            }
        }

        Q.assign(N, std::vector<T>(N));
        R.assign(N, std::vector<T>(N));
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                R[i][j] = R_data[i * N + j];
                Q[i][j] = Q_data[j * N + i];
            }
        }
    }
};
