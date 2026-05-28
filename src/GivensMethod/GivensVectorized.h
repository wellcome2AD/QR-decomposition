#pragma once

#include <vector>
#include <cmath>
#include <memory>
#include <immintrin.h>
#include <type_traits>

#include "../IQRSolver.h"

template <typename T>
class GivensVectorized : public IQRSolver<T> {
public:
    virtual void QR_decomposition(const std::vector<std::vector<T>>& A,
        std::vector<std::vector<T>>& Q,
        std::vector<std::vector<T>>& R) override
    {
        const int N = A.size();
        if (N == 0) return;

        size_t total = static_cast<size_t>(N) * N;
        auto R_data = std::make_unique<T[]>(total);
        auto Q_data = std::make_unique<T[]>(total);   // Q_T

        // инициализация: R = A, Q = E
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                R_data[i * N + j] = A[i][j];
                Q_data[i * N + j] = (i == j) ? T(1) : T(0);
            }
        }

        // прямой ход: вычисление R
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

                T* Rj_ptr = &R_data[j * N];
                T* Ri_ptr = &R_data[i * N];

                // векторизованное вращение строк R (ненулевые значения только в столбцах правее j)
                if constexpr (std::is_same_v<T, double>) {
                    size_t k = j;
                    for (; k + 3 < N; k += 4) {
                        __m256d Rj = _mm256_loadu_pd(&Rj_ptr[k]);
                        __m256d Ri = _mm256_loadu_pd(&Ri_ptr[k]);

                        __m256d c_vec = _mm256_set1_pd(c);
                        __m256d s_vec = _mm256_set1_pd(s);

                        __m256d new_Rj = _mm256_fmsub_pd(Rj, c_vec, _mm256_mul_pd(Ri, s_vec));
                        __m256d new_Ri = _mm256_fmadd_pd(Rj, s_vec, _mm256_mul_pd(Ri, c_vec));

                        _mm256_storeu_pd(&Rj_ptr[k], new_Rj);
                        _mm256_storeu_pd(&Ri_ptr[k], new_Ri);
                    }
                    for (; k < N; ++k) {
                        T Rj = Rj_ptr[k], Ri = Ri_ptr[k];
                        Rj_ptr[k] = Rj * c - Ri * s;
                        Ri_ptr[k] = Rj * s + Ri * c;
                    }
                }
                else if constexpr (std::is_same_v<T, float>) {
                    size_t k = j;
                    for (; k + 7 < N; k += 8) {
                        __m256 Rj = _mm256_loadu_ps(&Rj_ptr[k]);
                        __m256 Ri = _mm256_loadu_ps(&Ri_ptr[k]);

                        __m256 c_vec = _mm256_set1_ps(c);
                        __m256 s_vec = _mm256_set1_ps(s);

                        __m256 new_Rj = _mm256_fmsub_ps(Rj, c_vec, _mm256_mul_ps(Ri, s_vec));
                        __m256 new_Ri = _mm256_fmadd_ps(Rj, s_vec, _mm256_mul_ps(Ri, c_vec));

                        _mm256_storeu_ps(&Rj_ptr[k], new_Rj);
                        _mm256_storeu_ps(&Ri_ptr[k], new_Ri);
                    }
                    for (; k < N; ++k) {
                        T Rj = Rj_ptr[k], Ri = Ri_ptr[k];
                        Rj_ptr[k] = Rj * c - Ri * s;
                        Ri_ptr[k] = Rj * s + Ri * c;
                    }
                }

                // сохраняем c и s в виде tau
                double tau = s / (1.0 + c);
                R_data[i * N + j] = static_cast<T>(tau);
            }
        }

        // обратный ход: построение Q_T по сохранённым tau
        for (int j = 0; j < N - 1; ++j) {
            for (int i = j + 1; i < N; ++i) {
                if (std::abs(R_data[i * N + j]) < 1e-11) continue;

                double tau = static_cast<double>(R_data[i * N + j]);
                double tau2 = tau * tau;
                double denom = 1.0 + tau2;
                double c = (1.0 - tau2) / denom;
                double s = 2.0 * tau / denom;   // восстанавливаем c и s

                T* Qj_ptr = &Q_data[j * N];
                T* Qi_ptr = &Q_data[i * N];

                // применяем вращение строк Q_T
                if constexpr (std::is_same_v<T, double>) {
                    size_t k = 0;
                    for (; k + 3 < N; k += 4) {
                        __m256d Qj = _mm256_loadu_pd(&Qj_ptr[k]);
                        __m256d Qi = _mm256_loadu_pd(&Qi_ptr[k]);

                        __m256d c_vec = _mm256_set1_pd(c);
                        __m256d s_vec = _mm256_set1_pd(s);

                        __m256d new_Qj = _mm256_fmsub_pd(Qj, c_vec, _mm256_mul_pd(Qi, s_vec));
                        __m256d new_Qi = _mm256_fmadd_pd(Qj, s_vec, _mm256_mul_pd(Qi, c_vec));

                        _mm256_storeu_pd(&Qj_ptr[k], new_Qj);
                        _mm256_storeu_pd(&Qi_ptr[k], new_Qi);
                    }
                    for (; k < N; ++k) {
                        T Qj = Qj_ptr[k], Qi = Qi_ptr[k];
                        Qj_ptr[k] = Qj * c - Qi * s;
                        Qi_ptr[k] = Qj * s + Qi * c;
                    }
                }
                else if constexpr (std::is_same_v<T, float>) {
                    size_t k = 0;
                    for (; k + 7 < N; k += 8) {
                        __m256 Qj = _mm256_loadu_ps(&Qj_ptr[k]);
                        __m256 Qi = _mm256_loadu_ps(&Qi_ptr[k]);

                        __m256 c_vec = _mm256_set1_ps(c);
                        __m256 s_vec = _mm256_set1_ps(s);

                        __m256 new_Qj = _mm256_fmsub_ps(Qj, c_vec, _mm256_mul_ps(Qi, s_vec));
                        __m256 new_Qi = _mm256_fmadd_ps(Qj, s_vec, _mm256_mul_ps(Qi, c_vec));

                        _mm256_storeu_ps(&Qj_ptr[k], new_Qj);
                        _mm256_storeu_ps(&Qi_ptr[k], new_Qi);
                    }
                    for (; k < N; ++k) {
                        T Qj = Qj_ptr[k], Qi = Qi_ptr[k];
                        Qj_ptr[k] = Qj * c - Qi * s;
                        Qi_ptr[k] = Qj * s + Qi * c;
                    }
                }

                R_data[i * N + j] = T(0);   // восстанавливаем чистый ноль
            }
        }

        // копирование в выходные матрицы
        Q.assign(N, std::vector<T>(N));
        R.assign(N, std::vector<T>(N));
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                R[i][j] = R_data[i * N + j];
                Q[i][j] = Q_data[j * N + i]; // транспонирование
            }
        }
    }
};
