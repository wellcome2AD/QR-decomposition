#pragma once

#include "TMatrix.h"

inline static float sign(float x)
{
	if (x > 0.0)
		return 1.0;
	if (x < 0.0)
		return -1.0;
	return 1.0;
}

template <typename T>
static double count_beta(const TMatrix<T>& Ak, int k_step)
{
	const size_t N = Ak.Size();
	double sum_by_k_col = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(+:sum_by_k_col)
	for (int str = k_step; str < N; ++str) {
		sum_by_k_col += Ak[str][k_step] * Ak[str][k_step];
	}
	auto beta = sign(-Ak[k_step][k_step]) * sqrt(sum_by_k_col);
	return beta;
}

template <typename T>
static double count_mu(float beta, int k_step, const TMatrix<T>& Ak)
{
	double mu = 1 / sqrt(2.0 * beta * beta - 2 * beta * Ak[k_step][k_step]);
	return mu;
}

template <typename T>
static TVector<T> count_w(const TMatrix<T>& Ak, int k_step, float beta, float mu)
{
	const size_t N = Ak.Size();
	TVector<T> w(N);
#pragma omp parallel for num_threads(thread_num)
	for (int str = 0; str < N; ++str) {
		w[str] = str < k_step ? 0 : Ak[str][k_step] * mu;
	}
	w[k_step] -= beta * mu;
	return w;
}

template <typename T>
static TMatrix<T> count_H(float beta, float mu, TVector<T> w)
{
	size_t N = w.Size();
	TMatrix<T> H(N, N, 0);
#pragma omp parallel for num_threads(thread_num)
	for (int i = 0; i < N; ++i) {
		H[i][i] = 1.0;
		for (int j = 0; j < N; ++j) {
			H[i][j] -= w[i] * w[j] * 2.0;
		}
	}
	return H;
}

template <typename T>
void Household_QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R)
{
	const size_t N = A.Size();
	R = A;
	bool isFirst = true;
	for (size_t k = 0; k < N - 1; ++k) {
		auto beta = count_beta(R, k);
		auto mu = count_mu(beta, k, R);
		auto w = count_w(R, k, beta, mu);

		auto H = count_H(beta, mu, w);

		if (k == 0) {
			Q = H;
		}
		else {
			Q = Q * H;
		}

		R = H * R;
	}

	return;
}

template <typename T>
void Household_QR_decomposition_experimental(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R)
{
	const size_t N = A.Size();
	R = A;

	Q = TMatrix<T>(N, N, 0);
	for (auto i = 0; i < N; i++)
	{
		Q[i][i] = 1;
	}

	for (size_t k = 0; k < N - 1; ++k) {
		auto beta = count_beta(R, k);
		auto mu = count_mu(beta, k, R);
		auto w = count_w(R, k, beta, mu);

		// H * R = (I - 2 * w * wᵀ) * R = R - 2 * w * (wT * R)
		// формируем (wT * R)
		TMatrix<T> y(1, N, 0);
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				y[0][i] += w[j] * R[j][i];
			}
		}

		// y = 2 * y
		for (size_t i = 0; i < N; ++i)
		{
			y[0][i] *= 2.0;
		}

		// w * y = матрица размером n×m, где (i,j)-й элемент = w[i] * y[j]
		// w это вектор nx1, y это вектор 1xm
		TMatrix<T> wy(N, N);
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				wy[i][j] = w[i] * y[0][j];
			}
		}

		// R = R - w * y   (w * y — внешнее произведение, outer product)
		R = R - wy;
	}

	return;
}
//
//template<typename T>
//void build_q_inplace(TMatrix<T>& A, const std::vector<T>& tau, TMatrix<T>& Q) {
//	int n = A.Size();
//
//	// Инициализируем Q как единичную матрицу
//	Q = TMatrix<T>(n, n, 0);
//	for (int i = 0; i < n; i++) {
//		Q[i][i] = 1;
//	}
//
//	// Применяем отражения в обратном порядке
//	for (int k = n - 2; k >= 0; k--) {  // для n отражений (n-1 шагов)
//		if (std::abs(tau[k]) < 1e-12) continue;
//
//		// Для каждого столбца Q применяем отражение
//		for (int j = 0; j < n; j++) {
//			// Вычисляем vᵀ * Q[:,j]
//			T dot = Q[k][j];  // v[k] = 1
//
//			for (int i = k + 1; i < n; i++) {
//				dot += A[i][k] * Q[i][j];  // v[i] * Q[i][j]
//			}
//
//			dot *= tau[k];
//
//			// Обновляем Q[:,j] = Q[:,j] - v * dot
//			Q[k][j] -= dot;  // для v[k] = 1
//
//			for (int i = k + 1; i < n; i++) {
//				Q[i][j] -= dot * A[i][k];
//			}
//		}
//	}
//}
//
//template<typename T>
//void householder_qr_compact(TMatrix<T>& A, std::vector<T>& tau) {
//	int n = A.Size();
//	tau.resize(n - 1);
//
//	for (int k = 0; k < n - 1; k++) {
//		// 1. Вычисляем норму части столбца
//		T norm_x = 0;
//		for (int i = k; i < n; i++) {
//			norm_x += A[i][k] * A[i][k];
//		}
//		norm_x = std::sqrt(norm_x);
//
//		if (norm_x == 0) {
//			tau[k] = 0;
//			continue;
//		}
//
//		// 2. Вычисляем σ и первый элемент вектора v
//		T sign = (A[k][k] >= 0) ? 1 : -1;
//		T sigma = -sign * norm_x;
//		T v0 = A[k][k] - sigma;
//
//		// 3. Сохраняем τ (это и есть β из формулы H = I - τvvᵀ)
//		tau[k] = -sigma / (norm_x * norm_x - sigma * A[k][k]);
//
//		// 4. Сохраняем вектор v в столбце A (под диагональю)
//		A[k][k] = sigma;  // диагональный элемент R
//
//		// Масштабируем остальные элементы v
//		for (int i = k + 1; i < n; i++) {
//			A[i][k] /= v0;  // теперь A[i][k] = v[i] (при v[k]=1)
//		}
//
//		// 5. Применяем отражение к правой части матрицы A
//		//    A[k:n, k+1:n] = (I - τ v vᵀ) * A[k:n, k+1:n]
//		for (int j = k + 1; j < n; j++) {
//			// Вычисляем w = τ * (vᵀ * A[k:n, j])
//			T dot = A[k][j];  // v[k] = 1
//
//			for (int i = k + 1; i < n; i++) {
//				dot += A[i][k] * A[i][j];  // v[i] * A[i][j]
//			}
//
//			dot *= tau[k];
//
//			// Вычитаем: A[k:n, j] = A[k:n, j] - w * v
//			A[k][j] -= dot;  // для v[k] = 1
//
//			for (int i = k + 1; i < n; i++) {
//				A[i][j] -= dot * A[i][k];
//			}
//		}
//	}
//}
//
//// Функция для получения Q и R из компактного представления
//template<typename T>
//void extract_qr(TMatrix<T>& A_compact, const std::vector<T>& tau,
//	TMatrix<T>& Q, TMatrix<T>& R) {
//	int n = A_compact.Size();
//
//	// 1. Извлекаем R (верхняя треугольная часть)
//	R = TMatrix<T>(n, n, 0);
//	for (int i = 0; i < n; i++) {
//		for (int j = i; j < n; j++) {
//			R[i][j] = A_compact[i][j];
//		}
//	}
//
//	// 2. Строим Q из сохранённых векторов
//	build_q_inplace(A_compact, tau, Q);
//}
