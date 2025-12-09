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
		std::cout << "H" << k + 1 << ":" << H << std::endl;

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

	// вычисление R
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

		// w * y = матрица размером n×n, где (i,j)-й элемент = w[i] * y[j]
		// w это вектор nx1, y это вектор 1xn
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

		// сохраним в столбце k вектор wk с элемента k+1
		for (auto i = k + 1; i < N; ++i)
		{
			R[i][k] = w[i] / w[k];
		}
	}
	std::cout << "R:" << R << std::endl;
	Q = TMatrix<T>(N, N, 0);
	for (auto i = 0; i < N; i++)
	{
		Q[i][i] = 1;
	}

	/*
	// теперь для вычисления Q достанем вектора wk из R
	for (ptrdiff_t k = N - 2; k >= 0; --k) {
		// первые k элементов =0, w[k]=1;
		TVector<T> w(N, 0);
		w[k] = 1;
		double sum = 0.0;
		for (size_t i = k + 1; i < N; ++i) {
			w[i] = R[i][k];
			R[i][k] = 0.0; // занулим заодно эти элементы в R
			sum += w[i] * w[i];
		}
		// на диагонали хранится sigma, вычислим tau
		// tau = -sigma / (sigma^2 + sum_i( w[i]^2 )  где i = k+1..N-1
		// или эквивалентно beta = 2 * tau = -2sigma / (sigma^2 + sum_i( w[i]^2 )
		// double sigma = R[k][k];
		//double tau = -sigma / (sigma * sigma + sum); // todo: это надо вообще? у меня стоит коэффициент 2

		// Q = H_{N-2} * Q = Q - tau * w * (wT * Q)
		TMatrix<T> wTQ(1, N, 0);
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j) // todo: нужно оптимизировать, умножать не полностью; видимо, с k-го элемента
			{
				wTQ[0][i] += w[j] * Q[j][i];
			}
		}

		// w = τ * w
		for (size_t i = k; i < N; ++i)
		{
			w[i] *= 2; //tau;
		}

		// w * wTQ = матрица размером n×n, где (i,j)-й элемент = w[i] * y[j]
		// w это вектор nx1, wTQ это вектор 1xn
		TMatrix<T> wwTQ(N, N);
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				wwTQ[i][j] = w[i] * wTQ[0][j];
			}
		}

		// теперь вычитаем из Q полученную матрицу wwTQ
		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				Q[i][j] -= wwTQ[i][j];
			}
		}
	}
	*/

	// Строим Q из сохранённых в R векторов
	for (ptrdiff_t k = N - 2; k >= 0; --k) {
		// 1. Восстанавливаем вектор w
		//    w[k] = 1, w[k+1..N-1] берём из R
		TVector<T> w(N - k);  // только ненулевая часть!
		w[0] = 1.0;  // w[k] = 1

		double w_squared_sum = 0.0;
		for (size_t i = 0; i < N - k - 1; ++i) {
			w[i + 1] = R[k + 1 + i][k];  // берём из сохранённого
			R[k + 1 + i][k] = 0.0; // занулим заодно эти элементы в R
			w_squared_sum += w[i + 1] * w[i + 1];
		}

		// 2. Вычисляем τ = 2 / (wᵀw)
		//    wᵀw = 1² + sum(w[i]²) = 1 + w_squared_sum
		double tau = 2.0 / (1.0 + w_squared_sum);

		// 3. Вычисляем y = τ * (wᵀ * Q[k:N, :])
		//    y — вектор-строка размером N
		TVector<T> y(N, 0.0);

		for (size_t col = 0; col < N; ++col) {
			double dot = 0.0;
			// Умножаем только ненулевую часть w на соответствующие строки Q
			for (size_t i = 0; i < N - k; ++i) {
				dot += w[i] * Q[k + i][col];
			}
			y[col] = tau * dot;
		}

		// 4. Обновляем только Q[k:N, :] = Q[k:N, :] - w * y
		for (size_t i = 0; i < N - k; ++i) {
			for (size_t col = 0; col < N; ++col) {
				Q[k + i][col] -= w[i] * y[col];
			}
		}
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
