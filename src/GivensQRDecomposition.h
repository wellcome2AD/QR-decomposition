#pragma once

#include <iostream>
#include <cmath>

#include "TMatrix.h"

template <typename T>
void Givens_QR_decomposition(const TMatrix<T>& A, TMatrix<T>& Q, TMatrix<T>& R)
{
	auto n = A.Size(), m = A[0].Size();
	R = A;

	Q = TMatrix<T>(n, n, 0);
	for (auto i = 0; i < n; i++)
	{
		Q[i][i] = 1;
	}

	for (auto j = 0; j < m - 1; j++)
	{ // зануляем весь столбец j
//#pragma omp parallel for num_threads(thread_num)// if(n > 500)
		for (auto i = j + 1; i < n; i++)
		{ // c j+1 строк
			auto sqrt = std::sqrt(R[j][j] * R[j][j] + R[i][j] * R[i][j]);
			double c = R[j][j] / sqrt, s = -R[i][j] / sqrt;
			//#pragma omp parallel for num_threads(thread_num)// if(n > 500)
			for (auto k = j; k < m; k++)
			{ // вращаем только две строки -- i и j
				auto temp = R[j][k] * c - R[i][k] * s;
				R[i][k] = R[j][k] * s + R[i][k] * c;
				R[j][k] = temp;
			}

			//#pragma omp parallel for num_threads(thread_num)// if(n > 500)
			for (auto k = 0; k < n; k++)
			{ // меняются только две строки -- i и j
				auto temp = c * Q[k][j] - s * Q[k][i];
				Q[k][i] = s * Q[k][j] + c * Q[k][i];
				Q[k][j] = temp;
			}
		}
	}
}

