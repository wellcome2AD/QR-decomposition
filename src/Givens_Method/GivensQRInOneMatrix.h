#pragma once

#include <vector>
#include <iostream>
#include <cmath>

#include "../IQRSolver.h"
#include "math.h"


template <class T>
class QRMatrix
{
public:
	QRMatrix(const std::vector<std::vector<T>>& A)
		: _data(A.size()* A.size())
		, N(A.size())
	{
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i)
		{
			QAt(i, i) = 1.0;
			for (int j = 0; j < N; ++j)
			{
				RAt(i, j) = A[i][j];
			}
		}
	}

	T& RAt(size_t row, size_t col)
	{
		return _data[row * N + col].r;
	}
	T& QAt(size_t row, size_t col) {
		return _data[col * N + row].q;
	}

private:
	struct QRElement {
		T q;
		T r;
	};
	std::vector<QRElement> _data;
	size_t N;
};

template <typename T>
class GivensMethodQRInOneMatrix : public IQRSolver<T> {
public:
	// оптимизированная версия. уменьшены обращения к памяти. QR хранятся в одной структуре, каждый элемент попарно рядом в памяти. Q записана по столбцам, R -- по строкам.
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		auto N = A.size();
		QRMatrix<T> qr(A);

		for (auto j = 0; j < N - 1; j++)
		{ // зануляем весь столбец j
			for (auto i = j + 1; i < N; i++)
			{ // под главной диагональю
				double Rjj = qr.RAt(j, j);
				double Rij = qr.RAt(i, j);

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
						auto temp = qr.RAt(j, k) * c - qr.RAt(i, k) * s;
						qr.RAt(i, k) = qr.RAt(j, k) * s + qr.RAt(i, k) * c;
						qr.RAt(j, k) = temp;
					}

					// меняются только два столбца -- i и j
					auto temp = c * qr.QAt(k, j) - s * qr.QAt(k, i);
					qr.QAt(k, i) = s * qr.QAt(k, j) + c * qr.QAt(k, i);
					qr.QAt(k, j) = temp;
				}
			}
		}

		R = Q = std::vector<std::vector<T>>(N, std::vector<T>(N));
#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				R[i][j] = qr.RAt(i, j);
				Q[i][j] = qr.QAt(i, j);
			}
		}
	}
};
