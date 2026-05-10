#pragma once

#include <vector>
#include <set>
#include <iostream>
#include <cmath>

#include "../IQRSolver.h"
#include "math.h"

extern size_t thread_num;

template <typename T>
class GivensMethodConcurrentTasks : public IQRSolver<T> {
private:
	template <typename T>
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

		size_t size() const
		{
			return N;
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

public:
	// версия с параллельным выполнением задач на вращение непересекающихся строк
	virtual void QR_decomposition(const std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& Q, std::vector<std::vector<T>>& R) override
	{
		auto N = A.size();
		QRMatrix<T> qr(A);

		// expected number of rows: 190, real number: 192
		/*rows(16, 8) are excessively rotated 2 times
			rows(18, 9) are excessively rotated 2 times
			all needed rows are rotated*/
		size_t z_max = N - 1;
		size_t j_max = N - 1;
		for (size_t z = 2; z <= z_max; ++z)
		{
			if (z > thread_num)
			{
				// потоки могут одновременно вращать первые z-1 элементов диагонали...
				for (size_t j_begin = 0; j_begin <= j_max; j_begin += z - 1)
				{
					// ...но в этом нет смысла, если элементов осталось очень мало. Не параллелим, если одновременно работать будут меньше половины потоков.
//#pragma omp parallel for if (z_max - z < thread_num / 2)
					for (int j = j_begin; j + 1 < j_begin + z; ++j)
						if (j + z - 1 <= j_max)
							rotate(qr, j + z - 1, j);
				}
			}
			else
			{
				// Потоки могут одновременно вращать элементы диагонали, отстоящие друг от друга на "z" позиций...
				for (size_t j_begin = 0; j_begin < z; ++j_begin)
				{
					// ...но в этом нет смысла, если элементов осталось очень мало. Не параллелим, если одновременно работать будут меньше половины потоков.
					const int count = j_max - z + 1;
					// #pragma omp parallel for if ((j_max - j_begin) / z < thread_num / 2)
					for (int j = j_begin; j <= count; j += z)
						rotate(qr, j + z - 1, j);
				}
			}
		}

		rotate(qr, N - 1, 0);

		R = Q = std::vector<std::vector<T>>(N, std::vector<T>(N));
		// #pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				R[i][j] = qr.RAt(i, j);
				Q[i][j] = qr.QAt(i, j);
			}
		}
	}

private:
	void rotate(QRMatrix<T>& qr, size_t i, size_t j)
	{
		std::cout << i << " " << j << "\n";
		static std::vector<std::pair<size_t, size_t>> rotated;
		static std::set<std::pair<size_t, size_t>> rotated_set;
		rotated.emplace_back(i, j);
		rotated_set.insert(std::make_pair(i, j));
		size_t rotated_size = rotated.size();
		size_t rotated_set_size = rotated_set.size();
		auto N = qr.size();
		double Rjj = qr.RAt(j, j);
		double Rij = qr.RAt(i, j);

		if (std::abs(Rij) < 1e-11 * std::max(std::abs(Rjj), 1.0)) {
			return; // относительные числа
		}

		if (std::abs(Rij) < 1e-11 * std::abs(Rjj)) { // R[i][j] пренебрежимо мал по сравнению с R[j][j]					
			return;
		}

		if (std::abs(Rij) < 1e-11) { // уже 0
			return;
		}

		auto sqrt = std::sqrt(Rjj * Rjj + Rij * Rij);
		if (std::abs(sqrt) < 1e-11) { // для корректировки ошибки игнорируется малый знаменатель
			return;
		}
		double c = Rjj / sqrt, s = -Rij / sqrt;

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
		for (auto k = 0; k < N; k++)
		{
			// меняются только две строки -- i и j
			auto temp = qr.RAt(j, k) * c - qr.RAt(i, k) * s;
			qr.RAt(i, k) = qr.RAt(j, k) * s + qr.RAt(i, k) * c;
			qr.RAt(j, k) = temp;

			// меняются только два столбца -- i и j
			temp = c * qr.QAt(k, j) - s * qr.QAt(k, i);
			qr.QAt(k, i) = s * qr.QAt(k, j) + c * qr.QAt(k, i);
			qr.QAt(k, j) = temp;
		}
	}
};
