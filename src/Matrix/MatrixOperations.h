#pragma once

#include <iomanip>
#include <vector>
#include <string>
#include <cassert>

template <typename T>
inline double Fnorm(const std::vector<std::vector<T>>& m) {
	double res = 0.0;
	for (size_t i = 0; i < m.size(); ++i) {
		for (size_t j = 0; j < m[0].size(); ++j) {
			res += m[i][j] * m[i][j];
		}
	}
	return sqrt(res);
}

template <typename T>
inline void writeMatrixToFile(std::string fileName, std::vector<std::vector<T>> m) {
	std::ofstream file(fileName);
	for (size_t i = 0; i < m.size(); ++i) {
		for (size_t j = 0; j < m[0].size(); ++j) {
			file << m[i][j] << " ";
		}
		file << "; ";
	}
}

template <typename T>
inline void printMatrix(std::vector<std::vector<T>> m) {
	std::cout << "[\n  ";
	auto N = m.size();
	auto M = m[0].size();
	for (auto i = 0; i < N; ++i) {
		std::cout << std::left << std::setw(7);
		for (auto j = 0; j < M - 1; ++j) {
			std::cout << m[i][j] << ", ";
		}
		std::cout << m[i][M - 1];
		if (i == N - 1) {
			std::cout << "\n]";
		}
		std::cout << "\n  ";
	}
}

template <typename T>
inline void printFlatMatrix(std::vector<T> m, size_t N, size_t M) {
	std::cout << "[\n  ";
	for (auto i = 0; i < N; ++i) {
		std::cout << std::left << std::setw(7);
		for (auto j = 0; j < M - 1; ++j) {
			std::cout << m[i * M + j] << ", ";
		}
		std::cout << m[i * M + M - 1];
		if (i == N - 1) {
			std::cout << "\n]";
		}
		std::cout << "\n  ";
	}
}

template <typename T>
inline std::vector<std::vector<T>> generateMatrix(int N, int M) {
	std::vector<std::vector<T>> res(N, std::vector<T>(M));
	std::srand(std::time(0));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < M; ++j) {
			res[i][j] = std::rand() % 20 - 10;
		}
	}
	return res;
}

template <typename T>
inline std::vector<std::vector<T>> multiplyMatrix(const std::vector<std::vector<T>>& m1, const std::vector<std::vector<T>>& m2) {
	assert(m1.size() == m2[0].size());
	auto N = m1.size();
	std::vector<std::vector<T>> res(N, std::vector<T>(N, 0));
	int s = N / 4;
#pragma omp parallel for num_threads(thread_num)
	for (int i = 0; i < N; i++) {
		for (int jj = 0; jj < N; jj += s) {
			for (int kk = 0; kk < N; kk += s) {
				for (int j = jj; j < ((jj + s) > N ? N : (jj + s)); j++) {
					double temp = 0.0;
					for (int k = kk; k < ((kk + s) > N ? N : (kk + s)); k++) {
						temp += static_cast<double>(m1[i][k]) * static_cast<double>(m2[k][j]);
					}
					res[i][j] += static_cast<T>(temp);
				}
			}
		}
	}
	return res;
}

template <typename T>
std::vector<std::vector<T>> substractMatrix(const std::vector<std::vector<T>>& m1, const std::vector<std::vector<T>>& m2) {
	std::vector<std::vector<T>> res(m1);
	auto N = m1.size();
	for (auto i = 0; i < N; ++i) {
		for (auto j = 0; j < N; ++j) {
			res[i][j] -= m2[i][j];
		}
	}
	return res;
}

template <typename T>
T randomNotZeroValue()
{
	T value = 0;
	while (value == 0)
	{
		value = std::rand() % 20 - 10;
	}
	return value;
}

template <typename T>
std::vector<std::vector<T>> generateHessenbergMatrix(size_t size, bool upper, bool lower)
{
	std::srand(std::time(0));
	std::vector<std::vector<T>> res(size, std::vector<T>(size, 0));
	size_t N = size, M = size;
	if (upper == true && lower == true)
	{
		for (size_t i = 0; i < N; ++i) {
			if (i != 0) res[i][i - 1] = randomNotZeroValue<T>();
			if (i != N - 1) res[i][i + 1] = randomNotZeroValue<T>();
			res[i][i] = randomNotZeroValue<T>();
		}
	}
	else if (upper == true)
	{
		for (size_t i = 0; i < N; ++i) {
			if (i != 0) res[i][i - 1] = randomNotZeroValue<T>();
			for (size_t j = i; j < M; ++j) {
				res[i][j] = randomNotZeroValue<T>();
			}
		}
	}
	else if (lower == true)
	{
		for (size_t i = 0; i < N; ++i) {
			if (i != N - 1) res[i][i + 1] = randomNotZeroValue<T>();
			for (size_t j = 0; j < i + 1; ++j) {
				res[i][j] = randomNotZeroValue<T>();
			}
		}
	}
	return res;
}
