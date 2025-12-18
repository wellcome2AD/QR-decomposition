#pragma once

#include <vector>
#include <string>

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
inline std::vector<std::vector<T>> generate_matrix(int N, int M) {
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
	int s = 50;
#pragma omp parallel for num_threads(thread_num)
	for (int i = 0; i < N; i++) {
		for (int jj = 0; jj < N; jj += s) {
			for (int kk = 0; kk < N; kk += s) {
				for (int j = jj; j < ((jj + s) > N ? N : (jj + s)); j++) {
					T temp = 0;
					for (int k = kk; k < ((kk + s) > N ? N : (kk + s)); k++) {
						temp += m1[i][k] * m2[k][j];
					}
					res[i][j] += temp;
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