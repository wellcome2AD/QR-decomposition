#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <vector>

#include "QRDecomposition.h"
#include "TMatrix.h"

template <typename T>
double Fnorm(TMatrix<T> m) {
	double res = 0.0;
	for (size_t i = 0; i < m.Size(); ++i) {
		for (size_t j = 0; j < m[0].Size(); ++j) {
			res += m[i][j] * m[i][j];
		}
	}
	return sqrt(res);
}

template <typename T>
void writeMatrixToFile(std::string fileName, TMatrix<T> m) {
	std::ofstream file(fileName);
	for (size_t i = 0; i < m.Size(); ++i) {
		for (size_t j = 0; j < m[0].Size(); ++j) {
			file << m[i][j] << " ";
		}
		file << "; ";
	}
}

void test_with_known_system_1() {
	std::cout << "test size: " << 3 << std::endl;
	TMatrix<float> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} }, Q(3, 3), R(3, 3);
	
	writeMatrixToFile("test_data\\matrixA_" + std::to_string(3) + "_1.txt", A);

	double start = omp_get_wtime();
	QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(3) + "_1.txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(3) + "_1.txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

void test_with_known_system_2() {
	std::cout << "test size: " << 3 << std::endl;
	TMatrix<float> A = { {3, 2, -5}, {2, -1, 3}, {1, 2, -1} }, Q(3, 3), R(3, 3);
	
	writeMatrixToFile("test_data\\matrixA_" + std::to_string(3) + "_2.txt", A);

	double start = omp_get_wtime();
	QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(3) + "_2.txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(3) + "_2.txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

void test_with_known_system_3() {
	std::cout << "test size: " << 3 << std::endl;
	TMatrix<float> A = { { 4, 2, -1 }, { 5, 3, -2 }, { 3, 2, -3 } }, Q(3, 3), R(3, 3);
	
	writeMatrixToFile("test_data\\matrixA_" + std::to_string(3) + "_3.txt", A);

	double start = omp_get_wtime();
	QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(3) + "_3.txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(3) + "_3.txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

void test_with_generated_system(int N) {
	std::cout << "test size: " << N << std::endl;
	TMatrix<float> A = generate_matrix<float>(N, N), Q(N, N), R(N, N);
	writeMatrixToFile("test_data\\matrixA_" + std::to_string(N) + ".txt", A);

	double start = omp_get_wtime();
	QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(N) + ".txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

int main() {
	test_with_known_system_1();
	test_with_known_system_2();
	test_with_known_system_3();
	std::vector<int> sizes{100, 200, 300, 400, 500, 1000, 1500, 2000};
	for (auto&& size : sizes) {
		test_with_generated_system(size);
	}
	std::cout << "all tests passed";
	return 0;
}
