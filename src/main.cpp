#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <vector>

#include "HouseholdQRDecomposition.h"
#include "TMatrix.h"

typedef float currentType;

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

void HOUSEHOLD_test_with_generated_system(int N, bool writeToFile) {
	std::cout << "test size: " << N << std::endl;
	TMatrix<currentType> A = generate_matrix<currentType>(N, N), Q, R;

	double start = omp_get_wtime();
	Household_QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;
	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;

	if (writeToFile) {
		writeMatrixToFile("test_data\\matrixA_" + std::to_string(N) + ".txt", A);
		writeMatrixToFile("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
}

void HOUSEHOLD_OPTIMIZED_test_with_generated_system(int N, bool writeToFile) {
	std::cout << "test size: " << N << std::endl;
	TMatrix<currentType> A = generate_matrix<currentType>(N, N), Q, R;

	double start = omp_get_wtime();
	Household_QR_decomposition_experimental(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;
	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;

	if (writeToFile) {
		writeMatrixToFile("test_data\\matrixA_" + std::to_string(N) + ".txt", A);
		writeMatrixToFile("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
}

void defaultTest() {
	TMatrix<currentType> A, Q, R, expected_Q, expected_R;

	// 2х2 тест
	//{
	//	A = { {1, 2}, {3, 4} };
	//}
	// 3х3 тест
	{
		A = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
	}

	//// сгенерированная матрица указанного N
	//{
	//	int N = 4;
	//	A = generate_matrix<currentType>(N, N);
	//}

	Household_QR_decomposition(A, expected_Q, expected_R);
	Household_QR_decomposition_experimental(A, Q, R);

	//std::cout << "Q:" << Q << std::endl;
	//std::cout << "expec Q:" << expected_Q << std::endl;

	//std::cout << "R:" << R << std::endl;
	//std::cout << "expec R:" << expected_R << std::endl;

	auto QR = Q * R;
	std::cout << "Q * R:" << QR << std::endl;
	std::cout << "A:" << A << std::endl;
	
	std::cout << "abs error: " << Fnorm(QR - A) << std::endl;
	std::cout << "rel error: " << Fnorm(QR - A) / Fnorm(A) << std::endl << std::endl;
}

int main() {
	//{
	//	std::cout << "---------------Household method---------------\n";
	//	std::vector<int> sizes{ 100, 200, 300, 400, 500 };
	//	for (auto&& size : sizes) {
	//		HOUSEHOLD_test_with_generated_system(size, false);
	//	}
	//}
	{
		std::cout << "---------------Household optimized method---------------\n";
		std::vector<int> sizes{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500 };
		for (auto&& size : sizes) {
			HOUSEHOLD_OPTIMIZED_test_with_generated_system(size, false);
		}
	}
	std::cout << "all tests passed";

	return 0;
}
