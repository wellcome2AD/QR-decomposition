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

template<typename T>
bool nearly_equal_abs(T a, T b, T eps = std::numeric_limits<T>::epsilon()) {
	return std::abs(a - b) <= eps;
}

void HOUSEHOLD_test_with_known_system_1() {
	std::cout << "test size: " << 3 << std::endl;
	TMatrix<currentType> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} }, Q(3, 3), R(3, 3);

	writeMatrixToFile("test_data\\matrixA_" + std::to_string(3) + "_1.txt", A);

	double start = omp_get_wtime();
	Household_QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(3) + "_1.txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(3) + "_1.txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

void HOUSEHOLD_test_with_known_system_2() {
	std::cout << "test size: " << 3 << std::endl;
	TMatrix<currentType> A = { {3, 2, -5}, {2, -1, 3}, {1, 2, -1} }, Q(3, 3), R(3, 3);

	writeMatrixToFile("test_data\\matrixA_" + std::to_string(3) + "_2.txt", A);

	double start = omp_get_wtime();
	Household_QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(3) + "_2.txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(3) + "_2.txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

void HOUSEHOLD_test_with_known_system_3() {
	std::cout << "test size: " << 3 << std::endl;
	TMatrix<currentType> A = { { 4, 2, -1 }, { 5, 3, -2 }, { 3, 2, -3 } }, Q(3, 3), R(3, 3);

	writeMatrixToFile("test_data\\matrixA_" + std::to_string(3) + "_3.txt", A);

	double start = omp_get_wtime();
	Household_QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(3) + "_3.txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(3) + "_3.txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

void HOUSEHOLD_test_with_generated_system(int N) {
	std::cout << "test size: " << N << std::endl;
	TMatrix<currentType> A = generate_matrix<currentType>(N, N), Q(N, N), R(N, N);
	writeMatrixToFile("test_data\\matrixA_" + std::to_string(N) + ".txt", A);

	double start = omp_get_wtime();
	Household_QR_decomposition(A, Q, R);
	double end = omp_get_wtime();
	std::cout << "time: " << end - start << std::endl;

	writeMatrixToFile("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
	writeMatrixToFile("test_data\\matrixR_" + std::to_string(N) + ".txt", R);

	std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
	std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
}

int main() {
	///*{
	//	std::cout << "---------------Household method---------------\n";
	//	HOUSEHOLD_test_with_known_system_1();
	//	HOUSEHOLD_test_with_known_system_2();
	//	HOUSEHOLD_test_with_known_system_3();
	//	std::vector<int> sizes{ 100, 200, 300, 400, 500, 1000, 1500, 2000 };
	//	for (auto&& size : sizes) {
	//		HOUSEHOLD_test_with_generated_system(size);
	//	}
	//}*/
	//{
	//	std::cout << "---------------Givens method---------------\n";
	//	std::cout << "non-parallel execution\n";
	//	// GIVENS_test_with_known_system_1();
	//	// GIVENS_test_with_known_system_2();
	//	std::vector<int> sizes{ 100, 200, 300, 400, 500, 1000, 1500, 2000,2500, 3000, 3500};
	//	for (auto&& size : sizes) {
	//		GIVENS_test_with_generated_system(size);
	//	}
	//}
	//std::cout << "all tests passed";
	//return 0;

	TMatrix<currentType> A = { {1, 2}, {3, 4} }, Q, R, expected_Q, expected_R;
	int n = A.Size();
	//TMatrix<float> expected_R = { {-8.12404, -9.60114, -11.0782},
	//						      {-4.97338e-07, 0.904534, 1.80907},
	//						      {-1.92022e-07, 6.70552e-08, 8.9407e-08} },	
	//			   expected_Q = { {-0.123091, 0.904534, 0.408248},
	//							  {-0.492366, 0.301511, -0.816496},
	//							  {-0.86164, -0.301511, 0.408248} };
	
	Household_QR_decomposition(A, expected_Q, expected_R);
	Household_QR_decomposition_experimental(A, Q, R);

	std::cout << "Q:" << Q << std::endl;
	std::cout << "expec Q:" << expected_Q << std::endl;

	//std::cout << "R:" << R << std::endl;
	//std::cout << "expec R:" << expected_R << std::endl;
	//
	//auto QR = Q * R;
	//std::cout << "Q * R:" << QR << std::endl;
	//std::cout << "A:" << A << std::endl;
	//
	//std::cout << "abs error: " << Fnorm(QR - A) << std::endl;
	//std::cout << "rel error: " << Fnorm(QR - A) / Fnorm(A) << std::endl << std::endl;

	return 0;
}
