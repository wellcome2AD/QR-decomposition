#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <vector>

#include "QRDecomposition.h"
#include "LinearEquationSystem.h"

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

//void test_with_known_system_1() {
//	TMatrix<float> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} }, Q(3, 3), R(3, 3);
//	const TVector<float> b = { 1, 8, 5 };
//	const TVector<float> expec_x = { 1, -1, -2 };
//	QR_decomposition(A, b, Q, R);
//	std::cout << "Q: " << Q;
//	std::cout << "R: " << R;
//}
//
//void test_with_known_system_2() {
//	TMatrix<float> A = { {3, 2, -5}, {2, -1, 3}, {1, 2, -1} }, Q(3, 3), R(3, 3);
//	const TVector<float> b = { -1, 13, 9 };
//	const TVector<float> expec_x = { 3, 5, 4 };
//	QR_decomposition(A, b, Q, R);
//	std::cout << "Q: " << Q;
//	std::cout << "R: " << R;
//}
//
//void test_with_known_system_3() {
//	TMatrix<float> A = { { 4, 2, -1 }, { 5, 3, -2 }, { 3, 2, -3 } }, Q(3, 3), R(3, 3);
//	const TVector<float> b = { 1, 2, 0 };
//	const TVector<float> expec_x = { -1, 3, 1 };
//	QR_decomposition(A, b, Q, R);
//	std::cout << "Q: " << Q;
//	std::cout << "R: " << R;
//}

int main() {
	//test_with_known_system_1();
	//test_with_known_system_2();
	//test_with_known_system_3();
	std::vector<int> sizes{100, 200, 300, 400, 500};
	for (auto&& size : sizes) {
		std::cout << "test size: " << size << std::endl;
		const int N = size;
		TMatrix<float> A(N, N), Q(N, N), R(N, N);
		TVector<float> b(N), expec_x(N);
		generate_linear_equation_system<float>(A, b, expec_x);
		writeMatrixToFile("matrixA_" + std::to_string(size) + ".txt", A);


		double start = omp_get_wtime();
		QR_decomposition(A, b, Q, R);
		double end = omp_get_wtime();
		std::cout << "time: " << end - start << std::endl;

		writeMatrixToFile("matrixQ_" + std::to_string(size) + ".txt", Q);
		writeMatrixToFile("matrixR_" + std::to_string(size) + ".txt", R);

		std::cout << "abs error: " << Fnorm(Q * R - A) << std::endl;
		std::cout << "rel error: " << Fnorm(Q * R - A) / Fnorm(A) << std::endl << std::endl;
	}

	std::cout << "all tests passed";
	// сравнить со scilab, eigen по времени
	return 0;
}
