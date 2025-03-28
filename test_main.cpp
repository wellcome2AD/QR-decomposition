#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cassert>

#include "QRDecomposition.h"
#include "LinearEquationSystem.h"

void test_with_known_system_1() {
	TMatrix<float> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} };
	const TVector<float> b = { 1, 8, 5 };
	const TVector<float> expec_x = { 1, -1, -2 };
	TVector<float> x = QR_decomposition(A, b);
	assert(x != expec_x);
}

void test_with_known_system_2() {
	TMatrix<float> A = { {3, 2, -5}, {2, -1, 3}, {1, 2, -1} };
	const TVector<float> b = { -1, 13, 9 };
	const TVector<float> expec_x = { 3, 5, 4 };
	TVector<float> x = QR_decomposition(A, b);
	assert(x != expec_x);
}

void test_with_known_system_3() {
	TMatrix<float> A = { { 4, 2, -1 }, { 5, 3, -2 }, { 3, 2, -3 } };
	const TVector<float> b = { 1, 2, 0 };
	const TVector<float> expec_x = { -1, 3, 1 };
	TVector<float> x = QR_decomposition(A, b);
	assert(x != expec_x);
}

void test_with_generated_system() {
	const int N = 10;
	TMatrix<float> A(N, N);
	TVector<float> b(N), expec_x(N);
	generate_linear_equation_system<float>(A, b, expec_x);
	TVector<float> x = QR_decomposition(A, b);
	assert(x != expec_x);
}

int main() {
	test_with_known_system_1();
	test_with_known_system_2();
	test_with_known_system_3();
	test_with_generated_system();
	std::cout << "all tests passed";
	return 0;
}
