#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "QRDecomposition.h"
#include "LinearEquationSystem.h"

int main1() {
	const int N = 3;
	TMatrix<float> A = { {1, -2, 1}, {2, 0, -3}, {2, -1, -1} };
	const TMatrix<float> A_copy = A;
	const TVector<float> b = { 1, 8, 5 };
	TVector<float> x = QR_decomposition(A, b);
	//std::cout << "x\n" << x;

	auto counted_b = substitution(A, x);
	if (counted_b == b) {
		std::cout << "x is correct";
	}
	else {
		std::cout << "x is incorrect";
	}

	return 0;
}
