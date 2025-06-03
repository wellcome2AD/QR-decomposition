#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <vector>
#include <map>

#include "Householder_method/IHouseholderMethodSolver.h"
#include "Householder_method/HouseholderBasic.h"
#include "Householder_method/HouseholderWithoutMatrixMul.h"
#include "Householder_method/HouseholderWithNormW.h"
#include "TMatrix.h"

typedef double currentType;

template <typename T>
double Fnorm(TMatrix<T> m) {
	double res = 0.0;
	for (size_t i = 0; i < m.GetN(); ++i) {
		for (size_t j = 0; j < m.GetM(); ++j) {
			res += m.get(i, j) * m.get(i, j);
		}
	}
	return sqrt(res);
}

template <typename T>
void writeMatrixToFile(std::string fileName, TMatrix<T> m) {
	std::ofstream file(fileName);
	for (size_t i = 0; i < m.GetN(); ++i) {
		for (size_t j = 0; j < m.GetM(); ++j) {
			file << m.get(i,j) << " ";
		}
		file << "; ";
	}
}

template <typename T>
void HOUSEHOLDER_test_with_generated_system(IHouseholderMethodSolver<T>* solver, int N, bool writeToFile) {
	std::cout << "test size: " << N << std::endl;
	TMatrix<currentType> A = generate_matrix<currentType>(N, N), Q, R;

	double start = omp_get_wtime();
	solver->QR_decomposition(A, Q, R);
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

struct testParams {
	IHouseholderMethodSolver<currentType>* solver;
	std::string description;
	std::vector<int> sizes;
};

int main() {
	auto methods = std::map<int, testParams>{};

	methods[0] = {
		new HouseholderMethodBasic<currentType>(),
		"Householder basic version",
		{ 100, 200, 300, 400, 500 },
	};

	methods[1] = {
		new HouseholderMethodWithoutMatrixMults<currentType>(),
		"Householder without matrix multiplications",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500 },
	};

	methods[2] = {
		new HouseholderMethodWithNormW < currentType>,
		"Householder with normal w in-place",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500 },
	};

	for (const auto& method : methods) {
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes) {
			HOUSEHOLDER_test_with_generated_system(method.second.solver, size, false);
		}
	}
	std::cout << "all tests passed";

	return 0;
}
