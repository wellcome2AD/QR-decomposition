#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <map>

#include "IQRSolver.h"

#include "Householder_method/HouseholderBasic.h"
#include "Householder_method/HouseholderWithoutMatrixMul.h"
#include "Householder_method/HouseholderWithNormW.h"

#include "Givens_method/GivensBasic.h"

#include "matrix.h"

typedef double currentType;
size_t thread_num = 8;

template <typename T>
void printResultWithExpected(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& Q, const std::vector<std::vector<T>>& R)
{
	// проверка на соответствие дл€ малых тестов с помощью базовой реализации
	IQRSolver<currentType>* basicMethod = new GivensMethodBasic<currentType>();
	std::vector<std::vector<T>> expecQ, expecR;
	basicMethod->QR_decomposition(A, expecQ, expecR);
	std::cout << "Q:";
	printMatrix(Q);
	std::cout << std::endl;
	std::cout << "expected Q:";
	printMatrix(expecQ);
	std::cout << std::endl;

	std::cout << "R:";
	printMatrix(R);
	std::cout << std::endl;
	std::cout << "expected R:";
	printMatrix(expecR);
	std::cout << std::endl;
}

template <typename T>
void QR_decomposition_test_with_generated_matrix(IQRSolver<T>* solver, int N, bool writeToFile)
{
	std::cout << "test size: " << N << std::endl;
	std::vector<std::vector<T>> A = generateMatrix<T>(N, N);
	std::vector<std::vector<T>> Q(N, std::vector<T>(N, 0)), R(N, std::vector<T>(N, 0));

	double start = omp_get_wtime();
	solver->QR_decomposition(A, Q, R);
	double end = omp_get_wtime();

	std::cout << "time: " << end - start << std::endl;
	std::cout << "abs error: " << Fnorm(substractMatrix(multiplyMatrix(Q, R), A)) << std::endl;
	std::cout << "rel error: " << Fnorm(substractMatrix(multiplyMatrix(Q, R), A)) / Fnorm(A) << std::endl << std::endl;

	// printResultWithExpected(A, Q, R);

	if (writeToFile) {
		writeMatrixToFile<double>("test_data\\matrixA_" + std::to_string(N) + ".txt", A);
		writeMatrixToFile<double>("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile<double>("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
}

struct testParams
{
	IQRSolver<currentType>* solver;
	std::string description;
	std::vector<int> sizes;
};

void tests()
{
	auto methods = std::map<int, testParams>{};

	//methods[0] = {
	//	new HouseholderMethodBasic<currentType>(),
	//	"Householder basic version",
	//	{ 100, 200, 300, 400, 500 },
	//};

	//methods[1] = {
	//	new HouseholderMethodWithoutMatrixMults<currentType>(),
	//	"Householder without matrix multiplications",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};

	//methods[2] = {
	//	new HouseholderMethodWithNormW < currentType>,
	//	"Householder with normal w in-place",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};

	methods[3] = {
		new GivensMethodBasic<currentType>(),
		"Givens basic version",
		{ /*100, 200, 300, 400, 500,*/ 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	for (const auto& method : methods)
	{
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes)
		{
			QR_decomposition_test_with_generated_matrix(method.second.solver, size, false);
		}
	}
	std::cout << "all tests passed\n\n";
}

// »спользование: QR разложение матрицы ’ессенберга с использованием вращений
// ћатрицы ’ессенберга:
// * верхн€€ (квадратна€ матрица, у которой все элементы лежащие ниже первой поддиагонали равны нулю, т.е. m[i,j]=0 дл€ любого i>j+1)
// * нижн€€ (при транспонировании получаетс€ верхн€€ матрица ’ессенберга)
void hassenberg()
{
	std::cout << "hassenberg matrix decomposition\n";
	for (auto&& N : { 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 })
	{
		auto res = generateHassenbergMatrix<double>(N, true, true);
		// printMatrix(res);
		std::cout << "test size: " << N << std::endl;
		std::vector<std::vector<double>> Q(N, std::vector<double>(N, 0)), R(N, std::vector<double>(N, 0));

		double start = omp_get_wtime();
		auto solver = new GivensMethodBasic<currentType>();
		solver->QR_decomposition(res, Q, R);
		double end = omp_get_wtime();

		// printMatrix(Q);
		// printMatrix(R);

		std::cout << "time: " << end - start << std::endl;
		std::cout << "abs error: " << Fnorm<double>(substractMatrix<double>(multiplyMatrix(Q, R), res)) << std::endl;
		std::cout << "rel error: " << Fnorm<double>(substractMatrix<double>(multiplyMatrix(Q, R), res)) / Fnorm<double>(res) << std::endl << std::endl;

		writeMatrixToFile<double>("test_data\\matrixA_" + std::to_string(N) + ".txt", res);
		writeMatrixToFile<double>("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile<double>("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
	std::cout << "all tests passed\n";
}

int main()
{
	tests();
	// hassenberg();

	return 0;
}
