#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <map>

#include <sstream>

#include "IQRSolver.h"

#include "Householder_method/HouseholderBasic.h"
#include "Householder_method/HouseholderWithoutMatrixMul.h"
#include "Householder_method/HouseholderWithNormW.h"

#include "Givens_method/GivensBasic.h"
#include "Givens_Method/GivensQRInOneMatrix.h"
#include "Givens_Method/GivensConcurrentTasks.h"

#include "matrix.h"

typedef double currentType;
size_t thread_num = 8;

template <typename T>
void printResultWithExpected(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& Q, const std::vector<std::vector<T>>& R)
{
	// проверка на соответствие дл€ малых тестов с помощью базовой реализации
	IQRSolver<currentType>* basicMethod = new HouseholderMethodBasic<currentType>();
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

	std::cout << "Q * R:";
	printMatrix(multiplyMatrix(Q, R));
	std::cout << std::endl;
	std::cout << "A:";
	printMatrix(A);
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

	printResultWithExpected(A, Q, R);

	//std::cout << "Q * R:";
	//printMatrix(multiplyMatrix(Q, R));
	//std::cout << std::endl;
	//std::cout << "A:";
	//printMatrix(A);
	//std::cout << std::endl;

	if (writeToFile) {
		writeMatrixToFile("test_data\\matrixA_" + std::to_string(N) + ".txt", A);
		writeMatrixToFile("test_data\\matrixQ_" + std::to_string(N) + ".txt", Q);
		writeMatrixToFile("test_data\\matrixR_" + std::to_string(N) + ".txt", R);
	}
}

struct testParams
{
	IQRSolver<currentType>* solver;
	std::string description;
	std::vector<int> sizes;
};

void rotate(std::vector<std::vector<double>>& Q, std::vector<std::vector<double>>& R, size_t i, size_t j);

void test_rotations()
{
	std::vector<std::vector<double>> matrix = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };

	int N = matrix.size();

	std::cout << "matrix:";
	printMatrix(matrix);
	std::cout << std::endl;

	{
		std::vector<std::vector<double>> R(matrix), Q(N, std::vector<double>(N, 0));
		for (auto i = 0; i < N; i++)
			Q[i][i] = 1;

		rotate(Q, R, 1, 0);
		rotate(Q, R, 2, 0);
		rotate(Q, R, 2, 1);

		std::cout << "rotations: 10, 20, 21\n";

		std::cout << "Q:";
		printMatrix(Q);
		std::cout << std::endl;

		std::cout << "R:";
		printMatrix(R);
		std::cout << std::endl;
	}

	{
		std::vector<std::vector<double>> R(matrix), Q(N, std::vector<double>(N, 0));
		for (auto i = 0; i < N; i++)
			Q[i][i] = 1;

		rotate(Q, R, 2, 1);
		rotate(Q, R, 1, 0);
		rotate(Q, R, 2, 0);

		std::cout << "rotations: 21, 10, 20\n";

		std::cout << "Q:";
		printMatrix(Q);
		std::cout << std::endl;

		std::cout << "R:";
		printMatrix(R);
		std::cout << std::endl;
	}

	{
		std::vector<std::vector<double>> R(matrix), Q(N, std::vector<double>(N, 0));
		for (auto i = 0; i < N; i++)
			Q[i][i] = 1;

		rotate(Q, R, 1, 0);
		rotate(Q, R, 2, 1);
		rotate(Q, R, 2, 0);

		std::cout << "rotations: 10, 21, 20\n";

		std::cout << "Q:";
		printMatrix(Q);
		std::cout << std::endl;

		std::cout << "R:";
		printMatrix(R);
		std::cout << std::endl;
	}
}

int main()
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

	//methods[3] = {
	//	new GivensMethodBasic<currentType>(),
	//	"Givens basic version",
	//	{ 5/*100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000*/ },
	//};

	//methods[4] = {
	//	new GivensMethodQRInOneMatrix<currentType>(),
	//	"Givens with less memory accesses version",
	//	{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	//};
	/*
	methods[5] = {
		new GivensMethodConcurrentTasks<currentType>(),
		"Givens with concurrent execution of rotations",
		{ 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000 },
	};

	for (const auto& method : methods)
	{
		std::cout << method.second.description << std::endl;
		for (auto&& size : method.second.sizes)
		{
			QR_decomposition_test_with_generated_matrix(method.second.solver, size, false);
		}
	}
	*/

	test_rotations();

	std::cout << "all tests passed";

	return 0;
}

int check_pairs()
{
	// по парам чисел пон€ть, какие элементы под главной диагональю: не были занулены, были занулены несколько раз, какие пары строк были ошибочно вращены.
	// входные данные: длина диагонали N; пары чисел от 0 до N-1, обозначающие строки матрицы, которые вращались в ходе работы метода √ивенса » индекс элемента матрицы, который был занулЄн
	int diag_size = 20;
	std::string pairs = R"(
1 0
3 2
5 4
7 6
9 8
11 10
13 12
15 14
17 16
19 18
2 1
4 3
6 5
8 7
10 9
12 11
14 13
16 15
18 17
2 0
5 3
8 6
11 9
14 12
17 15
3 1
6 4
9 7
12 10
15 13
18 16
4 2
7 5
10 8
13 11
16 14
19 17
3 0
7 4
11 8
15 12
19 16
4 1
8 5
12 9
16 13
5 2
9 6
13 10
17 14
6 3
10 7
14 11
18 15
4 0
9 5
14 10
19 15
5 1
10 6
15 11
6 2
11 7
16 12
7 3
12 8
17 13
8 4
13 9
18 14
5 0
11 6
17 12
6 1
12 7
18 13
7 2
13 8
19 14
8 3
14 9
9 4
15 10
10 5
16 11
6 0
13 7
7 1
14 8
8 2
15 9
9 3
16 10
10 4
17 11
11 5
18 12
12 6
19 13
7 0
15 8
8 1
16 9
9 2
17 10
10 3
18 11
11 4
19 12
12 5
13 6
14 7
8 0
9 1
10 2
11 3
12 4
13 5
14 6
15 7
16 8
17 9
18 10
19 11
9 0
10 1
11 2
12 3
13 4
14 5
15 6
16 7
17 8
18 9
19 10
10 0
11 1
12 2
13 3
14 4
15 5
16 6
17 7
18 8
19 9
11 0
12 1
13 2
14 3
15 4
16 5
17 6
18 7
19 8
12 0
13 1
14 2
15 3
16 4
17 5
18 6
19 7
13 0
14 1
15 2
16 3
17 4
18 5
19 6
14 0
15 1
16 2
17 3
18 4
19 5
15 0
16 1
17 2
18 3
19 4
16 0
17 1
18 2
19 3
17 0
18 1
19 2
18 0
19 1
19 0)";
	std::stringstream ss(pairs);
	std::map<std::pair<int, int>, int> indexes;
	int a, b;
	int count = 0;
	while (ss >> a >> b)
	{
		indexes[{a, b}] += 1;
		count++;
	}

	std::cout << "expected number of rows: " << ((diag_size * diag_size - diag_size) / 2) << ", real number: " << count << std::endl;
	bool allRowsRotated = true;
	for (int i = 1; i < diag_size; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if (indexes[{i, j}] == 0)
			{
				allRowsRotated = false;
				std::cout << "rows (" << i << ", " << j << ") are not rotated\n";
			}

			if (indexes[{i, j}] > 1)
			{
				std::cout << "rows (" << i << ", " << j << ") are excessively rotated " << indexes[{i, j}] << " times\n";
			}

			indexes.erase({ i, j });
			// std::cout << indexes[i].first << " " << indexes[i].second << std::endl;
		}
	}

	for (auto&& [key, val] : indexes)
	{
		std::cout << "rows (" << key.first << ", " << key.second << ") don't to be rotated, they were rotated " << val << " times\n";
	}

	if (allRowsRotated)
	{
		std::cout << "all needed rows are rotated\n";
	}

	return 0;
}

void rotate(std::vector<std::vector<double>>& Q, std::vector<std::vector<double>>& R, size_t i, size_t j)
{
	auto N = R.size();
	double Rjj = R[j][j];
	double Rij = R[i][j];

	if (std::abs(Rij) < 1e-11 * std::max(std::abs(Rjj), 1.0)) {
		return; // относительные числа
	}

	if (std::abs(Rij) < 1e-11 * std::abs(Rjj)) { // R[i][j] пренебрежимо мал по сравнению с R[j][j]					
		return;
	}

	if (std::abs(Rij) < 1e-11) { // уже 0
		return;
	}

	auto sqrt = std::sqrt(Rjj * Rjj + Rij * Rij);
	if (std::abs(sqrt) < 1e-11) { // дл€ корректировки ошибки игнорируетс€ малый знаменатель
		return;
	}
	double c = Rjj / sqrt, s = -Rij / sqrt;

#pragma omp parallel for num_threads(thread_num) if (N >= 1000)
	for (auto k = 0; k < N; k++)
	{
		// мен€ютс€ только две строки -- i и j
		auto temp = R[j][k] * c - R[i][k] * s;
		R[i][k] = R[j][k] * s + R[i][k] * c;
		R[j][k] = temp;

		// мен€ютс€ только два столбца -- i и j
		temp = c * Q[k][j] - s * Q[k][i];
		Q[k][i] = s * Q[k][j] + c * Q[k][i];
		Q[k][j] = temp;
	}
}
