#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>

#include "invertible_matrix.h"

// Транспонирование матрицы
fmatrix Transpone(fmatrix mas) {
	ptrdiff_t rows = mas.size();
	ptrdiff_t cols = mas[0].size();
	fmatrix res(rows, fvector(cols, 0));
	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < rows; j++) {
			res[i][j] = mas[j][i];
		}
	}
	return res;
}

// Получение матрицы без i-й строки и j-го столбца
// (функция нужна для вычисления определителя и миноров)
fmatrix GetMatr(fmatrix mas, int row, int col) {
	int di, dj;
	ptrdiff_t rows = mas.size();
	ptrdiff_t cols = mas[0].size();
	fmatrix p(rows, fvector(cols, 0));
	di = 0;
	for (int i = 0; i < rows - 1; i++) { // проверка индекса строки
		if (i == row)  // строка совпала с вычеркиваемой
			di = 1;    // - дальше индексы на 1 больше
		dj = 0;
		for (int j = 0; j < cols - 1; j++) { // проверка индекса столбца
			if (j == col)  // столбец совпал с вычеркиваемым
				dj = 1;    // - дальше индексы на 1 больше
			p[i][j] = mas[i + di][j + dj];
		}
	}
	return p;
}

// Рекурсивное вычисление определителя
double Determinant(fmatrix mas, int m) {
	int k;
	ptrdiff_t rows = mas.size();
	fmatrix p(rows, fvector(rows, 0));
	double d = 0;
	k = 1; //(-1) в степени i
	if (m < 1) { printf("Определитель вычислить невозможно!"); return 0; }
	if (m == 1) { d = mas[0][0]; return(d); }
	if (m == 2) { d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]); return(d); }
	if (m > 2) {
		for (int i = 0; i < m; i++) {
			p = GetMatr(mas, i, 0);
			d = d + k * mas[i][0] * Determinant(p, m - 1);
			k = -k;
		}
	}
	return(d);
}

// Обратная матрица
fmatrix Mreverse(fmatrix mas) {
	ptrdiff_t rows = mas.size();
	fmatrix res(rows, fvector(rows, 0));
	double det = Determinant(mas, rows); // находим определитель исходной матрицы
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			res[i][j] = Determinant(GetMatr(mas, i, j), rows - 1);
			if ((i + j) % 2 == 1)       // если сумма индексов строки и столбца нечетная
				res[i][j] = -res[i][j];    // меняем знак минора
			res[i][j] = res[i][j] / det;
		}
	}
	return Transpone(res);
}
