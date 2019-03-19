#include <iostream>
#include <stdio.h>
#include <omp.h>
using namespace std;

void InitMatr(int rows,int cols, int** m)
{
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			m[i][j] = int(rand() % 10);
}


void InitKern(double kernel[3][3], int radius, double sigma)
{
	double norm = 0;

	for (int i = -radius; i <= radius; i++)
		for (int j = -radius; j <= radius; j++)
		{
			kernel[i + radius][j + radius] = (exp(-(i * i + j * j) / (2 * sigma * sigma)));
			norm += kernel[i + radius][j + radius];
		}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			kernel[i][j] /= norm;
}

void check(int rows, int cols,int **picture1, int **picture2)
{

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			if (picture1[i][j] - picture2[i][j] != 0)
			{
				cout << "       " << i << "  " << j << endl;
				//cout << result2[i*rows + j] << endl;
				//cout << image[i*rows + j] << endl;
				cout << "false" << endl;
			}
		}
}

void print(int rows, int cols, int **picture1)
{

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << picture1[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << endl;

}

inline int Clamp(int value, int min, int max)
{
	if (value < min)
		return min;

	if (value > max)
		return max;

	return value;
}

void Gauss_seq(int rows, int cols, double kernel[3][3], int **picture1, int **picture2 )
{
	double temp;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			temp = 0.0;
			for (int q = -1; q <= 1; q++)
				for (int l = -1; l <= 1; l++)
				{
					int idX = Clamp(i + q, 0, rows - 1);
					int idY = Clamp(j + l, 0, cols - 1);
					temp += picture1[idX][idY] * kernel[q + 1][l + 1];
				}
			picture2[i][j] = int(temp);
		}
}

void Gauss_par(int rows, int cols, double kernel[3][3], int **picture1, int **picture2)
{
	double temp;
	omp_set_num_threads(2);
#pragma omp parallel for firstprivate(temp)
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			temp = 0.0;
			for (int q = -1; q <= 1; q++)
				for (int l = -1; l <= 1; l++)
				{
					int idX = Clamp(i + q, 0, rows - 1);
					int idY = Clamp(j + l, 0, cols - 1);
					temp += picture1[idX][idY] * kernel[q + 1][l + 1];
				}
			picture2[i][j] = int(temp);
		}
}

int main()
{
	int rows ,cols;
	double st, end;
	double kernel[3][3];

	rows = 10000, cols = 5000;

	int** picture = new int*[rows];
	int** res_seq = new int*[rows];
	int** res_par = new int*[rows];

	for (int i = 0; i < rows; i++)
	{
		picture[i] = new int[cols];
		res_seq[i] = new int[cols];
		res_par[i] = new int[cols];
	}

	InitMatr(rows,cols, picture);
	InitKern(kernel, 1, 1.0);

	st = omp_get_wtime();
	Gauss_seq(rows, cols, kernel, picture, res_seq);
	end = omp_get_wtime();
	cout << "Time seq:" << end - st << std::endl;

	st = omp_get_wtime();
	Gauss_par(rows, cols, kernel, picture, res_par);
	end = omp_get_wtime();
	cout << "Time par:" << end - st << std::endl;

	check(rows, cols, res_seq, res_par);

	system("pause");

	for (int i = 0; i < rows; i++)
	{
		delete[] picture[i];
		delete[] res_seq[i];
		delete[] res_par[i];
	}
	delete[] picture;
	delete[] res_seq;
	delete[] res_par;

	return 0;
}