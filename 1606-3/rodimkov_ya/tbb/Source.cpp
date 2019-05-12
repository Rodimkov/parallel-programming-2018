#include <iostream>
#include <stdio.h>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range2d.h"
#include "tbb/task_scheduler_init.h"
#include <omp.h>
#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;
using namespace tbb;

void InitMatr(int rows, int cols, int** m)
{
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			m[i][j] = int(rand() % 255);
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

void ConvertPicture(int rows, int cols, int** m,Mat src)
{
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			m[i][j] = (src.at<Vec3b>(i, j)[0] + src.at<Vec3b>(i, j)[1] + src.at<Vec3b>(i, j)[2]) / 3;
}

void ConvertMat(int rows, int cols, int** m, Mat src)
{
	cout << rows << ' ' << cols << endl;

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
		
			src.at<Vec3b>(i, j)[0] = m[i][j];
			src.at<Vec3b>(i, j)[1] = m[i][j];
			src.at<Vec3b>(i, j)[2] = m[i][j];
		}
}

void check(int rows, int cols, int **picture1, int **picture2)
{

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			if (picture1[i][j] - picture2[i][j] != 0)
			{
				cout << "       " << i << "  " << j << endl;
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

void Gauss_seq(int rows, int cols, double kernel[3][3], int **picture1, int **picture2)
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

void par_mat(double kernel[3][3], int **picture1, int **picture2, int i , size_t size_rows, size_t size_cols)
{
	for (size_t j = 0; j < size_cols; j++)
	{
		double temp = 0;
		for (int q = -1; q <= 1; q++)
			for (int l = -1; l <= 1; l++)
			{
				int idX = Clamp(i + q, 0, size_rows - 1);
				int idY = Clamp(j + l, 0, size_cols - 1);
				temp += picture1[idX][idY] * kernel[q + 1][l + 1];
			}
		picture2[i][j] = int(temp);
	}
}

void parallel_matrix_multiply(double kernel[3][3], int **picture1, int **picture2, size_t size_rows, size_t size_cols)
{
	parallel_for(blocked_range<size_t>(0, size_rows), [=](const blocked_range<size_t>& r)
	{
		for (size_t i = r.begin(); i != r.end(); ++i)
		{
			par_mat(kernel, picture1, picture2, i, size_rows, size_cols);
		}
	});
}


int main()
{
	int rows, cols;
	double st, end;
	double kernel[3][3];

	Mat src;
	src = imread("1221.jpg");
	rows = src.rows;
	cols = src.cols;
	Mat res_mat = src.clone();

	int** picture = new int*[rows];
	int** res_seq = new int*[rows];
	int** res_par = new int*[rows];

	for (int i = 0; i < rows; i++)
	{
		picture[i] = new int[cols];
		res_seq[i] = new int[cols];
		res_par[i] = new int[cols];
	}

	ConvertPicture(rows, cols, picture, src);

	//ConvertMat(rows, cols, picture, src);

	InitKern(kernel, 1, 1.0);

	st = omp_get_wtime();
	Gauss_seq(rows, cols, kernel, picture, res_seq);
	end = omp_get_wtime();
	cout << "Time seq:" << end - st << std::endl;

	st = omp_get_wtime();
	Gauss_par(rows, cols, kernel, picture, res_seq);
	end = omp_get_wtime();
	cout << "Time par openmp:" << end - st << std::endl;

	st = omp_get_wtime();
	task_scheduler_init init(2);
	parallel_matrix_multiply(kernel, picture, res_par, rows, cols);
	end = omp_get_wtime();
	cout << "Time par TBB:" << end - st << std::endl;
	//print(rows, cols,res_seq);
	ConvertMat(rows, cols, res_seq, res_mat);

	//int WIDTH = 500, int HEIGHT = 500;

	namedWindow("image");
	imshow("image", src);
	namedWindow("image res");
	imshow("image res", res_mat);

	waitKey(0);

	/*
	Mat gray, edge, draw;
	cvtColor(src1, gray, CV_BGR2GRAY);
	Canny(gray, edge, 50, 150, 3);
	edge.convertTo(draw, CV_8U);
	namedWindow("image");
	imshow("image", draw);
	*/
	waitKey(0);

	

	//check(rows, cols, res_seq, res_par);


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