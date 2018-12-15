#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h> 
#include <mpi.h>


void kernel_generation(double kernel[3][3], int radius, double sigma)
{

	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			kernel[i + 1][j + 1] = (double)(exp(-(i*i + j*j) / (sigma*sigma))) / (sigma*sqrt(2 * 3.14));
		}
	}
}

int Clamp(int value, int min, int max)
{
	if (value < min)
		return min;

	if (value > max)
		return max;

	return value;
}

void validation_check(int * res, int * res2,int rows,int cols)
{
	bool flag = false;
	for (int i = 0; i < rows*cols; i++)
		if (res2[i] - res[i] != 0)
		{
			flag = !flag;
			break;
		}
	if (!flag)
		printf("the program worked correctly\n");
}

void sequential_launch(int * res,int * image, double kernel[3][3], int rows, int cols)
{
	double temp;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		{
			temp = 0;
			for (int q = -1; q <= 1; q++)
				for (int l = -1; l <= 1; l++)
				{
					int idX = Clamp(i + q, 0, rows - 1);
					int idY = Clamp(j + l, 0, cols - 1);
					temp += image[idX*rows + idY] * kernel[q + 1][l + 1];
				}
			res[i*cols + j] = Clamp(int(temp), 0, 255);
		}
}


int main(int argc, char* argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int diam = 3;
	int **img;
	double kernel[3][3];
	double sigma = 0;
	int rows = 0, cols = 0;

	if (argc == 4 && atoi(argv[1]) > 100 && atoi(argv[2]) > 100)
	{
		rows = atoi(argv[1]);
		cols = atoi(argv[2]);
		sigma = atof(argv[3]);

	}
	else
	{
		rows = 200;
		cols = 200;
		sigma = 1.0;
	}

	int *image = (int*)malloc(rows*cols * sizeof(int*));
	int *result2 = (int*)malloc(rows*cols * sizeof(int*));

	if (rank == 0)
	{
		img = (int**)malloc(rows * sizeof(int*));
		for (int i = 0; i < rows; ++i)
			img[i] = (int*)malloc(cols * sizeof(int));
		srand(time(NULL));
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				img[i][j] = rand() % 256;
			}

		int k = 0;

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				image[k] = img[i][j];
				k++;
			}
		kernel_generation(kernel, 1, sigma);
	}

	if (rank == 0)
	{
		sequential_launch(result2, image, kernel, rows, cols);
	}
	double start = MPI_Wtime();

	MPI_Bcast(kernel, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Status status;
	MPI_Datatype type, type2;

	int count = ((cols / size) + 2);

	MPI_Type_contiguous(count, MPI_INT, &type2);
	MPI_Type_commit(&type2);

	int stride = sizeof(int) * cols;

	MPI_Type_create_hvector(rows, 1, stride, type2, &type);
	MPI_Type_commit(&type);

	int tmp = cols - ((cols / size)*(size - 2) + (cols / size + 1)) - 1;

	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Send((image + ((tmp - i + 1) + (i - 1)*(cols / size + 1))), 1, type, i, NULL, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(image, 1, type, 0, NULL, MPI_COMM_WORLD, &status);
	}
	int temp = count;
	if (rank == 0)
	{
		temp = tmp + 2;
	}

	int *result = (int*)malloc(rows*temp * sizeof(int*));

	double qwerty;

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < temp; j++)
		{
			qwerty = 0;
			for (int k = -1; k <= 1; k++)
				for (int q = -1; q <= 1; q++)
				{
					int idX = Clamp(i + k, 0, rows - 1);
					int idY = Clamp(j + q, 0, temp - 1);

					qwerty += image[idX*rows + idY] * kernel[1 + k][1 + q];

				}
			result[i*temp + j] = Clamp(qwerty, 0, 255);
		}


	if (rank != 0)
	{
		MPI_Send(result, rows*count, MPI_INT, 0, NULL, MPI_COMM_WORLD);
	}
	if (rank == 0)
	{
		temp--;
		int *res = (int*)malloc(rows*count * sizeof(int*));
		for (int s = 1; s < size; s++)
		{

			MPI_Recv(res, rows*count, MPI_INT, s, NULL, MPI_COMM_WORLD, &status);

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < count; j++)
				{
					image[i*rows + ((j)+temp + (count - 2)*((s - 1)))] = res[i*count + j + 1];
				}
			}
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < temp; j++)
				{
					image[i*rows + j] = result[i*temp + j];
				}
		}
		int tempbus = temp;
		temp++;
		if (size == 1)
		{
			tempbus++;
		}

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < tempbus; j++)
			{
				image[i*rows + j] = result[i*(temp)+j];
			}

		printf("time = %f\n", (MPI_Wtime() - start));
	}
	if (rank == 0)
	{
		validation_check(image, result2, rows, cols);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}