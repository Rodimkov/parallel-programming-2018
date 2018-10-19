#include <ctime>
#include <iostream>
#include <cstdlib> 
#include <cmath>
#include <mpi.h>

using namespace std;

double function(double x)
{
	return exp(x);
}

int main()
{
	double res = 0, tmp = 0, result = 0, A = 1, B = 2;

	double lambda = double(rand()) / RAND_MAX;

	int ProcNum, ProcRank, N = 100000;


	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	srand(time(NULL) + ProcRank);

	for (int i = ProcRank*(int)(N / ProcNum); i < (ProcRank + 1)*(int)(N / ProcNum); i++)
	{
		double lambda = double(rand()) / RAND_MAX;
		tmp += function(A + lambda*(B - A));
	}

	res += (B - A) / double(N) *tmp;
	
	cout << "res = " << res << endl;

	MPI_Reduce(&res, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		tmp = 0;
		for (int i = (int)(N / ProcNum)*ProcNum; i < N; i++)
		{
			double lambda = double(rand()) / RAND_MAX;
			tmp += function(A + lambda*(B - A));
		}

		result += (B - A) / double(N) *tmp;

		cout << "result = " << result << endl;
		
	}

	MPI_Finalize();
	return 0;
}