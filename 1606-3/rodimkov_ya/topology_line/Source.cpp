#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "mpi.h"

using namespace std;

#if true

int main(int argc, char* argv[])
{
	int size, rank;
	int sending;
	int recieving;

	//-------------------------------------------------------------------------------------------
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc < 2)
	{
		sending = 1;
		recieving = size-1;
	}
	else
	{
		sending = atoi(argv[1]);
		recieving = atoi(argv[2]);
	}
	if (size > 1)
	{
		MPI_Status *status = new MPI_Status();
		int direction;
		if (sending < recieving) direction = 1;
		if (sending > recieving) direction = -1;
		if (sending == recieving) direction = 0;

		double t = 0;

		if (rank == sending)
		{
			t = MPI_Wtime();
			MPI_Send(&t, 1, MPI_DOUBLE, rank + direction, NULL, MPI_COMM_WORLD);
		}
		else if (rank == recieving)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, MPI_ANY_SOURCE, NULL, MPI_COMM_WORLD, status);
			std::cout << "my time =  " << MPI_Wtime() - t << endl;
		}
		else if ((rank > recieving && rank < sending) || (rank < recieving && rank > sending))
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, MPI_ANY_SOURCE, NULL, MPI_COMM_WORLD, status);
			MPI_Send(&t, 1, MPI_DOUBLE, rank + direction, NULL, MPI_COMM_WORLD);
		}
	}

	//------------------------------------------------------------------------------

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (size > 1)
	{
		MPI_Status *status = new MPI_Status();
		MPI_Comm comm;
		int dims[2];
		int coords[2];
		int period[2] = { 0,0 };
		dims[0] = 1;
		dims[1] = size;
		MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &comm);
		double t = 0;

		if (rank == sending)
		{
			t = MPI_Wtime();
			MPI_Send(&t, 1, MPI_DOUBLE , recieving, NULL, comm);
		}

		if (rank == recieving)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, MPI_ANY_SOURCE, NULL, comm, status);
			cout << "standart time =  " << MPI_Wtime() - t << endl;
		}
	}
	MPI_Finalize();
	return 0;
}
#elif false

int main(int argc, char* argv[])
{
	int size, rank;
	int sending;
	int recieving;

	//-------------------------------------------------------------------------------------------
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc < 2)
	{

		sending = 0;
		recieving = size - 1;
	}
	else
	{
		sending = atoi(argv[1]);
		recieving = atoi(argv[2]);
	}
	if (size > 1)
	{
		MPI_Status *status = new MPI_Status();
		int direction;
		if (sending < recieving) direction = 1;
		if (sending > recieving) direction = -1;
		if (sending == recieving) direction = 0;

		double t = 0;

		if (rank == sending)
		{
			t = MPI_Wtime();
			MPI_Send(&t, 1, MPI_DOUBLE, rank + direction, NULL, MPI_COMM_WORLD);
		}
		else if (rank == recieving)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, MPI_ANY_SOURCE, NULL, MPI_COMM_WORLD, status);
			std::cout << "my time =  " << MPI_Wtime() - t << endl;
		}
		else if ((rank > recieving && rank < sending) || (rank < recieving && rank > sending))
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, MPI_ANY_SOURCE, NULL, MPI_COMM_WORLD, status);
			MPI_Send(&t, 1, MPI_DOUBLE, rank + direction, NULL, MPI_COMM_WORLD);
		}
	}
}

	//------------------------------------------------------------------------------
#else
int main(int argc, char* argv[])
{
	int size, rank;
	int sending;
	int recieving;

	//-------------------------------------------------------------------------------------------
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc < 2)
	{

		sending = 0;
		recieving = size - 1;
	}
	else
	{
		sending = atoi(argv[1]);
		recieving = atoi(argv[2]);
	}
	if (size > 1)
	{
		MPI_Status *status = new MPI_Status();
		MPI_Comm comm;
		int dims[2];
		int coords[2];
		int period[2] = { 0,0 };
		dims[0] = 1;
		dims[1] = size;
		MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &comm);
		double t = 0;

		if (rank == sending)
		{
			t = MPI_Wtime();
			MPI_Send(&t, 1, MPI_DOUBLE, recieving, NULL, comm);
		}

		if (rank == recieving)
		{
			MPI_Recv(&t, 1, MPI_DOUBLE, MPI_ANY_SOURCE, NULL, comm, status);
			cout << "standart time =  " << MPI_Wtime() - t << endl;
		}
	}
	MPI_Finalize();
	return 0;
}
#endif
/*
int main(int argc, char* argv[])
{
	int size, rank;
	int sending;
	int recieving;

	//-------------------------------------------------------------------------------------------
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc < 2)
	{
		sending = 0;
		recieving = size - 1;
	}
	else
	{
		sending = atoi(argv[1]);
		recieving = atoi(argv[2]);
	}
	if (size > 1)
	{
		MPI_Status *status = new MPI_Status();
		MPI_Comm comm;
		int dims[2];
		int coords[2];
		int period[2] = { 0,0 };
		dims[0] = 1;
		dims[1] = size;

		if (size > 1)
		{
			MPI_Status *status = new MPI_Status();
			MPI_Comm comm;
			int dims[2];
			int coords[2];
			int period[2] = { 0,0 };
			dims[0] = 1;
			dims[1] = size;
			MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &comm);
			char *str = "";
			int length;

			if (rank == sending)
			{
				string message = "I'm happy to see you!";
				length = message.length();
				str = new char[length];
				for (int i = 0; i < length; i++)
				{
					str[i] = message[i];
				}
				MPI_Send(&length, 1, MPI_INT, recieving, 1, MPI_COMM_WORLD);
				MPI_Send(str, length, MPI_CHAR, recieving, 1, MPI_COMM_WORLD);
				cout << "rank = " << rank << " I sent" << endl;
			}

			if (rank == recieving)
			{
				MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
				str = new char[length];
				MPI_Recv(str, length, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
				cout << "rank = " << rank << "  Yes ! I received" << endl;

				for (int i = 0; i < length; i++)
				{
					cout << str[i];
				}
			}
		}
	}
	MPI_Finalize();
	return 0;
}

int main(int argc, char* argv[])
{
	int size, rank;
	int sending;
	int recieving;

	//-------------------------------------------------------------------------------------------
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc < 2)
	{
		sending = 0;
		recieving = size - 1;
	}
	else
	{
		sending = atoi(argv[1]);
		recieving = atoi(argv[2]);
	}
	if (size > 1)
	{
		MPI_Status *status = new MPI_Status();
		MPI_Comm comm;
		int dims[2];
		int coords[2];
		int period[2] = { 0,0 };
		dims[0] = 1;
		dims[1] = size;

		int direction;
		if (sending < recieving) direction = 1;
		if (sending > recieving) direction = -1;
		if (sending == recieving) direction = 0;

		char *str = "";
		int length;

		if (rank == sending)
		{
			string message = "I'm happy to see you!";
			length = message.length();
			str = new char[length];
			for (int i = 0; i < length; i++)
			{
				str[i] = message[i];
			}


			MPI_Send(&length, 1, MPI_INT, rank + direction, 1, MPI_COMM_WORLD);
			MPI_Send(str, length, MPI_CHAR, rank + direction, 1, MPI_COMM_WORLD);
			cout << "rank = " << rank << " I sent" << endl;

		}
		else if (rank == recieving)
		{
			MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
			str = new char[length];
			MPI_Recv(str, length, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
			cout << "rank = " << rank << "  Yes ! I received" << endl;

			for (int i = 0; i < length; i++)
			{
				cout << str[i];
			}
		}
		else if ((rank > recieving && rank < sending) || (rank < recieving && rank > sending))
		{
			MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
			str = new char[length];
			MPI_Recv(str, length, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
			int tmp;
			MPI_Send(&length, 1, MPI_INT, rank + direction, 1, MPI_COMM_WORLD);
			MPI_Send(str, length, MPI_CHAR, rank + direction, 1, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
	return 0;
}*/