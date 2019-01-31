#include "mpi.h"
#include <iostream>
#include <random>
#include <time.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono> 


// function to compute the 2-norm of a vector v of length n
double norm(double *v, int n) {
	double norm = 0.;

	for (int i = 0; i < n; i++)
		norm += v[i] * v[i];

	return sqrt(norm);
}

// initialise v to values between -10 and 10
void initialize(double *v, int n) {
	for (int i = 0; i < n; i++)
		v[i] = cos(double(i)) * 10.;
}


void normalize_vector(double *v, int n) {
	double norm = 0.;

	// compute the norm of v
	for (int i = 0; i < n; i++)
		norm += v[i] * v[i];
	norm = sqrt(norm);

	// normalize v
	for (int i = 0; i < n; i++)
		v[i] /= norm;
}


int main(int argc, char **argv) {
	const int N = 200000;
	double *v = (double*)malloc(N * sizeof(double));
	bool validated = false;


	//--------------------------------------------------------
	//MPI - Start
	//--------------------------------------------------------
	auto start_p = std::chrono::high_resolution_clock::now();
	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (rank == 0) {
		initialize(v, N);
		auto start_s = std::chrono::high_resolution_clock::now();
		normalize_vector(v, N);
		auto finish_s = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_s = finish_s - start_s;
		std::cout << "Serial time: " << elapsed_s.count() << " s\n";
		// chck the answer
		free(v);
	}


	const int nT = N / 4;

	int vT[N], v1[nT];
	

	if (rank == 0) {
		for (int i = 0; i < N; i++)
		{
			vT[i] = rand() % 10;
		}
	}
	start_p = std::chrono::high_resolution_clock::now();
	MPI_Scatter(&vT[0], nT, MPI_INT, &v1[0], nT, MPI_INT, 0, MPI_COMM_WORLD);

	double norm = 0.;

	// compute the norm of v
	for (int i = 0; i < nT; i++)
		norm += v1[i] * v1[i];

	int a, dest, source;
	for (int i = 1; i < size - 1; i *= 2)
	{
		a = 0;
		dest = rank + i; // получатель
		source = rank - i; // отправитель
		if (dest > size - 1) dest = MPI_PROC_NULL;
		if (source < 0)  source = MPI_PROC_NULL;

		MPI_Send(&norm, 1, MPI_INT, dest, i, MPI_COMM_WORLD);
		MPI_Recv(&a, 1, MPI_INT, source, i, MPI_COMM_WORLD, &status);
		norm = norm + a;
	}
	if (rank == size - 1) {
		norm = sqrt(norm);
	}

	MPI_Scatter(&norm, 1, MPI_INT, &norm, 1, MPI_INT, size - 1, MPI_COMM_WORLD);

	// normalize v
	for (int i = 0; i < nT; i++)
		v1[i] /= norm;

	MPI_Gather(&v1[0], nT, MPI_INT, &vT[0], nT, MPI_INT, 0, MPI_COMM_WORLD);


	if (rank == 0) {
		auto finish_p = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_p = finish_p - start_p;
		std::cout << "Parallel time: " << elapsed_p.count() << " s\n";
		
	}

	MPI_Finalize();



	//--------------------------------------------------------
	//MPI - Finish
	//--------------------------------------------------------



	return 0;
}
