#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "mpi.h"
using namespace std;

#define PI 3.14159265358979323846

int main(int argc, char** argv){
	double SUM = PI*PI/6;
	double S_n=0;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//processor 0 should generate vector, and partition and distribute elements to other processors
	//each processor is responsible for summing its own parts
	//at the end, all partial sums should be added together and made available for printout
	//report S-S_n for different values of n
	//I guess this implies that you don't need to do do it similarly to the first problem...?
	if(rank==0)//generate vector
	{
		for (int j=0, j<pow(2,14),j++)
		{ 
			vector_v[j] = 1/(j*j);
		}
	
	}

	for (int i=0, i<12,i++){ //partial sums
		if (i==0)
		{
			for (int j=0, j<pow(2,i+3),j++){
				S_n += vector_v[j];
			}
			printf("%d",SUM-S_n);
		}
		else
		{
			for (int j=pow(2,i+2), j<pow(2,i+3),j++)
			{
				S_n += vector_v[j];
			}
			printf("%d",SUM-S_n);
		}
	}
	return 0;
}
