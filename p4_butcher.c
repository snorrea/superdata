#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define PI 3.14159265358979323846

int ipow(int base, int power);

int main(int argc, char** argv){
	int k;
	double ALL_SUMS_SER[14-3];
	double ALL_SUMS_OMP[14-3];
	double ALL_SUMS_MPI[14-3];
	int rank, size, tag;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	for(k=3;k<=14;k++){
		unsigned long long int n_o_e = ipow(2,k),i;
		double SUM = PI*PI/6, vector_v[n_o_e],vector_mpi[n_o_e], S_n=0;
		if(rank==0){
		for (i=0; i< n_o_e; i++){
			vector_v[i] = 1./((i+1)*(i+1));
		}
		for(i=0; i< n_o_e; i++){
			ALL_SUMS_SER[k-3] = vector_v[i];
		}
		omp_set_num_threads(atoi(argv[1]));
		#pragma omp parallel for schedule(static) reduction(+:S_n)
		for(i=0; i< n_o_e; i++){
			ALL_SUMS_OMP[k-3] = vector_v[i];
		}
		}
		//MPI-bit
		
		unsigned long long int partition = n_o_e/(size);
		double my_vector[partition], my_sum=0;
		if(rank==0){
			for (i=0; i< n_o_e; i++){
				vector_mpi[i] = 1./((i+1)*(i+1));
			}
		}
		MPI_Scatter(vector_v,partition,MPI_DOUBLE,my_vector,partition,MPI_DOUBLE,0,MPI_COMM_WORLD);
		for(i=0;i<partition;i++){
			my_sum += my_vector[i];
		}
		MPI_Reduce(&my_sum,&S_n,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if(rank==0){
			ALL_SUMS_MPI[k-3] = S_n;
		}
	}
	if(rank==0){
		printf("Difference of serial and parallel values:\n");
		for(k=3;k<=14;k++){
			printf("For k=%d we get:\n",k);
			printf("\t OMP: %f\n",(ALL_SUMS_SER[k-3]-ALL_SUMS_OMP[k-3]));
			printf("\t MPI: %f\n",(ALL_SUMS_SER[k-3]-ALL_SUMS_MPI[k-3]));
		}
	}
	MPI_Finalize();
	return (EXIT_SUCCESS);
}

int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}
/* Write a program to compute the sum Sn using P processors where P is a
power of 2, and a distributed memory model (MPI). The program should
work as follows: Only processor 0 should be responsible for generating the
vector elements. Processor 0 should partition and distribute the vector elements
evenly among all the processors. Each processor should be responsible
for summing up its own part. At the end, all the partial sums should be
added together and made available on processor 0 for printout. Report the
difference S - Sn in double precision for different values of n.*/
