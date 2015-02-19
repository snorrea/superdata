#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

int ipow(int base, int power);

int main(int argc, char** argv){
	int rank, size, tag,length = atoi(argv[1]);
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	unsigned long long int n_o_e = ipow(2,length),i;
	double SUM = PI*PI/6, vector_v[n_o_e], S_n=0;
	unsigned long long int partition = n_o_e/(size);
	double my_vector[partition], my_sum=0, my_sums[size];
	clock_t start = clock(), end;
	if(rank==0){
		for (i=0; i< n_o_e; i++){
			vector_v[i] = 1./((i+1)*(i+1));
		}
		printf("# of elements: \t%llu\n",n_o_e);
		printf("# of elements per process: \t%llu\n",partition);
	}
	MPI_Scatter(vector_v,partition,MPI_DOUBLE,my_vector,partition,MPI_DOUBLE,0,MPI_COMM_WORLD);
	/*if(!rank==0){
		printf("Successful scatter from process %d!\n",rank);
	}*/
	for(i=0;i<partition;i++){
		my_sum += my_vector[i];
	}
	MPI_Reduce(&my_sum,&S_n,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(rank==0){
		printf("Difference: \t%f\n",SUM-S_n);
		printf("Ratio: \t\t%f\n",S_n/SUM);
		end = clock();
		double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
		printf("Time taken: \t%f\n",time_spent);
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
