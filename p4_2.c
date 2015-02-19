#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

int main(int argc, char** argv){
	omp_set_num_threads(atoi(argv[1]));
	printf("Number of processes running: %d\n",omp_get_num_procs());
	double SUM = PI*PI/6;
	double vector_v[(int) pow(2,14)];
	int i,j;
	for (j=0; j< (int) pow(2,14); j++){ //generate vector
		vector_v[j] = 1./((j+1)*(j+1));
		//printf("%d\n",vector_v[j]==0);
	}
	double S_n;
	clock_t start = clock(), end;
	for (i=0;i<12;i++){ //partial sums
		if (i==0)
		{
			#pragma omp parallel for schedule(static) reduction(+:S_n)
			for (j=0; j<(int) pow(2,i+3);j++){
				S_n += vector_v[j];
			}
			printf("Partial difference with k=3: %f \n",SUM-S_n);
		}
		else
		{
			#pragma omp parallel for schedule(static) reduction(+:S_n)
			for (j=(int)pow(2,i+2); j<(int)pow(2,i+3);j++)
			{
				S_n += vector_v[j];
			}
			//printf("%d",(int)pow(2,i+3)-(int)pow(2,i+2));
			printf("Partial difference with k=%d: %f \n",3+i,SUM-S_n);
		}
	}
	end = clock();
	double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
	//printf("Time taken %d seconds %d milliseconds\n",msec/1000, msec%1000);
	printf("Time taken %f\n",time_spent);
	return 0;
}


