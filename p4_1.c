#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

int main(int argc, char** argv){
	double SUM = PI*PI/6;
	double vector_v[(int) pow(2,14)];
	int i,j;
	for (j=0; j< (int) pow(2,14); j++){ //generate vector
		vector_v[j] = 1./((j+1)*(j+1));
		//printf("%d\n",vector_v[j]==0);
	}
	double S_n;
	for (i=0;i<12;i++){ //partial sums
		if (i==0)
		{
			for (j=0; j<(int) pow(2,i+3);j++){
				S_n += vector_v[j];
			}
			printf("%f \n",SUM-S_n);
		}
		else
		{
			for (j=(int)pow(2,i+2); j<(int)pow(2,i+3);j++)
			{
				S_n += vector_v[j];
			}
			//printf("%d",(int)pow(2,i+3)-(int)pow(2,i+2));
			printf("%f \n",SUM-S_n);
		}
	}
	printf("Total time taken: %f\n",time()-timestart);
	return 0;
}
