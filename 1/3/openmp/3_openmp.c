#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#if defined _OPENMP
#include <omp.h>
#endif

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

int main (int argc,char* argv[])
{ 
	int64_t *A;
	int T;
	int64_t N;

	if ((argc != 3) || ((N = atoll(argv[1])) <= 0) || ((T = atoi(argv[2])) < 1) )	
	{
	    printf("\nUsar: %s n\n  n: Dimension del vector\n T: cantidad threads\n", argv[0]);
	    exit(1);
	}

	A = (int64_t*)malloc(sizeof(int64_t)*N);

	#if defined _OPENMP
		omp_set_num_threads(T);
	#endif

	for (int64_t i = 0; i < N; i++) 
	{
		A[i] = i;
	}

	double timetick = dwalltime();

	int64_t pares = 0;
	#pragma omp parallel for reduction(+:pares)
	for (int64_t i = 0; i < N; i++)
	{
		if(A[i] % 2 == 0)
		{
			pares += 1;
		}
	}
	double tiempo = dwalltime() - timetick;	
	printf("Tiempo en segundos %f \n", tiempo);
	printf("Cantidad de pares: %ld \n", pares);
}
