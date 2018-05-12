#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#if defined _OPENMP
#include <omp.h>
#endif

// Cantidad de threads
int T = 1;

#define obtenerValorMatrizFila(M, F, C, N) (M[(F)*(N)+(C)])
#define asignarValorMatrizFila(M, F, C, N, VALUE) (M[(F)*(N)+(C)] = (VALUE))

#define obtenerValorMatrizColumna(M, F, C, N) (M[(F)+(N)*(C)])
#define asignarValorMatrizColumna(M, F, C, N, VALUE) (M[(F)+(N)*(C)] = (VALUE))


// C = AB, siendo A por filas, B por columnas, A por filas
void mulMatrices(double *A, double *B, double *C, int N)
{
	double sum;
	int i, j, k;
	#pragma omp parallel for private(sum, j, k)
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			sum = 0.0;
			for(k = 0; k < N; k++)
			{
				sum += obtenerValorMatrizFila(A, i, k, N) * obtenerValorMatrizColumna(B, k, j, N);
			}
			asignarValorMatrizFila(C, i, j, N, sum);
		}
	}  
}

void transponer(double *A, double *B, int N)
{
	int i, j;
	#pragma omp parallel for private (j)
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			B[i+j*N] = A[i*N+j];
		}
	}  
}

//Para calcular tiempo
double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

int main(int argc,char*argv[])
{
	double *A, *B, *C;
	double timetick;
	int N;

	if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0))
	{
		printf("\nUsar: %s n t\n n: Dimension de la matriz\n t: Cantidad de threads", argv[0]);
		exit(1);
	}

#if defined _OPENMP
	omp_set_num_threads(T);
#endif

	// Aloca memoria para las matrices
	// A sera la matriz a multiplicar por si mismo
	// B una copia de A almacenada por columnas
	A=(double*)malloc(sizeof(double)*N*N);
	B=(double*)malloc(sizeof(double)*N*N);
	C=(double*)malloc(sizeof(double)*N*N);

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			asignarValorMatrizFila(A, i, j, N, 1);
		}
	}   

	timetick = dwalltime();

	// Copia A en B pero ordenado por columnas
	transponer(A, B, N);
	
	// Multiplico A*B, siendo B = A ordenado por columnas
	mulMatrices(A, B, C, N);

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	// Verifica el resultado
	int correct = 1;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			if(obtenerValorMatrizFila(C, i, j, N) != N)
			{
				correct = 0;
				break;
			}
		}
	}   

	if(correct)
	{
		printf("Multiplicacion de matrices resultado correcto\n");
	}
	else
	{
		printf("Multiplicacion de matrices resultado erroneo\n");
	}

	free(A);
	free(B);
	free(C);
	return(0);
}
