#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1

// Cantidad de threads
int T, N;
double *A, *B, *C;
pthread_barrier_t threadBarrier;

#define obtenerValorMatrizFila(M, F, C, N) (M[(F)*(N)+(C)])
#define asignarValorMatrizFila(M, F, C, N, VALUE) (M[(F)*(N)+(C)] = (VALUE))

#define obtenerValorMatrizColumna(M, F, C, N) (M[(F)+(N)*(C)])
#define asignarValorMatrizColumna(M, F, C, N, VALUE) (M[(F)+(N)*(C)] = (VALUE))

typedef struct threads_args
{
	int start;
	int end;
}

void mulMatrices(double *A, double *B, double *C, int start, int end)
{
	for(int i = start; i <= end; i++)
	{
		for(int j = 0; j < N; j++)
		{
			asignarValorMatrizFila(C, i, j, N, 0);
			for(int k = 0; k < N; k++)
			{
				asignarValorMatrizFila(C, i, j, N, (obtenerValorMatrizFila(C, i, j, N) + obtenerValorMatrizFila(A, i, k, N)*obtenerValorMatrizColumna(B, k, j, N));
			}
		}
	}
}

void filasAColumnas(double *A, double *B, int start, int end)
{
	for(int i = start; i <= end; i++)
	{
		for(int j = 0; j < N; j++)
		{
			asignarValorMatrizColumna(B, i, j, N, obtenerValorMatrizFila(A, i, j, N));
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

void *ejercicioUno(void *args)
{
	threads_args *arg = (threads_args*) args;

	// Copia A en B pero ordenado por columnas
	filasAColumnas(A, B, arg->start, arg->end);

	pthread_barrier_wait(&barrier);
	
	// Realiza la multiplicacion
	mulMatrices(A, B, C, arg->start, arg->end);

}

int main(int argc,char*argv[])
{

	double timetick;


	//Controla los argumentos al programa
	if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0))
	{
		printf("\nUsar: %s n t\n n: Dimension de la matriz\n t: Cantidad de threads", argv[0]);
		exit(1);
	}

	// Aloca memoria para las matrices
	// A sera la matriz a multiplicar por si mismo
	// B una copia de A almacenada por columnas
	A=(double*)malloc(sizeof(double)*N*N);
	B=(double*)malloc(sizeof(double)*N*N);
	C=(double*)malloc(sizeof(double)*N*N);

	//Inicializa las matrices A 1, el resultado sera una matriz con todos sus valores en N
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			setValor(A, i, j, ORDENXFILAS, N, 1);
		}
	}   

	pthread_t threads[T];
	threads_args threadLimits[T];
	pthread_barrier_init(&threadBarrier, NULL, T);

	timetick = dwalltime();

	int extra = N % T;
	for(int i = 0; i < T; i++)
	{
		// Para corregir cuando N no es divisible por T
		if(i < extra)
		{
			threadLimits[i].start = i * (N/T) + i;
			threadLimits[i].end = (i + 1) * (N/T) - 1 + (i + 1); 
		}
		else
		{
			threadLimits[i].start = i * (N/T) + extra;
			threadLimits[i].end = (i + 1) * (N/T) - 1 + extra; 
		}
		pthread_create(&threads[i], NULL, &ejercicioUno, &threadLimits[i]);
	}

	for (int i = 0; i < t; i++)
        pthread_join(threads[i], NULL);

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	// Verifica el resultado
	int correct = 1;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			if(getValor(C, i, j, ORDENXFILAS, N) != N)
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
