#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1

// Cantidad de threads
int T = 1;

typedef struct mulMatricesT_args
{
	int start;
	int end;
	double *A;
	int ordenA;
	double *B;
	int ordenB;
	double *C;
	int ordenC;
	int N;
} mulMatricesT_args;

typedef struct filasAColumnasT_args
{
	int start;
	int end;
	double *A;
	double *B;
	int N;
} filasAColumnasT_args;

//Retorna el valor de la matriz en la posicion fila y columna segun el orden que este ordenada
double getValor(double *matriz, int fila, int columna, int orden, int N)
{
	if(orden == ORDENXFILAS)
	{
		return(matriz[fila*N+columna]);
	}
	else
	{
		return(matriz[fila+columna*N]);
	}
}

//Establece el valor de la matriz en la posicion fila y columna segun el orden que este ordenada
void setValor(double *matriz, int fila, int columna, int orden, int N, double valor)
{
	if(orden == ORDENXFILAS)
	{
		matriz[fila*N+columna] = valor;
	}
	else
	{
		matriz[fila+columna*N] = valor;
	}
}

void *mulMatrices_thread(void *arg)
{
	mulMatricesT_args *args = (mulMatricesT_args *)arg;
	for(int i = args->start; i <= args->end; i++)
	{
		for(int j = 0; j < args->N; j++)
		{
			setValor(args->C, i, j, args->ordenC, args->N, 0);
			for(int k = 0; k < args->N; k++)
			{
				setValor(args->C, i, j, args->ordenC, args->N, getValor(args->C, i, j, args->ordenC, args->N) + getValor(args->A, i, k, args->ordenA, args->N)*getValor(args->B, k, j, args->ordenB, args->N));
			}
		}
	}
	return NULL;
}

// C = AB
void mulMatrices(double *A, int ordenA, double *B, int ordenB, double *C, int ordenC, int lenght)
{
	pthread_t threads[T];
	mulMatricesT_args args[T];
	int extra = N % T;
	for(int i = 0; i < T; i++)
	{
		// Para corregir cuando N no es divisible por T
		if(i < extra)
		{
			args[i].start = i * (N/T) + i;
			args[i].end = (i + 1) * (N/T) - 1 + (i + 1); 
		}
		else
		{
			args[i].start = i * (N/T) + extra;
			args[i].end = (i + 1) * (N/T) - 1 + extra; 
		}

		args[i].A = A;
		args[i].B = B;
		args[i].C = C;

		args[i].ordenA = ordenA;
		args[i].ordenB = ordenB;
		args[i].ordenC = ordenC;

		args[i].N = N;

		pthread_create(&threads[i], NULL, &mulMatrices_thread, &args[i]);
	}

	for(int i = 0; i < T; i++)
	{
		pthread_join(threads[i], NULL);
	}
}

void *filasAColumnas_thread(void *arg)
{
	filasAColumnasT_args *args = (filasAColumnasT_args *)arg;
	for(int i = args->start; i <= args->end; i++)
	{
		for(int j = 0; j < args->N; j++)
		{
			setValor(args->B, i, j, ORDENXCOLUMNAS, args->N, getValor(args->A, i, j, ORDENXFILAS, args->N));
		}
	}

	return NULL;
}

void filasAColumnas(double *A, double *B, int N)
{
	pthread_t threads[T];
	filasAColumnasT_args args[T];
	int extra = N % T;
	for(int i = 0; i < T; i++)
	{
		// Para corregir cuando N no es divisible por T
		if(i < extra)
		{
			args[i].start = i * (N/T) + i;
			args[i].end = (i + 1) * (N/T) - 1 + (i + 1); 
		}
		else
		{
			args[i].start = i * (N/T) + extra;
			args[i].end = (i + 1) * (N/T) - 1 + extra; 
		}

		args[i].A = A;
		args[i].B = B;

		args[i].N = N;

		pthread_create(&threads[i], NULL, &filasAColumnas_thread, &args[i]);
	}

	for(int i = 0; i < T; i++)
	{
		pthread_join(threads[i], NULL);
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


	timetick = dwalltime();

	// Copia A en B pero ordenado por columnas
	filasAColumnas(A, B, N);
	
	// Realiza la multiplicacion
	mulMatrices(A, ORDENXFILAS, B, ORDENXCOLUMNAS, C, ORDENXFILAS, N);


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
