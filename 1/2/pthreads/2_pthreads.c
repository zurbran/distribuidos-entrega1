#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#define obtenerValorMatrizFila(M, F, C, N) (M[(F)*(N)+(C)])
#define asignarValorMatrizFila(M, F, C, N, VALOR) (M[(F)*(N)+(C)] = (VALOR))

#define obtenerValorMatrizColumna(M, F, C, N) (M[(F)+(N)*(C)])
#define asignarValorMatrizColumna(M, F, C, N, VALOR) (M[(F)+(N)*(C)] = (VALOR))

#define obtenerValorMatrizTriaSupColumna(M, F, C) (M[(F)+((C)*((C)+1))/2])
#define asignarValorMatrizTriaSupColumna(M, F, C, VALOR) (M[(F)+((C)*((C)+1))/2]= (VALOR))

#define obtenerValorMatrizTriaInfFila(M, F, C) (M[(C)+((F)*((F) + 1))/2])
#define asignarValorMatrizTriaInfFila(M, F, C, VALOR) (M[(C)+((F)*((F) + 1))/2]= (VALOR))

#define obtenerValorMatrizTriaSupFila(M,F,C,N) (M[(F)*(N)+(C)-((F)*((F)+1))/2])
#define asignarValorMatrizTriaSupFila(M,F,C,N,VALOR) (M[(F)*(N)+(C)-((F)*((F)+1))/2] = (VALOR))

int T,N;
double *A, *B, *C, *D, *E, *F, *U, *L, *M;
double *tC, *tE, *tF, *tTAC, *tTBE, *tTUF;
double *TAC, *TAAC, *ulTAAC, *TBE, *TLBE, *TDUF, *TUF, *TLBEDUF;
double ulAvg, bAvg, lAvg, uAvg = 0.0;
double sumaPromedio = 0.0;
pthread_barrier_t threadBarrier;
pthread_mutex_t sumMutex;

typedef struct threads_args
{
	int start;
	int end;
	int vectStart;
	int vectEnd;
	int	triaStart;
	int	triaEnd;
}threads_args;

void escalarPorMatriz(double *A, double *C, double escalar, int start, int end)
{
	for(int i= start; i < end; i++)
	{
		C[i]= escalar * A[i];
	}
}

void sumarMatrices(double *A, double *B, double *C, int start, int end)
{
	for(int i= start; i < end; i++)
	{
		C[i]= A[i] + B[i];
	}
}

void mulMatrices(double *A, double *B, double *C, int start, int end)
{
	double sum;
	for(int i = start; i < end; i++)
	{
		for(int j = 0; j < N; j++)
		{
			sum = 0.0;
			for(int k = 0; k < N; k++)
			{
				sum += obtenerValorMatrizFila(A, i, k, N)*obtenerValorMatrizColumna(B, k, j, N);
			}
			 asignarValorMatrizFila(C, i, j, N, sum);
		}
	}
}

void filasAColumnas(double *A, double *B, int start, int end)
{
	for(int i = start; i < end; i++)
	{
		for(int j = 0; j < N; j++)
		{
			asignarValorMatrizColumna(B, i, j, N, obtenerValorMatrizFila(A, i, j, N));
		}
	}
}

void sumarMatriz(double *M, int start, int end, pthread_mutex_t *mutex)
{
    double total = 0.0;

    for (int i = start; i < end; i++)
    {
        total += M[i];
    }

	pthread_mutex_lock(mutex);
    sumaPromedio += total;
    pthread_mutex_unlock(mutex);
}

void triangularInferiorPorCuadrada(double *L, double *A, double *C, int start, int end)
{
	double sum;
	for(int i = start; i < end; i++)
	{
		for(int j = 0; j < N; j++)
		{
			sum = 0.0;
			for(int k = 0; k <= i ; k++)
			{
				sum += obtenerValorMatrizTriaInfFila(L, i, k) * obtenerValorMatrizColumna(A, k, j, N);
			}
			asignarValorMatrizFila(C, i, j, N, sum);
		}
	}
}

void triangularSuperiorPorCuadrada(double *U, double *B, double *C,int start, int end)
{
	double sum;
	for(int i = start; i < end; i++)
	{
		for(int j = 0; j < N; j++)
		{
			sum = 0.0;
			for(int k = i; k < N; k++)
			{
				sum +=obtenerValorMatrizTriaSupFila(U, i, k, N)*obtenerValorMatrizColumna(B, k, j, N);
			}
			asignarValorMatrizFila(C, i, j, N, sum);
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

void *ejercicioDos(void *args){

	threads_args *arg = (threads_args*) args;

	// Calcular promedios
	sumarMatriz(U,arg->triaStart,arg->triaEnd,&sumMutex);

	if(pthread_barrier_wait(&threadBarrier) == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		uAvg = sumaPromedio / (N*N);
		sumaPromedio = 0;
	}

	pthread_barrier_wait(&threadBarrier);

	sumarMatriz(L,arg->triaStart,arg->triaEnd,&sumMutex);

	if(pthread_barrier_wait(&threadBarrier) == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		lAvg = sumaPromedio / (N*N);
		sumaPromedio = 0;
	}

	pthread_barrier_wait(&threadBarrier);

	ulAvg = uAvg * lAvg;

	sumarMatriz(B,arg->vectStart,arg->vectEnd,&sumMutex);

	if(pthread_barrier_wait(&threadBarrier) == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		bAvg = sumaPromedio / (N*N);
		sumaPromedio = 0;
	}
	// Fin de calcular promedios
	pthread_barrier_wait(&threadBarrier);

	// Paso matriz C a columnas para hacer A*C
	filasAColumnas(C, tC, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	// Multiplicar AtC y guardar en TAC
	mulMatrices(A, tC, TAC, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);
	// Paso a columna nuevamente
	filasAColumnas(TAC, tTAC, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	// Multiplicar A por tTAC y guardar en TAAC
	mulMatrices(A, tTAC, TAAC, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	// Multiplicar ulAvg por TAAC y almacenar en ulTAAC
	escalarPorMatriz(TAAC, ulTAAC, ulAvg, arg->vectStart, arg->vectEnd);
	// ulTAAC ahora contiene el primer termino ordenado por filas
	pthread_barrier_wait(&threadBarrier);
	// Ordeno a E por columnas
	filasAColumnas(E, tE, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	// Multiplicar B por tE y almacenar en TBE
	mulMatrices(B, tE, TBE, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	filasAColumnas(TBE, tTBE, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	// Multiplicar L por tTBE (BE) y almacenar en TLBE
	triangularInferiorPorCuadrada(L, tTBE, TLBE, arg->start, arg->end);
	/// TLBE ahora contiene el segundo termino sin el escalar multiplicado
	pthread_barrier_wait(&threadBarrier);

	// Preparo F para ser multiplicada pasandola a columnas
	filasAColumnas(F,tF, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	// Multiplicar U por tF y almacenar en TUF
	triangularSuperiorPorCuadrada(U, tF, TUF, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

	filasAColumnas(TUF,tTUF, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);

    // Multiplicar D por tTUF y almacenar en TDUF (ordenado por filas)
	mulMatrices(D, tTUF, TDUF, arg->start, arg->end);

	pthread_barrier_wait(&threadBarrier);
	// Dado que TLBE y TDUF estan ordenadas por filas se puede sumar como un vector
	sumarMatrices(TLBE,TDUF,TLBEDUF, arg->vectStart, arg->vectEnd);

	pthread_barrier_wait(&threadBarrier);
	// Multiplico el escalar (promedio de B, bAvg) a la matriz resultante de la suma (TLBEDUF)
	escalarPorMatriz(TLBEDUF, M, bAvg, arg->vectStart, arg->vectEnd);
	// M ahora contiene el segundo y ultimo termino
	pthread_barrier_wait(&threadBarrier);

	// Sumar el primer termino haciendo asi que M tenga el resultado final
	sumarMatrices(M, ulTAAC, M, arg->vectStart, arg->vectEnd);
	//M ahora contiene el resutlado

}

int main(int argc,char*argv[])
{
	double timetick;

	//Controla los argumentos al programa
	if ((argc != 3) || ((N = atoi(argv[1])) <= 0))
	{
		printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
		exit(1);
	}

    T= atoi(argv[2]);

	pthread_t threads[T];
	threads_args threadLimits[T];
	pthread_barrier_init(&threadBarrier, NULL, T);
	pthread_mutex_init(&sumMutex,NULL);

	// Aloca memoria para las matrices
	A = (double*)malloc(sizeof(double)*N*N);
	B = (double*)malloc(sizeof(double)*N*N);
	C = (double*)malloc(sizeof(double)*N*N);
	D = (double*)malloc(sizeof(double)*N*N);
	E = (double*)malloc(sizeof(double)*N*N);
	F = (double*)malloc(sizeof(double)*N*N);
	L = (double*)malloc(sizeof(double)*(N*(N+1))/2);
	U = (double*)malloc(sizeof(double)*(N*(N+1))/2);

	// Matrices transpuestas
	tC = (double*)malloc(sizeof(double)*N*N);
	tE = (double*)malloc(sizeof(double)*N*N);
	tF = (double*)malloc(sizeof(double)*N*N);
	tTAC = (double*)malloc(sizeof(double)*N*N);
	tTBE = (double*)malloc(sizeof(double)*N*N);
	tTUF = (double*)malloc(sizeof(double)*N*N);

	// Matrices intermedias
	TAAC = (double*)malloc(sizeof(double)*N*N);
	TAC = (double*)malloc(sizeof(double)*N*N);
    ulTAAC = (double*)malloc(sizeof(double)*N*N);
	TBE = (double*)malloc(sizeof(double)*N*N);
    TLBE = (double*)malloc(sizeof(double)*N*N);
	TUF = (double*)malloc(sizeof(double)*N*N);
    TDUF = (double*)malloc(sizeof(double)*N*N);
	TLBEDUF = (double*)malloc(sizeof(double)*N*N);

	// Resultado
	M = (double*)malloc(sizeof(double)*N*N);

    for(int i = 0; i < N*N ; i++)
    {
        A[i] = 1;
        B[i] = 1;
        C[i] = 1;
        D[i] = 1;
        E[i] = 1;
        F[i] = 1;
    }

    for(int i = 0; i < (N*(N+1))/2 ; i++)
    {
        U[i] = 1;
        L[i] = 1;
    }

	timetick = dwalltime();

	for(int i = 0; i < T; i++)
	{
		threadLimits[i].start = i * (N/T) ;
		threadLimits[i].end = (i + 1) * (N/T) ; 

		threadLimits[i].vectStart = i * ((N*N)/T) ;
		threadLimits[i].vectEnd = (i + 1) * ((N*N)/T) ;

		threadLimits[i].triaStart = i * (((N*(N+1))/2)/T) ;
		threadLimits[i].triaEnd = (i + 1) * (((N*(N+1))/2)/T) ; 

		pthread_create(&threads[i], NULL, &ejercicioDos, &threadLimits[i]);
	}
	
	for(int i = 0; i < T ; i++)
	{
		pthread_join(threads[i],NULL);
	}

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(L);
	free(U);
	free(tC);
	free(tE);
	free(tF);
	free(tTBE);
	free(tTUF);
	free(tTAC);
	free(TAAC);
	free(TAC);
	free(TBE);
	free(ulTAAC);
	free(TLBE);
	free(TDUF);
    free(TUF);
    free(TLBEDUF);
	free(M);
	return 0;
}