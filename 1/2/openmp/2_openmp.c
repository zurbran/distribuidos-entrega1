#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1
#define ORDENXSUPCOLUMNAS 2
#define ORDENXINFCOLUMNAS 3

#if defined _OPENMP
#include <omp.h>
#endif

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

int T = 1;
int N;

void mulMatrices(double *A, double *B, double *C)
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
				sum += obtenerValorMatrizColumna(A, i, k, N) * obtenerValorMatrizColumna(B, k, j, N);
			}
			asignarValorMatrizFila(C, i, j, N, sum);
		}
	}  
}

void triangularSuperiorPorCuadrada(double *U, double *B, double *C)
{
	double sum;
	int i, j, k;
	#pragma omp parallel for private(sum, j, k) schedule(dynamic, 64)
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			sum = 0.0;
			for(k = i; k < N ; k++)
			{
				sum += obtenerValorMatrizTriaSupFila(U, i, k, N) * obtenerValorMatrizColumna(B, k, j, N);
			}
			asignarValorMatrizFila(C, i, j, N, sum);
		}
	}  
}

void triangularInferiorPorCuadrada(double *L, double *A, double *C)
{
	double sum;
	int i, j, k;
	#pragma omp parallel for private(sum, j, k) schedule(dynamic, 64)
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			sum = 0.0;
			for(k = 0; k < i + 1; k++)
			{
				sum += obtenerValorMatrizTriaInfFila(L, i, k) * obtenerValorMatrizColumna(A, k, j, N);
			}
			asignarValorMatrizFila(C, i, j, N, sum);
		}
	}  
}

double sumarMatriz(double *A, int length)
{
	double sum = 0.0;
	#pragma omp parallel for reduction(+ : sum)
	for(int i = 0; i < length; i++)
	{
		sum += A[i];
	}

	return sum;
}

void escalarPorMatriz(double *A, double *C, double escalar, int length)
{
	#pragma omp parallel for
	for(int i = 0; i < length; i++)
	{
		C[i] = escalar * A[i];
	}
}

void sumarMatrices(double *A, double *B, double *C, int length)
{
	#pragma omp parallel for
	for(int i = 0; i < length; i++)
	{
		C[i] = A[i] + B[i];
	}
}

void filasAColumnas(double *A, double *B)
{
	#pragma omp parallel for
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			asignarValorMatrizColumna(B, i, j, N, obtenerValorMatrizFila(A, i, j, N));
		}
	}
}

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
	double timetick;

	if ((argc != 3) || ((N = atoi(argv[1])) <= 0))
	{
		printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
		exit(1);
	}

	#if defined _OPENMP
    	T= atoi(argv[2]);
		omp_set_num_threads(T);
	#endif

	// Aloca memoria para las matrices
	double *A = (double*)malloc(sizeof(double)*N*N);
	double *B = (double*)malloc(sizeof(double)*N*N);
	double *C = (double*)malloc(sizeof(double)*N*N);
	double *D = (double*)malloc(sizeof(double)*N*N);
	double *E = (double*)malloc(sizeof(double)*N*N);
	double *F = (double*)malloc(sizeof(double)*N*N);
	double *L = (double*)malloc(sizeof(double)*(N*(N+1))/2);
	double *U = (double*)malloc(sizeof(double)*(N*(N+1))/2);

	// Matrices transpuestas
	double *tC = (double*)malloc(sizeof(double)*N*N);
	double *tE = (double*)malloc(sizeof(double)*N*N);
	double *tF = (double*)malloc(sizeof(double)*N*N);
	double *tTAC = (double*)malloc(sizeof(double)*N*N);
	double *tTBE = (double*)malloc(sizeof(double)*N*N);
	double *tTUF = (double*)malloc(sizeof(double)*N*N);

	// Matrices intermedias
	double *TAAC = (double*)malloc(sizeof(double)*N*N);
	double *TAC = (double*)malloc(sizeof(double)*N*N);
    double *ulTAAC = (double*)malloc(sizeof(double)*N*N);
	double *TBE = (double*)malloc(sizeof(double)*N*N);
    double *TLBE = (double*)malloc(sizeof(double)*N*N);
	double *TUF = (double*)malloc(sizeof(double)*N*N);
    double *TDUF = (double*)malloc(sizeof(double)*N*N);
	double *TLBEDUF = (double*)malloc(sizeof(double)*N*N);

	// Resultado
	double *M = (double*)malloc(sizeof(double)*N*N);

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

	// Calcular promedios
	int squareLength = N*N;
	int triaLength = (N*(N+1))/2;
	double bAvg = sumarMatriz(B, squareLength) / squareLength;
	double uAvg = sumarMatriz(U, triaLength) / squareLength;
	double lAvg = sumarMatriz(L, triaLength) / squareLength;
	double ulAvg = uAvg * lAvg;
	// Fin de calcular promedios
	// Paso matriz C a columnas para hacer A*C
	filasAColumnas(C, tC);

	// Multiplicar AtC y guardar en TAC
	mulMatrices(A, tC, TAC);
	// Paso a columna nuevamente
	filasAColumnas(TAC, tTAC);

	// Multiplicar A por tTAC y guardar en TAAC
	mulMatrices(A, tTAC, TAAC);

	// Multiplicar ulAvg por TAAC y almacenar en ulTAAC
	escalarPorMatriz(TAAC, ulTAAC, ulAvg, squareLength);
	// ulTAAC ahora contiene el primer termino ordenado por filas
	// Ordeno a E por columnas
	filasAColumnas(E, tE);

	// Multiplicar B por tE y almacenar en TBE
	mulMatrices(B, tE, TBE);

	filasAColumnas(TBE, tTBE);

	// Multiplicar L por tTBE (BE) y almacenar en TLBE
	triangularInferiorPorCuadrada(L, tTBE, TLBE);
	/// TLBE ahora contiene el segundo termino sin el escalar multiplicado
	// Preparo F para ser multiplicada pasandola a columnas
	filasAColumnas(F,tF);
	// Multiplicar U por tF y almacenar en TUF
	triangularSuperiorPorCuadrada(U, tF, TUF);

	filasAColumnas(TUF,tTUF);

    // Multiplicar D por tTUF y almacenar en TDUF (ordenado por filas)
	mulMatrices(D, tTUF, TDUF);
	// Dado que TLBE y TDUF estan ordenadas por filas se puede sumar como un vector
	sumarMatrices(TLBE,TDUF,TLBEDUF, squareLength);
	// Multiplico el escalar (promedio de B, bAvg) a la matriz resultante de la suma (TLBEDUF)
	escalarPorMatriz(TLBEDUF, M, bAvg, squareLength);
	// M ahora contiene el segundo y ultimo termino
	// Sumar el primer termino haciendo asi que M tenga el resultado final
	sumarMatrices(M, ulTAAC, M, squareLength);
	//M ahora contiene el resutlado

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
