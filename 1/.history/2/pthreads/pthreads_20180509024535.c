#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1
#define ORDENXSUPCOLUMNAS 2
#define ORDENXINFCOLUMNAS 3


typedef struct vc_sum_arg
{
    int start;
    int end;
	double *A;
	double *B;
	double *C;
} sum;

int T,part = 0;
typedef struct avg_arg {
    double *matrix;
	int start;
	int end;
} average;

typedef struct mulTriangularInfCuadrada_args
{
	int start;
	int end;
	double *B;
	int ordenB;
	double *C;
	int ordenC;
	double *L;
	int ordenL;
	int N;
} mulTriangularInfCuadrada_args;

typedef struct filasAColumnasT_args
{
	int start;
	int end;
	double *A;
	double *B;
	int N;
} filasAColumnasT_args;

typedef struct sc_mul_arg
{
    int start;
    int end;
	double b;
    double *A;
	double *C;
} scalar;

void *calcularPromedio(void *args)
{
    average *argsPointer = args;
    double *matrix = argsPointer->matrix;
    int start = argsPointer->start;
	int end = argsPointer->end;
    double total = 0.0;

    for (int i = start; i <= end ; i++)
    {
        total += matrix[i];
    }

    total = total/(end-start+1);
    double *result = malloc(sizeof(*result));
    *result = total;
    pthread_exit(result);
}

double promedioVector(double *A,int length)
{
	pthread_t threads[T];
	average averageMatrix[T];
	int extra = length % T;
	for(int i = 0; i < T; i++)
	{
		// Para corregir cuando N no es divisible por T
		if(i < extra)
		{
			averageMatrix[i].start = i * (length/T) + i;
			averageMatrix[i].end = (i + 1) * (length/T) - 1 + (i + 1); 
		}
		else
		{
			averageMatrix[i].start = i * (length/T) + extra;
			averageMatrix[i].end = (i + 1) * (length/T) - 1 + extra; 
		}

		averageMatrix[i].matrix = A;

		pthread_create(&threads[i], NULL, calcularPromedio, &averageMatrix[i]);
	}

    double avg = 0.0;
    for(int i = 0; i < T; i++)
    {
        double *retAvg;
        pthread_join(threads[i],(void**)&retAvg);
        avg += *retAvg;
        free(retAvg);
    }
	return avg/(double)T;
}

double getValor(double *matriz, int fila, int columna, int orden, int N)
{
	if(orden == ORDENXFILAS)
	{
		return(matriz[fila*N+columna]);
	}
	else if(orden == ORDENXINFCOLUMNAS)
	{
		return(matriz[fila + columna * N - (columna * (columna + 1)) / 2]);
	}
	else if(orden == ORDENXSUPCOLUMNAS)
	{
		return(matriz[fila + (columna * (columna + 1)) / 2]);
	}
	else
	{
		return(matriz[fila+columna*N]);
	}
}

//Establece el valor de la matriz en la posicion fila y columna segun el orden que este ordenada
void setValor(double *matriz, int fila, int columna, int orden, int N, double valor)
{
	if(orden==ORDENXFILAS)
	{
		matriz[fila*N+columna] = valor;
	}
	else if(orden == ORDENXINFCOLUMNAS)
	{
		matriz[fila + columna * N - (columna * (columna + 1)) / 2] = valor;
	}
	else if(orden == ORDENXSUPCOLUMNAS)
	{
		matriz[fila + (columna * (columna + 1)) / 2] = valor;
	}
	else
	{
		matriz[fila+columna*N] = valor;
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

void *triangularInferiorPorCuadrada_thread(void *args)
{
	mulTriangularInfCuadrada_args *arg = (mulTriangularInfCuadrada_args *)args;
	int start = arg->start;
	int end = arg->end;
	int N = arg->N;
	double *C = arg->C;
	double *B = arg->B;
	double *L = arg->L;
	int ordenB = arg->ordenB;
	int ordenC = arg->ordenC;
	int ordenL = arg->ordenL;

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			setValor(C, i, j, ordenC, N, 0);
			for(int k = 0; k < i + 1; k++)
			{
				setValor(C, i, j, ordenC, N, getValor(C, i, j, ordenC, N) + getValor(L, i, k, ordenL, N)*getValor(B, k, j, ordenB, N));
			}
		}
	}  
}

void triangularInferiorPorCuadrada(double *L, int ordenL, double *B, int ordenB, double *C, int ordenC, int N)
{
	pthread_t threads[T];
	mulTriangularInfCuadrada_args args[T];
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

		args[i].B = B;
		args[i].L = L;
		args[i].C = C;

		args[i].ordenB = ordenB;
		args[i].ordenL = ordenL;
		args[i].ordenC = ordenC;

		args[i].N = N;

		pthread_create(&threads[i], NULL, &triangularInferiorPorCuadrada_thread, &args[i]);
	}
	for(int i = 0; i < T; i++)
	{
		pthread_join(threads[i], NULL);
	}
}

void *multiplyScalar( void *args )
{
    scalar *local = (scalar *) args;

	for(int i= local->start; i <= local->end; i++)
	{
		local->C[i]= local->b * local->A[i];
	}
}

void escalarPorVector(double b, double *A, double *C, int lenght)
{
	pthread_t threads[T];
	scalar sc_args[T];
	int extra = lenght % T;
    
	for(int i = 0; i < T; i++)
	{
		// Para corregir cuando N no es divisible por T
		if(i < extra)
		{
			sc_args[i].start = i * (lenght/T) + i;
			sc_args[i].end = (i + 1) * (lenght/T) - 1 + (i + 1); 
		}
		else
		{
			sc_args[i].start = i * (lenght/T) + extra;
			sc_args[i].end = (i + 1) * (lenght/T) - 1 + extra; 
		}

        sc_args[i].b= b;
	    sc_args[i].A= A;
	    sc_args[i].C= C;

		pthread_create(&threads[i], NULL, multiplyScalar, &sc_args[i]);
	}

	for(int i= 0; i < T; i++)
	{
		pthread_create(&threads[i],NULL, multiplyScalar,&sc_args);
	}

	for(int j = 0; j < T; j++)
    {
        if(pthread_join(threads[j], NULL))
        {
            printf("Error joining thread\n");
        }
    }
}

void *vectorSum( void *args )
{
    sum *local = (sum *) args;

	for(int i= local->start; i <= local->end; i++)
	{
		local->C[i]= local->A[i] + local->B[i];
	}
}

void sumaVectores(double *A, double *B, double *C, int lenght)
{
	pthread_t threads[T];
	int extra = lenght % T;
	sum vc_args[T];

    for(int i = 0; i < T; i++)
	{
		// Para corregir cuando N no es divisible por T
		if(i < extra)
		{
			vc_args[i].start = i * (lenght/T) + i;
			vc_args[i].end = (i + 1) * (lenght/T) - 1 + (i + 1); 
		}
		else
		{
			vc_args[i].start = i * (lenght/T) + extra;
			vc_args[i].end = (i + 1) * (lenght/T) - 1 + extra; 
		}

        vc_args[i].A= A;
	    vc_args[i].B= B;
	    vc_args[i].C= C;

		pthread_create(&threads[i], NULL, vectorSum, &vc_args[i]);
	}

	for(int j = 0; j < T; j++)
    {
        if(pthread_join(threads[j], NULL))
        {
            printf("Error joining thread\n");
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
	int N;
	double timetick;

	//Controla los argumentos al programa
	if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((atoi(argv[2]) != 2)&&(atoi(argv[2])!=4)))
	{
		printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
		exit(1);
	}

    T= atoi(argv[2]);

	// Aloca memoria para las matrices
	double *A = (double*)malloc(sizeof(double)*N*N);
	double *B = (double*)malloc(sizeof(double)*N*N);
	double *C = (double*)malloc(sizeof(double)*N*N);
	double *D = (double*)malloc(sizeof(double)*N*N);
	double *E = (double*)malloc(sizeof(double)*N*N);
	double *F = (double*)malloc(sizeof(double)*N*N);
	double *L = (double*)malloc(sizeof(double)*(N*(N+1))/2);
	double *U = (double*)malloc(sizeof(double)*(N*(N+1))/2);

	// Matrices intermedias
	double *T1A = (double*)malloc(sizeof(double)*N*N);
	double *T1B = (double*)malloc(sizeof(double)*N*N);
	double *T2A = (double*)malloc(sizeof(double)*N*N);
	double *T2B = (double*)malloc(sizeof(double)*N*N);
	double *T3A = (double*)malloc(sizeof(double)*N*N);
	double *T3B = (double*)malloc(sizeof(double)*N*N);

	// Resultado
	double *M = (double*)malloc(sizeof(double)*N*N);

	//Inicializa las matrices en 1
	memset(A, 1, sizeof(double) * N * N);
	memset(B, 1, sizeof(double) * N * N);
	memset(C, 1, sizeof(double) * N * N);
	memset(D, 1, sizeof(double) * N * N);
	memset(E, 1, sizeof(double) * N * N);
	memset(F, 1, sizeof(double) * N * N);
	memset(L, 1, sizeof(double) *(N*(N+1))/2);
	memset(U, 1, sizeof(double) *(N*(N+1))/2);

    for(int i = 0; i < (N*N); i++)
    {
       B[i] = i;
    }
    for(int i = 0; i < (N*(N+1))/2; i++)
    {
       U[i] = i;
       L[i] = i;
    }

	timetick = dwalltime();

	double uAvg = promedioVector(U, (N*(N+1))/2);
	double lAvg = promedioVector(L, (N*(N+1))/2);
	double uAvglAvg = uAvg * lAvg;
	double bAvg = promedioVector(B, N*N);
	

	sumaVectores(B, B, C, N*N);

	double cAvg = promedioVector(C, N*N);

    printf("El promedio de U es: %lf \n",uAvg);
	printf("El promedio de L es: %lf \n",lAvg);
	printf("El promedio de B es: %lf \n",bAvg);
	printf("El promedio de C es: %lf \n",cAvg);
	

	// cambiarOrdenMatriz(B,T1A,N,bAvgOffset,T,FilaAColumna);
	// cambiarOrdenMatriz(T1A,T2A,N,bAvgOffset,T,ColumnaAFila);

	// for(int i = 0; i < (N*N); i++)
    // {
    //    printf("Matriz B  %lf \n",B[i]);
    // }
	// 	for(int i = 0; i < (N*N); i++)
    // {
    //    printf("Matriz B flipeada %lf \n",T1A[i]);
    // }
	// for(int i = 0; i < (N*N); i++)
    // {
    //    printf("Matriz T1A flipeada %lf \n",T2A[i]);
    // }


	// double uAvg = promedioVector(U, (N*(N+1))/2);
	// double lAvg = promedioVector(L, (N*(N+1))/2);
	// double uAvglAvg = uAvg * lAvg;
	// double bAvg = promedioVector(B, N*N);

	// // Guardar A por columnas en T1A para multiplicar AA
	// filasAColumnas(T1A, A, N);

	// // Multiplicar AA y guardar en T1B
	// mulMatrices(A, ORDENXFILAS, T1A, ORDENXCOLUMNAS, T1B, ORDENXCOLUMNAS, N);

	// // Multiplicar T1B (AA) por C y guardar en T1A
	// mulMatrices(T1B, ORDENXCOLUMNAS, C, ORDENXFILAS, T1A, ORDENXFILAS, N);

	// // Multiplicar T1A por uAvglAvg y almacenar en T1A
	// escalarPorVector(uAvglAvg, T1A, T1A, N*N);

	// // T1A ahora contiene el primer termino (orden x filas)

	// // Multiplicar L por B y almacenar en T2A
	// triangularInferiorPorCuadrada(L, ORDENXINFCOLUMNAS, B, ORDENXFILAS, T2A, ORDENXCOLUMNAS, N);

	// // Multiplicar T2A (LB) por E y almacenar en T2B
	// mulMatrices(T2A, ORDENXCOLUMNAS, E, ORDENXFILAS, T2B, ORDENXFILAS, N);

	// // Multiplicar T2B (LBE) por bAvg y almacenar en T2B
	// escalarPorVector(bAvg, T2B, T2B, N*N);

	// /// T2B ahora contiene el segundo termino (orden x filas)

	// // Multiplicar D por U y almacenar en T3A
	// cuadradaPorTriangularSuperior(D, ORDENXFILAS, U, ORDENXSUPCOLUMNAS, T3A, ORDENXCOLUMNAS, N);

	// // Multiplicar T3A (DU) y almacenar en T3B
	// mulMatrices(T3A, ORDENXCOLUMNAS, F, ORDENXFILAS, T3B, ORDENXFILAS, N);

	// // Multiplicar T3B (DUF) por bAvg y almacenar en M (orden x filas)
	// escalarPorVector(bAvg, T3B, M, N*N);

	// // M ahora contiene el tercer termino

	// // Sumar como vector M y T2B (posible porque ambas son x filas)
	// sumaVectores(M, T2B, M, N*N);

	// // M ahora contiene la suma del tercer y segundo termino
	
	// // Sumar como vector M y T1A
	// sumaVectores(M, T1A, M, N*N);

	// // M ahora contiene el resutlado


	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	free(L);
	free(U);
	free(T1A);
	free(T1B);
	free(T2A);
	free(T2B);
	free(T3A);
	free(T3B);
	free(M);
	return 0;
}

