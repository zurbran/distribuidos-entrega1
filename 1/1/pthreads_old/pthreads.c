#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#define ORDENXFILAS 0
#define ORDENXCOLUMNAS 1 

typedef struct search
{
    int id;
    int offset;
	double *matrixA;
	double *matrixB;
	double *matrixR;
} search;

//Dimension por defecto de las matrices
int N=100;

//Retorna el valor de la matriz en la posicion fila y columna segun el orden que este ordenada
double getValor(double *matriz,int fila,int columna,int orden){
	if(orden==ORDENXFILAS){
		return(matriz[fila*N+columna]);
	}else{
		return(matriz[fila+columna*N]);
	}
}

//Establece el valor de la matriz en la posicion fila y columna segun el orden que este ordenada
void setValor(double *matriz,int fila,int columna,int orden,double valor){
	if(orden==ORDENXFILAS){
		matriz[fila*N+columna]=valor;
	}else{
		matriz[fila+columna*N]=valor;
	}
}

void *transpose(void *void_ptr)
{
    search *local = (search *) void_ptr;

    //Identifico la porcion que trabajara el thread
    int begin= local->id * local->offset;
    local->id++;
    
    for(int s= begin; s < (begin + local->offset); s++){
        for(int t= 0;t < N; t++){
            setValor(local->matrixB,s,t,ORDENXCOLUMNAS, getValor(local->matrixA, s, t, ORDENXFILAS));
        }
    }
}

void *multiply(void *void_ptr)
{
    search *local = (search *) void_ptr;

    //Identifico la porcion que trabajara el thread
    int begin= local->id * local->offset;
    local->id++;

    for(int s= begin; s < (begin + local->offset); s++)
    {
		for(int t= 0; t<N; t++){
			setValor(local->matrixR,s,t,ORDENXFILAS,0);
			for(int u= 0; u<N; u++){
				setValor(local->matrixR,s,t,ORDENXFILAS, getValor(local->matrixR,s,t,ORDENXFILAS) + getValor(local->matrixA,s,u,ORDENXFILAS)*getValor(local->matrixB,u,t,ORDENXCOLUMNAS));
			}
		}
    }
}

//Para calcular tiempo
double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

int main(int argc,char*argv[]){
	double *A,*B,*C;
	int i,j,k, T;
	int check=1;
	double timetick;
	search args;

	//Controla los argumentos al programa
	if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((atoi(argv[2]) != 2)&&(atoi(argv[2])!=4)))
	{
		printf("\nUsar: %s n\n  n: Dimension de la matriz (nxn X nxn)\n", argv[0]);
		exit(1);
	}

    T= atoi(argv[2]);

    pthread_t thread[T];

	// Aloca memoria para las matrices
	// A sera la matriz a multiplicar por si mismo
	// B una copia de A almacenada por columnas

	A=(double*)malloc(sizeof(double)*N*N);
	B=(double*)malloc(sizeof(double)*N*N);
	C=(double*)malloc(sizeof(double)*N*N);

	//Inicializa las matrices A 1, el resultado sera una matriz con todos sus valores en N
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			setValor(A,i,j,ORDENXFILAS,1);
            
		}
	}

	args.matrixA= A;
	args.matrixB= B;
	args.matrixR= C;
    args.id= 0;

	timetick = dwalltime();

	args.offset = N/T;

	// Copia A en B pero ordenado por columnas

	// for(i=0;i<N;i++){
	// 	for(j=0;j<N;j++){
	// 		setValor(B,i,j,ORDENXCOLUMNAS, getValor(A, i, j, ORDENXFILAS));
	// 	}
	// }

    for(int v = 0; v < T; v++)
    {
        if(pthread_create(&thread[v], NULL , transpose, (void *) &args)) 
        {
            printf("Error creating thread\n");
            return 3;
        }
    }

    for(int x = 0; x < T; x++)
    {
        printf("Joineando thread %d\n",x);
        if(pthread_join(thread[x], NULL))
        {
            printf("Error joining thread\n");
            return 4;
        }
        printf("Thread sumado\n");
    }

	// Realiza la multiplicacion

    args.id= 0;

	for(int q = 0; q < T; q++)
    {
        if(pthread_create(&thread[q], NULL , multiply, (void *) &args)) 
        {
            printf("Error creating thread\n");
            return 3;
        }
    }

    for(int w = 0; w < T; w++)
    {
        printf("Joineando thread %d\n",w);
        if(pthread_join(thread[w], NULL))
        {
            printf("Error joining thread\n");
            return 4;
        }
        printf("Thread sumado\n");
    }

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);

	// Verifica el resultado
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			check=check&&(getValor(args.matrixR,i,j,ORDENXFILAS)==N);
            //printf("%lf ",getValor(args.matrixR,i,j,ORDENXFILAS));
		}
	}

	if(check){
		printf("Multiplicacion de matrices resultado correcto\n");
	}else{
		printf("Multiplicacion de matrices resultado erroneo\n");
	}

	free(args.matrixA);
	free(args.matrixB);
	free(args.matrixR);
	return(0);
}
