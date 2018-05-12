#define obtenerValorMatrizFila(Matriz, Fila, Columna, N) (Matriz[(Fila)*(N)+(Columna)])
#define asignarValorMatrizFila(Matriz, Fila, Columna, N, X) (Matriz[(Fila)*(N)+(Columna)] = (X))

#define obtenerValorMatrizColumna(Matriz, Fila, Columna, N) (Matriz[(Fila)+(N)*(Columna)])
#define asignarValorMatrizColumna(Matriz, Fila, Columna, N, VALOR) (Matriz[(Fila)+(N)*(Columna)] = (VALOR))

#define obtenerValorMatrizTriaSupColumna(Matriz, Fila, Columna) (Matriz[(Fila)+((Columna)*((Columna)+1))/2])
#define asignarValorMatrizTriaSupColumna(Matriz, Fila, Columna, VALOR) (Matriz[(Fila)+((Columna)*((Columna)+1))/2]= (VALOR))

#define obtenerValorMatrizTriaInfFila(Matriz, Fila, Columna) (Matriz[(Columna)+((Fila)*((Fila) + 1))/2])
#define asignarValorMatrizTriaInfFila(Matriz, Fila, Columna, VALOR) (Matriz[(Columna)+((Fila)*((Fila) + 1))/2]= (VALOR))

#define obtenerValorMatrizTriaSupFila(Matriz,Fila,Columna,N) (Matriz[(Fila)*(N)+(Columna)-((Fila)*((Fila)+1))/2])
#define asignarValorMatrizTriaSupFila(Matriz,Fila,Columna,N,VALOR) (Matriz[(Fila)*(N)+(Columna)-((Fila)*((Fila)+1))/2] = (VALOR))