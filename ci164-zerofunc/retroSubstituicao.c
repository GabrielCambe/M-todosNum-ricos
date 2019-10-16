#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef REAL_DBL
typedef double* matriz;
typedef double* vetor;
#elseif
typedef float* matriz;
typedef float* vetor;
#endif

#define mymalloc(n,tipo) (tipo*)malloc(n*sizeof(tipo))
#define testmalloc(v) if(v == NULL) abort()
#define index(i,j,n) ((i*n)+j)
#define TAM_matriz 5
#define MAG_matriz 5
#define MAG_vetor 5
/*
float somaKahanSeq(){
}
*/


/*void imprimeMatrizQuadrada( matriz disp, unsigned int n ){
printf("Elementos da Matriz:\n");
}
*/

void gerarx( vetor x ){
  printf("x:\n");

  for(int i = 0; i < TAM_matriz; i++){
    x[i] = rand() % MAG_vetor;
    printf("%.0lf\n", x[i]);
  }
  printf("\n");  

  return;
}

void gerarM_sup( matriz M ){
  printf("M:\n");

  for(int i = 0; i < TAM_matriz; i++){
    for(int j = 0; j < i; j++){
      M[index(i,j,TAM_matriz)] = 0;
      printf("%.0lf ", M[index(i,j,TAM_matriz)]);
    }

    for(int j = i; j < TAM_matriz; j++){
      M[index(i,j,TAM_matriz)] = rand()  % MAG_matriz;
      printf("%.0lf ", M[index(i,j,TAM_matriz)]);
    }
    printf("\n");
  }
  printf("\n");

  return;
}

void gerarM_inf( matriz M ){
  printf("M:\n");

  for(int i = 0; i < TAM_matriz; i++){
    for(int j = 0; j <= i; j++){
      M[index(i,j,TAM_matriz)] = rand()  % MAG_matriz;;
      printf("%.0lf ", M[index(i,j,TAM_matriz)]);
    }

    for(int j = i+1; j < TAM_matriz; j++){
      M[index(i,j,TAM_matriz)] = 0;
      printf("%.0lf ", M[index(i,j,TAM_matriz)]);
    }
    printf("\n");
  }
  printf("\n");

  return;
}

void gerarb( matriz M, vetor x, vetor b){
  printf("b:\n");

  for(int i = 0; i < TAM_matriz; i++){
    b[i] = 0;
    for(int j = 0; j < TAM_matriz; j++){
      b[i] += M[index(i,j,TAM_matriz)] * x[j];
    }
    printf("%.0lf\n", b[i]);
  }

  printf("\n");
}

void gerarProblema_BckwrdSub( matriz* M, vetor* x, vetor* b ){
  *x = mymalloc(TAM_matriz, double);
  testmalloc(*x);
  gerarx( *x );

  *M = mymalloc(TAM_matriz*TAM_matriz, double);
  testmalloc(*M);
  gerarM_sup( *M );

  *b = mymalloc(TAM_matriz, double);
  testmalloc(*b);
  gerarb( *M, *x, *b);

}

void gerarProblema_FwrdSub( matriz* M, vetor* x, vetor* b ){
  *x = mymalloc(TAM_matriz, double);
  testmalloc(*x);
  gerarx( *x );

  *M = mymalloc(TAM_matriz*TAM_matriz, double);
  testmalloc(*M);
  gerarM_inf( *M );

  *b = mymalloc(TAM_matriz, double);
  testmalloc(*b);
  gerarb( *M, *x, *b);

}

void backwardSubstitution( matriz A, vetor b, unsigned int n, vetor* x_calc ){
  *x_calc = mymalloc(TAM_matriz, double);
  testmalloc(*x_calc);
  vetor x = *x_calc;

  printf("Operações BackwardSub:\n");

  int i,j;
  i = n-1; j = n-1;
  x[ i ] = b[ i ] / A[ index(i,j,n) ]; // termo isolado
  printf( "x[%d] = b[%d] / A[%d];\n", i, j, index(i,j,n) );
  for ( i = n-2; i >= 0; i-- ){
    x[ i ] = b[ i ];
    printf( "x[%d] = b[%d];\n", i, i );
    for ( j = i+1; j <= n-1; ++j ){
      x[ i ] -= A[ index(i, j, n) ] * x[ j ];
      printf( "x[%d] -= A[%d] * x[%d];\n", i, index(i,j,n), j );
    }
    x[ i ] /= A[ index(i, i, n) ];
    printf( "x[%d] /= A[%d];\n", i, index(i,i,n) );
  }

  printf("\nx Calculado:\n");
  for ( int i= 0; i < n; i++ )
    printf( "%f\n", x[ i ]);
}

void fowardSubstitution( matriz A, vetor b, unsigned int n, vetor* x_calc ){
  *x_calc = mymalloc(TAM_matriz, double);
  testmalloc(*x_calc);
  vetor x = *x_calc;

  printf("Operações FowardSub:\n");

  int i,j;
  i = 0; j = 0;
  x[ i ] = b[ i ] / A[ index(i,j,n) ]; // termo isolado
  printf( "x[%d] = b[%d] / A[%d];\n", i, j, index(i,j,n) );
  for ( i = 1; i <= (n-1); i++ ){
    x[ i ] = b[ i ];
    printf( "x[%d] = b[%d];\n", i, i );
    for ( j = 0; j < i; ++j ){
      x[ i ] -= A[ index(i, j, n) ] * x[ j ];
      printf( "x[%d] -= A[%d] * x[%d];\n", i, index(i,j,n), j );
    }
    x[ i ] /= A[ index(i, i, n) ];
    printf( "x[%d] /= A[%d];\n", i, index(i,i,n) );
  }

  printf("\nx Calculado:\n");
  for ( int i= 0; i < n; i++ )
    printf( "%f\n", x[ i ]);
}

int main(){
  srand( time(NULL) );
  vetor x; matriz M; vetor b; vetor x_calc;

#ifdef BACK_SUB  
  gerarProblema_BckwrdSub( &M, &x, &b );
  backwardSubstitution( M, b, TAM_matriz, &x_calc );
#endif

#ifdef FOWR_SUB
  gerarProblema_FwrdSub( &M, &x, &b );
  fowardSubstitution( M, b, TAM_matriz, &x_calc );
#endif

  free(M); free(x); free(b); free(x_calc);
  return 0;
}

//TESTES

//double M[16] ={ 2, 1, 4, 5, 0, 2, 1, 1, 0, 0, 3, 2, 0, 0, 0, 7 };
// X ESPERADO { 2.0, 3.0, 4.0, 1.0 }
//double b[4] = { 28, 11, 14, 7 };

/* INDEX TEST
  for(int i = 0; i < 4; i++)
   for(int j = 0; j < 4; j++)
    printf("[%d][%d] : [%d]\n", i, j, index(i,j,4));
*/
