#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

#define index(i,j,n) (i*n)+j
#define mymalloc(n,tipo) (tipo*)malloc(n*sizeof(tipo))
#define testmalloc(v) if(v == NULL) abort()


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  unsigned int tam = SL->n;
  real_t *res = mymalloc(tam, real_t);
  testmalloc(res);  
  
  //calcula resíduo
  for(int i = 0; i < tam; i++){
    res[i] = 0;
    for(int j = 0; j < tam; j++){
      res[i] += SL->A[index(i,j,tam)] * x[j];
    }
  }

#ifdef DEBUG_p
  printf("\nResiduo:");
  prnVetor(res, tam);
#endif

  real_t norma = 0;
  real_t diff = 0;
  //calcula norma do vetor ao residuo
  for(int i = 0; i < tam; i++){
    for(int j = 0; j < tam; j++){
      diff = abs(res[i] - SL->b[i]);
      norma += diff * diff;
    }
  }

  norma = sqrt(norma);

  free(res);
  return norma;
}

void backwardSubstitution( SistLinear_t* SL, real_t* x ){
  real_t* A = SL->A;
  real_t* b = SL->b;
  unsigned int n = SL->n;


#ifdef DEBUG_p
  printf("\nOperações BackwardSub:\n");
#endif

  int i,j;
  i = n-1; j = n-1;
  x[ i ] = b[ i ] / A[ index(i,j,n) ]; // termo isolado

#ifdef DEBUG_p
  printf( "x[%d] = b[%d] / A[%d];\n", i, j, index(i,j,n) );
#endif

  for ( i = n-2; i >= 0; i-- ){
    x[ i ] = b[ i ];

#ifdef DEBUG_p
    printf( "x[%d] = b[%d];\n", i, i );
#endif

    for ( j = i+1; j <= n-1; ++j ){
      x[ i ] -= A[ index(i, j, n) ] * x[ j ];

#ifdef DEBUG_p
      printf( "x[%d] -= A[%d] * x[%d];\n", i, index(i,j,n), j );
#endif
   
    }
    x[ i ] /= A[ index(i, i, n) ];

#ifdef DEBUG_p
    printf( "x[%d] /= A[%d];\n", i, index(i,i,n) );
#endif
  }

#ifdef DEBUG_p
  printf("\nx Calculado:\n");
  for ( int i= 0; i < n; i++ )
    printf( "%f\n", x[ i ]);
#endif
}

void fowardSubstitution( SistLinear_t* SL, real_t* x ){
  real_t* A = SL->A;
  real_t* b = SL->b;
  unsigned int n = SL->n;

#ifdef DEBUG_p
  printf("Operações FowardSub:\n");
#endif

  int i,j;
  i = 0; j = 0;
  x[ i ] = b[ i ] / A[ index(i,j,n) ]; // termo isolado

#ifdef DEBUG_p
  printf( "x[%d] = b[%d] / A[%d];\n", i, j, index(i,j,n) );
#endif

  for ( i = 1; i <= (n-1); i++ ){
    x[ i ] = b[ i ];

#ifdef DEBUG_p
    printf( "x[%d] = b[%d];\n", i, i );
#endif

    for ( j = 0; j < i; ++j ){
      x[ i ] -= A[ index(i, j, n) ] * x[ j ];

#ifdef DEBUG_p 
      printf( "x[%d] -= A[%d] * x[%d];\n", i, index(i,j,n), j );
#endif

    }
    x[ i ] /= A[ index(i, i, n) ];

#ifdef DEBUG_p
    printf( "x[%d] /= A[%d];\n", i, index(i,i,n) );
#endif
  }

#ifdef DEBUG_p
  printf("\nx Calculado:\n");
  for ( int i= 0; i < n; i++ )
    printf( "%f\n", x[ i ]);
#endif
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param pivotamento flag para indicar se o pivotamento parcial deve ser feito (!=0)

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, int pivotamento)
{
#ifdef DEBUG_p
  double tempo_exec = timestamp();
  printf("\nOperações ELiminação de Gauss:\n");
#endif
  unsigned int n = SL->n;

  //Cria tabela de pivotamento
  unsigned int* lut;
  lut = mymalloc(n,unsigned int);
  testmalloc(lut);
  for(int k = 0; k < n; k++) lut[k] = k;
  

  //Copia SistLin para modificação:
  SistLinear_t* SL_equiv = alocaSistLinear( n );
  for(int i=0; i < n; ++i)
    for(int j=0; j < n; ++j)
      SL_equiv->A[index(lut[i],lut[j],n)] = SL->A[index(lut[i],lut[j],n)];
  for(int i=0; i < n; ++i)
    SL_equiv->b[lut[i]] = SL->b[lut[i]];

  
  for(int k = 0; k < n-1; k++){ // k in [1..n-1]
    //Faz o pivomento
    //ENCONTRE i >= k tal que aik != 0
        //int i = 0;
        //while(SL->A[index(i,k,n)] == 0 && i <= n)
        //i++;
        //ABORTE se aii == 0, para todo i >= k
    //TROQUE a linha k com a linha i    
        //while(SL->A[index(i,k,n)] == 0 && i <= n)
        //i++;
        //if(i>n)

    //triangularização
    for(int i = k+1; i < n; i++){
      real_t m = SL_equiv->A[index(lut[i],lut[k],n)]/SL_equiv->A[index(lut[k],lut[k],n)];
      
#ifdef DEBUG_p
      printf("m = A_equiv[%d][%d]/A_equiv[%d][%d]\n", lut[i], lut[k], lut[k], lut[k]);
#endif
      // aik -= m * aik
      SL_equiv->A[index(lut[i],lut[k],n)] = 0;
#ifdef DEBUG_p
      printf("A_equiv[%d][%d] = 0;\n", lut[i], lut[k]);
#endif
      for(int j = k+1; j < n; j++){
	SL_equiv->A[index(lut[i],lut[j],n)] -= m*SL_equiv->A[index(lut[k],lut[j],n)];
#ifdef DEBUG_p
	printf("A_equiv[%d][%d] -= m * A_equiv[%d][%d]\n", lut[i], lut[j], lut[k], lut[j]);
#endif
      }
      SL_equiv->b[lut[i]] -= m*SL_equiv->b[lut[k]];
#ifdef DEBUG_p
      printf("b_equiv[%d] -= m * b_equiv[%d]\n", lut[i], lut[k]);
#endif
    }       
  }

  /*// ELIMINACAO GAUSS
    for(int k = 1; k < n-1; k++){
     int w = abs(akk);
     for(int j = k; j < n; j++){
      if(abs(ajk) > w){
       w = abs(ajk);
       r = j;
      }
      //Trocar linhas k e r
      for(int i = k+1; i < n; i++){
       //m = mik = aik/akk;
       //bi = bi - mik * bk;
       for(int j = k+1; j < n; j++){
        //aij = aij - mik*akj;
       }
      }
     }
    }
*/

#ifdef DEBUG_p
  printf("Sistema Triangular Equivalente:\n");
  prnSistLinear(SL_equiv);
#endif

  //RESOLUÇÃO DO SISTEMA TRIANGULAR
  backwardSubstitution(SL_equiv, x);

#ifdef DEBUG_p
  tempo_exec -= timestamp();
  printf("%Tempo de Execucao de eliminacaoGauss + backSubstitution: 10.10lf seg.\n", -(tempo_exec)*1000);
#endif
  
  free(lut);
  liberaSistLinear(SL_equiv);
  return 0;
}

/*!
  \brief Método de Gauss-Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
*/

int gaussJacobi (SistLinear_t *SL, real_t *x, real_t erro)
{
#ifdef DEBUG_p
  double tempo_exec = timestamp();
  printf("\nOperações Gauss-Jacobi:\n");
#endif

  unsigned int n = SL->n;
  
  for(int k = 0; k < MAXIT; ++k){
    for(int i = 0; i < n; ++i){
      x/*next*/[i] = SL->b[i];
      for(int j = 0; j < n; ++j){
	if(i != j){
	  x[i]/*next*/ -= SL->A[index(i,j,n)] * x[j]/*prev*/;
	}
      }
      x/*next*/[i] /= SL->A[index(i,i,n)];
    }
    if()// checar tolerancia
  }

  //checar maximo de iterações

#ifdef DEBUG_p
  tempo_exec -= timestamp();
  printf("%Tempo de Execucao de gaussJacobi: 10.10lf seg.\n", -(tempo_exec)*1000);
#endif
  
  return 0;
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */

/*int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro)
{
#ifdef DEBUG_p
  double tempo_exec = timestamp();
  printf("\nOperações Gauss-Seidel:\n");
#endif

  while(){
  }
  printf("Não houve convergência!");

#ifdef DEBUG_p
  tempo_exec -= timestamp();
  printf("%Tempo de Execucao de eliminacaoGauss + backSubstitution: 10.10lf seg.\n", -(tempo_exec)*1000);
#endif

return 0;
}
*/


// Alocaçao de memória
SistLinear_t* alocaSistLinear (unsigned int tam)
{
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  if ( SL ) {
    SL->A = (real_t *) malloc(tam * tam * sizeof(real_t));
    SL->b = (real_t *) malloc(tam * sizeof(real_t));

    if (!(SL->A) || !(SL->b))
      liberaSistLinear(SL);
  }
  
  SL->n = tam;

  return (SL);
}

// Liberacao de memória
void liberaSistLinear (SistLinear_t *SL)
{
  free(SL->A);
  free(SL->b);
  free(SL);
}

/*!
  \brief Cria coeficientes e termos independentes do SL
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser: comSolucao,
  eqNula, eqProporcional, eqCombLinear, hilbert 
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void inicializaSistLinear (SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max)
{
  unsigned int tam = SL->n;
  // para gerar valores no intervalo [0,coef_max]

 real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);

  // inicializa vetor b
  for (unsigned int i=0; i<tam; ++i) {
    SL->b[i] = (real_t)rand() * invRandMax;
  }
    
  if (tipo == hilbert) {
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[index(i,j,tam)] = 1.0 / (real_t)(i+j+1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[index(i,j,tam)] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % tam;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[nula*tam+j] = 0.0;
      }
      SL->b[nula] = 0.0;
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % tam;
      unsigned int propSrc = (propDst + 1) % tam;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[propDst*tam+j] = SL->A[propSrc*tam+j] * mult;
      }
      SL->b[propDst] = SL->b[propSrc] * mult;
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % tam;
      unsigned int combSrc1 = (combDst + 1) % tam;
      unsigned int combSrc2 = (combDst + 2) % tam;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[combDst*tam+j] = SL->A[combSrc1*tam+j] + SL->A[combSrc2*tam+j];
      }
      SL->b[combDst] = SL->b[combSrc1] + SL->b[combSrc2];
    }
    else if (tipo == diagDominante) {
      // aumenta o expoente dos termos da diagonal principal
      for (unsigned int i=0; i<tam; ++i) {
        SL->A[index(i,i,tam)] *= (real_t)tam;
      }
    }

  }
}


//MUDAR ESSAS FUNÇÕES CASO VOCÊ MUDE real_t!

SistLinear_t *lerSistLinear ()
{
  unsigned int n;
  SistLinear_t *SL;
  
  scanf("%d",&n);

  SL = alocaSistLinear (n);
  
  for(int i=0; i < n; ++i)
    for(int j=0; j < n; ++j)
      scanf ("%lg", &SL->A[index(i,j,n)]);

  for(int i=0; i < n; ++i)
    scanf ("%lg", &SL->b[i]);
  
  return SL;
}


void prnSistLinear (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n\t");
    for(int j=0; j < n; ++j)
      printf ("%10.5lg", SL->A[i*n+j]);
    printf ("   |   %.8lg", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  int i;

  printf ("\n");
  for(i=0; i < n; ++i)
      printf ("%10.10lg ", v[i]);
  printf ("\n\n");

}

