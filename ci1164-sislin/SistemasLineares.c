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

  real_t norma = 0;
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


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param pivotamento flag para indicar se o pivotamento parcial deve ser feito (!=0)

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, int pivotamento)
{
  unsigned int n = SL->n;
  for(int k = 0; k < n-1; k++){ // k in [1..n-1]
    //ENCONTRE i >= k tal que aik != 0
    //int i = 0;
    //while(SL->A[index(i,k,n)] == 0 && i <= n)
    //i++;

    //ABORTE se aii == 0, para todo i >= k
    
    //TROQUE a linha k com a linha i
    //while(SL->A[index(i,k,n)] == 0 && i <= n)
    //i++;
    //if(i>n)

    for(int i = k+1; i < n; i++){
      //m = mik = aik/akk
      //bi = bi - m*bk
      for(int j = k+1; j < n; j++){
	//aij = aij - makj
      }
    }       
  }

  /*
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
  

}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro)
{
  while(){
  }
  printf("Não houve convergência!");

}


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
	      SL->A[i*tam+j] = 1.0 / (real_t)(i+j+1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[i*tam+j] = (real_t)rand() * invRandMax;
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
        SL->A[i*tam+i] *= (real_t)tam;
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
      scanf ("%lg", &SL->A[i*n+j]);

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

