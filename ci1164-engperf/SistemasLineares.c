#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "SistemasLineares.h"

//-------------------------------------------------------------------
//
//  Não alterar esta parte do código.
//
/// Alocaçao de memória
SistLinear_t* alocaSistLinear (unsigned int n)
{
	SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
	if ( SL ) 
	{
		SL->A = (real_t *) malloc(n * n * sizeof(real_t));
		SL->b = (real_t *) malloc(n * sizeof(real_t));
		SL->x = (real_t *) malloc(n * sizeof(real_t));
		if (!(SL->A) || !(SL->b) || !(SL->x))
			liberaSistLinear(SL);
	}

	return (SL);
}

// Liberacao de memória
void liberaSistLinear (SistLinear_t *SL)
{
	free(SL->A);
	free(SL->b);
	free(SL->x);
	free(SL);
}
//-------------------------------------------------------------------

// Inicializa SL
void inicializaSistLinear (SistLinear_t *SL)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  // inicializa a matriz A
  for (unsigned int i=0; i<SL->n; ++i)
     for (unsigned int j=0; j<SL->n; ++j)
     {
        SL->A[i*SL->n+j] = ( (i==j) ? (double)(DIAG<<1) : 1.0 )  * (double)rand() * invRandMax;
     }

  // inicializa vetores b e x
  for (unsigned int i=0; i<SL->n; ++i)
  {
     SL->b[i] = (double)rand() * invRandMax;
     SL->x[i] = 0.0;
  }
}

// Calcula a normaL2 do resíduo
double normaL2Residuo(SistLinear_t *SL)
{

}


// Método de Gauss-Seidel
int gaussSeidel (SistLinear_t *SL, double eps, double *normaL2, unsigned int *iter)
{

	return (0);
}


