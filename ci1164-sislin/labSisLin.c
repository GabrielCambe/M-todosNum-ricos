#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"
#define index(i,j,n) (i*n)+j
#define mymalloc(n,tipo) (tipo*)malloc(n*sizeof(tipo))
#define testmalloc(v) if(v == NULL) abort()
#define SIST_tam 10

int main ()
{
    // inicializa gerador de nr aleatoreos
    srand(20192);
    
    SistLinear_t *SL = alocaSistLinear( SIST_tam );
    inicializaSistLinear(SL, comSolucao, 485);

    real_t *x = mymalloc(SIST_tam,real_t);
    testmalloc(x);

#ifdef DEBUG_p
    printf("Sistema Original:\n");    
    prnSistLinear( SL );
#endif


    //eliminacaoGauss(SL, x, 0);
    
    for(int t = 0; t < SIST_tam; ++t) x[t] = (real_t)rand(); 
    
    //gaussJacobi(SL, x, 0.5); // tol = 0.000061035

    gaussSeidel(SL, x, EPS);

#ifdef DEBUG_p
    real_t norma = normaL2Residuo(SL, x);
    printf("Norma: %10.10lg\n\n", norma);
#endif

    free(x);
    liberaSistLinear( SL );
}

