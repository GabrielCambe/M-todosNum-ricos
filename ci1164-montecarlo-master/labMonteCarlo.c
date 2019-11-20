#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "utils.h"

#define NRAND    ((double) rand() / RAND_MAX)  // drand48() 
#define SRAND(a) srand(a) // srand48(a)

double f(double* x, unsigned int n){
  double fx = 0;
  for(unsigned int i = 1; i <= n; i++){
    fx += ((x[i]*x[i]*x[i]*x[i]) - (16*x[i]*x[i]) + (5*x[i]))/2;
  }
  return fx;
}

// Integral Monte Carlo da função Styblinski-Tang de 2 variáveis
double styblinskiTang(double a, double b, int namostras)
{
  double resultado;
  double soma = 0.0;
  double x[2];
    
  printf("Metodo de Monte Carlo (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), variaveis = 2\n", a, b, namostras);
  
  double t_inicial = timestamp();
  ////////////////////////////
  for (unsigned int i = 1; i <= namostras; i++){
    x[0] = a + (( double ) rand () / RAND_MAX ) * ( b - a );
    x[1] = a + (( double ) rand () / RAND_MAX ) * ( b - a );
    
    soma += f(x,2);
  }
  resultado = (soma / namostras) * ( b - a ) * ( b - a );
  //////////////////////////
  double t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}


double retangulos_xy(double a, double b, int npontos) {
  const double h = (b-a)/npontos;
  double resultado;
  double soma = 0;
  double x[2];
  
  printf("Metodo dos Retangulos (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), h = (%lg)\n", a, b, npontos, h);
  
  double t_inicial = timestamp();
  
  ////////////////////
  for(double xi = a; xi <= b; xi+=h){
    for(double yj = a; yj <= b; yj+=h){
      x[0] = xi;
      x[1] = yj;
      soma += f(x,2);
    }
  }
  resultado = h*h*soma;
  //////////////////////
  
  double t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}


int main(int argc, char **argv) {
  if (argc < 5) {
    printf("Utilização: %s inicial final n_amostras n_variaveis\n", argv[0]);
    return 1;
  }

  // INICIAR VALOR DA SEMENTE
  time_t t;
  //SRAND((unsigned) time(&t));
  SRAND(1234);

  // CHAMAR FUNÇÕES DE INTEGRAÇÃO E EXIBIR RESULTADOS
  double volumeStyblinski = styblinskiTang(atof(argv[1]),atof(argv[2]),atoi(argv[3]));
  printf("%lf\n", volumeStyblinski);
  double volumeRetangulos = retangulos_xy(atof(argv[1]),atof(argv[2]),atoi(argv[4]));
  printf("%lf\n", volumeRetangulos);
  
  return 0;
}

