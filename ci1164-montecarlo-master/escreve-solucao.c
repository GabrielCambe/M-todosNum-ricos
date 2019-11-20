void escreve_solucao_gnuplot(char* arq_saida, double tempo_total_GaussSiedel, unsigned int num_iter, Real_t* residuo_iter, Sist_Lin* sist, Real_t* solucao){
  
    const Real_t hx = PI/((sist->nx)+1);
    const Real_t hy = PI/((sist->ny)+1);
    Real_t x = 0;
    Real_t y = 0;

    // Para o gnuplot
    FILE *arquivo_de_dados = fopen( "arqdat", "w+" );
    if( arquivo_de_dados != NULL ){
      fprintf(arquivo_de_dados, "%s", arq_saida);
      fclose(arquivo_de_dados);
    }
    
    FILE *num_pontos_x = fopen( "nx", "w+" );
    if( num_pontos_x != NULL ){
      fprintf(num_pontos_x, "%u", (sist->nx)+2);
      fclose(num_pontos_x);
    }
    
    FILE *num_pontos_y = fopen( "ny", "w+" );
    if( num_pontos_x != NULL ){
      fprintf(num_pontos_y, "%u", (sist->ny)+2);
      fclose(num_pontos_y);
    }
    
    // escreva os valores de x, y e u(x,y) no arquivo de saida
    for(unsigned int i = 0; i < ((sist->nx)+2); ++i, x += hx){
      y = 0;
      for(unsigned int j = 0; j < ((sist->ny)+2); ++j, y += hy){
	printf("%lf\t%lf\t%lf\n", x, y, solucao[ index(i,j,(sist->ny)+2) ]);
      }
      printf( "\n" );
    }
    
    return;
}
