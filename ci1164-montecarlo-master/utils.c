#include "utils.h"

// Retorna tempo em segundos
double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec + tp.tv_usec*1.0e-6));
}

