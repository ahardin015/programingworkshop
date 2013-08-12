/* Sequential Code for Steady-State Heat Problem*/

#include <math.h>
#include <time.h>
#define N      500
/*#define N      1000*/
#define EPSILON 0.01

int main (int argc, char *argv[])
{
  double diff;  //Change in value

  int i,j,k;
  double mean;  //Average boundary value
  double u[N][N];       //Old values
  double w[N][N];       //New Value
  double start,stop;
start = clock();
  /* set boundary values and compute mean boundary value */
  mean = 0.0;
  for(i = 0; i<N;i++){
    u[i][0] = u[i][N-1] = u[0][i] = 100.0;
    u[N-1][i] = 0.0;
    mean += u[i][0] + u[i][N-1] + u[0][i] + u[N-1][i];
  }
  mean /= (4.0 * N);

  /* Initialize interior values */

  for (i = 1; i < N-1; i++)
    for (j = 1; j < N-1; j++)
      u[i][j] = mean;

  
  /* compute steady-state solution */
  for (;;) {
    diff = 0.0;
    for (i = 1; i < N-1; i++)
      for (j = 1; j < N-1; j++) {
        w[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1])/4.0;
        if (fabs(w[i][j] - u[i][j]) > diff)
          diff = fabs(w[i][j] - u[i][j]);
        }
      if (diff <= EPSILON) break;
      for (i = 1; i < N-1; i++)
        for (j = 1; j < N-1; j++) u[i][j] = w[i][j];
}
stop=clock();
   /* Print solution */
  for (i = 0; i < N; i++) {
   for (j = 0; j < N;j++)
     printf ("%6.2f ", u[i][j]);
    putchar ('\n');
  }
  printf("Time %f",(stop-start)/1000000);
  return 0;
  }








