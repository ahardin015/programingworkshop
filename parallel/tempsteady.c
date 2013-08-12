/* Sequential Code for Steady-State Heat Problem*/

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EPSILON 0.01
#define BLOCK_LOW(id,p,n)   ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)  (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)  (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)


int main (int argc, char *argv[])
{
  int N = 10;
  int i,j,k;
  double mean = 0.0;  //Average boundary value
  double **u;  //Old values
  double *ustorage;  //Old values
  double globaldiff = 0.0;
  double b;
  int *count;
  int *disp;
  double **w;       //newa Value
  double *wstorage;       //newa Value
  double *nstorage;
  double **n;
  MPI_Status status, status2,status3,status4;
  MPI_Request request,request2,request3,request4;
  double *lowvalue;
  double *highvalue;
  
  int id;
  int p;
  int local_rows;

  double elapsed_time;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  MPI_Barrier(MPI_COMM_WORLD);
  elapsed_time = -MPI_Wtime();

 /* set boundary values and compute mean boundary value */
  local_rows = BLOCK_SIZE(id,p,N);
   double diff = 0.0;
  wstorage = (double *) malloc(local_rows*N*sizeof(double));
  w = (double **)malloc(local_rows*sizeof(double));
  for(i = 0; i < local_rows; i++)
   w[i] = &(wstorage[i*N]);
 
  nstorage = (double *) malloc(local_rows*N*sizeof(double));
  n = (double **)malloc(local_rows*sizeof(double));
    for(i = 0;i<local_rows;i++)
  	n[i]= &(nstorage[i*N]);
	if(n == NULL)
	{
	printf("NULL \n");
	return 0;
	}
	
  lowvalue = (double *)malloc(N*sizeof(double));
  highvalue = (double *)malloc(N*sizeof(double));

  if(id==p-1 ){
  ustorage = (double *)malloc(N*N*sizeof(double));
  u = (double **)malloc(N * sizeof(double *));
  for(i=0;i<N; i++){
  u[i] = &(ustorage[i*N]);
  }

  /*Initialize outervalues*/
      for(i = 0; i<N;i++){
      u[i][0] = u[i][N-1] = u[0][i] = 100.0;
      u[N-1][i] = 0.0;
      mean += u[i][0] + u[i][N-1] + u[0][i] + u[N-1][i];
    }
    mean /= (4.0 * N);
    /* Initialize interior values */

    for (i = 1; i < N-1; i++) {
      for (j = 1; j < N-1; j++){
        u[i][j] = mean;
       }}
// /*Get arrays ready for transfer items*/
 
   count = (int *) malloc(p*sizeof(int));
   disp = (int *) malloc(p*sizeof(int));
   disp[0]=0;
   count[0] = (BLOCK_SIZE(0,p,N));
     for (i=1; i <p; i++){
       disp[i] = (disp[i-1] + count[i-1]);
      count[i] = (BLOCK_SIZE(i,p,N));
}
 for(i=0;i<p;i++){
  disp[i] *= N;
  count[i] *= N;
 } }

for(i=0;i<p; i++)
   
    MPI_Scatterv (ustorage,(count),(disp), MPI_DOUBLE,
                  wstorage,local_rows*N,MPI_DOUBLE, 
                   p-1, MPI_COMM_WORLD );



/* Start of big loop*/
 for(diff=0.1;diff>EPSILON;){
//  for(k =0;k<2;k++){
  globaldiff = 0.0;
  diff = 0.0;

    if(id == 0){
    MPI_Isend(w[(local_rows-1)],N,MPI_DOUBLE,id+1,123,MPI_COMM_WORLD,&request);
    MPI_Irecv(highvalue,N,MPI_DOUBLE,id+1,222,MPI_COMM_WORLD,&request2);
    MPI_Wait(&request, &status);
    MPI_Wait(&request2, &status2);
    }
    if(id == p-1){
    MPI_Isend(w[0],N,MPI_DOUBLE,id-1,222,MPI_COMM_WORLD,&request);
    MPI_Irecv(lowvalue,N,MPI_DOUBLE,id-1,123,MPI_COMM_WORLD,&request2);
    MPI_Wait(&request, &status);
    MPI_Wait(&request2, &status2);
    }
    if(id > 0 && id < (p-1)){
      MPI_Isend(w[0],N,MPI_DOUBLE,id-1,222,MPI_COMM_WORLD,&request);
      MPI_Irecv(highvalue,N,MPI_DOUBLE,id+1,222,MPI_COMM_WORLD,&request2);
      
      MPI_Isend(w[(local_rows-1)],N,MPI_DOUBLE,id+1,123,MPI_COMM_WORLD,&request3);
      MPI_Irecv(lowvalue,N,MPI_DOUBLE,id-1,123,MPI_COMM_WORLD,&request4);
      
      MPI_Wait(&request, &status);
      MPI_Wait(&request2, &status2);      
      MPI_Wait(&request3, &status3);
      MPI_Wait(&request4, &status4);
      }
      
/*Calcualte values*/
  
 if(id == 0){
 	for(i=1;i<=(local_rows-1);i++){
	for(j=1;j<N-1;j++){
	if(local_rows == (i+1)){
	  n[i][j] = 0.25*(highvalue[j]+w[i][j-1]+w[i][j+1]+w[i-1][j]);
	   }
	  else{
	     n[i][j] = 0.25*(w[i-1][j]+w[i][j-1]+w[i][j+1]+w[i+1][j]);
		}
	if(fabs(w[i][j]-n[i][j]) > diff) 
  		diff = fabs(w[i][j]-n[i][j]);
  }}}
	   /* Processors in middle that need both low and high ghost values */
	   if(id < (p-1) && id > 0){
     for(i=0;i<=(local_rows-1);i++){
     for(j=1;j<N-1;j++){
      if(i == 0){
        n[i][j] = ((lowvalue[j]+w[i][j-1]+w[i][j+1]+w[i+1][j])*.25);
   }else if ( local_rows == (i+1)){
        n[i][j] = 0.25*(highvalue[j]+w[i][j-1]+w[i][j+1]+w[i-1][j]);
   }
   else{
    n[i][j] = 0.25*(w[i][j+1]+w[i][j-1]+w[i-1][j]+w[i+1][j]);          
   }
   if(fabs(w[i][j]-n[i][j]) > diff) 
  		diff = fabs(w[i][j]-n[i][j]);
	}}}
 	
   if(id == (p-1)){
        for(i=0;i<=(local_rows-2);i++){
        for(j=1;j<N-1;j++){
          if(i == 0){
         n[i][j] = 0.25*(lowvalue[j]+w[i][j-1]+w[i][j+1]+w[i+1][j]);
  	  
  	} else{
		n[i][j] = 0.25*(w[i][j+1]+w[i][j-1]+w[i-1][j]+w[i+1][j]);
	   }
	   if (fabs(w[i][j] - n[i][j]) > diff)
              	  diff = fabs(w[i][j] - n[i][j]);
	  }}}
 	  
	

	if(id == 0){
	for(i = 1;i<=(local_rows-1); i++){
	for(j=1;j<N-1;j++){
	w[i][j] = n[i][j];	
	}}}
	if(id < (p-1) && id>0){
	  for(i=0;i<=(local_rows-1);i++){
	  for(j=1;j<N-1;j++){
	  w[i][j] = n[i][j];
	}}}
	if(id == (p-1)){
	  for(i=0;i<=(local_rows-2);i++){
	  for(j=1;j<N-1;j++){
	  w[i][j] = n[i][j];
	  }}}
	  MPI_Allreduce(&diff,&globaldiff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
 	  diff = globaldiff;
	 
	  
}	  
// 	  
// 	  
// 	  
//   
// 
// 	/*Gatherv w because the iterations are finished*/
  for(i=0;i<p; i++)
   MPI_Gatherv(wstorage,N*local_rows,MPI_DOUBLE,ustorage,count,disp,MPI_DOUBLE,p-1,MPI_COMM_WORLD);
 
elapsed_time +=MPI_Wtime();
 
 if(id == p-1){
  /*Print solution*/
 for(i = 0;i<N;i++){
 for(j=0;j<N;j++)
 printf("%6.2f ", u[i][j]);
 putchar('\n');
 } 
 printf("Time %10.6f\n",elapsed_time);
 }
 MPI_Finalize();
   return 0;
  }








