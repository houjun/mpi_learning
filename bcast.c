#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
void main(int argc, char *argv[]) {
   int rank;
   double param=-1.0;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   if(rank==5) 
	param=23.0;
   if(rank < 3)
	MPI_Bcast(&param,1,MPI_DOUBLE,5,MPI_COMM_WORLD);
   printf("Proc:%d after broadcast parameter is %f \n",rank,param);
   MPI_Finalize();
}
