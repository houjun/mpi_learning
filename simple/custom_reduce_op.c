#include <mpi.h>
#include <stdio.h>
#include <math.h>
void onenorm(float *in, float *inout, int *len, MPI_Datatype *type) {  
    int i;
    for (i=0; i<*len; i++) {
        *inout = fabs(*in) + fabs(*inout);     
        in++;
        inout++;
    }
}

void main(int argc, char* argv[]) {

    int proc_num, my_rank;
    int root=0;
    float sendbuf, recvbuf;
    MPI_Op myop;

    int commute=0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
    MPI_Op_create(onenorm, commute, &myop);  
    
    sendbuf = my_rank*-1;
   
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Reduce(&sendbuf, &recvbuf, 1, MPI_FLOAT, myop, root, MPI_COMM_WORLD);   

    if(my_rank == root)
        printf("The operation yields %f\n", recvbuf);
    
    MPI_Finalize();
}
