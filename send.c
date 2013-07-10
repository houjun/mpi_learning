#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


void main(int argc, char *argv[]){

	int my_rank, proc_size;
	char buf[20];
	char test_100[20] = "HELLO WORLD 100";
	char test_200[20] = "HELLO WORLD 200";
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_size);


	if(my_rank == 0){

		MPI_Status status;
		MPI_Recv(buf, 20, MPI_CHAR, 1, 200, MPI_COMM_WORLD, &status);
                printf("%d:%s\n",my_rank,buf);

                MPI_Recv(buf, 20, MPI_CHAR, 1, 100, MPI_COMM_WORLD, &status);

		printf("%d:%s\n",my_rank,buf);
	}
	else{
		
		int sleep = 1;
		while(sleep);	

		MPI_Send(test_100, 20, MPI_CHAR, 0, 100, MPI_COMM_WORLD);
                MPI_Send(test_200, 20, MPI_CHAR, 0, 200, MPI_COMM_WORLD);



	}

	MPI_Finalize();


}
