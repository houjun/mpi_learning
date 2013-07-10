#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
/*
 * A guess number program, proc 0 has a number for all other procs to guess
 * output the proc id, the number, and the rank of the guesser.
 * Each time the guesser missed, proc 0 will tell whether the guessed number
 * is bigger or smaller.
 */
 
#define MAXNUM 1000000
void main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);

	int my_rank, nproc;
	int i, gotcha = 0;
	int low = 0, high = MAXNUM;
	int count = 0;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int* guess = malloc(nproc*sizeof(int));
	int* gol = malloc(nproc*sizeof(int));
	int number, guessed = -1, goled;


	srand(time(NULL)+(my_rank + 1));

	if(my_rank == 0){
		number = rand() % MAXNUM;
		while(gotcha == 0){

			MPI_Gather(&guessed, 1, MPI_INT, guess, 1, MPI_INT, 0, MPI_COMM_WORLD);

			for(i=1; i<nproc; i++){
				if(number == guess[i]){
					gotcha = 1;
					printf("Guess #%d: %d from proc %d, Gocha!\n",count,guess[i],i );
				}
				else{
					if(number > guess[i])
					gol[i] = 1;
					else
					gol[i] = 0;
			//		printf("Guess #%d: %d from proc %d is %s!\n",count, guess[i], i, gol[i]==1?"smaller":"larger");

				}
				

			}

			MPI_Bcast(&gotcha, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(gotcha == 1)
				break;
			MPI_Scatter(gol, 1, MPI_INT,&goled, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			count++;
		}

	}
	else{
		while(gotcha == 0){
			guessed = rand() % (high - low) + low;
			MPI_Gather(&guessed, 1, MPI_INT, guess, 1, MPI_INT, 0, MPI_COMM_WORLD);

			MPI_Bcast(&gotcha, 1, MPI_INT, 0, MPI_COMM_WORLD);

			if(gotcha == 1)
			break;

			MPI_Scatter(gol, 1, MPI_INT, &goled, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(goled == 1)
			low = guessed;
			else
			high = guessed;

			MPI_Barrier(MPI_COMM_WORLD); 
		}

	}


	MPI_Finalize();
}
