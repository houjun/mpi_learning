/*
 * "Game of Life" is a simple simulation developed by John Conway. In the "Game of Life" the domain is a 
 * 2-dimensional array of cells, and each cell can have one of two possible states,"alive" or "dead." 
 * The array is usually intialized using random numbers and the system then evolves in time. At each time step, 
 * each cell may or may not change its state, based on the number of adjacent alive cells, including diagonals. 
 *
 * There are three rules:
 * 1. If a cell has three neighbors that are alive, the cell will be alive. If it was already alive, it will remain so, 
 *    and if it was dead, it will become alive.
 * 2. If a cell has two neighbors that are alive, there is no change to the cell. If it was dead, it will remain dead, 
 *    and if it was alive, it will remain alive.
 * 3. In all other cases — the cell will be dead. If it was alive it becomes dead and if it was dead it remains dead.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define DEBUG 1
#define PERCENTAGE 60
#define SIMUTIME 10

void print_cell(int **cell, int arr_size_x, int arr_size_y){
    int i,j;
    for (i = 1; i < arr_size_x - 1; i++) {
        for (j = 1; j < arr_size_y - 1; j++) {
            printf("%d ",cell[i][j]);
        }
        printf("\n");
    }   

}

void update_ghostcell(int **cell, int arr_size_x, int arr_size_y){
    int i,j;
    for(i = 1; i < arr_size_x - 1; i++){
        cell[i][0] = cell[i][arr_size_y-2];
        cell[i][arr_size_y - 1] = cell[i][1];
    }
}

int count_neighbor(int **cell, int arr_size_x, int arr_size_y, int i, int j){
    int count = 0;
    if(i < 1 || i > arr_size_x - 2 || j < 1 || j > arr_size_y - 2)
        return -1;
    else          
        return cell[i-1][j-1] + cell[i-1][j] + cell[i-1][j+1] + cell[i][j-1] 
				+ cell[i][j+1] + cell[i+1][j-1] + cell[i+1][j] + cell[i+1][j+1];
   
}


void gen_life(int **cell, int arr_size_x, int arr_size_y, int rank){
	//randomly generate life
	srand(time(NULL)*rank);
	int i,j,tmp;
	for (i = 0; i < arr_size_x; i++) {
        for (j = 0; j < arr_size_y; j++) {
			if(i == 0 || j == 0 || i == arr_size_x - 1 || j == arr_size_y - 1){
				cell[i][j] = -1; //ghost cell init
				continue;			
			}		
			
            tmp = rand()%100;
            if(tmp >= PERCENTAGE)
                cell[i][j] = 1;
            else
                cell[i][j] = 0;
        }
    }

}

int polpulate(int **oldcell, int **newcell, int arr_size_x, int arr_size_y){
	int i,j,livecount = 0,count;
	for (i = 1; i < arr_size_x - 1; i++) {
        for (j = 1; j < arr_size_y - 1 ; j++) {
            count = count_neighbor(oldcell, arr_size_x, arr_size_y, i, j);
            if(count < 0)
                printf("ERROR!\n");
            if(count == 3)
                newcell[i][j] = 1;
            else if(count == 2)
                newcell[i][j] = oldcell[i][j];
            else
                newcell[i][j] = 0;
                
            if(newcell[i][j] == 1)
                livecount++;
        }
    }   
	return livecount;
}


void main(int argc, char* argv[]) {
    
    int my_rank, comm_size;
	int arr_size_x, arr_size_y;
    int i,j,tmp,count,iterationstep,livecount,pastlivecount;
    int** oldcell;
    int** newcell;
    int** tmpcell;
	int** complementcell;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if(my_rank==0){
	    printf("Enter array size:");
	    scanf("%d %d",&arr_size_x, &arr_size_y);
    }

	MPI_Status status;
    MPI_Bcast(&arr_size_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&arr_size_y, 1, MPI_INT, 0, MPI_COMM_WORLD);

#if 0
	printf("%d: x=%d, y=%d\n",my_rank,arr_size_x, arr_size_y);
#endif

    oldcell = malloc(sizeof(int*) * (arr_size_x + 2));
    newcell = malloc(sizeof(int*) * (arr_size_x + 2));
    complementcell = malloc(sizeof(int*) * ((arr_size_x+1)/2 + 2));
	

	if(my_rank == 0){
	//deal with upper half
		for (i = 0; i < arr_size_x/2 + 2; i++) {
			oldcell[i] = malloc(sizeof(int) * (arr_size_y + 2));
			newcell[i] = malloc(sizeof(int) * (arr_size_y + 2));
		}

		//store other half
		for (i = 0; i < (arr_size_x+1)/2 + 2; i++) {
           		complementcell[i] = malloc(sizeof(int) * (arr_size_y + 2));
	    	}
	
		gen_life(oldcell, arr_size_x/2 + 2, arr_size_y + 2, my_rank);

		iterationstep = 0;

		while(1){
    		MPI_Bcast(&iterationstep, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(iterationstep > SIMUTIME)
				break;
			printf("After %d unit of time:\n",iterationstep);


			update_ghostcell(oldcell, arr_size_x/2 + 2, arr_size_y + 2);

			//receive ghost cell info
			MPI_Recv(oldcell[0], arr_size_y + 2, MPI_INT, 1, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(oldcell[arr_size_x/2 + 1], arr_size_y + 2, MPI_INT, 1, 200, MPI_COMM_WORLD, &status);

			//send to ghost cell		
			MPI_Send(oldcell[1], arr_size_y + 2, MPI_INT, 1, 200, MPI_COMM_WORLD);
			MPI_Send(oldcell[arr_size_x/2], arr_size_y + 2, MPI_INT, 1, 100, MPI_COMM_WORLD);

			print_cell(oldcell, arr_size_x/2 + 2, arr_size_y + 2);

			for(i = 0; i < (arr_size_x+1)/2 + 2; i++)
				MPI_Recv(complementcell[i], arr_size_y + 2, MPI_INT, 1, i, MPI_COMM_WORLD, &status);

			//print other half		
			print_cell(complementcell, (arr_size_x+1)/2 + 2, arr_size_y + 2);

			//MPI_Barrier(MPI_COMM_WORLD);

			polpulate(oldcell, newcell, arr_size_x/2 + 2, arr_size_y + 2);

			tmpcell = oldcell;
			oldcell = newcell;
			newcell = tmpcell;	
			iterationstep ++;
		}

	}
	else if(my_rank == 1){
	//deal with lower half
		for(i = 0; i < (arr_size_x+1)/2 + 2; i++) {
        	oldcell[i] = malloc(sizeof(int) * (arr_size_y + 2));
        	newcell[i] = malloc(sizeof(int) * (arr_size_y + 2));
    	}

		gen_life(oldcell, (arr_size_x+1)/2 + 2, arr_size_y + 2, my_rank);

		while(1){
			MPI_Bcast(&iterationstep, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(iterationstep > SIMUTIME)
				break;
			update_ghostcell(oldcell, (arr_size_x+1)/2 + 2, arr_size_y + 2);

			//send to ghost cell		
			MPI_Send(oldcell[1], arr_size_y + 2, MPI_INT, 0, 200, MPI_COMM_WORLD);
			MPI_Send(oldcell[(arr_size_x+1)/2], arr_size_y + 2, MPI_INT, 0, 100, MPI_COMM_WORLD);

			//receive
			MPI_Recv(oldcell[0], arr_size_y + 2, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(oldcell[(arr_size_x+1)/2 + 1], arr_size_y + 2, MPI_INT, 0, 200, MPI_COMM_WORLD, &status);

			for(i = 0; i < (arr_size_x+1)/2 + 2; i++)
				MPI_Send(oldcell[i], arr_size_y + 2, MPI_INT, 0, i, MPI_COMM_WORLD);

			polpulate(oldcell, newcell, arr_size_x/2 + 2, arr_size_y + 2);		


			tmpcell = oldcell;
			oldcell = newcell;
			newcell = tmpcell;	
		}
	}

		


	MPI_Finalize();

}

