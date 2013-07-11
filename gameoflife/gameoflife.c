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
#include <time.h>
#define DEBUG 1
#define PERCENTAGE 60
#define SIMUTIME 10

void print_cell(int **cell, int arr_size_x, int arr_size_y, int print_ghostcell){
    int i, j, p;
    if(print_ghostcell == 1)
    	p = 0;
    else
    	p = -1;
    for (i = 0; i < arr_size_x + p; i++) {
        for (j = 0; j < arr_size_y + p; j++) {
            printf("%d ",cell[i][j]);
        }
        printf("\n");
    }   

}

void update_ghostcell(int **cell, int arr_size_x, int arr_size_y){
    int i,j;

    // four corner
    cell[0][0] 							= cell[arr_size_x - 2][arr_size_y - 2];
    cell[0][arr_size_y - 1] 			= cell[arr_size_x - 2][1];
    cell[arr_size_x - 1][0] 			= cell[1][arr_size_y - 2];
    cell[arr_size_x - 1][arr_size_y - 1] = cell[1][1];

    // four border
    for(i = 1; i < arr_size_x - 1; i++){
        cell[i][0] = cell[i][arr_size_y - 2];
        cell[i][arr_size_y - 1] = cell[i][1];
    }
    for(i = 1; i < arr_size_y - 1; i++){
        cell[0][i] = cell[arr_size_x - 2][i];
        cell[arr_size_x - 1][i] = cell[1][i];
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


void gen_life(int **cell, int arr_size_x, int arr_size_y){
	//r andomly generate life
	srand(time(NULL));
	int i,j,tmp;
	for (i = 0; i < arr_size_x; i++) {
        for (j = 0; j < arr_size_y; j++) {
			if(i == 0 || j == 0 || i == arr_size_x - 1 || j == arr_size_y - 1){
				cell[i][j] = -1; // ghost cell init
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

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    if(my_rank==0){
	    printf("Enter array size:");
	    scanf("%d %d",&arr_size_x, &arr_size_y);
    }

    // include ghost cell size to both
    arr_size_x += 2;
    arr_size_y += 2;

    MPI_Status status;
    // after broadcast each proc knows arr_size_x and arr_size_y
    MPI_Bcast(&arr_size_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&arr_size_y, 1, MPI_INT, 0, MPI_COMM_WORLD);

//    if(my_rank == 1){
//    	int sleep=1;
//    	while(sleep){;}
//    }

    // now we know how to divide the matrix vertically between ranks
    int rank_sub_arr_x = arr_size_x / comm_size;
    int rank_sub_arr_y = arr_size_y;
	//printf("PROC %d:%d %d\n",my_rank, rank_sub_arr_x, rank_sub_arr_y);

    /* Now only consider evenly distributed situation
     * i.e. 20 * 30 matrix size with 4 * 3 processes
     * TODO:possible residue
    int rank_sub_arr_x_res = arr_size_x + arr_size_x % comm_size;
    int rank_sub_arr_y_res = arr_size_y;
     */


    // allocate local storage
	int *loc_cell_1d = malloc((rank_sub_arr_x+2) * rank_sub_arr_y * sizeof(int));
	int **loc_cell = malloc((rank_sub_arr_x+2) * sizeof(int*));

	for(i = 0; i < (rank_sub_arr_x+2); i++){
		loc_cell[i] = &loc_cell_1d[i * rank_sub_arr_y];
	}


    if(my_rank == 0){

    	// proc 0 will generate and store the whole matrix
    	int *cell_1d = malloc(arr_size_x * arr_size_y * sizeof(int));
    	int **cell = malloc(arr_size_x * sizeof(int*));
    	for(i = 0; i < arr_size_x; i++){
    		cell[i] = &cell_1d[i * arr_size_y];
    	}

    	int print_ghostcell = 1;
    	gen_life(cell, arr_size_x, arr_size_y);
    	//print_cell(cell, arr_size_x, arr_size_y, 1);

    	update_ghostcell(cell, arr_size_x, arr_size_y);
    	print_cell(cell, arr_size_x, arr_size_y, print_ghostcell);

    	//printf("PROC 0 Comp Size:%d %d\n",arr_size_x, arr_size_y);

    	// distribute matrix
    	// special operation to first and last rank
    	for(i = 0; i < (rank_sub_arr_x + 1) * rank_sub_arr_y; i++){
    		loc_cell[i] = cell[i];
    	}

    	MPI_COMM_Send(cell[(rank_sub_arr_x * arr_size_x) - 1], (rank_sub_arr_x + 1) * rank_sub_arr_y, MPI_INT, comm_size - 1, 0, MPI_COMM_WORLD);

    	// all other ranks
    	for(i = 1; i < comm_size - 1; i++){
    		MPI_COMM_Send(cell[(rank_sub_arr_x * i) - 1], (rank_sub_arr_x + 2) * rank_sub_arr_y, MPI_INT, i, 0, MPI_COMM_WORLD);
    	}

    }
    else if(my_rank == comm_size - 1){


    	MPI_COMM_Recv(loc_cell, (rank_sub_arr_x + 1) * rank_sub_arr_y, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);



    }
    else{
    	MPI_COMM_Recv(loc_cell, (rank_sub_arr_x + 2) * rank_sub_arr_y, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);


    }

	MPI_Finalize();

}
























