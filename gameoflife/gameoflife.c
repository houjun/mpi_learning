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
 *
 */


/* 
 * For now the program divide the matrix(2-dimensional array of cells) by # of proc evenly.
 * Proc 0 is in charge of distribution, each proc get arr_size_x/comm_size rows of cells and calculate on their own.
 * After each time step, all proc send their data to proc 0 and then start over again.
 * TODO: When arr_size_x cannot be evenly divided by proc #
 * TODO: divide into submatrix instead of horizontally and measure the time
 *
 */
 
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define DEBUG 0
#define PERCENTAGE 60
#define SIMUTIME 10

void print_cell(int **cell, int arr_size_x, int arr_size_y, int print_ghostcell){
    int i, j, p;
    if(print_ghostcell == 1)
    	p = 0;
    else
    	p = -1;
    for (i = 0 - p; i < arr_size_x + p; i++) {
        for (j = 0 - p; j < arr_size_y + p; j++) {
            printf("%d ",cell[i][j]);
        }
        printf("\n");
    }   
    printf("\n");
    fflush(stdout);

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
//	srand(time(NULL));
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
    int i,j,count,livecount,totalcount;
    char intmp[10];
    double start_time, end_time, elapsed_time, all_time;

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
	int *old_loc_cell_1d = malloc((rank_sub_arr_x+2) * rank_sub_arr_y * sizeof(int));
	int **old_loc_cell = malloc((rank_sub_arr_x+2) * sizeof(int*));
	int *new_loc_cell_1d = malloc((rank_sub_arr_x+2) * rank_sub_arr_y * sizeof(int));
	int **new_loc_cell = malloc((rank_sub_arr_x+2) * sizeof(int*));
	int *tmp;

	int **cell;
	int print_ghostcell = 1;

	for(i = 0; i < (rank_sub_arr_x+2); i++){
		old_loc_cell[i] = &old_loc_cell_1d[i * rank_sub_arr_y];
		new_loc_cell[i] = &new_loc_cell_1d[i * rank_sub_arr_y];
	}


    // real rank sub-matrix size
    if(my_rank == 0 || my_rank == comm_size - 1)
    	rank_sub_arr_x++;
    else
    	rank_sub_arr_x += 2;

    if(my_rank == 0){

    			// proc 0 will generate and store the whole matrix
    			int *cell_1d = malloc(arr_size_x * arr_size_y * sizeof(int));
    			cell = malloc(arr_size_x * sizeof(int*));
    			for(i = 0; i < arr_size_x; i++){
    				cell[i] = &cell_1d[i * arr_size_y];
    			}

    			gen_life(cell, arr_size_x, arr_size_y);
    			//print_cell(cell, arr_size_x, arr_size_y, 1);

    			update_ghostcell(cell, arr_size_x, arr_size_y);
    			printf("Origin:\n");
    			print_cell(cell, arr_size_x, arr_size_y, print_ghostcell);
    }

    while(1){

    	start_time = MPI_Wtime();

		if(my_rank == 0){

			// distribute matrix
			// special operation to first: just copy, no need to send.
			for(i = 0; i < rank_sub_arr_x ; i++){
				for(j = 0; j < rank_sub_arr_y; j++)
					old_loc_cell[i][j] = cell[i][j];
			}
			// send to last rank
			MPI_Send(cell[(comm_size - 1) * (rank_sub_arr_x - 1) - 1], rank_sub_arr_x * rank_sub_arr_y
					, MPI_INT, comm_size - 1, 0, MPI_COMM_WORLD);

			// all other ranks
			for(i = 1; i < comm_size - 1; i++){
				MPI_Send(cell[((rank_sub_arr_x - 1) * i) - 1], (rank_sub_arr_x + 1) * rank_sub_arr_y, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

		}


		if(my_rank != 0)
			MPI_Recv(old_loc_cell[0], rank_sub_arr_x * rank_sub_arr_y, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

	#if DEBUG
		printf("From Proc %d:\n", my_rank);
		print_cell(old_loc_cell, rank_sub_arr_x, rank_sub_arr_y, print_ghostcell);
	#endif

		// populate
		livecount = polpulate(old_loc_cell, new_loc_cell, rank_sub_arr_x, rank_sub_arr_y);
		MPI_Reduce(&livecount, &totalcount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	#if DEBUG
		printf("After Popoulate from Proc %d:\n", my_rank);
		print_cell(old_loc_cell, rank_sub_arr_x, rank_sub_arr_y, print_ghostcell);
	#endif


		// collect results
		if(my_rank == 0){
			for(i = 1; i < comm_size; i++)
				MPI_Recv(cell[i * (rank_sub_arr_x - 1)], (rank_sub_arr_x - 1) * rank_sub_arr_y, MPI_INT, i, i, MPI_COMM_WORLD, &status);

			for(i = 1; i < rank_sub_arr_x - 1; i++){
				for(j = 0; j < rank_sub_arr_y; j++)
					cell[i][j] = new_loc_cell[i][j];
			}

			update_ghostcell(cell, arr_size_x, arr_size_y);
	    	printf("Living cells:%d\n", totalcount);
			print_cell(cell, arr_size_x, arr_size_y, print_ghostcell);

		}
		else{
			// only send row 1 ~ rank_sub_arr_x - 1
			MPI_Send(new_loc_cell[1], (rank_sub_arr_x - 2) * rank_sub_arr_y, MPI_INT, 0, my_rank, MPI_COMM_WORLD);
		}

    	end_time = MPI_Wtime();
    	elapsed_time = end_time - start_time;
		MPI_Reduce(&elapsed_time, &all_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    	if(my_rank == 0){

    		printf("Time needed: %f\n", all_time);

    		printf("Enter n to proceed, q to exit\n");
			scanf("%s",intmp);

			if(totalcount<=0){
				printf("All died out!\n");
				strcpy(intmp,"q"); //no need to go on once all died
			}
			MPI_Bcast(intmp, 10, MPI_CHAR, 0, MPI_COMM_WORLD);
			if(strcmp(intmp,"q")==0)
				break;
    	}
    	else{

	        MPI_Bcast(intmp, 10, MPI_CHAR, 0, MPI_COMM_WORLD);
	    	if(strcmp(intmp,"q")==0)
	    		break;
    	}

    }

	MPI_Finalize();

}



//	    if(1){
//	    	int sleep=1;
//	    	while(sleep){;}
//	    }
//	    MPI_Barrier(MPI_COMM_WORLD);





















