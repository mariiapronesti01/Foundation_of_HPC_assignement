#include "read_write.h"
#include <mpi.h>
#include <omp.h>

//------------------------STATIC-SERIAL------------------------------------------------------------------------------------------------------------------------------------------------------

//this subroutine is used to update the values of the cells according to the serial static evolution.

void update_static_serial(unsigned char *prevMatrix, unsigned char *newMatrix, long ysize){
	long size=ysize*(ysize+1);
	unsigned char max=255;
	unsigned char min=0;

	#pragma omp parallel for
	for(long i=ysize; i<size; i++){
		int count=0;

		//Get the values of row and column from i
		long r = (i/ysize);
		long c = (i%ysize);

		//Calculate the values of the row above/below and the left/right column while considering the boundary conditions
		long row_above =  r - 1;
		long row_below =  r + 1;
		long col_left = (c == 0) ? ysize-1 : c - 1;
		long col_right = (c == (ysize-1)) ? 0 : c + 1;

		//Get the number of alive neighbours
		count=prevMatrix[row_above*ysize+col_left]+prevMatrix[row_above*ysize+c]+prevMatrix[row_above*ysize+col_right]+prevMatrix[r*ysize+col_left]+prevMatrix[r*ysize+col_right]+prevMatrix[row_below*ysize+col_left]+prevMatrix[row_below*ysize+c]+prevMatrix[row_below*ysize+col_right];
		if((count==1530)||(count==1275))
			newMatrix[i-ysize]=min;
		else
			newMatrix[i-ysize]=max;}}
		
			
//this function is used to perform the static evolution for as many steps (t) as required and to write a snapshot of the evolved matrix every s steps.

void run_static(char * filename, int t, int s){
    
	long xsize = 0;
  	long ysize = 0;
  	int maxval = 0;    
  	unsigned char *prevMatrix;
	unsigned char *newMatrix;

	
	double start_time = MPI_Wtime();
 	
	//Initialize the new matrix by reading the pgm file where it's stored
  	read_pgm_image( (void**)&newMatrix, &maxval, &xsize, &ysize, filename);
  
  	long size=xsize*ysize;
	
	//Allocate the space needed to store the previous matrix
  	prevMatrix=(unsigned char*)malloc( (xsize+2)*ysize*sizeof(unsigned char) );
	

	printf("Performing the static evolution for %d times\n",t);
	//Perform the static evolution
 	char fout[26];
 	for(int i=1;i<=t;i++){
		#pragma omp parallel for
                for(long j=0; j< size; j++){
                        prevMatrix[j+xsize]=newMatrix[j];
                }
	
                #pragma omp parallel for
                for(long j=0; j< xsize; j++){
                        prevMatrix[j]=newMatrix[j+(xsize*(xsize-1))];
                        prevMatrix[((xsize+1)*xsize)+j]=newMatrix[j];
                }

   		update_static_serial(prevMatrix, newMatrix, xsize);

		if((s!=0)&&(i%s==0)){
			//Write the output image 
			sprintf(fout, "static_snapshot_%05d.pgm", i);
			write_pgm_image((void*)newMatrix, maxval, xsize, ysize, fout);
		}
	}

	if(s==0){
		//Write the output image 
		sprintf(fout, "static_snapshot_%05d.pgm", t);
		write_pgm_image((void*)newMatrix, maxval, xsize, ysize, fout);		
    	}	

	free(prevMatrix);
	free(newMatrix);
	
	double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        printf("%ld,%d,%d,%f\n",ysize,1,omp_get_max_threads(),elapsed_time);
   
}



//-----------------------------PARALLEL2----------------------------------------------------------------------------

//this subroutine is used to update the values of the cells according to the parallel static evolution.

void update_static_parallel(unsigned char *prevMatrix, unsigned char *newMatrix, long local_rows, long ysize){
	long size=(local_rows*ysize)-ysize;    //this is calculated here as to not have to recalculate for every iteration of the for loop
	unsigned char max=255;
	unsigned char min=0;

	#pragma omp parallel for
	for(long i=ysize; i<size; i++){
		int count=0;

		//Get the values of row and column from i
		long r = (i/ysize);
		long c = (i%ysize);

		//Calculate the values of the row above/below and the left/right column while considering the boundary conditions
		long row_above = r - 1;
		long row_below = r + 1;
		long col_left = (c == 0) ? ysize-1 : c - 1;
		long col_right = (c == (ysize-1)) ? 0 : c + 1;

		//Get the number of alive neighbours
		count=prevMatrix[row_above*ysize+col_left]+prevMatrix[row_above*ysize+c]+prevMatrix[row_above*ysize+col_right]+prevMatrix[r*ysize+col_left]+prevMatrix[r*ysize+col_right]+prevMatrix[row_below*ysize+col_left]+prevMatrix[row_below*ysize+c]+prevMatrix[row_below*ysize+col_right];
		if((count==1530)||(count==1275))
			newMatrix[i-ysize]=min;
		else
			newMatrix[i-ysize]=max;
	}
}
			

//this function is used to perform the static evolution in parallel for as many steps (t) as required and to write a snapshot of the evolved matrix every s steps.

void run_static_parallel(char * filename, int t, int s, int size, int rank, MPI_Status status, MPI_Request r){
    
  	long xsize = 0;
	long ysize = 0;
  	int maxval = 0;    
  	unsigned char *prevMLocal;
  	unsigned char *newMLocal;
  	unsigned char *completeMatrix;  
	double start_time, end_time;  
 
  	if(rank==0){
		start_time = MPI_Wtime();
		//Initialize the matrix by reading the pgm file where it's stored
		read_pgm_image( (void**)&completeMatrix, &maxval, &xsize, &ysize, filename);
  	}

  	MPI_Barrier(MPI_COMM_WORLD);

  	MPI_Bcast(&ysize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  	MPI_Bcast(&xsize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        
  	long rows_per_process = xsize / size;
  	long extra_rows = xsize % size;
  	long local_rows = rows_per_process + (rank < extra_rows ? 3 : 2);
	
	//Allocate the space for the 2 matrices
  	prevMLocal= (unsigned char*)malloc( local_rows*ysize*sizeof(unsigned char) );
  	newMLocal= (unsigned char*)malloc( (local_rows-2)*ysize*sizeof(unsigned char) ); 

	
    	// Calculate sendcounts and displacements
    	int sendcounts[size];
    	int displacements[size];
    	int total_rows = 0;

	#pragma omp parallel for
    	for (long i = 0; i < size; i++) {
        	int rows_to_send=(i<extra_rows)?(rows_per_process+1):(rows_per_process);
		sendcounts[i] = rows_to_send * ysize;
        	displacements[i] = total_rows * ysize;
        	total_rows += rows_to_send;
    	}

    	//Scatter data: divide and send the initial matrix to all local matrices
    	MPI_Scatterv(completeMatrix, sendcounts, displacements, MPI_UNSIGNED_CHAR, newMLocal, (local_rows-2) * ysize, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	
  	int rank_above= (rank==0)? size-1 : rank-1;
  	int rank_below= (rank==size-1)? 0: rank+1;

  	long Msize=(local_rows*ysize)-ysize;
  
  	char fout[26];

  	for(int i=1;i<=t;i++){
		//Initialize the previous local matrix
		//First copy the middle rows form newMLocal
	
		int tag_odd=2*i;
		int tag_even=2*i+1;

		#pragma omp parallel for
    		for(long i=ysize; i<Msize; i++){
    			prevMLocal[i]=newMLocal[i-ysize];
    		}

    		//Then exchange first and last rows with neighboring processes
    		MPI_Isend(&prevMLocal[ysize], ysize, MPI_UNSIGNED_CHAR, rank_above, tag_odd, MPI_COMM_WORLD,&r);
    		MPI_Recv(&prevMLocal[ysize*(local_rows-1)], ysize, MPI_UNSIGNED_CHAR, rank_below, tag_odd, MPI_COMM_WORLD, &status);
    
   		MPI_Isend(&prevMLocal[ysize*(local_rows-2)], ysize, MPI_UNSIGNED_CHAR, rank_below, tag_even, MPI_COMM_WORLD,&r);
  		MPI_Recv(prevMLocal, ysize, MPI_UNSIGNED_CHAR, rank_above, tag_even, MPI_COMM_WORLD, &status);			

		//Update the values of the cells
    		update_static_parallel(prevMLocal, newMLocal, local_rows, ysize);

    		if((s!=0)&&(i%s==0)){ 
			//Get the values of the final matrix
    			MPI_Gatherv(newMLocal, (local_rows-2) * ysize, MPI_UNSIGNED_CHAR, completeMatrix, sendcounts, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
   			MPI_Barrier(MPI_COMM_WORLD);

			if(rank==0){
				sprintf(fout, "static_snapshot_%05d.pgm", i);
				write_pgm_image((void*)completeMatrix, maxval, xsize, ysize, fout);
			}
		}
	}
	if(s==0){ 
		//Get the values of the final matrix
    		MPI_Gatherv(newMLocal, (local_rows-2) * ysize, MPI_UNSIGNED_CHAR, completeMatrix, sendcounts, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
   		MPI_Barrier(MPI_COMM_WORLD);

		if(rank==0){
			sprintf(fout, "static_snapshot_%05d.pgm", t);
			write_pgm_image((void*)completeMatrix, maxval, xsize, ysize, fout);
		}	
    	}

	//Free the memory
	free(prevMLocal);
	free(newMLocal);
	if(rank==0){
		free(completeMatrix);
		end_time = MPI_Wtime();
        	double elapsed_time = end_time - start_time;
        	printf("PARALLEL %ld,%d,%d,%f\n",ysize,size,omp_get_max_threads(),elapsed_time);
	}
}




