#include "read_write.h"
#include <mpi.h>
#include <omp.h>

//------------------------ORDERED-SERIAL--------------------------------------------------------------------------------------------------------------------------


//this subroutine is used to update the values of the cells according to the ordered evolution.

void update_ordered_serial(unsigned char *M, long xsize){
	long size=xsize*xsize;
	unsigned char max=255;
	unsigned char min=0;

	for(long i=0; i<size; i++){
		int count=0;

		//Get the values of row and column from i
		long r = (i/xsize);
		long c = (i%xsize);

		//Calculate the values of the row above/below and the left/right column while considering the boundary conditions
		long row_above = (r == 0) ? xsize-1 : r - 1;
		long row_below = (r == (xsize-1)) ? 0 : r + 1;
		long col_left = (c == 0) ? xsize-1 : c - 1;
		long col_right = (c == (xsize-1)) ? 0 : c + 1;

		//Get the number of alive neighbours
		count= M[row_above*xsize+col_left]+M[row_above*xsize+c]+M[row_above*xsize+col_right]+M[r*xsize+col_left]+M[r*xsize+col_right]+M[row_below*xsize+col_left]+M[row_below*xsize+c]+M[row_below*xsize+col_right];
		if((count==1530)||(count==1275))
			M[i]=min;
		else
			M[i]=max;
	}
}



//this function is used to perform the ordered evolution for as many steps (t) as required and to write a snapshot of the evolved matrix every s steps.

void run_ordered(char * filename, int t, int s){
    
  	long xsize = 0;
  	long ysize = 0;
  	int maxval = 0;    
  	unsigned char *M;

  	double start_time=MPI_Wtime();
	  
  	read_pgm_image((void**)&M, &maxval, &xsize, &ysize, filename);

  	char fout[27];
  	for(int i=1;i<=t;i++){
    		update_ordered_serial(M, xsize);
    		if((s!=0)&&(i%s==0)){ 
			//Write the output image 
			sprintf(fout, "ordered_snapshot_%05d.pgm", i);
			write_pgm_image( (void*)M, maxval, xsize, ysize, fout);		
    		}
  	}
	
	if(s==0){
		//Write the output image 
		sprintf(fout, "ordered_snapshot_%05d.pgm", t);
		write_pgm_image((void*)M, maxval, xsize, ysize, fout);
	}	

	free(M);

  	double end_time=MPI_Wtime();
  	double elapsed_time=end_time-start_time;
        printf("%ld,%d,%d,%f\n",ysize,1,omp_get_max_threads(),elapsed_time);
}






//------------------------ORDERED-PARALLEL--------------------------------------------------------------------------------------------------------------------------

//this subroutine is used to update the values of the cells according to the parallel ordered evolution.

void update_ordered_parallel(unsigned char *M, long local_rows, long ysize){
	long size=(local_rows*ysize)-ysize;  
	unsigned char max=255;
	unsigned char min=0;

	for(long i=ysize; i<size; i++){
		int count=0;

		//Get the values of row and column from i
		long r = (i/ysize);
		long c = (i%ysize);

		//Calculate the values of the row above/below and the left/right column while considering the boundary conditions
		long row_above =  r - 1;
		long row_below = r + 1;
		long col_left = (c == 0) ? ysize-1 : c - 1;
		long col_right = (c == (ysize-1)) ? 0 : c + 1;

		//Get the number of alive neighbours
		count=M[row_above*ysize+col_left]+M[row_above*ysize+c]+M[row_above*ysize+col_right]+M[r*ysize+col_left]+M[r*ysize+col_right]+M[row_below*ysize+col_left]+M[row_below*ysize+c]+M[row_below*ysize+col_right];
		if((count==1530)||(count==1275))
			M[i]=min;
		else
			M[i]=max;}}
		
		

//this function is used to perform the ordered evolution in parallel for as many steps (t) as required and to write a snapshot of the evolved matrix every s steps.

void run_ordered_parallel(char * filename, int t, int s, int size, int rank, MPI_Status status){
    
	long xsize = 0;
 	long ysize = 0;
 	int maxval = 0;    
 	unsigned char *M;
	unsigned char *MLocal;

 	double start_time, end_time;  


 	if(rank==0){
		start_time = MPI_Wtime();
		//Initialize the matrix by reading the pgm file where it's stored
		read_pgm_image( (void**)&M, &maxval, &xsize, &ysize, filename);
  	}


	MPI_Barrier(MPI_COMM_WORLD);

  	MPI_Bcast(&ysize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  	MPI_Bcast(&xsize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        
  	long rows_per_process = xsize / size;
  	long extra_rows = xsize % size;
  	long local_rows = rows_per_process + (rank < extra_rows ? 3 : 2);
  
	//Allocate the space for the 2 matrices
  	MLocal= (unsigned char*)malloc( local_rows*ysize*sizeof(unsigned char) );


    	// Calculate sendcounts and displacements
    	int sendcounts[size];
    	int displacements[size];
    	int total_rows = 0;
    	for (long i = 0; i < size; i++) {
        	int rows_to_send=(i<extra_rows)?(rows_per_process+1):(rows_per_process);
		sendcounts[i] = rows_to_send * ysize;
        	displacements[i] = total_rows * ysize;
        	total_rows += rows_to_send;
    	}

	

    	//Scatter data: divide and send the initial matrix to all local matrices
    	MPI_Scatterv(M, sendcounts, displacements, MPI_UNSIGNED_CHAR, &MLocal[ysize], (local_rows-2) * ysize, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);



	int rank_above= (rank==0)? size-1 : rank-1;
  	int rank_below= (rank==size-1)? 0: rank+1;

  	long Msize=(local_rows*ysize)-ysize;

	if(rank==size-1)
		MPI_Send(&MLocal[ysize*(local_rows-2)], ysize, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);

 	char fout[27];
  	for(int i=1;i<=t;i++){
		if(rank!=0){
			//Send the first line of the current process to the previous process as is (before the update, as the previous process will need it to perform its own updates)
			MPI_Send(&MLocal[ysize], ysize, MPI_UNSIGNED_CHAR, rank_above, 0, MPI_COMM_WORLD);
		}

		//Receive necessary lines
		MPI_Recv(&MLocal[ysize*(local_rows-1)], ysize, MPI_UNSIGNED_CHAR, rank_below, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(MLocal, ysize, MPI_UNSIGNED_CHAR, rank_above, 1, MPI_COMM_WORLD, &status);

		//Update the values of the cells
    		update_ordered_parallel(MLocal, local_rows, ysize);
		
		if(rank==0){
   			//Send the first line (updated) to the last process
                        MPI_Send(&MLocal[ysize], ysize, MPI_UNSIGNED_CHAR, rank_above, 0, MPI_COMM_WORLD);
                }   

		//Send necessary lines to the next process
                MPI_Send(&MLocal[ysize*(local_rows-2)], ysize, MPI_UNSIGNED_CHAR, rank_below, 1, MPI_COMM_WORLD);
 
    		if((s!=0)&&(i%s==0)){ 
			//Get the values of the final matrix
    			MPI_Gatherv(&MLocal[ysize], (local_rows-2) * ysize, MPI_UNSIGNED_CHAR, M, sendcounts, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
   			MPI_Barrier(MPI_COMM_WORLD);

			if(rank==0){
				sprintf(fout, "ordered_snapshot_%05d.pgm", i);
				write_pgm_image((void*)M, maxval, xsize, ysize, fout);
			}	
    		}
  	}

	if(rank==0)
		MPI_Recv(MLocal, ysize, MPI_UNSIGNED_CHAR, size-1, 1, MPI_COMM_WORLD, &status); //Needed to avoid an error at the end

	if(s==0){
		//Get the values of the final matrix
    		MPI_Gatherv(&MLocal[ysize], (local_rows-2) * ysize, MPI_UNSIGNED_CHAR, M, sendcounts, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
   		MPI_Barrier(MPI_COMM_WORLD);

		if(rank==0){
			sprintf(fout, "ordered_snapshot_%05d.pgm", t);
			write_pgm_image((void*)M, maxval, xsize, ysize, fout);
		}	
    	}
		

	//Free the memory
	free(MLocal);
	if(rank==0){
		free(M);
		end_time = MPI_Wtime();
     	   	double elapsed_time = end_time - start_time;
     	   	printf("%ld,%d,%d,%f\n",ysize,size,omp_get_max_threads(),elapsed_time);
	}
}







