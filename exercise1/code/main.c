#include <getopt.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include "read_write.h"
#include "static_evolution.h"
#include "ordered_evolution.h"

#define MAXVAL 255

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1





//-------------------------------------------RANDOM_PLAYGROUND---------------------------------------------------------------------------------------------------------------------------

//this function is used to generate the initial random playground.

void * generate_random( int maxval, long xsize, long ysize){
	unsigned char *cImage;  
	void *ptr;
 
	cImage= (unsigned char*)malloc( xsize*ysize*sizeof(unsigned char) );

	srand(time(NULL));

	long xy= xsize*ysize;
	int  half=  maxval/2;
	unsigned char minval= (unsigned char)0;
	unsigned char _maxval= (unsigned char)maxval;

	#pragma omp parallel for
	for(long i=0; i<xy; i++){
        	int random_number = (rand() % (maxval+1));
        	cImage[i]= (random_number > half)?_maxval:minval;
        }

 ptr= (void*)cImage;
 return ptr;
}


//--------------------PARALLEL VERSION----------------------------------------

//this function is used to generate the initial random playground in parallel.

void * generate_random_parallel( int maxval, long xsize, long ysize, int rank, int size){
	unsigned char *cImageLocal;  
	void *ptr;
	unsigned char *cImage=NULL;
  
 	
	long rows_per_process = xsize / size;
    	long extra_rows = xsize % size;
    	long local_rows = rows_per_process + (rank < extra_rows ? 1 : 0);
	
	cImageLocal= (unsigned char*)malloc( local_rows*ysize*sizeof(unsigned char) );
	
	srand(time(NULL)+rank); //seed with rank-specific value

	long xy= local_rows*ysize;
	int  half=  maxval/2;
	unsigned char minval= (unsigned char)0;
	unsigned char _maxval= (unsigned char)maxval;

	#pragma omp parallel for
	for(long i=0; i<xy; i++){
        	int random_number = (rand() % (maxval+1));
        	cImageLocal[i]= (random_number > half)?_maxval:minval;
        }
        
        if (rank == 0) {
        	cImage = (unsigned char*)malloc(xsize * ysize*sizeof(unsigned char));
	}

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

        
	MPI_Gatherv(cImageLocal, local_rows * ysize, MPI_UNSIGNED_CHAR, cImage, sendcounts, displacements, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
   	MPI_Barrier(MPI_COMM_WORLD);
	free(cImageLocal);
	
 	ptr= (void*)cImage;
 	return ptr;

}


//--------------------GET-ARGUMENTS----------------------------------------------------------------------------------------------------------------------------------------------



char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;



//----------------------MAIN---------------------------------------------------------------------------------------------------------------------------------------------------

int main( int argc, char **argv ) 
{   
    	//long xsize	   = XWIDTH;
   	//long ysize	   = YWIDTH;
    	int maxval     = MAXVAL;

    	MPI_Status status;
	MPI_Request req;

    	//Initialize MPI environment
	int mpi_provided_thread_level;
	MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
	if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) { 
		printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n"); 
		MPI_Finalize(); 
		exit( 1 );
	}


    	int rank, size;

    	// Get the rank of the current process
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    	// Get the total number of processes
    	MPI_Comm_size(MPI_COMM_WORLD, &size);

  

	//this part is needed get the values of the arguments from input
    	int action = 0;
    	char *optstring = "irk:e:f:n:s:";

    	int c;
    	while ((c = getopt(argc, argv, optstring)) != -1) {
    		switch(c) {
      
    		case 'i':
     	 		action = INIT; break;
      
    		case 'r':
      			action = RUN; break;
      
    		case 'k':
      			k = atoi(optarg); break;

    		case 'e':
      			e = atoi(optarg); break;

    		case 'f':
      			fname = (char*)malloc( sizeof(optarg)+1 );
      			sprintf(fname, "%s", optarg );
      			break;

    		case 'n':
      			n = atoi(optarg); break;

    		case 's':
      			s = atoi(optarg); break;

    		default :
      			printf("argument -%c not known\n", c ); break;
    		}
  	}

    
	//based on the values we read we can do different things:

    	//initialize playground
    	if(action==INIT){
		printf("Initializing random playground\n");
		if(fname==NULL){
			printf("No name was passed, default name will be used\n");
			fname=(char*)malloc(sizeof(fname_deflt)+1);
			sprintf(fname, "%s", fname_deflt);
		}
		if(size==1){
			void *ptr = generate_random( maxval, k, k );    
 			printf("The random image has been generated\n");
			write_pgm_image( ptr, maxval, k, k, fname);
        		printf("The random image has been written\n");
        		free(ptr);
        	} else {
        		void *ptr=generate_random_parallel(maxval,k,k,rank,size);
        		MPI_Barrier(MPI_COMM_WORLD);
        		if(rank==0){
				printf("The random image has been generated in parallel\n");
				write_pgm_image( ptr, maxval, k, k, fname);
        			printf("The random image has been written\n");
        		}
        		free(ptr);
        	}
    	}




    	//run the game
    	if(action==RUN){
		if (fname == NULL) {
                	printf("No name was passed, the program will try to read from %s\n\n", fname_deflt);
            	fname = (char*) malloc(sizeof(fname_deflt));
            	sprintf(fname, "%s", fname_deflt);
        	}
        	//run the ordered evolution
		if(e==ORDERED){
			if(size==1){
				run_ordered(fname, n, s);
			} else {
               	        	run_ordered_parallel(fname, n, s, size, rank, status );
			}
		//run the static evolution
		}else if(e==STATIC){
			if(size==1){
				run_static(fname, n, s );
			} else {
				run_static_parallel(fname, n, s, size, rank, status, req);
			}
		} else {
			printf("Invalid input: e can only take value in {0,1}\n"); 
			MPI_Finalize(); 
			exit( 1 );
		}
    	}
    
 
	MPI_Finalize();
    	return 0;
} 



