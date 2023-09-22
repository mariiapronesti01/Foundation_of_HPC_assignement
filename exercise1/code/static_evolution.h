#ifndef STATIC_EVOLUTION_H
#define STATIC_EVOLUTION_H

// Function to update the values of the cells in the static evolution.
void update_static_serial(unsigned char *prevMatrix, unsigned char *newMatrix, long xsize);

// Function to perform static evolution for the specified number of steps and save snapshots.
void run_static(char *filename, int t, int s);



// Function to update the values of the cells in the parallel static evolution. 
void update_static_parallel(unsigned char *prevMatrix, unsigned char *newMatrix, long local_rows, long ysize);

// Function to perform static evolution in parallel for the specified number of steps and save snapshots.
void run_static_parallel(char * filename, int t, int s, int size, int rank, MPI_Status status, MPI_Request r);

#endif // STATIC_EVOLUTION_H

