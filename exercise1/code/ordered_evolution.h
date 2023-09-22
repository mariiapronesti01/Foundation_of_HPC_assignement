#ifndef ORDERED_EVOLUTION_H
#define ORDERED_EVOLUTION_H

// Function to update the values of the cells in the ordered evolution.
void update_ordered_serial(unsigned char *M, long xsize);

// Function to perform ordered evolution for the specified number of steps and save snapshots.
void run_ordered(char *filename, int t, int s);



// Function to update the values of the cells in the ordered evolution in parallel.
void update_ordered_parallel(unsigned char *M, long local_rows, long ysize);

// Function to perform ordered evolution in parallel for the specified number of steps and save snapshots.
void run_ordered_parallel(char * filename, int t, int s, int size, int rank, MPI_Status status);

#endif // ORDERED_EVOLUTION_H

