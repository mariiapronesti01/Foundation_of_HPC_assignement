# Hybrid Parallel Version of Conway's Game of Life

This folder contains the source code and results for a hybrid parallel version of Conway's Game of Life, which utilizes both MPI and OpenMP for parallel processing. For more detailed information about the implementation, please consult the 'Report.pdf' file.

## Code and Compilation

Inside the 'code' folder, you will find source code and header files, along with a makefile for building the program. To compile the program, ensure that MPI is installed and relevant modules are loaded, and then run the following command:

```bash
make

This command will generate the executable 'my_program.x'.

## Running the Program

To execute the program, use 'mpirun' and specify the desired command-line arguments:

-f: Name of the file for reading/writing game data.
-k: Size of the matrix for initialization (default: 100).
-s: Specify how often to save snapshots of the game world (default: 0, saving only the final output).
-i: Initialize the game world.
-r: Run the game simulation.
-n: Number of simulation steps.
-e: Evolution policy (1 for static, 0 for ordered).


## Scalability Measurements

Detailed measurements related to three types of required scalability studies — OpenMP scalability, strong MPI, and weak MPI scalability — are available for both **EPYC** and **THIN** nodes of the *Orfeo* cluster. You can find these measurements in the corresponding folder, supplied with graphs. For further insights and analysis, please refer to the 'Report.pdf' file.
