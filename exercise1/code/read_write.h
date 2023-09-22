#ifndef READ_WRITE_H
#define READ_WRITE_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h> 

// Functions to work with pgm files
void write_pgm_image( void *image, int maxval, long xsize, long ysize, const char *image_name);
void read_pgm_image( void **image, int *maxval, long *xsize, long *ysize, const char *image_name);

#endif

