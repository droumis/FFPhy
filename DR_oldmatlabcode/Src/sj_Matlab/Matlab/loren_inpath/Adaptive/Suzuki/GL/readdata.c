#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include "adapt.h"


/*******************************************************************
  
	This function reads data about the rat's estimated path
	 filename - constant array containing path and filename
	 *x       - pointer to x-vector
	 samples  - how many samples to read in
	 skip	  - entries to skip per sample stored
    Return 0 for success, 1 for failure

*********************************************************************/
int readdata(char *filename, void *data, int nsamples, int size)
	/* reads in all of the binary data from the specified file */
{
	FILE 	*infile;
	int 	n;

    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", filename);
        exit(1);
    }
	/* read in the data */
	n = fread(data, size, nsamples, infile);
	fclose(infile);
	return n;
}	

int getnelements(char *filename, int esize)
	/* returns the number of elements of the specificied size in the binary file */
{
	int  nelem = 0;
	struct stat buf;

    if (stat(filename, &buf) != 0) {
        fprintf(stderr, "Error stating %s \n", filename);
        exit(1);
    }

	/* the size is the size of the file divided by the size of the element (esize */
	if (((nelem = buf.st_size / esize) * esize) != buf.st_size) {
        fprintf(stderr, "Warning: %d does not evenly divide size of file %s \n", esize, 
						filename);
	}
	return nelem;
}
