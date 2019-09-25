#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "hippoIO.h"

/* read vector of floats (of length nsamples) from into data, memory will not
   be allocated within this function 
   returns: the number of elements read in
*/
int readdata(char *filename, float *data, int nsamples)
{
	FILE 	*infile;
	int 	n;

    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", filename);
        throw(1);
    }
	/* read in the data */
	n = fread(data, sizeof(float), nsamples, infile);
	fclose(infile);
	return n;
}	

/* returns the number of elements of the specificied size in the binary file */
int getnelements(char *filename, int esize)
{
	int  nelem = 0;
	struct stat buf;

	if (stat(filename, &buf) != 0) {
	    fprintf(stderr, "Error stating %s \n", filename);
        throw(1);
	}

	/* the size is the size of the file divided by the size of the element (esize */
	if (((nelem = buf.st_size / esize) * esize) != buf.st_size) {
	    fprintf(stderr, "Warning: %d does not evenly divide size of file %s \n", esize, 
						filename);
	}
	return nelem;
}

double calcMatrixDet(double A[]) { return A[0]*A[3]-A[1]*A[2]; }

// return inverse of 2x2 matrix A in Ainv. A is vector with A{11,12,21,22}
void calcMatrixInv(double A[], double Ainv[])
{
    //  A={a b; c d}, |A|= ad-bc, Ainv= {d -b; -c a}/|A|
    
    double det= A[0]*A[3]-A[1]*A[2];
    hippo_Assert(det, "singular matrix");
    Ainv[0]= A[3]/det;
    Ainv[1]= -A[1]/det;
    Ainv[2]= -A[2]/det;
    Ainv[3]= A[0]/det;
}
