#include <stdio.h>
#include <stdlib.h>

main () {
	FILE 	*infile;
	int 	n;
	float	*data;
	int 	nsamples = 26400066;

    infile = fopen("thetahate.dat", "r");
	data = (float *) malloc(sizeof(float) * nsamples);
	/* read in the data */
	n = fread(data, sizeof(float), nsamples, infile);
	fprintf(stderr, "n = %d\n", n);
	fclose(infile);
}	
