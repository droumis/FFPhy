/* this function reads in the times from the .postimestamp file into a matlab
 * array and return the timstamps. */
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>

typedef unsigned int u32; 



void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{
    double      *timestamparray;
    double      *tptr;

    int 	i;
    u32 	timestamps[1000000];
    u32 	*tsptr;
    int 	ntimestamps;
    char	filename[200];
    char	tmpstring[200];


    gzFile 	timestampfile;

   /* Check numbers of arguments */
    if (nrhs != 1) {
	    mexErrMsgTxt("Usage: frametimes = getposframetimes(filename)" );
    }
    if (nlhs < 1) {
	    mexErrMsgTxt("Usage: frametimes = getposframetimes(filename)" );
    }


    /* get the name of the timestamp file */
    mxGetString(prhs[0], filename, mxGetNumberOfElements(prhs[0])+1);

    if ((timestampfile = gzopen(filename, "r")) == NULL) {
	mexPrintf("Error opening %s for reading\n", filename);
	return;
    }

    /* skip past the header */
    do {
        gzgets(timestampfile, tmpstring, 200);
    } while ((strncmp(tmpstring, "%%ENDHEADER", 10) != 0) &&
             (strncmp(tmpstring, "%%ENDCONFIG", 10) != 0));

    ntimestamps = 0;
    /* read in the timestamps */
    mexPrintf("reading\n");
    while (gzread(timestampfile, timestamps + ntimestamps++, sizeof(u32)) 
	    == sizeof(u32));
    mexPrintf("read %d timestamps\n", ntimestamps);

    /* create the mxArray */
    plhs[0] = mxCreateDoubleMatrix(ntimestamps, 1, mxREAL);
    timestamparray = mxGetPr(plhs[0]);
    
    tptr = timestamparray;
    tsptr = timestamps;
    /* copy the data into the timestamp array */
    for (i = 0; i < ntimestamps; i++, tptr++, tsptr++) {
	*tptr = (double) *tsptr;
    }
    return;
}





