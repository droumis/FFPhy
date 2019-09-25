#include "mexheader.h"
#include "mex.h"
#define MSCONVERT 1000

/*
	xcorrmain.c 

	Takes two one dimensional arrays of doubles, a bin, and a maximum time, and
	the start and end times for the two spike trains and returns the square 
	root normatilzed cross correlations, and the mean and confidence bounds for 
	the xcorr in a structure

*/


/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{
	const char 		*fieldnames[] = {"data1", "data2"};
	char 			tmpstring[80];
	mxArray			*datastruct;
	mxArray			*mxtmp;
	mxArray			*datatmp1;
	mxArray			*datatmp2;
	double			*data;
	double			*c1;
	double			*c2;
	double			*c1start, *c1tmp, *c2tmp, *c1end;
	double			*tmp;
	double			tmpnumspikes;
	double			intstart;
	double			bin;
	double			tmax;
	int				numspikes1;
	int				numspikes2;
	int 			dims[2]= {1,1};
	int 			ncontrec;
    int         	offset;
    int         	currenti;
    int         	i, j, k;
	int 			numbins;
	int 			*numeffectspikes;
	

  	/* Check numbers of arguments */
  	if (nrhs != 4) {
   		mexErrMsgTxt("Usage: xcorrstruct = spikexcorr(spikes1, spikes2, bin, tmax" );
	}
  	if (nlhs < 1) {
   		mexErrMsgTxt("Usage: xcorrstruct = spikexcorr(spikes1, spikes2, bin, tmax" );
  	}

	/* get the first list of times */
	c1 = mxGetPr(prhs[0]);
	numspikes1 = mxGetN(prhs[0]) * mxGetM(prhs[0]);
	/* get the second list of times */
	c2 = mxGetPr(prhs[1]);
	numspikes2 = mxGetN(prhs[1]) * mxGetM(prhs[1]);
	/* get the bin, convert to milliseconds */
	tmp = mxGetPr(prhs[2]);
	bin = *tmp / MSCONVERT;
	/* get tmax, convert to milliseconds */
	tmp = mxGetPr(prhs[3]) ;
	tmax = *tmp / MSCONVERT;

	/* create the datastructure */
	datastruct = mxCreateStructArray(2, dims, 2, fieldnames);

	/* compute the number of bins. note that there is a bin that is centered at
	   t = 0 */
	numbins = ceil((tmax - bin/2)  / bin) * 2 + 1;

	/* allocate space for the number of effective spikes 
		This variable substitues for Nb from Brillinger (1970) because for
		those spikes at the beginning of the spike train of B cannot
		be preceeded by a spike from A, and thus should not be counted 
	numeffectspikes = (int *) mxMalloc(sizeof(int) * numbins); */
	
	/* compute the correlation histograms */
	/* the following uninterpretable code does both the 1x2 and the 2x1
	 * correlations */

	/* for the moment, I am not going to worry about the number of
	effective spikes, and I am just going to assume the number of
	spikes is the number in the interval */
	
	for (i = 0; i < 2; i++) {
		/* compute both correlations */
		if (i == 0) {
			/* allocate space for the first list of cross correlation values */
			datatmp1 = mxCreateDoubleMatrix(1, numbins, mxREAL);
			data = mxGetPr(datatmp1);
			c1tmp = c1start = c1;
			c1end = c1 + numspikes1;
			c2tmp = c2;
		}
		else {
			/* allocate space for the second list of cross correlation values */
			datatmp2 = mxCreateDoubleMatrix(1, numbins, mxREAL); 
			data = mxGetPr(datatmp2);
			/* switch the numbers of spikes */
			c1tmp = c1start = c2;
			c1end = c2 + numspikes2;
			c2tmp = c1;
			tmpnumspikes = numspikes1;
			numspikes1 = numspikes2;
			numspikes2 = tmpnumspikes;
		}
		for (k = 0; k < numbins; k++) {
			data[k] = 0;
		}

		for (k = 0; k < numspikes2; k++) {
			/* fill up the histogram */
			/* start at the first time within tmax of c2 */
			while ((*c1start + tmax < *c2tmp)&&(c1start != c1end)) {
				c1start++;
			}	

			if ((*c1start + tmax < *c2tmp) && (c1start == c1end)) {
				/* we are at the end of the first spike train */
				break;
			}

			/* the beginning of the interval is tmax before c2tmp */
			intstart = *c2tmp - tmax;
			c1tmp = c1start;
			while ((*c1tmp - tmax < *c2tmp) && (c1tmp != c1end)) {
				/* add one to the relevant bin */
				data[round((*c1tmp - intstart) / bin)]++;
				c1tmp++;
			}
			/* move on to the next c2 spike */
			c2tmp++;
		}

		/* normalize and take the square root of the histogram bins */
		for (k = 0; k < numbins; k++) {
			data[k] = sqrt(data[k] / (bin * numspikes2));
		}
		if (i == 0) {
			/* set the data field */
			mxSetField(datastruct, 0, fieldnames[0], datatmp1);
		}
		else {
			/* set the data field */
			mxSetField(datastruct, 0, fieldnames[1], datatmp2);
		}
	}
	plhs[0] = datastruct;
	return;
}

