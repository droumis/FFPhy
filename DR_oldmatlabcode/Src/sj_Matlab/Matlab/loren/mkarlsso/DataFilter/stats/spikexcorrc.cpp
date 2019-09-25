#include "mexheader.h"
//#include "fileiomat.h"
#include "mex.h"
#define MSCONVERT 1000

/*
    xcorrmain.cpp 

    Takes two one dimensional arrays of doubles, a bin, and a maximum time, and
    the start and end times for the two spike trains and returns either the 
    square root normatilzed cross correlations or the standard normalized 
    cross correlations.a structure

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
    const char         *fieldnames[] = {"c1vsc2", "c2vsc1"};
    char             tmpstring[80];
    mxArray            *datastruct;
    mxArray            *mxtmp;
    mxArray            *datatmp1;
    mxArray            *datatmp2;
    double            *data;
    double            *c1;
    double            *c2;
    double            *c1start, *c1tmp, *c2tmp, *c1end, *c2end;
    double            *tmp;
    double            tmpnumspikes;
    double            intstart;
    double            bin;
    double            tmax;
    bool		sqrtnorm;
    int                numspikes1;
    int                numspikes2;
    int             dims[2]= {1,1};
    int             ncontrec;
    int             offset;
    int             currenti;
    int             i, j, k;
    int             numbins;
    int             *neffect;
    

      /* Check numbers of arguments */
      if (nrhs != 5) {
           mexErrMsgTxt("Usage: xcorrstruct = spikexcorr(spikes1, spikes2, bin, tmax, sqrtnorm" );
    }
      if (nlhs < 1) {
           mexErrMsgTxt("Usage: xcorrstruct = spikexcorr(spikes1, spikes2, bin, tmax, sqrtnorm" );
      }

    /* get the first list of times */
    c1 = mxGetPr(prhs[0]);
    numspikes1 = mxGetN(prhs[0]) * mxGetM(prhs[0]);
    /* get the second list of times */
    c2 = mxGetPr(prhs[1]);
    numspikes2 = mxGetN(prhs[1]) * mxGetM(prhs[1]);
    /* get the bin */
    tmp = mxGetPr(prhs[2]);
    bin = *tmp;
    /* get tmax */
    tmp = mxGetPr(prhs[3]) ;
    tmax = *tmp;
    /* get sqrtnorm, convert to int */
    tmp = mxGetPr(prhs[4]) ;
    sqrtnorm = (bool) *tmp;

    /* create the datastructure */
    datastruct = mxCreateStructArray(2, dims, 2, fieldnames);

    /* compute the number of bins.  */
    numbins = ceil(tmax / bin) * 2;

    /* compute the correlation histograms */
    /* the following uninterpretable code does both the 1x2 and the 2x1
     * correlations */

    for (i = 0; i < 2; i++) {
        /* compute both correlations */
        if (i == 0) {
            /* allocate space for the first list of cross correlation values */
            datatmp1 = mxCreateDoubleMatrix(1, numbins, mxREAL);
            data = mxGetPr(datatmp1);
            c1tmp = c1start = c1;
            c1end = c1 + numspikes1;
            c2end = c2 + numspikes2;
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
                data[(int) floor((*c1tmp - intstart) / bin)]++;
                c1tmp++;
            }
            /* move on to the next c2 spike */
            c2tmp++;
        }

	tmp = data;
        for (k = 0; k < numbins; k++, tmp++) {
	    if (sqrtnorm) {
		/* normalize and take the square root of the histogram bins */
		*tmp = sqrt(*tmp / (bin * numspikes2));
	    }
        }
        if (i == 0) {
            /* set the data field */
            mxSetField(datastruct, 0, fieldnames[0], datatmp1);
	    /* if we are not normalizing we can finish here */
	    if (!sqrtnorm) {
		break;
	    }
        }
        else {
            /* set the data field */
            mxSetField(datastruct, 0, fieldnames[1], datatmp2);
        }
    }
    plhs[0] = datastruct;
    return;
}

