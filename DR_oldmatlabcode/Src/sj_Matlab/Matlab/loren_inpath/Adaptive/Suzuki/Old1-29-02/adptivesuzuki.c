/*
    adaptivesuzuki.c
*/

#include "mex.h"
#include "mexheader.h"

#define NTHETAP 5 /* store 5 points at each timestep: csegt tval1 tval2 tval3 tval4 */

/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{
    double        *spiketimes;
    double        *timearray;
    double        *epst;
    double        *cpt;
    double        *spacing;
    double        *starttimes;
    double        *thetacurrt;
    double        timestep;

    int           nspikes;
    int           ntimesteps;
    int           i;
    int           startindex;

    /* Check numbers of arguments */
    if (nrhs != 5) {
           mexErrMsgTxt("Usage: [thetainit thetahat KSStat] = adaptivesuzuki(spiketimes, timearray, epst, cpt, starttimes)" );
    }
    if (nlhs < 3) {
           mexErrMsgTxt("Usage: [thetainit thetahat KSStat] = adaptivesuzuki(spiketimes, timearray, epst, cpt, starttimes)" );
    }

    /* get the list of spiketimes */
    spiketimes = mxGetPr(prhs[0]);
    nspikes = mxGetM(prhs[0]) * mxGetN(prhs[0]);

    /* get the list of times */
    timearray = mxGetPr(prhs[1]);
    timestep = timearray[1] - timearray[0];
    ntimesteps = mxGetN(prhs[1]) * mxGetM(prhs[1]);

    /* get the epsilons */
    epst = mxGetPr(prhs[2]);

    /* get the coordinates of the control points */
    cpt = mxGetPr(prhs[3]);
    ncpt = mxGetN(prhs[3]) * mxGetM(prhs[3]);
    spacing = (double *) mxCalloc(ncpt, sizeof(double));
    for ( csegt = 0 ; csegt < ncpt ; csegt++ ) {
	    spacing[csegt] = (cpt[csegt+1]-cpt[csegt]);
    }


    /* get the times for the first pass */
    starttimes = mxGetPr(prhs[4]);
    ntrials = mxGetN(prhs[4]) * mxGetM(prhs[4]);

    /* create the output datastructures */
    plhs[0] = mxCreateDoubleMatrix(ncpt, 1,  mxREAL);
    thetainit = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(NTHETAP, ntimesteps,  mxREAL);
    thetahat = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, nspikes,  mxREAL);
    KSStat = mxGetPr(plhs[2]);

    /* create the thetacurr array to hold the current value of thetahat t */
    thetacurrt = (double *) mxCalloc(ncpt, sizeof(double));

    /* set the initial thetahat variables */
    for (i = 0; i < ncpt; i++) {
        thetainit[i] = thetacurrt[i] = 1;
    }

    /* get the index of the first starttime */
    t = 0;
    while (timearray[t] > starttimes[0]) {
        t++;
    }
    startindex = t;

    /* find the first spike */
    nextspike = 0;
    while ((nextspike < nspikes) && (spiketimes[nextspike] <= starttimes[0])) {
        nextspike++;
    }
    startspike = nextspike;

    LambdaTInt = 0;
    lastspiketime = -1;

    for (trial = 0; trial < ntrials; trial++) {
        for (t = startindex; t < ntimesteps; t++) {
            /* find out how many spikes were present in the current timestep */
            dNt = 0;
            while ((nextspike < nspikes) && (spiketimes[nextspike] <= currenttime)){
                nextspike++;
                dNt++;
            }
	    

	}
    }

}


