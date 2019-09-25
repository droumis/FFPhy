/*
    adaptivesuzuki.c 
    2 D version for spatial and temporal spline

*/

#include "mex.h"
#include "mexheader.h"

#define NTHETAP 10 /* store 10 points at each timestep: csegx xval1 xval2 xval3 xval4
                        csegt tval1 tval2 tval3 tval4 */
#define MSEC_CONVERT  1000

/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{
    char             *tmpstring;
    mxArray            *datastruct;
    mxArray            *mxtmp;
    double            *spiketimes;
    double            *startindex;
    double            *endindex;
    double            *timestep, *ntimesteps;
    double            *fixtime;
    double            *testID;
    double            *epsx;
    double            *epst;
    double            *cpx;
    double            *cpt;
    double            *thetahat;
    double            *thetahatptr;
    double            *thetainit;
    double            *thetacurrx;
    double            *thetacurrt;
    double            *partialsx;
    double            *partialst;
    double            tmp;
    double            LambdaT;
    double            LambdaTInt;
    double            Delta;
    double            x1, x2, x3;
    double            *xtmp;
    double            t1, t2, t3;
    double            *ttmp;
    double            lastspiketime;
    double            isi;
    double            spaceval;
    double            timeval;
    double            *KSStat;
    double            *xspacing, *tspacing;
    double		*tmpptr;

    int            ntrials;
    int		   totaltimesteps = 0;
    int            maxspikes;
    int            trial;

    int            nspikes;
    int            nextspike;
    int            csegx;
    int            csegt;
    int            ncpx;
    int            ncpt;
    int            dNt;
    int            i, j, t;
    int            startt;

    
    /* initialized the necessary variables */
    tmpstring = (char *) mxCalloc(sizeof(char), 80);
    partialsx = (double *) mxCalloc(sizeof(double), 4);
    partialst = (double *) mxCalloc(sizeof(double), 4);
    xtmp = (double *) mxCalloc(sizeof(double), 4);
    ttmp = (double *) mxCalloc(sizeof(double), 4);

      /* Check numbers of arguments */
      if (nrhs != 11) {
           mexErrMsgTxt("Usage: [thetainit thetahat KSStat] = adaptivesuzuki(spiketimes, startindex, endindex, timestep, ntimesteps,  fixtime, testID, epsx, epst, cpx, cpt)" );
    }
      if (nlhs < 3) {
           mexErrMsgTxt("Usage: [thetainit thetahat KSStat] = adaptivesuzuki(spiketimes, startindex, endindex, timestep, ntimesteps,  fixtime, testID, epsx, epst, cpx, cpt)" );
    }

    /* get the list of spiketimes */
    spiketimes = mxGetPr(prhs[0]);
    nspikes = mxGetM(prhs[0])*mxGetN(prhs[0]);

    /* get the list of spikes indices */
    startindex = mxGetPr(prhs[1]);
    endindex = mxGetPr(prhs[2]);
    ntrials = mxGetM(prhs[2])*mxGetN(prhs[2]);

    /* get the list of times */
    timestep = mxGetPr(prhs[3]);
    ntimesteps = mxGetPr(prhs[4]);
    for (i = 0; i < ntrials; i++) {
        totaltimesteps += round(ntimesteps[i]);
    }

    /* Get the onset of fixation for each trial */
    fixtime = mxGetPr(prhs[5]);
    /* convert the fixation times to millisecond */
    //for (i = 0; i < ntrials; i++) {
    //    fixtime[i] /= MSEC_CONVERT;
    //}

    /* Get the valid trial array for this adaptive intensity function */
    testID = mxGetPr(prhs[6]);

    /* get the epsilons */
    epsx = mxGetPr(prhs[7]);
    epst = mxGetPr(prhs[8]);

    /* get the coordinates of the control points */
    cpx = mxGetPr(prhs[9]);
    ncpx = mxGetN(prhs[9]) * mxGetM(prhs[9]);
    cpt = mxGetPr(prhs[10]);
    ncpt = mxGetN(prhs[10]) * mxGetM(prhs[10]);
    xspacing = (double *) mxCalloc(ncpx, sizeof(double));
    tspacing = (double *) mxCalloc(ncpt, sizeof(double));

    /* create the output datastructures */
    plhs[0] = mxCreateDoubleMatrix(ncpx + ncpt, 1,  mxREAL);
    thetainit = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(NTHETAP, totaltimesteps,  mxREAL);
    thetahat = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, nspikes,  mxREAL);
    KSStat = mxGetPr(plhs[2]);

    /* set all of the KSStats to -1 */
    for (i = 0; i < nspikes; i++) {
        KSStat[i] = -1;
    }

    /* create the thetacurr array to hold the current value of thetahat x and t */
    thetacurrx = (double *) mxCalloc(ncpx, sizeof(double));
    thetacurrt = (double *) mxCalloc(ncpt, sizeof(double));

    /* set the initial thetahat variables */

    //for (i = 0; i < ncpx; i++) thetainit[i] = thetacurrx[i] = nspikes/(totaltimesteps*timestep[0]/1000);
    //set the initial x variable to be the averages, over the first valid trial, of the
    //firing rate for the various segments
    trial = 0;
    while (!testID[trial]) {
	trial++;
    }
    nextspike = round(startindex[trial]);
    /* advance to the first spike */
    while ((nextspike < round(endindex[trial])) && 
	   (spiketimes[nextspike] < 0)) {
	nextspike++;
    }
    /* get the number of spikes in the 0 - 300 ms delay period */
    dNt = 0;
    while ((nextspike < round(endindex[trial])) && 
	   (spiketimes[nextspike] < 300)) {
	   dNt++;
	   nextspike++;
    }
    for (i = 0; cpx[i] < 300; i++) {
        thetainit[i] = thetacurrx[i] = (float) dNt / 0.3;
    }
    /* get the number of spikes in the 300 - 800 ms stimulus presentation period */
    dNt = 0;
    while ((nextspike < round(endindex[trial])) && 
	   (spiketimes[nextspike]< 800)) {
	   dNt++;
	   nextspike++;
    }
    for (; cpx[i] < 800; i++) {
        thetainit[i] = thetacurrx[i] = (float) dNt / 0.5;
    }
    /* get the number of spikes in the 800 - 1500 ms stimulus presentation period */
    dNt = 0;
    while ((nextspike < round(endindex[trial])) && 
	   (spiketimes[nextspike] < 1500)) {
	   dNt++;
	   nextspike++;
    }
    for (; cpx[i] < 1500; i++) {
        thetainit[i] = thetacurrx[i] = (float) dNt / 0.7;
    }
    /* get the number of spikes in the 1500 - 2000 ms stimulus presentation period */
    dNt = 0;
    while ((nextspike < round(endindex[trial])) && 
	   (spiketimes[nextspike++] < 2000)) {
	   dNt++;
    }
    for (; i < ncpx; i++) {
        thetainit[i] = thetacurrx[i] = dNt / 0.5;
    }



    for (i = 0; i < ncpt; i++) thetainit[i+ncpx] = thetacurrt[i] = 1;

    /* Start the forward adaptive algorithm */
    thetahatptr = thetahat;
    for (trial = 0; trial < ntrials; trial++) {
	if (testID[trial]) {
	    lastspiketime = -1;
	    t = 0;
	    tmpptr = thetahatptr;
	    /* advance to the fixation time */
	    nextspike = round(startindex[trial]);
	    while ((nextspike < round(endindex[trial])) && 
	           (spiketimes[nextspike] < 0)){
		nextspike++;
            }
	    for (t = 0; t < ntimesteps[trial]; t++) {
	        /* find out how many spikes were present in the current timestep */
  	        dNt = 0;	  
    	        while ((nextspike < endindex[trial]) && 
		       (spiketimes[nextspike] <= t)){
		    nextspike++;
		    dNt++;
	        }
	        /* find the segment of the spline corresponding to the current 
	         * position. */
	        csegx = 1;
	        while (cpx[csegx+1] < t) {
                    csegx++;
	        }
	        x1 = (double) (t - cpx[csegx]) / (double) (cpx[csegx+1] - cpx[csegx]);
	        x2 = x1 * x1;
	        x3 = x2 * x1;
	        /* calculate x * Mc for the cardinal spline */
	        xtmp[0] = (-.5*x3 + x2 -.5*x1);
	        xtmp[1] = (1.5*x3 - 2.5*x2 + 1);
	        xtmp[2] = (-1.5*x3 + 2.0*x2 + .5*x1);
	        xtmp[3] = (.5*x3 - .5*x2); 
	        spaceval = (thetacurrx[csegx-1]*xtmp[0] + thetacurrx[csegx]*xtmp[1] +
	    		  thetacurrx[csegx+1]*xtmp[2] + thetacurrx[csegx+2]*xtmp[3]);
	        LambdaT = spaceval;
  
  
    	        /* if the last spike time was valid, add in the temporal spline element */
	        if (lastspiketime > 0) {
                    /* find the segment of the spline corresponding to the current 
		     * time. */
		    csegt = 1;
		    isi = t - lastspiketime;
		    while (cpt[csegt+1] < isi) {
                        csegt++;
		    }
		    t1 = (double) (isi - cpt[csegt]) / (double) (cpt[csegt+1] - cpt[csegt]);
		    t2 = t1 * t1;
		    t3 = t2 * t1;
		    /* calculate t * Mc for the cardinal spline */
		    ttmp[0] = (-.5*t3 + t2 -.5*t1);
		    ttmp[1] = (1.5*t3 - 2.5*t2 + 1);
		    ttmp[2] = (-1.5*t3 + 2.0*t2 + .5*t1);
		    ttmp[3] = (.5*t3 - .5*t2); 
		    timeval = thetacurrt[csegt-1]*ttmp[0] + thetacurrt[csegt]*ttmp[1] +
		              thetacurrt[csegt+1]*ttmp[2] + thetacurrt[csegt+2]*ttmp[3];
		    LambdaT = LambdaT * timeval;
		    /*LambdaTInt += LambdaT*timestep[0]; */
		    LambdaTInt += LambdaT*1e-3;
		    if (dNt) {
                        /* we need to subtract one because the nextspike index starts at 1
		         * for the first spike */
                        KSStat[nextspike-1] = 1-exp(-LambdaTInt); 
		        LambdaTInt = 0; 
		    }    
	        }

	        //Delta = dNt - (timestep[0] * LambdaT);
	        Delta = dNt - LambdaT * .001;

	        partialsx[0] = Delta * xtmp[0];
	        partialsx[1] = Delta * xtmp[1];
	        partialsx[2] = Delta * xtmp[2];
	        partialsx[3] = Delta * xtmp[3];


	        /* set the new values for the current segments */
		*(thetahatptr++) = csegx;
	        //thetahat[trial*(int) ntimesteps[trial]*NTHETAP+t*NTHETAP] = csegx;
	        for (i = 0; i < 4; i++) {
                    tmp = thetacurrx[csegx-1+i] + epsx[i] * partialsx[i];
		    //thetacurrx[csegx-1+i] = thetahat[trial*(int) ntimesteps[0]*NTHETAP+t*NTHETAP+i+1] = tmp; 
		    thetacurrx[csegx-1+i] = *(thetahatptr++) = tmp; 
                }

	        /* if the isi is valid, update thetacurrt */
	        if (lastspiketime > 0) {
		    partialst[0] = Delta * ttmp[0];
		    partialst[1] = Delta * ttmp[1];
		    partialst[2] = Delta * ttmp[2];
		    partialst[3] = Delta * ttmp[3];

		    /* set the new values for the current segments */
		    //thetahat[trial*(int) ntimesteps[0]*NTHETAP+t*NTHETAP+5] = csegt;
		    *(thetahatptr++) = csegt;
		    for (i = 0; i < 4; i++) {
		        tmp = thetacurrt[csegt-1+i] + epst[i] * partialst[i];
		        //thetacurrt[csegt-1+i] = thetahat[trial*(int) ntimesteps[0]*NTHETAP+t*NTHETAP+i+6] = tmp; 
		        thetacurrt[csegt-1+i] = *(thetahatptr++) = tmp; 
                    }
	        }
	        else {
		    /* ignore the current teporal spline point */
		    *thetahatptr = -1;
		    thetahatptr += 5;
	        }
	        if (dNt) {
	  	    lastspiketime = t;
	        }
	    }
	        
        }
	else {
	    thetahatptr += round(ntimesteps[trial]) * NTHETAP;
	    //mexPrintf("trial %d, advancing thetahatptr by %d, offset %d\n", trial, round(ntimesteps[trial]) * NTHETAP, thetahatptr - thetahat);
	}

    }
    mxFree(thetacurrx);
    mxFree(thetacurrt); 
    mxFree(tmpstring);
    mxFree(partialsx);
    mxFree(partialst);
    mxFree(xtmp);
    mxFree(ttmp); 
    return;
}
