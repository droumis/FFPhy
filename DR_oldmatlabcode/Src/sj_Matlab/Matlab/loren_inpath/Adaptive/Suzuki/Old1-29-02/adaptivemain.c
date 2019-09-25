/*
    adaptivemain.c 
    2 D version for spatial and temporal spline

    This version assumes that the control point list convers both the
    foward and backward directions of movement along the U track and, for
    movement along the second trajectory (from endpoint 2 to endpoint 1)
    adjusts the position to lie along the second set of points by adding the
    distance to the center control point

    To get the first pass estimate, we run the data backwards 

*/

#include "mex.h"
#include "mexheader.h"

#define NTHETAP 10 /* store 10 points at each timestep: csegx xval1 xval2 xval3 xval4
                        csegt tval1 tval2 tval3 tval4

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
    double            *timearray;
    double            *epsx;
    double            *epst;
    double            *cpx;
    double            *cpt;
    double            *thetahat;
    double            *thetainit;
    double            *thetacurrx;
    double            *thetacurrt;
    double            *partialsx;
    double            *partialst;
    double            tmp;
    double            LambdaT;
    double            LambdaTInt;
    double            Delta;
    double            timestep;
    double            x1, x2, x3;
    double            *xtmp;
    double            t1, t2, t3;
    double            *ttmp;
    double            lastspiketime;
    double            currenttime;
    double            isi;
    double            spaceval;
    double            timeval;
    double            *KSStat;
    double            *xspacing, *tspacing;

    int            ntrials;
    int            maxspikes;
    int            trial;

    int            nspikes;
    int        	   ntimesteps;
    int            nextspike;
    int            csegx;
    int            csegt;
    int            ncpx;
    int            ncpt;
    int            dNt;
    int            i, j, t;

    
    /* initialized the necessary variables */
    tmpstring = (char *) mxCalloc(sizeof(char), 80);
    partialsx = (double *) mxCalloc(sizeof(double), 4);
    partialst = (double *) mxCalloc(sizeof(double), 4);
    xtmp = (double *) mxCalloc(sizeof(double), 4);
    ttmp = (double *) mxCalloc(sizeof(double), 4);

      /* Check numbers of arguments */
      if (nrhs != 8) {
           mexErrMsgTxt("Usage: [thetainit thetahat KSStat] = adaptivemain(spiketimes, startindex, endindex, timearray,  epsx, epst, cpx, cpt)" );
    }
      if (nlhs < 3) {
           mexErrMsgTxt("Usage: [thetainit thetahat KSStat] = adaptivemain(spiketimes, startindex, endindex, timearray,  epsx, epst, cpx, cpt)" );
    }

    /* get the list of spiketimes */
    spiketimes = mxGetPr(prhs[0]);
    nspikes = mxGetM(prhs[0])*mxGetN(prhs[0]);

    /* get the list of spikes indices */
    startindex = mxGetPr(prhs[1]);
    endindex = mxGetPr(prhs[2]);
    ntrials = mxGetM(prhs[2])*mxGetN(prhs[2]);

    /* get the list of times */
    timearray = mxGetPr(prhs[3]);
    ntimesteps = mxGetN(prhs[3]) * mxGetM(prhs[3]);

    /* get the epsilons */
    epsx = mxGetPr(prhs[4]);
    epst = mxGetPr(prhs[5]);

    /* get the coordinates of the control points */
    cpx = mxGetPr(prhs[6]);
    ncpx = mxGetN(prhs[6]) * mxGetM(prhs[6]);
    cpt = mxGetPr(prhs[7]);
    ncpt = mxGetN(prhs[7]) * mxGetM(prhs[7]);
    xspacing = (double *) mxCalloc(ncpx, sizeof(double));
    tspacing = (double *) mxCalloc(ncpt, sizeof(double));

    /* create the output datastructures */
    plhs[0] = mxCreateDoubleMatrix(ncpx + ncpt, 1,  mxREAL);
    thetainit = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(NTHETAP, ntrials*ntimesteps,  mxREAL);
    thetahat = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, nspikes,  mxREAL);
    KSStat = mxGetPr(plhs[2]);

    /* create the thetacurr array to hold the current value of thetahat x and t */
    thetacurrx = (double *) mxCalloc(ncpx, sizeof(double));
    thetacurrt = (double *) mxCalloc(ncpt, sizeof(double));

    timestep = timearray[1] - timearray[0];

    /* set the initial thetahat variables */
    //mexPrintf("init = %f, ntimsteps = %d\n",nspikes/(ntimesteps*ntrials*timestep),ntimesteps);  

    for (i = 0; i < ncpx; i++) thetainit[i] = thetacurrx[i] = nspikes/(ntimesteps*ntrials*timestep);
    for (i = 0; i < ncpt; i++) thetainit[i+ncpx] = thetacurrt[i] = 1;

    /* Start the forward adaptive algorithm */
    for (trial = 0; trial < ntrials; trial++) {
      lastspiketime = -1;
      nextspike = (int) startindex[trial];
      for (t = 0; t < ntimesteps; t++) {
          currenttime = timearray[t];
          /* find out how many spikes were present in the current timestep */
          dNt = 0;	  
          while ((nextspike < endindex[trial]) && (spiketimes[nextspike]/1000 <= currenttime)){
              nextspike++;
              dNt++;
          }

          /* find the segment of the spline corresponding to the current 
          * position. */
          csegx = 1;
          while (cpx[csegx+1] < currenttime) {
              csegx++;
          }
          x1 = (double) (currenttime - cpx[csegx]) / (double) (cpx[csegx+1] - cpx[csegx]);
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
          //mexPrintf("csegx = %d, x1 = %f, LambdaT = %f\n",csegx,x1,LambdaT);  


          /* if the last spike time was valid, add in the temporal spline element */
          if (lastspiketime > 0) {
              /* find the segment of the spline corresponding to the current 
               * time. */
              csegt = 1;
              isi = currenttime - lastspiketime;
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
              LambdaTInt += LambdaT*timestep;
              if (dNt) {
                  /* we need to subtract one because the nextspike index starts at 1
                   * for the first spike */
                  KSStat[nextspike-1] = 1-exp(-LambdaTInt); 
                  LambdaTInt = 0; 
              }    
          }

          Delta = dNt - (timestep * LambdaT);

          partialsx[0] = Delta * xtmp[0];
          partialsx[1] = Delta * xtmp[1];
          partialsx[2] = Delta * xtmp[2];
          partialsx[3] = Delta * xtmp[3];

          //if (dNt != 0 ) mexPrintf("lambdaT = %f, Delta = %f\n",LambdaT, Delta);  

          /* set the new values for the current segments */
          thetahat[trial*ntimesteps*NTHETAP+t*NTHETAP] = csegx;
          for (i = 0; i < 4; i++) {
              tmp = thetacurrx[csegx-1+i] + epsx[i] * partialsx[i];
              thetacurrx[csegx-1+i] = thetahat[trial*ntimesteps*NTHETAP+t*NTHETAP+i+1] = tmp; 
          }
          /* if the isi is valid, update thetacurrt */
          if (lastspiketime > 0) {
              partialst[0] = Delta * ttmp[0];
              partialst[1] = Delta * ttmp[1];
              partialst[2] = Delta * ttmp[2];
              partialst[3] = Delta * ttmp[3];

              /* set the new values for the current segments */
              thetahat[trial*ntimesteps*NTHETAP+t*NTHETAP+5] = csegt;
              for (i = 0; i < 4; i++) {
                  tmp = thetacurrt[csegt-1+i] + epst[i] * partialst[i];
                  thetacurrt[csegt-1+i] = thetahat[trial*ntimesteps*NTHETAP+t*NTHETAP+i+6] = tmp; 
              }
          }
          else {
              /* ignore the current teporal spline point */
              thetahat[trial*ntimesteps*NTHETAP+t*NTHETAP+5] = -1;
          }
          if (dNt) {
              lastspiketime = currenttime;
          }
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
