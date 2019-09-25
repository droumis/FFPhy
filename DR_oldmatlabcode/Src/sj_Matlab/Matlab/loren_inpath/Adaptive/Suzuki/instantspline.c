/*
	adaptivestats.c 
	2 D version for spatial and temporal spline

*/

#include "mex.h"
#include "mexheader.h"

#define TIMESTEP 1.0 /* Input data sampled at 2 ms intervals */
#define NTHETAP 10 /* store 10 points at each timestep: csegx xval1 xval2 xval3 xval4
						csegt tval1 tval2 tval3 tval4

/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
      	int		   nlhs,		   /* number of expected outputs */
       	mxArray	   *plhs[],		/* array of pointers to output arguments */
	int		   nrhs,		   /* number of inputs */
	const mxArray *prhs[])		 /* array of pointers to input arguments */
{
	double			*cpx;
	double			*cpt;
	double			*thetahat;
	double			*thetainit;
	double			*thetainitx;
	double			*thetainitt;
        double                  *timeinto;
	double                  *xvals;
	double                  *tvals;

	int                     ncpx;
	int                     ncpt;
	int                     ntimesteps;
	int                     csegx;
        int                     csegt;
	int                     t;
	int			i;

  	/* Check numbers of arguments */
  	if (nrhs != 5) {
	        mexErrMsgTxt("Usage: [xvals tvals] = instantspline(cpx, cpt, thetainit, thetahat, timeinto)" );
	}
	if (nlhs < 2) {
	        mexErrMsgTxt("Usage: [xvals tvals] = instantspline(cpx, cpt, thetainit, thetahat, timeinto)" );
        }

	/* get the list of control point positions */
	cpx = mxGetPr(prhs[0]);
	ncpx = mxGetM(prhs[0]) * mxGetN(prhs[0]);
	cpt = mxGetPr(prhs[1]);
	ncpt = mxGetN(prhs[1]) * mxGetM(prhs[1]);

	/* get the initial values of the control points */
	thetainit = mxGetPr(prhs[2]);
	thetainitx = (double *) mxCalloc(ncpx, sizeof(double));
	for ( i=0; i<ncpx ; i++ ) thetainitx[i] = thetainit[i];
        thetainitt = (double *) mxCalloc(ncpt, sizeof(double));
	for ( i=0; i<ncpt ; i++ ) thetainitt[i] = thetainit[i+ncpx];

	/* get the update structure for the splines */
	thetahat = mxGetPr(prhs[3]);
	ntimesteps = floor(mxGetN(prhs[3]));

	/* get the extraction time */
	timeinto = mxGetPr(prhs[4]);
	if (timeinto[0]/TIMESTEP < ntimesteps) { ntimesteps = floor(timeinto[0]/TIMESTEP); }

	/* create the output datastructures */
	plhs[0] = mxCreateDoubleMatrix(1, ncpx, mxREAL);
	xvals = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, ncpt, mxREAL);
	tvals = mxGetPr(plhs[1]);

	for ( i=0 ; i<ncpx ; i++ ) xvals[i] = thetainitx[i]; 

	for ( i=0 ; i<ncpt ; i++ ) tvals[i] = thetainitt[i]; 
       
	
	//mexPrintf("ncpx = %d, ncpt = %d, ntimesteps = %d\n", ncpx, ncpt, ntimesteps);

	for ( t=0 ; t<ntimesteps ; t++ ) {
		/* update splines to next sample time */
		csegx = thetahat[t*NTHETAP];
	        for (i = 0; i < 4; i++) {
			if (csegx >= 1) xvals[csegx-1+i] = thetahat[ i+1 + t*NTHETAP ];
		}
		//xvals[0] = xvals[1];
		//xvals[ncpx-1] = xvals[ncpx-2];

		csegt = thetahat[5+t*NTHETAP];
		for (i = 0; i < 4; i++) {
			if (csegt >= 1) tvals[csegt-1+i] = thetahat[ i+1 + 5 + t*NTHETAP ];
		}

	}
	
	mxFree( thetainitx );
	mxFree( thetainitt );
	return;
}



