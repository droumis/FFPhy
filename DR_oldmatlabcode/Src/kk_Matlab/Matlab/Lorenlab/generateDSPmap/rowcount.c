

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>



void Usage(void) ;


/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{

      double 			  *dblptr;
      int               i,j,k,z;
      double      *values;
      double      *lookupvalues;
      double			  *output;
      double   tmpcount;
      int				nvalues  , nlookupvalues, nrows, ncols, nlookuprows, nlookupcols, issame;

      /* Check numbers of arguments */
      if (!(nrhs == 2) || !((nlhs == 1)||(nlhs == 0)) ){
         Usage();
      }
   

      values = mxGetData(prhs[0]);      
      nrows = mxGetM(prhs[0]);
      ncols = mxGetN(prhs[0]);
      lookupvalues = mxGetData(prhs[1]);      
      nlookuprows = mxGetM(prhs[1]);
      nlookupcols = mxGetN(prhs[1]);
      
      if (!(ncols == nlookupcols)&&(nrows > 0) && (nlookuprows > 0)){
         mexErrMsgTxt("The number of columns in both inputs must be the same.\n");
      }  
    
 
      output = (double *) mxCalloc(nrows, sizeof(double));

      for (i = 0; i < nrows; i++) {
         tmpcount = 0;
         for (j = 0; j < nlookuprows; j++) {
            issame = 1;
            for (k = 0; k < ncols; k++) { 
               if (values[i+(nrows*k)] != lookupvalues[j+(nlookuprows*k)]) {
                  issame = 0;
               }
            }
            if (issame) {
               tmpcount++;
               
            }   
         }
         output[i] = tmpcount;
         
      }

      plhs[0] = mxCreateDoubleMatrix(nrows,1, mxREAL);

      dblptr = mxGetPr(plhs[0]);

      for (i = 0; i < nrows; i++) {
               *dblptr = output[i];
               dblptr++;
      }


      mxFree(output);

      return;
}

void Usage(void)
{
      mexErrMsgTxt("Usage: one putput and two inputs\n");
}


