

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
      double      *directionpr;
      double			  *output;
      double   *compvalpr;
      double    compareval, direction; 
      int	nvalues  , nlookupvalues, nrows, ncols, nlookuprows, nlookupcols, foundindex, lastvalidindex, findnearest, foundindex2;

      /* Check numbers of arguments */
      if (!((nrhs == 2)||(nrhs == 3)||(nrhs == 4)) || !((nlhs == 1)||(nlhs == 0)) ) {
         Usage();
      }
   

      values = mxGetData(prhs[0]);      
      nvalues = mxGetM(prhs[0])*mxGetN(prhs[0]);
      nrows = mxGetM(prhs[0]);
      ncols = mxGetN(prhs[0]);
      
      compvalpr = mxGetData(prhs[1]);
            
      
      if (nrhs > 2) {
         directionpr = mxGetData(prhs[2]);
         if ((mxGetM(prhs[2]) * mxGetN(prhs[2])) != 1) {
            mexErrMsgTxt("Usage: DIRECTION must be a single value\n");
         }
         if (!((directionpr[0] == -1)||(directionpr[0] == 1)||(directionpr[0] == 0))) {
            mexErrMsgTxt("Usage: DIRECTION must be either -1 or 0 or 1\n");
         }
         direction = directionpr[0];
      }
      else {
         direction = 0;
      }
            
            
         
      if (nrhs == 4) {
         lookupvalues = mxGetData(prhs[3]);      
         nlookupvalues = mxGetM(prhs[3]) * mxGetN(prhs[3]);
      }
      /*indeces were not defined, so we will use all indeces as default*/
      else {
         nlookupvalues = nvalues;
         lookupvalues = (double *) mxCalloc(nvalues, sizeof(double));
         for (i=0;i < nvalues; i++) {
            lookupvalues[i] = i+1;
         }
      }
       

      if ((mxGetM(prhs[1]) * mxGetN(prhs[1])) != 1) {
         mexErrMsgTxt("Usage: COMPARE VALUE must be a single value\n");
      }
 
      
      lastvalidindex = -1;
      
      
      if (direction == 0) {
         findnearest = 1;
         direction = 1;
      }
      else {
         findnearest = 0;
      }

      compareval = compvalpr[0];
      output = (double *) mxCalloc(nvalues, sizeof(double));
      for (i = 0; i < nvalues; i++) {
         output[i] = values[i];
      }
      
      if (direction > 0) {
         i = 0;
         
      }
      else {
         i = nlookupvalues - 1;
         
      }
      
      while ((i >= 0) && (i < nlookupvalues)) {
         
         if ((lookupvalues[i] < 1) | (lookupvalues[i] > nvalues)) {
            mexErrMsgTxt("Usage: an index value was out of bounds\n");
         }
         
         
         foundindex = 0;
         j = lookupvalues[i]-1;
         
         if (output[j] == compareval) {
            while ((j >= 0) && (j < nvalues)) {
               if (values[j] != compareval) {
                  foundindex = j;
                  break;
               }
               else {
                  j = j+direction;
               }
            }
            if (foundindex) {
               j = lookupvalues[i]-1;
               while (j != foundindex) {
                  foundindex2 = 0;
                  for (z = 0; z < nlookupvalues;z++) {
                     if (lookupvalues[z]-1 == j) {
                        foundindex2 = 1;
                        break;
                     }
                  }
                  if (foundindex2) {
                     if (findnearest) {
                        if ((abs(j-foundindex) < abs(j-lastvalidindex))||(lastvalidindex == -1)) {
                           output[j] = values[foundindex];
                        }
                        else {
                           output[j] = values[lastvalidindex];
                        }
                     }
                     else {                     
                        output[j] = values[foundindex];
                     }
                  }
                  j = j+direction;
                  
               }
            }
            else if ((findnearest) && (lastvalidindex != -1)) {
               j = lookupvalues[i]-1;
               output[j] = values[lastvalidindex];
            }
               
         
         }
         else {
            lastvalidindex = i;
         }
           
         
         i = i + direction;
         
      }

      plhs[0] = mxCreateDoubleMatrix(nrows,ncols, mxREAL);

      dblptr = mxGetPr(plhs[0]);

      for (i = 0; i < nvalues; i++) {
               *dblptr = output[i];
               dblptr++;
      }


      /*mxFree(output);
      mxFree(lookupvalues);*/

      return;
}

void Usage(void)
{
      mexErrMsgTxt("Usage: one putput and 2-4 inputs\n");
}


