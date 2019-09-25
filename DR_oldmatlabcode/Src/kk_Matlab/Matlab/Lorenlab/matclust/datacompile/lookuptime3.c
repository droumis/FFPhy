

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
	short int		  tempmax;
	double      *spiketimes;
	double      *postimes;
        double         *pos;
	int			  tmpindex, bound1, bound2;
	double			  *output;
        double                    *output2;
        double                     timeratio;
	int				  nspiketimes, npostimes, lastnonzero, nextnonzero;

    /* Check numbers of arguments */
    if (!(nrhs == 3) || !(nlhs == 2)){
        Usage();
    }
   





	spiketimes = mxGetPr(prhs[0]);
	nspiketimes = mxGetM(prhs[0]) * mxGetN(prhs[0]);

	postimes = mxGetPr(prhs[1]);
	npostimes = mxGetM(prhs[1]) * mxGetN(prhs[1]);
         
        pos = mxGetPr(prhs[2]);
        
	
	

	output = (double *) mxCalloc(3* nspiketimes, sizeof(double));
        output2 = (double *) mxCalloc(nspiketimes, sizeof(double));
	

	

	for (i = 0; i < nspiketimes; i++) {
	  if ((spiketimes[i] >= postimes[0]) && (spiketimes[i] <= postimes[npostimes-1])) {
		
		tmpindex = (npostimes/2);
		bound1 = 0;
		bound2 = npostimes-1;
		
		while ((tmpindex != bound1) && (tmpindex != bound2)) {
			if (spiketimes[i] < postimes[tmpindex]) {
				bound2 = tmpindex;
				tmpindex = ((tmpindex-bound1)/2) + bound1;
			}
			else {
				bound1 = tmpindex;
				tmpindex = ((bound2-tmpindex)/2)+tmpindex;
			}
		}
	        
                
                output[3*i] = pos[3*tmpindex];
                output[3*i+1] = pos[3*tmpindex+1];
                output[3*i+2] = pos[3*tmpindex+2];
                output2[i] = tmpindex + 1;
                

                   
                   
          }
	  else {
		  tmpindex = 1;
                  output[3*i] = 0;
                  output[3*i+1] = 0;
                  output[3*i+2] = 0;
                  
                 
	  }
	  if (tmpindex == 0) {
             tmpindex = 1;
             output[3*i] = 0;
             output[3*i+1] = 0;
             output[3*i+2] = 0;
           
  
          }		
			




		
	}
		
		
					

	plhs[0] = mxCreateDoubleMatrix(3, nspiketimes, mxREAL);

	dblptr = mxGetPr(plhs[0]);

	for (i = 0; i < 3*nspiketimes; i++) {
		 *dblptr = output[i];
		 dblptr++;
	}
        plhs[1] = mxCreateDoubleMatrix(nspiketimes, 1, mxREAL);

	dblptr = mxGetPr(plhs[1]);

	for (i = 0; i < nspiketimes; i++) {
		 *dblptr = output2[i];
		 dblptr++;
	}

        
        

	mxFree(output);

	return;
}

void Usage(void)
{
	mexErrMsgTxt("wrong usage usage\n");
}


