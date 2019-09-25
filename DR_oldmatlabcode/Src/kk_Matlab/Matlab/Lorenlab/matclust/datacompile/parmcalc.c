

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
	int			      triggered,tempmaxloc,triggernum,detriggernum,triggerchannel, detriggered, retriggered, retriggernum;
	short int		  tempmax, tempmin, maxtime, mintime, tempmax2;
	unsigned int      *times;
	short int         *waves;
	double            *triggers;
	double			  *output;

	int                 maxchannel;
        int				  ntimepoints;
	int				  nwavepoints;

	int			      UnitsPerSec = 10000;
	int			      SamplingRate = 31200;
	float			  TimePerStep = (1/SamplingRate)*UnitsPerSec;
	float			  TimePerWindow = TimePerStep*40;

    /* Check numbers of arguments */
    if (!(nrhs == 3) || !(nlhs == 1)){
        Usage();
    }
   


/* function call: parameters = parmcalc(timestamps,waves,triggers)  
waves is 40 by 4 by n matrix of 16 bit int
triggers is a double vector of 4 trigger numbers  
timestamps is an unsigned 32bit int vector of n timestamps*/


	times = mxGetData(prhs[0]);
	ntimepoints = mxGetM(prhs[0]) * mxGetN(prhs[0]);

	waves = mxGetData(prhs[1]);
	nwavepoints = mxGetM(prhs[1]) * mxGetN(prhs[1]);

	triggers = mxGetPr(prhs[2]);
	
        
	output = (double *) mxCalloc(ntimepoints*6, sizeof(double));

	j = 0;
	for (i = 0; i < ntimepoints; i++) {
		triggernum = 40;
		detriggernum = 39;
		triggered = 0;
		detriggered = 0;
                tempmax = -1;
		tempmax2 = -1;
                tempmaxloc = 40;
	        retriggered = 0;
		retriggernum = 40;
			
		
		/*find the the earliest time when one channel comes above the threshhold, where before all were below threshold
                  and the time when all channels come back down below threshold, and also which channel had the highest peak
                */ 
			
		for (k=0; k<40; k++) {
			if (k>3) {
				for (j = 0; j<4; j++) {			
                                        
                                        /* record the maximum channel amplitude until all channels are back below threshold*/
                                        if ((waves[40*4*i + 40*j + k] > tempmax2)&&(!detriggered)) {
                                           tempmax2 = waves[40*4*i + 40*j + k];
                                           maxtime = k;
                                           maxchannel = j;
                                        }
                                        
                                        /* if any channel comes above threshhold for the first time ...  */
                                        if ((waves[40*4*i + 40*j + k]>triggers[j])&&(waves[40*4*i + 40*0 + k - 2]<triggers[0])&&(waves[40*4*i + 40*1 + k - 2]<triggers[1])&&(waves[40*4*i + 40*2 + k - 2]<triggers[2])&&(waves[40*4*i + 40*3 + k - 2]<triggers[3])&&(k<triggernum)) {
						

						triggered = 1;
						triggernum = k;
						triggerchannel = j;
					}	
                                        
                                        /* if all channels have come back below threshhold after triggering and a channel is now triggering again ...  */
                                        if ((waves[40*4*i + 40*j + k]>triggers[j])&&(waves[40*4*i + 40*0 + k - 2]<triggers[0])&&(waves[40*4*i + 40*1 + k - 2]<triggers[1])&&(waves[40*4*i + 40*2 + k - 2]<triggers[2])&&(waves[40*4*i + 40*3 + k - 2]<triggers[3])&&(detriggered)&&(!retriggered)) {
						
                                                retriggered = 1;
                                                retriggernum = k;
					}	
                                
                                }
                                /* when all channels come back below threshhold, the spike has detriggered */
                                if ((waves[40*4*i + 40*0 + k]<triggers[0])&&(waves[40*4*i + 40*1 + k]<triggers[1])&&(waves[40*4*i + 40*2 + k]<triggers[2])&&(waves[40*4*i + 40*3 + k ]<triggers[3])&&(triggered)&&(!detriggered)) {
						
                                        detriggered = 1;
				        detriggernum = k;
			        }	
				         
                        
		        }
			
		}
		
		/*between the trigger time to either the retrigger time or the end of the 40 points, find the minimum point of the channel that had the highest peak*/
		tempmin = waves[40*4*i+40*maxchannel+triggernum];
                mintime = triggernum;
                for (k=triggernum; k<retriggernum; k++) {
			if (waves[40*4*i+40*maxchannel+k] < tempmin) {
                           tempmin = waves[40*4*i+40*maxchannel+k];
                           mintime = k;
                        }
		}

		
		/*between the trigger and detrigger times, find the maximum of each channel*/
		for (j = 0; j<4; j++) {
			tempmax = -32767;
			
			
			
			for (k=triggernum; k<detriggernum+1; k++) {
				if (k>3) {
					
					
                                      if ((waves[40*4*i + 40*j + k]>tempmax)&&(triggered)) {	
                                              tempmax = waves[40*4*i + 40*j + k];
                                              
					}
							
				}
			}
			
		        if (tempmax == -32767) {
                           tempmax = 0;
                        }
                        output[6*i + j] = (double) tempmax;
                }
		
                output[6*i + 4] = (mintime-maxtime);
                output[6*i + 5] = (tempmax2-tempmin);

	
	}

			

	plhs[0] = mxCreateDoubleMatrix(6,ntimepoints, mxREAL);

	dblptr = mxGetPr(plhs[0]);

	for (i = 0; i < 6*ntimepoints; i++) {
		 *dblptr = output[i];
		 dblptr++;
	}


	mxFree(output);

	return;
}

void Usage(void)
{
	mexErrMsgTxt("wrong usage usage\n");
}


