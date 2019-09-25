/*
    readeeg.c 

    Reads in eeg data from a single channel eeg file between tstart and tend
    and puts it into a structure 

*/

#include "mexheader.h"
#include "mex.h"
#include "fileiomat.h"
#include <unistd.h>



#define NDIMS    1     /* number of dimensions in the data file */

/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{
    FILE             *efile;
    char             *eegfile;
    char            *timestring1;
    char            *timestring2;
    char            tmpstring[200];
    const char         *fieldnames[] = {"descript", "fields", "starttime", "samprate", "data"};
    u32     tstart, tend;
    unsigned int    stringlen;
    ContRec          tmpcontrec, emptycontrec;
    mxArray            *data;
    mxArray            *datatmp;
    mxArray            *datastruct;
    mxArray            *mxtmp;
    double            *tmpdata;
    u32   lasttime;
    double            *eegdata;
    double            starttime;
    double            startoffset;
    double            reclen;
    double            sampfreq;
    int             dims[2]= {1,1};
    int             ncontrec;
    int             offset;
    int             currenti;
    int             i, j, k, numinterp, nrec;
    

      /* Check numbers of arguments */
      if ((nrhs != 3)) {
           mexErrMsgTxt("Usage: eegrec = readeeg(eegfile, tstart, tend)");
    }
      if (nlhs < 1) {
        mexErrMsgTxt("readeeg must be called with one output argument");
      }


    /* get the file name */
    stringlen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
    eegfile = (char *) mxCalloc(stringlen, sizeof(char));
    mxGetString(prhs[0], eegfile, stringlen);

    /* get tstart */
    stringlen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
    timestring1 = (char *) mxCalloc(stringlen, sizeof(char));
    mxGetString(prhs[1], timestring1, stringlen);
    tstart = ParseTimestamp(timestring1);

    /* get tend */
    stringlen = (mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
    timestring2 = (char *) mxCalloc(stringlen, sizeof(char));
    mxGetString(prhs[2], timestring2, stringlen);
    tend = ParseTimestamp(timestring2);


    /* open the eeg file */
    if ((efile = fopen(eegfile, "r")) == NULL) {
        sprintf(tmpstring, "Error opening %s for reading\n", eegfile);
        mexErrMsgTxt(tmpstring);
    }

    /* get rid of the header */    
    do {
        fgets(tmpstring, 200, efile);
    } while ((strncmp(tmpstring, "%%ENDHEADER", 11) != 0));

    /* get the first record of the file */
    readcontrec(efile, &tmpcontrec);
    
    /* create an empty cont rec to fill in any holes in the data */
    emptycontrec.numsamples =  tmpcontrec.numsamples;
    emptycontrec.sampfreq =  tmpcontrec.sampfreq;
    mexPrintf("nsamp = %d, sampfreq = %f\n", tmpcontrec.numsamples, tmpcontrec.sampfreq);
    /* allocate space for the data */
    emptycontrec.data = (short *) mxCalloc(emptycontrec.numsamples, sizeof(short));
    for (i = 0; i < emptycontrec.numsamples; i++) {
	emptycontrec.data[i] = SHRT_MIN;
    }

    /* calculate the length (in time) of a record */
    reclen = (tmpcontrec.numsamples / tmpcontrec.sampfreq) * TIMESCALE;

    /* find the offsets of the first record to be included. */
    while ((tmpcontrec.timestamp < (tstart - reclen)) && (!feof(efile))) {
        readcontrec(efile, &tmpcontrec);
    mexPrintf("nsamp = %d, sampfreq = %f\n", tmpcontrec.numsamples, tmpcontrec.sampfreq);
    } 
    if (feof(efile)) {
        sprintf(tmpstring, "Error in start time %s\n", timestring1);
        mexErrMsgTxt(tmpstring);
    }
    starttime = (double) tmpcontrec.timestamp / TIMESCALE;
    lasttime = tmpcontrec.timestamp;

    /* put the eeg data2 into a Matlab structure with five fields */
    datastruct = mxCreateStructArray(2, dims, 5, fieldnames);

    /* set the first two fields */
    sprintf(tmpstring, "eeg data from %s: %s to %s\n", eegfile,
            timestring1, timestring2);
    mxtmp = mxCreateString((const char *)tmpstring); 
    mxSetField(datastruct, 0, fieldnames[0], mxtmp);
    sprintf(tmpstring, "eegamplitude");
    mxtmp = mxCreateString((const char *)tmpstring); 
    mxSetField(datastruct, 0, fieldnames[1], mxtmp);


    /* put the start time in the starttime field */
    datatmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    tmpdata = mxGetPr(datatmp);
    *tmpdata = (double) tmpcontrec.timestamp / TIMESCALE;
    starttime = *tmpdata;
    mxSetField(datastruct, 0, fieldnames[2], datatmp);

    /* calculate the number of records that need to be read in. Note that
     * because the clock cards run at slightly different frequencies than the
     * data acquisition cards, we add in a fudge factor of 1% to make sure we
     * have enough records */
    ncontrec =  (int) (ceil(((double) (tend - tstart)) / reclen) * 1.01);

    /* allocate space for the list of eeg samples. As the reclen is not exact
     * due to rounding error in the sampling frequency, round the sampling
     * frequency upward to get an upper bound on the number of samples */

    eegdata = (double *) mxCalloc(ncontrec * tmpcontrec.numsamples,
                                    sizeof(double));

    /* copy the first set of samples to eegdata */
    currenti = 0;
    for (i = 0 ; i < tmpcontrec.numsamples; i++) {
        eegdata[currenti++] = (double) tmpcontrec.data[i];
    }
    
    /* read in the rest of the records and put the eeg data into tmpdata */
    nrec = 1;
    while(!feof(efile)) {
        /* break out of the loop if the timestamp is > tend */
	if ((readcontrec(efile, &tmpcontrec) == 0) || 
	    (tmpcontrec.timestamp > tend)) {
            break;
        }
        /* check for "holes" in the data where the difference between two
         * adjacent timestamps is > 1.5 * reclen */
        if ((tmpcontrec.timestamp - lasttime) < (.5 * reclen)) {
	    mexPrintf("short time %ld\n", tmpcontrec.timestamp - lasttime);
	}
        if ((tmpcontrec.timestamp - lasttime) > (1.5 * reclen)) {
            sprintf(tmpstring, "Error: excessive difference between current and last timestamps: %ld vs %ld\nFilling in difference with SHRT_MIN\n", tmpcontrec.timestamp, lasttime);
	    mexPrintf(tmpstring);
	    numinterp = round((double)(tmpcontrec.timestamp - lasttime) / 
		              reclen) - 1;
	    for (k = 0; k < numinterp; k++) {
		emptycontrec.timestamp = lasttime + (u32) reclen;
		/* copy the samples to tmpdata */
		nrec++;
		for (i = 0 ; i < emptycontrec.numsamples; i++) {
		    eegdata[currenti++] = (double) emptycontrec.data[i];
		}
		lasttime += round(reclen);
	    }
        }
	/* copy the samples to tmpdata */
	nrec++;
	for (i = 0 ; i < tmpcontrec.numsamples; i++) {
	    eegdata[currenti++] = (double) tmpcontrec.data[i];
	}
	lasttime = tmpcontrec.timestamp;
    }    
    /* allocate space for the data */
    data = mxCreateDoubleMatrix(currenti, NDIMS, mxREAL);
    tmpdata = mxGetPr(data);

    /* copy the eeg data into tmpdata */
    for (i = 0 ; i < currenti; i++) {
        tmpdata[i] = eegdata[i];
    }

    /* set the data field */
    mxSetField(datastruct, 0, fieldnames[4], data);
    plhs[0] = datastruct;

    /* recalculate the sampling frequency based on the total time from starttime
     * to lasttime+reclen to correct for rounding errors */
    sampfreq = currenti / ((double)(lasttime + reclen) / TIMESCALE - starttime);
/*    mexPrintf("comparison of file versus computed sampfreq: %lf vs. %lf\n",
               tmpcontrec.sampfreq, sampfreq);
    mexPrintf("%d,  %lf, %lf, %lf\n", currenti, (double) lasttime/ TIMESCALE,
                (double)reclen/ TIMESCALE,  starttime); */
    /* put the sampling frequency in the sampfreq field */
    datatmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    tmpdata = mxGetPr(datatmp);
    *tmpdata = sampfreq; 
    mxSetField(datastruct, 0, fieldnames[3], datatmp);

    /* free up the allocated variables  */
    mxFree(eegfile);
    mxFree(timestring1);
    mxFree(timestring2);  
    fclose(efile); 
}







