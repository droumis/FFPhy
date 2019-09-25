/*
   readeeg.c 

   Reads in eeg data from a single channel eeg file between tstart and tend
   and puts it into a structure 

*/

#include "mexheader.h"
#include "mex.h"
#include "fileiomat.h"
#include <unistd.h>


#define NDIMS   1    /* number of dimensions in the data file */

 /******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
   int           nlhs,          /* number of expected outputs */
   mxArray       *plhs[],       /* array of pointers to output arguments */
   int           nrhs,          /* number of inputs */
   const mxArray *prhs[])       /* array of pointers to input arguments */
{
   FILE         *efile;
   char         *eegfile;
   char        *timestring1;
   char        *timestring2;
   char        tmpstring[200];
   const char        *fieldnames[] = {"descript", "fields", "starttime", "samprate", "data"};
   unsigned long    tstart, tend;
   unsigned int   stringlen;
   ContRec         tmpcontrec, emptycontrec;
   mxArray           *data;
   mxArray           *datatmp;
   mxArray           *datastruct;
   mxArray           *mxtmp;
   double           *tmpdata;
   unsigned long  lasttime;
   double           *eegdata;
   double           starttime;
   double           startoffset;
   double           *channelprop;
   double           reclen;
   double           samplingrate;
   int            dims[2]= {1,1};
   int            ncontrec;
   int            offset;
   int            currenti;
   int            i, j, k, numinterp, nrec;
   

     /* Check numbers of arguments */
     if ((nrhs != 4)) {
         mexErrMsgTxt("Usage: eegrec = readeeg(eegfile, tstart, tend, toffset)");
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
   /* get the proportion of the way through the complete list of channels so
   * that we can compute the correct start time */
   channelprop = mxGetPr(prhs[3]);

   /* open the eeg file */
   if ((efile = fopen(eegfile, "r")) == NULL) {
      sprintf(tmpstring, "Error opening %s for reading\n", eegfile);
      mexErrMsgTxt(tmpstring);
   }

   /* get rid of the header */     
   do {
      fgets(tmpstring, 200, efile);
   } while (strncmp(tmpstring, "%%ENDHEADER", 11) != 0);

   /* get the first record of the file */
   readcontrec(efile, &tmpcontrec);
   fprintf(stderr, "%lf\n", tmpcontrec.samplingrate);
   
   /* create an empty cont rec to fill in any holes in the data*/
   emptycontrec.numsamples =  tmpcontrec.numsamples;
   emptycontrec.samplingrate =  tmpcontrec.samplingrate;
   /* allocate space for the data */
   emptycontrec.data = (short *) mxCalloc(emptycontrec.numsamples, sizeof(short));
   for (i = 0; i < emptycontrec.numsamples; i++) {
      emptycontrec.data[i] = SHRT_MIN;
   }

   /* calculate the length (in time) of a record */
   reclen = (tmpcontrec.numsamples / tmpcontrec.samplingrate) * TIMESCALE;

   /* find the offsets of the first record to be included. */
   while ((tmpcontrec.timestamp < (tstart - reclen)) && (!feof(efile))) {
      readcontrec(efile, &tmpcontrec);
   } 
   if (feof(efile)) {
      sprintf(tmpstring, "Error in start time %s\n", timestring1);
      mexErrMsgTxt(tmpstring);
   }
   /* the actual start time should be offset for each channel, as the channels 
   * are sampled sequentially. The offset is computed by taking the relative
   * offset of the channel in the list of channels divided by the total
   * number of channels (given by channelprop) and multiplying that by the
    * time for a complete sample of all channels (1/ sampling rate),
   * as that will give the time for the first sample from this channel */
   starttime = (double) tmpcontrec.timestamp / TIMESCALE + 
                   *channelprop / tmpcontrec.samplingrate;
   mexPrintf("channelprob / samplingrate = %lf\n", *channelprop / tmpcontrec.samplingrate);
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
         sprintf(tmpstring, "Error: excessive difference between current and last timestamps: %ld vs %ld; Filling in difference with SHRT_MIN", tmpcontrec.timestamp, lasttime);
         /*mexErrMsgTxt(tmpstring);*/  
         mexPrintf(tmpstring);
         numinterp = round((double)(tmpcontrec.timestamp - lasttime) / 
                       reclen) - 1;
         sprintf(tmpstring, "Numinterp: %d; Reclen: %f\n", numinterp, reclen);
         mexPrintf(tmpstring);
         for (k = 0; k < numinterp; k++) {
            /* emptycontrec.timestamp = lasttime + (unsigned long) reclen; */
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

   /*mexPrintf("lasttime = %ld\n", lasttime);*/
   /* recalculate the sampling frequency based on the total time from starttime
    * to lasttime+reclen to correct for rounding errors */
   samplingrate = currenti / ((double)(lasttime + reclen) / TIMESCALE - starttime);
 /*    mexPrintf("comparison of file versus computed samplingrate: %lf vs. %lf\n",
            tmpcontrec.samplingrate, samplingrate);
   mexPrintf("%d, %lf, %lf, %lf\n", currenti, (double) lasttime/ TIMESCALE,
            (double)reclen/ TIMESCALE, starttime); */
   /* put the sampling frequency in the samplingrate field */
   datatmp = mxCreateDoubleMatrix(1, 1, mxREAL);
   tmpdata = mxGetPr(datatmp);
   *tmpdata = samplingrate; 
   mxSetField(datastruct, 0, fieldnames[3], datatmp);

   /* free up the allocated variables  */
   mxFree(eegfile);
   mxFree(timestring1);
   mxFree(timestring2);  
   fclose(efile); 
}



int readcontrec(FILE *file, ContRec *contrec)
   /* returns 0 on error */
{
   unsigned long time;

   if (fread(&(contrec->timestamp), sizeof(unsigned long), 1, file) != 1) 
      return 0;
   if (fread(&(contrec->numsamples), sizeof(int), 1, file) != 1) 
      return 0;
   if (fread(&(contrec->samplingrate), sizeof(double), 1, file) != 1) 
     return 0;

   /* allocate space for the data */
   contrec->data = (short *) mxCalloc(contrec->numsamples, sizeof(short));
   if (fread(contrec->data, sizeof(short), contrec->numsamples, file) != 
      contrec->numsamples) 
      return 0;
      
   return 1;
}


unsigned long ParseTimestamp(char *s)
{
char  *ptr;
char  *ptr2;
unsigned long   time;
unsigned long hour;
unsigned long min;
unsigned long sec;
float fracsec;
int    ncolons;
char  *fracptr;
char  timestr[100];   

   if(s == NULL){
   return(0);
   }
   /*
   ** copy the passed argument to the timestring for
   ** manipulation
   */
   strcpy(timestr,s);
   /*
   ** check for hr:min:sec.fracsec format vs min:sec
   */
   ncolons = strcount(timestr,':');
   fracsec = 0;
   if((fracptr = strchr(timestr,'.')) != NULL){
   sscanf(fracptr,"%f",&fracsec);
   *fracptr = '\0';
   };
   switch(ncolons){
   case 0:
      if(fracptr){
         sscanf(timestr,"%ld",&sec);
         time = (unsigned long) (sec*1e4 + (fracsec*1e4 + 0.5)); 
      } else {
         /*
         ** straight timestamp
         */
         sscanf(timestr,"%ld",&time);
      }
      break;
   case 1:
      /*
      ** find the colon
      */
      ptr = strchr(timestr,':');
      /*
      ** separate the minutes and the seconds into two strings
      */
      *ptr = '\0';
      /*
      ** read the minutes before the colon
      */
      sscanf(timestr,"%ld",&min);
      /*
      ** read the seconds after the colon
      */
      sscanf(ptr+1,"%ld",&sec);
      /*
      ** compute the timestamp
      */
      time = (unsigned long) (min*6e5 + sec*1e4 + (fracsec*1e4 + 0.5)); 
      break;
   case 2:
      /*
      ** find the first colon
      */
      ptr = strchr(timestr,':');
      /*
      ** find the second colon
      */
      ptr2 = strchr(ptr+1,':');
      /*
      ** separate the hours, minutes and the seconds into strings
      */
      *ptr = '\0';
      *ptr2 = '\0';
      /*
      ** read the hours before the first colon
      */
      sscanf(timestr,"%ld",&hour);
      /*
      ** read the minutes before the second colon
      */
      sscanf(ptr+1,"%ld",&min);
      /*
      ** read the seconds after the colon
      */
      sscanf(ptr2+1,"%ld",&sec);
      /*
      ** compute the timestamp
      */
      time = (unsigned long) (hour*36e6 + min*6e5 + sec*1e4 + (fracsec*1e4+0.5)); 
      break;
   default:
      fprintf(stderr,"unable to parse timestamp '%s'\n",timestr);
      return(0);
   }
   return(time);
}

int IsStringEmpty(char *str)
{
   if(str == NULL) return(1);
   /*
   ** scan the string to see if there are any non-white space
   ** characters in it
   */
   while(str && (*str != '\0')){
   if((*str != ' ') && (*str != '\t') && (*str != '\n')){
      /*
      ** found a non-white space character
      */
      return(0);
   }
   str++;
   }
   /*
   ** all white space
   */
   return(1);
}

int strcount(char *s, char c)
{
int    count;

   if(s == NULL) return(0);
   count = 0;
   while(*s != '\0'){
      if(*s == c) count++;
      s++;
   }
   return(count);
}



