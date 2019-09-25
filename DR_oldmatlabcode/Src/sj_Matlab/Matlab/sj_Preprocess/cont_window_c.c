/*
 * =============================================================
 *
 *  eegwindow_cpp.cpp
 *
 *  Caleb Kemere 2010-01-09
 *
 *  [data, timestamps]
 *    = eeg_window_cpp(filename, times, window_length);
 *
 *  [data, timestamps, samplingrate]
 *    = eeg_window_cpp(filename, times, window_length);
 *
 *       times, window_length in SECONDS
 *
 * =============================================================
 */

#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>

#include <time.h>

#include "mex.h"
#include "matrix.h"

#define MAX_FILENAME 200
#define HEADER_STRLEN 50
#define ERROR_STRLEN 500
#define MAX_NSAMPLES 5000

typedef uint32_t u32;

FILE *eegfile;

int nsamplesperbuf = 0;
int samplingrate = 0;
double actualSamplingRate;

double findRealTimestamp ();

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  char filename[MAX_FILENAME];
  int filenameLength=0;

  int nWindows;
  int windowLen; /* length of window in SAMPLES */
  double window[2];

  double *dEEGData, *dStartTime, *dStartActual;

  int16_t *data;

  char *header;
  char headerString[500];
  char tmpString[HEADER_STRLEN];

  fpos_t filePtr;
  size_t err;

  double offset;

  double recordTimestamps;
  double endRecord;

  uint32_t nextTimestamp;

  u32 dataStartTimestamp;

  int16_t *recordData;

  int nRecords;

  int i,j,odx;
  int temp;
  char errorString[ERROR_STRLEN];

  /* PROCESS INPUT ARGUMENTS */
  /* ----------------------------------------------------------*/
  /*  Check for proper number of arguments. */
  if (nrhs < 3) 
    mexErrMsgTxt("Missing arguments: filename, times, window-length required.");

  /* ================================================================
   *  Check to make sure file exists and is openable. 
   *   - Assume most other error checking taken care of by matlab
   *  wrapper, e.g., filename length, times is Nx1, verbose-flag
   *  provided.
   */
  
  mxGetString(prhs[0], filename, MAX_FILENAME-1);
  eegfile = fopen(filename, "rb");
  if (!eegfile) {
    sprintf(errorString,"Could not open file: %s; Return value %d\n", filename, (int) eegfile);
    mexErrMsgTxt(errorString);
  }

  /* = Read in header... ========================================== */
  fgets(tmpString,HEADER_STRLEN,eegfile);
  if (strcmp("%%BEGINHEADER\n",tmpString) != 0)
    mexErrMsgTxt("File has unknown header.");
  while (strcmp("%%ENDHEADER\n",tmpString) != 0) {
    fgets(tmpString,HEADER_STRLEN,eegfile);
    /* look for "% nsamplesperbuf %d" */
    if (strncmp(tmpString,"% nsamplesperbuf",16) == 0) {
      sscanf(tmpString,"%% nsamplesperbuf %d",&nsamplesperbuf);
    }
    if (strncmp(tmpString,"% samplingrate",14) == 0) {
      sscanf(tmpString,"%% samplingrate %d",&samplingrate);
    }
  }

  if ((samplingrate == 0) || (nsamplesperbuf == 0))
    mexErrMsgTxt("samplingrate or nsamplesperbuf not found.");

  /* ============================================================== */
  /* Initialization */
  /* mark start of file */
  fgetpos(eegfile, &filePtr);
  size_t dataStart = ftell(eegfile);

  nWindows = mxGetM(prhs[1]); /* length of times vector */

  err = fread(&nextTimestamp, sizeof(uint32_t), 1, eegfile);

  actualSamplingRate = 30000.0;
  i = 0;
  while (floor(actualSamplingRate/(i+1)) >= samplingrate)
    i = i + 1;
  actualSamplingRate = actualSamplingRate / i;

  /*fprintf(stderr,"Sampling rates: real = %f, actual = %f\n",
   * (double) samplingrate, actualSamplingRate);*/

  windowLen = ceil( mxGetScalar(prhs[2]) * (double) actualSamplingRate );
    /* IN SAMPLES!!! */

  /*   Now we can initialize our outputs...  */
  plhs[0] = mxCreateDoubleMatrix(windowLen,nWindows,mxREAL); /* output 1: eeg windows */
  plhs[1] = mxCreateDoubleMatrix(nWindows,1,mxREAL); /* output 2: window start timestamps */
  plhs[2] = mxCreateDoubleScalar(actualSamplingRate); /* output 3: sampling rate */
  recordTimestamps = (double) nsamplesperbuf * 10000.0 / actualSamplingRate; 

  dEEGData = (double *)mxGetData(plhs[0]);
  /* initialize to NaNs */
  double *dptr = dEEGData;
  for (j = 0; j < windowLen*nWindows; j++)
    *dptr++ = 0.0/0.0;

  dStartTime = (double *)mxGetData(prhs[1]);  /* These are what we're looking for... */
  dStartActual = (double *)mxGetData(plhs[1]); /* These are what we will find - should be similar.... */

  /* = Predict file records... ==================================== */
  struct stat filestat;
  stat(filename,&filestat);
  int recordLength = sizeof(u32) + sizeof(int16_t)*nsamplesperbuf;
  nRecords = (filestat.st_size - dataStart)/recordLength;
  if ( fmod(((double) filestat.st_size - dataStart), (double)recordLength) != 0) {
    mexWarnMsgTxt("File length does not match predicted.\n");
    fprintf(stdout,"%f %d %d %d %d\n",(double)filestat.st_size-dataStart,
        recordLength*nRecords, nRecords, recordLength,
        sizeof(int));
  }

  /* ============================================================== */
  /* Generate file index */

  endRecord = (nsamplesperbuf - 1) * 10000.0 / actualSamplingRate;
  recordData = (int16_t *)mxCalloc(nsamplesperbuf, sizeof(int16_t)); /* one record's data */

  fsetpos(eegfile,&filePtr); /* start of the file */

  /* ============================================================== */
  /* Ok, so we need to loop through the records and the window
     start times. Keep a look out for missing data at beginnings
     and ends, and for overlapping data - probably easiest to
     ftell and fseek for this.
     */
  for (i = 0; i < nWindows; i++) {
    fprintf(stderr,"+");
    fsetpos(eegfile, &filePtr); /* reset to the first record */
    dptr = dEEGData + i * windowLen;
    window[0] = dStartTime[i] * 10000.0; /* in 10kHz */
    window[1] = window[0] + windowLen * 10000.0/actualSamplingRate;

    err = fread(&nextTimestamp, sizeof(uint32_t), 1, eegfile);
    if (!err) mxErrMsgTxt("fread error");
    while ((nextTimestamp+endRecord) < window[0]) {
      err = fread(recordData, sizeof(int16_t), nsamplesperbuf, eegfile);
      if (!err) mxErrMsgTxt("fread error");
      err = fread(&nextTimestamp, sizeof(uint32_t), 1, eegfile);
      if (!err) mxErrMsgTxt("fread error");
    }

    if (nextTimestamp > window[1])
      continue; /* overshot */

    err = fread(recordData, sizeof(int16_t), nsamplesperbuf, eegfile);
    if (!err) mxErrMsgTxt("fread error");

    odx = round((window[0] - nextTimestamp) * actualSamplingRate/10000.0);
    dStartActual[i] = (nextTimestamp + odx*10000.0/actualSamplingRate) / 10000.0;

    dataStartTimestamp = nextTimestamp;
    for (j = 0; j < windowLen; j++) {
      dptr[j] = recordData[odx++];

      if (odx >= nsamplesperbuf) {
        err = fread(&nextTimestamp, sizeof(uint32_t), 1, eegfile);
        if (!err) mxErrMsgTxt("fread error");
        if (nextTimestamp > dataStartTimestamp + recordLength) {
          fprintf(stderr,"skip length: %d",(dataStartTimestamp + recordLength - nextTimestamp));
          if (nextTimestamp < window[1]) {
            odx = 0;
            j = j + (dataStartTimestamp + recordLength - nextTimestamp) * actualSamplingRate / 10000.0;
          }
        }
        else {
          err = fread(recordData, sizeof(int16_t), nsamplesperbuf, eegfile);
          if (!err) mxErrMsgTxt("fread error");
          odx = 0;
        }
        dataStartTimestamp = nextTimestamp;
      }
    }
  }
  fprintf(stderr,"\n");

  mxFree(data);
  fclose(eegfile);
}


/*
double findRealTimestamp () {
  static u32 lastTimestamps[2];
  int d1, d2;
  double offset;
  fpos_t filePos;
  size_t err;
  double recordTimestamps;
    
  recordTimestamps = (double) record->numSamples * 10000.0 / actualSamplingRate; 

  fgetpos(eegfile, &filePos);
  err = fread(&lastTimestamps[1], sizeof(u32), 1, eegfile);
  if (!err) mxErrMsgTxt("fread error");
  fseek(eegfile, sizeof(int) + sizeof(double) + record->numSamples*sizeof(int16_t), SEEK_CUR);
  err = fread(&lastTimestamps[0], sizeof(u32), 1, eegfile);
  if (!err) mxErrMsgTxt("fread error");
  fsetpos(eegfile, &filePos);

  offset = 0.0;
  d1 = round((lastTimestamps[0] - 2*recordTimestamps) - (record->timestamp + offset));
  d2 = round((lastTimestamps[1] - recordTimestamps) - (record->timestamp + offset));
  if ((d1 != 0) || (d2 != 0)) {
    offset = 1.0/3;
    d1 = round((lastTimestamps[0] - 2*recordTimestamps) - (record->timestamp + offset));
    d2 = round((lastTimestamps[1] - recordTimestamps) - (record->timestamp + offset));
    if ((d1 != 0) || (d2 != 0)) {
      offset = -1.0/3;
      d1 = round((lastTimestamps[0] - 2*recordTimestamps) - (record->timestamp + offset));
      d2 = round((lastTimestamps[1] - recordTimestamps) - (record->timestamp + offset));
    }
  }

  return offset;
}
*/
