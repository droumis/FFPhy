/*
 * =============================================================
 *
 *  epoch_extract_c.c
 *
 *  Caleb Kemere 2008-09-03
 *
 *  [epochs, samplingrate, numsamples]
 *    = epoch_extract(filename, verboseFlag);
 *
 * =============================================================
 */

#include <stdio.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

#define MAX_FILENAME 200
#define HEADER_STRLEN 50
typedef unsigned int u32;

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  FILE *eegfile;
  int verboseFlag = 0;
  double samplingrate, samplingrate_tmp;
  double actualSamplingrate;
  unsigned int numsamples, numsamples_tmp;
  u32 timestamp, timestamp_last;
  double tstamp_step;
  int f_step, c_step;
  double *tstamp_skip_ptr;
  short *data;
  short *matlabData;
  char *header;
  struct stat filestat;
  char filename[MAX_FILENAME];
  int filenameLength=0;
  char headerString[500];
  char tmpString[HEADER_STRLEN];
  size_t dataStart;
  int recordLength;
  long int recordsPerFile;
  long int i;
  int temp;
  char warnMsg[1000];
  
  /* PROCESS INPUT ARGUMENTS */
  /* ----------------------------------------------------------*/
  /*  Check for proper number of arguments. */
  if (nrhs < 1) 
    mexErrMsgTxt("Filename required.");
  
  /* Get the verbosity verboseFlag. */
  if (nrhs == 2)
    verboseFlag = mxGetScalar(prhs[1]);

  /* Check to make sure file exists and is openable. */
  if (mxIsChar(prhs[0])) {
    filenameLength = mxGetN(prhs[0]);
    if (mxGetM(prhs[0]) > filenameLength)
      filenameLength = mxGetM(prhs[0]);
    if (filenameLength > MAX_FILENAME)
      mexErrMsgTxt("Filename/path too long.");

    mxGetString(prhs[0], filename, filenameLength+1);
    eegfile = fopen(filename, "r");
    if (!eegfile) {
      printf("Filename %s\n", filename);
      printf("Return value %d\n", eegfile);
      mexErrMsgTxt("Could not open file.");
    }
  }

  /* ----------------------------------------------------------*/

  /* Get file size for time prediction. */
  stat(filename,&filestat);

  /* Read in header.... */
  fgets(tmpString,HEADER_STRLEN,eegfile);
  if (strcmp("%%BEGINHEADER\n",tmpString) != 0)
    mexErrMsgTxt("File has unknown header.");
  fgets(tmpString,HEADER_STRLEN,eegfile);
  headerString[0] = (char) 0;
  while (strcmp("%%ENDHEADER\n",tmpString) != 0) {
    strcat(headerString,tmpString);
    strcat(headerString,"\n");
    fgets(tmpString,HEADER_STRLEN,eegfile);
  }

  /* ----------------------------------------------------------*/

  /* PROCESS OUTPUT ARGUMENTS */
  if (nlhs != 3) 
    mexErrMsgTxt("Three outputs required.");

  /* ----------------------------------------------------------*/

  dataStart = ftell(eegfile);

  /* Read in first data record. */
  fread(&timestamp_last, sizeof(u32), 1, eegfile); 
  fread(&numsamples, sizeof(int), 1, eegfile);
  fread(&samplingrate, sizeof(double), 1, eegfile); 

  actualSamplingrate = 30000.0;
  i = 0;
  while (floor(actualSamplingrate/(i+1)) >= samplingrate)
    i = i + 1;
  actualSamplingrate = actualSamplingrate / i;

  plhs[1] = mxCreateDoubleScalar(actualSamplingrate); /* output 4: sampling rate */
  plhs[2] = mxCreateDoubleScalar(numsamples); /* output 3: numsamples */

  recordLength = sizeof(u32) + sizeof(int) + sizeof(double) + \
                 sizeof(short)*numsamples;
  recordsPerFile = (filestat.st_size - dataStart)/recordLength;
  if ( fmod(((double) filestat.st_size - dataStart), (double)recordLength) != 0) {
    sprintf(warnMsg,"File length does not match predicted: %d, %ld, %d %ld\n", dataStart, (long int) filestat.st_size, recordLength, recordsPerFile);
    mexWarnMsgTxt(warnMsg);
  }

  plhs[0] = mxCreateDoubleMatrix(2,recordsPerFile,mxREAL); /* epochs */
  tstamp_skip_ptr = mxGetPr(plhs[0]);
  *tstamp_skip_ptr++ = timestamp_last;
  tstamp_step = (numsamples * 10000) / actualSamplingrate;
  f_step = (int) floor(tstamp_step);
  c_step = (int) ceil(tstamp_step);

  data = malloc(sizeof(short)*numsamples);
  fread(data, sizeof(short), numsamples, eegfile);

  if (verboseFlag > 2)
    printf("Records per file: %d\n", recordsPerFile);

  for (i = 2; i <= recordsPerFile; i++) {

    temp = fread(&timestamp, sizeof(u32), 1, eegfile); 
    if (!temp) mxErrMsgTxt("fread error");
    if ( ((timestamp - timestamp_last) != f_step) &
       ((timestamp - timestamp_last) != c_step) ) {
      *tstamp_skip_ptr++ = timestamp_last;
      *tstamp_skip_ptr++ = timestamp;
    }
    timestamp_last = timestamp;

    temp = fread(&numsamples_tmp, sizeof(int), 1, eegfile);
    if (!temp) mxErrMsgTxt("fread error");
    if (numsamples_tmp != numsamples) {
      printf("Old: %X, new: %X\n", numsamples, numsamples_tmp);
      mxErrMsgTxt("Numsamples changed!!");
    }

    temp = fread(&samplingrate_tmp, sizeof(double), 1, eegfile); 
    if (!temp) mxErrMsgTxt("fread error");
    if (samplingrate_tmp != samplingrate) {
      sprintf(warnMsg,"Sampling rate changed: %e -> %e (timestamp %ld)\n",samplingrate, samplingrate_tmp, timestamp);
      mxErrMsgTxt(warnMsg);
    }

    temp = fread(data, sizeof(short), numsamples, eegfile);
    if (!temp) mxErrMsgTxt("fread error");

    if ((verboseFlag > 1) && (fmod(i,10000) == 0))
      printf(".",i);
  }
  *tstamp_skip_ptr = timestamp;
  
  fclose(eegfile);
}


