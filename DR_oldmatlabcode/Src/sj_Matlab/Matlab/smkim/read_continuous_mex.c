/*=============================================================================
%READ_CONTINUOUS_MEX Read single-channel continuous data from an NSpike *.eeg file
%
%   CONTINUOUS = READ_CONTINUOUS_MEX(FILENAME) opens the NSpike .eeg file
%   specified by FILENAME and reads all records in the file. CONTINUOUS is a
%   struct array in which each element corresponds to a contiguous (gapless)
%   chunk of continuous data. If there are N gaps in data acquisition, then the
%   length of CONTINUOUS will be (N+1). If the sampling rate is not the same
%   across all records, an error will be raised.
%
%   CONTINUOUS = READ_CONTINUOUS_MEX(FILENAME,START_TS,END_TS) opens the NSpike
%   .eeg file specified by FILENAME and reads records that cover the time
%   interval specified by the timestamp values START_TS, END_TS. START_TS and
%   END_TS must be uint32 scalars. CONTINUOUS is a scalar struct when START_TS
%   and END_TS are specified. The records must contiguously cover the time
%   interval between START_TS and END_TS; if there is a gap in the data, an
%   error will be raised.
%
%   The return struct CONTINUOUS has the following fields:
%     Fs: the nominal sampling rate (double scalar)
%     samples: a single-precision floating-point vector of voltages, expressed
%       in microvolts
%     timestamp: the uint32 timestamps of the samples
%
%   This MEX function depends on standard C libraries. It is meant to be called
%   through the wrapper M-function READ_CONTINUOUS.
%
%Written by SMK, 2009 May 29.
=============================================================================*/

#define _FILE_OFFSET_BITS 64
#define fopen fopen64
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "mex.h"

typedef struct continuous_record_t {
	uint32_t  timestamp;
	int			  num_samples;
  double    sampling_rate;
  int16_t*  samples;
} continuous_record_t;

/* NSpike timestamp clock ticks at 10kHz */
const double TS_PER_SEC = 10000.0;
/* criterion for detecting timing anomalies of the samples */
const int32_t TS_PRECISION = 2; /* 2 timestamp units = 0.2 ms */
/* criterion for comparing equality of sampling rates (because sampling_rate is
a double, the == binop can give misleading result */
const double SAMP_RATE_PRECISION = 0.000001;
/* maximum number of samples that we can store in a single chunk of contiguous
records. 5e7 samples is a generous allocation (over 9 hours at 1500 samples per
second). */
const int32_t MAX_NUM_SAMPLES = 5e7;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  char tmp_string[1000];
  mwSize alloc_sz;
  unsigned long i;
  unsigned long record_count = 0;
  continuous_record_t current_record, next_record;
  /* Very important: initialize the samples pointer to NULL! */
  current_record.samples = NULL;
  next_record.samples = NULL;
  
  /* Process arguments */
  if ( !((nrhs == 3) || (nrhs == 1)) || (nlhs > 1) ) {
    mexErrMsgTxt(
        "Usage: continuous = read_continuous(filename)\n"
        "   or  continuous = read_continuous(filename, tstart, tend)\n"
        "where tstart and tend are uint32 timestamps\n"); 
  }
  char* filename;
  FILE* file_handle;
  if (mxIsChar(prhs[0])) {
    alloc_sz = (mwSize) (1 + mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar));
    filename = (char*) mxMalloc(alloc_sz);
    mxGetString(prhs[0], filename, alloc_sz);
    if ((file_handle = fopen(filename, "r")) == NULL) {
      sprintf(tmp_string, "Could not open file %s\n", filename);
      mexErrMsgTxt(tmp_string);
    }
  } else {
    mexErrMsgTxt("First input argument must be a filename string\n");
  }
  uint32_t start_timestamp;
  uint32_t end_timestamp;
  if (nrhs == 1) {
    start_timestamp = 0;
    end_timestamp = UINT32_MAX;
  } else {
    /* START_TS and END_TS must be scalar uint32 */
    if ( !mxIsNumeric(prhs[1]) || !mxIsClass(prhs[1], "uint32") || 
        mxIsComplex(prhs[1]) || (mxGetNumberOfDimensions(prhs[1]) > 2) ||
        (mxGetM(prhs[1]) != 1) || (mxGetN(prhs[1]) != 1) ) {
      mexErrMsgTxt("Second argument must be real uint32 timestamp\n");
    } else {
      start_timestamp = (uint32_t) mxGetScalar(prhs[1]); 
    }
    if ( !mxIsNumeric(prhs[2]) || !mxIsClass(prhs[2], "uint32") || 
        mxIsComplex(prhs[2]) || (mxGetNumberOfDimensions(prhs[2]) > 2) ||
        (mxGetM(prhs[2]) != 1) || (mxGetN(prhs[2]) != 1) ) {
      mexErrMsgTxt("Third argument must be real uint32 timestamp\n");
    } else {
      end_timestamp = (uint32_t) mxGetScalar(prhs[2]); 
    }
    if (!(end_timestamp > start_timestamp)) {
      mexErrMsgTxt("END_TS must be later than START_TS\n");
    }
  } 
  /* Read newline-delimited chunks until the last line of the header is read */
  do {
    if (fgets(tmp_string, 1000, file_handle) == NULL) {
      sprintf(tmp_string, "Could not find end-of-header line in file %s\n", 
          filename);
      mexErrMsgTxt(tmp_string);
    }
  } while (strncmp(tmp_string, "%%ENDHEADER", 11) != 0);
  /* Read the first record */
  if (read_continuous_record(file_handle, &current_record) == 0) {
    sprintf(tmp_string, "Could not find valid continuous record at beginning "
        "of file %s\n", filename);
    mexErrMsgTxt(tmp_string);
  } else {
    record_count++;
  } 
  /* Grab the sampling rate of the first record and take this to be the nominal
  sampling rate. CAREFUL: WE ASSUME THAT THE SAMPLING RATE DOESN'T CHANGE! */ 
  double nominal_sampling_rate = current_record.sampling_rate;
  /* If time range is specified, make sure that the first record precedes the
  desired start time, and continue reading records until the samples fall
  within range */
  if (nrhs == 3) {
    if (current_record.timestamp >= start_timestamp) {
      sprintf(tmp_string, "First continuous record in file %s has timestamp "
          "%lu, which does not precede desired START_TS (%lu)\n", filename, 
          (unsigned long) current_record.timestamp, 
          (unsigned long) start_timestamp);
      mexErrMsgTxt(tmp_string);
    }
    while (current_record.timestamp + (uint32_t) ceil(TS_PER_SEC *
        current_record.num_samples / nominal_sampling_rate) <
        start_timestamp) {
      if (read_continuous_record(file_handle, &next_record) == 0) {
        sprintf(tmp_string, "Finished reading to end of file but did not "
            "reach desired START TS (%lu); last record has timestamp %lu\n", 
            (unsigned long) start_timestamp, 
            (unsigned long) current_record.timestamp);
        mexErrMsgTxt(tmp_string);
      } else {
        record_count++;
      }
      /* Check for anomalous fluctuation of sampling rate or out-or-order
      timestamps. It is very fucking important to use (double) FABS, not 
      (int) ABS */
      if (fabs(next_record.sampling_rate - nominal_sampling_rate) >
          SAMP_RATE_PRECISION) {
        sprintf(tmp_string, "Sampling rate of next record (%f) does not match "
            "nominal sampling rate (%f)", next_record.sampling_rate, 
            nominal_sampling_rate);
        mexErrMsgTxt(tmp_string);
      }
      if (next_record.timestamp <= current_record.timestamp) {
        sprintf(tmp_string, "Consecutive records have time stamps that are "
            "out of order: last %lu, next %lu", 
            (unsigned long) current_record.timestamp,
            (unsigned long) next_record.timestamp);
        mexErrMsgTxt(tmp_string);
      }
      /* Free memory allocated for current_record.samples */
      mxFree(current_record.samples);
      /* Copy next_record to current_record */
      current_record = next_record; 
      /* Now current_record.samples points to the same block of memory as
      next_record.samples, so we can unseat the next_record.samples pointer */
      next_record.samples = NULL;
    }
  } 

  /* Dynamically allocated arrays of data */
  uint32_t* num_samples;
  num_samples = mxMalloc((mwSize) sizeof(uint32_t));
  float** samples; /* jagged array of arrays */
  samples = mxMalloc((mwSize) sizeof(float*));
  uint32_t** timestamp; /* jagged array of arrays */
  timestamp = mxMalloc((mwSize) sizeof(uint32_t*));
  
  /* Assuming that subsequent records are about the same length as
  current_record, estimate the number of samples that will need to be read,
  with a 5% fudge factor to make sure that we allocate enough space. Allocate
  this number or MAX_NUM_SAMPLES, whichever is lesser. */
  alloc_sz = (mwSize) (1.05 * nominal_sampling_rate * (end_timestamp - 
      current_record.timestamp) / TS_PER_SEC * sizeof(double)); 
  if (alloc_sz > MAX_NUM_SAMPLES * sizeof(double)) {
    alloc_sz = (mwSize) (MAX_NUM_SAMPLES * sizeof(double));
  }
  /* Initialize the first chunk */
  num_samples[0] = 0;
  samples[0] = (float*) mxMalloc(alloc_sz);
  timestamp[0] = (uint32_t*) mxMalloc(alloc_sz);
  unsigned long int j = 1; /* length of the output struct array */

  /* Continue reading records until we reach END_TS or end of file */
  /* Remember that end_timestamp can equal UINT32_MAX! Be careful to avoid
  triggering an overflow wrap-around */
  while (current_record.timestamp - (uint32_t) 
      (0.50 + TS_PER_SEC / nominal_sampling_rate) <= end_timestamp) {
    /* Copy samples from current_record to output array, cast as double, and
    fill in uint32_t timestamps that are extrapolated from the start timestamp
    of the record given the sampling rate */
    for (i = 0; i < current_record.num_samples; i++) {
      samples[j-1][num_samples[j-1]] = (float) current_record.samples[i];
      /* the +0.50 is to force rounding to the nearest int */
      timestamp[j-1][num_samples[j-1]] = (uint32_t) (0.50 + 
          ((double) current_record.timestamp) + 
          ((double) i) * TS_PER_SEC / nominal_sampling_rate);
      num_samples[j-1]++;
    }

    if (read_continuous_record(file_handle, &next_record) == 0) {
      break;
    } else {
      record_count++;
      /* Check next record for anomalous fluctuation of sampling rate or
      inplausible timestamp. It is very fucking important to use (double) FABS,
      not (int) ABS */
      if (fabs(next_record.sampling_rate - nominal_sampling_rate) >
          SAMP_RATE_PRECISION) {
        sprintf(tmp_string, "Sampling rate of next record (%f) does not match "
            "nominal sampling rate (%f)", next_record.sampling_rate, 
            nominal_sampling_rate);
        mexErrMsgTxt(tmp_string);
      }
      if (next_record.timestamp <= current_record.timestamp) {
        sprintf(tmp_string, "Consecutive records have time stamps that are "
            "out of order: last %lu, next %lu", 
            (unsigned long) current_record.timestamp,
            (unsigned long) next_record.timestamp);
        mexErrMsgTxt(tmp_string);
      }

      uint32_t expected_next_timestamp = (uint32_t) (((double) 
          timestamp[j-1][0]) + ((double) num_samples[j-1]) * TS_PER_SEC / 
          nominal_sampling_rate + 0.50);
      int32_t time_gap = next_record.timestamp - expected_next_timestamp;  
      if (time_gap < -TS_PRECISION) {
        sprintf(tmp_string, "Next record is too early for sampling rate %f: " 
            "gap size is %lu, expected timetamp is %lu, recorded timestamp "
            "is %lu.\n", nominal_sampling_rate, (unsigned long) time_gap,
            (unsigned long) expected_next_timestamp, 
            (unsigned long) next_record.timestamp);
        mexErrMsgTxt(tmp_string);
      } else if (time_gap > TS_PRECISION) {
        if (nrhs == 3) {
          /* If START_TS and END_TS were specified, then any gap is a failure
          condition */
          sprintf(tmp_string, "Gap of %lu timestamp units detected: for "
              "sampling rate %f, expected timetamp is %lu, recorded timestamp "
              "is %lu.\n", (unsigned long) time_gap, nominal_sampling_rate, 
              (unsigned long) expected_next_timestamp, 
              (unsigned long) next_record.timestamp);
          mexErrMsgTxt(tmp_string);
        } else {
          /* Otherwise we just notify the user that we found a gap and advance to
          the next chunk */
          sprintf(tmp_string, "Gap of %lu timestamp units (%f seconds) "
              "detected after record number %lu\n", (unsigned long) 
              time_gap, ((double) time_gap)/TS_PER_SEC, record_count);
          mexWarnMsgTxt(tmp_string);
          /* to save memory, shrink the dynamic arrays corresponding to last
          finished chunk */
          samples[j-1] = (float*) mxRealloc(samples[j-1], (mwSize) 
              (num_samples[j-1]*sizeof(float)));
          timestamp[j-1] = (uint32_t*) mxRealloc(timestamp[j-1], (mwSize) 
              (num_samples[j-1]*sizeof(uint32_t)));
       /* LOOK!!!! j is now incremented to the next chunk
           |
           |
           V       (don't make an off-by-one error!) */
          j++; 
          /* Reallocate to accomodate next chunk */
          num_samples = (uint32_t*) mxRealloc(num_samples, 
              (mwSize) (j*sizeof(uint32_t)));
          samples = (float**) mxRealloc(samples, (mwSize) 
              (j*sizeof(float*)));
          timestamp = (uint32_t**) mxRealloc(timestamp, (mwSize) 
              (j*sizeof(uint32_t*)));
          alloc_sz = (mwSize) (1.05 * nominal_sampling_rate * (end_timestamp - 
              next_record.timestamp) / TS_PER_SEC * sizeof(double)); 
          if (alloc_sz > MAX_NUM_SAMPLES * sizeof(double)) {
            alloc_sz = (mwSize) (MAX_NUM_SAMPLES * sizeof(double));
          }
          num_samples[j-1] = 0;
          samples[j-1] = (float*) mxMalloc(alloc_sz);
          timestamp[j-1] = (uint32_t*) mxMalloc(alloc_sz);
        } 
      } 
      /* Free memory allocated for current_record, swap next_record to
      current_record, and unseat next_record->samples pointer */
      mxFree(current_record.samples);
      current_record = next_record;
      next_record.samples = NULL;
    }
  }

  /* If END_TS was specified, compare to the timestamp of the final sample that
  was read */
  if ((nrhs == 3) && (timestamp[j-1][num_samples[j-1]-1] <= end_timestamp)) {
    sprintf(tmp_string, "Failed to pass the desired END_TS (%lu); last "
        "sample was %lu\n", (unsigned long) end_timestamp, 
        (unsigned long) timestamp[j-1][num_samples[j-1]-1]);
    mexErrMsgTxt(tmp_string);
  }

  /* Allocate and assign output struct array. */
  mwSize dims[2] = {j, 1};
  const char* fieldnames[] = 
      {"Fs", "timestamp", "samples"};
  plhs[0] = mxCreateStructArray(2, dims, 3, fieldnames);
  if (plhs[0] == NULL) {
    mexErrMsgTxt("Unable to allocate output struct array\n");
  }
  for (i = 0; i < j; i++) {
    mxArray* out_timestamp = mxCreateNumericMatrix((mwSize) num_samples[i],
        1, mxUINT32_CLASS, mxREAL);
    mxArray* out_samples = mxCreateNumericMatrix((mwSize) num_samples[i],
        1, mxSINGLE_CLASS, mxREAL);

    memcpy((uint32_t*) mxGetData(out_timestamp), timestamp[i],
        num_samples[i]*sizeof(uint32_t));
    memcpy((float*) mxGetData(out_samples), samples[i], 
        num_samples[i]*sizeof(float));
    
    mxSetField(plhs[0], i, fieldnames[0], 
        mxCreateDoubleScalar(nominal_sampling_rate));
    mxSetField(plhs[0], i, fieldnames[1], out_timestamp);
    mxSetField(plhs[0], i, fieldnames[2], out_samples);
  }

  /* Clean up */
  if ((file_handle != NULL) && (fclose(file_handle) != 0)) {
    sprintf(tmp_string, "Could not close file %s\n", filename);
    mexErrMsgTxt(tmp_string);
  }
  mxFree(filename);
  if (current_record.samples != NULL) {
    mxFree(current_record.samples);
  }
  if (next_record.samples != NULL) {
    mxFree(next_record.samples);
  }
  mxFree(num_samples);
  for (i = 0; i < j; i++) {
    mxFree(samples[i]);
    mxFree(timestamp[i]);
  }
  mxFree(timestamp);
  mxFree(samples);
}

/* This function reads the next continuous record from the current position in
file_handle (and in doing so, advances the file read position) and seats the
pointer record_ptr to a struct that contains this data. If a valid record is
not found, the file read position is rewound to its original position, and 0 is
returned */
int read_continuous_record(FILE* file_handle, continuous_record_t* record_ptr) {
  char tmp_string[1000];
  off_t original_position;

  /* To prevent memory leak, we require that record_ptr->samples must be NULL
  (if it weren't NULL, then we might unseat it and orphan an allocated block of
  memory with no pointer) */
  if ((record_ptr != NULL) && (record_ptr->samples != NULL)) {
    sprintf(tmp_string,"Second input argument to read_continuous_record must "
        "point to a continuous_record_t with a NULL pointer in its samples "
        "field.\n");
    mexErrMsgTxt(tmp_string);
  }
  /* Grab the current file position */
  if ((original_position = ftello(file_handle)) == 1) {
    mexErrMsgTxt("Error in reading file position.\n");
  }
  /* The first three fields of the record are of fixed size and directly
  correspond to fields of the struct pointee of record_ptr */
  if ((fread(&(record_ptr->timestamp), sizeof(uint32_t), 1, file_handle) != 1) 
      || (fread(&(record_ptr->num_samples), sizeof(int), 1, file_handle) != 1) 
      || (fread(&(record_ptr->sampling_rate), sizeof(double), 1, file_handle) 
      != 1)) {
    if (fseeko(file_handle, original_position, SEEK_SET) != 0) {
      mexErrMsgTxt("Error in seeking to file position.\n");
    }
    return 0;
  }
  if (record_ptr->num_samples < 1) {
    sprintf(tmp_string,"Continuous record reports invalid number of samples %i\n",
        record_ptr->num_samples);
    mexErrMsgTxt(tmp_string);
  }
  if (record_ptr->sampling_rate < 0) {
    sprintf(tmp_string,"Continuous record reports invalid sampling rate %f\n",
        record_ptr->sampling_rate);
    mexErrMsgTxt(tmp_string);
  }
  /* The fourth field of the record is of variable size and contains the actual
  array data as packed binary 16-bit integers, not just a pointer; this is why 
  fread(record_ptr, sizeof(continuous_record_t), 1, file_handle) would be wrong
  */
  record_ptr->samples = (int16_t*) mxCalloc((mwSize) record_ptr->num_samples, 
      (mwSize) sizeof(int16_t));
  if (fread(record_ptr->samples, sizeof(int16_t), record_ptr->num_samples, 
      file_handle) != record_ptr->num_samples) { 
    /* Free the memory block that was just allocated */
    mxFree(record_ptr->samples);
    if (fseeko(file_handle, original_position, SEEK_SET) != 0) {
      mexErrMsgTxt("Error in seeking to file position.\n");
    }
    return 0;
  }     
  return 1;
}

