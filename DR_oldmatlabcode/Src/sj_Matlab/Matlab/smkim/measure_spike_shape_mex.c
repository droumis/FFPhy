/*=============================================================================
%MEASURE_SPIKE_SHAPE Measure shape of spike waveforms recorded in NSpike .tt format
%   SHAPE = MEASURE_SPIKE_SHAPE_MEX(X,T,I) finds the peaks, troughs, and 
%   zero-crossings of the thresold-triggered spike waveforms specified by 
%   X,T,I.
%   
%   X is a MxCxN array of single-precision floating-point values. C is the 
%   number of recording channels in the C-trode, N is the number of 
%   trigger events, and M is the number of samples acquired per window.
%   
%   T is a Mx1 vector of single-precision floating-point values, which 
%   gives the times of the M samples relative to the first sample in the 
%   window. T is assumed to be monotonically increasing. 
%   
%   I is an integer scalar between 1 and M, which indexes along T and along 
%   the first dimension of X. I specifies the time of the trigger within 
%   the M-sample acquisition window. I is expressed with MATLAB-style 
%   (i.e., the first element of T corresponds to I = 1).
%   
%   The output argument SHAPE is a struct whose fields are NxC arrays of 
%   single-precision floating-point values"
%
%   Troughs and peaks are found by parabolic interpolation, and
%   rising/falling/crossing points are found by linear interpolation. The
%   interpolation can be improved by pre-processing with cubic spline
%   interpolation, which is what is done in the m-file MEASURE_SPIKE_SHAPE.
%
%Written by SMK, 2009 August 26.
=============================================================================*/

#define _FILE_OFFSET_BITS 64
#define fopen fopen64
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include "mex.h"

/* this is the basic struct type for holding waveform shape parameters */
typedef struct waveform_shape_t {
  float   trough_amplitude;
  float   peak_amplitude;
  float   half_trough_time_before;
  float   trough_time;
  float   half_trough_time_after;
  float   zero_crossing_time;
  float   half_peak_time_before;
  float   peak_time;
  float   half_peak_time_after;
} waveform_shape_t;

void usage();
int is_valid_time_vector(float* time, unsigned int num_points,
    unsigned int trigger_idx);
int calculate_shape_parameters(float* samples, float* time, 
    unsigned int trigger_idx, unsigned int samples_per_window, 
    unsigned int num_channels, unsigned long event_idx, 
    waveform_shape_t* shape);
int parabolic_extremum_interpolation(float t1, float t2, float t3, float
    x1, float x2, float x3, float* t_extremum_ptr, float* x_extremum_ptr);
int linear_crossing_interpolation(float t1, float t2, float x1, float x2,
    float x_cross, float* t_cross_ptr);
void copy_to_output(float** output_data, waveform_shape_t* shape, 
    unsigned long num_events, unsigned int num_channels, 
    unsigned long event_idx);

/* platform-specific value of NaN as float; this value is used to indicate
that a shape parameter can not be computed because the waveform does not have
the expected shape */
float NAN_FLOAT; 

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  /* Obtain the platform-specific value of NaN as a single-precision float.
  mxGetNan() returns a double, which we downcast to a float(=mxSINGLE_CLASS).
  This downcast from (double) to (float) may yield unexpected behavior on some
  platforms! */
  NAN_FLOAT = (float) mxGetNaN();
  if (!isnanf(NAN_FLOAT)) {
    mexErrMsgTxt("Could not get float NaN value on this platform. "
        "This is a limitation of the MATLAB C API.\n");
  }

  char message[10000];

  /* Check input arguments */
  if ( !(nrhs == 3) || (nlhs > 1) ) {
    usage();
  }
  if ( !mxIsNumeric(prhs[0]) || (mxGetNumberOfDimensions(prhs[0]) != 3) || 
      !mxIsSingle(prhs[0]) || mxIsComplex(prhs[0]) || mxIsEmpty(prhs[0]) ) {
    mexErrMsgTxt("X must be a non-empty 3-dimensional array of real "
        "single-precision floating-point values\n");
  }
  const mwSize* x_dims = mxGetDimensions(prhs[0]);
  if ( !mxIsNumeric(prhs[1]) || (mxGetNumberOfDimensions(prhs[1]) != 2) || 
      !mxIsSingle(prhs[1]) || mxIsComplex(prhs[1]) || mxIsEmpty(prhs[1]) ||
      (mxGetM(prhs[1]) != x_dims[0]) ||
      (mxGetNumberOfElements(prhs[1]) != x_dims[0]) ) {
    /* ought to check for monotonic increasing values, but I'm lazy :( */
    mexErrMsgTxt("T must be a column vector of real single-precision "
        "floating-point values whose number of elements equals size(X,1)\n");
  }
  if ( !mxIsNumeric(prhs[2]) || (mxGetNumberOfDimensions(prhs[2]) != 2) || 
      mxIsComplex(prhs[2]) || (mxGetM(prhs[2]) != 1) || 
      (mxGetN(prhs[2]) != 1) || (mxGetScalar(prhs[2]) <= 1) || 
      (mxGetScalar(prhs[2]) >= x_dims[0]) ||
      (mxGetScalar(prhs[2]) != round(mxGetScalar(prhs[2]))) ) {
    mexErrMsgTxt(
        "I must be a real integer-valued scalar between 1 and size(X,1)\n");
  }

  /* Convert input from MATLAB to C data structures */
  unsigned int samples_per_window = (unsigned int) x_dims[0];
  unsigned int num_channels = (unsigned int) x_dims[1];
  unsigned long num_events = (unsigned long) x_dims[2];
  float* time = (float*) mxGetData(prhs[1]);
  /* subtract 1 to convert MATLAB-style index to C-style index */
  unsigned int trigger_idx = (unsigned int) (mxGetScalar(prhs[2]) - 1);
  /* check T against X and I */
  if ( !is_valid_time_vector(time, samples_per_window, trigger_idx) ) {
    mexErrMsgTxt("T must be a uniformly-spaced monotonically-increasing\n");
  }

  /* Allocate output MATLAB struct */
  const mwSize out_dims[2] = {1, 1};
  /* CAREFUL! Any changes to field_names must be accompanied by edits to the
  function copy_to_output */
  const char* field_names[] = {
      "trough_amplitude",
      "peak_amplitude",
      "half_trough_time_before",
      "trough_time",
      "half_trough_time_after",
      "zero_crossing_time",
      "half_peak_time_before",
      "peak_time",
      "half_peak_time_after" };
  /* count the number of elements in field_names */
  int num_fields = sizeof(field_names) / sizeof(*field_names);
  plhs[0] = mxCreateStructArray(2, out_dims, num_fields, field_names);
  if (plhs[0] == NULL) {
    mexErrMsgTxt("Unable to allocate output struct array\n");
  }
  /* array of pointers for populating the fields of plhs[0] */
  mxArray** field_values;
  field_values = mxMalloc((mwSize) num_fields*sizeof(mxArray*));
  /* array of pointers to the data parts of the elements of field_values */
  float** output_data;
  output_data = mxMalloc((mwSize) num_fields*sizeof(float*));
  unsigned int field_idx; /* iterator over struct fields */
  for (field_idx = 0; field_idx < num_fields; field_idx++) {
    /* conveniently, all fields are same size and data type */
    field_values[field_idx] = mxCreateNumericMatrix((mwSize) num_events, 
        (mwSize) num_channels, mxSINGLE_CLASS, mxREAL);
    if (field_values[field_idx] == NULL) {
      sprintf(message, "Unable to allocate mxArray for %uth field of the "
          "output MATLAB struct (%s)\n", field_idx, field_names[field_idx]);
      mexErrMsgTxt(message);
    }
    output_data[field_idx] = (float*) mxGetData(field_values[field_idx]);
    mxSetFieldByNumber(plhs[0], 0, (mwIndex) field_idx, 
        field_values[field_idx]); 
  }

  /* Construct an array waveform_shape_t[num_channels] to store waveform
  parameters recorded on all channels for one threshold-trigger event */
  waveform_shape_t* shape;
  shape = (waveform_shape_t*) mxMalloc((mwSize) num_channels *
      sizeof(waveform_shape_t));

  /* Process the threshold-trigger events one by one */
  unsigned long event_idx;
  for (event_idx = 0; event_idx < num_events; event_idx++) {
    /* Calculate shape parameters for the ith threshold-trigger event */
    if ( !calculate_shape_parameters((float*) mxGetData(prhs[0]), time, 
        trigger_idx, samples_per_window, num_channels, event_idx, shape) ) {
      sprintf(message, "Error processing trigger event number %lu\n", 
          event_idx);
      mexErrMsgTxt(message);
    }
    /* Copy results to the output MATLAB struct */
    copy_to_output(output_data, shape, num_events, num_channels, event_idx);
  }

  /* Cleanup */
  mxFree(shape);
  /* it's okay to destroy these arrays of pointers; plhs[0] retains pointers to
  the pointees and will automatically cleanup upon exit */
  mxFree(field_values); 
  mxFree(output_data);

} 


/******************************************************************************
is_valid_time_vector(time, samples_per_window, trigger_idx)

Check if the vector time is monotonically-increasing.
*******************************************************************************/
int is_valid_time_vector(float* time, unsigned int samples_per_window,
    unsigned int trigger_idx) {

  if (samples_per_window <= trigger_idx) {
    mexPrintf("number of elements in time vector must exceed trigger_idx\n");
    return 0;
  }
  unsigned int i;
  float time_diff = time[1] - time[0];
  if (time_diff <= 0) {
    mexPrintf("time vector must be monotonic increasing\n");
    return 0;
  }
  const float tol = 1e-8; /* 10 nanoseconds */
  /* this loop start from index 1! */
  for (i = 1; i < samples_per_window; i++) {
    if (fabs((time[i] - time[i-1]) - time_diff) >= tol) {
      mexPrintf("first-order differences between elements of time vector are "
          "discrepant: difference is %g, tolerance is %g\n", 
          time[i] - time[i-1], tol);
      return 0;
    }
  }
  return 1;

}

/******************************************************************************
calculate_shape_parameters(samples, time, trigger_idx, samples_per_window,
num_channels, event_idx, shape)

Read (SAMPLES_PER_WINDOW x NUM_CHANNELS) floats from contiguous memory in
SAMPLES, corresponding to the EVENT_INDEXth threshold-trigger event, and
calculate waveform parameters (using TIME and TRIGGER_IDX) to pass to SHAPE.

NSpike threshold trigger rules: if *any* one of the channels goes more negative
than its respective channel-specific threshold, then acquisition is triggered.
There is no guarantee that *every* channel in the group will have a discernible
spike waveform, but there is a guarantee that *at least one* channel in the
group will exceed its threshold at trigger_idx.

I designed this algorithm to be robust to noisy channels. Noisy channels are
ignored and only minimal dependence between channels is assumed, so that a
noisy channels do not interfere with the estimation of waveform shape on the
other channels.

Step 1: If a channel has a non-negative sample value at trigger_idx, then it is
probably noisy. All waveform shape parameters are assigned NaN.

Step 2: For all channels that are negative at the trigger time, find the first
sample greater than zero *after* the trigger. shape.zero_crossing_time is found
by linear interpolation between this first positive sample and the immediately
preceding (negative or zero-valued) sample. If no post-trigger positive sample
can be found, then shape.zero_crossing_time is assigned NaN.

Step 3: For all channels that are negative at trigger and have a defined
shape.zero_crossing_time, search for the minimum sample between these two
timepoints. If this minimum is found after the trigger, then it is the trough
and we can calculate shape.trough_time and shape.trough_amplitude by parabolic
interpolation over three samples: the minimum, the sample before, and the
sample after. However, if the minimum is found *at* the trigger, then we must
search earlier to find the real minimum. (This can happen when the trough
arrives at slightly different times on different channels; papers by the
Buzsaki group show examples recorded with multisite silicon probes.) To bound
this pre-trigger search, search backwards to find the last positive-valued
sample *before* the trigger. If this pre-trigger positive sample exists, then
we find the minimum (trough) between this timepoint and the first
positive-valued sample after the trigger. If no such pre-trigger positive
sample can be found, then we look for the minimum between the beginning of the
window to the first positive-valued sample after the trigger. If the minimum is
found at the beginning of the window, then shape.trough_amplitude and
shape.trough_time are NaN. Note that this second search before the trigger can
introduce problems of its own; for example, when the spike is the second of a
floatt, then the preceding spike trough can be misidentified as the trough of
this spike. Hence we only resort to this pre-trigger search for a minimum when
a plausible trough can not be found after the trigger.

Step 4: For all channels that have a defined shape.trough_time and
shape.trough_amplitude, search forward and backward from the minimum to find
the first/last samples which fall below 0.5*shape.trough_amplitude (be careful
with sign). Calculate shape.half_trough_time_before and
shape.half_trough_time_after by linear interpolation. The half-amplitude point
before the trough may not exist if the initial negativity of the waveform is
truncated at the beginning of the window.

Step 5: For all channels that have a defined shape.zero_crossing_time, look
ahead to find the first negative-valued sample after the zero-crossing, which
helps to establish a time bound on the peak time. Again, we need to be careful.
Sometimes the second zero-crossing (from positive to negative) is just a
spurious noisy wiggle in the neighborhood of shape.zero_crossing_time,
especially if the spike has a low amplitude on this channel; we don't want to
misidentify this wiggle as the peak on that channel. Or a negative-valued
sample may not exist within the window after the zero-crossing, because the
spike has a broad positive peak which is truncated before it returns to zero.
For each channel which has a valid shape.zero_crossing_time, we find the
earliest return to negativity that occurs after shape.zero_crossing_time. We
then pick among these channels the *latest* of these second-negativity
crossings. (This is sort of a maximin approach). Now search for the maximum
between the post-trigger zero-crossing and this post-peak bound (or the end of
the window, if the bound is undefined). If this maximum is found at the end of
the window, then shape.peak_time and shape.peak_amplitude are NaN; otherwise,
calculate the peak with parabolic interpolation.

Step 6: For all channels that have a defined shape.peak_time and
shape.peak_amplitude, search forward and backward from the maximum to find the
first/last samples which fall below 0.5*shape.peak_amplitude, and calculate
shape.half_peak_time_before and shape.half_peak_time_after by linear
interpolation. The half-amplitude point after the peak may not exist if the
positive phase of the waveform is truncated at the end of the window.
*******************************************************************************/
int calculate_shape_parameters(float* samples, float* time, 
    unsigned int trigger_idx, unsigned int samples_per_window, 
    unsigned int num_channels, unsigned long event_idx, 
    waveform_shape_t* shape) {

  unsigned int c; /* channel index */
  unsigned int i; /* index for sample within window */
  unsigned long offset[num_channels];
  for (c = 0; c < num_channels; c++) {
    /* pointers to the correct location in the samples[] array */
    offset[c] = event_idx*num_channels*samples_per_window + 
        c*samples_per_window;
    /* fill in default NaN values */
    shape[c].trough_amplitude = NAN_FLOAT;
    shape[c].peak_amplitude = NAN_FLOAT;
    shape[c].half_trough_time_before = NAN_FLOAT;
    shape[c].trough_time = NAN_FLOAT;
    shape[c].half_trough_time_after = NAN_FLOAT;
    shape[c].zero_crossing_time = NAN_FLOAT;
    shape[c].half_peak_time_before = NAN_FLOAT;
    shape[c].peak_time = NAN_FLOAT;
    shape[c].half_peak_time_after = NAN_FLOAT;
  }

  unsigned int first_negative_after_trigger_idx[num_channels];
  unsigned int first_positive_after_trough_idx[num_channels];
  unsigned int minimum_idx[num_channels];
  unsigned int half_trough_before_idx[num_channels];
  unsigned int half_trough_after_idx[num_channels];
  unsigned int first_negative_after_zero_crossing_idx[num_channels];
  unsigned int maximum_idx[num_channels];
  unsigned int half_peak_before_idx[num_channels];
  unsigned int half_peak_after_idx[num_channels];
  float current_val;
  float best_val;
  /* unlike the other *idx variables, this one is not an array */
  unsigned int end_search_for_peak_idx;

  /* find the post-trough zero-crossing */
  for (c = 0; c < num_channels; c++) {
    /* check that the waveform is negative at trigger_idx; if it isn't, then
    scan forwards until a negativity is reached; otherwise, we start at
    trigger_idx */
    first_negative_after_trigger_idx[c] = samples_per_window - 1;
    for (i = trigger_idx; i < samples_per_window; i++) {
      if (samples[offset[c] + i] < 0) {
        first_negative_after_trigger_idx[c] = i;
        break;
      }
    }
    if (first_negative_after_trigger_idx[c] == samples_per_window - 1) {
      /* this is not a valid spike */
      continue;
    }
    /* find the zero-crossing after the trough and before the peak */
    first_positive_after_trough_idx[c] = first_negative_after_trigger_idx[c];
    for (i = first_negative_after_trigger_idx[c]; i < samples_per_window; 
        i++) {
      if (samples[offset[c] + i] > 0) {
        first_positive_after_trough_idx[c] = i;
        break;
      }
    }    
    if (first_positive_after_trough_idx[c] == 
        first_negative_after_trigger_idx[c]) {
      /* no positive sample found after trigger; zero-crossing is undefined, as
      are other waveform features */
      continue;
    } else if ( 
        (samples[offset[c] + first_positive_after_trough_idx[c] - 1] > 0) ||
        (samples[offset[c] + first_positive_after_trough_idx[c]] <= 0) ) {
      mexPrintf("expected negative-to-posibive zero-crossing, but instead "
          "consecutive samples share the same sign\n");
      return 0;
    } else if ( !linear_crossing_interpolation( 
        time[first_positive_after_trough_idx[c] - 1], 
        time[first_positive_after_trough_idx[c]], 
        samples[offset[c] + first_positive_after_trough_idx[c] - 1],
        samples[offset[c] + first_positive_after_trough_idx[c]], 
        0, &(shape[c].zero_crossing_time)) ) {
      mexPrintf("error while computing zero-crossing\n");
      return 0;
    }
  }

  /* for each channel, find the earliest negative-valued sample that comes
  after the zero-crossing */
  for (c = 0; c < num_channels; c++) {
    if (!isnanf(shape[c].zero_crossing_time)) {
      first_negative_after_zero_crossing_idx[c] = samples_per_window - 1;
      for (i = first_positive_after_trough_idx[c]; i < samples_per_window; 
          i++) {
        if (samples[offset[c] + i] < 0) {
          first_negative_after_zero_crossing_idx[c] = i;
          break;
        }
      }
    } else {
      /* if waveform does not have a valid zero-crossing after the trough, set
      to an impossible sentinel value that is guaranteed to be smaller than any
      valid value */
      first_negative_after_zero_crossing_idx[c] = 0;
    }
  }

  /* find the trough and its half-width */
  for (c = 0; c < num_channels; c++) {
    if (!isnanf(shape[c].zero_crossing_time)) {
      minimum_idx[c] = first_negative_after_trigger_idx[c];
      best_val = samples[offset[c] + minimum_idx[c]];
      for (i = first_negative_after_trigger_idx[c]; 
          i < first_positive_after_trough_idx[c]; i++) {
        current_val = samples[offset[c] + i];
        /* remember that the trough is negative */
        if (current_val < best_val) {
          best_val = current_val;
          minimum_idx[c] = i;
        }
      } 
      if (minimum_idx[c] == first_negative_after_trigger_idx[c]) {
        /* redo the search: look before trigger_idx in case the trigger came
        after the trough on this channel */
        minimum_idx[c] = 0;
        best_val = samples[offset[c] + minimum_idx[c]];
        for (i = first_positive_after_trough_idx[c] - 1; (i >= 0) && 
            (i < samples_per_window); i--) {
          current_val = samples[offset[c] + i];
          if (current_val > 0) {
            /* if we enter a preceding positive phase, terminate search and do
            not assign new value to minimum_idx[c] */
            break;
          } else if (current_val < best_val) {
            best_val = current_val;
            minimum_idx[c] = i;
          }
        }
      } 
      if (minimum_idx[c] == 0) {
        /* couldn't find the trough; give up */
        continue;
      } else if ( !parabolic_extremum_interpolation(
          time[minimum_idx[c] - 1], time[minimum_idx[c]], 
          time[minimum_idx[c] + 1], samples[offset[c] + minimum_idx[c] - 1], 
          samples[offset[c] + minimum_idx[c]], 
          samples[offset[c] + minimum_idx[c] + 1], &(shape[c].trough_time), 
          &(shape[c].trough_amplitude)) ) {
        mexPrintf("error while computing trough\n");
        return 0;
      }
    }
    if (!isnanf(shape[c].trough_time)) {
      /* find the falling half-width of the trough */
      half_trough_after_idx[c] = minimum_idx[c];
      for (i = minimum_idx[c]; 
          i <= first_positive_after_trough_idx[c]; i++) {
        /* be careful with sign! trough is negative */
        if (samples[offset[c] + i] > 0.5*shape[c].trough_amplitude) {
          half_trough_after_idx[c] = i;
          break;
        }
      }
      if (half_trough_after_idx[c] == minimum_idx[c]) {
        /* if the trough is so narrow that it only appears on one sample, then we
        can't measure it's width (and it's not a real spike anyway) so just
        give up */
        continue;
      } else if ( !linear_crossing_interpolation(
          time[half_trough_after_idx[c] - 1],
          time[half_trough_after_idx[c]], 
          samples[offset[c] + half_trough_after_idx[c] - 1],
          samples[offset[c] + half_trough_after_idx[c]], 
          0.5*shape[c].trough_amplitude, 
          &(shape[c].half_trough_time_after)) ) {
        mexPrintf("error while computing half-trough fall\n");
        return 0;
      }
      /* find the rising half-width of the trough */
      half_trough_before_idx[c] = minimum_idx[c];  
      for (i = minimum_idx[c]; (i >= 0) && (i < samples_per_window); i--) {
        /* be careful with sign! trough is negative */
        if (samples[offset[c] + i] > 0.5*shape[c].trough_amplitude) {
          half_trough_before_idx[c] = i;
          break;
        }
      }
      if (half_trough_before_idx[c] == minimum_idx[c]) {
        /* the rising half-width may be truncated at the beginning of the
        window if the trough is wide */
        shape[c].half_trough_time_before = NAN_FLOAT;
      } else if ( !linear_crossing_interpolation(
          time[half_trough_before_idx[c]],
          time[half_trough_before_idx[c] + 1], 
          samples[offset[c] + half_trough_before_idx[c]],
          samples[offset[c] + half_trough_before_idx[c] + 1], 
          0.5*shape[c].trough_amplitude, 
          &(shape[c].half_trough_time_before)) ) {
        mexPrintf("error while computing half-trough rise\n");
        return 0;
      }
    }
  }

  /* among those channels which have well-defined zero-crossings, find the
  latest occurence of the first negative-going zero-crossing after the
  post-trough zero-crossing. if none of the channels has a valid
  zero_crossing_time, then end_search_for_peak_idx will retain its sentinel
  initialization value of zero */
  end_search_for_peak_idx = 0;
  for (c = 0; c < num_channels; c++) {
    if (first_negative_after_zero_crossing_idx[c] > end_search_for_peak_idx) {
      end_search_for_peak_idx = first_negative_after_zero_crossing_idx[c];
    }
  }

  /* find the peak and its halfwidth */
  for (c = 0; c < num_channels; c++) {
    /* Check that the waveform has a defined zero-crossing after trough */
    if (isnanf(shape[c].zero_crossing_time) || 
        (end_search_for_peak_idx == 0)) { 
      continue;
    }
    maximum_idx[c] = end_search_for_peak_idx;
    best_val = samples[offset[c] + maximum_idx[c]];
    for (i = first_positive_after_trough_idx[c]; 
        i <= end_search_for_peak_idx; i++) {
      current_val = samples[offset[c] + i];
      if (current_val > best_val) {
        best_val = current_val;
        maximum_idx[c] = i;
      }
    } 
    if (maximum_idx[c] == end_search_for_peak_idx) {
      /* couldn't find the peak; give up */
      continue;
    } else if ( !parabolic_extremum_interpolation(
        time[maximum_idx[c] - 1], time[maximum_idx[c]], 
        time[maximum_idx[c] + 1], samples[offset[c] + maximum_idx[c] - 1], 
        samples[offset[c] + maximum_idx[c]], 
        samples[offset[c] + maximum_idx[c] + 1], &(shape[c].peak_time), 
        &(shape[c].peak_amplitude)) ) {
      mexPrintf("error while computing trough\n");
      return 0;
    }
    if (!isnanf(shape[c].peak_time)) {
      /* find the rising half-width of the peak */
      half_peak_before_idx[c] = maximum_idx[c];  
      for (i = maximum_idx[c]; i >= first_positive_after_trough_idx[c] - 1; 
          i--) {
        /* be careful with sign! peak is positive */
        if (samples[offset[c] + i] < 0.5*shape[c].peak_amplitude) {
          half_peak_before_idx[c] = i;
          break;
        }
      }
      if (half_peak_before_idx[c] == maximum_idx[c]) {
        /* if the peak is so narrow that it only appears on one sample, then we
        can't measure it's width (and it's not a real spike anyway) so just
        give up */
      } else if ( !linear_crossing_interpolation(
          time[half_peak_before_idx[c]],
          time[half_peak_before_idx[c] + 1], 
          samples[offset[c] + half_peak_before_idx[c]],
          samples[offset[c] + half_peak_before_idx[c] + 1], 
          0.5*shape[c].peak_amplitude, 
          &(shape[c].half_peak_time_before)) ) {
        mexPrintf("error while computing half-peak rise\n");
        return 0;
      }
      /* find the falling half-width of the peak */
      half_peak_after_idx[c] = maximum_idx[c];
      for (i = maximum_idx[c]; i < samples_per_window; i++) {
        /* be careful with sign! peak is positive */
        if (samples[offset[c] + i] < 0.5*shape[c].peak_amplitude) {
          half_peak_after_idx[c] = i;
          break;
        }
      }
      if (half_peak_after_idx[c] == maximum_idx[c]) {
        /* the falling half-width of the peak may be truncated at the end of
        the window if the peak is broad */
        shape[c].half_peak_time_after = NAN_FLOAT;
      } else if ( !linear_crossing_interpolation(
          time[half_peak_after_idx[c] - 1],
          time[half_peak_after_idx[c]], 
          samples[offset[c] + half_peak_after_idx[c] - 1],
          samples[offset[c] + half_peak_after_idx[c]], 
          0.5*shape[c].peak_amplitude, 
          &(shape[c].half_peak_time_after)) ) {
        mexPrintf("error while computing half-peak fall\n");
        return 0;
      }
    }
  }

  /* if we made it this far, then return success */
  return 1;

}


/******************************************************************************
parabolic_extremum_interpolation(t1, t2, t3, x1, x2, x3, t_out, x_out)

Given three points (t1,x1), (t2,x2), (t3,x3), such that t1 < t2 < t3 and
either ((x1 <= x2) && (x2 >= x3)) or ((x1 >= x2) && (x2 <= x3)), fit a parabola
to the three points and find the vertex.
*******************************************************************************/
int parabolic_extremum_interpolation(float t1, float t2, float t3, float x1, 
    float x2, float x3, float* t_extremum_ptr, float* x_extremum_ptr) {

  float a, b, c, d;
  if ( (t1 < t2) && (t2 < t3) && 
      ( ((x1 <= x2) && (x2 >= x3)) || ((x1 >= x2) && (x2 <= x3)) ) ) {
    d = -(t1 - t2)*(t2 - t3)*(t3 - t1);
    a = (t3*(x2 - x1) + t2*(x1 - x3) + t1*(x3 - x2))/d;
    b = (t3*t3*(x1 - x2) + t2*t2*(x3 - x1) + t1*t1*(x2 - x3))/d;
    c = (x1*t2*t3*(t2-t3) + x2*t1*t3*(t3-t1) + x3*t1*t2*(t1-t2))/d;
    *t_extremum_ptr = -b/(2*a);
    *x_extremum_ptr = c - b*b/(4*a);
    return 1;
  } else {
    mexPrintf("parabolic_extremum_interpolation can not handle these inputs: "
        "t1=%.9f, t2=%.9f, t3=%.9f, x1=%.9f, x2=%.9f, x3=%.9f\n", 
        t1,t2,t3,x1,x2,x3);
    return 0;
  }

}


/******************************************************************************
linear_crossing_interpolation(t1, t2, x1, x2, x_crossing_value, t_out)

Given two points (t1,x1) and (t2,x2) and a crossing value x_cross, such that
t1 < t2 and either (x1 <= x_cross <= x2) or (x1 >= x_cross >= x2), find the
value t_cross at which the line connecting the points crosses x_cross
*******************************************************************************/
int linear_crossing_interpolation(float t1, float t2, float x1, float x2, 
    float x_cross, float* t_cross_ptr) {

  /* check that x_cross is between x1 and x2 */
  if ( (t1 < t2) && (x1 != x2) && 
      ( ((x1 <= x_cross) && (x_cross <= x2)) || 
      ((x1 >= x_cross) && (x_cross >= x2)) ) ) {
    *t_cross_ptr = (t2*(x_cross - x1) + t1*(x2 - x_cross)) / (x2 - x1);
    return 1;
  } else {
    mexPrintf("linear_crossing_interpolation can not handle these inputs: "
        "t1=%.9f, t2=%.9f, x1=%.9f, x2=%.9f, x_cross=%.9f\n",t1,t2,x1,x2,x_cross);
    return 0;
  }

}


/******************************************************************************
copy_to_output(output_data, shape, num_events, num_channels, event_idx)

copy values from the NUM_CHANNELS elements of the struct array SHAPE to the
corresponding array elements of OUTPUT_DATA 
*******************************************************************************/
void copy_to_output(float** output_data, waveform_shape_t* shape, 
    unsigned long num_events, unsigned int num_channels, 
    unsigned long event_idx) {

  unsigned int c; /* channel index */
  unsigned long offset;

  for (c = 0; c < num_channels; c++) {
    offset = c*num_events + event_idx;
    /* implementation here must match the fields of the waveform_shape_t
    struct type to field_names. output_data[f] points to the data array that
    corresponds to field_names[f] */
    *(output_data[0] + offset) = shape[c].trough_amplitude;
    *(output_data[1] + offset) = shape[c].peak_amplitude;
    *(output_data[2] + offset) = shape[c].half_trough_time_before;
    *(output_data[3] + offset) = shape[c].trough_time;
    *(output_data[4] + offset) = shape[c].half_trough_time_after;
    *(output_data[5] + offset) = shape[c].zero_crossing_time;
    *(output_data[6] + offset) = shape[c].half_peak_time_before;
    *(output_data[7] + offset) = shape[c].peak_time;
    *(output_data[8] + offset) = shape[c].half_peak_time_after;
  }

}


/******************************************************************************
usage()

Tell the user how to call this function correctly
*******************************************************************************/
void usage() {

  mexErrMsgTxt(
"MEASURE_SPIKE_SHAPE Measure shape of spike waveforms recorded in NSpike .tt format\n"
"   SHAPE = MEASURE_SPIKE_SHAPE_MEX(X,T,I) finds the peaks, troughs, and \n"
"   zero-crossings of the thresold-triggered spike waveforms specified by "
"   X,T,I.\n"
"   \n"
"   X is a MxCxN array of single-precision floating-point values. C is the \n"
"   number of recording channels in the polytrode, N is the number of \n"
"   trigger events, and M is the number of samples acquired per window.\n"
"   \n"
"   T is a Mx1 vector of single-precision floating-point values, which \n"
"   gives the times of the M samples relative to the first sample in the \n"
"   window. T is assumed to be monotonically increasing.\n"
"   \n"
"   I is an integer scalar between 1 and M, which indexes along T and along \n"
"   the first dimension of X. I specifies the time of the trigger within \n"
"   the M-sample acquisition window. I is expressed with MATLAB-style \n"
"   (i.e., the first element of T corresponds to I = 1).\n"
"   \n"
"   The output argument SHAPE is a struct whose fields are NxC arrays of \n"
"   single-precision floating-point values\n"
  );

}

