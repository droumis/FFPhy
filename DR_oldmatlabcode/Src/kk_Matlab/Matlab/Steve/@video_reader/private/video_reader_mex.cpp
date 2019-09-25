/*=============================================================================
%VIDEO_READER_MEX Interface for reading NSpike video data.
%
%   VIDEO_READER_MEX is an interface for manipulating a C++ VideoReader object.
%   The memory used by VIDEO_READER_MEX is not visible in the MATLAB workspace
%   and may not be cleared as expected when an error is encountered. The
%   video_reader MATLAB class (also written by SMK) is a friendly wrapper
%   around this function.
%
%   VIDEO_READER_MEX('open',mpeg_stream_filename,mpeg_index_filename,
%   timestamps_filename) initializes the interface. There is no return value.
%   An error will be raised and the mex file will exit if the arguments are
%   invalid.
%
%   VIDEO_READER_MEX('close') closes the interface and clears memory. There is
%   no return value.
%
%   VIDEO_READER_MEX('width'), VIDEO_READER_MEX('height'),
%   VIDEO_READER_MEX('aspect_ratio_numerator'),
%   VIDEO_READER_MEX('aspect_ratio_denominator') return the corresponding
%   attributes of the mpeg picture as scalar doubles.
%
%   VIDEO_READER_MEX('number_of_frames') returns the number of video frames in
%   the mpeg stream, which must equal the number of file offsets in the mpeg
%   index as well as the number of video timestamps.
%
%   VIDEO_READER_MEX('all_timestamps') returns a vector of uint32 timestamps
%   read from timetamps_filanem.
%
%   VIDEO_READER_MEX('seek',offset,origin) updates the internal state of the
%   underlying VideoReader object so that it is populated with data
%   corresponding to the video frame in the stream position specified by
%   (offset, origin). The offset argument follows the C++ convention for
%   file i/o streams, and origin argument can assume the following string or
%   integer values:
%     'beg' or -1: beginning of stream
%     'cur' or  0: current position in stream
%     'end' or +1: end of stream
%   
%   VIDEO_READER_MEX('current_index') returns the current stream position,
%   counting from the first video frame according to MATLAB convention. Thus,
%   an index of 1 corresponds to the first video frame, and an index of
%   VIDEO_READER_MEX('number_of_frames') corresponds to the last video frame.
%
%   VIDEO_READER_MEX('current_timestamp') returns the uint32 timestamp
%   corresponding to the current stream position.
%
%   VIDEO_READER_MEX('current_picture') returns the uint8 pixel values of the
%   picture at the current stream position. This return value is a
%   3-dimensional array of size 
%     [VIDEO_READER_MEX('height') VIDEO_READER_MEX('width') 3]
%   oriented to match MATLAB's image(...) function
%
%Depends on:
%   VideoReader C++ class (written by SMK) and standard C++ libraries
%
%   libmpeg2-4-dev
%   mpeg2dec
%
%  compile as mex -lmpeg2 -lmpeg2convert video_reader_mex.cpp VideoReader.cpp
%
%Written by SMK, 2009 January 1.
=============================================================================*/

#include "mex.h"
#include <cstring>
#include "VideoReader.h"

/* NSpike timestamp clock ticks at 10kHz */
const double TS_PER_SEC = 10000.0;

/* this static object behaves as a singleton */
static VideoReader video_reader;

const char* usage = "Usage:\n"
    "video_reader_mex('open',mpeg_stream_filename, "
    "mpeg_offset_filename, timestamps_filename)\n"
    "video_reader_mex('close')\n"
    "video_reader_mex('seek', offset, origin)\n"
    "lhs = video_reader_mex('width')\n"
    "lhs = video_reader_mex('height')\n"
    "lhs = video_reader_mex('aspect_ratio_numerator')\n"
    "lhs = video_reader_mex('aspect_ratio_denominator')\n"
    "lhs = video_reader_mex('number_of_frames')\n"
    "lhs = video_reader_mex('current_picture')\n"
    "lhs = video_reader_mex('current_index')\n"
    "lhs = video_reader_mex('current_timestamp')\n"
    "lhs = video_reader_mex('all_timestamps')\n";

/* method wrappers */
static void cleanup();
void open(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);
void close(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);
void seek(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);
void width(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);
void height(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);
void aspect_ratio_numerator(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);
void aspect_ratio_denominator(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);
void number_of_frames(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);
void current_picture(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);
void current_index(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);
void current_timestamp(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);
void all_timestamps(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]);

static void cleanup() {
  // do I need to do anything to delete video_reader ?!
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  mexAtExit(cleanup);

  mwSize char_length;
  char tmp_char[100];

  if ((nrhs < 1) || (!mxIsChar(prhs[0]))) {
    mexErrMsgTxt(usage);
  }
  char_length = (mwSize) (1 + mxGetM(prhs[0]) * mxGetN(prhs[0]) * 
      sizeof(mxChar));
  if (mxGetString(prhs[0], tmp_char, char_length)) {
    mexErrMsgTxt(usage);
  } else if (strcmp(tmp_char, "open") == 0) {
    open(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "close") == 0) {
    close(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "seek") == 0) {
    seek(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "width") == 0) {
    width(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "height") == 0) {
    height(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "aspect_ratio_numerator") == 0) {
    aspect_ratio_numerator(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "aspect_ratio_denominator") == 0) {
    aspect_ratio_denominator(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "number_of_frames") == 0) {
    number_of_frames(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "current_picture") == 0) {
    current_picture(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "current_index") == 0) {
    current_index(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "current_timestamp") == 0) {
    current_timestamp(nlhs, plhs, nrhs, prhs);
  } else if (strcmp(tmp_char, "all_timestamps") == 0) {
    all_timestamps(nlhs, plhs, nrhs, prhs);
  } else {
    mexErrMsgTxt(usage);
  }

}


void open(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  mwSize char_length;
  char tmp_char[1000];
  std::string mpeg_stream_filename;
  std::string mpeg_index_filename;
  std::string timestamps_filename;

  if ((nrhs != 4) || (nlhs != 0)) { 
    mexErrMsgTxt(usage);
  }
  if (video_reader.ready()) {
    mexErrMsgTxt("object is already initialized\n");
  }
  if (!mxIsChar(prhs[1])) {
    mexErrMsgTxt("Second argument must be .mpeg file name\n");
  }
  char_length = (mwSize) (1 + 
      mxGetM(prhs[1]) * mxGetN(prhs[1]) * sizeof(mxChar));
  if (mxGetString(prhs[1], tmp_char, char_length)) {
    mexErrMsgTxt("Could not read input argument into char array\n");
  }
  mpeg_stream_filename.assign(tmp_char, (std::string::size_type) char_length);

  if (!mxIsChar(prhs[2])) {
    mexErrMsgTxt("Third argument must be .mpegoffset file name\n");
  }
  char_length = (mwSize) (1 + 
      mxGetM(prhs[2]) * mxGetN(prhs[2]) * sizeof(mxChar));
  if (mxGetString(prhs[2], tmp_char, char_length)) {
    mexErrMsgTxt("Could not read input argument into char array\n");
  }
  mpeg_index_filename.assign(tmp_char, (std::string::size_type) char_length);

  if (!mxIsChar(prhs[3])) {
    mexErrMsgTxt("Third input argument to constructor must be string\n");
  }
  char_length = (mwSize) (1 + 
      mxGetM(prhs[3]) * mxGetN(prhs[3]) * sizeof(mxChar));
  if (mxGetString(prhs[3], tmp_char, char_length)) {
    mexErrMsgTxt("Could not read input argument into char array\n");
  }
  timestamps_filename.assign(tmp_char, (std::string::size_type) char_length);

  video_reader.open(mpeg_stream_filename, mpeg_index_filename, 
      timestamps_filename);
  if (video_reader.ready()) {
    mexLock(); // keep this mex-file locked until the close method is called
  }
  // set to beginning
  video_reader.seek_frame(0, std::ios::beg);

} 


void close(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  if ((nrhs != 1) || (nlhs != 0)) { 

    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  video_reader.close();
  if (!video_reader.ready()) {
    mexUnlock(); // unlock after closing
  }
}


void seek(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  
  frameoff_t offset;
  mwSize char_length;
  char tmp_char[10];
  std::ios::seekdir origin;

  if ((nrhs != 3) || (nlhs != 0)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }

  if (mxIsNumeric(prhs[1]) && (mxGetNumberOfElements(prhs[1]) == 1)) {
    offset = (frameoff_t) mxGetScalar(prhs[1]);
  }
  /* origin can be specified as string whose legal values are
    'beg' -1: beginning of stream
    'cur'  0: current position in stream
    'end'  1: end of stream
    Note that these use the C++ convention of starting from 0 and ending with (n-1)
  */
  if (mxIsChar(prhs[2])) {
    char_length = (mwSize) (1 + mxGetM(prhs[2]) * mxGetN(prhs[2]) * 
        sizeof(mxChar));
    if (mxGetString(prhs[2], tmp_char, char_length)) {
      mexErrMsgTxt("Could not read input argument into char array\n");
    }
    if (strcmp(tmp_char, "beg") == 0) {
      origin = std::ios::beg;
    } else if (strcmp(tmp_char, "cur") == 0) {
      origin = std::ios::cur;
    } else if (strcmp(tmp_char, "end") == 0) {
      origin = std::ios::end;
    } else {
      mexErrMsgTxt(usage);
    }
  } else if (mxIsNumeric(prhs[2]) && (mxGetNumberOfElements(prhs[2]) == 1)) {
    switch ((int) mxGetScalar(prhs[2])) {
      case -1: {
        origin = std::ios::beg;
        break;
      }
      case 0: {
        origin = std::ios::cur;
        break;
      }
      case 1: {
        origin = std::ios::end;
        break;
      }
    }
  } else {
    mexErrMsgTxt(usage);
  }
  if (!video_reader.seek_frame(offset, origin)) {
    mexErrMsgTxt("out of bounds\n");
  }
}


void width(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  if ((nrhs != 1) || (nlhs > 1)) { 

    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  plhs[0] = mxCreateDoubleScalar((double) video_reader.width());

}


void height(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  if ((nrhs != 1) || (nlhs > 1)) { 

    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  plhs[0] = mxCreateDoubleScalar((double) video_reader.height());

}


void aspect_ratio_numerator(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  plhs[0] = mxCreateDoubleScalar(
      (double) video_reader.aspect_ratio_numerator());

}


void aspect_ratio_denominator(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  plhs[0] = mxCreateDoubleScalar(
      (double) video_reader.aspect_ratio_denominator());

}


void number_of_frames(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  plhs[0] = mxCreateDoubleScalar((double) video_reader.number_of_frames());

}


void current_picture(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  const uint8_t* data;
  data = video_reader.current_picture();
  mwSize dims[3];

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  } 
  if (data != 0) {
    dims[0] = (mwSize) video_reader.height();
    dims[1] = (mwSize) video_reader.width();
    dims[2] = (mwSize) 3; // each pixel has a RGB triplet
    plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    // pointer to the data field of CData
    uint8_t* out_ptr;
    out_ptr = (uint8_t*) mxGetData(plhs[0]);
    // copy values
    unsigned long int i, j, k; 
    for (i = 0; i < dims[0]; i++) {
      for (j = 0; j < dims[1]; j++) {
        for (k = 0; k < dims[2]; k++) {
          out_ptr[i + dims[0]*(j + dims[1]*k)] = 
              data[dims[2]*(dims[1]*i + j) + k];
        }
      }
    }
  } else {
    dims[0] = (mwSize) 0;
    dims[1] = (mwSize) 0;
    dims[2] = (mwSize) 3;
    plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
  }

}
  

void current_index(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  }
  // +1.0 is to make indexing consistent with MATLAB convention
  plhs[0] = mxCreateDoubleScalar(1.0 + (double) video_reader.current_index());

}


void current_timestamp(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  } 
  plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
  memcpy(mxGetData(plhs[0]), &video_reader.current_timestamp(),
      1*sizeof(uint32_t));

}


void all_timestamps(int nlhs, mxArray* plhs[], int nrhs, 
    const mxArray* prhs[]) {

  if ((nrhs != 1) || (nlhs > 1)) { 
    mexErrMsgTxt(usage);
  }
  if (!video_reader.ready()) {
    mexErrMsgTxt("object is not ready\n");
  } 
  plhs[0] = mxCreateNumericMatrix(video_reader.number_of_frames(), 1, 
    mxUINT32_CLASS, mxREAL);
  uint32_t* out_ptr;
  out_ptr = (uint32_t*) mxGetData(plhs[0]);
  /* copy data from C++ STL vector to matlab vector */
  uint32_t i;
  for (i = 0; i < video_reader.number_of_frames(); i++) {
    out_ptr[i] = video_reader.all_timestamps()[i];
  }

}


