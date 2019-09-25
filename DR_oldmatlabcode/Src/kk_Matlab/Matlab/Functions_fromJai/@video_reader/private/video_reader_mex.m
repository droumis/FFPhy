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
%Written by SMK, 2009 January 1.
%

