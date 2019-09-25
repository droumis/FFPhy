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
%
