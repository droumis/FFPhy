function [eeg, times,samplingRate] = eeg_window(filename, eegTimes, windowLength, varargin)
% function [eeg, times] = eeg_window(filename, eegTimes, windowLength, varargin)

[otherArgs] = procOptions(varargin);


% Error checking
MAX_FILENAME = 200;
if length(filename) > (MAX_FILENAME-1)
  error('eeg_window:maxFilename','Filename: %s too long.',filename);
end

eegTimes = sort(eegTimes(:)); % make times a Nx1

[eeg, times, samplingRate] = eeg_window_c(filename, eegTimes, windowLength);
eeg = transpose(eeg);
