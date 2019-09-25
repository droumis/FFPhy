function [epochs, samplingrate, numsamples] ...
  = epoch_extract(filename, verboseFlag);

% Caleb Kemere 2008-09-03
% 
% function [epochs, samplingrate, numsamples] ...
%   = epoch_extract(filename, verboseFlag);
% 
% Function loadeeg:
%  Gets epochs from an EEG file
%
%  Input arguments:
%    "filename" - Filename of .eeg file. Should be fully
%      qualified for fopen.
%    "verboseFlag" (optional) - Set to 1 to report number of
%      records in file without opening. Set to 2 to report progress
%      (10000's of records read).
%
%  Output arguments:
%    "epochs" - start and end times for epochs
%    "samplingrate" (optional) - EEG sampling rate (Hz)
%    "numsamples" (optional) - number of data points per EEG block
%

[epochs, samplingrate, numsamples] = epoch_extract_c(filename, verboseFlag);

first_zero = find(epochs(1,:)==0,1);
epochs_new = epochs(:,1:first_zero-1);
clear('epochs');
epochs = transpose(epochs_new);
