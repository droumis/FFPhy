% peaks = FINDEEGPEAKS(eegstruct, index, maxfreq)
%	Returns the locations of the peaks in the eeg by finding all
%	locations where the slope of the the eeg signals changes sign from
%	positive to negative
%   maxfreq gives the minimum frequencies of the peaks in Hz

function [peaks] = findapeaks(eegstruct, index, maxfreq)


data = eegstruct{index(1)}{index(2)}{index(3)}.data;

% calculate an estimate of the first derivative
deriv = diff(data(1:end-1));
derivtmp = diff(data(2:end));

peaks = (find((deriv > 0) & (derivtmp <= 0))+1);


if (nargin < 3)
	return
else
	% get rid of the all but the first of each set of high frequeny peaks
	minpdiff = eegstruct{index(1)}{index(2)}{index(3)}.samprate / maxfreq;
	ok = find(diff(peaks) > minpdiff);
	peaks = [peaks(1) peaks(ok+1)];
end
