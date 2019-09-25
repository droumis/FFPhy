% peaks = FINDEEGPEAKS(data)
%	Returns the locations of the peaks in the data by finding all
%	locations where the slope of the data changes sign from
%	positive to negative

function [peaks] = findpeaks(data);

% calculate an estimate of the first derivative
deriv = diff(data(1:end-1));
derivtmp = diff(data(2:end));

peaks = (find((deriv >= 0) & (derivtmp <= 0))+1);
