% neweegstruct = APPLYFILTER(eegstruct, index, filter, samprate)
%	Applies the filter to the data in the indexed element of eeg struct. 
%	If samprate is not specificied, the sampling rate is preserved from the
%	input data, otherwise the data is resampled at the specificed rate.
%	Returns a structure with starttime, samprate, and data fields

function [ns] = applyfilter(eegstruct, index, filter, samprate)

eeg = eegstruct{index(1)}{index(2)}{index(3)};
data = eeg.data;
datasamprate = eeg.samprate;

resample = 1;
% check to see if samprate was specified or if the specified sampling rate is the same as
% that for the original data
if ((nargin ~= 4) | (samprate == datasamprate))
	resample = 0;
end

% copy the input structure to the output structure
ns = eeg;
ns.data = [];

% apply the filter
newdata = filtfilt(filter.tf.num, filter.tf.den, data);

if (resample)
	% only do 2e6 points at a time
	maxdatalen = 2e6;
	ndatapoints = length(newdata);
	newstarttime = 1;
	for i = 1:maxdatalen:ndatapoints
		dataendtime = min(i + maxdatalen, ndatapoints);
		datatimes = ((i:dataendtime) - 1) / datasamprate;

		newendtime = floor(dataendtime * samprate / datasamprate);
		newtimes = ((newstarttime:newendtime) - 1) / samprate;
		ns.data = [ns.data  ...
					round(interp1(datatimes, newdata(i:dataendtime), ...
						 newtimes, 'linear'))];
		newstarttime = newendtime + 1;
	end
	ns.samprate = samprate;
else
	ns.data = newdata;
end
