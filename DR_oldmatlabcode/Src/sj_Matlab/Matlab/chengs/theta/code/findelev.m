function elev = findelev(feeg,nstd,times)
%
%       nstd:   threshold this many std above baseline 
%
% Find times when the envelop of the filtered signal exceeds mean+n*std.

% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious 
% flucutations in the envelope)
SMOOTHING_WIDTH = 0.002; % 2 ms

% take the Hilbert-transform amplitude envelope of the ripple trace and
% smooth with a Gaussian with standard deviation of SMOOTHING_WIDTH
% which spans 4 SD (note that this uses the gausswin function
% of the Signal Processing toolbox)
gaussian_kernel = gausswin(ceil(8*SMOOTHING_WIDTH*feeg.samprate),4);
% normalize
gaussian_kernel = gaussian_kernel/sum(gaussian_kernel);
% filter the ripple-filtered trace
smoothed_envelope = filtfilt(gaussian_kernel,[1],feeg.data(:,3));


baseline= mean(smoothed_envelope);
dev= std(smoothed_envelope);
threshold= baseline+nstd*dev;

if threshold<baseline; warning('threshold<baseline'); end

fprintf(1, 'baseline= %.2f, threshold= %.2f, std= %.2f\n', baseline, threshold, dev);

flags = smoothed_envelope > threshold;

% figure out the actual times (in seconds) which correspond to these first and
% last eeg samples for each ripple snippet
timestamps = feeg.starttime + ((1:size(feeg.data,1))' - 1 )/feeg.samprate;

elev= interp1(timestamps, flags, times)>0.1;



