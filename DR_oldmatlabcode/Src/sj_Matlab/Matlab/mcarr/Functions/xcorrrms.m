function rms = xcorrrms(time, cc, rmstmax, rmsmin)
%function rms = xcorrrms(time, cc, rmstmax, rmsmin)
% 
% Computes the root mean squre of the unnormalized crosscorrelation histogram 
% from -rmstmax to rmstmax.
%
% if the histogram has less then rmsmin counts in it, rms is returned as NaN;
%
if (sum(cc) < rmsmin)
    rms = NaN;
    return;
end

% find the bins between +/- rmstmax
bins = find(abs(time) <= rmstmax);

% squre the times, mulitple by the bin counts, take the mean and then the sqrt

rms = sqrt(mean(cc(bins) .* (time(bins).^2)));


