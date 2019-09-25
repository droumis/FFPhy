function [ec, cc, cc1, cc2] = excesscorr(time, cc, ns1, ns2, sw1, sw2)
%function ec = excesscorr(time, xcorrhist, nspikes1, nspikes2, smoothingwidth1, smoothingwidth2)
% 
% Output - cc = Normalized cc, cc1 = smooth cc, cc2 = baseline cc
%
% Computes the excess correlation of the crosscorrelation histogram as follows:
% The histogram is normalized by sqrt(nspikes1 * nspikes2).  It is then 
% smoothed, first using a gaussian with stdev smoothingwidth1, and separately
% using a gaussian with stdev smoothingwidth2.
% The excess correlation is defined as the difference the two smoothed
% normalized histogram (width1 smoothed - width2 smoothed) at zero lag.

% normalize the cross correlegram
if (isempty(cc))
    % this occurs when there is no data, so we return NaN
    ec = NaN;
    cc = NaN;
    cc1 = NaN;
    cc2 = NaN;
    return;
end

cc = cc ./ sqrt(ns1 * ns2);

% smooth it using sw1
nstd=round(sw1/(time(2) - time(1)));
g1 = gaussian(nstd, 5*nstd+1);
cc1 = smoothvect(cc, g1);

nstd=round(sw2/(time(2) - time(1)));
g2 = gaussian(nstd, 6*nstd+1);
cc2 = smoothvect(cc, g2);

% find the minium bin value and then the index or indices of the closes bins
mn = min(abs(time));

zerobins = find(abs(time) == mn);

ec = mean(cc1(zerobins) - cc2(zerobins));

