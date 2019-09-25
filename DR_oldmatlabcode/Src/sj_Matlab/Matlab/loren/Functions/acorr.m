%	a = ACORR(spikestruct, bin, tmax)
%       returns the normalized autocorrelation of the spiketrain in spikestruct
%       bin and tmax should be in units of msec

function [ac] = acorr(spikestruct, bin, tmax)

descript = sprintf('Autocorrelation of %s, bin = %lf, tmax = %lf',  ...
					inputname(1), bin, tmax);
ac.descript = char(descript, spikestruct.descript);

nbins = floor(tmax / bin);
a = zeros(2 * nbins + 1, 2);
s = spikestruct.data(:,1) * 1e3; % convert to ms
nspikes = length(s);

if (nspikes < 10) 
	% too few spikes for a decent acorr
	ac = [];
	return
end

%bin the spikes
numsbins = (s(nspikes) - s(1)) / bin;

bins = bin/2 : bin : (numsbins + 1/2) * bin;

%[sbinned tmp] = hist(s - s(1), numsbins);
[sbinned tmp] = hist(s - s(1), bins);


a(:,1) = [-1 * nbins * bin:bin:nbins*bin]'; 

len = length(sbinned);
for i = 1 : nbins
	a(nbins+1+i, 2) = sbinned(1:len-i) * sbinned(1+i:len)';
	a(nbins+1-i, 2) = a(nbins+1+i,2);
end
a(nbins+1,2) = nspikes;
a(:,2) = a(:,2) ./ nspikes;
ac.data = a;
