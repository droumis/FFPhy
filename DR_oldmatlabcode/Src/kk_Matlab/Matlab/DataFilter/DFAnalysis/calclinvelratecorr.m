function out = calclinvelratecor(index, excludetimes, spikes, linpos, varargin)
% out = calclinvelratecor(index, excludetimes, spikes, linpos, options)
% Calculates the correlation between velocity and firing rate in 2 second bins,
% considering only non excluded times and bins with at least 1 spike
%
% Options:
%   'binsize', binsize -- set binsize in seconds
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%

appendindex = 0;
binsize = 2;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'binsize'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if (isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data))
    out = [-10 -10];
    if (appendindex)
	out = [out index];
    end
    return;
end

s = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);

vtimes = linpos{index(1)}{index(2)}.statematrix.time;
rws = 1:length(vtimes)';

tmpind = sub2ind(size(linpos{index(1)}{index(2)}.statematrix.linearVelocity), rws', linpos{index(1)}{index(2)}.statematrix.referenceWell);

v = abs(linpos{index(1)}{index(2)}.statematrix.linearVelocity(tmpind));

% get rid of excluded spikes
s = s(find(~isExcluded(s, excludetimes)));
if (isempty(s))
    out = [-10 -10];
    if (appendindex)
	out = [out index];
    end
    return;
end

v(find(isExcluded(vtimes, excludetimes))) = Inf;

% create a set of time bins for the spikes and the velocity
tbins = vtimes(1):binsize:vtimes(end);

% calculate the average velocity for each bin with a boxcar filter
nf = round(binsize / (vtimes(2) - vtimes(1)));
boxcar = ones(nf, 1) ./ nf;
vsmooth = smoothvect(v, boxcar);
% the velocity for each bin is the velocity in the center of the bin
vbin = vsmooth((1+nf/2):nf:end); 

% get the spike rate in each bin
bins = vtimes((1+nf/2):nf:end); 
sc = hist(s, bins);
sc = sc ./ binsize;


valid = find((sc > 0) & isfinite(vbin'));
% make sure there are at least 25 points
if (length(valid) < 25)
    out = [-10 -10];
    if (appendindex)
	out = [out index];
    end
    return;
end


%scatter(sc(valid), vbin(valid))
%pause
[c p] = corrcoef(sc(valid), vbin(valid));
if (~isfinite(c(1,2)))
    out = [-10 -10]
else
    out = [c(1,2) p(1,2)];
end
if (appendindex)
    out = [out index];
end
