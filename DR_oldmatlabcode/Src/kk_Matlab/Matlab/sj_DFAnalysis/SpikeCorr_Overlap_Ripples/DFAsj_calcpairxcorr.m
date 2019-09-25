function out = DFAsj_calcpairxcorr(ind, excludetimes, spikes, varargin)
% Shantanu - Renamed from calcpairxcorr
% Changing to return only correlation
% peak distance returned by peakdistance_traj or peakdistance_linpos

%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
%
% Shantanu:
% Changing from working for multicellanal to singlecellanal
% Relative Spike Timing vs Peak Distance
% This is also ideal for sequence compression index. See Sens paper - Fig 7
% Need to limit C to 100ms and add a condition of 30 spikes
% Also, could use CoM rather than peak place fields. There must be another function
% calcxcorrmeasure is a better function with more control which also uses
% linfields (trajdata) to calculate overlap rather than getting place fields from scratch
% Only returns time lag - not the excess correlation
% Does not call xcorrms for RMS time lag, rather it is time of max cross-corr
% Place field comparison is distance between peaks
%
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.

bin = 0.0001;
%bin = 0.001;
sw1 = 0.005; % 5ms smoothing
sw2 = 0.1; % 100ms smoothing. Defunct for this code 
tmax = .5;

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            
            case 'bin'
                bin = varargin{option+1};
            case 'tmax'
                tmax = varargin{option+1};
            case 'plotxcorr'
                plotxcorr = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

if (nargin < 5)
    binsize = 2;
end

warning('OFF','MATLAB:divideByZero');

try
    spikes1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
    spikes1 = spikes1(find(~isExcluded(spikes1,excludetimes)));
    spikes2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data(:,1);
    spikes2 = spikes2(find(~isExcluded(spikes2,excludetimes)));
catch
    spikes1 = [];
    spikes2 = [];
end

tet1=ind(3); cell1=ind(4); tet2=ind(5); cell2=ind(6);

if (~isempty(spikes1) && ~isempty(spikes2))
    
    tmpcorr = spikexcorr(spikes1, spikes2, bin, tmax);
    maxval = max(tmpcorr.c1vsc2); % This is an integer: no of spikes in bin
    out.time = [];
    for addup = 1:maxval
        out.time = [out.time tmpcorr.time(find(tmpcorr.c1vsc2 == addup))];
    end
    
    % Get smooth correlation
    [~, ~, normsmoothcorr, basecorr] = ...
    excesscorr(tmpcorr.time, tmpcorr.c1vsc2, tmpcorr.nspikes1, tmpcorr.nspikes2, sw1, sw2);
    maxval = max(normsmoothcorr);
    out.maxtime = [];
    for addup = 1:length(maxval)
        out.maxtime = [out.maxtime tmpcorr.time(find(normsmoothcorr == maxval(addup)))];
    end
    
    out.c1vsc2 = tmpcorr.c1vsc2;
    out.normsmoothcorr = normsmoothcorr;
    out.timebase = tmpcorr.time;
    out.nspikes1 = length(spikes1);
    out.nspikes2 = length(spikes2);
    %out.ind = ind(i,:);
    out.index = ind;
    if isempty(out.time),
        out.time = NaN;
    end
    if isempty(out.maxtime),
        out.maxtime = NaN;
    end
else
    out.time = NaN;
    out.maxtime = NaN;
    out.c1vsc2 = NaN;
    out.normsmoothcorr = NaN;
    out.timebase = NaN;
    out.nspikes1 = NaN;
    out.nspikes2 = NaN;
    out.index = ind;
end

