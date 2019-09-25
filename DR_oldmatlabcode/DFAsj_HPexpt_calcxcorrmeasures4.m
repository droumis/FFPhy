function out = DFAsj_HPexpt_calcxcorrmeasures4(ind, excludetimes, spikes, varargin)

% SJ - 08/2012.
% Shantanu - Renamed from calcxcorrmeasures
% function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
%
% Shantanu:
% Better than calcpairxcorr -
% More control and gets overlap from linfields-trajdata rather than calculate from scratch
% Uses xcorrms for RMS time lag
% Uses excesscorr for excess correlation
% Place field comparison is overlap - needs linfields
%
%
% Shantanu - What I am doing in sj version
% - For sleep - there is no overlap measure. So get rid of that. No need for linfields.
%   Overlap during run is calculated by DFAsj_calcoverlap called by main script
%   DFSsj_xcorrmeasure..
% - I want to return out as structure with fields
% - excesscorr is not the main output I am looking for here. I want to
% return the normalized correlation in the -100ms to +100ms window



% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.
%   'edgespikes', 0 or 1
%		  0 will enforce the exclude periods and will thus
%   	    	    exclude all spikes outside of the exclude windows. (default)
%		  1 will include all spikes where at least one of the spikes is
%		    not excluded
%   'sw1', n	The smaller smoothing width for the excess correlation measure
%		Default 0.005 (5 ms)
%   'sw2', n	The larger smoothing width for the excess correlation measure
%		Default 0.250 (250 ms)
%   'rmstmax', n The maximum absolute time to be used for the rms calculation
%		Default 0.1 (100 ms)
%   'rmsmincounts', n The minimum number of counts in the histogram for the rms
%  		measurement to be carried out.
%		Default 10.
%   'plotxcorr', 1 or 0, default 0
%       plot xcorrelation if plotxcorr = 1
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%

%load /data/loren/mkarlsso/Ten/tenlinfields
%linfields = trajdata;

appendindex = 1; % always appendindex
forripples = 1; % Default is forripples
forbeh = 0; % Beh time scale correlations

% These do not change for different conditions. For Rip Dis. bin was 2ms, and sw1 was 5 ms
bin = 0.01; % 0.01 = 10ms bins    0.002=2ms bins
sw1 = bin*3; % for smoothing     %sw1 = bin*3; % 6 ms smoothing
edgespikes = 0;
% For ripple correlation, do
tmax = 0.5;

% This is legacy - for baseline correlation like Sen's paper
sw2 = 0.1; % For ripples. 100ms smoothing in excesscorr - big window for baseline
rmstmax = 0.05; % For ripples. Look for concentration within 50 ms
rmsmincounts = 10;


for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'forripples'
                forripples = varargin{option+1};
            case 'forbeh'
                forbeh = varargin{option+1};
            case 'bin'
                bin = varargin{option+1};
            case 'tmax'
                tmax = varargin{option+1};
            case 'edgespikes'
                edgespikes = varargin{option+1};
            case 'sw1'
                sw1 = varargin{option+1};
            case 'sw2'
                sw2 = varargin{option+1};
            case 'rmstmax'
                rmstmax = varargin{option+1};
            case 'plotxcorr'
                plotxcorr = varargin{option+1};
            case 'rmsmincounts'
                rmsmincounts = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end


if forripples==0
    tmax = 1;
    sw2 = 0.25;
    rmstmax = 0.1;
    rmsmincounts = 10;
end

if forbeh==1
    tmax = 5;
    sw2 = 1;
    rmstmax = 0.5;
    rmsmincounts = 10;
end


% Get total time for data
totaleptime = diff(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange)./10000; % in secs
excltime = sum(diff(excludetimes'));
T = totaleptime - excltime;

% for each cell we calculate the cross correlation
try
    t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
    t2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data(:,1);
catch
    % if either of those produced an error (meaning no spikes in epoch), we return NaN
    out.index = ind;
    out.T = T;
    out.p1p2 = NaN;
    out.corrstruct = NaN;
    out.rawcorr = NaN;
    out.rawcorr_sm = NaN;
    out.normcorr = NaN;
    out.normcorr_sm = NaN;
    out.Neventscorr = 0;
    % Legacy
    %out.basecorr = NaN;
    %out.ec = NaN;
    %out.rms = NaN;
    %out.probcoa = NaN;
    %out.coactivez = NaN;
    return;
end

tet1=ind(3);
cell1=ind(4);
tet2=ind(5);
cell2=ind(6);

%apply the exclude rule
t1inc = t1(find(~isExcluded(t1, excludetimes)));
t2inc = t2(find(~isExcluded(t2, excludetimes)));

% Get raw-cross-correlation
% --------------------------

if ~isempty(t1inc) && ~isempty(t2inc)
    
    xc = spikexcorr(t1inc, t2inc, bin, tmax);
    % if we want to include edge spikes, we need to add in the correlation of the
    % excluded t1 spikes with the included t2spikes
    if ((edgespikes) & (~isempty(xc.time)))
        t1ex = t1(find(isExcluded(t1, excludetimes)));
        if (~isempty(t1ex))
            tmpxc = spikexcorr(t1ex, t2inc_forcorr, bin, tmax);
            % add these values to the original histogram
            xc.c1vsc2 = xc.c1vsc2 + tmpxc.c1vsc2;
        end
    end
    
    % Expected probability
    p1 = xc.nspikes1/T; p2 = xc.nspikes2/T; % Fir rate in Hz
    exp_p = p1*p2; % per sec
    
    %     % Crosscov
    %     crosscov = (xc.c1vsc2 ./ (bin*T))-exp_p;
    %     % Convert to Z-score
    %     factor = sqrt((bin*T) / exp_p);
    %     Zcrosscov = crosscov .* (factor);
    
    %Raw Corr  Normalize by geometric mean: Units will be coincidences/ spike
    rawcorr = xc.c1vsc2;
    normcorr = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);
    
    % Smooth values as well
    try
        nstd=round(sw1/(xc.time(2) - xc.time(1))); % will be 3 std
        %g1 = gaussian(nstd, 2*nstd+1);    %sigma=3, npoints=7 (+/- 3 bins
        % Can also do a 3 point gaussian with sigma=3. Almost a boxcar.
        g1 = gaussian(nstd, 2*nstd+1);
        
    catch
        keyboard;
    end
    rawcorr_sm = smoothvect(rawcorr, g1);
    normcorr_sm = smoothvect(normcorr, g1);
    
    % Plot to check
    %figure; hold on; plot(xc.time,normcorr,'r'); plot(xc.time, normcorr_sm, 'k','LineWidth',2);
    
    % Get Nevents in raw correlation from -200ms to 200ms
    bins = find(abs(xc.time)<=0.2);
    Neventscorr = sum(xc.c1vsc2(bins));
    
    
    % % Legacy - COmment
    % % compute the excess correlation at 0 lag, and get normalized crosscorrln
    % [ec, normrawcorr, normsmoothcorr, basecorr] = ...
    %     excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
    % rms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);
    %
    % % get the probability that both were coactive in the included intervals
    % n1 = nInInterval(t1inc, excludetimes, 0);
    % n2 = nInInterval(t2inc, excludetimes, 0);
    % probcoa = mean(logical(n1) & logical(n2));
    % % get the zscore for this probability of coactivity
    % coactivez = coactivezscore(n1, n2);
    
    
    % Output
    % -------
    out.index = ind;
    % Corrln
    % ------
    out.T = T;
    out.p1p2 = exp_p;
    out.corrstruct = xc;
    out.rawcorr = rawcorr;
    out.rawcorr_sm = rawcorr_sm;
    out.normcorr = normcorr;
    out.normcorr_sm = normcorr_sm;
    out.Neventscorr = Neventscorr;
    
else
    out.index = ind;
    out.T = T;
    out.p1p2 = NaN;
    out.corrstruct = NaN;
    out.rawcorr = NaN;
    out.rawcorr_sm = NaN;
    out.normcorr = NaN;
    out.normcorr_sm = NaN;
    out.Neventscorr = 0;
end



