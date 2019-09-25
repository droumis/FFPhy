function out = DFAsj_getthetacrosscovLAG_timecondition4(ind, excludetimes, spikes, varargin)
% Ver4 : Starting 10Feb2014 - Sync codes with everyone
% Changed defns of normcorr, etc a bit
% Can change binsize conditions for getting correct lag

% Shantanu, Dec 2013. Add a condition for minimum length of include time segments

% Shantanu Oct 2013 - Similar to DFAsj_getthetacovariogram/ DFAsj_calcxcorrmeasures
% Cross-Cov as in Siapas (2005), Wierzynski (2009);
% Firing rate adjusted. Standardized measure (Z-scored) also returned
% -----------------------------------------------

% Here, the focus is just on the cross-correlation. Usual way is to return the normalized xcorr.

% Here, a shuffled cc for continuous data will be calculated, and raw-shuffled will be returned
% Both the "raw" and "shuffled" will be normalized by the geometric mean of the no of spikes (sqrt(N1*N2)),
% so flank normalization will not be need (normalization by mean rate)
% Alternative is to normalize by subtracting psth1.psth2. But this is not trial-based data
% For continuous data, spiktimes will have to be pushed appropriate amounts - try 5sec+1sec jitter, and 10sec+1sec jitter

% I can also do a cross-covariance (Siapas / Wierzynski) by binning spikes to get histogram, then subtracting mean
% fir rate in each bin for both cells, and then doing xcorr. Siapas uses a diff method, of getting cross-corr,
% then normalizing by:
% Divide by (binsize * Total time of recording) - (expected co-occurence in bin by product of fir rates)

% Calculating shuffle in trial-based data is easier  since you shuffle identity of trials.




% -----------------------------------------------
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

% Values for theta. For ripple corrln, use lower binsize.
% Also, for ripples, you will treat it as trials and get covariogram by <S1.S2> - <P1.P2> (See Brody, 1999)
bin = 0.01; % 10ms bins           %bin = 0.004;    %4 ms bins
sw1 = bin*3; %    % bin = bin*2.5  % 10ms smoothing             % for smoothing corrln. Necessary?
sw2 = 0.2; % Old - used in excesscorr to get baseline corrln - might be useful approach. Not used currently.

edgespikes = 0;
tmax = 0.5; % +/- 500ms for corrln

thrstime = 1; % Minimum length of includetime segments in sec

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
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
            case 'thrstime'
                thrstime = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end




% for each cell we calculate the cross correlation
try
    t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
    t2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data(:,1);
catch
    % if either of those produced an error (meaning no spikes in epoch), we return NaN
    out.index = ind;
    out.T = NaN;
    out.p1p2 = NaN;
    out.corrstruct = NaN;
    out.rawcorr = NaN;
    out.rawcorr_sm = NaN;
    out.normcorr = NaN;
    out.normcorr_sm = NaN;
    out.Zcrosscov = NaN;
    out.crosscov = NaN;
    out.Zcrosscov_sm = NaN;
    out.crosscov_sm = NaN;
    out.Neventscorr = 0;
end

tet1=ind(3);
cell1=ind(4);
tet2=ind(5);
cell2=ind(6);

% % If using exclude periods, looks like this
% % ---------------------------------------
% %apply the exclude rule
% t1inc = t1(find(~isExcluded(t1, excludetimes)));
% t2inc = t2(find(~isExcluded(t2, excludetimes)));
% % Get total time for data
% totaleptime = diff(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange)./10000; % in secs
% excltime = sum(diff(excludetimes'));
% T = totaleptime - excltime;

totaleptime = diff(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange)./10000; % in secs
excltime = sum(diff(excludetimes'));

% Use Include periods instead
% --------------------------
% Get IncludeTimes from flanking edges in excludetimes and epoch start-end
epstend = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange./10000;
incl=[];
incl(:,1) = excludetimes(1:end-1,2);
incl(:,2) = excludetimes(2:end,1);
incl = [epstend(1),incl(1,1) ;incl];
incl = [incl; incl(end,2),epstend(2)];

% Length of include periods
incl_lths = diff(incl')';
% Discard anything < thrstime
discard = find(incl_lths<thrstime);
incl(discard,:)=[];

% Now, use isIncluded instead of a isExcluded
% -----------------------------------------------
%apply the include rule
t1inc = t1(find(isIncluded(t1, incl)));
t2inc = t2(find(isIncluded(t2, incl)));
% Get total time for data
incl_lths = diff(incl')';
T = sum(incl_lths);



% Get raw-cross-correlation
% --------------------------

if ~isempty(t1inc) && ~isempty(t2inc)
    
    xc = spikexcorr(t1inc, t2inc, bin, tmax);
    % if we want to include edge spikes, we need to add in the correlation of the
    % excluded t1 spikes with the included t2spikes
    if ((edgespikes) & (~isempty(xc.time)))
        t1ex = t1(find(isExcluded(t1, excludetimes)));
        if (~isempty(t1ex))
            tmpxc = spikexcorr(t1ex, t2inc, bin, tmax);
            % add these values to the original histogram
            xc.c1vsc2 = xc.c1vsc2 + tmpxc.c1vsc2;
        end
    end
    
    
    % Expected probability
    p1 = xc.nspikes1/T; p2 = xc.nspikes2/T; % Fir rate in Hz
    exp_p = p1*p2; % per sec
    
    % Crosscov
    crosscov = (xc.c1vsc2 ./ (bin*T))-exp_p;
    % Convert to Z-score
    factor = sqrt((bin*T) / exp_p);
    Zcrosscov = crosscov .* (factor);
    
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
    Zcrosscov_sm = smoothvect(Zcrosscov, g1);
    crosscov_sm = smoothvect(crosscov, g1);
    
    % Plot to check
    % figure; hold on; plot(xc.time,corr_sm,'k');
    % figure; hold on; plot(xc.time, rawcorr_sm, 'k');
    % figure; hold on; plot(xc.time,Zcrosscov ,'r'); plot(xc.time,Zcrosscov_sm,'k','LineWidth',2);
    % figure; hold on; plot(xc.time,crosscov_sm,'b');
    
    % Get Nevents in raw correlation from -200ms to 200ms
    bins = find(abs(xc.time)<=0.2);
    Neventscorr = sum(xc.c1vsc2(bins));
    
    % Output
    out.index = ind;
    out.T = T;
    out.p1p2 = exp_p;
    out.corrstruct = xc;
    out.rawcorr = rawcorr;
    out.rawcorr_sm = rawcorr_sm;
    out.normcorr = normcorr;
    out.normcorr_sm = normcorr_sm;
    out.Zcrosscov = Zcrosscov;
    out.crosscov = crosscov;
    out.Zcrosscov_sm = Zcrosscov_sm;
    out.crosscov_sm = crosscov_sm;
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
    out.Zcrosscov = NaN;
    out.crosscov = NaN;
    out.Zcrosscov_sm = NaN;
    out.crosscov_sm = NaN;
    out.Neventscorr = 0;
    return;
end


% Can get PEAK and LATENCY TO PEAK OF ZCrosscov AS WELL. Do in calling function.

% To compute a statistic for Zcrosscov
% ------------------------------------
% As in Siapas (2005), sqrt(2)*erfinv(1-alpha).
% Alpha is 0.01, but should be reduced by no. of bins in -200 to 200ms, = 40
% Strictly, this comparison should be to non-smoothened Zcrosscov, and for Gaussian approximation
% Valus is 3.66 for alpha=0.01/40; and 2.58 for alpha=0.01

% As in Wierzynski, generate 40-dim random unit vectors, say 1000 times.
% Each shuffle is a putative Zcrosscov_shuf. Use the 1000 peak values as a test statistic
% Note that this is the same for allpairs, as long as normal approx is valid.
% Normal approx valid when lambda > 10 (lambda = T*bin*pi*pj). For T=20min, pi*pj > 0.83
% If not, geerate Poisson vectors (J-lambda)/sqrt(lambda), where J is Poisson with intensity lambda
% THis is just the standardized measure for Poisson process, since mean and variance is lambda


% if (appendindex)
%     out.index = ind;
% end



