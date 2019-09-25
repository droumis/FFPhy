function out = DFAsj_getthetacovariogram(ind, excludetimes, spikes, varargin)

% Shantanu Oct 2013 - Similar to DFAsj_calcxcorrmeasures
% Here, the focus is just on the corss-correlation. Usual way is to return the normalized xcorr.
% Attempt to get covaiogram by shuffling. See various refs, including Brody, 1999

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
bin = 0.01; % 10 ms bins
sw1 = bin*3; % for smoothing corrln. Necessary?
sw2 = 0.2; % Old - used in excesscorr to get baseline corrln - might be useful approach. Not used currently.

edgespikes = 0;
tmax = 0.5; % +/- 500ms for corrln

shufflelag1 = 1; jitter1 = 1; 
shufflelag2 = 5; jitter2 = 2; 

lags = [1, 2, 5];

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
            case 'shufflelag1'
                shufflelag1 = varargin{option+1};
            case 'shufflelag2'
                shufflelag2 = varargin{option+1};
            case 'jitter1'
                jitter1 = varargin{option+1};
            case 'jitter2'
                jitter2 = varargin{option+1};
            case 'lags'
                lags = varargin{option+1};
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
    out.corr = NaN;
    out.rawcorr = NaN;
    out.shufcorr = NaN;
    out.covgm = NaN;
    out.rawcorr_sm = NaN;
    out.shufcorr_sm = NaN;
    out.covgm_sm = NaN;
    out.Neventscorr = 0;
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
% Normalize by geometric mean: Units will be coincidences/ spike
rawcorr = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);
% Smooth 
nstd=round(sw1/(xc.time(2) - xc.time(1))); g1 = gaussian(nstd, 5*nstd+1);
rawcorr_sm = smoothvect(rawcorr, g1);

% Now get the shuffle corr for both lags. This will take time. Do only 1 lag for now
% -------------------------------------------------------------------------------------------------
% Shuffling with jitter does not make sense. You will flatten corrln. Do only once, with multiple lags

for i = 1:length(lags)
    push = lags(i);
    currt2 = t2inc + push;
    sxc = spikexcorr(t1inc, currt2, bin, tmax);    
    if ((edgespikes) & (~isempty(sxc.time)))
        t1ex = t1(find(isExcluded(t1, excludetimes)));
        if (~isempty(t1ex))
            tmpxc = spikexcorr(t1ex, currt2, bin, tmax);
            % add these values to the original histogram
            sxc.c1vsc2 = sxc.c1vsc2 + tmpxc.c1vsc2;
        end
    end
    % Normalize
    shufcorr(i,:) = sxc.c1vsc2 ./ sqrt(sxc.nspikes1 * sxc.nspikes2);  
    % Now get the covariogram
    covgm(i,:) = rawcorr - shufcorr(i,:);
    % Smooth 
    shufcorr_sm(i,:) = smoothvect(shufcorr(i,:), g1);
    covgm_sm(i,:) = smoothvect(covgm(i,:), g1);  
end    

% %Plotting to check
clr = ['r','b','c'];
figure; hold on; plot(xc.time, rawcorr_sm,'k-'); 
for i=1:length(lags)    
    plot(xc.time, shufcorr_sm(i,:),[clr(i) '--']); plot(xc.time, covgm_sm(i,:),[clr(i) '-']);
end

% Shuffling multiple times
% -----------------------
shufflelag1 = lags(1); jitter1 = lags(1)/2;
shufcorrt = [];
for i=1:100    
    r = randperm(100);
    push = shufflelag1 + (r(1)/100)*jitter1; %To push 5sec + random time in 0.01-1sec
    currt2 = t2inc + push;
    sxc = spikexcorr(t1inc, currt2, bin, tmax);    
    if ((edgespikes) & (~isempty(sxc.time)))
        t1ex = t1(find(isExcluded(t1, excludetimes)));
        if (~isempty(t1ex))
            tmpxc = spikexcorr(t1ex, currt2, bin, tmax);
            % add these values to the original histogram
            sxc.c1vsc2 = sxc.c1vsc2 + tmpxc.c1vsc2;
        end
    end
    % Normalize
    shufcorrt(i,:) = sxc.c1vsc2 ./ sqrt(sxc.nspikes1 * sxc.nspikes2);   
end
shufcorrm = mean(shufcorrt,1); % take mean of all shuffled correlations
covgmm = rawcorr - shufcorrm;
shufcorr_smm = smoothvect(shufcorrm, g1);
covgm_smm = smoothvect(covgmm, g1);  
    
% %Plot to check
figure; hold on; plot(xc.time, rawcorr_sm,'k-');  
plot(xc.time, shufcorr_smm,'b--'); plot(xc.time, covgm_smm, 'r-');

figure; hold on;
plot(xc.time, covgm_smm, 'r-');
plot(xc.time,covgm_sm(1,:),'r--');

% Get Nevents in raw correlation from -200ms to 200ms
bins = find(abs(xc.time)<=0.2);
Neventscorr = sum(xc.c1vsc2(bins));

% Can get PEAK and LATENCY TO PEAK OF COVGM AS WELL







% if (appendindex)
%     out.index = ind;
% end
out.index = ind;
out.corr = xc;
out.rawcorr = rawcorr;
out.shufcorr = shufcorr;
out.covgm = covgm;
out.rawcorr_sm = rawcorr_sm;
out.shufcorr_sm = shufcorrm_sm;
out.covgm_sm = covgm_sm;
out.Neventscorr = Neventscorr;



