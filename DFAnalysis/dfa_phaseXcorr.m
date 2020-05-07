function out = dfa_phaseXcorr(idx, timeFilter, varargin)
% filterframework analysis function 
% return pairwise phase cross-correlation
%
% Input:
%   idx : [day ntA clustA ntB clustB ep1..epN]
%
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.

bin = 0.01; % seconds
sw1 = bin*3; % seconds
sw2 = 0.250;
rmstmax = 0.1;
rmsmincounts = 50;
tmax = 1;
try assign(varargin{:}); catch; end

% check for required data in varargin
reqData = {'spikes', 'lick'};
for s = 1:length(reqData)
    if ~any(cell2mat(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 0)))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

day = idx(1);
ntA = idx(2);
clustA = idx(3);
ntB = idx(4);
clustB = idx(5);
eps = idx(6:end);

out = init_out();
out.index = idx;

% for each cell we calculate the cross correlation 
t1 = [];
t2 = [];
totalEpsDur = 0;
for iep = 1:length(eps)
    totalEpsDur = totalEpsDur + ...
        diff(spikes{day}{eps(iep)}{ntA}{clustA}.timerange)./10000; % in secs
    try
        t1 = [t1; spikes{day}{eps(iep)}{ntA}{clustA}.data];
        t2 = [t2; spikes{day}{eps(iep)}{ntB}{clustB}.data];
    catch
        warning('couldnt load spikes %s d%d [ntA clA ntB clB]:: %d %d %d %d\n', ...
            animal, day, nyA, clustA, ntB, clustB);
        return
        % if either of those produced an error, we return NaNs 
    end
end
%apply the exclude rule
t1inc = [];
t2inc = [];
timeFilterCat = cell2mat([{timeFilter{:}}']); % combine epoch timefilters
t1cnd = t1(~logical(isExcluded(t1, timeFilterCat)));
t2cnd = t2(~logical(isExcluded(t2, timeFilterCat)));

% Get total time for included data
totalExclDur = nansum(diff(timeFilterCat,[],2));
T = totalEpsDur-totalExclDur;

%% get lickbout licks
[intraBoutXP, ~] = getLickBoutLicks(animal, ...
    [repmat(day,length(eps),1) eps'], varargin);
intraBoutXPvec = intraBoutXP{day}{eps};
licks = cell2mat([arrayfun(@(x) lick{day}{x}.starttime, eps, 'un', 0)']);
eventTimes = intraBoutXPvec(ismember(intraBoutXPvec, licks));

numEventsPerEp = [];
for iep = 1:length(eps)
    epStartTime = spikes{day}{eps(iep)}{ntA}{clustA}.timerange(1);
    epEndTime = spikes{day}{eps(iep)}{ntA}{clustA}.timerange(2);
    numEventsPerEp = [numEventsPerEp; sum(logical(isExcluded(eventTimes, ...
        [epStartTime epEndTime])))];
end

%% transform spikes into a phase vector
% how? each spike is a spike time.. i need each spike time relative to
% preceeding lick then ratio to interlickinterval duration to get %.. then
% radians.. then that + (2pi*[preceding lick number])

[~, ~, spikebin] = histcounts(t1cnd, eventTimes);
spikesILImask = spikebin~=0;
t1LB = t1cnd(spikesILImask);
binTimeStart = eventTimes(spikebin(spikesILImask));
binTimeEnd = eventTimes(spikebin(spikesILImask)+1);

ILI = (binTimeEnd - binTimeStart);
t1FromLick = (t1LB - binTimeStart);

pctSpikeILI = t1FromLick./ILI;

% for each valid ILI, increment its spikes to be pct of ILI + the ILI index
t1SpikePctILIbin = spikebin(spikesILImask)+pctSpikeILI;

% spiketrain B
[~, ~, spikebin] = histcounts(t2cnd, eventTimes);
spikesILImask = spikebin~=0;
t2LB = t2cnd(spikesILImask);
binTimeStart = eventTimes(spikebin(spikesILImask));
binTimeEnd = eventTimes(spikebin(spikesILImask)+1);

ILI = (binTimeEnd - binTimeStart);
t2FromLick = (t2LB - binTimeStart);

pctSpikeILI = t2FromLick./ILI;

% for each valid ILI, increment its spikes to be pct of ILI + the ILI index
t2SpikePctILIbin = spikebin(spikesILImask)+pctSpikeILI;

%% Phase xcorr

% get xc
xc = spikexcorr(t1SpikePctILIbin,t2SpikePctILIbin, bin, tmax);
% compute the excess correlation at 0 lag
exc = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
% compute RMS
xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);

% % if we want to include edge spikes, we need to add in the correlation of the
% % excluded t1 spikes with the included t2spikes
% if ((edgespikes) & (~isempty(xc.time)))
%     t1ex = t1(find(isExcluded(t1, excludetimes)));
%     if (~isempty(t1ex))
%         tmpxc = spikexcorr(t1ex, t2inc, bin, tmax);
%         % add these values to the original histogram
%         xc.c1vsc2 = xc.c1vsc2 + tmpxc.c1vsc2;
%     end
% end

% Expected probability
p1 = xc.nspikes1/T; p2 = xc.nspikes2/T; % Fir rate in Hz
exp_p = p1*p2; % per sec

% Crosscov
crosscov = (xc.c1vsc2 ./ (bin*T))-exp_p;
% Convert to Z-score
factor = sqrt((bin*T) / exp_p);
Zcrosscov = crosscov .* (factor);

% normxc:: Normalize by geometric mean: Units will be coincidences/ spike
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);

% smooth
try
    nstd=round(sw1/(xc.time(2) - xc.time(1)));
    g1 = gaussian(nstd, 2*nstd+1); %g1 = gaussian(nstd, 2*nstd+1);
    % Can also do a 3 point gaussian with sigma=3. Almost a boxcar.
    %     g1 = gaussian(nstd, nstd);
catch
    fprintf('missing spikes \n')
    return
end


% output
out.normxc = normxc;

out.crosscov = crosscov;
out.Zcrosscov = Zcrosscov;

out.corr_sm = smoothvect(xc.c1vsc2, g1);
out.normxc_sm = smoothvect(normxc, g1);
out.Zcrosscov_sm = smoothvect(Zcrosscov, g1);
out.crosscov_sm = smoothvect(crosscov, g1);

out.T = T;
out.bin = bin;
out.tmax = tmax;
out.xc = xc;
out.normxc = normxc;
% out.expProb = expProb;
% out.xcShfmean = xcShfmean;
% out.xcShfLag0 = xcShfLag0;
out.excesscorr = nanmean(exc); % nanmean in case there are two bins equally near lag zero
out.xcrms = xcrms;

end

function out = init_out()
out.corr_sm = NaN;

out.Zcrosscov = NaN;
out.crosscov = NaN;
out.crosscov_sm = NaN;
out.Zcrosscov_sm = NaN;
out.T = NaN;
out.tmax = nan;
out.bin = nan;
out.xc = nan;
out.normxc = nan;
out.normxc_sm = nan;
% out.expProb = expProb;
out.excesscorr = nan;
out.xcrms = nan;
end