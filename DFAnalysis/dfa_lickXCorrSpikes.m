

function out = dfa_lickXCorrSpikes(idx, excludeIntervals, varargin)

%{
run with filterframework via singleDayCellAnal
- this should really be done after epochs are combined... 
- gets spikes (burst) x lick (burst) xcorr and phase clustering
- example script: lickXcorrSU_20190916.m

idx: [day ntrode cluster eps:]
excludeIntervals: [start end; ...] timefilter (applied to spiketimes)
varargin required data: {'spikes', 'lick'}

- phase modulation score as used in karalis sirota 2018/19
- uses the circ stats toolbox's circ_rtest to get the pval and z which is (R^2/n)..
log(R^2/n) where R is the normed mean resultant vector and n is num samples
circ_r gets you the mean resultant vector, r
- note, karalis and sirota don't use any point process with less than 200 points in
the times of interest..
$DR19
%}

fprintf('%d %d %d %d %d\n', idx)
reqData = {'lick', 'spikes', 'cellinfo'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error('missing required data')
    end
end


bin = .01;
tmax = 1;
excShortBin = bin*2;
excLongBin = .250;
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
minBoutLicks = 3;
rmsmincounts = 1; % min bin count within rmstamax. otherwise nan
rmstmax = .25; % seconds
computeShuf = 1;
numShufs = 100;
maxShift = 500; % ms

% eventType = 'lick';
% psthbin = .001;
% TF = 1; % legacy tetfilter
% plotfigs = 1;
% displayplots = 0;
% saveplots = 1;
% savefigas = {'png', 'mfig'};
if ~isempty(varargin)
    assign(varargin{:});
end

% init output
out = init_out();
out.animal = animal;
out.index = idx;
out.bin = bin;
out.tmax = tmax;

day = idx(1);
eps = idx(4:5);
nt = idx(2);
clust = idx(3);
out.info = cellinfo{day}{eps(1)}{nt}{clust};

%% spiketimes
spikeTimes = [];
for e = eps
    try
        spikeTimes = [spikeTimes; spikes{day}{e}{nt}{clust}.data(:,1)];
    catch
        continue
    end
end

if isempty(spikeTimes)
    fprintf('spiketimes empty\n');
    return
end

includetimes = ~isExcluded(spikeTimes, excludeIntervals);
spikeTimes = spikeTimes(includetimes);
if isempty(spikeTimes)
    fprintf('spikeimes empty\n');
    return
end

%% Get lick-burst intervals
licksVecs = getLickBout([], animal, [day eps(1); day eps(2)], 'lick', lick, 'maxIntraBurstILI', ...
    maxILIthresh, 'minBoutLicks', minBoutLicks);
burstIntvs = [];
for e = eps
try
    burstIntvs = [burstIntvs; vec2list(licksVecs{day}{e}.lickBout, licksVecs{day}{e}.time)];
catch
    fprintf('error defining lick bouts for %d %d \n', day, ep)
    return
end
end
%% get licks within lick-burst intervals
lbLicks = [];
idlicks = [];
for e = eps
    licks = lick{day}{e}.starttime;
    idlicks = [idlicks; lick{day}{e}.id];
    lbLicks = [lbLicks; licks(logical(isExcluded(licks, burstIntvs)))]; % isIncluded
end
fprintf('%d licks within %d bursts \n',numel(lbLicks), size(burstIntvs,1))
%% time mod
%% run spike x lick xcorr
time = (-tmax+bin/2):bin:(tmax-bin/2);
xc = spikexcorr(spikeTimes, lbLicks, bin, tmax);
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
nstd = round(excShortBin/bin);
g1 = gaussian(nstd, 2*nstd+1);
smthxc = smoothvect(normxc, g1); % smooth xc
excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
    excLongBin));
% T = abs(diff(swr{day}{ep}{1}.timerange)); %total ep time
% xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
% p1 = xc.nspikes1/T;
% p2 = xc.nspikes2/T; % fr in Hz
% expProb = p1*p2; % per sec. Expected probability

out.time = time;
out.xc = xc;
out.normxc = normxc;
out.smthxc = smthxc;
out.excorr = excorr;
% out.xcrms = xcrms;
% out.p1 = p1;
% out.p2 = p2;
% out.expProb = expProb;

%% phasemod stuff
%% get spike events within lick-burst intervals
burstSpikeTime = spikeTimes(logical(isExcluded(spikeTimes, burstIntvs)));
burstSpikeTime = sort(burstSpikeTime(~isnan(burstSpikeTime)));
fprintf('spikes within %d lickbouts:  %d (%.02f pct) \n', size(burstIntvs,1), ...
    numel(burstSpikeTime), length(burstSpikeTime)/length(spikeTimes)*100)
if isempty(burstSpikeTime)
    fprintf('no spikes in lick bouts for %d %d %d %d %d skipping\n', idx)
    return
end

%% Get the containing licks for each spike
% histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin 
[~,~,spikeLickBinidx] = histcounts(burstSpikeTime, lbLicks);
% exclude swr's in 0 bin, before first lick edge
spikeLickBinidx = spikeLickBinidx(spikeLickBinidx > 0); 
burstSpikeTime = burstSpikeTime(spikeLickBinidx > 0);
% time since pre-containing lick
spikeTimeSinceLick = burstSpikeTime - lbLicks(spikeLickBinidx);
% get inter lick interval containing each spike
spikeBinILI = lbLicks(spikeLickBinidx+1) - lbLicks(spikeLickBinidx); % +1 is post-containing lick
% exclude swr-containing ili out of lick burst range 60 : 250 ms
spikeValidILI = all([spikeBinILI > minILIthresh spikeBinILI < maxILIthresh],2);
fprintf('lBurst spikes in valid ILI: %d (%.02f pct) \n', sum(spikeValidILI), ...
    sum(spikeValidILI)/length(spikeTimes)*100)
burstSpikeTime = burstSpikeTime(spikeValidILI);
if isempty(burstSpikeTime)
    return
end
spikeBinILI = spikeBinILI(spikeValidILI);
spikeTimeSinceLick = spikeTimeSinceLick(spikeValidILI);
% save spike-containing burst id
% [~,~,spikeLickBurstIdx] = histcounts(burstSpikeTime, [burstIntvs(:,1); inf]);
% burstContainSpike = burstIntvs(spikeLickBurstIdx,:);

%% phasemod
spikePctSinceLick = spikeTimeSinceLick ./ spikeBinILI;
spikeLickPhase = 2*pi*spikePctSinceLick; % pct of full cycle
meanvec = mean(exp(1i*spikeLickPhase)); % get mean resultant vector
meanMRVmag = abs(meanvec); % vector magnitude
vecang = angle(meanvec); % for plotting mrv
[~, z] = circ_rtest(spikeLickPhase); % z is mean res vec
phasemod = log(z); % log variance normalizes (karalis,sirota)

% out.spikeTimes = spikeTimes;
% out.licks = licks;
% out.idlicks = idlicks;
% out.boutIntvs = burstIntvs;
% out.swrTimeSinceLick = spikeTimeSinceLick;
% out.swrBinILI = spikeBinILI;
% out.burstSWRStart = burstSpikeTime;
% out.burstContainSpike = burstContainSpike;
% % out.dayEpoch = repmat([day ep], length(burstSpikeTime),1);
% out.spikePctSinceLick = spikePctSinceLick;
out.spikeLickPhase = spikeLickPhase;
out.meanMRVmag = meanMRVmag;
out.vecang = vecang;
out.phasemod = phasemod;


%% Compute shuffled lickDin x swr
if computeShuf
    tic
    for i = 1:numShufs
        % shuffle spike times locally 
        r = randi([-maxShift maxShift],length(spikeTimes(:,1)),1)/1e3;
        spikeTimesShuf = sort(spikeTimes+r);
       %% time mod
       %% run spike x lick xcorr
        xc = spikexcorr(spikeTimesShuf, lbLicks, bin, tmax);
        normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
        try
            smthxc = smoothvect(normxc, g1);
        catch % if normxc is empty
            continue
        end
        excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
            excLongBin));
%         xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
%         p1 = xc.nspikes1/T;
%         p2 = xc.nspikes2/T; % fr in Hz
%         expProb = p1*p2; % per sec. Expected probability
        
        out.xcShuf{end+1} = xc;
        out.normxcShuf{end+1} = normxc;
        out.smthxcShuf{end+1} = smthxc;
        out.excorrShuf{end+1} = excorr;
%         out.xcrmsShuf{end+1} = xcrms;
%         out.p1Shuf{end+1} = p1;
%         out.p2Shuf{end+1} = p2;
%         out.expProbShuf{end+1} = expProb;

       %% phasemod stuff
       %% get spike events within lick-burst intervals
        burstSpikeTime = spikeTimesShuf(logical(isExcluded(spikeTimesShuf, burstIntvs)));
        burstSpikeTime = sort(burstSpikeTime(~isnan(burstSpikeTime)));
%         fprintf('spikes within %d lickbouts:  %d (%.02f pct swrs) \n', size(burstIntvs,1), ...
%             numel(burstSpikeTime), length(burstSpikeTime)/length(spikeTimesShuf)*100)
        if isempty(burstSpikeTime)
            continue
        end
        
        % histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin 
        [~,~,spikeLickBinidx] = histcounts(burstSpikeTime, lbLicks);
        % exclude spikes in 0 bin, before first lick edge
        spikeLickBinidx = spikeLickBinidx(spikeLickBinidx > 0);
        burstSpikeTime = burstSpikeTime(spikeLickBinidx > 0);
        % time since pre-containing lick
        spikeTimeSinceLick = burstSpikeTime - lbLicks(spikeLickBinidx);
        % get inter lick interval containing each spike
        spikeBinILI = lbLicks(spikeLickBinidx+1) - lbLicks(spikeLickBinidx);
        % exclude spike-containing ili out of lick burst range 60 : 250 ms
        spikeValidILI = all([spikeBinILI > minILIthresh spikeBinILI < maxILIthresh],2);
        burstSpikeTime = burstSpikeTime(spikeValidILI);
        spikeBinILI = spikeBinILI(spikeValidILI);
        spikeTimeSinceLick = spikeTimeSinceLick(spikeValidILI);
%         [~,~,spikeLickBurstIdx] = histcounts(burstSpikeTime, [burstIntvs(:,1); inf]);
%         burstContainSpike = burstIntvs(spikeLickBurstIdx,:);
        if isempty(spikeTimeSinceLick)
            return
        end
       % phasemod
        spikePctSinceLick = spikeTimeSinceLick ./ spikeBinILI;
        spikeLickPhase = 2*pi*(spikePctSinceLick);
        meanvec = mean(exp(1i*spikeLickPhase));
        meanMRVmag = abs(meanvec);
        vecang = angle(meanvec);
        [~, z] = circ_rtest(spikeLickPhase); % pval is the stat rayleigh test. z is mean res vec
        phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)
        
%         out.burstContainSpikeShuf{end+1} = burstContainSpike;
%         out.burstSpikeTimeShuf{end+1} = burstSpikeTime;
%         out.swrTimeSinceLickShuf{end+1} = spikeTimeSinceLick;
%         out.spikePctSinceLickShuf{end+1} = spikePctSinceLick;
        out.spikeLickPhaseShuf{end+1} = spikeLickPhase;
        out.meanMRVmagShuf{end+1} = meanMRVmag;
        out.vecangShuf{end+1} = vecang;
        out.phasemodShuf{end+1} = phasemod;
%         
%         out.swrBinILIShuf{end+1} = spikeBinILI;
%         out.swrInBurstStartShuf{end+1} = burstSpikeTime;
%         out.swrTimeSinceLickShuf{end+1} = spikeTimeSinceLick;
%         out.swrBurstIntervalShuf{end+1} = burstContainSwr;
        
    end
    fprintf('shuffle took %.02f s\n', toc);
end
% varargin{end+1} = 'byDay';
% varargin{end+1} = 1;
% out.evTrigSpike = dfa_eventTrigSpiking(idx, excludeIntervals, 'byDay', 1, 'eventType', ...
%     eventType, 'lick', lick, 'spikes', spikes);
end

function out = init_out()
out.index = [];
out.animal = [];
out.info = [];
out.time = [];
% out.boutIntvs = [];
out.xc = [];
out.normxc = [];
out.smthxc = [];
out.excorr = [];

% iliPhase
% out.burstContainSpike = [];
% out.spikePctSinceLick = [];
out.spikeLickPhase = [];
out.meanMRVmag = [];
out.vecang = [];
out.phasemod = [];

% shuffle
% xcTime
out.xcShuf = {};
out.normxcShuf = {};
out.smthxcShuf = {};
out.excorrShuf = {};
% iliPhase
out.spikeLickPhaseShuf = {};
out.meanMRVmagShuf = {};
out.vecangShuf = {};
out.phasemodShuf = {};

% evTrigSpike
out.evTrigSpike = [];
end
