

function out = dfa_lickswrcorr(idx, excludeIntervals, varargin)
%{
run with filterframework via singleepochanal
- gets swr (burst) x lick (burst) xcorr and phase clustering
- plot epoch result (move this out to script).

idx: [day epoch]
excludeIntervals: [start end; ...] timefilter
varargin required data: {'ca1rippleskons','task', 'lick', 'DIO'}

- phase modulation score as used in karalis sirota 2018/19
- uses the circ stats toolbox's circ_rtest to get the pval and z which is (R^2/n)..
log(R^2/n) where R is the normed mean resultant vector and n is num samples
circ_r gets you the mean resultant vector, r
- note, karalis and sirota don't use any point process with less than 200 points in
the times of interest..
- DR 19
%}

fprintf('%d %d\n',idx)
reqData = {'ca1rippleskons','task', 'lick', 'DIO'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

pconf = paramconfig;
bin = .01;
tmax = 1;
excShortBin = bin*2;
excLongBin = .250;
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
minBoutLicks = 3;
rmsmincounts = 1; % min bin count within rmstamax. otherwise nan
rmstmax = .25; % seconds
compute_shuffle = 1;
numshuffs = 100;
shuffOffset = 250; %ms
% sigpct = .975;
eventType = 'ca1rippleskons';
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1,1);
ep = idx(1,2);
out = init_out(); % init output
out.index = idx;
out.animal = animal;

%% Get SWR events
evid = find(contains(varargin(1:2:end), eventType));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};
eventTimes = [];
try
    eventTimes = events{day}{ep}{1}.starttime;
catch
    try
        eventTimes = events{day}{ep}.starttime;
    catch
        fprintf('no events detected for day%d ep%d\n', day,ep)
        return
    end
end
ecbefore = size(eventTimes,1);
eventTimes = eventTimes(~isExcluded(eventTimes(:,1),excludeIntervals),:);
ecafter = size(eventTimes,1);
fprintf('%d of %d events discarded bc excluded periods in timefilter: d%d e%d\n',...
    ecbefore-ecafter, ecbefore, day,ep)
if isempty(eventTimes)
    return
end

%% Get SWR, Lick events within lick-burst intervals
licksVecs = getLickBout([], animal, [day ep], 'lick', lick, 'maxIntraBurstILI', ...
    maxILIthresh, 'minBoutLicks', minBoutLicks);
try
    burstIntvs = vec2list(licksVecs{day}{ep}.lickBout, licksVecs{day}{ep}.time);
catch
    fprintf('error defining lick bouts for %d %d \n', day, ep)
    return
end
% filter for swr within lick-bursts (this should already be done with
% timefilter.. but whatevs)
swrInBurstStart = eventTimes(logical(isExcluded(eventTimes, burstIntvs)));
swrInBurstStart = sort(swrInBurstStart(~isnan(swrInBurstStart)));
fprintf('swrs within %d lickbouts:  %d (%.02f pct swrs) \n', size(burstIntvs,1), ...
    numel(swrInBurstStart), length(swrInBurstStart)/length(eventTimes)*100)
if isempty(swrInBurstStart)
    fprintf('no swrs in lick bouts for %d %d.. skipping\n', day, ep)
    return
end
% filter for licks within lick-bursts
licks = lick{day}{ep}.starttime;
idlicks = lick{day}{ep}.id;
lbLicks = licks(logical(isExcluded(licks, burstIntvs)));
fprintf('licks within %d bursts: %d (%.02f pct) \n',size(burstIntvs,1), ...
    numel(lbLicks), length(lbLicks)/length(licks)*100)

[~,~,swrLickBinidx] = histcounts(swrInBurstStart, lbLicks);
% exclude swr's in 0 bin. how does this happen? maybe if same time as first lick
swrLickBinidx = swrLickBinidx(swrLickBinidx > 0);
swrInBurstStart = swrInBurstStart(swrLickBinidx > 0);

%% time since last lick within valid ILI
swrTimeSinceLick = swrInBurstStart - lbLicks(swrLickBinidx);
if any(swrTimeSinceLick) < 0
    error('swr rel time must be pos')
end
% get inter lick interval for each swr
swrBinILI = lbLicks(swrLickBinidx+1) - lbLicks(swrLickBinidx);
% exclude ili out of lick burst range
swrValidILI = all([swrBinILI > minILIthresh swrBinILI < maxILIthresh],2);
fprintf('lBurst swrs in valid ILI: %d (%.02f pct swrs) \n', sum(swrValidILI), ...
    sum(swrValidILI)/length(eventTimes)*100)
swrInBurstStart = swrInBurstStart(swrValidILI);
if isempty(swrInBurstStart)
    return
end
swrBinILI = swrBinILI(swrValidILI);
swrTimeSinceLick = swrTimeSinceLick(swrValidILI);
% save burst start/end time for each enclosed and valid swr
[~,~,swrLickBurstIdx] = histcounts(swrInBurstStart, [burstIntvs(:,1); inf]);
swrBurstInterval = burstIntvs(swrLickBurstIdx,:);
out.swrStart = eventTimes;
out.licks = licks;
out.idlicks = idlicks;
out.boutIntvs = burstIntvs;
out.swrBinILI = swrBinILI;
out.swrInBurstStart = swrInBurstStart;
out.swrTimeSinceLick = swrTimeSinceLick;
out.swrBurstInterval = swrBurstInterval;
out.dayEpoch = repmat([day ep], length(swrInBurstStart),1);
%% xcorr
time = (-tmax+bin/2):bin:(tmax-bin/2);
nstd = round(excShortBin/bin);
g1 = gaussian(nstd, 2*nstd+1);
xc = spikexcorr(swrInBurstStart, lbLicks, bin, tmax);
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
smthxc = smoothvect(normxc, g1); % smooth xc
excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
    excLongBin));
T = abs(diff(events{day}{ep}{1}.timerange)); %total ep time
xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
p1 = xc.nspikes1/T;
p2 = xc.nspikes2/T; % fr in Hz
expProb = p1*p2; % per sec. Expected probability

out.time = time;
out.xc = xc;
out.normxc = normxc;
out.smthxc = smthxc;
out.excorr = excorr;
out.xcrms = xcrms;
out.p1 = p1;
out.p2 = p2;
out.expProb = expProb;

%% phasemod
swrPctSinceLick = swrTimeSinceLick ./ swrBinILI;
swrLickPhase = 2*pi*swrPctSinceLick; % pct of full cycle.. cart to polar
meanvec = mean(exp(1i*swrLickPhase));
meanMRVmag = abs(meanvec);
vecang = angle(meanvec);
[pval, z] = circ_rtest(swrLickPhase); % pval is the stat rayleigh test. z is mean res vec
phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)

out.vecang = vecang;
out.swrLickPhase = swrLickPhase;
out.swrPctSinceLick = swrPctSinceLick;
out.phasemod = phasemod;
out.meanMRVmag = meanMRVmag;

%% Compute shuffled lickDin x swr
if compute_shuffle
    tic
    meanMRVmag = [];
    smthxc = [];
    excorr = [];
    phasemod = [];
    % shuffle swr start time
    r = randi([-shuffOffset shuffOffset],length(eventTimes(:,1)),numshuffs)/1e3;
    swrStartShuf = bsxfun(@plus,eventTimes(:,1),r); %, make shuffle array
    out.swrStartShuf = swrStartShuf;
    %%
    for i = 1:numshuffs
        swrInBurstStart = swrStartShuf(logical(isExcluded(swrStartShuf(:,i), burstIntvs)),i);
        swrInBurstStart = sort(swrInBurstStart(~isnan(swrInBurstStart)));
        if isempty(swrInBurstStart)
            continue
        end
        [~,~,swrLickBinidx] = histcounts(swrInBurstStart, lbLicks);
        swrLickBinidx = swrLickBinidx(swrLickBinidx > 0);
        swrInBurstStart = swrInBurstStart(swrLickBinidx > 0);
        swrTimeSinceLick = swrInBurstStart - lbLicks(swrLickBinidx);
        if any(swrTimeSinceLick) < 0
            error('swr rel time must be pos')
        end
        swrBinILI = lbLicks(swrLickBinidx+1) - lbLicks(swrLickBinidx);
        swrValidILI = all([swrBinILI > minILIthresh swrBinILI < maxILIthresh],2);
        swrInBurstStart = swrInBurstStart(swrValidILI);
        swrBinILI = swrBinILI(swrValidILI);
        swrTimeSinceLick = swrTimeSinceLick(swrValidILI);
        [~,~,swrLickBurstIdx] = histcounts(swrInBurstStart, [burstIntvs(:,1); inf]);
        swrBurstInterval = burstIntvs(swrLickBurstIdx,:);
        out.swrBinILIShuf{end+1} = swrBinILI;
        out.swrInBurstStartShuf{end+1} = swrInBurstStart;
        out.swrTimeSinceLickShuf{end+1} = swrTimeSinceLick;
        out.swrBurstIntervalShuf{end+1} = swrBurstInterval;
        
        %% xcorr
        xc = spikexcorr(sort(swrInBurstStart), lbLicks, bin, tmax);
        normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
        try
            smthxc = smoothvect(normxc, g1);
        catch % if normxc is empty
            continue
        end
        excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
            excLongBin));
        xcrms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts); % RMS
        p1 = xc.nspikes1/T;
        p2 = xc.nspikes2/T; % fr in Hz
        expProb = p1*p2; % per sec. Expected probability
        
        out.xcShuf{end+1} = xc;
        out.normxcShuf{end+1} = normxc;
        out.smthxcShuf{end+1} = smthxc;
        out.excorrShuf{end+1} = excorr;
        out.xcrmsShuf{end+1} = xcrms;
        out.p1Shuf{end+1} = p1;
        out.p2Shuf{end+1} = p2;
        out.expProbShuf{end+1} = expProb;
        
        %% phasemod
        swrPctSinceLick = swrTimeSinceLick ./ swrBinILI;
        swrLickPhase = 2*pi*(swrPctSinceLick);
        meanvec = mean(exp(1i*swrLickPhase));
        meanMRVmag = abs(meanvec);
        vecang = angle(meanvec);
        [~, z] = circ_rtest(swrLickPhase); % pval is the stat rayleigh test. z is mean res vec
        phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)
        
        out.swrBurstIntervalShuf{end+1} = swrBurstInterval;
        out.swrInBurstStartShuf{end+1} = swrInBurstStart;
        out.swrTimeSinceLickShuf{end+1} = swrTimeSinceLick;
        out.meanMRVmagShuf{end+1} = meanMRVmag;
        out.swrLickPhaseShuf{end+1} = swrLickPhase;
        out.swrPctSinceLickShuf{end+1} = swrPctSinceLick;
        out.vecangShuf{end+1} = vecang;
        out.phasemodShuf{end+1} = phasemod;
    end
    fprintf('shuffle took %.02f s\n', toc);
end
end

function out = init_out()
out.index = [];
out.animal = [];
out.time = [];
out.licks = [];
out.idlicks = [];
out.boutIntvs = [];
out.swrStart = [];
out.dayEpoch = [];
% xcTime
out.xc = [];
out.normxc = [];
out.smthxc = [];
out.excorr = [];
out.xcrms = [];
out.p1 = [];
out.p2 = [];
out.expProb = [];
% iliPhase
out.swrBinILI = [];
out.swrBurstInterval = [];
out.swrInBurstStart = [];
out.swrTimeSinceLick = [];
out.meanMRVmag = [];
out.swrLickPhase = [];
out.swrPctSinceLick = [];
out.vecang = [];
out.phasemod = [];

%% shuffle
out.swrStartShuf = [];
% xcTime
out.xcShuf = {};
out.normxcShuf = {};
out.smthxcShuf = {};
out.excorrShuf = {};
out.xcrmsShuf = {};
out.p1Shuf = {};
out.p2Shuf = {};
out.expProbShuf = {};
% iliPhase
out.swrBinILIShuf = {};
out.swrBurstIntervalShuf = {};
out.swrInBurstStartShuf = {};
out.swrTimeSinceLickShuf = {};
out.meanMRVmagShuf = {};
out.swrLickPhaseShuf = {};
out.swrPctSinceLickShuf = {};
out.vecangShuf = {};
out.phasemodShuf = {};

end