

function out = dfa_lickswrcorr(idx, excludeIntervals, varargin)

%{
run with filterframework via singleepochanal
- gets swr (burst) x lick (burst) xcorr and phase clustering
- example script: swrlickxcorr_20191031.m

idx: [day epoch]
excludeIntervals: [start end; ...] timefilter (applied to recover swrs)
varargin required data: {'ca1rippleskons', 'lick'}

- phase modulation score as used in karalis sirota 2018/19
- uses the circ stats toolbox's circ_rtest to get the pval and z which is (R^2/n)..
log(R^2/n) where R is the normed mean resultant vector and n is num samples
circ_r gets you the mean resultant vector, r
- note, karalis and sirota don't use any point process with less than 200 points in
the times of interest..
- DR 19
%}

fprintf('%d %d\n',idx)
reqData = {'ca1rippleskons', 'lick'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
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
compute_shuffle = 1;
numshuffs = 100;
maxShift = 500; %ms
eventType = 'ca1rippleskons';
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1,1);
ep = idx(1,2);
out = init_out(); % init output
out.index = idx;
out.animal = animal;

%% Recover SWR events

evid = find(contains(varargin(1:2:end), eventType));
o = [1:2:length(varargin)]+1;
swr = varargin{o(evid)};
swrTimes = [];
try
    swrTimes = swr{day}{ep}{1}.starttime;
catch
    try
        swrTimes = swr{day}{ep}.starttime;
    catch
        fprintf('no swr events detected for day%d ep%d\n', day, ep)
        return
    end
end
ecbefore = size(swrTimes,1);
swrTimes = swrTimes(~isExcluded(swrTimes(:,1), excludeIntervals),:); % isIncluded
ecafter = size(swrTimes,1);
fprintf('%d of %d swr events discarded bc excluded periods in timefilter: d%d e%d\n',...
    ecbefore-ecafter, ecbefore, day,ep)
if isempty(swrTimes)
    return
end

%% Get lick-burst intervals
licksVecs = getLickBout([], animal, [day ep], 'lick', lick, 'maxIntraBurstILI', ...
    maxILIthresh, 'minBoutLicks', minBoutLicks);
try
    burstIntvs = vec2list(licksVecs{day}{ep}.lickBout, licksVecs{day}{ep}.time);
catch
    fprintf('error defining lick bouts for %d %d \n', day, ep)
    return
end

%% get licks within lick-burst intervals
licks = lick{day}{ep}.starttime;
idlicks = lick{day}{ep}.id;
lbLicks = licks(logical(isExcluded(licks, burstIntvs))); % isIncluded
fprintf('licks within %d bursts: %d (%.02f pct) \n',size(burstIntvs,1), ...
    numel(lbLicks), length(lbLicks)/length(licks)*100)

%% calc swr vs lick xcorr
time = (-tmax+bin/2):bin:(tmax-bin/2);
xc = spikexcorr(swrTimes, lbLicks, bin, tmax);
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
nstd = round(excShortBin/bin);
g1 = gaussian(nstd, 2*nstd+1);
smthxc = smoothvect(normxc, g1); % smooth xc
excorr = nanmean(excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, excShortBin, ...
    excLongBin));
T = abs(diff(swr{day}{ep}{1}.timerange)); %total ep time
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

%% phasemod stuff
%% get swr events within lick-burst intervals
burstSWRStart = swrTimes(logical(isExcluded(swrTimes, burstIntvs)));
burstSWRStart = sort(burstSWRStart(~isnan(burstSWRStart)));
fprintf('swrs within %d lickbouts:  %d (%.02f pct swrs) \n', size(burstIntvs,1), ...
    numel(burstSWRStart), length(burstSWRStart)/length(swrTimes)*100)
if isempty(burstSWRStart)
    fprintf('no swrs in lick bouts for %d %d.. skipping\n', day, ep)
    return
end

%% Get the containing licks for each swr
% histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin 
[~,~,swrLickBinidx] = histcounts(burstSWRStart, lbLicks);
% exclude swr's in 0 bin, before first lick edge
swrLickBinidx = swrLickBinidx(swrLickBinidx > 0); 
burstSWRStart = burstSWRStart(swrLickBinidx > 0);

% time since pre-containing lick
swrTimeSinceLick = burstSWRStart - lbLicks(swrLickBinidx);

% get inter lick interval containing each swr
swrBinILI = lbLicks(swrLickBinidx+1) - lbLicks(swrLickBinidx); % +1 is post-containing lick

% exclude swr-containing ili out of lick burst range 60 : 250 ms
swrValidILI = all([swrBinILI > minILIthresh swrBinILI < maxILIthresh],2);
fprintf('lBurst swrs in valid ILI: %d (%.02f pct swrs) \n', sum(swrValidILI), ...
    sum(swrValidILI)/length(swrTimes)*100)
burstSWRStart = burstSWRStart(swrValidILI);
if isempty(burstSWRStart)
    return
end
swrBinILI = swrBinILI(swrValidILI);
swrTimeSinceLick = swrTimeSinceLick(swrValidILI);

% save swr-containing burst id
[~,~,swrLickBurstIdx] = histcounts(burstSWRStart, [burstIntvs(:,1); inf]);
burstContainSwr = burstIntvs(swrLickBurstIdx,:);

% phasemod
swrPctSinceLick = swrTimeSinceLick ./ swrBinILI;
swrLickPhase = 2*pi*swrPctSinceLick; % pct of full cycle
meanvec = mean(exp(1i*swrLickPhase)); % get mean resultant vector
meanMRVmag = abs(meanvec); % vector magnitude
vecang = angle(meanvec);
Z = meanMRVmag^2/n;
[~, z] = circ_rtest(swrLickPhase); % z is mean res vec
if Z ~= z
    error('-----wut?')
end
phasemod = log(z); % log variance normalizes (karalis,sirota)

out.swrStart = swrTimes;
out.licks = licks;
out.idlicks = idlicks;
out.boutIntvs = burstIntvs;
out.swrTimeSinceLick = swrTimeSinceLick;
out.swrBinILI = swrBinILI;
out.burstSWRStart = burstSWRStart;
out.burstContainSwr = burstContainSwr;
out.dayEpoch = repmat([day ep], length(burstSWRStart),1);
out.swrPctSinceLick = swrPctSinceLick;
out.swrLickPhase = swrLickPhase;
out.meanMRVmag = meanMRVmag;
out.vecang = vecang; % for polar scatter plots
out.phasemod = phasemod;

%% Compute shuffled lickDin x swr
if compute_shuffle
    tic
    meanMRVmag = [];
    smthxc = [];
    excorr = [];
    phasemod = [];
    % shuffle swr start time
    r = randi([-maxShift maxShift],length(swrTimes(:,1)),numshuffs)/1e3;
    swrStartShuf = sort(bsxfun(@plus,swrTimes(:,1),r)); %, make shuffle array
%     out.swrStartShuf = swrStartShuf;
    %%
    for i = 1:numshuffs
                %% xcorr
        xc = spikexcorr(swrStartShuf(:,i), lbLicks, bin, tmax);
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
        
        %% phasemod stuff
       %% get spike events within lick-burst intervals
        burstSWRStart = swrStartShuf(logical(isExcluded(swrStartShuf(:,i), burstIntvs)),i);
        burstSWRStart = burstSWRStart(~isnan(burstSWRStart));
        if isempty(burstSWRStart)
            continue
        end
        
        % histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin 
        [~,~,swrLickBinidx] = histcounts(burstSWRStart, lbLicks);
        % exclude swrs in 0 bin, before first lick edge
        swrLickBinidx = swrLickBinidx(swrLickBinidx > 0);
        burstSWRStart = burstSWRStart(swrLickBinidx > 0);
        % time since pre-containing lick
        swrTimeSinceLick = burstSWRStart - lbLicks(swrLickBinidx);
        % get inter lick interval containing each swr
        swrBinILI = lbLicks(swrLickBinidx+1) - lbLicks(swrLickBinidx);
        % exclude swr-containing ili out of lick burst range 60 : 250 ms
        swrValidILI = all([swrBinILI > minILIthresh swrBinILI < maxILIthresh],2);
        burstSWRStart = burstSWRStart(swrValidILI);
        swrBinILI = swrBinILI(swrValidILI);
        swrTimeSinceLick = swrTimeSinceLick(swrValidILI);
        [~,~,swrLickBurstIdx] = histcounts(burstSWRStart, [burstIntvs(:,1); inf]);
        burstContainSwr = burstIntvs(swrLickBurstIdx,:);
        
        % phasemod
        swrPctSinceLick = swrTimeSinceLick ./ swrBinILI;
        swrLickPhase = 2*pi*(swrPctSinceLick);
        meanvec = mean(exp(1i*swrLickPhase));
        meanMRVmag = abs(meanvec);
        vecang = angle(meanvec);
        [~, z] = circ_rtest(swrLickPhase); % pval is the stat rayleigh test. z is mean res vec
        phasemod = log(z); % i think log makes it 'variance normalized' (karalis,sirota)
        
%         out.swrBinILIShuf{end+1} = swrBinILI;
%         out.swrInBurstStartShuf{end+1} = burstSWRStart;
%         out.swrTimeSinceLickShuf{end+1} = swrTimeSinceLick;
%         out.swrBurstIntervalShuf{end+1} = burstContainSwr;

%         out.swrBurstIntervalShuf{end+1} = burstContainSwr;
%         out.swrInBurstStartShuf{end+1} = burstSWRStart;
%         out.swrTimeSinceLickShuf{end+1} = swrTimeSinceLick;
%         out.swrPctSinceLickShuf{end+1} = swrPctSinceLick;
        out.swrLickPhaseShuf{end+1} = swrLickPhase;
        out.meanMRVmagShuf{end+1} = meanMRVmag;
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
% out.swrStartShuf = [];
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