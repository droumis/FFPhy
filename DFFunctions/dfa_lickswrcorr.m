

function out = dfa_lickswrcorr(idx, timeFilter, varargin)
% out = dfa_lickswrcorr(idx, timeFilter, varargin)
% gets swr x xp xcorr and phase clustering
%
% Iterator:
% - singleDayAnal (city)
% - no longer: singleepochanal (space)
% 
% args:
% - idx [day epochs]
% - timeFilter: applied to swr
%               .. not applied to licks because they need a different
%               timefilter.. currently 'within XPburst intervals'
% 
% varargs:
% - data (required) {'ca1rippleskons', 'lick'}
% - win:
% - eventType:
% - LFPTypes: 
%{

Notes:
- example script: swrlickxcorr_20191031.m
- phase modulation score as used in karalis sirota 2018/19
- uses the circ stats toolbox's circ_rtest to get the pval and z which is (R^2/n)..
log(R^2/n) where R is the normed mean resultant vector and n is num samples
circ_r gets you the mean resultant vector, r
- note, karalis and sirota don't use any point process with less than 200 points in
the times of interest..

%}

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'lick'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end
byDay = 1;
bin = .01;
tmax = 1;
eventType = 'ca1rippleskons';
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
minBoutLicks = 3;
% xcorr / excorr
excShortBin = bin*2;
excLongBin = .250;
rmsmincounts = 1; % min bin count within rmstamax. otherwise nan
rmstmax = .25; % seconds
% shuf
compute_shuffle = 1;
numshuffs = 100;
maxShift = 500; %ms
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1);
if byDay
    eps = idx(2:end);
else
    eps = idx(2);
end

out = init_out(); % init output
out.index = idx;
out.animal = animal;

%% Get SWRs in timeFilter

evid = find(contains(varargin(1:2:end), eventType));
o = [1:2:length(varargin)]+1;
swr = varargin{o(evid)};
swrTimes = [];
maxthresh = [];
swrEnd = [];
for e = 1:length(eps)
    ep = eps(e);
    try
        swrTimes = [swrTimes; swr{day}{ep}{1}.starttime];
        maxthresh = [maxthresh; swr{day}{ep}{1}.maxthresh];
        swrEnd = [swrEnd; swr{day}{ep}{1}.endtime];
    catch
        try
            swrTimes = [swrTimes; swr{day}{ep}.starttime];
            maxthresh = [maxthresh; swr{day}{ep}.maxthresh];
            swrEnd = [swrEnd; swr{day}{ep}.endtime];
        catch
            fprintf('no swr events detected for day%d ep%d\n', day, ep)
            return
        end
    end
end
% apply timefilter to swrs
incSWRs = ~isIncluded(swrTimes(:,1), timeFilter);
fprintf('%d of %d swr events discarded bc excluded periods in timefilter: d%d\n',...
    sum(incSWRs), length(incSWRs), day)
swrTimes = swrTimes(incSWRs,:);
out.maxthresh = maxthresh;
out.swrEnd = swrEnd;

if isempty(swrTimes)
    return
end

%% Get licks in lick-burst intervals
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:});
boutTimes = cell2mat({boutTimes{day}{eps}}');
intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
intraBoutXP = intraBoutXP(:,1);
fprintf('%d XP within %d bursts \n', numel(intraBoutXP), size(boutTimes,1))

% intraBoutXP = lick{day}{ep}.starttime;
% idlicks = lick{day}{ep}.id;
% lbLicks = licks(logical(isExcluded(licks, burstIntvs))); % isIncluded
% keep valid lickburst licks
%
%% calc swr vs lick xcorr
time = (-tmax+bin/2):bin:(tmax-bin/2);
xc = spikexcorr(swrTimes, intraBoutXP, bin, tmax); % arg 2 is the referencec signal (should be lick)
normxc = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2); % normalize xc
nstd = 1; %round(excShortBin/bin);
g1 = gaussian(nstd, 2*nstd+1);
smthxc = smoothvect(normxc, g1); % smooth xc
plot(smthxc)
hold on
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

%% timemod shuff.. shuffle each swr - / + .5 s
if compute_shuffle
    % shuffle swr start time
    r = randi([-maxShift maxShift],length(swrTimes(:,1)), numshuffs)/1e3;
    swrStartSh = sort(swrTimes+r,1);
    c = 0;
    for i = 1:numshuffs
       %% xcorr
        xc = spikexcorr(swrStartSh(:,i), intraBoutXP, bin, tmax);
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
        
        c = c+1;
        out.normxcSh(c,:) = normxc;
        out.xcSh{c,1} = xc;
%         out.normxcShuf{end+1} = normxc;
        out.smthxcSh(c,:) = smthxc;
        out.excorrSh(c) = excorr;
%         out.xcrmsShuf{end+1} = xcrms;
%         out.p1Shuf{end+1} = p1;
%         out.p2Shuf{end+1} = p2;
%         out.expProbShuf{end+1} = expProb;
    end
end 

%% phasemod stuff
% get swr events within lick-burst intervals
burstSWRTimes = swrTimes(isIncluded(swrTimes, boutTimes));
burstSWRTimes = sort(burstSWRTimes(~isnan(burstSWRTimes)));
fprintf(' %d swrs within %d lickbouts. (%.02f pct swrs) \n', size(boutTimes,1), ...
    numel(burstSWRTimes), length(burstSWRTimes)/length(swrTimes)*100)
if isempty(burstSWRTimes)
    fprintf('no swrs in lick bouts for %d %d.. skipping\n', day, ep)
    return
end
% Get the swr-containing licks
% histc given licks as edges: get bin idx. lickedges(binidx) is the prior edge of bin, 
% aka the swr-preceding lick 
[~,~,swrLickBinidx] = histcounts(burstSWRTimes, intraBoutXP);
% exclude swr's in 0 bin, before first lick edge
% swrLickBinidx = swrLickBinidx(swrLickBinidx > 0); 
% burstSWRTimes = burstSWRTimes(swrLickBinidx > 0);
% for r = 1:length(burstSWRTimes)
%     
% end
% k = bsxfun(@minus, intraBoutXP(:), burstSWRTimes(:)'); 
% k(k >= 0) = -inf;
% [~,swrLickBinidx] = max(k);
% swrLickBinidx(all(k == -inf)) = 0;

% time since last lick
swrTimeSinceLick = burstSWRTimes - intraBoutXP(swrLickBinidx);
if any(swrTimeSinceLick > .250)
    error('wut')
end
% swr-containing ILI
for b = 1:length(burstSWRTimes)
    pre = max(intraBoutXP(intraBoutXP<=burstSWRTimes(b)));
    pst = min(intraBoutXP(intraBoutXP>=burstSWRTimes(b)));
    ili(b) = pst-pre;
end

%%
% 
% figure
% % clf
% xs = boutTimes(:,1);
% xe = boutTimes(:,2);
% yl = ylim;
% patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
%     length(xe)),'r', 'FaceAlpha',.2, 'edgecolor','none');
% 
% hold on
% l = line(intraBoutXP(:,[1 1])',ylim, 'linewidth',1, 'color', 'k');
% l = line(burstSWRTimes(:,[1 1])',ylim, 'linewidth',1, 'color', 'r');
% hold off
% 
% l = line(burstSWRTimes(:,[1 1])',ylim, 'linewidth',1, 'color', 'k');
%%

swrBinILI = intraBoutXP(swrLickBinidx+1) - intraBoutXP(swrLickBinidx); % +1 is post-containing lick
if any(swrBinILI > .250)
    error('wut')
end
% exclude swr-containing ili out of lick burst range 60 : 250 ms
swrValidILI = all([swrBinILI > minILIthresh swrBinILI < maxILIthresh],2);
fprintf('lBurst swrs in valid ILI: %d (%.02f pct swrs) \n', sum(swrValidILI), ...
    sum(swrValidILI)/length(swrTimes)*100)
burstSWRTimes = burstSWRTimes(swrValidILI);
if isempty(burstSWRTimes)
    return
end
swrBinILI = swrBinILI(swrValidILI);
swrTimeSinceLick = swrTimeSinceLick(swrValidILI);

% save swr-containing burst id
[~,~,swrLickBurstIdx] = histcounts(burstSWRTimes, [boutTimes(:,1); inf]);
burstContainSwr = boutTimes(swrLickBurstIdx,:);

% save the spike-containing burst lick order
% ultimately I want the count, relative to burst, for each spike-containing ILI
% so group the licks by lickburst ID, then get a per-burst cumcount with
% splitapply.. then, bc it's same dim as LBlicks, index from cumcount for each spike-containing ILI idx
[~, ~, lickBurstBinIdx] = histcounts(intraBoutXP, [boutTimes(:,1); inf]);
Y = cell2mat(splitapply(@(x) {[1:length(x)]'}, lickBurstBinIdx, lickBurstBinIdx));
swrBurstLickNum = Y(swrLickBinidx(swrValidILI));

% phasemod
swrPctSinceLick = swrTimeSinceLick ./ swrBinILI;
swrLickPhase = 2*pi*swrPctSinceLick; % pct of full cycle
meanvec = mean(exp(1i*swrLickPhase)); % get mean resultant vector
meanMRVmag = abs(meanvec); % vector magnitude
vecang = angle(meanvec);
% Z = meanMRVmag^2/n;
[pv, z] = circ_rtest(swrLickPhase); % z is mean res vec
% if Z ~= z
%     error('-----wut?')
% end
phasemod = log(z); % log variance normalizes (karalis,sirota)

out.swrStart = swrTimes;
out.licks = intraBoutXP;
% out.idlicks = idlicks;
out.boutTimes = boutTimes;
out.swrTimeSinceLick = swrTimeSinceLick;
out.swrBinILI = swrBinILI;
out.burstSWRStart = burstSWRTimes;
out.burstContainSwr = burstContainSwr;
out.dayEpoch = repmat([day ep], length(burstSWRTimes),1);

out.lickBurstBinIdx = lickBurstBinIdx;
out.swrBurstLickNum = swrBurstLickNum;
out.swrPctSinceLick = swrPctSinceLick;
out.swrLickPhase = swrLickPhase;
out.meanMRVmag = meanMRVmag;
out.vecang = vecang; % for polar scatter plots
out.phasemod = phasemod;

if compute_shuffle
    %% phasemod shuff
    out.swrLickphaseSh = 2*pi*rand(length(swrPctSinceLick), numshuffs);
    for s = 1:size(out.swrLickphaseSh,2)
        [~, z] = circ_rtest(out.swrLickphaseSh(:,s)); % z is mean res vec
        out.phasemodSh(s) = log(z); % log variance normalizes (karalis,sirota)
    end
%     fprintf('shuffle took %.02f s\n', toc);
end
end

function out = init_out()
out.index = [];
out.animal = [];
out.maxthresh = [];
out.swrEnd = [];
% xcTime
out.time = [];
out.xc = [];
out.normxc = [];
out.smthxc = [];
out.excorr = [];
out.xcrms = [];
out.p1 = [];
out.p2 = [];
out.expProb = [];

out.swrStart = [];
out.licks = [];
% out.idlicks = [];
out.boutTimes = [];
% iliPhase
out.swrTimeSinceLick = [];
out.swrBinILI = [];
out.burstSWRStart = [];
out.burstContainSwr = [];
out.dayEpoch = [];

out.lickBurstBinIdx = [];
out.swrBurstLickNum = [];
out.swrPctSinceLick = [];
out.swrLickPhase = [];
out.meanMRVmag = [];
out.vecang = [];
out.phasemod = [];

%% shuffle
% out.swrStartShuf = [];
% xcTime
out.xcSh = [];
out.normxcSh = [];
out.smthxcSh = [];
out.excorrSh = [];
% out.xcrmsShuf = {};
% out.p1Shuf = {};
% out.p2Shuf = {};
% out.expProbShuf = {};
% iliPhase
% out.swrBinILIShuf = {};
% out.swrBurstIntervalShuf = {};
% out.swrInBurstStartShuf = {};
% out.swrTimeSinceLickShuf = {};
% out.meanMRVmagShuf = {};
% out.swrLickPhaseShuf = {};
% out.swrPctSinceLickShuf = {};
% out.vecangShuf = {};
out.phasemodSh = [];
out.swrLickphaseSh = [];
end