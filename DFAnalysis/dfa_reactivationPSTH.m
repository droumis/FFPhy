
function out = dfa_reactivationPSTH(idx, timeFilter, varargin)
% out = dfa_reactivationPSTH(idx, timeFilter, varargin)
% computes 'reactivation strength' a la Peyrache 2009, 2010..
% required vargs : {'cellinfo', 'templateFilter', 'rippleFilter', 'burstFilter', 'cellfilter'};
% match times are currently hardcoded
%
%                 Salubrious Squirrel
%                                              ,.-'''     '-.
%                                           ,'     _        '-
%                                         ,'     -'           \
%                                        /                     \
%                                        '                      |
%              _                        |              ..       \
%            _(C)                       /                        |
%          ,'_   `.                     |       -''-.__        _'
%         '(o `    `.                   /            __',,--:,-'
%       ,'      .    `._.--....._      |      __.--''  '
%      @<<       ;               '-.   |     /'
%       |/------'     .              `\|    ,'
% ___________ ,`; ,'   ;                    | ____________________
%            /'  ;   ,'        ,-         _,'
%          ,'  ,'  ,'        .'          |
%   _.(._,' _,'  ,' `._    .'            /
%  ((_O)))''_, -'      `-..'           ,'
%   |   |//'              M         ,''
%    `-'                   J.     ,'
%                   ,------'P     /
%                  ((.__________.'
% %
% Iterator:
% - singleDayAnal (city)
% - no longer: singleepochanal (space)
% 
% args:
% - idx [day epochs]
% - timeFilter: applied to swr
%               .. not applied to licks because they need a different
%               timefilter.. currently 'within XPburst intervals' is
%               hardcorded
% 
% varargs:
% - data (required) {'ca1rippleskons', 'lick'}
% - 
% - win:
% - eventType:
% - LFPTypes: 
%
%
%
%
%
%{

Notes:
- city:fireworks
- used to use the timeFilter as intervals for the TEMPLATE spikes
        % i.e. timefilter = {'get2dstate', '$velocity > 4'};

FFPhy V0.1
@DKR

%}

% check for required data in varargin
reqData = {'cellinfo'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error(sprintf('missing varg: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end
% assign params
% eventType = 'ca1rippleskons';
byDay = 1; % currently requires use of singleDayAnal (city)
minILIthresh = .06; % seconds
maxILIthresh = .250; % seconds
minBoutLicks = 3;

numCellsThresh = 5;
numSpikesThresh = 100;
minLickTimeFromSwr = 1;
numShufs = 100;
bin = 0.01; % seconds
fullModel = 1;
perPC = 1;
consensus_numtets = 2;   % minimum # of tets for consensus event detection
minstdthresh = 2;        % STD. how big your ripples are
exclusion_dur = .5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
rippleFilter = {{'getconstimes', '($cons == 1)', ...
    'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
    'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
    'minvelocity', minvelocity,'maxvelocity',maxvelocity}};

templateFilter = {{'get2dstate','($velocity>4)'}};

if ~isempty(varargin)
    assign(varargin{:})
end

day = idx(1);
if byDay
    eps = idx(2:end);
else
    eps = idx(2);
end

% init output
out = init_out(); % init output
out.index = idx;
out.animal = animal;
andef = animaldef(animal);

%% load lick burst XP
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:}, 'minILIthresh', minILIthresh, 'maxILIthresh', maxILIthresh, 'minBoutLicks', ...
    minBoutLicks);
boutTimes = cell2mat({boutTimes{day}{eps}}');
intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
intraBoutXP = intraBoutXP(:,1);
fprintf('%d XP within %d bursts \n', numel(intraBoutXP), size(boutTimes,1))

%% Get SWRs
swr = loaddatastruct(andef{2}, animal, 'ca1rippleskons', day); 
% evid = find(contains(varargin(1:2:end), 'ca1rippleskons'));
% o = [1:2:length(varargin)]+1;
% swr = varargin{o(evid)};
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

% apply RIPPLE FILTER
f.animal = animaldef(animal);
f.epochs{1} = [repmat(day, length(eps), 1), eps'];
rF = setfiltertime(f, rippleFilter);
rippleFilterOUT = cell2mat(rF.excludetime{1}');
incSWRs = ~isIncluded(swrTimes(:,1), rippleFilterOUT);
fprintf('%d of %d candidate swr events selected: d%d\n', sum(incSWRs), length(incSWRs), day)
swrTimes = swrTimes(incSWRs,:);
out.maxthresh = maxthresh;
out.swrEnd = swrEnd;

if isempty(swrTimes)
    return
end

%% Filter for SU/cells
su = evaluatefilter(cellinfo, cellfilter);
numCells=size(su,1);
if numCells<numCellsThresh
    fprintf('%s %d %d not enough cells: %d < %d\n', animal, day, ep, numCells, numCellsThresh);
    return
end

%% load spikes for day
spikes = loaddatastruct(andef{2}, animal, 'spikes', day);
epTimeRange = spikes{day}{ep}{su(1,3)}{su(1,4)}.timerange;
times = epTimeRange(1):.001:epTimeRange(2); % ms rez ep time


%% Get Template spikes
suKept = [];
c = 0;
timeBinEdges = epTimeRange(1):bin:epTimeRange(2);

f.animal = animaldef(animal);
f.epochs{1} = [repmat(day, length(eps), 1), eps'];
tF = setfiltertime(f, templateFilter);
tFilterOUT = cell2mat(tF.excludetime{1}');

templMask = ~isIncluded(timeBinEdges, tFilterOUT);
fprintf('run template pct eptime %.02f (%.03f s)\n',sum(templMask)/length(templMask)*100, ...
    sum(templMask)*bin)
templateSpikes = [];
for s = 1:length(su(:,1))
    iSpks = spikes{day}{ep}{su(s,3)}{su(s,4)}.data;
    if length(iSpks) < numSpikesThresh
        continue
    end
    suKept = [suKept; s];
    c = c+1;
    allSpikes(c,:) = zscore(histc(iSpks,timeBinEdges)')';
    hmmm...
    iSpks = iSpks(~isIncluded(iSpks, tFilterOUT));
    templateSpikes(c,:) = histc(iSpks,timeBinEdges);
end
templateBinMask = ~isIncluded(timeBinEdges, tFilterOUT);

% g = timeBinEdges(templMask);
% spikesTemplateTime = g(1:end-1)+bin/2;
su = su(suKept,:);
out.su = su;
% out.binSpikes = binSpikes;
% out.spikes = cell2mat(outSpikes);
fprintf('%d cells \n', length(suKept));
% get spike bin times. Each time in reactTime will match a reactivation
% strength score

%% get Match Spikes
% match spikes are binned a little differently so that they are centered on
% event onset
winCenterTimes = swrTimes;
idxWin = win(1)/bin:win(2)/bin;
Win = win(1):bin:win(2);
matchBinStackedTimes = bsxfun(@plus, winCenterTimes, Win);
for ev = 1:size(matchBinStackedTimes,1)
    for s = 1:length(su(:,1))
        iSpks = spikes{day}{ep}{su(s,3)}{su(s,4)}.data;      
        matchSpikes{ev}(s,:) = histc(iSpks, matchBinStackedTimes(e,:));
    end
end

out.reactTime = reactTime;
%% Z-scored binned spikes, this is the 'Q' matrix in peyrache
zSpikesTemplate=zscore(templateSpikes')';

zSpikesMatch=zscore(binSpikesMatch')';
numMatchBins=size(binSpikesMatch,2);

%% get swr conditions times, idx

% lag0idx = ceil(length(idxWin)/2);
etaTime = idxWin*bin;
% % this is NOT swr times!
% swrTime = vec2list(~isExcluded(times, timeFilter), times);


% get nearest bin less than each swr start time
[~, swrTimeIdx] = histc(swrTimes(:,1), timeBinEdges);
% exclude swr too close to beginning or end to get a full window of data
swrTimeIdx = swrTimeIdx(swrTimeIdx+min(idxWin)>0);
swrTimeIdx = swrTimeIdx(swrTimeIdx+max(idxWin)<size(zSpikesMatch,2));
if isempty(swrTimeIdx)
    fprintf('no swr\n');
    return
end
% % burst intervals, and swr's, licks within them
% l = setfiltertime(f, burstFilter);
% 
% swrBurstTime = swrTime(~isExcluded(swrTime(:,1), l.excludetime{1}{1}),:);
% 
% burstIntervals = vec2list(~isExcluded(times, l.excludetime{1}{1}), times);
% lickBurstTime = lick{day}{ep}.starttime(~isExcluded(lick{day}{ep}.starttime, l.excludetime{1}{1}));
% lickBurstTime = lickBurstTime(lickBurstTime+min(idxWin)>0);
% lickBurstTime = lickBurstTime(lickBurstTime+max(idxWin)<size(zSpikesMatch,2));
% % exclude licks near swrs
% [~, d] = knnsearch(swrBurstTime(:,1), lickBurstTime);
% lickBurstTime = lickBurstTime(d>minLickTimeFromSwr,:);
% 
% % exclude events, intervals too close to ep edges
% [~, swrBurstTimeIdx] = histc(swrBurstTime(:,1), matchTimeEdges);
% swrBurstTimeIdx = swrBurstTimeIdx(swrBurstTimeIdx+min(idxWin)>0);
% swrBurstTimeIdx = swrBurstTimeIdx(swrBurstTimeIdx+max(idxWin)<size(zSpikesMatch,2));
% 
% [~, lickBurstTimeIdx] = histc(lickBurstTime, matchTimeEdges);
% lickBurstTimeIdx = lickBurstTimeIdx(lickBurstTimeIdx+min(idxWin)>0);
% lickBurstTimeIdx = lickBurstTimeIdx(lickBurstTimeIdx+max(idxWin)<size(zSpikesMatch,2));
% 
% % out
% out.etaTime = etaTime;
% out.swrTime = swrTime;
% out.swrBurstTime = swrBurstTime;
% out.burstIntervals = burstIntervals;
% out.lickBurstTime = lickBurstTime;

%% create diagless-cc template mat (common to full and perpc)
% create pairwise cell activity correlation matrix
zCCtemplate=(zSpikesTemplate * zSpikesTemplate') / numBinsTemplate; %effectively same as corrcoef(ZspikeBins')
zCCtemplateNoDiag = zCCtemplate - diag(diag(zCCtemplate)); % cancel out the diag i=j terms

%% Full Model version of Reactivation Strength
if fullModel
    reactFull=diag(zSpikesMatch' * zCCtemplateNoDiag * zSpikesMatch)';
    out.reactFull = reactFull;
    %% SWR rxn psth
    [out.swrReactPSTHfull, out.swrReactETAfull, out.swrReactETAfullShufs] = ...
        stackPSTH(reactFull, swrTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
    %% SWR-Burst rxn psth
    if ~isempty(swrBurstTimeIdx)
        [out.swrBurstReactPSTHfull, out.swrBurstReactETAfull, out.swrBurstReactETAfullShufs] = ...
            stackPSTH(reactFull, swrBurstTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
    end
    %% lick rxn psth
    [out.lickReactPSTHfull, out.lickReactETAfull, out.lickReactETAfullShufs] = ...
        stackPSTH(reactFull, lickBurstTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
end

%% PER PC
% matrix is symmetric so divide by 2 and number of bins for normalization
if perPC
    %% decomposing correlation matrix into eigenvalues and eigenvectors
    [eigVec, eigVal] = eig(zCCtemplate); % now C=V*D*V'
    % Sort, so that the first explains most variance
    [eigValSort, eigValSortIdx] = sort(diag(eigVal),'descend');
    eigVecSort = eigVec(:,eigValSortIdx);
    
    % get significant eigenvalues
    eigValSigThresh = (1+sqrt(1/numBinsTemplate/numCells))^2;
    eigValSortSig = eigValSort(eigValSort>eigValSigThresh);
    eigVecSortSig = eigVecSort(:,eigValSort>eigValSigThresh);
    numEigValSig = length(eigValSortSig);
    
    for iEig = 1:numEigValSig
        PC = eigVecSortSig(:,iEig)*eigVecSortSig(:,iEig)'*eigValSortSig(iEig);
        PCsNoDiag(iEig,:,:) = PC - diag(diag(PC));
    end
    
    out.eigVecSortSig = eigVecSortSig;
    out.eigValSortSig = eigValSortSig;
    %% compute reactivation trace per sig PC
    for jEig=1:numEigValSig
        for it=1:numMatchBins
            % outer product of the instantaneous ensemble firing
            % in each i,j is the product of firing of cell i and cell j (in this time bin).
            ensOutProd=zSpikesMatch(:,it)*zSpikesMatch(:,it)';
            % dot-product ensembleOuterProd with corr matrix (whose diagonal has been 0'd)
            out.reactPerPC{jEig}(it)=sum(sum(ensOutProd.*squeeze(PCsNoDiag(jEig,:,:))));
        end
        %% SWR rxn psth per PC
        [out.swrReactPSTHPerPC{jEig}, out.swrReactETAPerPC{jEig}, out.swrReactETAPerPCShufs{jEig}] = ...
            stackPSTH(out.reactPerPC{jEig}, swrTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
        %% SWR-Burst rxn psth
        if ~isempty(swrBurstTimeIdx)
            [out.swrBurstReactPSTHPerPC{jEig}, out.swrBurstReactETAPerPC{jEig}, out.swrBurstReactETAPerPCShufs{jEig}] = ...
                stackPSTH(out.reactPerPC{jEig}, swrBurstTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
        end
        %% lick rxn psth
        [out.lickReactPSTHPerPC{jEig}, out.lickReactETAPerPC{jEig}, out.lickReactETAPerPCShufs{jEig}] = ...
            stackPSTH(out.reactPerPC{jEig}, lickBurstTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
    end
end

end

function out = init_out()
out.idx = [];
out.animal = [];
out.su = [];
% out.binSpikes = [];
% out.spikes = [];
out.swrTime = [];
out.swrBurstTime = [];
out.etaTime = [];

out.reactFull = []; % full epoch 'match' intervals reactivation strength series
out.reactTime = []; % full epoch 'match' time.. aka reactivation

out.swrReactPSTHfull = [];
out.swrReactETAfull = [];
out.swrReactETAfullShufs = [];

out.swrBurstReactPSTHfull = [];
out.swrBurstReactETAfull = [];
out.swrBurstReactETAfullShufs = [];

out.burstIntervals = [];
out.lickBurstTime = [];
out.lickReactPSTHfull = [];
out.lickReactETAfull = [];
out.lickReactETAfullShufs = [];

% per pc
out.eigVecSortSig = [];
out.eigValSortSig = [];
out.reactPerPC = [];

out.swrReactPSTHPerPC = [];
out.swrReactETAPerPC = [];
out.swrReactETAPerPCShufs = [];

out.swrBurstReactPSTHPerPC = [];
out.swrBurstReactETAPerPC = [];
out.swrBurstReactETAPerPCShufs = [];

out.lickReactPSTHPerPC = [];
out.lickReactETAPerPC = [];
out.lickReactETAPerPCShufs = [];
end