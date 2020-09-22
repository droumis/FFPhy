
function out = dfa_reactivationPSTH(idx, timeFilter, varargin)
% out = dfa_reactivationPSTH(idx, timeFilter, varargin)
% computes 'reactivation strength' a la Peyrache 2009, 2010..
% required vargs : {'cellinfo', 'templateFilter', 'rippleFilter', 'burstFilter', 'cellfilter'};
% match times are currently hardcoded
%
%
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
cellFilter = '(isequal($area, ''ca1'')) && ($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
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
andef = animaldef(animal);

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
swrRXN.animal = animaldef(animal);
swrRXN.epochs{1} = [repmat(day, length(eps), 1), eps'];
rF = setfiltertime(swrRXN, rippleFilter);
rippleFilterOUT = cell2mat(rF.excludetime{1}');
incSWRs = ~isIncluded(swrTimes(:,1), rippleFilterOUT);
fprintf('%d of %d candidate swr events selected: d%d\n', sum(incSWRs), length(incSWRs), day)
swrTimes = swrTimes(incSWRs,:);

if isempty(swrTimes)
    return
end

%% load lick burst XP
[intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
    varargin{:}, 'minILIthresh', minILIthresh, 'maxILIthresh', maxILIthresh, 'minBoutLicks', ...
    minBoutLicks);
boutTimes = cell2mat({boutTimes{day}{eps}}');
intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
intraBoutXP = intraBoutXP(:,1);


d = bsxfun(@minus, swrTimes', intraBoutXP);
d(d < 0) = inf;
minDistFromXPtoNextSWR = min(d, [], 2);
XPnotNearSWRidx = minDistFromXPtoNextSWR > .5;

fprintf('%d XP within %d bouts and not near swr\n', sum(XPnotNearSWRidx), size(boutTimes,1))

% deal with the much larger # of xp than swrs in later days
% otherwise can't compute bc array too large
if sum(XPnotNearSWRidx) > sum(incSWRs)
    % of the xp indices, choose as many as there are swrs, randomly
    xpi = find(XPnotNearSWRidx);
    keepxp = sort(datasample(xpi,sum(incSWRs),'Replace',false));
    dropxp = xpi(~ismember(xpi, keepxp));
    XPnotNearSWRidx(dropxp) = 0;
    fprintf('***TOO MANY XP! using %d XP\n', sum(XPnotNearSWRidx))

end
% subsample to match # of swrs
% XPnotNearSWRidx
    
intraBoutXPnotSWR = intraBoutXP(XPnotNearSWRidx);

%% swr intra burst and extra burst idx
swrIntraBurstIdx = isIncluded(swrTimes, boutTimes);
swrExtraBurstIdx = ~isIncluded(swrTimes, boutTimes);
%% load spikes apply cellFilter
su = evaluatefilter(cellinfo, cellFilter);
numCells=size(su,1);
if numCells<numCellsThresh
    fprintf('%s %d %d not enough cells: %d < %d\n', animal, day, ep, numCells, numCellsThresh);
    return
end

spikes = loaddatastruct(andef{2}, animal, 'spikes', day);
epTimeRange = spikes{day}{ep}{su(1,3)}{su(1,4)}.timerange;
times = epTimeRange(1):.001:epTimeRange(2); % ms rez ep time


%% Get Template and Match Q (binned, z-scored spikes)
suKept = [];
c = 0;
timeBinEdges = epTimeRange(1):bin:epTimeRange(2);

swrRXN.animal = animaldef(animal);
swrRXN.epochs{1} = [repmat(day, length(eps), 1), eps'];
tF = setfiltertime(swrRXN, templateFilter);
templateExcludeTimes = cell2mat(tF.excludetime{1}');

templMask = ~isIncluded(timeBinEdges, templateExcludeTimes);
fprintf('run template pct eptime %.02f (%.03fpct s)\n',sum(templMask)/length(templMask)*100, ...
    sum(templMask)*bin)

idxWin = win(1)/bin:win(2)/bin;

% match spikes bin centers at event onset
Win = win(1):bin:win(2);
SWRmatchBinStackedTimes = bsxfun(@plus, swrTimes, Win);
XPmatchBinStackedTimes = bsxfun(@plus, intraBoutXPnotSWR, Win);

XPmatchSpikesZ = [];
SWRmatchSpikesZ = [];
for s = 1:length(su(:,1)) % for each su, 
    iSpks = spikes{day}{ep}{su(s,3)}{su(s,4)}.data;
    if length(iSpks) < numSpikesThresh
        continue
    end
    suKept = [suKept; s];
    c = c+1;
    spikesFRBin = histc(iSpks,timeBinEdges)*(1/bin); % FR
    MspikesFRBin(c) = nanmean(spikesFRBin);
    STDspikesFRBin(c) = nanstd(spikesFRBin);
    ZspikesFRBin(c,:) = (spikesFRBin - MspikesFRBin(c)) / STDspikesFRBin(c);
    
    % zscore matchSpikes using the mean and std from allspikes
    for ev = 1:size(SWRmatchBinStackedTimes,1)
        SWRmatchSpikesZ(c,:,ev) = (histc(iSpks, SWRmatchBinStackedTimes(ev,:))*(1/bin) - ...
            MspikesFRBin(c)) / STDspikesFRBin(c);
    end
    for ev = 1:size(XPmatchBinStackedTimes,1)
        XPmatchSpikesZ(c,:,ev) = (histc(iSpks, XPmatchBinStackedTimes(ev,:))*(1/bin) - ...
            MspikesFRBin(c)) / STDspikesFRBin(c);
    end
    
end
templateBinMask = ~isIncluded(timeBinEdges, templateExcludeTimes);
tmpSpkZBin = ZspikesFRBin(:, templateBinMask);
allMatchSpkZBin = ZspikesFRBin(:, ~templateBinMask); % use for zscoring match results
% stack horizontally
SWRmatchSpikesZ2D = SWRmatchSpikesZ(:,:); % cell x time x event -> cell x time concat
XPmatchSpikesZ2D = XPmatchSpikesZ(:,:);

% g = timeBinEdges(templMask);
% spikesTemplateTime = g(1:end-1)+bin/2;
su = su(suKept,:);
% out.binSpikes = binSpikes;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
% out.spikes = cell2mat(outSpikes);
fprintf('%d cells \n', length(suKept));
% get spike bin times. Each time in reactTime will match a reactivation
% strength score

%% get Match Spikes
% match spikes are binned a little differently so that they are centered on
% event onset

% for ev = 1:size(matchBinStackedTimes,1)
%     for s = 1:length(su(:,1))
%         iSpks = spikes{day}{ep}{su(s,3)}{su(s,4)}.data;
%         matchSpikes{ev}(s,:) = histc(iSpks, matchBinStackedTimes(e,:));
%     end
% end

% out.reactTime = reactTime;
%% Z-scored binned spikes, this is the 'Q' matrix in peyrache
% zSpikesTemplate=zscore(templateZspikesFRBin')';

% zSpikesMatch=zscore(binSpikesMatch')';
% numMatchBins=size(binSpikesMatch,2);

%% get swr conditions times, idx

% lag0idx = ceil(length(idxWin)/2);
etaTime = idxWin*bin;
% % this is NOT swr times!
% swrTime = vec2list(~isExcluded(times, timeFilter), times);

%
% % get nearest bin less than each swr start time
% [~, swrTimeIdx] = histc(swrTimes(:,1), timeBinEdges);
% % exclude swr too close to beginning or end to get a full window of data
% swrTimeIdx = swrTimeIdx(swrTimeIdx+min(idxWin)>0);
% swrTimeIdx = swrTimeIdx(swrTimeIdx+max(idxWin)<size(zSpikesMatch,2));
% if isempty(swrTimeIdx)
%     fprintf('no swr\n');
%     return
% end
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
nBinsTemplate = size(tmpSpkZBin, 2);
zCCtemplate=(tmpSpkZBin * tmpSpkZBin') / nBinsTemplate; %effectively same as corrcoef(ZspikeBins')
zCCtemplateNoDiag = zCCtemplate - diag(diag(zCCtemplate)); % cancel out the diag i=j terms
% ceil(size(matchSpikesZ,2)/2)
%% Full Model version of Reactivation Strength
if fullModel
    %% SWR rxn psth
    swrRXN2D = diag(SWRmatchSpikesZ2D' * zCCtemplateNoDiag * SWRmatchSpikesZ2D)';
    swrRXN = reshape(swrRXN2D, [size(SWRmatchSpikesZ,2) size(SWRmatchSpikesZ,3)])';
    
        
    % SWR-intraBurst rxn psth
    swrIntraBurstRXN = swrRXN(swrIntraBurstIdx,:);
    swrIntraBurstRXN_mean = nanmean(swrIntraBurstRXN);
    swrIntraBurstRXN_sem = sem(swrIntraBurstRXN,1);
    
    % SWR-extraBurst rxn psth
    swrExtraBurstRXN = swrRXN(swrExtraBurstIdx, :);
    swrExtraBurstRXN_mean = nanmean(swrExtraBurstRXN);
    swrExtraBurstRXN_sem = sem(swrExtraBurstRXN,1);
    
    
    % Zscore the swr RXN based on MEAN, STD RXN of ALL NON-TEMPLATE TIMES
    RXNallMatchSpkZBin = diag(allMatchSpkZBin' * zCCtemplateNoDiag * allMatchSpkZBin);
    meanRXNallMatchSpkZBin = nanmean(RXNallMatchSpkZBin);
    stdRXNallMatchSpkZBin = nanstd(RXNallMatchSpkZBin);
    swrRXNZ = (swrRXN - meanRXNallMatchSpkZBin) / stdRXNallMatchSpkZBin;

        
    % SWR-intraBurst rxn psth
    swrIntraBurstRXNZ = swrRXNZ(swrIntraBurstIdx,:);
    swrIntraBurstRXNZ_mean = nanmean(swrIntraBurstRXNZ);
    swrIntraBurstRXNZ_sem = sem(swrIntraBurstRXNZ,1);
    
    % SWR-extraBurst rxn psth
    swrExtraBurstRXNZ = swrRXNZ(swrExtraBurstIdx, :);
    swrExtraBurstRXNZ_mean = nanmean(swrExtraBurstRXNZ);
    swrExtraBurstRXNZ_sem = sem(swrExtraBurstRXNZ,1);
    
    %% lick rxn psth, no SWR
    xpRXN2D=diag(XPmatchSpikesZ2D' * zCCtemplateNoDiag * XPmatchSpikesZ2D)';
    XPnoswrRXN = reshape(xpRXN2D, [size(XPmatchSpikesZ,2) size(XPmatchSpikesZ,3)])';
    XPnoswrRXN_mean = nanmean(XPnoswrRXN);
    XPnoswrRXN_sem = sem(XPnoswrRXN,1);
    
    % Zscore the XP RXN based on MEAN, STD RXN of ALL NON-TEMPLATE TIMES
    XPnoswrRXNZ =  (XPnoswrRXN - meanRXNallMatchSpkZBin) / stdRXNallMatchSpkZBin;
    XPnoswrRXNZ_mean = nanmean(XPnoswrRXNZ);
    XPnoswrRXNZ_sem = sem(XPnoswrRXNZ,1);
end

% plot([swrIntraBurstRXN_mean' swrExtraBurstRXN_mean' XPnoswrRXN_mean'])

%% PER PC
% matrix is symmetric so divide by 2 and number of bins for normalization
if perPC
    %% decomposing correlation matrix into eigenvalues and eigenvectors
    [eigVec, eigVal] = eig(zCCtemplate); % now C=V*D*V'
    % Sort, so that the first explains most variance
    [eigValSort, eigValSortIdx] = sort(diag(eigVal),'descend');
    eigVecSort = eigVec(:,eigValSortIdx);
    
    % get significant eigenvalues
    eigValSigThresh = (1+sqrt(1/nBinsTemplate/numCells))^2;
    eigValSortSig = eigValSort(eigValSort>eigValSigThresh);
    eigVecSortSig = eigVecSort(:,eigValSort>eigValSigThresh);
    numEigValSig = length(eigValSortSig);
    
    for iEig = 1:numEigValSig
        PC = eigVecSortSig(:,iEig)*eigVecSortSig(:,iEig)'*eigValSortSig(iEig);
        PCsNoDiag(iEig,:,:) = PC - diag(diag(PC));
    end

    %% compute reactivation trace per sig PC
    rxnPerPCSWR = [];
    reactPerPCXP = [];
    for jEig=1:numEigValSig % for each PC
        
        numMatchBinsSWR = size(SWRmatchSpikesZ2D,2);
        for itb=1:numMatchBinsSWR % for each time bin in the SWR match psth times
            % outer product of the instantaneous ensemble firing
            % in each i,j is the product of firing of cell i and cell j (in this time bin).
            ensOutProd=SWRmatchSpikesZ2D(:,itb)*SWRmatchSpikesZ2D(:,itb)';
            % dot-product ensembleOuterProd with corr matrix (whose diagonal has been 0'd)
            rxnPerPCSWR(itb) = sum(sum(ensOutProd.*squeeze(PCsNoDiag(jEig,:,:))));
        end
        swrRXNperPC = reshape(rxnPerPCSWR, [size(SWRmatchSpikesZ,2) size(SWRmatchSpikesZ,3)])';
        
        swrIntraBurstRXNperPC(:,:,jEig) = swrRXNperPC(swrIntraBurstIdx,:);
        swrExtraBurstRXNperPC(:,:,jEig) = swrRXNperPC(swrExtraBurstIdx,:);
        
        numMatchBinsXP = size(XPmatchSpikesZ2D,2);
        for itb=1:numMatchBinsXP
            % outer product of the instantaneous ensemble firing
            % in each i,j is the product of firing of cell i and cell j (in this time bin).
            ensOutProd=XPmatchSpikesZ2D(:,itb)*XPmatchSpikesZ2D(:,itb)';
            % dot-product ensembleOuterProd with corr matrix (whose diagonal has been 0'd)
            reactPerPCXP(itb)=sum(sum(ensOutProd.*squeeze(PCsNoDiag(jEig,:,:))));
        end
        xpRXNperPC(:,:,jEig) = reshape(reactPerPCXP, [size(XPmatchSpikesZ,2) size(XPmatchSpikesZ,3)])';
    end
    swrIntraBurstRXNperPC_mean = squeeze(nanmean(swrIntraBurstRXNperPC, 1))';
    swrIntraBurstRXNperPC_sem = squeeze(sem(swrIntraBurstRXNperPC, 1))';
    
    swrExtraBurstRXNperPC_mean = squeeze(nanmean(swrExtraBurstRXNperPC, 1))';
    swrExtraBurstRXNperPC_sem = squeeze(sem(swrIntraBurstRXNperPC, 1))';
    
    xpRXNperPC_mean = squeeze(nanmean(xpRXNperPC, 1))';
    xpRXNperPC_sem  = squeeze(sem(xpRXNperPC, 1))';
    
    %         %% SWR rxn psth per PC
    %         [out.swrReactPSTHPerPC{jEig}, out.swrReactETAPerPC{jEig}, out.swrReactETAPerPCShufs{jEig}] = ...
    %             stackPSTH(out.reactPerPC{jEig}, swrTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
    %         %% SWR-Burst rxn psth
    %         if ~isempty(swrBurstTimeIdx)
    %             [out.swrBurstReactPSTHPerPC{jEig}, out.swrBurstReactETAPerPC{jEig}, out.swrBurstReactETAPerPCShufs{jEig}] = ...
    %                 stackPSTH(out.reactPerPC{jEig}, swrBurstTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
    %         end
    %         %% lick rxn psth
    %         [out.lickReactPSTHPerPC{jEig}, out.lickReactETAPerPC{jEig}, out.lickReactETAPerPCShufs{jEig}] = ...
    %             stackPSTH(out.reactPerPC{jEig}, lickBurstTimeIdx, 'idxWin', idxWin, 'numShufs', numShufs);
    %     end
    
    out.index = idx;
    out.animal = animal;
    out.maxthresh = maxthresh;
    out.swrEnd = swrEnd;
    out.su = su;
    out.eigVecSortSig = eigVecSortSig;
    out.eigValSortSig = eigValSortSig;
    out.etaTime = etaTime;
    
    % full
    out.swrIntraBurstRXN = swrIntraBurstRXN;
    out.swrIntraBurstRXN_mean = swrIntraBurstRXN_mean;
    out.swrIntraBurstRXN_sem = swrIntraBurstRXN_sem;    
    
    out.swrIntraBurstRXNZ = swrIntraBurstRXNZ;
    out.swrIntraBurstRXNZ_mean = swrIntraBurstRXNZ_mean;
    out.swrIntraBurstRXNZ_sem = swrIntraBurstRXNZ_sem;    
    
    out.swrExtraBurstRXN = swrExtraBurstRXN;
    out.swrExtraBurstRXN_mean = swrExtraBurstRXN_mean;
    out.swrExtraBurstRXN_sem = swrExtraBurstRXN_sem;
    
    out.swrExtraBurstRXNZ = swrExtraBurstRXNZ;
    out.swrExtraBurstRXNZ_mean = swrExtraBurstRXNZ_mean;
    out.swrExtraBurstRXNZ_sem = swrExtraBurstRXNZ_sem;
    
    out.XPnoswrRXN = XPnoswrRXN;
    out.XPnoswrRXN_mean = XPnoswrRXN_mean;
    out.XPnoswrRXN_sem = XPnoswrRXN_sem;
    
    out.XPnoswrRXNZ = XPnoswrRXNZ;
    out.XPnoswrRXNZ_mean = XPnoswrRXNZ_mean;
    out.XPnoswrRXNZ_sem = XPnoswrRXNZ_sem;
    
    % per pc
    out.swrIntraBurstRXNperPC_mean = swrIntraBurstRXNperPC_mean ;
    out.swrIntraBurstRXNperPC_sem = swrIntraBurstRXNperPC_sem;
    
    out.swrExtraBurstRXNperPC_mean = swrExtraBurstRXNperPC_mean;
    out.swrExtraBurstRXNperPC_sem = swrExtraBurstRXNperPC_sem;
    
    out.xpRXNperPC_mean = xpRXNperPC_mean;
    out.xpRXNperPC_sem  = xpRXNperPC_sem;
    
    
end

end

function out = init_out()
out.index = [];
out.animal = [];
out.maxthresh = [];
out.swrEnd = [];
out.su = [];
out.eigVecSortSig = [];
out.eigValSortSig = [];
out.etaTime = [];

out.swrIntraBurstRXN = [];
out.swrIntraBurstRXN_mean = [];
out.swrIntraBurstRXN_sem = [];    

out.swrIntraBurstRXNZ = [];
out.swrIntraBurstRXNZ_mean = [];
out.swrIntraBurstRXNZ_sem = [];    

out.swrExtraBurstRXN = [];
out.swrExtraBurstRXN_mean = [];
out.swrExtraBurstRXN_sem = [];

out.swrExtraBurstRXNZ = [];
out.swrExtraBurstRXNZ_mean = [];
out.swrExtraBurstRXNZ_sem = [];

out.XPnoswrRXN = [];
out.XPnoswrRXN_mean = [];
out.XPnoswrRXN_sem = [];

out.XPnoswrRXNZ = [];
out.XPnoswrRXNZ_mean = [];
out.XPnoswrRXNZ_sem = [];

% per pc
out.swrIntraBurstRXNperPC_mean = [];
out.swrIntraBurstRXNperPC_sem = [];

out.swrExtraBurstRXNperPC_mean = [];
out.swrExtraBurstRXNperPC_sem = [];

out.xpRXNperPC_mean = [];
out.xpRXNperPC_sem  = [];
end