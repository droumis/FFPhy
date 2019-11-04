
% uses the excludeIntervals as intervals for the TEMPLATE spikes
% this would probably be high velocity times, created with the timefilter:
% timefilter = {'get2dstate', '$velocity > 4'};
% uses the lick data to generate lickBurstIntervals


function out = dfa_reactivationPSTH(idx, nonSwrIntervals, varargin)

% check for required data in varargin
reqData = {'spikes', 'lick', 'ca1rippleskons', 'cellinfo'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error(sprintf('missing data: %s ', reqData{~ismember(reqData,varargin(1:2:end))}));
    end
end

% assign params
numCellsThresh = 5;
numSpikesThresh = 50;
minLickTimeFromSwr = 1;
numShufs = 100;
bin = 0.01; % seconds
fullModel = 1;
perPC = 0;
if ~isempty(varargin)
    assign(varargin{:})
end

% init output
out = init_out(idx, animal);
day = idx(1);
ep = idx(2);

%% Filter for SU
su = evaluatefilter(cellinfo, cellfilter);
numCells=size(su,1);
if numCells<numCellsThresh
    fprintf('%s %d %d not enough cells: %d < %d\n', animal, day, ep, numCells, numCellsThresh);
    return
end
epTimeRange = spikes{day}{ep}{su(1,3)}{su(1,4)}.timerange;
times = epTimeRange(1):.001:epTimeRange(2); % ms rez ep time

%% Bin SU spikes
timeBinEdges = epTimeRange(1):bin:epTimeRange(2);
for s = 1:length(su(:,1))
    iSpks = spikes{day}{ep}{su(s,3)}{su(s,4)}.data;
    if length(iSpks) < numSpikesThresh
        fprintf('%s %d %d not enough spikes: %d < %d\n', animal, day, ep, length(iSpks), numSpikesThresh);
        return
    end
    outSpikes{s,1} = [repmat(s,length(iSpks),1) iSpks];
    binSpikes(s,:) = histc(iSpks,timeBinEdges);
end
out.su = su;
out.binSpikes = binSpikes;
out.spikes = cell2mat(outSpikes);

%% Get Template- and Match- binned zSpikes
f.animal = animaldef(animal);
f.epochs{1} = idx;
tF = setfiltertime(f, templateFilter);
templMask = ~isExcluded(timeBinEdges, tF.excludetime{1}{1});
fprintf('template pct eptime %.02f (%.03f s)\n',sum(templMask)/length(templMask)*100, ...
    sum(templMask)*bin)
% get spike bin times
g = timeBinEdges(templMask);
spikesTemplateTime = g(1:end-1)+bin/2;
matchTimeEdges = timeBinEdges(~templMask);
reactTime = matchTimeEdges(1:end-1)+bin/2;
% bin template, match spikes
binSpikesTemplate = binSpikes(:,templMask);
numBinsTemplate=size(binSpikesTemplate,2);
binSpikesMatch = binSpikes(:,~templMask);
% Z-transform binned firing of each cell
zSpikesTemplate=zscore(binSpikesTemplate')';
zSpikesMatch=zscore(binSpikesMatch')';
numMatchBins=size(binSpikesMatch,2);
% out
out.reactTime = reactTime;

%% get swr conditions times, idx
idxWin = win(1)/bin:win(2)/bin;
lag0idx = ceil(length(idxWin)/2);
etaTime = idxWin*bin;
swrTime = vec2list(~isExcluded(times, nonSwrIntervals), times);
% get nearest bin less than each swr start time
[~, swrTimeIdx] = histc(swrTime(:,1), matchTimeEdges);
% exclude swr too close to beginning or end to get a full window of data
swrTimeIdx = swrTimeIdx(swrTimeIdx+min(idxWin)>0);
swrTimeIdx = swrTimeIdx(swrTimeIdx+max(idxWin)<size(zSpikesMatch,2));
if isempty(swrTimeIdx)
    fprintf('no swr\n');
    return
end
% burst intervals, and swr's, licks within them
l = setfiltertime(f, burstTimeFilter);

swrBurstTime = swrTime(~isExcluded(swrTime(:,1), l.excludetime{1}{1}),:);

burstIntervals = vec2list(~isExcluded(times, l.excludetime{1}{1}), times);
lickBurstTime = lick{day}{ep}.starttime(~isExcluded(lick{day}{ep}.starttime, l.excludetime{1}{1}));
lickBurstTime = lickBurstTime(lickBurstTime+min(idxWin)>0);
lickBurstTime = lickBurstTime(lickBurstTime+max(idxWin)<size(zSpikesMatch,2));
% exclude licks near swrs
[~, d] = knnsearch(swrBurstTime(:,1), lickBurstTime);
lickBurstTime = lickBurstTime(d>minLickTimeFromSwr,:);

% exclude events, intervals too close to ep edges
[~, swrBurstTimeIdx] = histc(swrBurstTime(:,1), matchTimeEdges);
swrBurstTimeIdx = swrBurstTimeIdx(swrBurstTimeIdx+min(idxWin)>0);
swrBurstTimeIdx = swrBurstTimeIdx(swrBurstTimeIdx+max(idxWin)<size(zSpikesMatch,2));

[~, lickBurstTimeIdx] = histc(lickBurstTime, matchTimeEdges);
lickBurstTimeIdx = lickBurstTimeIdx(lickBurstTimeIdx+min(idxWin)>0);
lickBurstTimeIdx = lickBurstTimeIdx(lickBurstTimeIdx+max(idxWin)<size(zSpikesMatch,2));

% out
out.etaTime = etaTime;
out.swrTime = swrTime;
out.swrBurstTime = swrBurstTime;
out.burstIntervals = burstIntervals;
out.lickBurstTime = lickBurstTime;

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

function out = init_out(idx, animal, varargin)
out.idx = idx;
out.animal = animal;
out.su = [];
out.binSpikes = [];
out.spikes = [];
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