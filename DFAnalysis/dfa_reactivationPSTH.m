
% uses the excludeIntervals as intervals for the TEMPLATE spikes
% this would probably be high velocity times, created with the timefilter:
% timefilter = {'get2dstate', '$velocity > 4'};
% uses the lick data to generate lickBurstIntervals


function out = dfa_reactivationPSTH(idx, templateIntervals, varargin)

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
numShuffs = 100;
bin = 0.01; % seconds
fullModel = 1;
perPC = 0;
if ~isempty(varargin)
    assign(varargin{:})
end

% init output
out = init_out(idx, animal);

%% Filter and bin All su spikes
su = evaluatefilter(cellinfo, cellfilter);
numCells=size(su,1);
if numCells<numCellsThresh
    fprintf('%s %d %d not enough cells: %d < %d\n', animal, idx(1), idx(2), numCells, numCellsThresh);
    return
end

epTimeRange = spikes{idx(1)}{idx(2)}{su(1,3)}{su(1,4)}.timerange;
times = epTimeRange(1):.001:epTimeRange(2); % ms rez ep time
andef = animaldef(animal);
f.animal = animaldef(animal);
f.epochs{1} = idx;

timeBinEdges = epTimeRange(1):bin:epTimeRange(2);
for s = 1:length(su(:,1))
    iSpks = spikes{idx(1)}{idx(2)}{su(s,3)}{su(s,4)}.data;
    if length(iSpks) < numSpikesThresh
        fprintf('%s %d %d not enough spikes: %d < %d\n', animal, idx(1), idx(2), length(iSpks), numSpikesThresh);
        return
    end
    outSpikes{s,1} = [repmat(s,length(iSpks),1) iSpks];
    binSpikes(s,:) = histc(iSpks,timeBinEdges);
end
out.su = su;
out.binSpikes = binSpikes;
out.spikes = cell2mat(outSpikes);
%% Get Template and Match binned zSpikes
templMask = logical(isExcluded(timeBinEdges, templateIntervals));
fprintf('template pct eptime %.02f (%.03f s)\n',sum(templMask)/length(templMask)*100, ...
    sum(templMask)*bin)
% get time bins
tmpTimeEdges = timeBinEdges(templMask);
spikesTemplateTime = tmpTimeEdges(1:end-1)+bin/2;
matchTimeEdges = timeBinEdges(~templMask);
reactTime = matchTimeEdges(1:end-1)+bin/2;
out.reactTime = reactTime;
% bin template spikes
spikesTemplate = binSpikes(:,templMask);
numBinsTemplate=size(spikesTemplate,2);
spikesMatch = binSpikes(:,~templMask);

% Z-transform binned firing of each cell
zSpikesTemplate=zscore(spikesTemplate')';
zSpikesMatch=zscore(spikesMatch')';
numMatchBins=size(spikesMatch,2);

%% create no-diag cc template mat
% create pairwise cell activity correlation matrix
zCCtemplate=(zSpikesTemplate*zSpikesTemplate')/numBinsTemplate; %effectively same as corrcoef(ZspikeBins')
zCCtemplateNoDiag = zCCtemplate - diag(diag(zCCtemplate)); % cancel out the diag i=j terms

%% Full Model version of Reactivation Strength
if fullModel
    reactFull=diag(zSpikesMatch'*zCCtemplateNoDiag*zSpikesMatch)';
%     reactFullNorm=(1/(2*numMatchBins))*sum(reactFull);
    out.reactFull = reactFull;
    %% all SWR psth
    d = setfiltertime(f, swrTimeFilter);
    swrTime = vec2list(~isExcluded(times, d.excludetime{1}{1}), times);
    out.swrTime = swrTime;
    [~, swrTimeIdx] = histc(swrTime(:,1), matchTimeEdges);
    idxWin = win(1)/bin:win(2)/bin;
    lag0idx = ceil(length(idxWin)/2);
    swrReactPSTHtime = idxWin*bin;
    swrTimeIdx = swrTimeIdx(swrTimeIdx+min(idxWin)>0);
    swrTimeIdx = swrTimeIdx(swrTimeIdx+max(idxWin)<length(reactFull));
    if ~isempty(swrTimeIdx)
        swrReactPSTHfull = cell2mat(arrayfun(@(r) reactFull(r+idxWin), swrTimeIdx, 'un', 0));
        swrReactETAfull = nanmean(swrReactPSTHfull);
        % all swr shuff
        swrReactETAfullShufs = [];
        for sh=1:numShuffs
            swrReactPSTHfullShuf=[];
            for qq=1:length(swrTimeIdx)
                shiftBy=round(rand(1)*size(swrReactPSTHfull,2));
                swrReactPSTHfullShuf(qq,:)=circshift(swrReactPSTHfull(qq,:),shiftBy,2);
            end
            swrReactETAfullShufs(sh,:) = nanmean(swrReactPSTHfullShuf);
        end
        swrReactETAfullShufsMean = nanmean(swrReactETAfullShufs);
        swrReactETAfullShufsSEM = nanstd(swrReactETAfullShufs)/...
            sqrt(size(swrReactETAfullShufs,1));
        %     swrP = ranksum(swrReactETAfullShufs(:,lag0idx), swrReactETAfull(:,lag0idx),'alpha',0.05,'tail','left');

    out.swrReactPSTHfull = swrReactPSTHfull;
    out.swrReactETAfull = swrReactETAfull;
    out.reactPSTHtime = swrReactPSTHtime;
    out.swrReactETAfullShufs = swrReactETAfullShufs;
    end
    
    %% SWR-Burst psth
    l = setfiltertime(f, burstTimeFilter);
    swrBurstTime = swrTime(~isExcluded(swrTime(:,1), l.excludetime{1}{1}),:);
    out.swrBurstTime = swrBurstTime;
    [~, swrBurstTimeIdx] = histc(swrBurstTime(:,1), matchTimeEdges);
    swrBurstTimeIdx = swrBurstTimeIdx(swrBurstTimeIdx+min(idxWin)>0);
    swrBurstTimeIdx = swrBurstTimeIdx(swrBurstTimeIdx+max(idxWin)<length(reactFull));
    if ~isempty(swrBurstTimeIdx)
        swrBurstReactPSTHfull = cell2mat(arrayfun(@(r) reactFull(r+idxWin), swrBurstTimeIdx, 'un', 0));
        swrBurstReactETAfull = nanmean(swrBurstReactPSTHfull);
        % swr-burst shuff
        swrBurstReactETAfullShufs = [];
        for sh=1:numShuffs
            swrBurstReactPSTHfullShuf=[];
            for qq=1:length(swrBurstTimeIdx)
                shiftBy=round(rand(1)*size(swrBurstReactPSTHfull,2));
                swrBurstReactPSTHfullShuf(qq,:)=circshift(swrBurstReactPSTHfull(qq,:),shiftBy,2);
            end
            swrBurstReactETAfullShufs(sh,:) = nanmean(swrBurstReactPSTHfullShuf);
        end
        swrBurstReactETAfullShufsMean = nanmean(swrBurstReactETAfullShufs);
        swrBurstReactETAfullShufsSEM = nanstd(swrBurstReactETAfullShufs)/...
            sqrt(size(swrBurstReactETAfullShufs,1));
        swrBurstP = ranksum(swrBurstReactETAfullShufs(:,lag0idx), swrBurstReactETAfull(:,lag0idx),'alpha',0.05,'tail','left');

        out.swrBurstReactPSTHfull = swrBurstReactPSTHfull;
        out.swrBurstReactETAfull = swrBurstReactETAfull;
        out.swrBurstReactETAfullShufs = swrBurstReactETAfullShufs;
    end
    
    %% Full Model lick reactivation psth, eta
    burstIntervals = vec2list(~isExcluded(times, l.excludetime{1}{1}), times);
    allLicks = lick{idx(1)}{idx(2)}.starttime;
    lickBurstTime = allLicks(logical(isExcluded(allLicks, burstIntervals)));
    % exclude licks within 1 second of a swr
    [~, d] = knnsearch(swrBurstTime(:), lickBurstTime);
    lickBurstTime = lickBurstTime(d > 1);
    lickBurstTime = lickBurstTime(lickBurstTime+min(idxWin)>0);
    lickBurstTime = lickBurstTime(lickBurstTime+max(idxWin)<length(reactFull));
    out.lickBurstTime = lickBurstTime;
    [~, lickBurstTimeIdx] = histc(lickBurstTime, matchTimeEdges);
    lickBurstTimeIdx = lickBurstTimeIdx(lickBurstTimeIdx+min(idxWin)>0);
    lickBurstTimeIdx = lickBurstTimeIdx(lickBurstTimeIdx+max(idxWin)<length(reactFull));
    if ~isempty(lickBurstTimeIdx)
        lickReactPSTHfull = cell2mat(arrayfun(@(r) reactFull(r+idxWin), lickBurstTimeIdx, 'un', 0));
        lickReactETAfull = nanmean(lickReactPSTHfull);
        % lick shuff
        lickReactETAfullShufs = [];
        for sh=1:numShuffs
            lickReactPSTHfullShuf=[];
            for qq=1:length(lickBurstTimeIdx)
                shiftBy=round(rand(1)*size(lickReactPSTHfull,2));
                lickReactPSTHfullShuf(qq,:)=circshift(lickReactPSTHfull(qq,:),shiftBy,2);
            end
            lickReactETAfullShufs(sh,:) = nanmean(lickReactPSTHfullShuf);
        end
        lickReactETAfullShufsMean = nanmean(lickReactETAfullShufs);
        lickReactETAfullShufsSEM = nanstd(lickReactETAfullShufs)/...
            sqrt(size(lickReactETAfullShufs,1));
        %         lickP = ranksum(lickReactETAfullShufs(:,lag0idx), lickReactETAfull(:,lag0idx),'alpha',0.05,'tail','left');
        out.lickReactPSTHfull = lickReactPSTHfull;
        out.lickReactETAfull = lickReactETAfull;
        out.lickReactETAfullShufs = lickReactETAfullShufs;
    end
end

%% PER PC
% matrix is symmetric so divide by 2 and number of bins for normalization
if perPC
    % decomposing correlation matrix into eigenvalues and eigenvectors
    [eigVec, eigVal] = eig(zCCtemplate); % now C=V*D*V'
    % Sort, so that the first explains most variance
    [eigValSort, eigValSortIdx] = sort(diag(eigVal),'descend');
    eigVecSort = eigVec(:,eigValSortIdx);
    
    % determining significant eigenvalues (as in Peyrache)
    eigValSigThresh = (1+sqrt(1/numBinsTemplate/numCells))^2;
    eigValSig = eigValSort(eigValSort>eigValSigThresh);
    eigVecSig = eigVecSort(:,eigValSort>eigValSigThresh);
    numEigValSig = length(eigValSig);
    
    for iEig = 1:numEigValSig
        PC = eigVecSig(:,iEig)*eigVecSig(:,iEig)'*eigValSig(iEig);
        PCsNoDiag(iEig,:,:) = PC - diag(diag(PC));
    end
    
    reactPerPC = [];
    for it=1:numMatchBins
        % outer product of the instantaneous ensemble firing
        % in each i,j is the product of firing of cell i and cell j (in this time bin).
        ensOutProd=zSpikesMatch(:,it)*zSpikesMatch(:,it)';
        for jEig=1:numEigValSig
            % dot-product ensembleOuterProd with corr matrix (whose diagonal has been 0'd)
            reactPerPC(jEig,it)=sum(sum(ensOutProd.*squeeze(PCsNoDiag(jEig,:,:))));
        end
    end
    out.reactPerPC = reactPerPC; % per pc full epoch 'match' intervals reactivation strength series
    %% per PC shuffle and test psth
    
    reactSignPCs = [];
    for iEig=1:numEigValSig
        % comparing reactivation strength to per PC circshift match zbinspikes
        reactPerPCShufmean=[];
        reactPerPCShuf=[];
        for qq=1:numShuffs
            zSpikesMatchShuf=[];
            for tr=1:numCells
                shiftBy=round(rand(1)*numMatchBins);
                zSpikesMatchShuf(tr,:)=circshift(zSpikesMatch(tr,:),shiftBy,2);
            end
            
            for i2=1:numMatchBins
                ensOutProdShuf=zSpikesMatchShuf(:,i2)*zSpikesMatchShuf(:,i2)';
                reactPerPCShuf(i2)=sum(sum(ensOutProdShuf.*squeeze(PCsNoDiag(iEig,:,:))));
            end
            reactPerPCShufmean =[reactPerPCShufmean mean(reactPerPCShuf)];
        end
        reactSignPC=mean(mean(reactPerPC(iEig,:))<reactPerPCShufmean);
        reactSignPCs{iEig}=[reactSignPCs reactSignPC];
        
        %% per PC psth and circshift test
        
        swrReactPSTHperPC = cell2mat(arrayfun(@(r) reactPerPC(iEig,r+idxWin), swrTimeIdx, 'un', 0));
        swrReactETAperPC = nanmean(swrReactPSTHperPC);
        % circshift test psth
        swrReactETAperPCShufsTrim=[];
        swrReactETAperPCShufs=[];
        for sh=1:numShuffs
            swrReactPSTHperPCShuf=[];
            for qq=1:length(swrTimeIdx)
                shiftBy=round(rand(1)*size(swrReactPSTHperPC,2));
                swrReactPSTHperPCShuf(qq,:)=circshift(swrReactPSTHperPC(qq,:),shiftBy,2);
                
            end
            swrReactETAperPCShuf=nanmean(swrReactPSTHperPCShuf);
            swrReactETAperPCShufs=[swrReactETAperPCShufs; swrReactETAperPCShuf];
            swrReactETAperPCShufsTrim=[swrReactETAperPCShufsTrim nanmean(swrReactETAperPCShuf(trimIdx:end-trimIdx))];
        end
        
        % determining significance of each PC's react PSTH with shuffling
        if nanmean(nanmean(swrReactETAperPC(trimIdx:end-trimIdx))>swrReactETAperPCShufsTrim)>0.95
            sigxcorr=1;
        else
            sigxcorr=0;
        end
        
        % need to add burst swr and lick to this psth,shuffle
        %         swrBurstReactPSTHperPC =
        %         swrBurstReactETAperPC =
        %         swrBurstReactPSTHperPCshuff =
        %         swrBurstReactETAperPCshuff =
        %
        %         lickReactPSTHperPC =
        %         lickReactETAperPC =
        %         lickReactPSTHperPCshuff =
        %         lickReactETAperPCshuff =
        
    end
end

end

function out = init_out(idx, animal, varargin)
out.idx = idx;
out.animal = animal;
out.su = [];
out.binSpikes = [];
out.spikes = [];
out.reactFull = []; % full epoch 'match' intervals reactivation strength series
out.reactTime = []; % full epoch 'match' time.. aka reactivation

out.swrTime = [];
out.swrReactPSTHfull = [];
out.swrReactETAfull = [];
out.reactPSTHtime = [];
out.swrReactETAfullShufs = [];

out.swrBurstTime = [];
out.swrBurstReactPSTHfull = [];
out.swrBurstReactETAfull = [];
out.swrBurstReactETAfullShufs = [];

out.lickBurstTime = [];
out.lickReactPSTHfull = [];
out.lickReactETAfull = [];
out.lickReactETAfullShufs = [];
end


