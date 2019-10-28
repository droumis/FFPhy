
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
numCellsThresh = 10;
savefigas = {'png'};
numShuffs = 100;
bin = 0.100; % seconds
fullModel = 1;
perPC = 1;
plotfigs = 1;
savefigs = 1;
trimIdx = 4;
if ~isempty(varargin)
    assign(varargin{:})
end

% init output
out = init_out(idx, animal);

%% Filter and bin All su spikes
su = evaluatefilter(cellinfo, cellfilter);
numCells=size(su,1);
if numCells<numCellsThresh
    fprintf('not enough cells:(%d) %d %d\n', numCells, idx(1), idx(2));
    return
end
epTimeRange = spikes{idx(1)}{idx(2)}{su(1,3)}{su(s,4)}.timerange;
timeBinEdges = epTimeRange(1):bin:epTimeRange(2);
for s = 1:length(su(:,1))
    iSpks = spikes{idx(1)}{idx(2)}{su(s,1)}{su(s,2)}.data;
    binSpikes(s,:) = histc(iSpks,timeBinEdges);
end
%% Get Template and Match binned zSpikes
templMask = logical(isExcluded(timeBinEdges, templateIntervals));
fprintf('template pct eptime %.02f (%.03f s)\n',sum(templMask)/length(templMask)*100, ...
    sum(templMask)*bin)

tmpTimeEdges = timeBinEdges(templMask);
spikesTemplateTime = tmpTimeEdges(1:end-1)+bin/2;
matchTimeEdges = timeBinEdges(~templMask);
spikesMatchTime = matchTimeEdges(1:end-1)+bin/2;

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
reactFull=diag(zSpikesMatch'*zCCtemplateNoDiag*zSpikesMatch)';
Mtemplatematch=(1/(2*numMatchBins))*sum(reactFull);
%  = (1/(2*numMatchBins))*sum(reactSeries);
% tempMatchTS = MtempMatch;
% tempMatch=[Mtemplatematch];

%% full model time-shuffle
if fullModel
    % Test significance of reactivation by comparing reactivation strength to that from
    % shuffled Qmatch matrices
    %     MtemplatematchShufs=[];
    reactFullShuf = [];
    MtemplatematchShufs = [];
    for qq = 1:numShuffs
        zSpikesMatchShuf = [];
        for ic = 1:numCells
            shiftamount = round(rand(1)*numMatchBins);
            zSpikesMatchShuf(ic,:) = circshift(zSpikesMatch(ic,:),shiftamount,2);
        end
        MtemplatematchTimeSeriesShuf =diag(zSpikesMatchShuf'*zCCtemplateNoDiag*zSpikesMatchShuf);
        MtemplatematchShuf=(1/(2*numMatchBins))*sum(MtemplatematchTimeSeriesShuf);
        MtemplatematchShufs=[MtemplatematchShufs MtemplatematchShuf];
    end
    reactSignFull=nanmean(Mtemplatematch<MtemplatematchShufs);
    
    %% Full Model swr reactivation psth, eta
    andef = animaldef(animal);
    f.animal = animaldef(animal);
    f.epochs{1} = idx;
    d = setfiltertime(f, swrTimeFilter);
    times = epTimeRange(1):.001:epTimeRange(2);
    swrIntervals = vec2list(~isExcluded(times, d.excludetime{1}{1}), times);
    % use swr starts to check for inclusion into matchTimeEdges
    [h, swrTimeIdx] = histc(swrIntervals(:,1), matchTimeEdges);
    idxWin = win(1)/bin:win(2)/bin;
    swrReactPSTHtime = idxWin*.1;
    swrTimeIdx = swrTimeIdx(swrTimeIdx+min(idxWin)>0);
    swrReactPSTH = cell2mat(arrayfun(@(r) reactFull(r+idxWin), swrTimeIdx, 'un', 0));
    swrReactETA = nanmean(swrReactPSTH);
    % plot(idxWin, swrReactETA)
    
    % swr psth shuff
    swrReactPSTHshuff =
    swrReactETAshuff =
    
    %% Burst swr psth, eta
    swrBurstReactPSTH =
    swrBurstReactETA =
    swrBurstReactPSTHshuff =
    swrBurstReactETAshuff =
    
    %% Full Model lick reactivation psth, eta
    l = setfiltertime(f, lickTimeFilter);
    burstIntervals = vec2list(~isExcluded(times, l.excludetime{1}{1}), times);
    allLicks = lick{idx(1)}{idx(2)}.starttime;
    licksFilt = allLicks(logical(isExcluded(allLicks, burstIntervals)));
    
    [h, burstTimeIdx] = histc(licksFilt, matchTimeEdges);
    idxWin = win(1)/bin:win(2)/bin;
    burstTimeIdx = burstTimeIdx(burstTimeIdx+min(idxWin)>0);
    lickReactPSTH = cell2mat(arrayfun(@(r) reactFull(r+idxWin), burstTimeIdx, 'un', 0));
    lickReactETA = mean(lickReactPSTH);
    % plot(idxWin, lickReactETA)
    
    lickReactPSTHshuff =
    lickReactETAshuff =
    %% stat testing
    % p-value at lag 0
    lag0idx = [];
    swrP = ranksum(swrReactETAshuff(:,lag0idx), swrReactETA(:,lag0idx),'alpha',0.05,'tail','left');
    swrBurstP = ranksum(swrBurstReactETAshuff(:,lag0idx), swrBurstReactETA(:,lag0idx),'alpha',0.05,'tail','left');
    lickP = ranksum(lickReactETAshuff(:,lag0idx), lickReactETA(:,lag0idx),'alpha',0.05,'tail','left');
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
    %% per PC shuffle and test psth
    
    reactSignPCs = [];
    for iEig=1:numEigValSig
        % comparing reactivation strength to per PC circshift match zbinspikes
        reactPerPCShufmean=[];
        reactPerPCShuf=[];
        for qq=1:numShuffs
            zSpikesMatchShuf=[];
            for tr=1:numCells
                shiftamount=round(rand(1)*numMatchBins);
                zSpikesMatchShuf(tr,:)=circshift(zSpikesMatch(tr,:),shiftamount,2);
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
        for tt=1:numShuffs
            swrReactPSTHperPCShuf=[];
            for qq=1:length(swrTimeIdx)
                shiftamount=round(rand(1)*size(swrReactPSTHperPC,2));
                swrReactPSTHperPCShuf(qq,:)=circshift(swrReactPSTHperPC(qq,:),shiftamount,2);
                
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
if plotfigs
    %% plot Full ETA w sem errorbars
    % plot swr ETA with ebars, shuff
    
    subplot(3,1,1)
    plot(swrReactPSTHtime, swrReactETA,'k','linewidth',3);
    hold on;
    plot(swrReactPSTHtime,swrReactETAshuff,'r','linewidth',3);
    errorbar(swrReactPSTHtime,swrReactETA,nanstd(swrReactPSTH)/sqrt(size(swrReactPSTH,1)),'k','linewidth',0.5)
    errorbar(swrReactPSTHtime,swrReactETAshuff,nanstd(swrReactETAshuff)/sqrt(size(swrReactETAshuff,1)),'r','linewidth',0.5)
    ylabel('Reactivation strength')
    xlabel('Time (s)')
    title('swrReactETA')
    
    % plot swr burst ETA with ebars, shuff
    subplot(3,1,2)
    plot(swrReactPSTHtime, swrBurstReactETA,'k','linewidth',3);
    hold on;
    plot(swrReactPSTHtime,swrBurstReactETAshuff,'r','linewidth',3);
    errorbar(swrReactPSTHtime,swrBurstReactETA, ...
        nanstd(swrBurstReactPSTH)/sqrt(size(swrBurstReactPSTH,1)),'k','linewidth',0.5)
    errorbar(swrReactPSTHtime,swrBurstReactETAshuff, ...
        nanstd(swrBurstReactPSTHshuff)/sqrt(size(swrBurstReactPSTHshuff,1)),'r','linewidth',0.5)
    ylabel('Reactivation strength')
    xlabel('Time (s)')
    title('swrBurstReactETA')
    
    % plot lick ETA with ebars, shuff
    subplot(3,1,3)
    plot(swrReactPSTHtime, lickReactETA,'k','linewidth',3);
    hold on;
    plot(swrReactPSTHtime,lickReactETAshuff,'r','linewidth',3);
    errorbar(swrReactPSTHtime,lickReactETA,nanstd(lickReactPSTH)/sqrt(size(lickReactPSTH,1)),'k','linewidth',0.5)
    errorbar(swrReactPSTHtime,lickReactETAshuff,nanstd(lickReactETAshuff)/sqrt(size(lickReactETAshuff,1)),'r','linewidth',0.5)
    ylabel('Reactivation strength')
    xlabel('Time (s)')
    title('lickReactETA')
    
    %% plot time snippet of spikes + zbinSpikes traces + reactivation result (full,pc)
    exWin = [];
    % plot spikes in time window
    
    % plot zbinspikes traces in time window
    % plot like lfp
    
    % plot reactivation result in time window
    % reactPerPC (1 trace per sig pc)
    % reactFull
    
    %% plot per PC method Demo
    % plot the steps and examples of the main parts
    % templateCC Eigvec, val, sig PC's, match z spikes, per pc react
    if 1
        figure
        subplot(2,3,1); imagesc(zCCtemplateNoDiag); ax = gca; ax.YDir = 'normal'; title('zCCtemplateNoDiag')
        ylabel('cell')
        subplot(2,3,2); imagesc(eigVal); ax = gca; ax.YDir = 'normal'; title('eigval')
        subplot(2,3,3); imagesc(eigVec); ax = gca; ax.YDir = 'normal'; title('eigvec')
        for iEig = 1:numEigValSig
            subplot(2,3,iEig+3); imagesc(squeeze(PCsNoDiag(iEig,:,:))); ax = gca; ax.YDir = 'normal';
            title(sprintf('zCCtemplate SigPC%dnoDiag', iEig))
        end
    end
    %% Plot reactStrength swr-psth, lick-psth, burst-swr-psth
    % plot burst swr-psth, burst lick-psth, nonburst swr-psth
    % underneath heatrasters, plot ETA, stdfill, shuffle
    
    %%
    % plot the spatial representation of the significant PC's, a la Peyrache
end
%% return output
out.reactSeries = reactFull; % full epoch 'match' intervals reactivation strength series
out.reactSeriesPerPC = reactPerPC; % per pc full epoch 'match' intervals reactivation strength series
out.timeBinEdges = timeBinEdges(templMask); % full epoch 'match' intervals time

out.swrReactPSTH = swrReactPSTH;
out.swrReactETA = swrReactETA;
out.swrReactPSTHtime = swrReactPSTHtime;

out.swrReactPSTHshuff = swrReactPSTHshuff;
out.swrReactETAshuff = swrReactETAshuff;

out.swrBurstReactPSTH = swrBurstReactPSTH; % nBurstSWRDin x periSWRTimeWin
out.swrNoBurstReactPSTH = swrNoBurstReactPSTH; % nNonBurstSWRDin x periSWRTimeWin

out.singleepochReactPSTH = lickReactPSTH;% nLickDin x periLickTimeWin
out.lickReactPSTHtime = [];
end

function out = init_out(idx, animal, varargin)
out.idx = idx;
out.idx = animal;
out.reactSeries = [];
out.reactSeriesPerPC = [];
out.timeBinEdges = [];

out.swrReactPSTH = [];
out.swrReactETA = [];
out.swrReactPSTHtime = [];

out.swrReactPSTHshuff = [];
out.swrReactETAshuff = [];

out.swrBurstReactPSTH = [];
out.swrNoBurstReactPSTH = [];

out.lickReactPSTH = [];
out.lickReactPSTHtime = [];
end


