
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
numShuffS = 1000;
bin = 0.100; % seconds
if ~isempty(varargin)
    assign(varargin{:})
end
% init output
out = init_out(idx);
%%
su = evaluatefilter(cellinfo, cellfilter);
% Bin the spikes
epTimeRange = spikes{idx(1)}{idx(2)}{su(1,3)}{su(s,4)}.timerange;
timeBinEdges = epTimeRange(1):bin:epTimeRange(2);
for s = 1:length(su(:,1))
    iSpks = spikes{idx(1)}{idx(2)}{su(s,1)}{su(s,2)}.data;
    binSpikes(s,:) = histc(iSpks,timeBinEdges);
end

runBool = logical(isExcluded(timeBinEdges, templateIntervals));
fprintf('%.02f pct of ep time (%.03f s) as template\n',sum(runBool)/length(runBool)*100, ...
    sum(runBool)*bin)
tmpTimeEdges = timeBinEdges(runBool);
spikesTemplateTime = tmpTimeEdges(1:end-1)+bin/2;
matchTimeEdges = timeBinEdges(~runBool);
spikesMatchTime = matchTimeEdges(1:end-1)+bin/2;

spikesTemplate = binSpikes(:,runBool);
spikesMatch = binSpikes(:,~runBool);
%% now compute reactivation

% Z-transform firing of each cell and rotate so that rows are cells:
numBinsTemplate=size(spikesTemplate,2);
numCells=size(spikesTemplate,1);
if numCells<numCellsThresh
    fprintf('not enough cells:(%d) %d %d\n', numCells, idx(1), idx(2));
    return
end
zSpikesTemplate=zscore(spikesTemplate')';
% pairwise cell activity correlation matrix
zCorrCoef=(zSpikesTemplate*zSpikesTemplate')/numBinsTemplate; %effectively same as corrcoef(ZspikeBins')

% Put Cii=0 for all i will cancel out the i=j terms:
zCorrCoefNoDiag = zCorrCoef - diag(diag(zCorrCoef));
zSpikesMatch=zscore(spikesMatch')';
numMatchBins=size(spikesMatch,2);

% 'FULL MODEL': PRE-PRINCIPAL-COMPONENT DECOMPOSITION
reactSeries=diag(zSpikesMatch'*zCorrCoefNoDiag*zSpikesMatch)';
%  = (1/(2*numMatchBins))*sum(reactSeries);
% tempMatchTS = MtempMatch;
% tempMatch=[Mtemplatematch];
%% per swr, lick reactSeries psth, eta
andef = animaldef(animal);
f.animal = animaldef(animal);
f.epochs{1} = idx;
d = setfiltertime(f, swrTimeFilter);
times = epTimeRange(1):.001:epTimeRange(2);
swrIntervals = vec2list(~isExcluded(times, d.excludetime{1}{1}), times);
% use swr starts to check for inclusion into matchTimeEdges
[h, swrTimeIdx] = histc(swrIntervals(:,1), matchTimeEdges);
idxWin = win(1)/bin:win(2)/bin;
swrTimeIdx = swrTimeIdx(swrTimeIdx+min(idxWin)>0);
swrReactPSTH = cell2mat(arrayfun(@(r) reactSeries(r+idxWin), swrTimeIdx, 'un', 0));
swrReactETA = mean(swrReactPSTH);
% plot(idxWin, swrReactETA)
%% per lick reactSeries psth, eta
l = setfiltertime(f, lickTimeFilter);
burstIntervals = vec2list(~isExcluded(times, l.excludetime{1}{1}), times);
allLicks = lick{idx(1)}{idx(2)}.starttime;
licksFilt = allLicks(logical(isExcluded(allLicks, burstIntervals)));

[h, burstTimeIdx] = histc(licksFilt, matchTimeEdges);
idxWin = win(1)/bin:win(2)/bin;
burstTimeIdx = burstTimeIdx(burstTimeIdx+min(idxWin)>0);
lickReactPSTH = cell2mat(arrayfun(@(r) reactSeries(r+idxWin), burstTimeIdx, 'un', 0));
lickReactETA = mean(lickReactPSTH);
% hold on;
% plot(idxWin, lickReactETA)
% hold off;
%% full model shuffle
if 0
    % Test significance of reactivation by comparing reactivation strength to that from
    % shuffled Qmatch matrices
%     MtemplatematchShufs=[];
    reactSeriesShuf = [];
    for qq=1:numShuffS
        zSpikesMatchShuf=[];
        for ic=1:numCells
            shiftamount=round(rand(1)*numMatchBins);
            zSpikesMatchShuf(ic,:)=circshift(zSpikesMatchShuf(ic,:),shiftamount,2);
        end
        reactSeriesShuf(qq,:) =diag(zSpikesMatchShuf'*C0*zSpikesMatchShuf);
%         MtemplatematchShuf=(1/(2*numMatchBins))*sum(reactSeriesShuf);
%         MtemplatematchShufs=[MtemplatematchShufs MtemplatematchShuf];
    end
    % =mean(Mtemplatematch<MtemplatematchShufs);
    % =[reactSeries reactSeriesFullModel];
end

%% PER PC PRE-PRINCIPAL-COMPONENT DECOMPOSITION
% For each time point, I do outer product of the instantaneous ensemble firing
% vector, resulting in a matrix where in each i,j is the product of firing
% of cell i and cell j (in this time bin). If we dot-product this with the
% correlation matrix (whose diagonal has been 0'd), we get for each pair of
% cells i,j, their firing product, times their correlation, which is what
% we want. Since the matric is symmetric we have to divide by 2, and also
% by the number of bins for normalization (this last one not entirely clear
% to me but that's just a constant).

% decomposing correlation matrix into eigenvalues and eigenvectors
[eigVec, eigVal] = eig(zCorrCoef); % now C=V*D*V'
% Sorting in descending order, so that now the first is the one explaining most variance
[eigValSort, eigValSortIdx] = sort(diag(eigVal),'descend');
eigVecSort = eigVec(:,eigValSortIdx);

% determining significant eigenvalues (as in Peyrache)
eigValSigThresh = (1+sqrt(1/numBins/numCells))^2;
eigValSig = eigValSort(eigValSort>eigValSigThresh);
eigVecSig = eigVecSort(:,eigValSort>eigValSigThresh);
numEigValSig = length(eigValSig);

for i=1:numEigValSig
    curP=eigVecSig(:,i)*eigVecSig(:,i)'*eigValSig(i);
    % again, I think putting 0's on the diagonal cancels out the i=j terms
    % doing it for all P matrices
    curP = curP - diag(diag(curP));
    allPs(i,:,:)=curP;

end
% Creating Mtemplatematches, which holds for each time bin (column) the
% projection on the i'th component (rows)
reactSeriesPerPC = [];
for i=1:numMatchBins
    ensembleOuterProd=zSpikesMatch(:,i)*zSpikesMatch(:,i)';
    for j=1:numEigValSig
        reactSeriesPerPC(j,i)=sum(sum(ensembleOuterProd.*squeeze(allPs(j,:,:))));
    end
end

%% for each sig PC, test its significance, then test mean riptrig react significance

allcrosscorrs=[];
for ii=1:numsigeigenvals

    % Test significance of reactivation by comparing reactivation strength to that from
    % shuffled Qmatch matrices

    MtemplatematchShufsPC1=[];
    MtemplatematchTimeSeriesPCShuf=[];
    for qq=1:numShuffleRuns
        QmatchShuf=[];
        for tr=1:numCells
            shiftamount=round(rand(1)*numMatchBins);
            QmatchShuf(tr,:)=circshift(Qmatch(tr,:),shiftamount,2);
        end

        for i2=1:numMatchBins

            ensembleOuterProdShuf=QmatchShuf(:,i2)*QmatchShuf(:,i2)';

            MtemplatematchTimeSeriesPCShuf(i2)=sum(sum(ensembleOuterProdShuf.*squeeze(allPs(ii,:,:))));
        end

        MtemplatematchShufsPC1=[MtemplatematchShufsPC1 mean(MtemplatematchTimeSeriesPCShuf)];

    end

    reactivationSignifPC=mean(mean(MtemplatematchesTimeSeriesPC(ii,:))<MtemplatematchShufsPC1);
    allallreactivationSignifPC{ii}=[allallreactivationSignifPC{ii} reactivationSignifPC];

  % creating SWR-triggered reactivation
    ripindsPC=ripinds(ripinds-maxlaginbins>0&ripinds+maxlaginbins<size(MtemplatematchesTimeSeriesPC(ii,:),2));
    riptrigreactmat=[];
    numrips1=length(ripindsPC);
    for qq=1:numrips1
        riptrigreactmat=[riptrigreactmat;MtemplatematchesTimeSeriesPC(ii,ripindsPC(qq)-maxlaginbins:ripindsPC(qq)+maxlaginbins)];
    end

    currawcrosscorr=mean(riptrigreactmat)';


    % determining significant cross-corr by shuffling
    limit1=45;
    allshufmaxes=[];
    allmeanshufs=[];

    for tt=1:numShuffleRuns
        riptrigreactmatshuf=[];
        for qq=1:numrips1
            shiftamount=round(rand(1)*size(riptrigreactmat,2));
            riptrigreactmatshuf(qq,:)=circshift(riptrigreactmat(qq,:),shiftamount,2);

        end
        meanshuf=mean(riptrigreactmatshuf);
        allmeanshufs=[allmeanshufs;meanshuf];
        allshufmaxes=[allshufmaxes mean(meanshuf(limit1:end-limit1))];
    end
    if mean(mean(currawcrosscorr(limit1:end-limit1))>allshufmaxes)>0.95
        sigxcorr=1;
    else
        sigxcorr=0;
    end


%                 figure;errorbar(1:size(riptrigreactmat,2),mean(riptrigreactmat),std(riptrigreactmat)/sqrt(size(riptrigreactmat,1)))
%                 axis([41 61 0 0.001])
%                 title(['ind= ' num2str(curSleepCellInds(1,1:3)) ' numcells= ' num2str(numCells) ' eigvalnum= ' num2str(ii) ' sig= ' num2str(sigxcorr)])
%                 saveas(gcf, [savedirimgsX 'reactrip-'   num2str(curSleepCellInds(1,1:3)) 'eigval' num2str(ii) '.jpg'])

    allallcrosscorrs{ii}=[allallcrosscorrs{ii};currawcrosscorr'];
    allallcrosscorrsShufs{ii}=[allallcrosscorrsShufs{ii};allmeanshufs];

    allallanimdays{ii}=[allallanimdays{ii};curSleepCellInds(1,1:3)];
    allallsigxcorr{ii}=[allallsigxcorr{ii};sigxcorr];
    allallnumcells{ii}=[allallnumcells{ii}; numCells];
    allallnumbins{ii}=[allallnumbins{ii}; numBins];
    allallnummatchbins{ii}=[allallnummatchbins{ii};numMatchBins];
    allalltemplatematchesTimeSeriesPC{ii}{end+1}=MtemplatematchesTimeSeriesPC(ii,:);

end


%% collect lick and swr react psth from full series
burstReactSWRPSTH = [];
nonBurstReactSWRPSTH = [];
timeSWRPSTH = [];
burstReactLickPSTH = [];
timelickPSTH = [];

%% plot epoch
% plot burst swr-psth, burst lick-psth, nonburst swr-psth

% also plot mean, std, shuffle sign bounds
% plot binned distr a la peyrach, gideon

% plot the representation of the significant PC's, a la Peyrache

%% return output
out.reactSeries = reactSeries; % full epoch 'match' intervals reactivation strength series
out.reactSeriesPerPC = reactSeriesPerPC; % per pc full epoch 'match' intervals reactivation strength series
out.timeReactSeries = timeBinEdges(runBool); % full epoch 'match' intervals time

out.burstReactSWRPSTH = burstReactSWRPSTH; % nBurstSWRDin x periSWRTimeWin
out.nonBurstReactSWRPSTH = nonBurstReactSWRPSTH; % nNonBurstSWRDin x periSWRTimeWin
out.timeSWRPSTH = timeSWRPSTH; % time

out.burstReactLickPSTH = burstReactLickPSTH;% nLickDin x periLickTimeWin
out.timelickPSTH = timelickPSTH;
end

function out = init_out(idx, varargin)
out.reactSeries = []; % full epoch 'match' intervals reactivation strength series
out.reactSeriesPerPC = []; % per pc full epoch 'match' intervals reactivation strength series
out.timeReactSeries = []; % full epoch 'match' intervals time

out.burstReactSWRPSTH = []; % nBurstSWRDin x periSWRTimeWin
out.nonBurstReactSWRPSTH = []; % nNonBurstSWRDin x periSWRTimeWin
out.timeSWRPSTH = []; % time

out.burstReactLickPSTH = [];% nLickDin x periLickTimeWin
out.timelickPSTH = [];
end


