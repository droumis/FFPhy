


%{

- bin >4cm/s spikes into 100 ms intervals

%}

animal = 'D10';
day = 6;
epoch = 2;

numShuffleRuns=100;
numcellthresh=4;
velThresh = 4; %cm/s
%% load data
andef = animaldef(animal);
spikes = loaddatastruct(andef{2}, animal, 'spikes', day);
pos = loaddatastruct(andef{2}, animal, 'pos', day);
lick = loaddatastruct(andef{2}, animal, 'lick', day);
ripple = loaddatastruct(andef{2}, animal, 'ca1rippleskons', day);
%%
% get a list of nonmu spike indices.
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'reactivationPLTH';
Fp.params = {'>4cm/s', 'ca1SU', 'wtrackdays', 'exemplar_wepochs', Fp.filtfunction};
Fp = load_filter_params(Fp);
%%
F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
    'excludetime', Fp.timefilter,'cells', Fp.cellfilter, 'iterator',Fp.iterator);
%%
Fr.params = {'ripples'};
Frip = load_filter_params(Fr);
Frip = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
    'excludetime', Frip.timefilter,'cells', Fp.cellfilter, 'iterator',Fp.iterator);
%%
ani = 1;
su = F(ani).data{1}{1};
spikesBinned = [];
s = 1;
epTimeRange = spikes{day}{epoch}{su(s,1)}{su(s,2)}.timerange;
binnedTimeEdges = epTimeRange(1):0.100:epTimeRange(2);
for s = 1:length(su(:,1))
    iSpks = spikes{day}{epoch}{su(s,1)}{su(s,2)}.data;
    spikesBinned(s,:) = histc(iSpks,binnedTimeEdges);
end

%% filter spikes: template (>4cm/s), match (<4cm/s)
andef = animaldef(animal);
out = get2dstate(andef{2},animal,[day epoch]);
vel = out{day}{epoch}.velocity;
velBool = vel > velThresh;
runIntervals = vec2list(~velBool, out{day}{epoch}.time);

runBool = ~isExcluded(binnedTimeEdges, runIntervals);

spikesTemplateTime = binnedTimeEdges(runBool);
spikesMatchTime = binnedTimeEdges(~runBool);

spikesTemplate = spikesBinned(:,runBool);
spikesMatch = spikesBinned(:,~runBool);
%%

% Z-transform firing of each cell and rotate so that rows are cells:
numBins=size(spikesTemplate,2);
numCells=size(spikesTemplate,1);
if numCells>=numcellthresh
    ZspikeBins=zscore(spikesTemplate')';
    % correlation matrix
    C=(ZspikeBins*ZspikeBins')/numBins; %effectively same as corrcoef(ZspikeBins')
    
    % Put Cii=0 for all i will cancel out the i=j terms:
    C0 = C - diag(diag(C));
    
    Qmatch=zscore(spikesMatch')';
    numMatchBins=size(spikesMatch,2);
    
    %===== 'FULL MODEL': PRE-PRINCIPAL-COMPONENT DECOMPOSITION
    % I do the same analyses below per- principal component
    MtempMatch=diag(Qmatch'*C0*Qmatch)';
    Mtemplatematch=(1/(2*numMatchBins))*sum(MtempMatch);
    
    tempMatch{end+1}=MtempMatch;
    tempMatch=[tempMatch Mtemplatematch];
    
    % Test significance of reactivation by comparing reactivation strength to that from
    % shuffled Qmatch matrices
    
    MtemplatematchShufs=[];
    for qq=1:numShuffleRuns
        QmatchShuf=[];
        for tr=1:numCells
            shiftamount=round(rand(1)*numMatchBins);
            QmatchShuf(tr,:)=circshift(Qmatch(tr,:),shiftamount,2);
        end
        
        MtemplatematchTimeSeriesShuf=diag(QmatchShuf'*C0*QmatchShuf);
        MtemplatematchShuf=(1/(2*numMatchBins))*sum(MtemplatematchTimeSeriesShuf);
        MtemplatematchShufs=[MtemplatematchShufs MtemplatematchShuf];
    end
    
    reactivationSignifFullModel=mean(Mtemplatematch<MtemplatematchShufs);
    allallreactivationSignificanceFullModel=[allallreactivationSignificanceFullModel reactivationSignifFullModel];
end
%% Collect lick and swr reactivation PSTH

% use the lick and swr as psth events instead of as filters on spikes
lickvec = getLickBout([], animal, epochs, varargin);
lickBurstIntervals = vec2list(lickvec{day}{epoch}.lickBout, lickvec{day}{epoch}.time);
spikesSWRMatch =

%  ripple-triggered psth of reactivation.

ripvec=zeros(1,numMatchBins);
ripvec(allsleepspikebins{sleepind1(cursleepep)}.rippleBins)=1;
ripinds=find(ripvec);
ripindsFullModel=ripinds(ripinds-maxlaginbins>0&ripinds+maxlaginbins<numMatchBins);
riptrigreactmat=[];
numrips1=length(ripindsFullModel);
for qq=1:numrips1
    riptrigreactmat=[riptrigreactmat;MtemplatematchTimeSeries(ripindsFullModel(qq)-maxlaginbins:ripindsFullModel(qq)+maxlaginbins)];
end

%% plot
ripbounds = Frip.excludetime{1}{1}(:,2);
plot(spikesMatchTime, MtempMatch)
hold on
line([ripbounds'; ripbounds'], [ones(1,length(ripbounds)); 100*ones(1,length(ripbounds))], 'color', [.5 .5 .5 .1])
hold off

