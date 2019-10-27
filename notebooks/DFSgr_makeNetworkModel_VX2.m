
% previously called DFSgr_makeNetworkModel_V8groupAwakes.m (changed Nov. 6th, 2015)

%clear; %close all;
runscript = 0;
savedata = 1; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells
plotGraphs=1;
savedirX = '/opt/data15/gideon/GR_ProcessedData/';
savedirimgsX=['~/Code/Matlab/analysis/stage2analysis/Dropbox/reactrip/'];

%val=1; savefile = [savedirX 'ACnetwork_NoMove'];
%val=2; savefile = [savedirX 'ACnetwork_Move'];
%val=3; savefile = [savedirX 'ACnetwork_FastMove'];
%val=4; savefile = [savedirX 'ACnetwork_SlowMove'];
val=5; savefile = [savedirX 'ACnetwork_AllEpoch'];
%val=6; savefile = [savedirX 'PFCnetwork_AllEpoch'];
%val=7; savefile = [savedirX 'ACnetwork_AllEpochSS'];



savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    animals = {'Parker','Rosenthal','Nadal','Borg'};
    
    
    runepochfilter = 'isequal($type, ''run'')';
    sleepepochfilter ='isequal($type, ''sleep'') && isequal($audprot, ''0'')';
    
    if val<=5
    cellfilter = 'strcmp($area, ''AC'') && ($numspikes > 100)'; %
    else
    cellfilter = 'strcmp($area, ''PFC'') && ($numspikes > 100)'; %
        
    end
    
    
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    
    timefilter_sleep = {{'kk_getsleep', '(($sleep == 1))'},{'DFTFsj_getvelpos', '(($absvel < 10000))'}};
    timefilter_nrem = {{'kk_getnrem', '(($nrem == 1))'},{'DFTFsj_getvelpos', '(($absvel < 10000))'}};
    timefilter_rem = {{'kk_getrem', '(($rem == 1))'},{'DFTFsj_getvelpos', '(($absvel < 10000))'}};
    if val==1
        timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel < 4))'},...
            {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    elseif val==2
        timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel >= 4))'},...
            {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    elseif val==3
        timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel > 12))'},...
            {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    elseif val==4
        timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel >=4 & $absvel < 12))'},...
            {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    elseif val==5|val==6
        timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel >=0))'},...
            {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
        
    end
    
    
    
    % Iterator
    % --------
    iterator = 'multicellanal';
    
    % Filter creation
    % ----------------
    modg = createfilter('animal',animals,'epochs',sleepepochfilter, 'cells',...
        cellfilter, 'excludetime',timefilter_nrem,'iterator', iterator);
    %
    
    modf = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
        cellfilter,  'excludetime', timefilter_place_new, 'iterator', iterator);
    
    
    disp('Done Filter Creation');
    
    
    if val==7
    modg = setfilterfunction(modg,'DFAgr_getBinnedDataSeparatedSWRs',{'spikes','ripples','sound','cellinfo','tetinfo','pos'},'thrstime',1); % With includetime condition
    modf = setfilterfunction(modf,'DFAgr_getBinnedDataSeparatedSWRs',{'spikes','ripples','sound','cellinfo','tetinfo','pos'},'thrstime',1); % With includetime condition
    else
    modg = setfilterfunction(modg,'DFAgr_getBinnedData',{'spikes','ripples','sound','cellinfo','tetinfo','pos'},'thrstime',1); % With includetime condition
    modf = setfilterfunction(modf,'DFAgr_getBinnedData',{'spikes','ripples','sound','cellinfo','tetinfo','pos'},'thrstime',1); % With includetime condition
        
    end
    
    % Run analysis
    % ------------

    modg = runfilter(modg);
    modf = runfilter(modf);
    
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save(savefile);
    end
    
else
    
    load(savefile);
    
end % end runscript

if ~exist('savedata')
    return
end


% -------------------------  Filter Format Done -------------------------

% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------
gatherdata = 0; savegatherdata = 1;
switch val
    case 1
        gatherdatafile = [savedirX 'ACnetwork_NoMove_gather'];
    case 2
        gatherdatafile = [savedirX 'ACnetwork_Move_gather'];
    case 3
        gatherdatafile = [savedirX 'ACnetwork_FastMove_gather'];
    case 4
        gatherdatafile = [savedirX 'ACnetwork_SlowMove_gather'];
    case 5
        gatherdatafile = [savedirX 'ACnetwork_AllEpoch_gather'];
    case 6
        gatherdatafile = [savedirX 'PFCnetwork_AllEpoch_gather'];
    case 7
        gatherdatafile = [savedirX 'ACnetwork_AllEpochSS_gather'];
 
end

numShuffleRuns=100;
numcellthresh=4;
tic

if gatherdata
    plot1=0;
    % Gathering for AWAKE
    cnt=0;
    allwakeinds=[];
    allwakespikebins={};
    allwakesoundbins={};
    allwakespeedbins={};
    allwakeanimdayep=[];
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            curSpikeBinsTmp=modf(an).output{1}(i).spikeBins;
            if ~isempty(curSpikeBinsTmp)
                cnt=cnt+1;
                
                indNoAnim=modf(an).output{1}(i).indices;
                indWAnim=[an*ones(size(indNoAnim,1),1) indNoAnim];
                allwakeinds=[allwakeinds;indWAnim];
                allwakeanimdayep=[allwakeanimdayep;indWAnim(1,1:3)];
                
                allwakespikebins{cnt}.inds=indWAnim;
                allwakespikebins{cnt}.spikeBins=modf(an).output{1}(i).spikeBins;
                allwakesoundbins{cnt}=modf(an).output{1}(i).soundBins;
                allwakespeedbins{cnt}=modf(an).output{1}(i).speedBins;
            end
        end
        
    end
    
    % Gathering for SLEEP
    
    cnt=0;
    allsleepinds=[];
    allsleepspikebins={};
    allsleepanimdayep=[];
    for an = 1:length(modg)
        for i=1:length(modg(an).output{1})
            curSpikeBinsTmp=modg(an).output{1}(i).spikeBins;
            if ~isempty(curSpikeBinsTmp)
                
                cnt=cnt+1;
                
                indNoAnim=modg(an).output{1}(i).indices;
                indWAnim=[an*ones(size(indNoAnim,1),1) indNoAnim];
                
                allsleepinds=[allsleepinds;indWAnim];
                allsleepanimdayep=[allsleepanimdayep;indWAnim(1,1:3)];
                
                allsleepspikebins{cnt}.inds=indWAnim;
                allsleepspikebins{cnt}.spikeBins=modg(an).output{1}(i).spikeBins;
                allsleepspikebins{cnt}.rippleBins=modg(an).output{1}(i).rippleBins;
            end
            
        end
        
    end
    allallcrosscorrs={};for k=1:20,allallcrosscorrs{k}=[];end
    allallcrosscorrsShufs={};for k=1:20,allallcrosscorrsShufs{k}=[];end
    allallanimdays={};for k=1:20,allallanimdays{k}=[];end
    allallsigxcorr={};for k=1:20,allallsigxcorr{k}=[];end
    allallnumcells={};for k=1:20,allallnumcells{k}=[];end
    allallnumbins={};for k=1:20,allallnumbins{k}=[];end
    allallnummatchbins={};for k=1:20,allallnummatchbins{k}=[];end
    allallmeantemplatematches={};for k=1:20,allallmeantemplatematches{k}=[];end
    allalltemplatematchesTimeSeriesPC={};for k=1:20,allalltemplatematchesTimeSeriesPC{k}=[];end
    allallanimdayspertime={};for k=1:20,allallanimdayspertime{k}=[];end
    allallreactivationSignifPC={};for k=1:20,allallreactivationSignifPC{k}=[];end
    
    
    allalltemplatematchTimeSeries={};
    allallripplebins={};
    allalltemplatematch=[];
    allallreactivationSignificanceFullModel=[];
    allallcrosscorrsFullModel=[];
    allallcrosscorrsFullModelShuf=[];
    
    allallanimdaysFullModel=[];
    allallsigxcorrFullModel=[];
    allallnumcellsFullModel=[];
    allallnumbinsFullModel=[];
    allallnummatchbinsFullModel=[];
    
    % in the previous version of this code, I looked at consecutive pairs
    % of awake and sleep epochs. Here I concatenate awake epochs, calculate
    % correlation models on them together, and observe how reactivation
    % across the rest sessions develops throughout the day
    uniqueawakeanimdays=unique(allwakeanimdayep(:,1:2),'rows');
    binsize=0.2;
    xcorrmaxlaginS=10;
    maxlaginbins=xcorrmaxlaginS/binsize;
    task2eps=[1 13;1 14;1 15;2 10;2 11;2 12];
    
    for wakeind1=1:size(uniqueawakeanimdays,1)
        wakeind1
        curanimday=uniqueawakeanimdays(wakeind1,:);
        if isempty(find(ismember(task2eps,curanimday,'rows')))
            
            allawakedailyepochsinds=find(ismember(allwakeanimdayep(:,1:2),curanimday,'rows'));
            for curep=1:length(allawakedailyepochsinds)
                curepcellinds=allwakespikebins{allawakedailyepochsinds(curep)}.inds(:,[1 2 4 5]);
                curepspikebins=allwakespikebins{allawakedailyepochsinds(curep)}.spikeBins;
                
                if curep==1
                    consensuscellinds=curepcellinds;
                    consensusspikebins=curepspikebins;
                else
                    [consensuscellinds iCurcells iConsensuscells]=intersect(curepcellinds,consensuscellinds,'rows');
                    consensusspikebins=[consensusspikebins(:,iConsensuscells); curepspikebins(:,iCurcells)];
                end
            end
            
            sleepind1=find(ismember(allsleepanimdayep(:,1:2),curanimday,'rows'));
            
            for cursleepep=1:length(sleepind1)
                
                cursleepcellinds=allsleepspikebins{sleepind1(cursleepep)}.inds(:,[1 2 4 5]);
                
                curSleepSpikeBins=allsleepspikebins{sleepind1(cursleepep)}.spikeBins;
                
                [consensuscellindsSleepWake iCurcells iConsensuscells]=intersect(cursleepcellinds,consensuscellinds,'rows');
                
                curSleepSpikeBins=curSleepSpikeBins(:,iCurcells)';
                curSleepCellInds=allsleepspikebins{sleepind1(cursleepep)}.inds(iCurcells,:);
                
                curWakeSpikeBins=consensusspikebins(:,iConsensuscells)';
                % sanity check
                if size(curSleepSpikeBins,1)~=size(curSleepCellInds,1)| size(curWakeSpikeBins,1)~=size(curSleepCellInds,1)
                    keyboard
                end
                
                % By Peyrache et al. 2009: Replay of rule-learning related neural patterns in the prefrontal cortex during sleep
                
                % Z-transform firing of each cell and rotate so that rows are cells:
                numBins=size(curWakeSpikeBins,2);
                numCells=size(curWakeSpikeBins,1);
                if numCells>=numcellthresh
                
                ZspikeBins=zscore(curWakeSpikeBins')';
                % correlation matrix
                C=(ZspikeBins*ZspikeBins')/numBins; %effectively same as corrcoef(ZspikeBins')
                
                % Put Cii=0 for all i will cancel out the i=j terms:
                C0 = C - diag(diag(C));
                
                
                Qmatch=zscore(curSleepSpikeBins')';
                numMatchBins=size(curSleepSpikeBins,2);
                
                %===== 'FULL MODEL': PRE-PRINCIPAL-COMPONENT DECOMPOSITION
                % I do the same analyses below per- principal component
                MtemplatematchTimeSeries=diag(Qmatch'*C0*Qmatch)';
                Mtemplatematch=(1/(2*numMatchBins))*sum(MtemplatematchTimeSeries);
                
                allalltemplatematchTimeSeries{end+1}=MtemplatematchTimeSeries;
                allalltemplatematch=[allalltemplatematch Mtemplatematch];
                
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
                 % THIS IS WHERE I SUBSTITUTE THE 2 MEASURES ABOVE,
                % CHANGE IF YOU WANT
                currawcrosscorrFullModel=mean(riptrigreactmat)';
                allallcrosscorrsFullModel=[allallcrosscorrsFullModel;currawcrosscorrFullModel'];
                allallripplebins{end+1}=ripvec;
                
                % determining significant cross-corr by shuffling
                limit1=45;
                allshufmaxes=[];
                allmeanshufs=[];
                
                for tt=1:numShuffleRuns
                    riptrigreactmatshuf=[];
                    for qq=1:numrips1
                        shiftamount=round(rand(1)*size(riptrigreactmat,2));
                        riptrigreactmatshuf(qq,:)=circshift(riptrigreactmat(qq,:),shiftamount,2);
                        
                        % riptrigreactmatshuf(qq,:)=riptrigreactmat(qq,randperm(size(riptrigreactmat,2)));
                    end
                    meanshuf=mean(riptrigreactmatshuf);
                    allmeanshufs=[allmeanshufs;meanshuf];
                    %   shufcrosscorr=xcorr(MtemplatematchesTimeSeriesPC(ii,:),ripvec(randperm(length(ripvec))),maxlaginbins,'coeff')';
                    allshufmaxes=[allshufmaxes mean(meanshuf(limit1:end-limit1))];
                end
                
                allallcrosscorrsFullModelShuf=[allallcrosscorrsFullModelShuf;allmeanshufs];

                
                if mean(mean(currawcrosscorrFullModel(limit1:end-limit1))>allshufmaxes)>0.95
                    sigxcorr=1;
                else
                    sigxcorr=0;
                end
                try
                allallanimdaysFullModel=[allallanimdaysFullModel;curSleepCellInds(1,1:3)];
                allallsigxcorrFullModel=[allallsigxcorrFullModel;sigxcorr];
                allallnumcellsFullModel=[allallnumcellsFullModel; numCells];
                allallnumbinsFullModel=[allallnumbinsFullModel; numBins];
                allallnummatchbinsFullModel=[allallnummatchbinsFullModel;numMatchBins];
                catch
                    keyboard
                end
                
                % Alternative loop way, which is less efficient but easier to understand. The result is
                % identical.
                % For each time point, I do outer product of the instantaneous ensemble firing
                % vector, resulting in a matrix where in each i,j is the product of firing
                % of cell i and cell j (in this time bin). If we dot-product this with the
                % correlation matrix (whose diagonal has been 0'd), we get for each pair of
                % cells i,j, their firing product, times their correlation, which is what
                % we want. Since the matric is symmetric we have to divide by 2, and also
                % by the number of bins for normalization (this last one not entirely clear
                % to me but that's just a constant).
                %                 M=[];
                %                 for i=1:numMatchBins
                %                     ensembleOuterProd=Qmatch(:,i)*Qmatch(:,i)';
                %                     M(i)=(1/(2*numMatchBins))*sum(sum(ensembleOuterProd.*C0));
                %                 end
                % % these are already identical
                % update!
                % MtemplatematchTimeSeries=M;
                
                
                % decomposing correlation matrix into eigenvalues and eigenvectors
                [V, D]=eig(C); % now C=V*D*V'
                
                
                % Sorting in descending order, so that now the first is the one explaining
                % most variance
                eigenvals1=diag(D);
                [s, sortind]=sort(eigenvals1,'descend');
                eigenvals=eigenvals1(sortind);
                eigenvectors=V(:,sortind);
                
                % determining significant eigenvalues (as in Peyrache)
                q=numBins/numCells;
                eigvalsigthresh=(1+sqrt(1/q))^2;
                
                sigeigenvals=eigenvals(eigenvals>eigvalsigthresh);
                sigeigenvectors=eigenvectors(:,eigenvals>eigvalsigthresh);
                numsigeigenvals=length(sigeigenvals);
                % These are N matrices, each is the outer product of the i'th eigenvector with itself
                % times the corresponding eigenvalue
                %                 allPs=[];
                %                 for i=1:numCells
                %                     allPs(i,:,:)=eigenvectors(:,i)*eigenvectors(:,i)'*eigenvals(i);
                %                 end
                % as a confirmation, it is also that C is equal to adding up the outer
                % projections of its eigenvectors weighted by the corresponding eigenvalue:
                % cc=squeeze(sum(allPs,1))
                % (now C=cc)
                
                % but here I only look at significant ones:
                allPs=[];
                for i=1:numsigeigenvals
                    curP=sigeigenvectors(:,i)*sigeigenvectors(:,i)'*sigeigenvals(i);
                    % again, I think putting 0's on the diagonal cancels out the i=j terms
                    % doing it for all P matrices
                    curP = curP - diag(diag(curP));
                    allPs(i,:,:)=curP;
                    
                end
                
                
                
                
                % Creating Mtemplatematches, which holds for each time bin (column) the
                % projection on the i'th component (rows)
                MtemplatematchesTimeSeriesPC=[];
                for i=1:numMatchBins
                    ensembleOuterProd=Qmatch(:,i)*Qmatch(:,i)';
                    for j=1:numsigeigenvals
                        MtemplatematchesTimeSeriesPC(j,i)=sum(sum(ensembleOuterProd.*squeeze(allPs(j,:,:))));
                    end
                end
                
              
                
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
                    %
                    
                    
                    allallcrosscorrs{ii}=[allallcrosscorrs{ii};currawcrosscorr'];
                    allallcrosscorrsShufs{ii}=[allallcrosscorrsShufs{ii};allmeanshufs];
                    
                    allallanimdays{ii}=[allallanimdays{ii};curSleepCellInds(1,1:3)];
                    allallsigxcorr{ii}=[allallsigxcorr{ii};sigxcorr];
                    allallnumcells{ii}=[allallnumcells{ii}; numCells];
                    allallnumbins{ii}=[allallnumbins{ii}; numBins];
                    allallnummatchbins{ii}=[allallnummatchbins{ii};numMatchBins];
                    allalltemplatematchesTimeSeriesPC{ii}{end+1}=MtemplatematchesTimeSeriesPC(ii,:);
                    
                end
                
                if plot1
                    % plotting the time series for the global template match and for the
                    % projection on each of the components
                    figure('Position',[100,100,800 1000])
                    subplot(2,1,1)
                    timeseries=0.2:0.2:length(MtemplatematchTimeSeries)*0.2;
                    plot(timeseries,MtemplatematchTimeSeries,'k','linewidth',2)
                    hold on
                    %                 plot(MtemplatematchTimeSeriesesTimeSeries')
                    %    plot(allsleepspikebins{sleepind1}.rippleBins,-0.002,'ro')
                    hold on;plot(timeseries(ripvec==1),max(MtemplatematchTimeSeries),'ro')
                    xlabel('Time (S)')
                    ylabel('reactivation strength')
                    title(num2str(curSleepCellInds(1,1:3)))
                    subplot(2,1,2)
                    crosscorrtimeseries=-4:0.2:4;
                    %                plot(allcrosscorrs)
                    %   plot(crosscorrtimeseries,allcrosscorrs(:,1),'k','linewidth',2)
                    if length(currawcrosscorrFullModel>0)
                        plot(currawcrosscorrFullModel,'k','linewidth',2)
                    end
                    xlabel('Time (S)')
                    ylabel('Cross correlation')
                    
                   keyboard
                    close all
                end
                
                
                
                
            end
            
            end
        end
    end
    
    
    
    % day and epoch correction
    % this maps the sleep epochs to their number in the daily sequence
    % columns: animal,day,epoch,epochnum
    % epochnum is 1/2/3
    animdaysleepepmapping=[1,2,5,1;1,2,7,2;1,2,9,3;1,3,4,1;1,3,6,2;1,3,8,3;1,4,4,1;1,4,6,2;1,4,8,3;1,5,4,1;1,5,6,2;1,5,8,3;1,6,4,1;1,6,6,2;1,6,8,3;1,7,4,1;1,7,6,2;1,7,8,3;1,8,4,1;1,8,6,2;1,8,8,3;1,9,4,1;1,9,6,2;1,9,8,3;1,10,4,1;1,10,6,2;1,10,8,3;1,11,4,1;1,11,6,2;1,11,8,3;1,12,6,2;1,12,8,3;1,13,4,1;1,13,6,2;1,13,8,3;1,14,4,1;1,14,6,2;1,14,8,3;1,15,4,1;1,15,6,2;1,15,8,3;2,1,4,1;2,1,6,2;2,1,8,3;2,2,4,1;2,2,6,2;2,2,8,3;2,3,4,1;2,3,6,2;2,3,8,3;2,4,4,1;2,4,6,2;2,4,8,3;2,5,4,1;2,5,6,2;2,5,8,3;2,6,4,1;2,6,6,2;2,6,8,3;2,7,4,1;2,7,6,2;2,7,8,3;2,8,4,1;2,8,6,2;2,8,8,3;2,9,4,1;2,9,6,2;2,9,8,3;2,10,4,1;2,10,6,2;2,10,8,3;2,11,4,1;2,11,6,2;2,11,8,3;2,12,4,1;2,12,6,2;2,12,8,3;3,8,3,1;3,8,5,2;3,9,3,1;3,9,4,2;3,10,3,1;3,10,5,2;3,11,3,1;3,11,5,2;3,11,7,3;3,12,3,1;3,12,5,2;3,12,7,3;3 13 3 1;3 13 5 2;3 13 7 3;3,14,3,1;3,14,5,2;3,14,7,3;3,15,4,1;3,15,6,2;3,15,8,3;3,16,3,1;3,16,5,2;3,16,7,3;3,17,3,1;3,17,5,2;3,17,7,3;4 1 3 1;4 1 5 2;4 2 3 1;4 2 5 2;4,3,3,1;4,4,3,1;4,4,5,2;4,5,3,1;4,5,5,2;4,6,3,1;4,6,5,2;4,7,3,1;4,7,5,2;4,8,3,1;4,8,5,2;4,8,7,3;4,9,3,1;4,9,5,2;4,10,3,1;4,10,5,2;4,11,3,1;4,11,5,2;4,11,7,3];
    
    % correcting Full Model
    
    tmp1ind=allallanimdaysFullModel;
    correctedind=[];
    for j=1:size(tmp1ind,1)
        whichrow=find(ismember(animdaysleepepmapping(:,1:3),tmp1ind(j,:),'rows'));
        correctedind(j,:)=animdaysleepepmapping(whichrow,[1 2 4]);
        
        % correcting day for Nadal
        if correctedind(j,1)==3
            correctedind(j,2)=correctedind(j,2)-7;
        end
        
    end
    correctedindsFullModel=correctedind;
    correctedepsFullModel=correctedind(:,3);
    
    
    % correcting per-PC
    correctedinds={};
    correctedeps={};
    for i=1:1:3
        tmp1ind=allallanimdays{i};
        correctedind=[];
        for j=1:size(tmp1ind,1)
            whichrow=find(ismember(animdaysleepepmapping(:,1:3),tmp1ind(j,:),'rows'));
            correctedind(j,:)=animdaysleepepmapping(whichrow,[1 2 4]);
            
            % correcting day for Nadal
            if correctedind(j,1)==3
                correctedind(j,2)=correctedind(j,2)-7;
            end
            
        end
        correctedinds{i}=correctedind;
        correctedeps{i}=correctedind(:,3);
        
    end
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
end


toc
return


%%
savedirX = '/opt/data15/gideon/GR_ProcessedData/';

%AllEpochgatherdatafile = [savedirX 'ACnetwork_AllEpochSS_gather'];
AllEpochgatherdatafile = [savedirX 'ACnetwork_AllEpoch_gather'];


% ++ All Epoch
load(AllEpochgatherdatafile)
%load(SlowMovegatherdatafile)

% Full Model Variables
AllEpCorrectedIndsFullModel=correctedindsFullModel;
AllEpCorrectedEpsFullModel=correctedepsFullModel;
AllEpallalltemplatematchTimeSeriesFullModel=allalltemplatematchTimeSeries;
AllEpallalltemplatematchFullModel=allalltemplatematch;
AllEpallallnumcellsFullModel=allallnumcellsFullModel;
AllEpallallcrosscorrsFullModel=allallcrosscorrsFullModel;
AllEpallallcrosscorrsFullModelShuf=allallcrosscorrsFullModelShuf;
AllEpallallsigxcorrFullModel=allallsigxcorrFullModel;
AllEpallallreactivationSignificanceFullModel=allallreactivationSignificanceFullModel;
% PC Model Variables
AllEpCorrectedIndsPC=correctedinds;
AllEpCorrectedEpsPC=correctedeps;
AllEpallalltemplatematchTimeSeriesPC=allalltemplatematchesTimeSeriesPC;
AllEpallallnumcellsPC=allallnumcells;
AllEpallallcrosscorrsPC=allallcrosscorrs;
AllEpallallcrosscorrsPCShuf=allallcrosscorrsShufs;
AllEpallallsigxcorrPC=allallsigxcorr;
AllEpallallreactivationSignifPC=allallreactivationSignifPC;





%%  Fraction of sleep epochs with significant reactivation, according to FULL MODEL
figure;bar([mean(AllEpallallreactivationSignificanceFullModel<0.05) ],'w','linewidth',2)
ylim([0 1])
title('Fraction of sleep epochs with significant reactivation')
set(gca,'XTickLabel',{'All Ep'})
%%  Fraction of sleep epochs with significant reactivation, according to PC MODEL
figure;bar([mean(AllEpallallreactivationSignifPC{1}<0.05) ],'w','linewidth',2)
axis([-1 3 0 1])
title('Fraction of sleep epochs with significant reactivation')
set(gca,'XTickLabel',{'All Ep'})





%%
% Here I create a structure of reactivation significance for the shuffled
% data. As each epoch is shuffled numShuffleRuns times (currently 100), I
% take the significance for the non-shuffled data, and "stretch" it by a
% factor of 100. Eventually this yields
% AllEpallallreactivationSignifPCShuf, which is in the same format as AllEpallallreactivationSignifPC
% but for the shuffled data
AllEpallallreactivationSignifPCShuf={};
for ii=1:5
    tmpsigvecshuf=zeros(size(AllEpallallcrosscorrsPCShuf{ii},1),1);
    for jj=1:length(AllEpallallreactivationSignifPC{ii})
        tmpsigvecshuf((jj-1)*numShuffleRuns+1:jj*numShuffleRuns)=AllEpallallreactivationSignifPC{ii}(jj);
    end
    AllEpallallreactivationSignifPCShuf{ii}=tmpsigvecshuf;

end
%% For Full Model

tmpsigvecshuf=zeros(size(AllEpallallcrosscorrsFullModelShuf,1),1);
for jj=1:length(AllEpallallreactivationSignificanceFullModel)
    tmpsigvecshuf((jj-1)*numShuffleRuns+1:jj*numShuffleRuns)=AllEpallallreactivationSignificanceFullModel(jj);    
end
AllEpallallreactivationSignificanceFullModelShuf=tmpsigvecshuf;

%%
clear a1 e1 s1;
% reactivation-SWR corr by epoch, BY PC
% Choose:
xaxis1=linspace(-10,10,101);
% whether there is significant reactivation strength overall
doonlysigreact=1;
% Currently not using: whether there is a significant cross-corr with SWRs
doonlysig=0;
domove=3;
pc1=1;
tozscore=1;
if doonlysigreact
    if doonlysig
        if domove==1
            % Move, only significantly SWR-corr
            e1=MoveCorrectedEpsPC{pc1}(MoveallallsigxcorrPC{pc1}==1&MoveallallreactivationSignifPC{pc1}'<0.05);
            a1=MoveallallcrosscorrsPC{pc1}(MoveallallsigxcorrPC{pc1}==1&MoveallallreactivationSignifPC{pc1}'<0.05,:);
        elseif domove==2
            % No Move, only significantly SWR-corr
            e1=NoMoveCorrectedEpsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1&NoMoveallallreactivationSignifPC{pc1}'<0.05);
            a1=NoMoveallallcrosscorrsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1&NoMoveallallreactivationSignifPC{pc1}'<0.05,:);
        elseif domove==3
            e1=AllEpCorrectedEpsPC{pc1}(AllEpallallsigxcorrPC{pc1}==1&AllEpallallreactivationSignifPC{pc1}'<0.05);
            a1=AllEpallallcrosscorrsPC{pc1}(AllEpallallsigxcorrPC{pc1}==1&AllEpallallreactivationSignifPC{pc1}'<0.05,:);
            
        end
    else
        if domove==1
            % Move, only significantly SWR-corr
            e1=MoveCorrectedEpsPC{pc1}(MoveallallreactivationSignifPC{pc1}'<0.05);
            a1=MoveallallcrosscorrsPC{pc1}(MoveallallreactivationSignifPC{pc1}'<0.05,:);
        elseif domove==2
            % No Move, only significantly SWR-corr
            e1=NoMoveCorrectedEpsPC{pc1}(NoMoveallallreactivationSignifPC{pc1}'<0.05);
            a1=NoMoveallallcrosscorrsPC{pc1}(NoMoveallallreactivationSignifPC{pc1}'<0.05,:);
        elseif domove==3
            e1=AllEpCorrectedEpsPC{pc1}(AllEpallallreactivationSignifPC{pc1}'<0.05);
            a1=AllEpallallcrosscorrsPC{pc1}(AllEpallallreactivationSignifPC{pc1}'<0.05,:);
            s1=AllEpallallcrosscorrsPCShuf{pc1}(AllEpallallreactivationSignifPCShuf{pc1}'<0.05,:);
        end
        
    end
else
    if doonlysig
        if domove==1
            % Move, only significantly SWR-corr
            e1=MoveCorrectedEpsPC{pc1}(MoveallallsigxcorrPC{pc1}==1);
            a1=MoveallallcrosscorrsPC{pc1}(MoveallallsigxcorrPC{pc1}==1,:);
        elseif domove==2
            % No Move, only significantly SWR-corr
            e1=NoMoveCorrectedEpsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1);
            a1=NoMoveallallcrosscorrsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1,:);
        elseif domove==3
            e1=AllEpCorrectedEpsPC{pc1}(AllEpallallsigxcorrPC{pc1}==1);
            a1=AllEpallallcrosscorrsPC{pc1}(AllEpallallsigxcorrPC{pc1}==1,:);
            
        end
    else
        if domove==1
            % Move, only significantly SWR-corr
            e1=MoveCorrectedEpsPC{pc1};
            a1=MoveallallcrosscorrsPC{pc1};
        elseif domove==2
            % No Move, only significantly SWR-corr
            e1=NoMoveCorrectedEpsPC{pc1};
            a1=NoMoveallallcrosscorrsPC{pc1};
        elseif domove==3
            e1=AllEpCorrectedEpsPC{pc1};
            a1=AllEpallallcrosscorrsPC{pc1};
            s1=AllEpallallcrosscorrsPCShuf{pc1};

            
        end
        
    end
    
end
if tozscore
    a1=zscore(a1')';
    s1=zscore(s1')';
end

ep1=a1(e1==1,:);
ep2=a1(e1==2,:);
ep3=a1(e1==3,:);

% zscored:
% ep1Z=a1Z(e1==1,:);
% ep2Z=a1Z(e1==2,:);
% ep3Z=a1Z(e1==3,:);
% not separating epochs

% subplot(2,2,1)
% plot(xaxis1,nanmean(a1))
% subplot(2,2,3)
% plot(xaxis1,nanmean(a1Z))
%subplot(2,1,1)
figure('Position',[900 100 400 900]);
subplot(2,1,1)
plot(xaxis1,nanmean(a1),'k','linewidth',3)
hold on
plot(xaxis1,nanmean(s1),'r','linewidth',3)
errorbar(xaxis1,nanmean(a1),nanstd(a1)/sqrt(size(a1,1)),'k','linewidth',0.5)
errorbar(xaxis1,nanmean(s1),nanstd(s1)/sqrt(size(s1,1)),'r','linewidth',0.5)

xlim([-4 4])
ylabel('Z-scored reactivation strength')
xlabel('Time (s)')
subplot(2,1,2)
plot(xaxis1,nanmean(a1),'k','linewidth',3)
hold on
plot(xaxis1,nanmean(s1),'r','linewidth',3)
errorbar(xaxis1,nanmean(a1),nanstd(a1)/sqrt(size(a1,1)),'k','linewidth',0.5)
errorbar(xaxis1,nanmean(s1),nanstd(s1)/sqrt(size(s1,1)),'r','linewidth',0.5)

xlim([-1 1])
ylabel('Z-scored reactivation strength')
xlabel('Time (s)')

%
allwin=48:54;
prewin=48:50;
postwin=51:54;

a1zoomed=a1(:,allwin);
s1zoomed=s1(:,allwin);
xaxis1zoomed=xaxis1(allwin);

a1zoomedPre=a1(:,prewin);
s1zoomedPre=s1(:,prewin);
xaxis1zoomedPre=xaxis1(prewin);

a1zoomedPost=a1(:,postwin);
s1zoomedPost=s1(:,postwin);
xaxis1zoomedPost=xaxis1(postwin);
format bank
%%
figure('Position',[100,100,1200 1000])
subplot(3,1,1)
plot(xaxis1zoomed,nanmean(a1zoomed),'k','linewidth',3);
hold on;
plot(xaxis1zoomed,nanmean(s1zoomed),'r','linewidth',3);
errorbar(xaxis1zoomed,nanmean(a1zoomed),nanstd(a1zoomed)/sqrt(size(a1zoomed,1)),'k','linewidth',0.5)
errorbar(xaxis1zoomed,nanmean(s1zoomed),nanstd(s1zoomed)/sqrt(size(s1zoomed,1)),'r','linewidth',0.5)
ylabel('Reactivation strength')
xlabel('Time (s)')
for curbin1=allwin
%    p=ranksum(s1(:,curbin1), a1(:,curbin1));
    p=ranksum(s1(:,curbin1), a1(:,curbin1),'alpha',0.05,'tail','left')

    text(xaxis1(curbin1),1.2,sprintf('%0.4f',p),'fontsize',10 )
end

subplot(3,1,2)
plot(xaxis1zoomedPre,nanmean(a1zoomedPre),'k','linewidth',3)
hold on;
errorbar(xaxis1zoomedPre,nanmean(a1zoomedPre),nanstd(a1zoomedPre)/sqrt(size(a1zoomedPre,1)),'k','linewidth',0.5)
errorbar(xaxis1zoomedPre,nanmean(s1zoomedPre),nanstd(s1zoomedPre)/sqrt(size(s1zoomedPre,1)),'r','linewidth',0.5)

ylabel('Reactivation strength')
xlabel('Time (s)')
subplot(3,1,3)
plot(xaxis1zoomedPost,nanmean(a1zoomedPost),'k','linewidth',3)
hold on;
errorbar(xaxis1zoomedPost,nanmean(a1zoomedPost),nanstd(a1zoomedPost)/sqrt(size(a1zoomedPost,1)),'k','linewidth',0.5)
errorbar(xaxis1zoomedPost,nanmean(s1zoomedPost),nanstd(s1zoomedPost)/sqrt(size(s1zoomedPost,1)),'r','linewidth',0.5)

ylabel('Reactivation strength')
xlabel('Time (s)')
%% p-value at lag 0
p=ranksum(s1(:,51), a1(:,51),'alpha',0.05,'tail','left')
%% per epoch
figure('Position',[900 100 400 900]);
plot(xaxis1,nanmean(ep1),'g','linewidth',3)
hold on
plot(xaxis1,nanmean(ep2),'b','linewidth',3)
plot(xaxis1,nanmean(ep3),'m','linewidth',3)

plot(xaxis1,nanmean(s1),'r','linewidth',3)
errorbar(xaxis1,nanmean(ep1),nanstd(ep1)/sqrt(size(ep1,1)),'g','linewidth',0.5)
errorbar(xaxis1,nanmean(ep2),nanstd(ep2)/sqrt(size(ep2,1)),'b','linewidth',0.5)
errorbar(xaxis1,nanmean(ep3),nanstd(ep3)/sqrt(size(ep3,1)),'m','linewidth',0.5)

errorbar(xaxis1,nanmean(s1),nanstd(s1)/sqrt(size(s1,1)),'r','linewidth',0.5)

xlim([-4 4])
ylabel('Reactivation strength')
xlabel('Time (s)')


%% FULL MODEL
%
clear a1 e1 s1;
% reactivation-SWR corr by epoch, BY PC
% Choose:
xaxis1=linspace(-10,10,101);
% whether there is significant reactivation strength overall
doonlysigreact=1;
% Currently not using: whether there is a significant cross-corr with SWRs
doonlysig=0;
domove=3;
pc1=1;
tozscore=1;
if doonlysigreact
    if doonlysig
%         if domove==1
%             % Move, only significantly SWR-corr
%             e1=MoveCorrectedEpsPC{pc1}(MoveallallsigxcorrPC{pc1}==1&MoveallallreactivationSignifPC{pc1}'<0.05);
%             a1=MoveallallcrosscorrsPC{pc1}(MoveallallsigxcorrPC{pc1}==1&MoveallallreactivationSignifPC{pc1}'<0.05,:);
%         elseif domove==2
%             % No Move, only significantly SWR-corr
%             e1=NoMoveCorrectedEpsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1&NoMoveallallreactivationSignifPC{pc1}'<0.05);
%             a1=NoMoveallallcrosscorrsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1&NoMoveallallreactivationSignifPC{pc1}'<0.05,:);
%         elseif domove==3
        %    e1=AllEpCorrectedEpsFullModel(AllEpallallsigxcorrFullModel==1&AllEpallallreactivationSignificanceFullModel'<0.05);
        %    a1=AllEpallallcrosscorrsFullModel(AllEpallallsigxcorrFullModel==1&AllEpallallreactivationSignificanceFullModel'<0.05,:);
            
        %end
    else
%         if domove==1
%             % Move, only significantly SWR-corr
%             e1=MoveCorrectedEpsPC{pc1}(MoveallallreactivationSignifPC{pc1}'<0.05);
%             a1=MoveallallcrosscorrsPC{pc1}(MoveallallreactivationSignifPC{pc1}'<0.05,:);
%         elseif domove==2
%             % No Move, only significantly SWR-corr
%             e1=NoMoveCorrectedEpsPC{pc1}(NoMoveallallreactivationSignifPC{pc1}'<0.05);
%             a1=NoMoveallallcrosscorrsPC{pc1}(NoMoveallallreactivationSignifPC{pc1}'<0.05,:);
%         elseif domove==3
            e1=AllEpCorrectedEpsFullModel(AllEpallallreactivationSignificanceFullModel'<0.05);
            a1=AllEpallallcrosscorrsFullModel(AllEpallallreactivationSignificanceFullModel'<0.05,:);
            s1=AllEpallallcrosscorrsFullModelShuf(AllEpallallreactivationSignificanceFullModelShuf'<0.05,:);
            
      %  end
        
    end
else
    if doonlysig
%         if domove==1
%             % Move, only significantly SWR-corr
%             e1=MoveCorrectedEpsPC{pc1}(MoveallallsigxcorrPC{pc1}==1);
%             a1=MoveallallcrosscorrsPC{pc1}(MoveallallsigxcorrPC{pc1}==1,:);
%         elseif domove==2
%             % No Move, only significantly SWR-corr
%             e1=NoMoveCorrectedEpsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1);
%             a1=NoMoveallallcrosscorrsPC{pc1}(NoMoveallallsigxcorrPC{pc1}==1,:);
%         elseif domove==3
    %        e1=AllEpCorrectedEpsPC{pc1}(AllEpallallsigxcorrPC{pc1}==1);
     %       a1=AllEpallallcrosscorrsPC{pc1}(AllEpallallsigxcorrPC{pc1}==1,:);
            
      %  end
    else
%         if domove==1
%             % Move, only significantly SWR-corr
%             e1=MoveCorrectedEpsPC{pc1};
%             a1=MoveallallcrosscorrsPC{pc1};
%         elseif domove==2
%             % No Move, only significantly SWR-corr
%             e1=NoMoveCorrectedEpsPC{pc1};
%             a1=NoMoveallallcrosscorrsPC{pc1};
%         elseif domove==3
       

            e1=AllEpCorrectedEpsFullModel;
            a1=AllEpallallcrosscorrsFullModel;
            s1=AllEpallallcrosscorrsFullModelShuf;
    %    end
        
    end
    
end
if tozscore
    a1=zscore(a1')';
    s1=zscore(s1')';
end

ep1=a1(e1==1,:);
ep2=a1(e1==2,:);
ep3=a1(e1==3,:);

% zscored:
% ep1Z=a1Z(e1==1,:);
% ep2Z=a1Z(e1==2,:);
% ep3Z=a1Z(e1==3,:);
% not separating epochs

% subplot(2,2,1)
% plot(xaxis1,nanmean(a1))
% subplot(2,2,3)
% plot(xaxis1,nanmean(a1Z))
%subplot(2,1,1)
figure
subplot(2,1,1)
plot(xaxis1,nanmean(a1),'k','linewidth',3)
hold on
plot(xaxis1,nanmean(s1),'r','linewidth',3)
errorbar(xaxis1,nanmean(a1),nanstd(a1)/sqrt(size(a1,1)),'k','linewidth',0.5)
errorbar(xaxis1,nanmean(s1),nanstd(s1)/sqrt(size(s1,1)),'r','linewidth',0.5)

xlim([-4 4])
ylabel('Z-scored reactivation strength')
xlabel('Time (s)')
subplot(2,1,2)
plot(xaxis1,nanmean(a1),'k','linewidth',3)
hold on
plot(xaxis1,nanmean(s1),'r','linewidth',3)
errorbar(xaxis1,nanmean(a1),nanstd(a1)/sqrt(size(a1,1)),'k','linewidth',0.5)
errorbar(xaxis1,nanmean(s1),nanstd(s1)/sqrt(size(s1,1)),'r','linewidth',0.5)

xlim([-1 1])
ylabel('Z-scored reactivation strength')
xlabel('Time (s)')

%
allwin=48:54;
prewin=48:50;
postwin=51:54;

a1zoomed=a1(:,allwin);
s1zoomed=s1(:,allwin);
xaxis1zoomed=xaxis1(allwin);

a1zoomedPre=a1(:,prewin);
s1zoomedPre=s1(:,prewin);
xaxis1zoomedPre=xaxis1(prewin);

a1zoomedPost=a1(:,postwin);
s1zoomedPost=s1(:,postwin);
xaxis1zoomedPost=xaxis1(postwin);
format bank
figure('Position',[100,100,1200 1000])
subplot(3,1,1)
plot(xaxis1zoomed,nanmean(a1zoomed),'k','linewidth',3);
hold on;
plot(xaxis1zoomed,nanmean(s1zoomed),'r','linewidth',3);
errorbar(xaxis1zoomed,nanmean(a1zoomed),nanstd(a1zoomed)/sqrt(size(a1zoomed,1)),'k','linewidth',0.5)
errorbar(xaxis1zoomed,nanmean(s1zoomed),nanstd(s1zoomed)/sqrt(size(s1zoomed,1)),'r','linewidth',0.5)
ylabel('Reactivation strength')
xlabel('Time (s)')
for curbin1=allwin
    p=ranksum(s1(:,curbin1), a1(:,curbin1));
    text(xaxis1(curbin1),1.2,sprintf('%0.4f',p),'fontsize',10 )
end

subplot(3,1,2)
plot(xaxis1zoomedPre,nanmean(a1zoomedPre),'k','linewidth',3)
hold on;
errorbar(xaxis1zoomedPre,nanmean(a1zoomedPre),nanstd(a1zoomedPre)/sqrt(size(a1zoomedPre,1)),'k','linewidth',0.5)
errorbar(xaxis1zoomedPre,nanmean(s1zoomedPre),nanstd(s1zoomedPre)/sqrt(size(s1zoomedPre,1)),'r','linewidth',0.5)

ylabel('Reactivation strength')
xlabel('Time (s)')
subplot(3,1,3)
plot(xaxis1zoomedPost,nanmean(a1zoomedPost),'k','linewidth',3)
hold on;
errorbar(xaxis1zoomedPost,nanmean(a1zoomedPost),nanstd(a1zoomedPost)/sqrt(size(a1zoomedPost,1)),'k','linewidth',0.5)
errorbar(xaxis1zoomedPost,nanmean(s1zoomedPost),nanstd(s1zoomedPost)/sqrt(size(s1zoomedPost,1)),'r','linewidth',0.5)

ylabel('Reactivation strength')
xlabel('Time (s)')

%% plot example
xaxis1=linspace(-10,10,101);

%bb=gaussian(1,3);
% examples: 66,39,47,65,44
ensemblenum=39;

figure('Position',[200 100 800 400]);
curtrace=allalltemplatematchesTimeSeriesPC{1}{ensemblenum};
timeseries1=0.2:0.2:length(curtrace)*0.2;
subplot(1,5,1:4);
plot(timeseries1,curtrace,'k','linewidth',2);
hold on
xval=[timeseries1(allallripplebins{ensemblenum}>0);timeseries1(allallripplebins{ensemblenum}>0)]
yval=[50*ones(1,length(xval));52*ones(1,length(xval))]
plot(xval,yval,'r','linewidth',2)

%plot(timeseries1(allallripplebins{ensemblenum}>0),60,'ro')
ylim([-10 180])
ylabel('Reactivation strength')
xlabel('Time (S)')

subplot(1,5,5);
%plot(xaxis1,filtfilt(bb,1,AllEpallallcrosscorrsPC{1}(ensemblenum,:)),'r','linewidth',2)
plot(xaxis1,zscore(AllEpallallcrosscorrsPC{1}(ensemblenum,:)')','r','linewidth',2)

xlim([-4 4])
ylabel('Z-scored reactivation')
xlabel('Time rel. to SWRs')

%% hist
[n1 x1]=hist(curtrace,100);
figure
bar(x1,n1/sum(n1),'k')
set(gca,'YScale','log')
xlabel('Reactivation strength')
title(['Skewness= ' num2str(skewness(curtrace))])
%% skewness values across experiments
allskewnessvals=[];
for ensemblenum1=1:length(allalltemplatematchesTimeSeriesPC{1})
    curtrace=allalltemplatematchesTimeSeriesPC{1}{ensemblenum1};
    if AllEpallallreactivationSignifPC{1}(ensemblenum1)<0.05
        allskewnessvals=[allskewnessvals skewness(curtrace)];
    end
end
figure;
hist(allskewnessvals)
xlabel('Skewness')
ylabel('Count')
